from __future__ import annotations

import gzip
import struct
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import requests


@dataclass(frozen=True)
class TabixReference:
    bins: dict[int, list[tuple[int, int]]]
    linear_offsets: list[int]


def _read_i32(data: bytes) -> int:
    return struct.unpack("<i", data)[0]


def _read_u32(data: bytes) -> int:
    return struct.unpack("<I", data)[0]


def _read_u64(data: bytes) -> int:
    return struct.unpack("<Q", data)[0]


def fetch_bytes(url: str, start: int, end: int, timeout: int = 120) -> bytes:
    response = requests.get(
        url,
        headers={"Range": f"bytes={start}-{end}"},
        timeout=timeout,
    )
    if response.status_code not in {200, 206}:
        response.raise_for_status()
    return response.content


def load_tabix_index(index_path: Path, index_url: str) -> tuple[list[str], list[TabixReference]]:
    if index_path.exists():
        raw = index_path.read_bytes()
    else:
        response = requests.get(index_url, timeout=120)
        response.raise_for_status()
        raw = response.content
        index_path.parent.mkdir(parents=True, exist_ok=True)
        index_path.write_bytes(raw)
    data = gzip.decompress(raw)
    if data[:4] != b"TBI\x01":
        raise ValueError(f"Not a TBI index: {index_url}")
    offset = 4
    n_ref = _read_i32(data[offset : offset + 4])
    offset += 4
    offset += 6 * 4
    l_nm = _read_i32(data[offset : offset + 4])
    offset += 4
    names_blob = data[offset : offset + l_nm]
    offset += l_nm
    names = names_blob.rstrip(b"\x00").decode("utf-8").split("\x00")
    references = []
    for _ in range(n_ref):
        n_bin = _read_i32(data[offset : offset + 4])
        offset += 4
        bins = {}
        for _ in range(n_bin):
            bin_id = _read_u32(data[offset : offset + 4])
            offset += 4
            n_chunk = _read_i32(data[offset : offset + 4])
            offset += 4
            chunks = []
            for _ in range(n_chunk):
                beginning = _read_u64(data[offset : offset + 8])
                offset += 8
                ending = _read_u64(data[offset : offset + 8])
                offset += 8
                chunks.append((beginning, ending))
            bins[bin_id] = chunks
        n_intv = _read_i32(data[offset : offset + 4])
        offset += 4
        linear_offsets = []
        for _ in range(n_intv):
            linear_offsets.append(_read_u64(data[offset : offset + 8]))
            offset += 8
        references.append(TabixReference(bins, linear_offsets))
    return names, references


def reg2bins(begin: int, end: int) -> list[int]:
    end -= 1
    bins = [0]
    for shift, offset in [
        (26, 1),
        (23, 9),
        (20, 73),
        (17, 585),
        (14, 4681),
    ]:
        start_bin = offset + (begin >> shift)
        end_bin = offset + (end >> shift)
        bins.extend(range(start_bin, end_bin + 1))
    return bins


def merge_chunks(chunks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not chunks:
        return []
    chunks = sorted(chunks)
    merged = [chunks[0]]
    for beginning, ending in chunks[1:]:
        old_beginning, old_ending = merged[-1]
        if beginning <= old_ending:
            merged[-1] = (old_beginning, max(old_ending, ending))
        else:
            merged.append((beginning, ending))
    return merged


def bgzf_blocks(data: bytes) -> list[tuple[int, bytes]]:
    blocks = []
    for line in data.splitlines():
        if not line:
            continue
        fields = line.split(b"\t")
        if len(fields) >= 5:
            blocks.append((int(fields[0]), fields[4]))
    return blocks


def read_bgzf_range(
    data_url: str,
    start_virtual: int,
    end_virtual: int,
) -> bytes:
    start_coffset = start_virtual >> 16
    start_uoffset = start_virtual & 0xFFFF
    end_coffset = end_virtual >> 16
    end_uoffset = end_virtual & 0xFFFF
    compressed = fetch_bytes(
        data_url,
        start_coffset,
        end_coffset + 65535,
    )
    position = 0
    current_coffset = start_coffset
    pieces: list[bytes] = []
    while current_coffset <= end_coffset:
        if position >= len(compressed):
            compressed += fetch_bytes(
                data_url,
                start_coffset + position,
                end_coffset + 65535,
            )
            position = 0
            start_coffset += len(compressed)
        if position + 18 > len(compressed):
            more = fetch_bytes(
                data_url,
                current_coffset,
                end_coffset + 65535,
            )
            compressed = more
            position = 0
        block_size = struct.unpack("<H", compressed[position + 16 : position + 18])[0] + 1
        if position + block_size > len(compressed):
            compressed = fetch_bytes(
                data_url,
                current_coffset,
                current_coffset + block_size - 1,
            )
            position = 0
        block = compressed[position : position + block_size]
        uncompressed = gzip.decompress(block)
        low = start_uoffset if current_coffset == (start_virtual >> 16) else 0
        high = end_uoffset if current_coffset == end_coffset else len(uncompressed)
        pieces.append(uncompressed[low:high])
        position += block_size
        current_coffset += block_size
    return b"".join(pieces)


def query_region(
    data_url: str,
    index_url: str,
    chromosome: str,
    start: int,
    end: int,
    cache_dir: Path,
) -> list[str]:
    host = urlparse(index_url).netloc.replace(":", "_")
    index_path = cache_dir / f"{host}_{Path(urlparse(index_url).path).name}"
    names, references = load_tabix_index(index_path, index_url)
    normalized = chromosome.removeprefix("chr")
    candidates = [normalized, f"chr{normalized}"]
    reference_id = next(
        (index for index, name in enumerate(names) if name in candidates),
        None,
    )
    if reference_id is None:
        return []
    reference = references[reference_id]
    begin = max(0, start - 1)
    end = max(begin + 1, end)
    # Point SNV queries only need the deepest 16-kb bins. Including all
    # parent bins can pull chromosome-scale chunks and defeat HTTP range access.
    deep_bins = range(
        4681 + (begin >> 14),
        4681 + ((end - 1) >> 14) + 1,
    )
    chunks = []
    for bin_id in deep_bins:
        chunks.extend(reference.bins.get(bin_id, []))
    if not chunks:
        for bin_id in reg2bins(begin, end):
            chunks.extend(reference.bins.get(bin_id, []))
    linear_index = begin >> 14
    minimum_offset = (
        reference.linear_offsets[linear_index]
        if linear_index < len(reference.linear_offsets)
        else 0
    )
    chunks = [chunk for chunk in chunks if chunk[1] > minimum_offset]
    lines = []
    for chunk_begin, chunk_end in merge_chunks(chunks):
        raw = read_bgzf_range(data_url, chunk_begin, chunk_end)
        for line in raw.decode("utf-8", "replace").splitlines():
            fields = line.split("\t")
            chromosome_index = next(
                (
                    index
                    for index, value in enumerate(fields[:5])
                    if value.removeprefix("chr") == normalized
                    and index + 1 < len(fields)
                ),
                None,
            )
            if chromosome_index is None:
                continue
            try:
                line_position = int(fields[chromosome_index + 1])
            except ValueError:
                continue
            if begin <= line_position <= end:
                lines.append(line)
    return lines


def query_point(
    data_url: str,
    index_url: str,
    chromosome: str,
    position: int,
    cache_dir: Path,
) -> list[str]:
    return query_region(
        data_url,
        index_url,
        chromosome,
        position,
        position,
        cache_dir,
    )
