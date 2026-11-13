from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key


CHUNK_SIZE = 40
INTESTINE_PATTERN = re.compile(
    r"colon|intest|ileum|caecum|rectum|colonic|mucosa",
    re.IGNORECASE,
)
MACROPHAGE_PATTERN = re.compile(r"macrophage|monocyte|myeloid", re.IGNORECASE)


def aggregate_score_rows(
    adata: ad.AnnData,
    values: np.ndarray,
    quantiles: np.ndarray,
) -> list[dict]:
    track = adata.var.copy()
    names = track["biosample_name"].fillna("").astype(str).to_numpy()
    intestine_mask = np.array([bool(INTESTINE_PATTERN.search(x)) for x in names])
    myeloid_mask = np.array([bool(MACROPHAGE_PATTERN.search(x)) for x in names])
    obs = adata.obs.reset_index(drop=True)
    records: list[dict] = []

    for row_index, obs_row in obs.iterrows():
        raw = np.asarray(values[row_index], dtype=float)
        q = np.asarray(quantiles[row_index], dtype=float)
        finite = np.isfinite(raw) & np.isfinite(q)
        if not finite.any():
            continue

        abs_raw = np.where(finite, np.abs(raw), -np.inf)
        raw_idx = int(np.argmax(abs_raw))
        abs_q = np.where(finite, np.abs(q), -np.inf)
        q_idx = int(np.argmax(abs_q))
        record = {
            "variant_id": str(obs_row.get("variant")),
            "gene_id": str(obs_row.get("gene_id", "")).split(".")[0],
            "gene_name": obs_row.get("gene_name"),
            "gene_strand": obs_row.get("strand"),
            "raw_at_max_abs_raw": float(raw[raw_idx]),
            "q_at_max_abs_raw": float(q[raw_idx]),
            "track_at_max_abs_raw": track.iloc[raw_idx].get("name"),
            "biosample_at_max_abs_raw": names[raw_idx],
            "max_abs_q": float(abs_q[q_idx]),
            "q_at_max_abs_q": float(q[q_idx]),
            "track_at_max_abs_q": track.iloc[q_idx].get("name"),
            "biosample_at_max_abs_q": names[q_idx],
            "colon_intestine_max_abs_q": np.nan,
            "colon_intestine_q_at_max_abs_q": np.nan,
            "colon_intestine_biosample_at_max_abs_q": None,
            "myeloid_max_abs_q": np.nan,
            "myeloid_q_at_max_abs_q": np.nan,
            "myeloid_biosample_at_max_abs_q": None,
        }

        for prefix, mask in (
            ("colon_intestine", intestine_mask),
            ("myeloid", myeloid_mask),
        ):
            subset = finite & mask
            if not subset.any():
                continue
            subset_abs_q = np.where(subset, np.abs(q), -np.inf)
            idx = int(np.argmax(subset_abs_q))
            record[f"{prefix}_max_abs_q"] = float(abs(q[idx]))
            record[f"{prefix}_q_at_max_abs_q"] = float(q[idx])
            record[f"{prefix}_biosample_at_max_abs_q"] = names[idx]
        records.append(record)
    return records


def query_one_chunk(client, variants: list, scratch_path: Path) -> tuple[list[dict], list[str]]:
    scores = client.query_variants(
        variants,
        requested_scorers=["RNA_SEQ"],
        progress_bar=False,
        max_workers=4,
    )
    adata = scores["RNA_SEQ"]
    values = np.asarray(adata.X)
    quantiles = np.asarray(adata.layers["quantiles"])
    records = aggregate_score_rows(adata, values, quantiles)
    pd.DataFrame(records).to_csv(scratch_path, index=False)
    return records, []


def query_with_isolation(client, variants: list, scratch_dir: Path) -> tuple[list[dict], list[str]]:
    records: list[dict] = []
    invalid: list[str] = []
    for offset in range(0, len(variants), CHUNK_SIZE):
        chunk = variants[offset : offset + CHUNK_SIZE]
        chunk_hash = hashlib.sha1(
            "|".join(str(value) for value in chunk).encode("ascii")
        ).hexdigest()[:12]
        chunk_path = scratch_dir / f"atlas_{chunk_hash}_{len(chunk)}.csv"
        if chunk_path.exists():
            records.extend(pd.read_csv(chunk_path).to_dict("records"))
            continue
        try:
            chunk_records, _ = query_one_chunk(client, chunk, chunk_path)
            records.extend(chunk_records)
            print(f"Atlas: completed {min(offset + len(chunk), len(variants))}/{len(variants)}")
        except Exception:
            if len(chunk) == 1:
                invalid.append(str(chunk[0]))
                continue
            midpoint = len(chunk) // 2
            left, left_invalid = query_with_isolation(
                client, chunk[:midpoint], scratch_dir
            )
            right, right_invalid = query_with_isolation(
                client, chunk[midpoint:], scratch_dir
            )
            records.extend(left + right)
            invalid.extend(left_invalid + right_invalid)
    return records, invalid


def main() -> None:
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    configure_proxy()
    key = load_api_key()
    client = atlas.create(key, timeout=30)
    input_path = ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    frame = pd.read_csv(input_path)
    primary = frame[frame["study_role"] == "primary_ibd"].drop_duplicates(
        "canonical_variant"
    )
    variant_strings = primary["canonical_variant"].astype(str).tolist()
    variants = [genome.Variant.from_str(value) for value in variant_strings]

    scratch_dir = ROOT / "work" / "atlas_chunks"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    records, invalid = query_with_isolation(client, variants, scratch_dir)

    output_dir = ROOT / "results" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.DataFrame(records)
    scores.to_csv(output_dir / "atlas_rna_seq_gene_quantiles.csv", index=False)
    pd.DataFrame({"variant_id": invalid}).to_csv(
        output_dir / "atlas_invalid_variants.csv", index=False
    )

    top_columns = [
        "variant_id",
        "gene_id",
        "gene_name",
        "raw_at_max_abs_raw",
        "q_at_max_abs_raw",
        "max_abs_q",
        "q_at_max_abs_q",
        "biosample_at_max_abs_q",
        "colon_intestine_max_abs_q",
        "colon_intestine_q_at_max_abs_q",
        "colon_intestine_biosample_at_max_abs_q",
    ]
    scores["abs_signed_q"] = scores["q_at_max_abs_q"].abs()
    top = (
        scores.sort_values(["variant_id", "abs_signed_q"], ascending=[True, False])
        .groupby("variant_id", as_index=False)
        .head(5)[top_columns]
    )
    top.to_csv(output_dir / "atlas_variant_top_genes.csv", index=False)

    print(f"Atlas scored variants: {len(primary) - len(invalid)}")
    print(f"Atlas invalid variants: {len(invalid)}")
    print(f"Gene-level rows: {len(scores)}")


if __name__ == "__main__":
    main()
