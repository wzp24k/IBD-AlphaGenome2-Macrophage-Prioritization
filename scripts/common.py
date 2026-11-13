from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

import requests


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "analysis.json"


def load_config() -> dict[str, Any]:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_api_key() -> str:
    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if key:
        return key
    key_path = ROOT / "API.txt"
    if not key_path.exists():
        raise FileNotFoundError(
            "Set ALPHAGENOME_API_KEY or place the key in API.txt."
        )
    key = key_path.read_text(encoding="utf-8").strip()
    if not key:
        raise ValueError("API.txt is empty.")
    return key


def configure_proxy() -> None:
    if os.name != "nt":
        return
    if any(os.environ.get(name) for name in ("HTTPS_PROXY", "https_proxy")):
        return
    try:
        import winreg

        location = r"Software\Microsoft\Windows\CurrentVersion\Internet Settings"
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, location) as key:
            enabled, _ = winreg.QueryValueEx(key, "ProxyEnable")
            if not enabled:
                return
            proxy, _ = winreg.QueryValueEx(key, "ProxyServer")
    except (FileNotFoundError, OSError):
        return
    if not proxy:
        return
    if "=" in proxy:
        parts = dict(item.split("=", 1) for item in proxy.split(";") if "=" in item)
        proxy = parts.get("https") or parts.get("http") or ""
    if proxy and "://" not in proxy:
        proxy = f"http://{proxy}"
    if proxy:
        os.environ.setdefault("HTTPS_PROXY", proxy)
        os.environ.setdefault("HTTP_PROXY", proxy)
        os.environ.setdefault("https_proxy", proxy)
        os.environ.setdefault("http_proxy", proxy)


def post_graphql(
    session: requests.Session,
    endpoint: str,
    query: str,
    variables: dict[str, Any] | None = None,
    *,
    attempts: int = 4,
    timeout: int = 180,
) -> dict[str, Any]:
    payload: dict[str, Any] = {"query": query}
    if variables is not None:
        payload["variables"] = variables
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = session.post(endpoint, json=payload, timeout=timeout)
            response.raise_for_status()
            data = response.json()
            if data.get("errors"):
                raise RuntimeError("; ".join(e.get("message", "") for e in data["errors"]))
            return data
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"GraphQL request failed: {last_error}")
