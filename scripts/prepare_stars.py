#!/usr/bin/env python3
"""One-off script: download HYG v4.1 CSV, filter to naked-eye stars, write data/stars.json.

Usage: python scripts/prepare_stars.py
"""
import csv
import json
import os
import urllib.request
from pathlib import Path

HYG_URL = "https://raw.githubusercontent.com/astronexus/HYG-Database/refs/heads/main/hyg/CURRENT/hygdata_v41.csv"
RAW_CSV = Path(__file__).parent / "hygdata_v41.csv"
OUTPUT_JSON = Path(__file__).resolve().parent.parent / "data" / "stars.json"
MAG_LIMIT = 6.5


def download_csv() -> Path:
    """Download HYG CSV if not already cached locally. Returns path."""
    if not RAW_CSV.exists():
        print(f"Downloading {HYG_URL}...")
        urllib.request.urlretrieve(HYG_URL, RAW_CSV)
        print(f"Saved to {RAW_CSV}")
    else:
        print(f"Using cached {RAW_CSV}")
    return RAW_CSV


def parse_and_filter(csv_path: Path) -> list[dict]:
    """Read CSV, filter mag <= 6.5, extract fields.

    CRITICAL: HYG 'ra' column is in HOURS (0-24).
    Convert to DEGREES by multiplying by 15.
    Store as 'ra' in degrees in the output JSON.

    Returns list of dicts with keys: ra, dec, mag, hip, proper, con
    """
    stars = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mag_str = row.get("mag", "")
            if not mag_str:
                continue
            mag = float(mag_str)
            if mag > MAG_LIMIT:
                continue

            ra_hours = float(row["ra"])
            ra_deg = ra_hours * 15.0       # CRITICAL CONVERSION

            hip_str = row.get("hip", "")
            star = {
                "ra": round(ra_deg, 4),
                "dec": round(float(row["dec"]), 4),
                "mag": round(mag, 2),
                "hip": int(hip_str) if hip_str else None,
                "proper": row.get("proper") or None,
                "con": row.get("con") or None,
            }
            stars.append(star)

    stars.sort(key=lambda s: s["mag"])  # brightest first
    return stars


def write_json(stars: list[dict], output_path: Path) -> None:
    """Write star list to compact JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(stars, f, separators=(",", ":"))
    size_kb = output_path.stat().st_size / 1024
    print(f"Wrote {len(stars)} stars to {output_path} ({size_kb:.0f} KB)")


def main() -> None:
    csv_path = download_csv()
    stars = parse_and_filter(csv_path)
    write_json(stars, OUTPUT_JSON)


if __name__ == "__main__":
    main()
