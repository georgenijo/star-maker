"""Load and cache the star catalog from data/stars.json."""
from pathlib import Path
import json
from typing import TypedDict


class Star(TypedDict):
    ra: float       # Right ascension in degrees (0-360)
    dec: float      # Declination in degrees (-90 to +90)
    mag: float      # Apparent visual magnitude
    hip: int | None
    proper: str | None
    con: str | None


_DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "stars.json"
_cache: list[Star] | None = None


def load_stars() -> list[Star]:
    """Load stars from data/stars.json. Cached after first call.

    Returns list of Star dicts, sorted by magnitude (brightest first).
    Raises FileNotFoundError if data/stars.json does not exist.
    """
    global _cache
    if _cache is None:
        with open(_DATA_PATH) as f:
            raw = json.load(f)
        _cache = [s for s in raw if s.get("proper") != "Sol"]
    return _cache
