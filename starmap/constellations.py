"""Load and parse constellation line data from d3-celestial GeoJSON."""
from pathlib import Path
import json
from typing import Iterator

_DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "constellations.json"
_cache: dict | None = None


def load_constellations() -> dict:
    """Load raw GeoJSON FeatureCollection. Cached after first call.

    Returns the full GeoJSON dict with 'type' and 'features' keys.
    """
    global _cache
    if _cache is None:
        with open(_DATA_PATH) as f:
            _cache = json.load(f)
    return _cache


def iter_line_segments() -> Iterator[tuple[float, float, float, float]]:
    """Yield all constellation line segments as (ra1, dec1, ra2, dec2).

    Each segment is a pair of consecutive points within a polyline.
    A constellation's MultiLineString contains multiple polylines,
    each with N points producing N-1 segments.

    Coordinates are degrees. RA may be negative (-180 to +180).
    """
    data = load_constellations()
    for feature in data["features"]:
        for polyline in feature["geometry"]["coordinates"]:
            for i in range(len(polyline) - 1):
                ra1, dec1 = polyline[i]
                ra2, dec2 = polyline[i + 1]
                yield (ra1, dec1, ra2, dec2)
