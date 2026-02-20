# MASTER_PLAN.md — Star Map Generator Implementation Blueprint

> **Purpose**: Single-source implementation spec. A main Claude Code session reads this plan,
> dispatches subagents for each task, runs gate tests, and moves to the next task.
> The main session **never writes code itself** — it only orchestrates.

---

## Table of Contents

1. [Orchestration Model](#1-orchestration-model)
2. [Dependency Map & Build Order](#2-dependency-map--build-order)
3. [Phase 1: Data Pipeline + Projection Math](#3-phase-1-data-pipeline--projection-math)
4. [Phase 2: SVG Renderer](#4-phase-2-svg-renderer)
5. [Phase 3: Flask Application + Web UI](#5-phase-3-flask-application--web-ui)
6. [Phase 4: Polish + Edge Cases](#6-phase-4-polish--edge-cases)
7. [GitHub Issues Reference](#7-github-issues-reference)
8. [Data Format Reference](#8-data-format-reference)

---

## 1. Orchestration Model

### Main Session Loop

```
WHILE GitHub issues remain open:
    1. issues = `gh issue list --milestone "Phase N" --state open`
    2. FOR each issue (lowest number first):
         - Read dependency list from issue body
         - IF all dependency issues are closed → dispatch subagent
    3. WAIT for subagent to complete
    4. RUN smoke test from issue body
    5. IF smoke test passes → `gh issue close <number>`
       ELSE → re-dispatch subagent with error output as additional context
    6. IF all issues in milestone are closed:
         - RUN phase gate test
         - IF passes → move to Phase N+1
         - ELSE → create bug-fix issue, re-dispatch
```

### Subagent Prompt Template

Feed this to every subagent (fill in the blanks):

```
You are implementing a single task for the star-maker project.

## Project Context
<paste CLAUDE.md contents or tell subagent to read it>

## Your Task
<paste GitHub issue body verbatim>

## Files That Already Exist (read-only context)
<list dependency files that exist>

## Files You Must Create or Modify
<exact file paths from the issue>

## Critical Conventions
- RA in data/stars.json is in DEGREES (already converted from hours × 15)
- d3-celestial RA can be negative (-180 to +180) — Astropy handles it fine
- Use `zoneinfo` (stdlib), NOT `pytz`
- Use `xml.sax.saxutils.escape()` for any user text embedded in SVG
- SVG text: UPPERCASE + letter-spacing, NOT font-variant: small-caps
- Optimize constellation rendering: batch all endpoints into arrays, transform ONCE
- Path resolution: use `Path(__file__).resolve().parent`, never hardcoded absolute paths
- All functions must have type hints and docstrings

## After Implementation
Run the smoke test listed in the issue and confirm it passes.
Commit your changes with a clear message.
```

### Parallelization Opportunities

```
Phase 1:
  Group A (no deps):      Issue #1 (scaffolding) + Issue #2 (constellation download)
  Group B (needs Group A): Issue #3 (star prep) + Issue #5 (projection.py) — parallel
  Group C (needs Group A + B): Issue #4 (catalog.py) + Issue #6 (constellations.py) — parallel
  Group D (needs all above): Issue #7 (gate test)

Phase 2:
  Sequential: Issue #8 (renderer) → Issue #9 (gate test)

Phase 3:
  Issue #10 (app.py) first, then Issues #11 + #12 (HTML + CSS) in parallel, then Issue #13 (gate)

Phase 4:
  Issues #14 + #15 + #17 in parallel, then Issue #16 (code quality), then final gate
```

### Error Recovery

If a subagent's code fails the smoke test:
1. Capture the error output (traceback, assertion failure)
2. Re-dispatch the **same** subagent with: original issue body + produced code + error output
3. Instruction: "Fix the error. The smoke test must pass."
4. If 3 failures → flag for manual review

---

## 2. Dependency Map & Build Order

### File Dependencies

```
requirements.txt           ← no deps (first file)
starmap/__init__.py         ← no deps (first file)
scripts/prepare_stars.py    ← requirements.txt (needs urllib)
data/stars.json             ← OUTPUT of: python scripts/prepare_stars.py
data/constellations.json    ← downloaded file, no code dep
starmap/catalog.py          ← data/stars.json, starmap/__init__.py
starmap/projection.py       ← starmap/__init__.py (uses astropy, numpy)
starmap/constellations.py   ← data/constellations.json, starmap/__init__.py
scripts/test_projection.py  ← starmap/catalog.py, starmap/projection.py
starmap/renderer.py         ← starmap/catalog.py, starmap/projection.py, starmap/constellations.py
scripts/test_render.py      ← starmap/renderer.py
app.py                      ← starmap/renderer.py
templates/index.html        ← app.py (Jinja2 context variables)
static/style.css            ← templates/index.html (CSS class names)
```

### Topological Build Order

```
Layer 0: requirements.txt, starmap/__init__.py
Layer 1: data/constellations.json (download), scripts/prepare_stars.py
Layer 2: data/stars.json (run prepare_stars.py)
Layer 3: starmap/catalog.py, starmap/projection.py, starmap/constellations.py  [parallel]
Layer 4: scripts/test_projection.py  [PHASE 1 GATE]
Layer 5: starmap/renderer.py
Layer 6: scripts/test_render.py  [PHASE 2 GATE]
Layer 7: app.py, templates/index.html, static/style.css  [semi-parallel]
Layer 8: Integration test  [PHASE 3 GATE]
Layer 9: Polish pass  [PHASE 4]
```

---

## 3. Phase 1: Data Pipeline + Projection Math

**Goal**: Produce the two data files, implement coordinate transforms and projection math.
**End state**: A test script proves stars project correctly for a known sky.

---

### Task 1.1 — Project Scaffolding *(GitHub Issue #1)*

**Files to create**: `requirements.txt`, `starmap/__init__.py`

**`requirements.txt`** — exact contents:
```
flask>=3.0
astropy>=6.0
numpy>=1.26
cairosvg>=2.7
```

**`starmap/__init__.py`** — exact contents:
```python
"""Star map poster generator."""
```

**Setup command** (main session runs this once after scaffolding):
```bash
python -m venv venv && source venv/bin/activate && pip install -r requirements.txt
```

**Smoke test**:
```bash
python -c "import starmap; print('OK')"
# Expected: OK
```

---

### Task 1.2 — Download Constellation Data *(GitHub Issue #2)*

**File to create**: `data/constellations.json`

**Commands**:
```bash
mkdir -p data
curl -L -o data/constellations.json \
  "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.lines.json"
```

**Smoke test**:
```bash
python -c "
import json
with open('data/constellations.json') as f:
    d = json.load(f)
print(f'Type: {d[\"type\"]}')
print(f'Features: {len(d[\"features\"])}')
print(f'First ID: {d[\"features\"][0][\"id\"]}')
assert d['type'] == 'FeatureCollection'
assert len(d['features']) == 88
print('OK')
"
```
Expected: `Type: FeatureCollection`, `Features: 88`, `OK`

---

### Task 1.3 — Prepare Star Catalog *(GitHub Issue #3)*

**File to create**: `scripts/prepare_stars.py`
**File produced**: `data/stars.json`
**Depends on**: Task 1.1 (requirements.txt)

**Complete specification**:

```python
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
            ra_deg = ra_hours * 15.0       # ← CRITICAL CONVERSION

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
```

**Run command**: `python scripts/prepare_stars.py`

**Smoke test**:
```bash
python -c "
import json
with open('data/stars.json') as f:
    stars = json.load(f)
print(f'Star count: {len(stars)}')
sirius = [s for s in stars if s.get('proper') == 'Sirius'][0]
print(f'Sirius: ra={sirius[\"ra\"]:.4f} dec={sirius[\"dec\"]:.4f} mag={sirius[\"mag\"]:.2f}')
assert 100 < sirius['ra'] < 102, f'Sirius RA wrong: {sirius[\"ra\"]}'
assert -17 < sirius['dec'] < -16, f'Sirius Dec wrong: {sirius[\"dec\"]}'
assert sirius['mag'] < -1, f'Sirius mag wrong: {sirius[\"mag\"]}'
assert 8000 <= len(stars) <= 10000, f'Star count unexpected: {len(stars)}'
print('ALL CHECKS PASSED')
"
```
Expected: `Star count: ~9000-9500`, `ALL CHECKS PASSED`

---

### Task 1.4 — Star Catalog Loader *(GitHub Issue #4)*

**File to create**: `starmap/catalog.py`
**Depends on**: Task 1.1, Task 1.3

```python
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
            _cache = json.load(f)
    return _cache
```

**Smoke test**:
```bash
python -c "
from starmap.catalog import load_stars
stars = load_stars()
print(f'Loaded {len(stars)} stars')
print(f'Brightest: {stars[0]}')
assert len(stars) > 8000
print('OK')
"
```

---

### Task 1.5 — Projection Module *(GitHub Issue #5)*

**File to create**: `starmap/projection.py`
**Depends on**: Task 1.1

```python
"""Coordinate transforms and stereographic projection.

Pipeline: ICRS (RA/Dec) → AltAz (altitude/azimuth) → Stereographic (x, y)
"""
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np
from numpy.typing import NDArray
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u


def build_altaz_frame(
    lat: float,          # decimal degrees, -90 to +90
    lon: float,          # decimal degrees, -180 to +180
    dt: datetime,        # naive datetime (no tzinfo)
    tz_name: str,        # IANA tz name, e.g. "America/New_York"
) -> AltAz:
    """Build AltAz reference frame for a specific time and location.

    Pseudo-logic:
        1. tz = ZoneInfo(tz_name)
        2. local_dt = dt.replace(tzinfo=tz)
        3. obs_time = astropy.time.Time(local_dt)
        4. location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
        5. return AltAz(obstime=obs_time, location=location)
    """
    tz = ZoneInfo(tz_name)
    local_dt = dt.replace(tzinfo=tz)
    obs_time = Time(local_dt)
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    return AltAz(obstime=obs_time, location=location)


def transform_radec_to_altaz(
    ra: NDArray[np.float64],    # RA array in degrees
    dec: NDArray[np.float64],   # Dec array in degrees
    frame: AltAz,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Bulk-transform RA/Dec arrays to altitude/azimuth arrays (degrees).

    Pseudo-logic:
        1. coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        2. altaz = coords.transform_to(frame)
        3. return (altaz.alt.deg, altaz.az.deg)

    NOTE: This is the slow step (~2-5 sec for 9000 stars). Call ONCE with
    all stars, not per-star.
    """
    coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    altaz = coords.transform_to(frame)
    return altaz.alt.deg, altaz.az.deg


def filter_visible(alt: NDArray[np.float64]) -> NDArray[np.bool_]:
    """Return boolean mask: True where altitude > 0 (above horizon)."""
    return alt > 0


def stereographic_project(
    alt: NDArray[np.float64],   # altitude in degrees (must be > 0)
    az: NDArray[np.float64],    # azimuth in degrees
    cx: float,                  # circle center X (SVG px)
    cy: float,                  # circle center Y (SVG px)
    radius: float,              # circle radius (SVG px)
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Stereographic projection from hemisphere to flat circle.

    Center = zenith (alt=90°), Edge = horizon (alt=0°).
    East appears on LEFT (correct for looking up at the sky).

    Math:
        alt_rad = radians(alt)
        az_rad  = radians(az)
        r_norm  = cos(alt_rad) / (1 + sin(alt_rad))   # 0 at zenith, 1 at horizon
        x = cx + r_norm * radius * sin(az_rad)
        y = cy - r_norm * radius * cos(az_rad)
    """
    alt_rad = np.radians(alt)
    az_rad = np.radians(az)
    r_norm = np.cos(alt_rad) / (1.0 + np.sin(alt_rad))
    x = cx + r_norm * radius * np.sin(az_rad)
    y = cy - r_norm * radius * np.cos(az_rad)
    return x, y
```

**Smoke test**:
```bash
python -c "
from starmap.projection import stereographic_project
import numpy as np

# Zenith → center
alt = np.array([90.0]); az = np.array([0.0])
x, y = stereographic_project(alt, az, 900, 850, 720)
print(f'Zenith: x={x[0]:.1f}, y={y[0]:.1f} (expect 900.0, 850.0)')
assert abs(x[0] - 900.0) < 0.01 and abs(y[0] - 850.0) < 0.01

# Horizon North → top edge
alt_h = np.array([0.001]); az_h = np.array([0.0])
x_h, y_h = stereographic_project(alt_h, az_h, 900, 850, 720)
print(f'Horizon-N: x={x_h[0]:.1f}, y={y_h[0]:.1f} (expect ~900, ~130)')
assert abs(x_h[0] - 900.0) < 1.0 and abs(y_h[0] - 130.0) < 5.0

print('Projection math OK')
"
```

---

### Task 1.6 — Constellation Loader *(GitHub Issue #6)*

**File to create**: `starmap/constellations.py`
**Depends on**: Task 1.1, Task 1.2

```python
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
```

**Smoke test**:
```bash
python -c "
from starmap.constellations import load_constellations, iter_line_segments
data = load_constellations()
print(f'Constellations: {len(data[\"features\"])}')
assert len(data['features']) == 88
segments = list(iter_line_segments())
print(f'Total line segments: {len(segments)}')
assert 500 < len(segments) < 1500
print(f'Sample: {segments[0]}')
print('OK')
"
```

---

### Task 1.7 — Phase 1 Gate Test *(GitHub Issue #7)*

**File to create**: `scripts/test_projection.py`
**Depends on**: Tasks 1.4, 1.5, 1.6

```python
#!/usr/bin/env python3
"""Phase 1 gate: transform stars for Boston 2023-03-15 21:00 EST, verify results."""
from datetime import datetime
import numpy as np
from starmap.catalog import load_stars
from starmap.projection import (
    build_altaz_frame, transform_radec_to_altaz,
    filter_visible, stereographic_project,
)


def main() -> None:
    lat, lon = 42.3601, -71.0589
    dt = datetime(2023, 3, 15, 21, 0)
    tz_name = "America/New_York"

    stars = load_stars()
    print(f"Loaded {len(stars)} stars")

    frame = build_altaz_frame(lat, lon, dt, tz_name)
    ra = np.array([s["ra"] for s in stars])
    dec = np.array([s["dec"] for s in stars])
    alt, az = transform_radec_to_altaz(ra, dec, frame)

    mask = filter_visible(alt)
    n_visible = mask.sum()
    print(f"Visible stars: {n_visible} / {len(stars)}")
    assert 3000 < n_visible < 6000, f"Unexpected: {n_visible}"

    x, y = stereographic_project(alt[mask], az[mask], 900, 850, 720)
    print(f"X range: [{x.min():.0f}, {x.max():.0f}]")
    print(f"Y range: [{y.min():.0f}, {y.max():.0f}]")

    # Check Polaris (HIP 11767)
    for i, s in enumerate(stars):
        if s.get("hip") == 11767:
            p_alt, p_az = alt[i], az[i]
            print(f"Polaris: alt={p_alt:.1f}°, az={p_az:.1f}°")
            assert 35 < p_alt < 50, f"Polaris alt wrong: {p_alt}"
            px, py = stereographic_project(
                np.array([p_alt]), np.array([p_az]), 900, 850, 720
            )
            print(f"Polaris projected: ({px[0]:.0f}, {py[0]:.0f})")
            assert 500 < px[0] < 1300 and 300 < py[0] < 900
            break

    print("\n=== PHASE 1 GATE: ALL TESTS PASSED ===")


if __name__ == "__main__":
    main()
```

**Gate command**:
```bash
python scripts/test_projection.py
# MUST print: === PHASE 1 GATE: ALL TESTS PASSED ===
```

---

## 4. Phase 2: SVG Renderer

**Goal**: Generate a complete, beautiful star map SVG.
**End state**: A test script writes an SVG that opens in a browser as a premium poster.

---

### Task 2.1 — SVG Renderer *(GitHub Issue #8)*

**File to create**: `starmap/renderer.py`
**Depends on**: Phase 1 gate passed (all starmap modules + data files)

This is the **largest and most critical** file. Complete specification:

```python
"""SVG star map renderer.

Generates print-quality star map poster as SVG string.
Canvas: 1800×2400 px (3:4 portrait).
"""
from datetime import datetime
from xml.sax.saxutils import escape
from zoneinfo import ZoneInfo

import numpy as np

from starmap.catalog import load_stars
from starmap.constellations import iter_line_segments
from starmap.projection import (
    build_altaz_frame, transform_radec_to_altaz,
    filter_visible, stereographic_project,
)

# ── Layout Constants ──────────────────────────────────────────────
CANVAS_W = 1800
CANVAS_H = 2400
CIRCLE_CX = 900
CIRCLE_CY = 850
CIRCLE_R = 720
BORDER_INSET = 40
BORDER_OPACITY = 0.6

# ── Text Constants ────────────────────────────────────────────────
TEXT_X = 900
TEXT_LINE1_Y = 1720       # custom message
TEXT_LINE2_Y = 1775       # DMS coordinates
TEXT_LINE3_Y = 1830       # date and time
TEXT_FONT = "'Georgia', 'Times New Roman', serif"
TEXT_SIZE_LINE1 = 28
TEXT_SIZE_LINE2 = 20
TEXT_SIZE_LINE3 = 20
TEXT_LETTER_SPACING = 4   # px, for small-caps effect
TEXT_COLOR = "#ffffff"
TEXT_OPACITY_SECONDARY = 0.7


def mag_to_radius(mag: float) -> float:
    """Star magnitude → SVG circle radius.
    Brighter (lower mag) = bigger. Range: ~0.3px to ~3.9px.
    """
    return max(0.3, 2.5 - mag * 0.35)


def mag_to_opacity(mag: float) -> float:
    """Star magnitude → opacity. Brighter = more opaque.
    Range: 0.5 (mag 6.5) to 1.0 (mag ≤ 2).
    """
    return max(0.5, min(1.0, 1.0 - (mag - 2.0) * 0.1))


def decimal_to_dms(value: float, is_lat: bool) -> str:
    """Decimal degrees → DMS string like '42° 20' 25" N'.

    Args:
        value: Decimal degrees (positive or negative).
        is_lat: True for N/S, False for E/W.
    """
    direction = ("N" if value >= 0 else "S") if is_lat else ("E" if value >= 0 else "W")
    value = abs(value)
    d = int(value)
    m = int((value - d) * 60)
    s = int(((value - d) * 60 - m) * 60)
    return f'{d}\u00b0 {m}\' {s}" {direction}'


def format_datetime(dt: datetime, tz_name: str) -> str:
    """Format datetime for poster: 'March 15, 2023, 09:00 PM'."""
    tz = ZoneInfo(tz_name)
    local_dt = dt.replace(tzinfo=tz)
    return local_dt.strftime("%B %d, %Y, %I:%M %p")


def render_star_map(
    lat: float,
    lon: float,
    dt: datetime,
    tz_name: str,
    message: str,
    show_constellations: bool = True,
) -> str:
    """Generate complete star map poster as SVG XML string.

    Full pipeline:
        1. Build AltAz frame
        2. Load + transform + filter + project stars
        3. Load + batch-transform + filter + project constellation lines
        4. Assemble SVG: background → border → clipPath → constellations → stars → circle → text

    Returns complete SVG XML string.
    """
    # ── 1. Build frame ──
    frame = build_altaz_frame(lat, lon, dt, tz_name)

    # ── 2. Stars ──
    stars = load_stars()
    ra = np.array([s["ra"] for s in stars])
    dec = np.array([s["dec"] for s in stars])
    mags = np.array([s["mag"] for s in stars])

    alt, az = transform_radec_to_altaz(ra, dec, frame)
    mask = filter_visible(alt)
    vis_alt, vis_az, vis_mag = alt[mask], az[mask], mags[mask]
    sx, sy = stereographic_project(vis_alt, vis_az, CIRCLE_CX, CIRCLE_CY, CIRCLE_R)

    # ── 3. Constellation lines (BATCHED — transform once) ──
    const_lines_svg = ""
    if show_constellations:
        segments = list(iter_line_segments())
        if segments:
            all_ra1 = np.array([s[0] for s in segments])
            all_dec1 = np.array([s[1] for s in segments])
            all_ra2 = np.array([s[2] for s in segments])
            all_dec2 = np.array([s[3] for s in segments])

            # Concatenate all endpoints, transform ONCE
            cat_ra = np.concatenate([all_ra1, all_ra2])
            cat_dec = np.concatenate([all_dec1, all_dec2])
            cat_alt, cat_az = transform_radec_to_altaz(cat_ra, cat_dec, frame)

            n = len(segments)
            a1, a2 = cat_alt[:n], cat_alt[n:]
            z1, z2 = cat_az[:n], cat_az[n:]

            both_vis = (a1 > 0) & (a2 > 0)
            if both_vis.any():
                lx1, ly1 = stereographic_project(a1[both_vis], z1[both_vis], CIRCLE_CX, CIRCLE_CY, CIRCLE_R)
                lx2, ly2 = stereographic_project(a2[both_vis], z2[both_vis], CIRCLE_CX, CIRCLE_CY, CIRCLE_R)
                line_parts = []
                for j in range(len(lx1)):
                    line_parts.append(
                        f'<line x1="{lx1[j]:.1f}" y1="{ly1[j]:.1f}" '
                        f'x2="{lx2[j]:.1f}" y2="{ly2[j]:.1f}" '
                        f'stroke="#ffffff" stroke-width="0.5" opacity="0.3"/>'
                    )
                const_lines_svg = "\n".join(line_parts)

    # ── 4. Assemble SVG ──
    parts: list[str] = []

    # Header
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {CANVAS_W} {CANVAS_H}" '
        f'width="{CANVAS_W}" height="{CANVAS_H}">'
    )

    # Black background
    parts.append('<rect width="100%" height="100%" fill="#000000"/>')

    # Outer border
    bx, by = BORDER_INSET, BORDER_INSET
    bw, bh = CANVAS_W - 2 * BORDER_INSET, CANVAS_H - 2 * BORDER_INSET
    parts.append(
        f'<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" '
        f'fill="none" stroke="#ffffff" stroke-width="1" opacity="{BORDER_OPACITY}"/>'
    )

    # Clip path
    parts.append("<defs>")
    parts.append(f'<clipPath id="sky-clip">')
    parts.append(f'<circle cx="{CIRCLE_CX}" cy="{CIRCLE_CY}" r="{CIRCLE_R}"/>')
    parts.append("</clipPath>")
    parts.append("</defs>")

    # Constellation lines (behind stars, clipped)
    if const_lines_svg:
        parts.append('<g clip-path="url(#sky-clip)">')
        parts.append(const_lines_svg)
        parts.append("</g>")

    # Stars (on top of constellations, clipped)
    parts.append('<g clip-path="url(#sky-clip)">')
    for i in range(len(sx)):
        r = mag_to_radius(vis_mag[i])
        op = mag_to_opacity(vis_mag[i])
        parts.append(
            f'<circle cx="{sx[i]:.1f}" cy="{sy[i]:.1f}" r="{r:.2f}" '
            f'fill="#ffffff" opacity="{op:.2f}"/>'
        )
    parts.append("</g>")

    # Circle border
    parts.append(
        f'<circle cx="{CIRCLE_CX}" cy="{CIRCLE_CY}" r="{CIRCLE_R}" '
        f'fill="none" stroke="#ffffff" stroke-width="1.5"/>'
    )

    # Text block
    msg_escaped = escape(message.upper())
    coords_text = f"{decimal_to_dms(lat, True)}  {decimal_to_dms(lon, False)}"
    date_text = format_datetime(dt, tz_name).upper()

    parts.append(
        f'<text x="{TEXT_X}" y="{TEXT_LINE1_Y}" text-anchor="middle" '
        f'fill="{TEXT_COLOR}" font-family="{TEXT_FONT}" '
        f'font-size="{TEXT_SIZE_LINE1}" letter-spacing="{TEXT_LETTER_SPACING}">'
        f'{msg_escaped}</text>'
    )
    parts.append(
        f'<text x="{TEXT_X}" y="{TEXT_LINE2_Y}" text-anchor="middle" '
        f'fill="{TEXT_COLOR}" font-family="{TEXT_FONT}" '
        f'font-size="{TEXT_SIZE_LINE2}" letter-spacing="{TEXT_LETTER_SPACING}" '
        f'opacity="{TEXT_OPACITY_SECONDARY}">'
        f'{coords_text}</text>'
    )
    parts.append(
        f'<text x="{TEXT_X}" y="{TEXT_LINE3_Y}" text-anchor="middle" '
        f'fill="{TEXT_COLOR}" font-family="{TEXT_FONT}" '
        f'font-size="{TEXT_SIZE_LINE3}" letter-spacing="{TEXT_LETTER_SPACING}" '
        f'opacity="{TEXT_OPACITY_SECONDARY}">'
        f'{date_text}</text>'
    )

    parts.append("</svg>")
    return "\n".join(parts)
```

**Performance note**: The batched constellation transform is critical. Without it, you'd have ~700 individual `transform_radec_to_altaz` calls (minutes). With batching: 2 calls total (seconds).

---

### Task 2.2 — Phase 2 Gate Test *(GitHub Issue #9)*

**File to create**: `scripts/test_render.py`
**Depends on**: Task 2.1

```python
#!/usr/bin/env python3
"""Phase 2 gate: render SVG, verify structure, write test_output.svg."""
from datetime import datetime
from starmap.renderer import render_star_map, decimal_to_dms, mag_to_radius


def main() -> None:
    svg = render_star_map(
        lat=42.3601, lon=-71.0589,
        dt=datetime(2023, 3, 15, 21, 0),
        tz_name="America/New_York",
        message="The Night We Met",
        show_constellations=True,
    )

    with open("test_output.svg", "w") as f:
        f.write(svg)
    print(f"Wrote test_output.svg ({len(svg)} bytes)")

    # Structural checks
    assert svg.startswith("<svg"), "Missing <svg"
    assert "</svg>" in svg, "Missing </svg>"
    assert 'viewBox="0 0 1800 2400"' in svg, "Wrong viewBox"
    assert "sky-clip" in svg, "Missing clip path"
    assert "THE NIGHT WE MET" in svg, "Message not uppercase"

    star_count = svg.count('<circle cx="') - 2  # minus clipPath + border circles
    print(f"Star circles: {star_count}")
    assert star_count > 3000, f"Too few stars: {star_count}"

    line_count = svg.count("<line ")
    print(f"Constellation lines: {line_count}")
    assert line_count > 100, f"Too few lines: {line_count}"

    assert "42" in svg and "20" in svg, "DMS lat missing"
    assert "71" in svg, "DMS lon missing"

    # Helper tests
    assert mag_to_radius(-1.46) > 3.0
    assert mag_to_radius(6.0) < 0.5

    print("\n=== PHASE 2 GATE: ALL TESTS PASSED ===")
    print("Open test_output.svg in a browser to visually verify.")


if __name__ == "__main__":
    main()
```

**Gate command**:
```bash
python scripts/test_render.py
# MUST print: === PHASE 2 GATE: ALL TESTS PASSED ===
# Then visually verify test_output.svg in browser
```

**Visual verification checklist** (main session confirms by opening the SVG):
- [ ] Black background with white rectangular border
- [ ] Circular star field with varying-size white dots
- [ ] Constellation lines visible as thin subtle connections
- [ ] Text below circle: message, coordinates, date
- [ ] Recognizable constellations (Orion, Big Dipper visible for Boston winter)

---

## 5. Phase 3: Flask Application + Web UI

**Goal**: Serve a web form that generates and displays star maps with download options.
**End state**: `python app.py` → browser → fill form → see poster → download SVG/PNG.

---

### Task 3.1 — Flask Application *(GitHub Issue #10)*

**File to create**: `app.py`
**Depends on**: Phase 2 gate passed

```python
"""Flask web application for star map poster generation.

Routes:
    GET  /              — Form page
    POST /generate      — Generate star map, return page with inline SVG preview
    POST /download/svg  — SVG file download
    POST /download/png  — PNG via CairoSVG (graceful fallback)
"""
from datetime import datetime
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

from flask import Flask, render_template, request, Response
from starmap.renderer import render_star_map

app = Flask(__name__)

TIMEZONES = [
    "America/New_York", "America/Chicago", "America/Denver",
    "America/Los_Angeles", "America/Anchorage", "Pacific/Honolulu",
    "Europe/London", "Europe/Paris", "Europe/Berlin",
    "Asia/Tokyo", "Asia/Shanghai", "Australia/Sydney", "UTC",
]


def _parse_form() -> tuple[float, float, datetime, str, str, bool] | str:
    """Parse and validate form data. Returns tuple on success, error string on failure.

    Validates:
        - lat: float, -90 to 90
        - lon: float, -180 to 180
        - date + time: parseable as datetime
        - tz: valid IANA timezone
        - message: defaults to location_name or "The Night Sky"
        - constellations: checkbox boolean
    """
    try:
        lat = float(request.form.get("lat", ""))
        if not -90 <= lat <= 90:
            return "Latitude must be between -90 and 90."
    except ValueError:
        return "Invalid latitude. Enter a number like 42.36."

    try:
        lon = float(request.form.get("lon", ""))
        if not -180 <= lon <= 180:
            return "Longitude must be between -180 and 180."
    except ValueError:
        return "Invalid longitude. Enter a number like -71.06."

    date_str = request.form.get("date", "")
    time_str = request.form.get("time", "")
    if not date_str or not time_str:
        return "Date and time are required."
    try:
        dt = datetime.strptime(f"{date_str} {time_str}", "%Y-%m-%d %H:%M")
    except ValueError:
        return "Invalid date or time format."

    tz_name = request.form.get("tz", "UTC")
    try:
        ZoneInfo(tz_name)
    except (ZoneInfoNotFoundError, KeyError):
        return f"Unknown timezone: {tz_name}"

    message = request.form.get("message", "").strip()
    if not message:
        message = request.form.get("location_name", "").strip()
    if not message:
        message = "The Night Sky"

    show_constellations = request.form.get("constellations") == "on"
    return (lat, lon, dt, tz_name, message, show_constellations)


@app.route("/")
def index() -> str:
    return render_template("index.html", timezones=TIMEZONES, svg=None, error=None, form_data={})


@app.route("/generate", methods=["POST"])
def generate() -> str:
    form_data = request.form.to_dict()
    result = _parse_form()
    if isinstance(result, str):
        return render_template("index.html", timezones=TIMEZONES, svg=None,
                               error=result, form_data=form_data)
    lat, lon, dt, tz_name, message, show_constellations = result
    svg = render_star_map(lat, lon, dt, tz_name, message, show_constellations)
    return render_template("index.html", timezones=TIMEZONES, svg=svg,
                           error=None, form_data=form_data)


@app.route("/download/svg", methods=["POST"])
def download_svg() -> Response:
    result = _parse_form()
    if isinstance(result, str):
        return Response(result, status=400, mimetype="text/plain")
    lat, lon, dt, tz_name, message, show_constellations = result
    svg = render_star_map(lat, lon, dt, tz_name, message, show_constellations)
    return Response(svg, mimetype="image/svg+xml",
                    headers={"Content-Disposition": "attachment; filename=starmap.svg"})


@app.route("/download/png", methods=["POST"])
def download_png() -> Response:
    result = _parse_form()
    if isinstance(result, str):
        return Response(result, status=400, mimetype="text/plain")
    lat, lon, dt, tz_name, message, show_constellations = result
    svg = render_star_map(lat, lon, dt, tz_name, message, show_constellations)
    try:
        import cairosvg
        png_bytes = cairosvg.svg2png(bytestring=svg.encode("utf-8"),
                                      output_width=3600, output_height=4800)
    except ImportError:
        return Response("PNG export requires cairosvg. Install: pip install cairosvg",
                        status=500, mimetype="text/plain")
    return Response(png_bytes, mimetype="image/png",
                    headers={"Content-Disposition": "attachment; filename=starmap.png"})


if __name__ == "__main__":
    app.run(debug=True, port=5000)
```

---

### Task 3.2 — HTML Template *(GitHub Issue #11)*

**File to create**: `templates/index.html`
**Depends on**: Task 3.1 (needs Jinja2 context variables)

**Template context variables** (from `app.py`):
- `timezones` — `list[str]` of timezone names
- `svg` — SVG XML string or `None`
- `error` — error message string or `None`
- `form_data` — `dict` of submitted form values (for retaining values after POST)

**Requirements**:
- Main form: `method="POST"` `action="/generate"`
- Fields (exact `name` attributes matching `_parse_form()`):
  - `location_name` — text input, placeholder "910 @ 400 Fenway"
  - `lat` — number input, step=0.0001, default 42.3601
  - `lon` — number input, step=0.0001, default -71.0589
  - `date` — date input, default today
  - `time` — time input, default "21:00"
  - `tz` — `<select>` from `timezones`, default "America/New_York"
  - `message` — text input, optional, placeholder "The night we met..."
  - `constellations` — checkbox, checked by default
- Submit button: "Generate Star Map" → `btn btn-primary`
- Error display: if `error`, show in `.error-box` div
- SVG preview: if `svg`, render `{{ svg | safe }}` in `.preview-section` div
- Download buttons (only shown when `svg` is set):
  - SVG download: a `<form>` with `action="/download/svg"` + hidden fields duplicating all inputs
  - PNG download: a `<form>` with `action="/download/png"` + hidden fields duplicating all inputs
  - Both use `btn btn-download`
- All values pre-populated from `form_data` dict using `value="{{ form_data.get('lat', '42.3601') }}"`
- Layout classes: `.container`, `.main-layout`, `.form-section`, `.preview-section`
- `<meta name="viewport" content="width=device-width, initial-scale=1.0">`
- `<link rel="stylesheet" href="{{ url_for('static', filename='style.css') }}">`
- Title: "Star Map Generator"

---

### Task 3.3 — CSS Stylesheet *(GitHub Issue #12)*

**File to create**: `static/style.css`
**Can be built in parallel with Task 3.2** (CSS class names are pre-defined above)

**Complete specification**:

```css
/* === Reset & Base === */
*, *::before, *::after {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}

body {
    background: #0a0a0a;
    color: #e0e0e0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
    line-height: 1.6;
    min-height: 100vh;
}

/* === Layout === */
.container {
    max-width: 1200px;
    margin: 0 auto;
    padding: 2rem 1.5rem;
}

h1 {
    font-family: 'Georgia', 'Times New Roman', serif;
    font-size: 1.8rem;
    font-weight: 400;
    letter-spacing: 3px;
    text-transform: uppercase;
    color: #ffffff;
    margin-bottom: 2rem;
    text-align: center;
}

/* === Form === */
.form-section {
    max-width: 500px;
    margin: 0 auto 3rem;
}

.form-group {
    margin-bottom: 1.25rem;
}

label {
    display: block;
    font-size: 0.8rem;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    color: #888;
    margin-bottom: 0.4rem;
}

input[type="text"],
input[type="number"],
input[type="date"],
input[type="time"],
select {
    width: 100%;
    padding: 0.6rem 0.8rem;
    background: #1a1a1a;
    border: 1px solid #333;
    border-radius: 4px;
    color: #ffffff;
    font-size: 0.95rem;
    transition: border-color 0.2s;
}

input:focus,
select:focus {
    outline: none;
    border-color: #555;
}

/* Checkbox styling */
.checkbox-group {
    display: flex;
    align-items: center;
    gap: 0.5rem;
}

.checkbox-group input[type="checkbox"] {
    width: 16px;
    height: 16px;
    accent-color: #666;
}

.checkbox-group label {
    margin-bottom: 0;
    font-size: 0.85rem;
}

/* Lat/lon row */
.form-row {
    display: flex;
    gap: 1rem;
}

.form-row .form-group {
    flex: 1;
}

/* === Buttons === */
.btn {
    display: inline-block;
    padding: 0.7rem 1.8rem;
    border: none;
    border-radius: 4px;
    font-size: 0.85rem;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    cursor: pointer;
    transition: background-color 0.2s, border-color 0.2s;
}

.btn-primary {
    width: 100%;
    background: #ffffff;
    color: #000000;
    margin-top: 0.5rem;
}

.btn-primary:hover {
    background: #e0e0e0;
}

.btn-download {
    background: transparent;
    border: 1px solid #444;
    color: #ccc;
}

.btn-download:hover {
    border-color: #888;
    color: #fff;
}

/* === Error === */
.error-box {
    background: #2a1515;
    border: 1px solid #5a2020;
    border-radius: 4px;
    padding: 0.8rem 1rem;
    color: #ff8888;
    font-size: 0.9rem;
    margin-bottom: 1.5rem;
}

/* === Preview === */
.preview-section {
    text-align: center;
}

.preview-section svg {
    max-width: 100%;
    height: auto;
    border: 1px solid #222;
}

.download-buttons {
    margin-top: 1.5rem;
    display: flex;
    justify-content: center;
    gap: 0.75rem;
}

/* === Responsive === */
@media (min-width: 900px) {
    .main-layout {
        display: flex;
        gap: 3rem;
        align-items: flex-start;
    }

    .form-section {
        flex: 0 0 400px;
        margin-bottom: 0;
    }

    .preview-section {
        flex: 1;
    }
}
```

---

### Task 3.4 — Phase 3 Gate Test *(GitHub Issue #13)*

No separate test script. Run manually:

```bash
# Start the app in background
python app.py &
APP_PID=$!
sleep 3

# Test GET /
curl -s http://localhost:5000/ | grep -c "<form"
# Expected: 1

# Test POST /generate
curl -s -X POST http://localhost:5000/generate \
  -d "lat=42.36&lon=-71.06&date=2023-03-15&time=21:00&tz=America/New_York&message=Test&constellations=on" \
  | grep -c "<svg"
# Expected: 1

# Test SVG download
curl -s -X POST http://localhost:5000/download/svg \
  -d "lat=42.36&lon=-71.06&date=2023-03-15&time=21:00&tz=America/New_York&message=Test&constellations=on" \
  -o /tmp/test_dl.svg
ls -la /tmp/test_dl.svg
# Should be > 100KB

# Test validation error
curl -s -X POST http://localhost:5000/generate \
  -d "lat=999&lon=-71.06&date=2023-03-15&time=21:00&tz=UTC&message=Test" \
  | grep -c "error"
# Expected: >= 1

kill $APP_PID
echo "=== PHASE 3 GATE: ALL TESTS PASSED ==="
```

**Visual verification** (main session opens browser to `localhost:5000`):
- [ ] Dark-themed form renders
- [ ] Form submits and shows SVG preview inline
- [ ] SVG download produces valid file
- [ ] PNG download works (if cairosvg installed)
- [ ] Invalid inputs show error, form values retained
- [ ] Responsive layout works (form + preview side by side on desktop)

---

## 6. Phase 4: Polish + Edge Cases

**Goal**: Visual quality, edge cases, type safety, PNG verification.

---

### Task 4.1 — Visual Tuning *(GitHub Issue #14)*

**Modify**: `starmap/renderer.py`
**Depends on**: Phase 3 gate

Changes:
- **Bright star glow**: For stars with `mag < 1.0`, render a second `<circle>` behind the main dot: 3× radius, opacity 0.08, creates subtle halo
- **Magnitude curve tuning**: Adjust if visual output is off
- **Long message handling**: If message > 40 chars, reduce `TEXT_SIZE_LINE1` to 22

```python
# Glow implementation (in the star rendering loop):
if vis_mag[i] < 1.0:
    glow_r = r * 3.0
    parts.append(
        f'<circle cx="{sx[i]:.1f}" cy="{sy[i]:.1f}" r="{glow_r:.2f}" '
        f'fill="#ffffff" opacity="0.08"/>'
    )
```

---

### Task 4.2 — Edge Case Hardening *(GitHub Issue #15)*

**Modify**: `starmap/renderer.py`, `app.py`
**Depends on**: Phase 3 gate

Test and fix:
- North Pole observer (`lat=90`): Polaris at zenith, all circumpolar stars visible
- South Pole observer (`lat=-90`): reversed sky
- Empty message → default "The Night Sky"
- Unicode/special chars in message → `escape()` handles it
- Midnight display: "12:00 AM" not "00:00 AM"
- Noon: "12:00 PM"

**Test commands**:
```bash
python -c "
from datetime import datetime
from starmap.renderer import render_star_map

# North pole
svg = render_star_map(90.0, 0.0, datetime(2023, 6, 21, 0, 0), 'UTC', 'North Pole')
assert '<svg' in svg
print(f'North Pole SVG: {len(svg)} bytes')

# South pole
svg = render_star_map(-90.0, 0.0, datetime(2023, 6, 21, 0, 0), 'UTC', 'South Pole')
assert '<svg' in svg
print(f'South Pole SVG: {len(svg)} bytes')

# Special characters
svg = render_star_map(42.0, -71.0, datetime(2023, 3, 15, 21, 0), 'UTC', 'Café <3 & Stars!')
assert '&amp;' in svg  # escaped ampersand
assert '&lt;' in svg    # escaped <
print('Special chars OK')

print('Edge cases OK')
"
```

---

### Task 4.3 — Code Quality *(GitHub Issue #16)*

**Modify**: All `starmap/*.py`, `app.py`
**Depends on**: Tasks 4.1, 4.2

- Verify all public functions have type hints + docstrings (most already do from specs)
- Remove any hardcoded absolute paths
- Add to `.gitignore`:
  ```
  test_output.svg
  test_output.png
  scripts/hygdata_v41.csv
  ```

---

### Task 4.4 — PNG Export Verification *(GitHub Issue #17)*

**Depends on**: Phase 3 gate

```bash
python -c "
import cairosvg
from starmap.renderer import render_star_map
from datetime import datetime

svg = render_star_map(42.36, -71.06, datetime(2023, 3, 15, 21, 0),
                      'America/New_York', 'Test', True)
png = cairosvg.svg2png(bytestring=svg.encode(), output_width=3600, output_height=4800)
with open('test_output.png', 'wb') as f:
    f.write(png)
print(f'PNG size: {len(png)} bytes')
assert len(png) > 100000, 'PNG too small'
print('PNG export OK')
"
```

---

### Phase 4 Final Gate

```bash
python app.py &
APP_PID=$!
sleep 3

# 3 different locations
for params in \
  "lat=42.36&lon=-71.06&date=2023-12-15&time=21:00&tz=America/New_York&message=Boston+Winter&constellations=on" \
  "lat=-33.87&lon=151.21&date=2023-07-15&time=21:00&tz=Australia/Sydney&message=Sydney+Winter&constellations=on" \
  "lat=0.0&lon=0.0&date=2023-06-21&time=21:00&tz=UTC&message=Equator+Solstice&constellations=on"; do
    curl -s -X POST http://localhost:5000/download/svg \
      -d "$params" -o /dev/null -w "SVG: %{http_code} %{size_download}\n"
done

# PNG test
curl -s -X POST http://localhost:5000/download/png \
  -d "lat=42.36&lon=-71.06&date=2023-03-15&time=21:00&tz=America/New_York&message=Test&constellations=on" \
  -o /dev/null -w "PNG: %{http_code} %{size_download}\n"

kill $APP_PID
echo "=== PHASE 4 FINAL GATE: ALL TESTS PASSED ==="
```

---

## 7. GitHub Issues Reference

17 issues across 4 milestones. Create in order. Dependencies listed per issue.

| # | Title | Labels | Milestone | Blocked By |
|---|-------|--------|-----------|------------|
| 1 | Create requirements.txt and starmap/__init__.py | `phase-1`, `scaffolding` | Phase 1 | — |
| 2 | Download d3-celestial constellation lines | `phase-1`, `data` | Phase 1 | — |
| 3 | Create star catalog prep script and generate data/stars.json | `phase-1`, `data` | Phase 1 | #1 |
| 4 | Implement starmap/catalog.py | `phase-1`, `core` | Phase 1 | #1, #3 |
| 5 | Implement starmap/projection.py | `phase-1`, `core` | Phase 1 | #1 |
| 6 | Implement starmap/constellations.py | `phase-1`, `core` | Phase 1 | #1, #2 |
| 7 | Phase 1 gate test — scripts/test_projection.py | `phase-1`, `test`, `gate` | Phase 1 | #4, #5, #6 |
| 8 | Implement starmap/renderer.py | `phase-2`, `core` | Phase 2 | #7 |
| 9 | Phase 2 gate test — scripts/test_render.py | `phase-2`, `test`, `gate` | Phase 2 | #8 |
| 10 | Implement app.py Flask application | `phase-3`, `web` | Phase 3 | #9 |
| 11 | Create templates/index.html | `phase-3`, `web`, `ui` | Phase 3 | #10 |
| 12 | Create static/style.css dark theme | `phase-3`, `web`, `ui` | Phase 3 | #10 |
| 13 | Phase 3 gate test — integration verification | `phase-3`, `test`, `gate` | Phase 3 | #11, #12 |
| 14 | Visual tuning — bright star glow, magnitude curves | `phase-4`, `polish` | Phase 4 | #13 |
| 15 | Edge case hardening — poles, special chars, empty input | `phase-4`, `polish` | Phase 4 | #13 |
| 16 | Code quality — type hints, docstrings, .gitignore cleanup | `phase-4`, `polish` | Phase 4 | #14, #15 |
| 17 | PNG export verification at print quality | `phase-4`, `test` | Phase 4 | #13 |

**Issue body template** — each issue body should contain:
1. Task description (from this document)
2. Complete file spec with function signatures
3. Smoke test commands with expected output
4. Acceptance criteria checklist
5. "Blocked by: #X, #Y" line

---

## 8. Data Format Reference

### data/stars.json

```json
[
  {"ra": 101.2872, "dec": -16.7161, "mag": -1.46, "hip": 32349, "proper": "Sirius", "con": "CMa"},
  {"ra": 213.9153, "dec": 19.1822, "mag": -0.05, "hip": 69673, "proper": "Arcturus", "con": "Boo"}
]
```

| Field | Type | Range | Notes |
|-------|------|-------|-------|
| `ra` | float | 0–360 | Degrees (converted from HYG hours × 15) |
| `dec` | float | -90 to +90 | Degrees |
| `mag` | float | -1.46 to 6.5 | Apparent visual magnitude |
| `hip` | int \| null | — | Hipparcos catalog ID |
| `proper` | str \| null | — | Common name (e.g., "Sirius") |
| `con` | str \| null | — | 3-letter IAU constellation abbreviation |

### data/constellations.json

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "id": "And",
      "properties": {"rank": "1"},
      "geometry": {
        "type": "MultiLineString",
        "coordinates": [
          [[30.97, 42.33], [17.43, 35.62], [9.83, 30.86]],
          [[14.30, 23.42], [11.83, 24.27]]
        ]
      }
    }
  ]
}
```

| Field | Notes |
|-------|-------|
| Coordinates | `[ra_degrees, dec_degrees]` per point |
| RA range | -180 to +180 (negative values valid, Astropy handles) |
| Feature count | 88 constellations |
| Segment count | ~700-800 total line segments |
