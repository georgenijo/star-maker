# CLAUDE.md -- star-maker

Star map poster generator. Takes a date, time, and location; renders an astronomically accurate SVG star chart suitable for printing.

## Tech Stack

- Python 3.11+
- Flask >= 3.0 -- web UI and routes
- Astropy >= 6.0 -- ICRS-to-AltAz coordinate transforms, time handling
- NumPy >= 1.26 -- vectorized math
- CairoSVG >= 2.7 -- optional PNG export
- `zoneinfo` (stdlib) -- timezone handling (not pytz)

## Project Structure

```
star-maker/
├── CLAUDE.md              # This file -- project context
├── PLAN.md                # Phased build roadmap
├── requirements.txt
├── app.py                 # Flask app -- routes, form, SVG/PNG endpoints
├── starmap/
│   ├── __init__.py
│   ├── catalog.py         # Load + filter stars.json
│   ├── projection.py      # ICRS -> AltAz -> stereographic (x,y)
│   ├── constellations.py  # Load + transform constellation lines
│   └── renderer.py        # SVG XML string generation
├── data/
│   ├── stars.json         # Filtered HYG catalog (~9k stars, mag <= 6.5)
│   └── constellations.json # d3-celestial constellation lines (RA/Dec GeoJSON)
├── static/
│   └── style.css          # Dark theme form styling
└── templates/
    └── index.html         # Flask template -- form + inline SVG preview
```

## Architecture Decisions

### Star Catalog: HYG Database v4.2
- Source CSV has ~120k stars. We preprocess offline to ~9k stars (mag <= 6.5) and bundle as JSON.
- Fields needed: `ra` (deg), `dec` (deg), `mag`, `hip` (Hipparcos ID), `proper` (name), `con` (constellation abbr).
- Why HYG: free (CC BY-SA 4.0), includes HIP + HR cross-references, comprehensive to naked-eye limit.
- Source: https://www.astronexus.com/hyg

### Constellation Lines: d3-celestial
- File: `constellations.lines.json` from github.com/ofrohn/d3-celestial (MIT license).
- Format: GeoJSON FeatureCollection. Each feature is a constellation with `MultiLineString` geometry.
- Coordinates are RA/Dec in degrees. These go through the same ICRS -> AltAz -> stereographic pipeline as catalog stars. No star-ID lookup needed.

### Projection: Stereographic
- Formula: `r = cos(alt) / (1 + sin(alt))`, normalized to circle radius.
- Center = zenith, edge = horizon.
- `x = cx + r * R * sin(az)`, `y = cy - r * R * cos(az)`.
- East appears on the LEFT (correct for looking up at the sky).

### Coordinate Pipeline
- Astropy `SkyCoord` (ICRS) bulk-transformed to `AltAz` for given time/location.
- Vectorized: pass all star coords as arrays, not one-by-one.
- Time via `astropy.time.Time` + `zoneinfo`. Location via `EarthLocation`.

### SVG Rendering
- Programmatic XML string building (no matplotlib, no template engine for SVG).
- Canvas: 1800 x 2400 px (3:4 portrait). Star circle: centered at (900, 850), radius 720px.
- Star sizing: `radius = max(0.3, 2.5 - mag * 0.35)`.
- Font: Georgia with serif fallback. Uppercase + `letter-spacing` for small-caps effect.
- Constellation lines: white, opacity 0.3, stroke-width 0.5.

### PNG Export
- CairoSVG converts SVG string to PNG. For print quality, scale up (e.g., 7200x10800 for 24x36" at 300dpi).

## Commands

```bash
# Setup
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Run
python app.py
# Serves at http://localhost:5000

# PNG export requires system Cairo library
# macOS: brew install cairo
# Ubuntu: sudo apt install libcairo2-dev
```

## Gotchas

1. **CairoSVG does not support `feGaussianBlur`**. SVG filter effects are silently ignored in PNG output. Primary reason Milky Way is deferred to v2.
2. **CairoSVG font rendering** requires system-installed fonts. Georgia must be present or fallback serif will be used.
3. **SVG `font-variant: small-caps`** is unreliable across renderers. Use uppercase + `letter-spacing` instead.
4. **Stereographic distortion** grows near the horizon. Acceptable for decorative poster use.
5. **d3-celestial RA convention**: RA stored as longitude (negative values exist). Astropy handles negative RA in degrees fine.
6. **Astropy first-run**: may download IERS data on first coordinate transform. Needs internet on first run.
7. **Coordinate display**: convert decimal lat/lon to DMS for poster text (e.g., `42deg 20' 25" N`).

## Scope

### v1
- Star field (magnitude-scaled dots with opacity variation)
- Constellation lines (toggleable)
- Stereographic projection with correct sky orientation
- Text block: custom message, DMS coordinates, formatted date/time
- Circular star field + rectangular outer border
- Flask form with all inputs
- SVG and PNG download

### Deferred (v2+)
- Milky Way band rendering
- Heart-shaped clip mask
- Color themes
- Moon/planet positions
- Grid overlays
- Geocoding (city name -> lat/lon)
- Font selection
