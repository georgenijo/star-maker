# PLAN.md -- star-maker build roadmap

Phased build plan. Each phase produces a testable artifact. Phases are sequential.

---

## Phase 1: Data Pipeline

**Goal**: Download, process, and validate the two data files. Implement coordinate transform and projection math.

### Tasks

1. **Download and preprocess HYG v4.2**
   - Download CSV from https://www.astronexus.com/hyg
   - Write a one-off script (`scripts/prepare_stars.py`) that:
     - Reads the CSV, filters to `mag <= 6.5`
     - Extracts: `ra`, `dec`, `mag`, `hip`, `proper`, `con`
     - Writes `data/stars.json`
   - Commit the output JSON (not the raw CSV)

2. **Obtain constellation lines**
   - Download `constellations.lines.json` from d3-celestial GitHub
   - Save as `data/constellations.json`

3. **Implement `starmap/catalog.py`**
   - `load_stars()` -> list of star dicts from `data/stars.json`
   - Cache on first call

4. **Implement `starmap/constellations.py`**
   - `load_constellations()` -> parsed GeoJSON features
   - Helper to iterate line segments as `(ra1, dec1, ra2, dec2)` pairs

5. **Implement `starmap/projection.py`**
   - `build_altaz_frame(lat, lon, dt, tz_name)` -> AltAz frame
   - `transform_radec_to_altaz(ra_arr, dec_arr, altaz_frame)` -> (alt_arr, az_arr)
   - `stereographic_project(alt_arr, az_arr, cx, cy, radius)` -> (x_arr, y_arr)
   - `filter_visible(alt_arr)` -> boolean mask (alt > 0)

6. **Validate**: test script (`scripts/test_projection.py`) that transforms stars for a known location/time and spot-checks Polaris position

### Files
- `scripts/prepare_stars.py`, `scripts/test_projection.py`
- `data/stars.json`, `data/constellations.json`
- `starmap/__init__.py`, `starmap/catalog.py`, `starmap/constellations.py`, `starmap/projection.py`
- `requirements.txt`

### Done When
- `data/stars.json` has ~8,000-9,500 entries
- `data/constellations.json` has 88 constellation features
- Test script runs, prints plausible visible star count, Polaris near center for northern hemisphere

---

## Phase 2: SVG Renderer

**Goal**: Generate a complete star map SVG from inputs. Output a file that opens in a browser.

### Tasks

1. **Implement `starmap/renderer.py`**
   - `render_star_map(params) -> str` (SVG XML string)
   - Params: `lat`, `lon`, `date`, `time`, `tz_name`, `message`, `show_constellations`
   - Pipeline: load data -> transform coords -> project -> build SVG

2. **SVG structure**:
   - Black background, outer border (40px inset)
   - `<clipPath>` circle for star field
   - Star dots (clipped): `<circle>` per star, sized by magnitude
   - Constellation lines (clipped): `<line>` per segment
   - Circle border stroke
   - Text block: message, DMS coords, formatted date

3. **Helpers**: `mag_to_radius()`, `mag_to_opacity()`, `decimal_to_dms()`, `format_datetime()`

4. **Validate**: `scripts/test_render.py` writes `test_output.svg`, visually verify in browser

### Files
- `starmap/renderer.py`
- `scripts/test_render.py`

### Done When
- SVG opens in browser showing black poster with circular star field
- Stars visible as varying-size white dots
- Constellation lines connect stars when enabled
- Text block shows message, coordinates, date below the circle
- Recognizable constellations visible (e.g., Orion in winter sky)

---

## Phase 3: Flask Application

**Goal**: Serve a web form that generates and displays the star map. Enable downloads.

### Tasks

1. **Implement `app.py`**
   - `GET /` -- form page
   - `POST /generate` -- render SVG, return page with inline preview
   - `POST /download/svg` -- SVG file download
   - `POST /download/png` -- PNG via CairoSVG (graceful fallback if unavailable)

2. **Implement `templates/index.html`**
   - Form: location name, lat, lon, date, time, timezone dropdown, constellation toggle
   - Preview area with inline SVG
   - Download buttons (appear after generation)
   - Form retains values after submission

3. **Implement `static/style.css`**
   - Dark theme (#0a0a0a background, light text)
   - Clean form inputs, subtle borders
   - Responsive: form top, preview below

4. **Validation**: lat/lon range, date/time parsing, timezone name check

### Files
- `app.py`
- `templates/index.html`
- `static/style.css`

### Done When
- `python app.py` serves on localhost:5000
- Form submits and shows SVG preview
- SVG and PNG download buttons work
- Invalid inputs show clear errors, form values retained
- Dark theme looks clean

---

## Phase 4: Polish

**Goal**: Visual quality, edge cases, UX.

### Tasks

1. **Visuals**: tune magnitude curves, add glow to brightest stars, render lines behind stars
2. **Constellation lines**: clip at horizon boundary, handle coord misalignment
3. **Text**: spacing tuning, long message handling, special character rendering
4. **Form UX**: loading indicator, sensible defaults, browser geolocation button
5. **PNG**: graceful CairoSVG fallback, test at multiple scales, verify font rendering
6. **Edge cases**: pole observer, daytime (empty sky), extreme dates
7. **Code quality**: type hints, docstrings, no hardcoded paths

### Files
- All `starmap/*.py` files, `app.py`, `templates/index.html`, `static/style.css`

### Done When
- Looks good for 3 test cases: northern winter, southern summer, equatorial
- PNG export clean at 2x scale with correct fonts
- No crashes on edge-case inputs
- Type hints and docstrings on all public functions
