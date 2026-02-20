# CLAUDE.md — Star Map Generator

Local Python web app that generates astronomically accurate, print-quality star map posters. User enters a date, time, location, and custom message → app renders a beautiful SVG/PNG star map they can download.

## Reference Style

The output should look like a premium star map poster:
- **Black background** with a large **circular star field** in the upper portion
- **White/cream stars** as dots, sized proportionally to brightness (magnitude)
- **Constellation lines** connecting stars (thin, subtle white/gray lines)
- **Milky Way band** rendered as a soft, semi-transparent glow/cloud across the star field
- **Thin white circular border** around the star field
- **Thin white rectangular border** around the entire poster (with padding)
- **Text below the circle** in elegant small-caps serif font:
  - Line 1: Custom message (e.g. location name, "The night we met", etc.)
  - Line 2: Coordinates in DMS format (e.g. `42° 20' 25" N 71° 6' 12" W`)
  - Line 3: Date and time (e.g. `March 05, 2023, 12:00 PM`)
- Overall poster aspect ratio: roughly **3:4** (portrait, like 18×24 or 24×36 print)

## Tech Stack

- **Python 3.11+**
- **Flask** — localhost web server with a simple form UI
- **Astropy** — coordinate transforms (RA/Dec → Alt/Az), time handling, star catalog access
- **NumPy** — math
- **SVG output** (generated programmatically as XML strings, no matplotlib) — vector = infinite scale for print
- **Optional PNG export** via `cairosvg` (pip install cairosvg)

## How It Works

### Astronomy Pipeline

1. **Input**: date, time, timezone, latitude, longitude
2. **Convert to observation time**: Use `astropy.time.Time` with timezone → UTC
3. **Load star catalog**: Use the Yale Bright Star Catalog (BSC5, ~9,110 stars) or Hipparcos via `astropy.coordinates`. Stars down to about magnitude 6.5 (naked-eye visible).
4. **Coordinate transform**: For each star, convert RA/Dec (ICRS) → AltAz for the given time and location using `astropy.coordinates.AltAz`
5. **Filter**: Keep only stars above the horizon (altitude > 0°)
6. **Project**: Stereographic projection from the hemisphere onto a flat circle:
   - Center of circle = zenith (directly overhead)
   - Edge of circle = horizon
   - `r = cos(alt) / (1 + sin(alt))` normalized to circle radius
   - `theta = azimuth` (with North up, East left for proper sky orientation)
   - `x = center_x + r * sin(az)`, `y = center_y - r * cos(az)`
7. **Star sizing**: Map magnitude to dot radius. Brighter = bigger. Suggested: `radius = max(0.3, 2.5 - mag * 0.35)` (tweak as needed)
8. **Constellation lines**: Use IAU constellation stick figures (a known dataset of star pairs). Transform both endpoints the same way, draw a line if both are above horizon.
9. **Milky Way**: Approximate with a set of predefined galactic plane points transformed to the local sky, rendered as a subtle semi-transparent gradient or blurred path. A simpler approach: use a precomputed Milky Way outline (galactic latitude ≈ 0° band) and render it as a soft, wide, semi-transparent white stroke or filled region.

### SVG Rendering

Generate the SVG as a string (no external rendering library). Structure:

```xml
<svg xmlns="..." viewBox="0 0 1800 2400" width="1800" height="2400">
  <!-- Black background -->
  <rect width="100%" height="100%" fill="#000000"/>
  
  <!-- Outer border (thin white rectangle with margin) -->
  <rect x="40" y="40" width="1720" height="2320" fill="none" stroke="#ffffff" stroke-width="1" opacity="0.6"/>
  
  <!-- Star field clipping circle -->
  <defs>
    <clipPath id="sky-clip">
      <circle cx="900" cy="850" r="720"/>
    </clipPath>
  </defs>
  
  <!-- Milky Way (behind stars, inside clip) -->
  <g clip-path="url(#sky-clip)">
    <!-- Milky Way path with gaussian blur filter -->
  </g>
  
  <!-- Stars (inside clip) -->
  <g clip-path="url(#sky-clip)">
    <circle cx="..." cy="..." r="..." fill="#ffffff" opacity="..."/>
    <!-- ... more stars ... -->
  </g>
  
  <!-- Constellation lines (inside clip) -->
  <g clip-path="url(#sky-clip)">
    <line x1="..." y1="..." x2="..." y2="..." stroke="#ffffff" stroke-width="0.5" opacity="0.3"/>
  </g>
  
  <!-- Circle border -->
  <circle cx="900" cy="850" r="720" fill="none" stroke="#ffffff" stroke-width="1.5"/>
  
  <!-- Text below circle -->
  <text x="900" y="1750" text-anchor="middle" fill="#ffffff" font-family="'Georgia', 'Times New Roman', serif" font-size="28" letter-spacing="4">
    910 @ 400 FENWAY
  </text>
  <!-- ... more text lines ... -->
</svg>
```

Key SVG dimensions:
- Canvas: **1800 × 2400 px** (3:4 ratio, good for 12×16" at 150dpi or 6×8" at 300dpi)
- Star circle: centered horizontally, upper portion, radius ~720px
- Text area: below the circle, centered
- All text in small caps style (use CSS `font-variant: small-caps` or manually uppercase with smaller font for lowercase)

### Constellation Data

Use the standard IAU constellation line data. This is a well-known dataset that maps pairs of star identifiers (usually HIP numbers or HR numbers) to draw "stick figure" constellation outlines. 

A good source: the `constellationship.fab` file from Stellarium, which lists pairs of Hipparcos IDs for each constellation line segment. Alternatively, use a simplified JSON file with ~700 line segments mapping HR or HIP star pairs.

Bundle this as a JSON or CSV data file in the project:
```json
{
  "Orion": [[26727, 27366], [27366, 26311], ...],
  "Ursa Major": [[...], ...],
  ...
}
```

### Milky Way Rendering

The Milky Way is the hardest part visually. Approaches (from simplest to best):

**Simple (recommended to start):** Skip it. Constellation lines + dense star field already looks great.

**Medium:** Render the galactic plane as a wide, blurred, semi-transparent white band. Take a set of points along galactic latitude 0°, transform them to the observer's sky, and draw a thick (~100px), heavily blurred, low-opacity white path.

**Advanced:** Use a precomputed Milky Way outline (a polygon of the bright regions). Stellarium's `milkyway.dat` or similar datasets provide this. Render as a filled region with SVG gaussian blur and low opacity.

Start with "Simple" — get the stars and constellations right first, then add Milky Way as an enhancement.

## Project Structure

```
starmap-generator/
├── CLAUDE.md              # This file
├── requirements.txt       # Python dependencies
├── app.py                 # Flask app — routes, form handling, SVG generation entry point
├── starmap/
│   ├── __init__.py
│   ├── catalog.py         # Star catalog loading (BSC/Hipparcos via astropy)
│   ├── projection.py      # Coordinate transforms + stereographic projection
│   ├── constellations.py  # Constellation line data + lookup
│   ├── renderer.py        # SVG string generation
│   └── milkyway.py        # Milky Way rendering (optional, can be empty initially)
├── data/
│   └── constellations.json  # Constellation line pairs (HIP or HR IDs)
├── static/
│   └── style.css          # Minimal form styling
└── templates/
    └── index.html         # Flask template — form + preview
```

## Web UI (localhost)

Simple single-page Flask app:

### Form Fields
- **Location name** (text) — displayed on poster (e.g. "910 @ 400 Fenway")
- **Latitude** (number, decimal degrees) — e.g. 42.3404
- **Longitude** (number, decimal degrees) — e.g. -71.1028
- **Date** (date picker)
- **Time** (time picker, 12h or 24h)
- **Timezone** (dropdown — common US timezones + UTC, or free text like "America/New_York")
- **Custom message** (optional text — if provided, replaces location name as first line)
- **Show constellation lines** (checkbox, default on)
- **Show Milky Way** (checkbox, default on if implemented)

### Output
- "Generate" button → renders SVG inline on the page as a preview
- "Download SVG" button → serves the SVG as a file download
- "Download PNG" button → converts SVG to high-res PNG via cairosvg (if available)

### Styling
- Dark theme (matches the poster vibe)
- Form on the left or top, preview on the right or bottom
- Keep it simple — this is a tool, not a product

## Commands

```bash
# Setup
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Run
python app.py
# → http://localhost:5000

# Dependencies
pip install flask astropy numpy cairosvg
```

## Requirements.txt

```
flask>=3.0
astropy>=6.0
numpy>=1.26
cairosvg>=2.7
```

## Key Implementation Notes

1. **Star catalog**: `astropy` can query the Hipparcos catalog via `astropy.coordinates.SkyCoord` and vizier, but for offline use, download the BSC5 catalog as a CSV/JSON and bundle it. Alternatively, use `astropy.coordinates.get_constellation` with a hardcoded bright star list. The simplest reliable approach: bundle a JSON of the ~9,000 BSC stars with RA, Dec, magnitude, and HR number.

2. **Performance**: The coordinate transform for ~9,000 stars takes a few seconds with astropy. This is fine for a one-shot generator. Cache the catalog load.

3. **Text styling**: SVG `font-variant: small-caps` may not render in all viewers. Safer approach: manually format text as uppercase with CSS `letter-spacing` for the small-caps look.

4. **Coordinate display**: Convert decimal lat/lon to DMS for the poster text: `42° 20' 25" N 71° 6' 12" W`

5. **Sky orientation**: Traditional star maps show East on the left (mirror of ground maps) because you're looking UP. Make sure azimuth mapping reflects this: `x = center + r * sin(az)` with az measured clockwise from North.

6. **Star opacity**: Optionally vary opacity with magnitude — dimmer stars slightly more transparent. Gives a more natural depth feel.

7. **PNG export size**: When using cairosvg, scale up for print quality. For a 24×36" print at 300dpi: 7200×10800px. The SVG viewBox stays at 1800×2400, cairosvg handles the scaling.

## Stretch Goals (Not Required for v1)

- Heart-shaped mask option (clip stars to a heart path instead of a circle)
- Color themes (dark on light, navy blue, deep purple)
- Moon position rendered as a crescent
- Planet positions labeled
- Grid lines (RA/Dec or Alt/Az)
- Geocoding (type a city name → auto lat/lon lookup)
- Font selection
