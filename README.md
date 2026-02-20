# Star Map Generator

Generate astronomically accurate star map posters for any date, time, and location. Renders a print-quality SVG (or PNG) showing the night sky with constellation lines, styled as a minimal poster.

## Features

- Accurate star positions using the HYG star catalog (~9,000 naked-eye stars)
- 88 constellation line overlays (toggleable)
- Stereographic projection centered on the zenith
- Custom message, coordinates in DMS format, and formatted date/time
- SVG and PNG download (PNG at 3600x4800px for print)
- Dark-themed responsive web UI

## Setup

Requires Python 3.11+ and (optionally) the Cairo C library for PNG export.

```bash
# Clone and set up
git clone https://github.com/georgenijo/star-maker.git
cd star-maker
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Optional: install Cairo for PNG export
# macOS:
brew install cairo
# Ubuntu:
sudo apt install libcairo2-dev
```

## Usage

```bash
source venv/bin/activate
python app.py
```

Open `http://localhost:8080` in your browser. Fill in a location (lat/lon), date, time, timezone, and an optional message, then click **Generate Star Map**.

Download the result as SVG (vector, scalable to any print size) or PNG.

## Project Structure

```
starmap/
  catalog.py         # Star catalog loader (HYG database)
  projection.py      # ICRS -> AltAz -> stereographic projection
  constellations.py  # Constellation line data loader
  renderer.py        # SVG poster generation
app.py               # Flask web application
templates/
  index.html         # Form and preview template
static/
  style.css          # Dark theme styles
data/
  stars.json         # Pre-processed star catalog (~9k stars, mag <= 6.5)
  constellations.json # d3-celestial constellation lines
scripts/
  prepare_stars.py   # One-off script to regenerate stars.json from HYG CSV
```

## Data Sources

- **Stars**: [HYG Database v4.1](https://www.astronexus.com/hyg) (CC BY-SA 4.0)
- **Constellation lines**: [d3-celestial](https://github.com/ofrohn/d3-celestial) (MIT)

## License

MIT
