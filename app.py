"""Flask web application for star map poster generation.

Routes:
    GET  /              -- Form page
    POST /generate      -- Generate star map, return page with inline SVG preview
    POST /download/svg  -- SVG file download
    POST /download/png  -- PNG via CairoSVG (graceful fallback)
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
    """Render the form page with no star map."""
    return render_template("index.html", timezones=TIMEZONES, svg=None, error=None, form_data={})


@app.route("/generate", methods=["POST"])
def generate() -> str:
    """Generate star map from form data, return page with inline SVG preview."""
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
    """Generate and download star map as SVG file."""
    result = _parse_form()
    if isinstance(result, str):
        return Response(result, status=400, mimetype="text/plain")
    lat, lon, dt, tz_name, message, show_constellations = result
    svg = render_star_map(lat, lon, dt, tz_name, message, show_constellations)
    return Response(svg, mimetype="image/svg+xml",
                    headers={"Content-Disposition": "attachment; filename=starmap.svg"})


@app.route("/download/png", methods=["POST"])
def download_png() -> Response:
    """Generate and download star map as PNG file via CairoSVG."""
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
