"""SVG star map renderer.

Generates print-quality star map poster as SVG string.
Canvas: 1800x2400 px (3:4 portrait).
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

# -- Layout Constants --
CANVAS_W = 1800
CANVAS_H = 2400
CIRCLE_CX = 900
CIRCLE_CY = 850
CIRCLE_R = 720
BORDER_INSET = 40
BORDER_OPACITY = 0.6

# -- Text Constants --
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
    """Star magnitude -> SVG circle radius.

    Brighter (lower mag) = bigger. Range: ~0.3px to ~3.9px.
    """
    return max(0.3, 2.5 - mag * 0.35)


def mag_to_opacity(mag: float) -> float:
    """Star magnitude -> opacity. Brighter = more opaque.

    Range: 0.5 (mag 6.5) to 1.0 (mag <= 2).
    """
    return max(0.5, min(1.0, 1.0 - (mag - 2.0) * 0.1))


def decimal_to_dms(value: float, is_lat: bool) -> str:
    """Decimal degrees -> DMS string like '42 deg 20' 25" N'.

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
    """Format datetime for poster: 'March 15, 2023, 09:00 PM'.

    Args:
        dt: Naive datetime (no tzinfo).
        tz_name: IANA timezone name, e.g. 'America/New_York'.

    Returns:
        Formatted date/time string.
    """
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
        4. Assemble SVG: background -> border -> clipPath -> constellations -> stars -> circle -> text

    Args:
        lat: Observer latitude in decimal degrees (-90 to +90).
        lon: Observer longitude in decimal degrees (-180 to +180).
        dt: Naive datetime for the observation.
        tz_name: IANA timezone name, e.g. 'America/New_York'.
        message: Custom message text for the poster.
        show_constellations: Whether to render constellation lines.

    Returns:
        Complete SVG XML string.
    """
    # -- 1. Build frame --
    frame = build_altaz_frame(lat, lon, dt, tz_name)

    # -- 2. Stars --
    stars = load_stars()
    ra = np.array([s["ra"] for s in stars])
    dec = np.array([s["dec"] for s in stars])
    mags = np.array([s["mag"] for s in stars])

    alt, az = transform_radec_to_altaz(ra, dec, frame)
    mask = filter_visible(alt)
    vis_alt, vis_az, vis_mag = alt[mask], az[mask], mags[mask]
    sx, sy = stereographic_project(vis_alt, vis_az, CIRCLE_CX, CIRCLE_CY, CIRCLE_R)

    # -- 3. Constellation lines (BATCHED -- transform once) --
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

    # -- 4. Assemble SVG --
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
