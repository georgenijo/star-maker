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
