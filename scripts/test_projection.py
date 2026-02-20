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
            print(f"Polaris: alt={p_alt:.1f} deg, az={p_az:.1f} deg")
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
