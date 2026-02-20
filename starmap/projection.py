"""Coordinate transforms and stereographic projection.

Pipeline: ICRS (RA/Dec) -> AltAz (altitude/azimuth) -> Stereographic (x, y)
"""
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np
from numpy.typing import NDArray
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u


def build_altaz_frame(
    lat: float,
    lon: float,
    dt: datetime,
    tz_name: str,
) -> AltAz:
    """Build AltAz reference frame for a specific time and location.

    Args:
        lat: decimal degrees, -90 to +90
        lon: decimal degrees, -180 to +180
        dt: naive datetime (no tzinfo)
        tz_name: IANA tz name, e.g. "America/New_York"
    """
    tz = ZoneInfo(tz_name)
    local_dt = dt.replace(tzinfo=tz)
    obs_time = Time(local_dt)
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    return AltAz(obstime=obs_time, location=location)


def transform_radec_to_altaz(
    ra: NDArray[np.float64],
    dec: NDArray[np.float64],
    frame: AltAz,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Bulk-transform RA/Dec arrays to altitude/azimuth arrays (degrees).

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
    alt: NDArray[np.float64],
    az: NDArray[np.float64],
    cx: float,
    cy: float,
    radius: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Stereographic projection from hemisphere to flat circle.

    Center = zenith (alt=90 deg), Edge = horizon (alt=0 deg).
    East appears on LEFT (correct for looking up at the sky).

    Math:
        alt_rad = radians(alt)
        az_rad  = radians(az)
        r_norm  = cos(alt_rad) / (1 + sin(alt_rad))
        x = cx + r_norm * radius * sin(az_rad)
        y = cy - r_norm * radius * cos(az_rad)
    """
    alt_rad = np.radians(alt)
    az_rad = np.radians(az)
    r_norm = np.cos(alt_rad) / (1.0 + np.sin(alt_rad))
    x = cx + r_norm * radius * np.sin(az_rad)
    y = cy - r_norm * radius * np.cos(az_rad)
    return x, y
