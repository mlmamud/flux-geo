#!/usr/bin/env python3
"""Plot and animate PFLOTRAN multi-well outputs.

This script is intentionally lightweight (h5py + numpy + matplotlib) and works
straight from PFLOTRAN's XMF/HDF5 outputs.

What it can do:
- List available variables in `pflotran.h5`
- Plot 2D slices for any cell-centered variable:
    - XY slice at constant Z
    - XZ slice at constant Y (vertical cross section)
    - YZ slice at constant X (vertical cross section)
- Animate that slice over time (GIF or MP4)
- Plot a time-series for a variable at the nearest cell to a requested (x,y,z)

Typical usage (from this folder):

  # list variables
  python plot_outputs.py --list-vars

  # plot gas saturation on a Z-slice near 105 m
  python plot_outputs.py --var "Gas Saturation" --z 105 --outdir figures

  # animate gas saturation (every 2nd output)
  python plot_outputs.py --var "Gas Saturation" --z 105 --every 2 --animate --gif

    # animation with a STATIC colorbar (global min/max across all frames)
    python plot_outputs.py --var "Gas Saturation" --z 105 --animate --gif --scale global

    # plot ALL available variables (last-time snapshot for each)
    python plot_outputs.py --all-vars --z 105 --outdir figures_allvars

    # vertical cross sections
    python plot_outputs.py --plane xz --y 775 --var "Gas Saturation" --outdir figures_xz
    python plot_outputs.py --plane yz --x 775 --var "Gas Saturation" --outdir figures_yz

    # batch across ALL scenarios under a folder (finds all pflotran.h5):
    # (useful once you have multiple cases run)
    python plot_outputs.py --all-scenarios --search-root ../ --all-vars --outdir figures_batch

  # time-series at a point
  python plot_outputs.py --var "Liquid Pressure [Pa]" --point 775 775 105 --timeseries


cd /home/lal/wsl_codes/Carbon/co2/sco2/multi-well && rm -rf figures_gif_xy_slow_jet_test && python3 plot_outputs.py --plane xy --var "Gas Pressure [Pa]" --animate --gif --outdir figures_gif_xy_slow_jet_test --method grid --every 1
Notes:
- Coordinates (XC, YC, ZC) are cell centers.
- Z is in meters (consistent with the mesh).
"""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np


def _lazy_import_matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def _lazy_import_mpl_tri():
    import matplotlib.tri as mtri

    return mtri


def default_cmap(var: str, delta: bool) -> str:
    if delta:
        return "RdBu_r"
    v = var.lower()
    if "saturation" in v:
        return "turbo"
    if "pressure" in v:
        return "jet"
    if "temperature" in v:
        return "inferno"
    if "permeability" in v:
        return "turbo"
    if "density" in v:
        return "turbo"
    return "turbo"


@dataclass(frozen=True)
class TimeGroup:
    name: str
    index: int
    time_years: float


# PFLOTRAN HDF5 has a couple common top-level time group naming conventions.
# 1) (older) "  0 Time  1.00000E+00 y"
# 2) (common) "Time:  1.00000E+00 y"
_TIME_RE_OLD = re.compile(r"^\s*(\d+)\s+Time\s+([0-9.Ee+-]+)\s+y\s*$")
_TIME_RE_NEW = re.compile(r"^\s*Time:\s*([0-9.Ee+-]+)\s+y\s*$")


def parse_time_group_name(name: str) -> Optional[TimeGroup]:
    match_old = _TIME_RE_OLD.match(name)
    if match_old:
        return TimeGroup(name=name, index=int(match_old.group(1)), time_years=float(match_old.group(2)))
    match_new = _TIME_RE_NEW.match(name)
    if match_new:
        # index will be assigned in list_time_groups based on sorted time
        return TimeGroup(name=name, index=-1, time_years=float(match_new.group(1)))
    return None


def list_time_groups(h5) -> list[TimeGroup]:
    groups: list[TimeGroup] = []
    for key in h5.keys():
        tg = parse_time_group_name(key)
        if tg is not None:
            groups.append(tg)
    # Prefer sorting by time (more robust across naming conventions).
    groups.sort(key=lambda t: (float(t.time_years), str(t.name)))
    # Re-index sequentially so filenames/time selection remain stable.
    groups = [TimeGroup(name=t.name, index=i, time_years=t.time_years) for i, t in enumerate(groups)]
    return groups


def read_cell_centers(h5, sample_shape: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return cell center coordinates as flat arrays (ncell,).

    Supports two common PFLOTRAN layouts:
    - Unstructured: Domain/XC, Domain/YC, Domain/ZC are (ncell,)
    - Structured: Coordinates/X [m], Y [m], Z [m] are node coordinates (nx+1, ny+1, nz+1)
      and fields are (nx, ny, nz)
    """

    if "Domain" in h5 and all(k in h5["Domain"] for k in ("XC", "YC", "ZC")):
        xc = np.asarray(h5["Domain"]["XC"], dtype=float).ravel()
        yc = np.asarray(h5["Domain"]["YC"], dtype=float).ravel()
        zc = np.asarray(h5["Domain"]["ZC"], dtype=float).ravel()
        return xc, yc, zc

    if "Coordinates" not in h5:
        raise KeyError(
            "Could not find coordinates in HDF5. Expected either Domain/(XC,YC,ZC) or Coordinates/(X [m],Y [m],Z [m])."
        )

    coords = h5["Coordinates"]
    if not all(k in coords for k in ("X [m]", "Y [m]", "Z [m]")):
        raise KeyError(
            "Coordinates group exists, but missing one of: 'X [m]', 'Y [m]', 'Z [m]'."
        )

    x_edges = np.asarray(coords["X [m]"], dtype=float).ravel()
    y_edges = np.asarray(coords["Y [m]"], dtype=float).ravel()
    z_edges = np.asarray(coords["Z [m]"], dtype=float).ravel()

    if len(sample_shape) != 3:
        raise ValueError(
            f"Structured grid detected via Coordinates, but sample field has shape {sample_shape}; expected 3D (nx, ny, nz)."
        )
    nx, ny, nz = (int(sample_shape[0]), int(sample_shape[1]), int(sample_shape[2]))
    if x_edges.size != nx + 1 or y_edges.size != ny + 1 or z_edges.size != nz + 1:
        raise ValueError(
            "Coordinates edge arrays do not match field shape. "
            f"Got X:{x_edges.size},Y:{y_edges.size},Z:{z_edges.size} but field shape (nx,ny,nz)=({nx},{ny},{nz})."
        )

    x_cent = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_cent = 0.5 * (y_edges[:-1] + y_edges[1:])
    z_cent = 0.5 * (z_edges[:-1] + z_edges[1:])
    Xc, Yc, Zc = np.meshgrid(x_cent, y_cent, z_cent, indexing="ij")
    return Xc.ravel(), Yc.ravel(), Zc.ravel()


def plot_timeseries_multi(
    *,
    times: np.ndarray,
    series: dict[str, np.ndarray],
    title: str,
    ylabel: str,
    outpath: str,
):
    plt = _lazy_import_matplotlib()

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    for label, y in series.items():
        ax.plot(times, y, lw=2, label=label)
    ax.set_xlabel("Time [years]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def estimate_default_slice(zc: np.ndarray) -> tuple[float, float]:
    """Pick a reasonable default Z slice and tolerance.

    - z0 = median ZC
    - tol = ~half of the minimum spacing between unique Z levels
    """

    z_unique = np.unique(np.asarray(zc, dtype=float))
    z_unique.sort()

    z0 = float(np.median(z_unique))

    dz = np.diff(z_unique)
    dz = dz[dz > 0]
    if dz.size == 0:
        return z0, 1.0

    tol = float(0.51 * np.min(dz))
    return z0, tol


def estimate_default_plane(coord: np.ndarray) -> tuple[float, float]:
    """Pick a reasonable default plane location and tolerance for X/Y/Z.

    - c0 = median(coord)
    - tol = ~half of minimum spacing between unique coordinate levels
    """

    c_unique = np.unique(np.asarray(coord, dtype=float))
    c_unique.sort()

    c0 = float(np.median(c_unique))

    dc = np.diff(c_unique)
    dc = dc[dc > 0]
    if dc.size == 0:
        return c0, 1.0

    tol = float(0.51 * np.min(dc))
    return c0, tol


def select_cells_by_z(zc: np.ndarray, z: float, ztol: float) -> np.ndarray:
    return np.abs(zc - z) <= ztol


def select_cells_by_x(xc: np.ndarray, x: float, xtol: float) -> np.ndarray:
    return np.abs(xc - x) <= xtol


def select_cells_by_y(yc: np.ndarray, y: float, ytol: float) -> np.ndarray:
    return np.abs(yc - y) <= ytol


def maybe_downsample(mask: np.ndarray, max_points: int, rng: np.random.Generator) -> np.ndarray:
    idx = np.flatnonzero(mask)
    if max_points is None or idx.size <= max_points:
        return mask
    keep = rng.choice(idx, size=max_points, replace=False)
    new_mask = np.zeros_like(mask, dtype=bool)
    new_mask[keep] = True
    return new_mask


def get_available_vars(h5) -> list[str]:
    groups = list_time_groups(h5)
    if not groups:
        raise RuntimeError("No time groups found in HDF5 (expected keys like '  0 Time  ... y').")
    first = groups[0].name
    return sorted(list(h5[first].keys()))


def find_h5_files(search_root: str, filename: str = "pflotran.h5") -> list[str]:
    paths: list[str] = []
    for root, _, files in os.walk(search_root):
        if filename in files:
            paths.append(os.path.join(root, filename))
    paths.sort()
    return paths


def read_field(h5, group_name: str, var: str) -> np.ndarray:
    if group_name not in h5:
        raise KeyError(f"Time group not found: {group_name!r}")
    g = h5[group_name]
    if var not in g:
        available = sorted(g.keys())
        raise KeyError(
            f"Variable {var!r} not found in {group_name!r}. Available: {available}"
        )
    return np.asarray(g[var])


def nearest_cell_index(xc: np.ndarray, yc: np.ndarray, zc: np.ndarray, point: tuple[float, float, float]) -> int:
    x0, y0, z0 = point
    dx = xc - x0
    dy = yc - y0
    dz = zc - z0
    # distance in meters
    dist2 = dx * dx + dy * dy + dz * dz
    return int(np.argmin(dist2))


def read_well_xy_from_well_file(path: str) -> Optional[tuple[float, float]]:
    """Parse (x,y) from PFLOTRAN .well file header.

    Expected line like:
      Top of hole (x,y,z) [m]:  7.750E+02 7.750E+02 2.000E+02
    """

    if not os.path.exists(path):
        return None
    top_re = re.compile(r"^\s*Top of hole\s*\(x,y,z\)\s*\[m\]:\s*([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s*$")
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for _ in range(200):
                line = f.readline()
                if not line:
                    break
                m = top_re.match(line)
                if m:
                    return float(m.group(1)), float(m.group(2))
    except OSError:
        return None
    return None


def read_well_path_from_well_file(path: str) -> Optional[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Parse a polyline of well segment centers (x,y,z) from a PFLOTRAN .well file.

    Expected lines like:
      Segment center (x,y,z) [m]:  7.750E+02 7.750E+02 2.500E+01

    Returns (x, y, z) arrays, or None if not found.
    """

    if not os.path.exists(path):
        return None
    seg_re = re.compile(
        r"^\s*Segment center\s*\(x,y,z\)\s*\[m\]:\s*([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s*$"
    )
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for _ in range(400):
                line = f.readline()
                if not line:
                    break
                m = seg_re.match(line)
                if m:
                    xs.append(float(m.group(1)))
                    ys.append(float(m.group(2)))
                    zs.append(float(m.group(3)))
    except OSError:
        return None
    if not xs:
        return None
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float), np.asarray(zs, dtype=float)


def _try_build_grid(x: np.ndarray, y: np.ndarray):
    """Try to interpret scattered (x,y) as a rectilinear grid.

    Returns (X2, Y2, ij) where:
      - X2, Y2 are 2D arrays from meshgrid
      - ij is an (n,2) array with integer indices mapping each point to (i,j)
    Returns None if the points do not lie on a full grid.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    x_unique.sort()
    y_unique.sort()

    nx = x_unique.size
    ny = y_unique.size
    if nx * ny != x.size:
        return None

    # Map each point to grid indices
    ix = np.searchsorted(x_unique, x)
    iy = np.searchsorted(y_unique, y)

    # Ensure each (iy,ix) appears exactly once
    flat = iy * nx + ix
    if np.unique(flat).size != flat.size:
        return None

    X2, Y2 = np.meshgrid(x_unique, y_unique)
    ij = np.vstack([iy, ix]).T
    return X2, Y2, ij


def _auto_vrange(
    values: np.ndarray,
    mask: np.ndarray,
    robust: bool,
    robust_percentiles: tuple[float, float] = (2.0, 98.0),
) -> tuple[float, float]:
    v = np.asarray(values)[mask]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0, 1.0
    if robust:
        p_lo, p_hi = robust_percentiles
        lo, hi = np.percentile(v, [p_lo, p_hi])
        # If the field is sparse (e.g., mostly zeros) these percentiles can collapse.
        # In that case, recompute on nonzero values to keep the colorbar responsive.
        if np.isclose(lo, hi):
            eps = 1e-30
            v_nz = v[np.abs(v) > eps]
            if v_nz.size >= max(50, int(0.001 * v.size)):
                lo, hi = np.percentile(v_nz, [p_lo, p_hi])
            if np.isclose(lo, hi):
                lo, hi = float(np.min(v)), float(np.max(v))
    else:
        lo, hi = float(np.min(v)), float(np.max(v))
    if np.isclose(lo, hi):
        hi = lo + 1.0
    return float(lo), float(hi)


def plot_slice(
    *,
    xc: np.ndarray,
    yc: np.ndarray,
    values: np.ndarray,
    mask: np.ndarray,
    title: str,
    outpath: str,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    s: float,
    method: str,
    levels: int,
    wells_xy: Optional[list[tuple[str, float, float]]],
    cbar_label: str,
):
    plt = _lazy_import_matplotlib()

    fig, ax = plt.subplots(figsize=(9, 6.8), constrained_layout=True)

    xm = xc[mask]
    ym = yc[mask]
    vm = values[mask]

    mappable = None
    if method == "grid":
        grid = _try_build_grid(xm, ym)
        if grid is None:
            method = "tri"
        else:
            X2, Y2, ij = grid
            V2 = np.full(X2.shape, np.nan, dtype=float)
            V2[ij[:, 0], ij[:, 1]] = vm
            mappable = ax.pcolormesh(X2, Y2, V2, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)

    if mappable is None and method == "tri":
        mtri = _lazy_import_mpl_tri()
        tri = mtri.Triangulation(xm, ym)
        mappable = ax.tricontourf(tri, vm, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)

    if mappable is None:
        mappable = ax.scatter(xm, ym, c=vm, s=s, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)

    if wells_xy:
        for name, wx, wy in wells_xy:
            ax.plot(wx, wy, marker="*", markersize=12, color="white", markeredgecolor="black", markeredgewidth=1.0)
            ax.text(wx + 15, wy + 15, name, color="black", fontsize=10, weight="bold")

    ax.set_title(title)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_aspect("equal", adjustable="box")

    cb = fig.colorbar(mappable, ax=ax)
    cb.set_label(cbar_label)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_slice_2d(
    *,
    u: np.ndarray,
    v: np.ndarray,
    values: np.ndarray,
    mask: np.ndarray,
    title: str,
    outpath: str,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    s: float,
    method: str,
    levels: int,
    overlays: Optional[list[tuple[str, np.ndarray, np.ndarray, dict]]],
    xlabel: str,
    ylabel: str,
    cbar_label: str,
    equal_aspect: bool,
):
    plt = _lazy_import_matplotlib()

    fig, ax = plt.subplots(figsize=(9, 6.8), constrained_layout=True)

    um = u[mask]
    vm = v[mask]
    valm = values[mask]

    mappable = None
    if method == "grid":
        grid = _try_build_grid(um, vm)
        if grid is None:
            method = "tri"
        else:
            U2, V2, ij = grid
            W2 = np.full(U2.shape, np.nan, dtype=float)
            W2[ij[:, 0], ij[:, 1]] = valm
            mappable = ax.pcolormesh(U2, V2, W2, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)

    if mappable is None and method == "tri":
        mtri = _lazy_import_mpl_tri()
        tri = mtri.Triangulation(um, vm)
        mappable = ax.tricontourf(tri, valm, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)

    if mappable is None:
        mappable = ax.scatter(um, vm, c=valm, s=s, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)

    if overlays:
        for name, ou, ov, style in overlays:
            style = dict(style)
            ax.plot(ou, ov, **style)
            if ou.size and ov.size:
                ax.text(float(ou[-1]) + 10, float(ov[-1]) + 10, name, color="black", fontsize=10, weight="bold")

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")

    cb = fig.colorbar(mappable, ax=ax)
    cb.set_label(cbar_label)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_timeseries(
    *,
    times: np.ndarray,
    series: np.ndarray,
    title: str,
    ylabel: str,
    outpath: str,
):
    plt = _lazy_import_matplotlib()

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(times, series, lw=2)
    ax.set_xlabel("Time [years]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def animate_slice(
    *,
    xc: np.ndarray,
    yc: np.ndarray,
    mask: np.ndarray,
    time_groups: list[TimeGroup],
    load_values,
    title_prefix: str,
    outpath: str,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    s: float,
    fps: int,
    gif: bool,
    method: str,
    levels: int,
    wells_xy: Optional[list[tuple[str, float, float]]],
    cbar_label: str,
    scale: str,
    robust: bool,
    scale_alpha: float,
):
    plt = _lazy_import_matplotlib()
    import matplotlib.animation as animation

    xm = xc[mask]
    ym = yc[mask]

    # Prefer grid animation for speed/quality when possible
    grid = _try_build_grid(xm, ym) if method == "grid" else None
    if method == "grid" and grid is None:
        method = "tri"

    mtri = _lazy_import_mpl_tri() if method == "tri" else None
    tri = mtri.Triangulation(xm, ym) if method == "tri" else None
    contour_artists: list = []

    def _levels_for(vmin_use: Optional[float], vmax_use: Optional[float]):
        if not isinstance(levels, int):
            return levels
        if vmin_use is None or vmax_use is None:
            return levels
        lo = float(vmin_use)
        hi = float(vmax_use)
        if not np.isfinite(lo) or not np.isfinite(hi):
            return levels
        if np.isclose(lo, hi):
            hi = lo + 1.0
        return np.linspace(lo, hi, int(levels))

    fig, ax = plt.subplots(figsize=(8, 6))
    values0, t0 = load_values(time_groups[0])

    # If we're auto-scaling per-frame, don't let the initial draw pick a tiny
    # +/-eps range when early frames are nearly constant.
    init_vmin = vmin
    init_vmax = vmax
    current_vmin = vmin
    current_vmax = vmax
    if scale in {"frame", "smooth"} and (vmin is None or vmax is None):
        init_vmin, init_vmax = _auto_vrange(values0, mask, robust=robust)
        current_vmin, current_vmax = init_vmin, init_vmax

    vm0 = values0[mask]

    mappable = None
    if method == "grid" and grid is not None:
        X2, Y2, ij = grid
        V2 = np.full(X2.shape, np.nan, dtype=float)
        V2[ij[:, 0], ij[:, 1]] = vm0
        mappable = ax.pcolormesh(X2, Y2, V2, shading="auto", cmap=cmap, vmin=init_vmin, vmax=init_vmax)
        grid_state = (X2, Y2, ij, V2)
    elif method == "tri" and tri is not None:
        before = list(ax.collections)
        mappable = ax.tricontourf(
            tri,
            vm0,
            levels=_levels_for(init_vmin, init_vmax),
            cmap=cmap,
            vmin=init_vmin,
            vmax=init_vmax,
        )
        after = list(ax.collections)
        # store only the artists created by tricontourf
        contour_artists = [c for c in after if c not in before]
        grid_state = None
    else:
        mappable = ax.scatter(xm, ym, c=vm0, s=s, cmap=cmap, vmin=init_vmin, vmax=init_vmax, linewidths=0)
        grid_state = None

    if wells_xy:
        for name, wx, wy in wells_xy:
            ax.plot(wx, wy, marker="*", markersize=12, color="white", markeredgecolor="black", markeredgewidth=1.0)
            ax.text(wx + 15, wy + 15, name, color="black", fontsize=10, weight="bold")

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(mappable, ax=ax)
    cb.set_label(cbar_label)
    # Stable tick formatting for wide-range fields (e.g., pressure).
    # IMPORTANT: don't force explicit tick arrays per-frame; that can produce
    # overlapping/duplicated labels with scientific offset text.
    try:
        import matplotlib.ticker as mticker

        cb.locator = mticker.MaxNLocator(nbins=7)
        cb.formatter = mticker.ScalarFormatter(useMathText=True)
        cb.formatter.set_powerlimits((-3, 3))
        cb.update_ticks()
    except Exception:
        pass

    # Two-line title prevents overlap with colorbar scientific offset text.
    title = ax.set_title(f"{title_prefix}\nt = {t0:.3g} y")

    def update(frame_idx: int):
        nonlocal current_vmin, current_vmax
        tg = time_groups[frame_idx]
        values, t = load_values(tg)

        vm = values[mask]
        # Optional rescaling of the colormap per-frame.
        if scale in {"frame", "smooth"}:
            new_vmin, new_vmax = _auto_vrange(values, mask, robust=robust)
            if scale == "smooth":
                a = float(np.clip(scale_alpha, 0.0, 1.0))
                if current_vmin is None:
                    current_vmin = new_vmin
                else:
                    current_vmin = a * current_vmin + (1.0 - a) * new_vmin
                if current_vmax is None:
                    current_vmax = new_vmax
                else:
                    current_vmax = a * current_vmax + (1.0 - a) * new_vmax
            else:
                current_vmin, current_vmax = new_vmin, new_vmax

            try:
                mappable.set_clim(current_vmin, current_vmax)
                cb.update_normal(mappable)
                cb.update_ticks()
            except Exception:
                pass
        if method == "grid" and grid_state is not None:
            _, _, ij, V2 = grid_state
            V2[:] = np.nan
            V2[ij[:, 0], ij[:, 1]] = vm
            mappable.set_array(V2.ravel())
        elif method == "tri" and tri is not None:
            # tricontourf needs redraw; remove only the previous contour artists
            nonlocal contour_artists
            for artist in contour_artists:
                try:
                    artist.remove()
                except Exception:
                    pass
            before = list(ax.collections)
            # If we're scaling per-frame, use the updated clim.
            vmin_use = current_vmin if scale in {"frame", "smooth"} else vmin
            vmax_use = current_vmax if scale in {"frame", "smooth"} else vmax
            new_contour = ax.tricontourf(
                tri,
                vm,
                levels=_levels_for(vmin_use, vmax_use),
                cmap=cmap,
                vmin=vmin_use,
                vmax=vmax_use,
            )
            after = list(ax.collections)
            contour_artists = [c for c in after if c not in before]
            try:
                cb.update_normal(new_contour)
                cb.update_ticks()
            except Exception:
                pass
        else:
            mappable.set_array(vm)

        if scale in {"frame", "smooth"} and current_vmin is not None and current_vmax is not None:
            title.set_text(
                f"{title_prefix}\n"
                f"t = {t:.3g} y | range=[{current_vmin:.3g}, {current_vmax:.3g}]"
            )
        else:
            title.set_text(f"{title_prefix}\nt = {t:.3g} y")
        return (title,)

    ani = animation.FuncAnimation(fig, update, frames=len(time_groups), interval=1000 / max(fps, 1), blit=False)

    if gif:
        writer = animation.PillowWriter(fps=fps)
    else:
        if animation.writers.is_available("ffmpeg"):
            writer = animation.FFMpegWriter(fps=fps)
        else:
            # fallback to GIF if ffmpeg is missing
            writer = animation.PillowWriter(fps=fps)
            if not outpath.lower().endswith(".gif"):
                outpath = os.path.splitext(outpath)[0] + ".gif"

    ani.save(outpath, writer=writer, dpi=150)
    plt.close(fig)


def animate_slice_2d(
    *,
    u: np.ndarray,
    v: np.ndarray,
    mask: np.ndarray,
    time_groups: list[TimeGroup],
    load_values,
    title_prefix: str,
    outpath: str,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    s: float,
    fps: int,
    gif: bool,
    method: str,
    levels: int,
    overlays: Optional[list[tuple[str, np.ndarray, np.ndarray, dict]]],
    xlabel: str,
    ylabel: str,
    cbar_label: str,
    scale: str,
    robust: bool,
    scale_alpha: float,
    equal_aspect: bool,
):
    plt = _lazy_import_matplotlib()
    import matplotlib.animation as animation

    um = u[mask]
    vm = v[mask]

    grid = _try_build_grid(um, vm) if method == "grid" else None
    if method == "grid" and grid is None:
        method = "tri"

    mtri = _lazy_import_mpl_tri() if method == "tri" else None
    tri = mtri.Triangulation(um, vm) if method == "tri" else None
    contour_artists: list = []

    def _levels_for(vmin_use: Optional[float], vmax_use: Optional[float]):
        if not isinstance(levels, int):
            return levels
        if vmin_use is None or vmax_use is None:
            return levels
        lo = float(vmin_use)
        hi = float(vmax_use)
        if not np.isfinite(lo) or not np.isfinite(hi):
            return levels
        if np.isclose(lo, hi):
            hi = lo + 1.0
        return np.linspace(lo, hi, int(levels))

    fig, ax = plt.subplots(figsize=(8, 6))
    values0, t0 = load_values(time_groups[0])
    w0 = values0[mask]

    # If we're auto-scaling per-frame, don't let the initial draw pick a tiny
    # +/-eps range when early frames are nearly constant.
    init_vmin = vmin
    init_vmax = vmax
    current_vmin = vmin
    current_vmax = vmax
    if scale in {"frame", "smooth"} and (vmin is None or vmax is None):
        init_vmin, init_vmax = _auto_vrange(values0, mask, robust=robust)
        current_vmin, current_vmax = init_vmin, init_vmax

    mappable = None
    if method == "grid" and grid is not None:
        U2, V2, ij = grid
        W2 = np.full(U2.shape, np.nan, dtype=float)
        W2[ij[:, 0], ij[:, 1]] = w0
        mappable = ax.pcolormesh(U2, V2, W2, shading="auto", cmap=cmap, vmin=init_vmin, vmax=init_vmax)
        grid_state = (U2, V2, ij, W2)
    elif method == "tri" and tri is not None:
        before = list(ax.collections)
        mappable = ax.tricontourf(
            tri,
            w0,
            levels=_levels_for(init_vmin, init_vmax),
            cmap=cmap,
            vmin=init_vmin,
            vmax=init_vmax,
        )
        after = list(ax.collections)
        contour_artists = [c for c in after if c not in before]
        grid_state = None
    else:
        mappable = ax.scatter(um, vm, c=w0, s=s, cmap=cmap, vmin=init_vmin, vmax=init_vmax, linewidths=0)
        grid_state = None

    if overlays:
        for name, ou, ov, style in overlays:
            style = dict(style)
            ax.plot(ou, ov, **style)
            if ou.size and ov.size:
                ax.text(float(ou[-1]) + 10, float(ov[-1]) + 10, name, color="black", fontsize=10, weight="bold")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(mappable, ax=ax)
    cb.set_label(cbar_label)
    try:
        import matplotlib.ticker as mticker

        cb.locator = mticker.MaxNLocator(nbins=7)
        cb.formatter = mticker.ScalarFormatter(useMathText=True)
        cb.formatter.set_powerlimits((-3, 3))
        cb.update_ticks()
    except Exception:
        pass

    title = ax.set_title(f"{title_prefix}\nt = {t0:.3g} y")

    def update(frame_idx: int):
        nonlocal current_vmin, current_vmax
        tg = time_groups[frame_idx]
        values, t = load_values(tg)
        w = values[mask]

        if scale in {"frame", "smooth"}:
            new_vmin, new_vmax = _auto_vrange(values, mask, robust=robust)
            if scale == "smooth":
                a = float(np.clip(scale_alpha, 0.0, 1.0))
                current_vmin = new_vmin if current_vmin is None else a * current_vmin + (1.0 - a) * new_vmin
                current_vmax = new_vmax if current_vmax is None else a * current_vmax + (1.0 - a) * new_vmax
            else:
                current_vmin, current_vmax = new_vmin, new_vmax
            try:
                mappable.set_clim(current_vmin, current_vmax)
                cb.update_normal(mappable)
                cb.update_ticks()
            except Exception:
                pass

        if method == "grid" and grid_state is not None:
            _, _, ij, W2 = grid_state
            W2[:] = np.nan
            W2[ij[:, 0], ij[:, 1]] = w
            mappable.set_array(W2.ravel())
        elif method == "tri" and tri is not None:
            nonlocal contour_artists
            for artist in contour_artists:
                try:
                    artist.remove()
                except Exception:
                    pass
            before = list(ax.collections)
            vmin_use = current_vmin if scale in {"frame", "smooth"} else vmin
            vmax_use = current_vmax if scale in {"frame", "smooth"} else vmax
            new_contour = ax.tricontourf(
                tri,
                w,
                levels=_levels_for(vmin_use, vmax_use),
                cmap=cmap,
                vmin=vmin_use,
                vmax=vmax_use,
            )
            after = list(ax.collections)
            contour_artists = [c for c in after if c not in before]
            try:
                cb.update_normal(new_contour)
                cb.update_ticks()
            except Exception:
                pass
        else:
            mappable.set_array(w)

        if scale in {"frame", "smooth"} and current_vmin is not None and current_vmax is not None:
            title.set_text(
                f"{title_prefix}\n"
                f"t = {t:.3g} y | range=[{current_vmin:.3g}, {current_vmax:.3g}]"
            )
        else:
            title.set_text(f"{title_prefix}\nt = {t:.3g} y")
        return (title,)

    ani = animation.FuncAnimation(fig, update, frames=len(time_groups), interval=1000 / max(fps, 1), blit=False)

    if gif:
        writer = animation.PillowWriter(fps=fps)
    else:
        if animation.writers.is_available("ffmpeg"):
            writer = animation.FFMpegWriter(fps=fps)
        else:
            writer = animation.PillowWriter(fps=fps)
            if not outpath.lower().endswith(".gif"):
                outpath = os.path.splitext(outpath)[0] + ".gif"

    ani.save(outpath, writer=writer, dpi=150)
    plt.close(fig)


def run_one_case(*, h5_path: str, outdir: str, args) -> None:
    import h5py

    os.makedirs(outdir, exist_ok=True)

    if bool(getattr(args, "clean", False)):
        import glob
        import shutil

        for pat in ("*.png", "*.gif", "*.mp4"):
            for p in glob.glob(os.path.join(outdir, pat)):
                try:
                    os.remove(p)
                except OSError:
                    pass
        # also remove any curated output subfolder
        showcase_dir = os.path.join(outdir, "showcase")
        if os.path.isdir(showcase_dir):
            shutil.rmtree(showcase_dir, ignore_errors=True)

    case_dir = os.path.dirname(os.path.abspath(h5_path))

    with h5py.File(h5_path, "r") as h5:
        time_groups = list_time_groups(h5)
        if not time_groups:
            raise RuntimeError(f"No PFLOTRAN time groups found in HDF5 output: {h5_path}")

        if args.every > 1:
            time_groups = time_groups[:: args.every]

        # Determine sample field shape so we can interpret Coordinates-based structured grids.
        sample_group = h5[time_groups[0].name]
        try:
            sample_var = next(iter(sample_group.keys()))
        except StopIteration:
            raise RuntimeError(f"No variables found in time group: {time_groups[0].name}")
        sample_shape = tuple(np.asarray(sample_group[sample_var]).shape)

        xc, yc, zc = read_cell_centers(h5, sample_shape)

        plane = args.plane
        if plane == "xy":
            z, ztol = args.z, args.ztol
            if z is None or ztol is None:
                z_auto, ztol_auto = estimate_default_slice(zc)
                z = z_auto if z is None else z
                ztol = ztol_auto if ztol is None else ztol

            z_unique = np.unique(np.asarray(zc, dtype=float))
            z_unique.sort()
            if z_unique.size > 0:
                z_min = float(z_unique[0])
                z_max = float(z_unique[-1])
                if float(z) < z_min - float(ztol) or float(z) > z_max + float(ztol):
                    z_near = float(z_unique[int(np.argmin(np.abs(z_unique - float(z))))])
                    print(
                        f"[warn] Requested --z={float(z):.3g} m is outside mesh Z range "
                        f"[{z_min:.3g}, {z_max:.3g}] m; using nearest Z={z_near:.3g} m."
                    )
                    z = z_near
            mask = select_cells_by_z(zc, z, ztol)
            u = xc
            v = yc
            xlabel, ylabel = "X [m]", "Y [m]"
            equal_aspect = True
            slice_desc = f"z≈{z:.2f}±{ztol:.2f} m"
            file_tag = f"xy_z{z:.2f}"
            plane_value = z
            plane_tol = ztol
        elif plane == "xz":
            y0, ytol = args.y, args.ytol
            if y0 is None or ytol is None:
                y_auto, ytol_auto = estimate_default_plane(yc)
                y0 = y_auto if y0 is None else y0
                ytol = ytol_auto if ytol is None else ytol
            mask = select_cells_by_y(yc, y0, ytol)
            u = xc
            v = zc
            xlabel, ylabel = "X [m]", "Z [m]"
            equal_aspect = False
            slice_desc = f"y≈{y0:.2f}±{ytol:.2f} m"
            file_tag = f"xz_y{y0:.2f}"
            plane_value = y0
            plane_tol = ytol
        else:  # yz
            x0, xtol = args.x, args.xtol
            if x0 is None or xtol is None:
                x_auto, xtol_auto = estimate_default_plane(xc)
                x0 = x_auto if x0 is None else x0
                xtol = xtol_auto if xtol is None else xtol
            mask = select_cells_by_x(xc, x0, xtol)
            u = yc
            v = zc
            xlabel, ylabel = "Y [m]", "Z [m]"
            equal_aspect = False
            slice_desc = f"x≈{x0:.2f}±{xtol:.2f} m"
            file_tag = f"yz_x{x0:.2f}"
            plane_value = x0
            plane_tol = xtol

        if args.method != "grid":
            rng = np.random.default_rng(args.seed)
            mask = maybe_downsample(mask, args.max_points, rng)

        wells_xy: list[tuple[str, float, float]] = []
        overlays_2d: list[tuple[str, np.ndarray, np.ndarray, dict]] = []

        def _add_well(name: str, well_filename: str):
            well_path = os.path.join(case_dir, well_filename)
            xy = read_well_xy_from_well_file(well_path)
            if xy is not None:
                wells_xy.append((name, xy[0], xy[1]))
            wpath = read_well_path_from_well_file(well_path)
            if wpath is None:
                return
            wx, wy, wz = wpath
            if plane == "xy":
                return
            if plane == "xz":
                if xy is None:
                    return
                if abs(xy[1] - float(plane_value)) <= float(plane_tol):
                    overlays_2d.append(
                        (
                            name,
                            np.asarray(wx, dtype=float),
                            np.asarray(wz, dtype=float),
                            {
                                "color": "white",
                                "lw": 2.0,
                                "marker": "o",
                                "ms": 3.0,
                                "mec": "black",
                                "mew": 0.6,
                            },
                        )
                    )
            else:  # yz
                if xy is None:
                    return
                if abs(xy[0] - float(plane_value)) <= float(plane_tol):
                    overlays_2d.append(
                        (
                            name,
                            np.asarray(wy, dtype=float),
                            np.asarray(wz, dtype=float),
                            {"color": "white", "lw": 2.0, "marker": "o", "ms": 3.0, "mec": "black", "mew": 0.6},
                        )
                    )

        _add_well("Injector", "injector.well")
        _add_well("Producer", "producer.well")

        available_vars = get_available_vars(h5)
        if args.all_vars:
            vars_to_plot = available_vars
        elif args.var:
            vars_to_plot = args.var
        elif bool(getattr(args, "important6", False)):
            preferred = [
                "Gas_Saturation",
                "Liquid_Saturation",
                "Gas_Pressure [Pa]",
                "Liquid_Pressure [Pa]",
                "Gas_Mole_Fraction_CO2",
                "Temperature [C]",
            ]
            picked = [v for v in preferred if v in available_vars]
            if len(picked) < 1:
                picked = available_vars[:6]
            vars_to_plot = picked
        else:
            vars_to_plot = ["Gas_Saturation"]

        def read1d(group_name: str, var: str) -> np.ndarray:
            return np.asarray(read_field(h5, group_name, var), dtype=float).ravel()

        # If the user didn't specify an XY z-slice but asked for the "important" set,
        # try to pick a Z level that actually changes for Gas_Saturation.
        if plane == "xy" and bool(getattr(args, "important6", False)) and args.z is None and "Gas_Saturation" in vars_to_plot:
            # Choose the Z level with the largest change in max saturation over time.
            z_unique = np.unique(np.asarray(zc, dtype=float))
            z_unique.sort()
            if z_unique.size > 0:
                max_by_z = np.full((z_unique.size, len(time_groups)), np.nan, dtype=float)
                masks_by_z = [select_cells_by_z(zc, float(zz), float(plane_tol)) for zz in z_unique]
                for it, tg in enumerate(time_groups):
                    vv = read1d(tg.name, "Gas_Saturation")
                    for iz, mz in enumerate(masks_by_z):
                        w = vv[mz]
                        w = w[np.isfinite(w)]
                        if w.size:
                            max_by_z[iz, it] = float(np.max(w))
                deltas = np.nanmax(max_by_z, axis=1) - np.nanmin(max_by_z, axis=1)
                if np.any(np.isfinite(deltas)):
                    iz_best = int(np.nanargmax(deltas))
                    z_best = float(z_unique[iz_best])
                    mask = masks_by_z[iz_best]
                    slice_desc = f"z≈{z_best:.2f}±{float(plane_tol):.2f} m (auto)"
                    file_tag = f"xy_z{z_best:.2f}"
                    plane_value = z_best
                    print(f"[important6] XY: auto-picked z={z_best:.3g} m (Δmax sat={float(deltas[iz_best]):.3g}).")

        # Curated outputs: a couple of vivid animations + a few summary plots.
        if getattr(args, "showcase", False):
            showcase_dir = os.path.join(outdir, "showcase")
            os.makedirs(showcase_dir, exist_ok=True)

            def _pick_var(preferred: list[str], contains: list[str] | None = None) -> Optional[str]:
                for p in preferred:
                    if p in available_vars:
                        return p
                if contains:
                    for vname in available_vars:
                        lv = vname.lower()
                        if all(c.lower() in lv for c in contains):
                            return vname
                return None

            sat_var = _pick_var(["Gas_Saturation"], contains=["gas", "saturation"])
            pres_var = _pick_var(["Gas_Pressure [Pa]", "Liquid_Pressure [Pa]"], contains=["pressure", "pa"])
            temp_var = _pick_var(["Temperature [C]"], contains=["temperature"])

            chosen = [v for v in (sat_var, pres_var, temp_var) if v is not None]
            if not chosen:
                raise RuntimeError("--showcase: could not identify any variables to plot.")

            times = np.array([tg.time_years for tg in time_groups], dtype=float)

            # For XY showcase, auto-pick a Z slice that actually changes for saturation
            # when the requested z is missing or out-of-range.
            if plane == "xy":
                z_unique = np.unique(np.asarray(zc, dtype=float))
                z_unique.sort()
                if z_unique.size > 0:
                    z_min = float(z_unique[0])
                    z_max = float(z_unique[-1])
                    z_req = args.z
                    z_req_ok = z_req is not None and (z_min - float(plane_tol) <= float(z_req) <= z_max + float(plane_tol))

                    if sat_var is not None and (z_req is None or not z_req_ok):
                        # Choose the Z level with the largest change in max saturation over time.
                        # (Fast + robust, and tends to produce visually informative animations.)
                        max_by_z = np.full((z_unique.size, len(time_groups)), np.nan, dtype=float)
                        masks_by_z = [select_cells_by_z(zc, float(zz), float(plane_tol)) for zz in z_unique]
                        for it, tg in enumerate(time_groups):
                            vv = read1d(tg.name, sat_var)
                            for iz, mz in enumerate(masks_by_z):
                                w = vv[mz]
                                w = w[np.isfinite(w)]
                                if w.size:
                                    max_by_z[iz, it] = float(np.max(w))

                        deltas = np.nanmax(max_by_z, axis=1) - np.nanmin(max_by_z, axis=1)
                        if np.any(np.isfinite(deltas)):
                            iz_best = int(np.nanargmax(deltas))
                            z_best = float(z_unique[iz_best])
                            mask = masks_by_z[iz_best]
                            slice_desc = f"z≈{z_best:.2f}±{float(plane_tol):.2f} m (auto)"
                            file_tag = f"xy_z{z_best:.2f}"
                            print(
                                f"[showcase] XY: auto-picked z={z_best:.3g} m for {sat_var} "
                                f"(max change Δ={float(deltas[iz_best]):.3g})."
                            )

            # 1) Saturation animation + stats
            if sat_var is not None:
                baseline_sat = read1d(time_groups[0].name, sat_var) if args.delta else None

                def sat_transform(arr: np.ndarray) -> np.ndarray:
                    return arr if baseline_sat is None else (arr - baseline_sat)

                def sat_loader(tg: TimeGroup):
                    return sat_transform(read1d(tg.name, sat_var)), tg.time_years

                safe_sat = sat_var.replace(" ", "_").replace("/", "_").replace("[", "").replace("]", "")
                out_anim = os.path.join(showcase_dir, f"anim_{safe_sat}_{file_tag}.gif")

                # Global colorbar across time for a clean animation.
                lo = np.inf
                hi = -np.inf
                for tg in time_groups:
                    vv = sat_transform(read1d(tg.name, sat_var))
                    v_lo, v_hi = _auto_vrange(
                        vv,
                        mask,
                        robust=True,
                        robust_percentiles=(float(args.global_pct[0]), float(args.global_pct[1])),
                    )
                    lo = min(lo, v_lo)
                    hi = max(hi, v_hi)

                animate_slice_2d(
                    u=u,
                    v=v,
                    mask=mask,
                    time_groups=time_groups,
                    load_values=sat_loader,
                    title_prefix=f"{sat_var} | {slice_desc}",
                    outpath=out_anim,
                    cmap="turbo",
                    vmin=lo,
                    vmax=hi,
                    s=args.marker_size,
                    fps=max(int(args.fps), 6),
                    gif=True,
                    method=args.method,
                    levels=args.levels,
                    overlays=overlays_2d if overlays_2d else None,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    cbar_label=sat_var if not args.delta else f"Δ {sat_var}",
                    scale="global",
                    robust=True,
                    scale_alpha=args.scale_alpha,
                    equal_aspect=equal_aspect,
                )

                # Saturation stats on the chosen slice.
                mean_v = np.empty(len(time_groups), dtype=float)
                p95_v = np.empty(len(time_groups), dtype=float)
                max_v = np.empty(len(time_groups), dtype=float)
                for i, tg in enumerate(time_groups):
                    vv = sat_transform(read1d(tg.name, sat_var))
                    w = vv[mask]
                    w = w[np.isfinite(w)]
                    if w.size == 0:
                        mean_v[i] = np.nan
                        p95_v[i] = np.nan
                        max_v[i] = np.nan
                    else:
                        mean_v[i] = float(np.mean(w))
                        p95_v[i] = float(np.percentile(w, 95))
                        max_v[i] = float(np.max(w))
                out_stats = os.path.join(showcase_dir, f"stats_{safe_sat}_{file_tag}.png")
                plot_timeseries_multi(
                    times=times,
                    series={"mean": mean_v, "p95": p95_v, "max": max_v},
                    title=f"{sat_var} on slice | {slice_desc}",
                    ylabel=sat_var if not args.delta else f"Δ {sat_var}",
                    outpath=out_stats,
                )

            # 2) Pressure delta animation (often visually interesting)
            if pres_var is not None:
                baseline_p = read1d(time_groups[0].name, pres_var)

                def pres_loader(tg: TimeGroup):
                    return read1d(tg.name, pres_var) - baseline_p, tg.time_years

                safe_p = pres_var.replace(" ", "_").replace("/", "_").replace("[", "").replace("]", "")
                out_anim_p = os.path.join(showcase_dir, f"anim_d{safe_p}_{file_tag}.gif")

                lo = np.inf
                hi = -np.inf
                for tg in time_groups:
                    dv = read1d(tg.name, pres_var) - baseline_p
                    v_lo, v_hi = _auto_vrange(
                        dv,
                        mask,
                        robust=True,
                        robust_percentiles=(float(args.global_pct[0]), float(args.global_pct[1])),
                    )
                    lo = min(lo, v_lo)
                    hi = max(hi, v_hi)
                # Make the diverging colormap symmetric when possible.
                if np.isfinite(lo) and np.isfinite(hi):
                    m = max(abs(float(lo)), abs(float(hi)))
                    lo, hi = -m, m

                animate_slice_2d(
                    u=u,
                    v=v,
                    mask=mask,
                    time_groups=time_groups,
                    load_values=pres_loader,
                    title_prefix=f"Δ {pres_var} (t - t0) | {slice_desc}",
                    outpath=out_anim_p,
                    cmap="RdBu_r",
                    vmin=lo,
                    vmax=hi,
                    s=args.marker_size,
                    fps=max(int(args.fps), 6),
                    gif=True,
                    method=args.method,
                    levels=args.levels,
                    overlays=overlays_2d if overlays_2d else None,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    cbar_label=f"Δ {pres_var}",
                    scale="global",
                    robust=True,
                    scale_alpha=args.scale_alpha,
                    equal_aspect=equal_aspect,
                )

            print(f"[showcase] Wrote curated plots to: {showcase_dir}")
            return

        for var in vars_to_plot:
            baseline = read1d(time_groups[0].name, var) if args.delta else None

            def transform(arr: np.ndarray) -> np.ndarray:
                if baseline is None:
                    return arr
                return arr - baseline

            cbar_label = var if not args.delta else f"Δ {var}"
            cmap = args.cmap or default_cmap(var, delta=args.delta)

            last = time_groups[-1]
            values_last = transform(read1d(last.name, var))

            vmin = args.vmin
            vmax = args.vmax
            if vmin is None or vmax is None:
                vmin0, vmax0 = _auto_vrange(values_last, mask, robust=args.robust)
                vmin = vmin0 if vmin is None else vmin
                vmax = vmax0 if vmax is None else vmax

            safe_var = var.replace(" ", "_").replace("/", "_").replace("[", "").replace("]", "")
            out_png = os.path.join(outdir, f"slice_{safe_var}_{file_tag}_t{last.index:03d}.png")
            if plane == "xy":
                plot_slice(
                    xc=xc,
                    yc=yc,
                    values=values_last,
                    mask=mask,
                    title=f"{var} | {slice_desc} | t={last.time_years:.3g} y",
                    outpath=out_png,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    s=args.marker_size,
                    method=args.method,
                    levels=args.levels,
                    wells_xy=wells_xy if wells_xy else None,
                    cbar_label=cbar_label,
                )
            else:
                plot_slice_2d(
                    u=u,
                    v=v,
                    values=values_last,
                    mask=mask,
                    title=f"{var} | {slice_desc} | t={last.time_years:.3g} y",
                    outpath=out_png,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    s=args.marker_size,
                    method=args.method,
                    levels=args.levels,
                    overlays=overlays_2d if overlays_2d else None,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    cbar_label=cbar_label,
                    equal_aspect=equal_aspect,
                )

            if args.timeseries:
                if args.point is None:
                    raise SystemExit("--timeseries requires --point X Y Z")
                point = (float(args.point[0]), float(args.point[1]), float(args.point[2]))
                idx = nearest_cell_index(xc, yc, zc, point)

                times = np.array([tg.time_years for tg in time_groups], dtype=float)
                series = np.empty(len(time_groups), dtype=float)
                for i, tg in enumerate(time_groups):
                    series[i] = float(transform(read1d(tg.name, var))[idx])

                out_ts = os.path.join(
                    outdir,
                    f"timeseries_{safe_var}_at_{point[0]:.1f}_{point[1]:.1f}_{point[2]:.1f}.png",
                )
                plot_timeseries(
                    times=times,
                    series=series,
                    title=f"{var} at nearest cell to ({point[0]:.1f}, {point[1]:.1f}, {point[2]:.1f}) m (cell #{idx})",
                    ylabel=cbar_label,
                    outpath=out_ts,
                )

            if args.animate:
                anim_vmin = vmin
                anim_vmax = vmax
                if (args.vmin is None or args.vmax is None) and args.scale == "global":
                    lo = np.inf
                    hi = -np.inf
                    for tg in time_groups:
                        vv = transform(read1d(tg.name, var))
                        # For a static (global) colorbar, exclude extreme outliers by default
                        # so a single bad cell doesn't destroy the contrast.
                        clip_extremes = (not bool(args.global_absolute)) or bool(args.global_robust)
                        if clip_extremes:
                            v_lo, v_hi = _auto_vrange(
                                vv,
                                mask,
                                robust=True,
                                robust_percentiles=(float(args.global_pct[0]), float(args.global_pct[1])),
                            )
                        else:
                            v_lo, v_hi = _auto_vrange(vv, mask, robust=False)
                        lo = min(lo, v_lo)
                        hi = max(hi, v_hi)
                    anim_vmin = lo if args.vmin is None else args.vmin
                    anim_vmax = hi if args.vmax is None else args.vmax
                elif (args.vmin is None or args.vmax is None) and args.scale in {"frame", "smooth"}:
                    vv0 = transform(read1d(time_groups[0].name, var))
                    lo0, hi0 = _auto_vrange(vv0, mask, robust=args.robust)
                    anim_vmin = lo0 if args.vmin is None else args.vmin
                    anim_vmax = hi0 if args.vmax is None else args.vmax

                def loader(tg: TimeGroup):
                    vv = transform(read1d(tg.name, var))
                    return vv, tg.time_years

                ext = ".gif" if args.gif else ".mp4"
                out_anim = os.path.join(outdir, f"anim_{safe_var}_{file_tag}{ext}")

                snap_n = int(getattr(args, "snapshot_count", 0) or 0)
                if snap_n > 0:
                    n_frames = len(time_groups)
                    if n_frames > 0:
                        if snap_n >= n_frames:
                            snap_idx = list(range(n_frames))
                        else:
                            snap_idx = np.linspace(0, n_frames - 1, snap_n, dtype=int).tolist()
                            # de-dup while preserving order
                            seen: set[int] = set()
                            snap_idx = [i for i in snap_idx if not (int(i) in seen or seen.add(int(i)))]

                        for j, fi in enumerate(snap_idx, start=1):
                            tg = time_groups[int(fi)]
                            vv = transform(read1d(tg.name, var))
                            out_snap = os.path.join(
                                outdir,
                                f"snap{j:02d}_{safe_var}_{file_tag}_t{tg.index:03d}.png",
                            )
                            if plane == "xy":
                                plot_slice(
                                    xc=xc,
                                    yc=yc,
                                    values=vv,
                                    mask=mask,
                                    title=f"{var} | {slice_desc} | t={tg.time_years:.3g} y",
                                    outpath=out_snap,
                                    cmap=cmap,
                                    vmin=anim_vmin,
                                    vmax=anim_vmax,
                                    s=args.marker_size,
                                    method=args.method,
                                    levels=args.levels,
                                    wells_xy=wells_xy if wells_xy else None,
                                    cbar_label=cbar_label,
                                )
                            else:
                                plot_slice_2d(
                                    u=u,
                                    v=v,
                                    values=vv,
                                    mask=mask,
                                    title=f"{var} | {slice_desc} | t={tg.time_years:.3g} y",
                                    outpath=out_snap,
                                    cmap=cmap,
                                    vmin=anim_vmin,
                                    vmax=anim_vmax,
                                    s=args.marker_size,
                                    method=args.method,
                                    levels=args.levels,
                                    overlays=overlays_2d if overlays_2d else None,
                                    xlabel=xlabel,
                                    ylabel=ylabel,
                                    cbar_label=cbar_label,
                                    equal_aspect=equal_aspect,
                                )
                if plane == "xy":
                    animate_slice(
                        xc=xc,
                        yc=yc,
                        mask=mask,
                        time_groups=time_groups,
                        load_values=loader,
                        title_prefix=f"{var} | {slice_desc}",
                        outpath=out_anim,
                        cmap=cmap,
                        vmin=anim_vmin,
                        vmax=anim_vmax,
                        s=args.marker_size,
                        fps=args.fps,
                        gif=args.gif,
                        method=args.method,
                        levels=args.levels,
                        wells_xy=wells_xy if wells_xy else None,
                        cbar_label=cbar_label,
                        scale=args.scale,
                        robust=args.robust,
                        scale_alpha=args.scale_alpha,
                    )
                else:
                    animate_slice_2d(
                        u=u,
                        v=v,
                        mask=mask,
                        time_groups=time_groups,
                        load_values=loader,
                        title_prefix=f"{var} | {slice_desc}",
                        outpath=out_anim,
                        cmap=cmap,
                        vmin=anim_vmin,
                        vmax=anim_vmax,
                        s=args.marker_size,
                        fps=args.fps,
                        gif=args.gif,
                        method=args.method,
                        levels=args.levels,
                        overlays=overlays_2d if overlays_2d else None,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        cbar_label=cbar_label,
                        scale=args.scale,
                        robust=args.robust,
                        scale_alpha=args.scale_alpha,
                        equal_aspect=equal_aspect,
                    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", default="pflotran.h5", help="PFLOTRAN HDF5 output (default: pflotran.h5)")
    parser.add_argument("--outdir", default="figures", help="Output directory for plots")

    parser.add_argument(
        "--all-scenarios",
        action="store_true",
        help="Find and process every pflotran.h5 under --search-root (ignores --h5)",
    )
    parser.add_argument(
        "--search-root",
        default=".",
        help="Root folder to search for scenarios when using --all-scenarios (default: .)",
    )

    parser.add_argument("--list-vars", action="store_true", help="Print available variables and exit")

    parser.add_argument(
        "--showcase",
        action="store_true",
        help=(
            "Generate a curated set of 'interesting' plots/animations (e.g., Gas_Saturation animation, "
            "pressure-delta animation, and saturation slice stats) into <outdir>/showcase."
        ),
    )

    parser.add_argument(
        "--var",
        action="append",
        default=None,
        help=(
            "Variable name to plot (must match HDF5 dataset name). "
            "You can repeat --var multiple times to generate multiple plots/animations in one run."
        ),
    )
    parser.add_argument(
        "--all-vars",
        action="store_true",
        help="Plot all variables found in the HDF5 (overrides --var)",
    )

    parser.add_argument(
        "--important6",
        action="store_true",
        help=(
            "Plot/animate a curated set of 6 key variables (Gas/Liquid Saturation, Gas/Liquid Pressure, "
            "Gas CO2 mole fraction, Temperature) when --var is not provided."
        ),
    )

    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete existing *.png/*.gif/*.mp4 in --outdir before writing new outputs.",
    )

    parser.add_argument(
        "--plane",
        choices=["xy", "xz", "yz"],
        default="xy",
        help="Slice plane: xy (const Z), xz (const Y), yz (const X)",
    )

    parser.add_argument("--z", type=float, default=None, help="Z slice center [m] for --plane xy (default: auto)")
    parser.add_argument("--ztol", type=float, default=None, help="Z slice half-thickness [m] for --plane xy (default: auto)")

    parser.add_argument("--y", type=float, default=None, help="Y plane [m] for --plane xz (default: auto)")
    parser.add_argument("--ytol", type=float, default=None, help="Y plane tolerance [m] for --plane xz (default: auto)")

    parser.add_argument("--x", type=float, default=None, help="X plane [m] for --plane yz (default: auto)")
    parser.add_argument("--xtol", type=float, default=None, help="X plane tolerance [m] for --plane yz (default: auto)")

    parser.add_argument("--every", type=int, default=1, help="Use every Nth time step (default: 1)")
    parser.add_argument("--max-points", type=int, default=40000, help="Max points to plot per frame (default: 40000)")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for downsampling")

    parser.add_argument(
        "--cmap",
        default=None,
        help="Matplotlib colormap (default chooses a sensible one based on variable name)",
    )
    parser.add_argument("--vmin", type=float, default=None, help="Colorbar min")
    parser.add_argument("--vmax", type=float, default=None, help="Colorbar max")
    parser.add_argument("--marker-size", type=float, default=4.0, help="Scatter marker size")

    parser.add_argument(
        "--method",
        choices=["grid", "tri", "scatter"],
        default="grid",
        help="Slice rendering method: grid=pcolormesh (best), tri=contourf on triangulation, scatter=dots",
    )
    parser.add_argument("--levels", type=int, default=30, help="Contour levels for tri method")
    parser.add_argument(
        "--delta",
        action="store_true",
        help="Plot change from initial time (value - value_at_t0) for slices/animations/timeseries",
    )
    parser.add_argument(
        "--scale",
        choices=["global", "frame", "smooth"],
        default="smooth",
        help="Color scaling mode for animations: global (fixed), frame (auto each frame), smooth (auto but stabilized)",
    )
    parser.add_argument(
        "--scale-alpha",
        type=float,
        default=0.75,
        help="Smoothing factor for --scale smooth (higher = less jitter)",
    )
    parser.add_argument(
        "--robust",
        action="store_true",
        help="Use robust (2-98%%) color limits when auto-scaling",
    )
    parser.add_argument(
        "--global-robust",
        action="store_true",
        help=(
            "(Deprecated) Kept for compatibility. For --scale global, prefer --global-pct or --global-absolute."
        ),
    )
    parser.add_argument(
        "--global-absolute",
        action="store_true",
        help="For --scale global, use true global min/max (no outlier clipping)",
    )
    parser.add_argument(
        "--global-pct",
        nargs=2,
        type=float,
        metavar=("PLOW", "PHIGH"),
        default=(2.0, 98.0),
        help="For --scale global, clip extremes using percentiles (default: 2 98)",
    )

    parser.add_argument("--animate", action="store_true", help="Create an animation over time")
    parser.add_argument("--fps", type=int, default=4, help="Animation frames per second")
    parser.add_argument("--gif", action="store_true", help="Force GIF output (default tries MP4 if ffmpeg is available)")

    parser.add_argument(
        "--snapshot-count",
        type=int,
        default=0,
        help=(
            "When used with --animate, also write N evenly-spaced PNG snapshots "
            "(snap01_*, snap02_*, ...) using the SAME color limits as the animation."
        ),
    )

    parser.add_argument(
        "--point",
        nargs=3,
        type=float,
        default=None,
        metavar=("X", "Y", "Z"),
        help="Point (x y z) [m] for nearest-cell time series",
    )
    parser.add_argument("--timeseries", action="store_true", help="Generate point time series (requires --point)")

    args = parser.parse_args()

    # List-vars mode only makes sense for a single file.
    if args.list_vars:
        import h5py

        with h5py.File(args.h5, "r") as h5:
            for v in get_available_vars(h5):
                print(v)
        return 0

    if args.all_scenarios:
        h5_paths = find_h5_files(args.search_root, filename=os.path.basename(args.h5))
        if not h5_paths:
            raise SystemExit(f"No {os.path.basename(args.h5)} found under: {args.search_root}")
        for p in h5_paths:
            scenario_dir = os.path.dirname(os.path.abspath(p))
            rel = os.path.relpath(scenario_dir, os.path.abspath(args.search_root))
            outdir = os.path.join(args.outdir, rel)
            print(f"[case] {scenario_dir} -> {outdir}")
            run_one_case(h5_path=p, outdir=outdir, args=args)
        print(f"Wrote plots to: {args.outdir}")
        return 0

    run_one_case(h5_path=args.h5, outdir=args.outdir, args=args)
    print(f"Wrote plots to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
