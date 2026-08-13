#!/opt/miniconda3/bin/python3
"""Precompute specnm normal-mode solutions for the bodies the specnm.jl notebook
visualizes, at a fixed fmax=0.01 Hz, and write one JSON file per model to
src/assets/data/specnm_precomputed/<model>.json.

Run this once offline (it's slow -- mesh construction alone is 40s-2min per
family per model); the notebook itself never calls specnm, it only reads the
JSON this script produces. Re-run only when a model file under
src/assets/data/specnm_models/ changes or FMAX needs to change.

Usage:
    /opt/miniconda3/bin/python3 scripts/precompute_specnm.py [model ...]
    (defaults to prem_ani mars vpremoon europa if no models given)
"""
import itertools
import json
import os
import sys
import time

import numpy as np
import specnm

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODELS_DIR = os.path.join(REPO_ROOT, "src", "assets", "data", "specnm_models")
OUT_DIR = os.path.join(REPO_ROOT, "src", "assets", "data", "specnm_precomputed")
FMAX = 0.01
DEFAULT_MODELS = ["prem_ani", "mars", "vpremoon", "europa"]


def overtones(l_arr, l1_start):
    """Overtone number n for each mode. specnm returns modes already grouped
    into contiguous runs of matching angular order l (verified empirically),
    so a single groupby pass numbers each run 0,1,2,... -- except l==1, which
    starts at l1_start instead (spheroidal: 2, skipping translation/rotation;
    toroidal: 1, skipping the nonexistent 0T1)."""
    out = np.empty(len(l_arr), dtype=int)
    i = 0
    for l, group in itertools.groupby(l_arr):
        start = l1_start if l == 1 else 0
        for k in range(len(list(group))):
            out[i] = start + k
            i += 1
    return out


def component_arrays(mode, sph_type, ef):
    """Split specnm's interleaved-by-component eigenfunction matrix (nmodes,
    ncols) into one (nmodes, nradial) array per physical component."""
    if mode == "spheroidal":
        if sph_type == 3:
            return ["U", "V", "P"], [ef[:, 0::3], ef[:, 1::3], ef[:, 2::3]]
        return ["U", "V"], [ef[:, 0::2], ef[:, 1::2]]
    return ["W"], [ef]


def round_list(v, sig=6):
    """Round to significant figures, not decimal places. Eigenfunctions here come
    back from specnm normalized to a physical-units scale as small as ~1e-10 --
    round(x, digits=5) (fixed decimal places) silently zeroed every single one of
    them out, which is a data bug, not a plotting one (confirmed by comparing the
    precomputed JSON against a fresh live specnm call: identical near-zero raw
    values, so the widget was faithfully plotting already-destroyed data)."""
    out = []
    for x in v:
        x = float(x)
        out.append(float(f"{x:.{sig}g}") if np.isfinite(x) and x != 0 else 0.0)
    return out


def solve_family(solver_fn, problem_method, model_fname, attenuation_mode, l1_start):
    obj = solver_fn(model_fname, fmax=FMAX)
    out = getattr(obj, problem_method)(attenuation_mode=attenuation_mode, fmax=FMAX / 2)
    l = np.array(out["angular orders"], dtype=int)
    f = np.array(out["frequencies"], dtype=float) * 1000.0  # mHz
    n = overtones(l, l1_start)
    ef = np.array(out["eigenfunctions"], dtype=float)
    mode = str(obj.mode)
    sph_type = int(obj.sph_type) if mode == "spheroidal" else None
    radius = float(obj.radius) / 1000.0
    r = np.array(obj.r, dtype=float) / 1000.0
    depth = radius - r
    labels, comps = component_arrays(mode, sph_type, ef)
    # mode-major flat array per component: ef[c][mode*nradial + i] -- numpy's
    # default C (row-major) flatten of a (nmodes, nradial) array already is this.
    ef_flat = [round_list(c.flatten()) for c in comps]
    payload = {
        "l": [int(x) for x in l],
        "f": round_list(f),
        "n": [int(x) for x in n],
        "labels": labels,
        "depth": round_list(depth),
        "nradial": len(depth),
        "ef": ef_flat,
    }
    return payload, obj


def precompute_one(model):
    t0 = time.time()
    model_fname = os.path.join(MODELS_DIR, model)
    if not os.path.isfile(model_fname):
        print(f"[{model}] SKIP: no such model file at {model_fname}", flush=True)
        return

    print(f"[{model}] rayleigh: meshing + solving (fmax={FMAX})...", flush=True)
    sph_payload, ray = solve_family(specnm.rayleigh, "rayleigh_problem", model_fname, "elastic", 2)
    print(f"[{model}] rayleigh done in {time.time()-t0:.1f}s ({len(sph_payload['l'])} modes)", flush=True)

    t1 = time.time()
    print(f"[{model}] love: meshing + solving...", flush=True)
    # attenuation_mode="full" was the notebook's original live-solve setting for
    # Love (asymmetric with Rayleigh's "elastic", already flagged there as
    # unconfirmed-intentional) -- empirically it does not converge: on europa it
    # ran 2000+ shift-invert retries ("warning: not enough eigenvalues found")
    # without ever finishing. "elastic" solves cleanly in well under a second.
    # Using "elastic" for both families here, not just because it's faster --
    # because "full" appears to hang outright.
    tor_payload, lov = solve_family(specnm.love, "love_problem", model_fname, "elastic", 1)
    print(f"[{model}] love done in {time.time()-t1:.1f}s ({len(tor_payload['l'])} modes)", flush=True)

    earth_model = {
        "depth": round_list(float(ray.radius) / 1000.0 - np.array(ray.r) / 1000.0),
        "rho": round_list(np.array(ray.rho)),
        "vp": round_list(np.array(ray.vp)),
        "vs": round_list(np.array(ray.vs)),
    }
    result = {
        "spheroidal": sph_payload,
        "toroidal": tor_payload,
        "earth_model": earth_model,
        "model_name": model,
        "fmax": FMAX,
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, f"{model}.json")
    with open(out_path, "w") as fh:
        json.dump(result, fh, separators=(",", ":"))
    size_kb = os.path.getsize(out_path) / 1024
    print(f"[{model}] wrote {out_path} ({size_kb:.0f} KB) total {time.time()-t0:.1f}s", flush=True)


def main():
    models = sys.argv[1:] or DEFAULT_MODELS
    print(f"precomputing specnm output for: {models} (fmax={FMAX})", flush=True)
    failures = []
    for m in models:
        try:
            precompute_one(m)
        except Exception as e:
            print(f"[{m}] FAILED: {type(e).__name__}: {e}", flush=True)
            failures.append(m)
    if failures:
        print(f"DONE with failures: {failures}", flush=True)
        sys.exit(1)
    print("DONE, all models succeeded", flush=True)


if __name__ == "__main__":
    main()
