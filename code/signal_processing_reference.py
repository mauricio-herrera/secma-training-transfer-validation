"""Reference 6-DoF signal-processing functions used in the SECMA analysis.

Input dataframe columns:
Timestamp, PosX, PosY, PosZ, RotX, RotY, RotZ, RotW

This public reference implementation intentionally contains no participant-name parsing
or source-file identifiers. It reproduces the kinematic definitions used in the study
when supplied with appropriately formatted de-identified raw telemetry.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.signal import butter, sosfiltfilt

GAP_THRESHOLD_S = 0.05
LOWPASS_CUTOFF_HZ = 10.0
FILTER_ORDER = 4


def enforce_quaternion_sign_continuity(q: np.ndarray):
    q = np.asarray(q, dtype=float).copy()
    norms = np.linalg.norm(q, axis=1)
    valid = norms > 0
    q[valid] /= norms[valid, None]
    flips = 0
    for i in range(1, len(q)):
        if np.dot(q[i - 1], q[i]) < 0:
            q[i] *= -1.0
            flips += 1
    return q, flips


def lowpass_uniform(x: np.ndarray, fs: float, cutoff: float = LOWPASS_CUTOFF_HZ,
                    order: int = FILTER_ORDER):
    if len(x) < 20 or fs <= 2 * cutoff:
        return x
    sos = butter(order, cutoff, btype="low", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=0)


def compute_trial_metrics(df: pd.DataFrame, gap_threshold_s: float = GAP_THRESHOLD_S,
                          cutoff_hz: float = LOWPASS_CUTOFF_HZ):
    required = ["Timestamp", "PosX", "PosY", "PosZ", "RotX", "RotY", "RotZ", "RotW"]
    x = df.copy()
    x["Timestamp"] = pd.to_datetime(x["Timestamp"], errors="coerce")
    for c in required[1:]:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    x = x.dropna(subset=required).sort_values("Timestamp")
    x = x.loc[~x["Timestamp"].duplicated(keep="first")]

    t = (x.Timestamp - x.Timestamp.iloc[0]).dt.total_seconds().to_numpy(float)
    p = x[["PosX", "PosY", "PosZ"]].to_numpy(float)
    q, flips = enforce_quaternion_sign_continuity(
        x[["RotX", "RotY", "RotZ", "RotW"]].to_numpy(float)
    )
    dt = np.diff(t)
    cuts = np.where((dt > gap_threshold_s) | (dt <= 0))[0] + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(t)]
    valid_dt = dt[(dt > 0) & (dt <= gap_threshold_s)]
    med_dt = float(np.median(valid_dt))
    fs = 1.0 / med_dt

    active_time = 0.0
    path_length = 0.0
    angular_speed = []
    linear_speed = []

    for a, b in zip(starts, ends):
        if b - a < 5:
            continue
        ts, ps, qs = t[a:b], p[a:b], q[a:b]
        local_dt = np.diff(ts)
        local_dt = local_dt[(local_dt > 0) & (local_dt <= gap_threshold_s)]
        if len(local_dt) == 0:
            continue
        h = float(np.median(local_dt))
        active_time += ts[-1] - ts[0]
        tu = np.arange(ts[0], ts[-1] + 0.5 * h, h)
        pu = np.column_stack([np.interp(tu, ts, ps[:, j]) for j in range(3)])
        pu = lowpass_uniform(pu, 1.0 / h, cutoff_hz)
        steps = np.linalg.norm(np.diff(pu, axis=0), axis=1)
        path_length += float(steps.sum())
        linear_speed.append(steps / h)

        r = Rotation.from_quat(qs)
        rotvec = (r[:-1].inv() * r[1:]).as_rotvec()
        dtn = np.diff(ts)
        good = (dtn > 0) & (dtn <= gap_threshold_s)
        if good.sum() >= 4:
            omega = rotvec[good] / dtn[good, None]
            tm = ((ts[:-1] + ts[1:]) / 2.0)[good]
            tg = np.arange(tm[0], tm[-1] + 0.5 * h, h)
            omega_u = np.column_stack([np.interp(tg, tm, omega[:, j]) for j in range(3)])
            omega_u = lowpass_uniform(omega_u, 1.0 / h, cutoff_hz)
            angular_speed.append(np.linalg.norm(omega_u, axis=1))

    speed = np.concatenate(linear_speed) if linear_speed else np.array([])
    ang = np.concatenate(angular_speed) if angular_speed else np.array([])
    return {
        "n_samples": len(x),
        "wall_time_s": float(t[-1] - t[0]),
        "active_time_s": active_time,
        "total_path_m": path_length,
        "mean_speed_m_s": float(np.mean(speed)),
        "p95_speed_m_s": float(np.percentile(speed, 95)),
        "p95_angspeed_rad_s": float(np.percentile(ang, 95)),
        "median_dt_s": med_dt,
        "effective_hz": fs,
        "jitter_sd_ms": float(np.std(valid_dt, ddof=1) * 1000),
        "p99_dt_ms": float(np.quantile(valid_dt, 0.99) * 1000),
        "n_gaps_gt_50ms": int((dt > gap_threshold_s).sum()),
        "max_gap_s": float(dt.max()),
        "q_sign_flips": flips,
    }
