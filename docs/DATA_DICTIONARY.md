# Data dictionary

## Public identifiers
- `participant`: pseudonymous participant code (`P001`--`P018`). No direct identifiers are included.
- `trial_id`: pseudonymous trial identifier. Original acquisition filenames are excluded.
- `group`: `BOX` or `SECMA`.
- `phase`: `PRE` or `POST`.
- `rep`: repeated trial number within participant and phase.

## Kinematic variables
- `active_time_s`: duration of valid contiguous task segments in seconds; gaps >50 ms are excluded.
- `total_path_m`: cumulative filtered 3-D instrument path length in metres.
- `p95_angspeed_rad_s`: 95th percentile of quaternion-derived angular-speed magnitude in rad/s.
- `mean_speed_m_s`, `p95_speed_m_s`: descriptive translational-speed metrics.

## Tracking-QC variables
- `n_samples`: number of valid samples in the trial.
- `wall_time_s`: first-to-last timestamp duration.
- `median_dt_s`: median valid sample interval.
- `effective_hz`: inverse median sample interval.
- `jitter_sd_ms`: SD of valid inter-sample intervals in ms.
- `p99_dt_ms`: 99th percentile of valid inter-sample intervals in ms.
- `n_gaps_gt_50ms`: count of gaps >50 ms.
- `max_gap_s`: maximum inter-sample gap.
- `q_sign_flips`: number of quaternion sign-continuity corrections.

## Privacy
Original filenames, participant names, email-like strings, and the internal participant key are deliberately excluded from the public/repository package.
