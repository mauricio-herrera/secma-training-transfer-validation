# SECMA training-transfer validation: reproducibility package

This package supports the revised manuscript **“Training-transfer validation of a portable VR laparoscopic simulator: 6-DoF kinematic assessment on a common physical reference task.”**

## Scope
The observed-data analysis contains **18 participants (9 BOX, 9 SECMA)** and **66 unique PRE/POST assessment trials**. All empirical estimates use observed participants only. Posterior-predictive simulation is used solely to study prospective sample-size/equivalence operating characteristics. Synthetic draws are never counted as experimental participants.

## Public data
`data/subject_phase_metrics_v2.csv` contains the participant-phase outcomes used for the primary analysis. `data/trial_metrics_anonymized.csv` contains trial-level QC and kinematic outcomes with pseudonymous trial and participant identifiers. Original filenames and the internal identity key are excluded.

## Reproduction
Create a Python environment and run:

```bash
pip install -r requirements.txt
python code/reproduce_statistics.py
```

The script writes regenerated statistical tables to `reproduced_results/`.

`code/signal_processing_reference.py` provides the generic quaternion/timestamp processing functions used to define the kinematic outcomes. Raw participant telemetry is not included in the public package because source acquisition records may contain identifiers; release of raw telemetry should follow the applicable ethics/consent conditions.

## Analysis principles
- primary between-arm model: POST ~ PRE + training group with HC3 SE;
- repeated trials are not treated as independent participants; participant-phase medians define the subject-level dataset;
- equivalence is evaluated only for the designated primary ANCOVA;
- the 0.2 x SD_PRE margin is a standardized reference margin, not a clinical threshold;
- posterior-predictive cohorts are design simulations, not augmented experimental observations.

## Repository release
This directory is prepared for GitHub publication and Zenodo archival release. Before public release, authors should add the final manuscript citation, repository URL/DOI, ethics-approved data-sharing statement, and chosen data license.
