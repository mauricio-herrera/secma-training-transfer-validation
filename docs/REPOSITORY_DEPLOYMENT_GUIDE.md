# GitHub + Zenodo deployment guide

## Recommended strategy
Publishing this de-identified package is strongly recommended for the revision because it makes the empirical N=18 analysis, the separation between observed and simulated data, and the exact statistical pipeline independently inspectable. The public repository should contain **derived anonymized data and code only**; do not upload the raw acquisition files or the internal participant key unless institutional ethics/consent explicitly permits it.

## GitHub
1. Create a new public repository, e.g. `secma-training-transfer-validation`.
2. Copy the contents of this folder to the repository root.
3. Replace `TO_BE_ADDED` in `CITATION.cff` with the GitHub repository URL.
4. Confirm the data-sharing statement and choose the data license.
5. Commit and push.
6. Create a tagged release `v1.0.0` corresponding exactly to the resubmission analysis. Do not change files after tagging; create a new release for later changes.

## Zenodo
1. Connect the GitHub repository to Zenodo.
2. Enable archival for the repository.
3. Publish the GitHub `v1.0.0` release; Zenodo will archive that release and mint a DOI.
4. Add the version DOI to the manuscript Data/Code Availability statement and the response to reviewers.
5. Update `CITATION.cff` and repository README with the DOI in a subsequent metadata-only release if needed.

## Recommended public contents
- anonymized participant-phase outcomes;
- anonymized trial-level QC/outcome table;
- primary ANCOVA, TOST, ICC, sensitivity and posterior-predictive result tables;
- fully reproducible statistical script;
- generic signal-processing reference implementation;
- data dictionary and environment requirements.

## Keep private / working only
- participant-name key;
- original filenames containing names/email strings;
- raw 6-DoF files unless ethics/consent explicitly permits public sharing;
- earlier synthetic development datasets;
- superseded analysis outputs.

## Metadata note
Zenodo supports both `CITATION.cff` and `.zenodo.json` for GitHub releases. If both are present, Zenodo uses `.zenodo.json` for the archived release metadata, so keep the two files synchronized before creating `v1.0.0`.
