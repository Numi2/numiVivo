# Cross-study temporal prediction: GSE226572 to Kang

**FAIL: 0/8 donors meet the predefined biological gate.** The unchanged native temporal interpolation model trained on 18 observed responses from GSE226572's three donors and predicted Kang's six-hour IFNB response for eight donors, using their controls. Each donor had to improve RMSE by at least 5% against both mean response and no change.

Candidate improvements were 0.39–1.88% versus mean response and 15.47–28.19% versus no change. All eight improved numerically over both baselines, but none met the required margin over mean response. These are reused development cohorts; this is not independent validation or a promoted product default.

The endpoint is whole-PBMC aggregate expression across 12,993 shared features. Counts are normalized to one million using the complete source feature denominator, then log1p transformed and projected. The dense Kang matrix has 16 aggregate rows, not individual cells. Composition and study differences remain part of this endpoint.

The native model is the unchanged executable from [TemporalInterpolation](../TemporalInterpolation/README.md). Six hours uses equal weights on the four- and eight-hour response knots. The mean baseline weights donors equally and each donor's observed times equally. Protocol and input identities are retained in `protocol.json`; all eight input/output hashes were frozen before scoring in `freeze.json`.

`run.py` and `run.log` preserve the first checker failure: its independent reference mistakenly used 0.75/0.25 interpolation weights. `score_frozen.py` corrects this to 0.5/0.5, verifies every frozen input/output hash, independently checks both native candidate and baseline arrays to 1e-12, and scores the existing outputs without refitting or regenerating predictions. `results.json` retains every donor, RMSE, gate and numerical discrepancy. Full inputs and outputs remain at `/Users/n/numivivo-duration-crossstudy-20260912`; source identities and binary hash are recorded, so rerunning requires those retained artifacts.

This result limits the preceding within-study success: adding duration alone does not establish useful transfer beyond a mean-response baseline. Further work must address biological context and validate on held-out studies without choosing models from their outcomes.
