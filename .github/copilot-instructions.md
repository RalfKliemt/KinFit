# Project Guidelines

## Code Style
- Preserve the existing split between public headers in include/ and implementations in source/.
- Keep ROOT reflection/IO hooks intact for TObject-based classes:
  - Keep ClassDef(...) in headers.
  - Keep matching ClassImp(...) in source files.
  - Keep include/FitterLinkDef.h in sync when classes are added/removed.
- Follow the existing naming style: private data members use the f prefix (for example fChi2, fCands).
- Prefer small, local edits over broad refactors in matrix-heavy fitter logic.

## Architecture
- Core fit engine: KinFitter (iterative constrained kinematic fit, Jacobians, covariance updates).
- Physics object container: KFitParticle (track parameters + covariance matrix + ROOT serialization).
- Analysis orchestration: KFitAnalyzer (event loop and output writing), KFitDecayBuilder (candidate combinations and best-fit selection).
- Utility/reconstruction helpers: KFitDecayCandFinder, KFitVertexFinder, CoordinateConversion.
- Typical flow: build KFitParticle objects -> configure constraint in KinFitter/KFitAnalyzer -> run fit -> read fitted daughters/mother/probability.

## Build and Test
- Build/install (from repo root):
  - mkdir -p build
  - cd build
  - cmake .. -DCMAKE_INSTALL_PREFIX=<install-prefix>
  - make -j
  - make install
- There is no automated unit test suite in this repository.
- Validate behavior with the provided ROOT macros and examples:
  - analysis_user.C
  - QA/fit_toyMC_fromPluto.C
- For usage and integration notes, link to README.md rather than duplicating setup text.

## Conventions and Pitfalls
- Do not hand-edit generated ROOT dictionary files (for example KinFitDict.cxx in build outputs).
- In KinFitter, constraints are implemented as mutually exclusive mode flags; do not silently combine constraints unless explicitly redesigning that logic.
- Preserve covariance parameter ordering documented in KFitParticle: 1/p, theta, phi, R, Z.
- Changes to serialized ROOT classes can affect compatibility. Treat ClassDef version changes as deliberate API/IO changes, not incidental edits.
