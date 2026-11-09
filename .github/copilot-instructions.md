# Project Guidelines

## Code Style
- Headers and sources live side by side: core/ holds the fitter (KinFitter, KinFitParticle,
  vertex finders, KFitDecayCandFinder, CoordinateConversion), analysis/ holds the
  event-loop layer (KFitAnalyzer, KFitDecayBuilder). Include headers by bare file name.
- Keep ROOT reflection/IO hooks intact for TObject-based classes:
  - Keep ClassDef(...) in headers.
  - Keep matching ClassImp(...) in source files.
  - Keep FitterLinkDef.h (repo root) and the class lists in CMakeLists.txt in sync when
    classes are added/removed.
- Follow the existing naming style: private data members use the f prefix (for example fChi2, fCands).
- Prefer small, local edits over broad refactors in matrix-heavy fitter logic.

## Architecture
- Core fit engine: KinFitter (iterative constrained kinematic fit, Jacobians, covariance updates).
- Physics object container: KinFitParticle (track parameters + covariance matrix + ROOT serialization).
- Analysis orchestration: KFitAnalyzer (event loop and output writing), KFitDecayBuilder (candidate combinations and best-fit selection).
- Utility/reconstruction helpers: KFitDecayCandFinder, KFitVertexFinder, CoordinateConversion.
- Typical flow: build KinFitParticle objects -> configure constraint in KinFitter/KFitAnalyzer -> run fit -> read fitted daughters/mother/probability.

## Build and Test
- Build/install (from repo root):
  - mkdir -p build
  - cd build
  - cmake .. -DCMAKE_INSTALL_PREFIX=<install-prefix>
  - make -j
  - make install
- There is no automated unit test suite yet (planned in QA/ with CTest, see TODO.md).
- Validate behavior with the examples in docs/examples/:
  - toy_lambda_fit.cxx and vertex_constB.cxx are self-contained; build them with
    -DKINFIT_BUILD_EXAMPLES=ON and run ./build/toy_lambda_fit, ./build/vertex_constB.
  - analysis_user.C and fit_toyMC_fromPluto.C need external input files.
- User documentation lives in docs/ (index: docs/README.md). Keep it in sync when
  changing public behaviour; link to it rather than duplicating text.

## Conventions and Pitfalls
- Do not hand-edit generated ROOT dictionary files (for example KinFitDict.cxx in build outputs).
- In KinFitter, constraints are implemented as mutually exclusive mode flags; do not silently combine constraints unless explicitly redesigning that logic.
- Preserve covariance parameter ordering documented in KinFitParticle: 1/p, theta, phi, R, Z.
- Changes to serialized ROOT classes can affect compatibility. Treat ClassDef version changes as deliberate API/IO changes, not incidental edits.
