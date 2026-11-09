# KinFit review & extension plan

Working list for the full-project review. Tick items off as they land; each
phase should end in its own commit(s).

## 0. Housekeeping

- [x] Commit the pending `KFitVertexFinderConstB` work (staged files plus the
      `fit()`/verbosity/projection fixes in the working tree)
- [x] Remove the stale `build/` cache (ROOT 6.40.00 → 6.40.04 breaks configure)
- [x] Unify file headers (tabs/whitespace in `analysis_user.C`,
      `QA/fit_toyMC_fromPluto.C`)

## 1. Documentation (`docs/`)

- [x] `docs/README.md` – index / overview of the package
- [x] `docs/installation.md` – build, install, `find_package`, ROOT versions
- [x] `docs/conventions.md` – units (MeV/c vs GeV/c, mm), track parametrisation
      (1/p, θ, φ, R, Z), sign of R, covariance ordering, coordinate system, Bz sign
- [x] `docs/kinfitter.md` – all constraints (3C, 4C, Vertex, MissingParticle,
      Mass, MissingMass, MassVtx): what they constrain, required inputs, NDF,
      outputs (chi², prob, pulls, fitted daughters/mother)
- [x] `docs/vertexing.md` – `KFitVertexFinder` (straight lines) vs
      `KFitVertexFinderConstB` (helix, constant Bz), `KFitDecayCandFinder`
- [x] `docs/analyzer.md` – `KFitAnalyzer` / `KFitDecayBuilder` workflow,
      input/output tree layout
- [x] Move examples to `docs/examples/`:
  - [x] `analysis_user.C`
  - [x] `QA/fit_toyMC_fromPluto.C`
  - [x] new minimal standalone example (compiled, no input file needed)
        → `toy_lambda_fit.cxx`, `vertex_constB.cxx`, CMake option
        `KINFIT_BUILD_EXAMPLES`
- [x] Fix README: dangling "see" link, missing `QA/fit_toyMC.C` and
      `QA/test_files/` references, point to `docs/`
- [x] Update `.github/copilot-instructions.md` (test suite, example locations)

## 2. Tests (`QA/`, CTest)

- [ ] CMake: `option(KINFIT_BUILD_TESTS)`, `enable_testing()`, `QA/CMakeLists.txt`
- [ ] Small self-contained test helper (no external framework, or decide on
      GoogleTest/Catch2) + toy generator for tracks with known truth
- [ ] `KinFitParticle`: construction, getters/setters, R/Z ↔ XYZ conversion
- [ ] `CoordinateConversion`: Cartesian ↔ spherical round-trip, error propagation
      vs numeric Jacobian
- [ ] `KFitVertexFinder`: exact intersection, parallel tracks
- [ ] `KFitVertexFinderConstB`: helix toy (from scratchpad), loopers, B=0
      equivalence, `fit()` after setters, singular case
- [ ] `KinFitter` per constraint: unsmeared input → chi² ≈ 0 and constraint
      satisfied; smeared input → pull mean ≈ 0, σ ≈ 1, prob flat
- [ ] Analytic Jacobians vs finite differences for every constraint
- [ ] `KFitDecayCandFinder`: decay candidate & covariance vs toy truth
- [ ] I/O round-trip: write/read `KinFitParticle` in `TClonesArray`
- [ ] `KFitAnalyzer` end-to-end on a small generated input file
- [ ] Wire into CI (GitHub Actions with ROOT container)

## 3. Math review

Found while writing the docs (toy MC, see docs/kinfitter.md):

- [ ] **MissingMass broken**: `f_eval` uses P_init **+** Σpᵢ, `Feta_eval`
      uses P_init − Σpᵢ → 0/500 converged
- [ ] **MassVtx** converges in ~0.2 % of toy events
- [ ] **Convergence test divides by measured values**: any parameter exactly 0
      (R = 0 / Z = 0 from MC) → NaN → `fit()` never returns true
- [ ] MissingParticle with one measured track: D V Dᵀ singular (rank ≤ 3) —
      document or reformulate
- [ ] Multiple mass constraints: `addMassConstraint` clears
      `fFlexiParticlesInFit` before `f_eval` checks it → `setUseCandInFit`
      path unreachable; NDF incremented per call
- [ ] Combining constraint types silently overwrites `d` (exclusive flags)
- [ ] `fit()` calls `abort()` if no constraint is set
- [ ] 3C: daughters use only diagonal covariance; `fMother.setCovariance()`
      before `fMother = mother` is dead code; not yet toy-verified
- [ ] Fitted covariance returned only as diagonal
- [ ] 4C ⟨prob⟩ = 0.40 in a 500-event toy (expect 0.50) — investigate
- [ ] `KFitDecayCandFinder`: flight line anchored at (decay − primary) instead
      of decay vertex; sign of R from primary y; covariance uses √ of
      correlations and inconsistent signs; dθ/dx uses z/r² instead of z²/r²
- [ ] `KinFitParticle::setPhiRad/Deg` call `cand.SetTheta`;
      `setTOF/setRapidity` take `int`; `fSector` is `double`
- [ ] `CoordinateConversion`: defaults `9 * 10 ^ (-999)` (XOR, not power);
      `theta = atan(pT/pz)` wrong for pz < 0; `using namespace std` in header;
      not in dictionary

- [ ] `KinFitter` core: Lagrange-multiplier iteration, covariance update,
      chi²/NDF bookkeeping, convergence criteria
- [ ] Each constraint function f and its Jacobians (Fη, Fξ) — derivation in
      `docs/`, check against finite differences
- [ ] Unknown parameters (ξ) handling: vertex, missing particle
- [ ] Pull definition (and the `-10000` sentinel in `KinFitter.cxx`)
- [ ] Unit consistency: MeV vs GeV (`analysis_user.C` uses GeV values)
- [ ] `CoordinateConversion` error propagation formulas
- [ ] `KFitDecayCandFinder` covariance of the reconstructed decay candidate
- [ ] `KFitVertexFinder*`: unweighted LSQ → covariance-weighted option?
- [ ] Numerical robustness: explicit inverses vs decompositions, near-singular
      matrices, angle wrapping of φ

## 4. I/O review

Found while writing the docs:

- [ ] `KFitAnalyzer` event loop starts at 1 → first event never analysed
- [ ] `KFitDecayBuilder::fillFitCands()` writes `particleCounter[1e6]`
      after the last combination (out of bounds)
- [ ] Unknown task in `KFitDecayBuilder::doFit` → `abort()` via `fit()`
- [ ] Split-mode warning "custom Streamer" when writing `KinFitParticle`
      (inherits `TLorentzVector`)

- [ ] `KFitAnalyzer`: unchecked `TFile`/`TTree` opens, leaked input file,
      ownership of `TClonesArray`, behaviour at `nEvents > entries`
- [ ] Input file list / `TChain` support (constructor says "file(list)")
- [ ] `KinFitParticle` streamer: `ClassDef` version, transient members,
      schema evolution from v1
- [ ] `KinFitter` has `ClassDef` — does it need to be persistent?
- [ ] Missing `ClassDef`/dictionary entries vs `FitterLinkDef.h`
- [ ] Output tree layout documented and stable

## 5. Extensions (to be decided)

- [ ] Covariance-weighted / χ²-based vertex fit in constant B, returning a
      vertex covariance
- [ ] Helix (constant-B) track model inside `KinFitter` vertex constraints
- [ ] Combined constraints (e.g. 4C + mass) instead of exclusive modes
- [ ] Python bindings / sync with PyKinFit
- [ ] Combinatorics & Lists
- [ ] Parametrized B field
- [ ] …
