# KinFit: A Kinematic Fitting Package for Hadron Physics Experiments

KinFit is a C++/ROOT library for kinematic fitting of charged-particle tracks.
It adjusts measured tracks within their uncertainties to fulfil physical
constraints and returns the improved tracks, χ², fit probability and pulls.
Constraints include four-momentum conservation, invariant mass, a missing
particle and a common vertex.

**Full documentation: [docs/](docs/README.md)**

## Installation

```bash
git clone https://github.com/KinFit/KinFit.git
cd KinFit
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/your/install/prefix
cmake --build build -j
cmake --install build
```

KinFit needs ROOT 6 (Core, Matrix, Physics, MathCore, RIO, Tree). Use it from
CMake with

```cmake
find_package(KinFit REQUIRED)
target_link_libraries(your_executable PRIVATE KinFit::KinFit)
```

or in ROOT with `gSystem->Load("/your/install/prefix/lib/libKinFit")`.
See [docs/installation.md](docs/installation.md) for details.

## Usage in a nutshell

```cpp
KinFitParticle proton(lvProton, R1, Z1);   // TLorentzVector + R, Z of the track
proton.setCovariance(cov1);                // 5x5 in (1/p, theta, phi, R, Z)
KinFitParticle pion(lvPion, R2, Z2);
pion.setCovariance(cov2);

KinFitter fitter({proton, pion});
fitter.addMassConstraint(1115.683);        // MeV/c^2
if (fitter.fit()) {
  KinFitParticle fittedProton = fitter.getDaughter(0);
  double prob = fitter.getProb();
}
```

- Units, track parametrisation and covariance layout: [docs/conventions.md](docs/conventions.md)
- All constraints and fit results: [docs/kinfitter.md](docs/kinfitter.md)
- Vertex finding (straight lines and constant Bz) and decay candidates:
  [docs/vertexing.md](docs/vertexing.md)
- Automated event-by-event fitting of ROOT files with `KFitAnalyzer`:
  [docs/analyzer.md](docs/analyzer.md)
- Runnable examples: [docs/examples/](docs/examples/). Build them with
  `-DKINFIT_BUILD_EXAMPLES=ON`.

## Python implementation

[PyKinFit](https://github.com/KinFit/PyKinFit)
