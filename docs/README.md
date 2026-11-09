# KinFit documentation

KinFit is a C++/ROOT library for **kinematic fitting** of charged-particle
tracks in hadron-physics experiments. Measured track parameters are adjusted
within their uncertainties so that they fulfil physical constraints, such as
four-momentum conservation, an invariant mass or a common vertex. The fit
minimises a χ² using Lagrange multipliers and returns the improved tracks,
their covariances, χ², fit probability and pulls.

## Contents

| Page | What it covers |
| --- | --- |
| [installation.md](installation.md) | Building, installing, using KinFit from CMake or a ROOT session |
| [conventions.md](conventions.md) | Units, track parametrisation (1/p, θ, φ, R, Z), sign of R, covariance layout, magnetic field |
| [kinfitter.md](kinfitter.md) | The `KinFitter` class: every constraint, its inputs and outputs, convergence, pulls |
| [vertexing.md](vertexing.md) | Vertex finding with straight lines or helices in constant Bz; building a decay candidate |
| [analyzer.md](analyzer.md) | Automated event-by-event fitting of ROOT files with `KFitAnalyzer` |
| [examples/](examples/) | Runnable examples |

## Classes at a glance

| Class | Role |
| --- | --- |
| `KinFitParticle` | A track: `TLorentzVector` + (R, Z) + 5×5 covariance + bookkeeping (PID, charge, …) |
| `KinFitter` | The fitter. One constraint type per fitter instance |
| `KFitVertexFinder` | Least-squares common point of straight tracks |
| `KFitVertexFinderConstB` | Common point of helix tracks in a constant field Bz |
| `KFitDecayCandFinder` | Builds a decaying "mother" candidate (e.g. Λ) from its daughters and two vertices, as input for a 3C fit |
| `KFitDecayBuilder` | Loops over all particle combinations of an event and keeps the best fit |
| `KFitAnalyzer` | Reads a tree of `KinFitParticle`s, runs `KFitDecayBuilder` per event, writes the best fit |
| `CoordinateConversion` | Helper to convert between (p, θ, φ, R, Z) and Cartesian track parameters |

## Quick start

```cpp
#include "KinFitter.h"

// 1. Build particles from measured tracks (units: MeV/c, mm, rad)
TLorentzVector lvProton, lvPion;             // momentum + mass hypothesis
lvProton.SetXYZM(px1, py1, pz1, 938.272);
lvPion.SetXYZM(px2, py2, pz2, 139.570);

KinFitParticle proton(lvProton, R1, Z1);     // R, Z: closest approach to beam
KinFitParticle pion(lvPion, R2, Z2);
proton.setCovariance(covProton);              // 5x5 in (1/p, θ, φ, R, Z)
pion.setCovariance(covPion);

// 2. Fit
KinFitter fitter({proton, pion});
fitter.addMassConstraint(1115.683);           // M(p π-) = M(Λ)
if (fitter.fit()) {
  KinFitParticle pFit = fitter.getDaughter(0);
  double prob = fitter.getProb();
  double pullInvP = fitter.getPull(0);        // pull of 1/p of particle 0
}
```

[examples/toy_lambda_fit.cxx](examples/toy_lambda_fit.cxx) is a complete
version of this that generates its own input. It runs with no data file.

## Examples

| File | Shows | Input |
| --- | --- | --- |
| [toy_lambda_fit.cxx](examples/toy_lambda_fit.cxx) | Mass and vertex fit of generated Λ → pπ⁻, pull and probability checks | none |
| [vertex_constB.cxx](examples/vertex_constB.cxx) | Helix vertex finding in Bz = 1 T compared with straight lines | none |
| [fit_toyMC_fromPluto.C](examples/fit_toyMC_fromPluto.C) | Vertex fit on a Pluto toy MC tree, filling QA histograms | Pluto tree (not shipped) |
| [analysis_user.C](examples/analysis_user.C) | `KFitAnalyzer` driver macro | `KinFitParticle` tree (not shipped) |

Run the two self-contained ones as compiled programs:

```bash
cmake -S . -B build -DKINFIT_BUILD_EXAMPLES=ON
cmake --build build
./build/toy_lambda_fit 2000
./build/vertex_constB 10
```

or interpreted in ROOT:

```bash
cd docs/examples
root -l -e 'gSystem->Load("<build-or-install>/libKinFit")' 'toy_lambda_fit.cxx(2000)'
```

Expected output of `toy_lambda_fit`:

```text
Mass fit:   converged 2000/2000
  M(p pi-) before: mean  1115.75  rms   1.46 MeV/c^2
  M(p pi-) after : mean  1115.68  rms   0.00 MeV/c^2
  <prob> = 0.502 (expect 0.5 for a flat distribution)
  proton pulls (mean / width):  1/p -0.04/0.98  theta +0.01/0.98  phi -0.03/0.98

Vertex fit: converged 2000/2000
  <prob> = 0.489
  proton pulls (mean / width):  theta +0.01/1.00  phi -0.03/1.00  R -0.02/1.00  Z -0.01/1.00
```

Pulls centred at 0 with width 1 and a mean probability of 0.5 show that the
input covariances are consistent with the smearing.

## Related

- Python implementation: <https://github.com/KinFit/PyKinFit>
