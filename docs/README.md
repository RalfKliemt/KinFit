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
| [vertexing.md](vertexing.md) | Vertex finding with straight lines, helices in constant Bz or any field; building a decay candidate |
| [analyzer.md](analyzer.md) | Automated event-by-event fitting of ROOT files with `KinFitAnalyzer` |
| [examples/](examples/) | Runnable examples |
| [../QA/README.md](../QA/README.md) | The automated test suite |

## Classes at a glance

| Class | Role |
| --- | --- |
| `KinFitParticle` | A track: `TLorentzVector` + (R, Z) + 5×5 covariance + bookkeeping (PID, charge, …) |
| `KinFitter` | The fitter: one kinematic constraint (4C, 3C, missing mass/particle) plus mass and common-vertex constraints |
| `KinFitField`, `KinFitFieldUniform`, `KinFitFieldByZ`, `KinFitFieldByXZ`, `KinFitFieldFunction` | Magnetic field models: B = 0, constant Bz/By, By(z) / By(x, z) histograms, user function |
| `KinFitTrajectory` | Track propagation in a field (exact helix if uniform, Runge-Kutta otherwise) |
| `KinFitVertexFinder` | Least-squares common point of straight tracks |
| `KinFitVertexFinderConstB` | Common point of helix tracks in a constant field Bz |
| `KinFitVertexFinderField` | Common point of tracks in any field model |
| `KinFitDecayCandFinder` | Builds a decaying "mother" candidate (e.g. Λ) from its daughters and two vertices, as input for a 3C fit |
| `KinFitDecayBuilder` | Loops over all particle combinations of an event and keeps the best fit |
| `KinFitAnalyzer` | Reads a tree of `KinFitParticle`s, runs `KinFitDecayBuilder` per event, writes the best fit |
| `KinFitCoordinateConversion` | Helper to convert between (p, θ, φ, R, Z) and Cartesian track parameters |

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
| [toy_lambda_fit.cxx](examples/toy_lambda_fit.cxx) | Mass and common-vertex fit of generated Λ → pπ⁻, pulls, probability, fitted vertex | none |
| [vertex_constB.cxx](examples/vertex_constB.cxx) | Vertex finding in Bz = 1 T: helix finder, general-field finder and straight lines compared | none |
| [vertex_field_fit.cxx](examples/vertex_field_fit.cxx) | Common-vertex + 4C fit of p π⁻ γ in a By(z) dipole field, photon measured by energy and hit | none |
| [fit_toyMC_fromPluto.C](examples/fit_toyMC_fromPluto.C) | Vertex fit on a Pluto toy MC tree, filling QA histograms | Pluto tree (not shipped) |
| [analysis_user.C](examples/analysis_user.C) | `KinFitAnalyzer` driver macro | `KinFitParticle` tree (not shipped) |

Run the three self-contained ones as compiled programs:

```bash
cmake -S . -B build -DKINFIT_BUILD_EXAMPLES=ON
cmake --build build
./build/toy_lambda_fit 2000
./build/vertex_constB 10
./build/vertex_field_fit 500
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

Vertex fit: converged 2000/2000 (NDF = 2 tracks x 2 - 3 = 1)
  <prob> = 0.489
  proton pulls (mean / width):  theta +0.01/1.00  phi -0.03/1.00  R -0.02/1.00  Z -0.01/1.00
  |fitted vertex - truth|: mean 1.80 mm
  vertex pulls (mean / width):  x +0.04/0.97  y -0.01/1.01  z -0.02/1.00
```

Expected output of `vertex_field_fit`:

```text
converged 500/500, NDF = 5 (4 + 2 tracks x 2 - 3)
<prob> = 0.488 (expect 0.5 for a flat distribution)
M(sum of fitted daughters) - M(mother): mean 8.1e-08, rms 4.5e-07 MeV/c^2

vertex resolution [mm, 68% of |estimate - truth|]
        start (finder)      fit     fit pull mean/width
  x           0.365          0.347      -0.07/0.99
  y           0.372          0.359      +0.03/1.00
  z           1.072          1.054      +0.02/0.98

pulls (mean/width): photon E +0.05/1.08, proton theta -0.03/1.01
```

Pulls centred at 0 with width 1 and a mean probability of 0.5 show that the
input covariances are consistent with the smearing.

## Related

- Python implementation: <https://github.com/KinFit/PyKinFit>
