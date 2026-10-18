# Round 10b — The secondary fold has a verified saddle-node energy law

**18 October 2026. Model-only theoretical prediction, not experimental
validation or a dynamical simulation.**

## Resuming the interrupted round

The [full four-variable Round 10 audit](ROUND10_TRUE_BIFURCATION_SEQUENCE.md)
already established that, for reflection-symmetric channel profiles

\[
w_\pm(x)=1+0.05\cos(kx)
       \pm0.0004[\cos(2kx)-\cos(4kx)]
\quad (\mathrm{mm}),\qquad k=2\pi/(0.3\mathrm{\ mm}),
\]

the **positive-phase** channel at an assumed sidewall Young angle
36.5°, and depth 0.38 mm, develops a **stable asymmetric
equilibrium** when liquid volume falls below the symmetric critical
volume. A second saddle-node ends that stable branch.

The full free-contact-angle model uses
\(z=(l,r,\theta_l,\theta_r)\), fixed area \(A(z)\), and
surface-tension-normalized energy \(E(z)\); equilibrium has
\(\nabla E+\lambda\nabla A=0\).
Stability is assessed using the Hessian projected onto the
**three-dimensional fixed-volume tangent space**.

### Decisive independently reproduced coordinates

| Local event | Predicted volume | Centre displacement |
|---|---:|---:|
| Symmetric pitchfork at positive phase | 0.3159443695 µL | 0 µm |
| Secondary fold (asymmetric stable + unstable collide) | **0.3133490185 µL** | **23.600694 µm** |

The interval is
\[
\boxed{\Delta V=2.595351\,\mathrm{nL}},
\]
just **0.8215%** of the liquid volume at onset. In a
quasistatic evaporative experiment, the model therefore
predicts *continuous* symmetry breaking followed by a later
discontinuity if the evolving branch remains selected.

At the secondary fold, an independent numerical continuation
gave fixed-volume eigenvalues approximately
\[
(0,\;0.3530304,\;3.5563264).
\]
The zero mode is not inferred from a plotted turning point
alone; the other two modes remain positive.

## Independent check of the 7 pJ saddle

At the *initial* symmetric critical volume, a separate
full four-variable KKT solve finds an **index-one saddle**
at centre displacement
\[
c_s=\boxed{33.6298313\,\mu\mathrm{m}}.
\]

The stationary equation residual is below
\(7\times10^{-16}\), and the projected Hessian spectrum is
\[
(-0.281833,\;0.217396,\;4.173446).
\]

The normalized energy difference from the symmetric critical
state is
\[
\boxed{\Delta E/\gamma d=0.000255112434\ \mathrm{mm}}.
\]

For **illustrative**, not measured, surface tension
\(\gamma=0.072\,\mathrm{N/m}\) and depth
\(d=0.38\,\mathrm{mm}\), this gives a candidate local
barrier of **6.979876 pJ**, approximately
**1.6965 billion** \(k_BT\) at 298 K. This is not a
global minimax activation barrier until competing pathways
and all relevant contact-line modes have been excluded.
Thermally driven hops are not a plausible *sole* explanation
for transitions across this modeled barrier far from a
fold; evaporation and nonequilibrium forcing must be considered.

The shallow actual sidewall apparent angles at the saddle are
approximately **9.29° and 21.76°**, despite the
assumed intrinsic Young angle 36.5°. This strongly limits
direct extrapolation to experiments without independent wetting
and meniscus imaging.

## New result: near-fold energy-barrier scaling

Track the stable and index-one unstable asymmetric stationary
solutions at the **same total volume** by solving
equal meniscus pressure, contact Young angles and total area,
then project the full Hessian. As the second fold is approached
from above,

\[
\boxed{\Delta E(V)\;\propto\;
  (V-V_{\mathrm{fold}})^{3/2}}.
\]

A direct numerical fit at volume distances 0.005–0.15 nL
gives a log-log exponent **1.5000918**, consistent
with the standard saddle-node catastrophe normal form.
This is a confirmation of **classical** catastrophe theory,
not discovery of a new critical exponent.

Using the illustrative \(\gamma\) and depth above:

| Above second fold (nL) | Stable centre (µm) | Saddle centre (µm) | Barrier (pJ) |
|---:|---:|---:|---:|
| 0.02 | 22.5359 | 24.6216 | 0.004703 |
| 0.08 | 21.4202 | 25.6044 | 0.037627 |
| 0.40 | 18.3802 | 27.8988 | 0.420907 |
| 1.00 | 14.5199 | 30.1482 | 1.665387 |
| 2.00 | 8.2363 | 32.5366 | 4.717625 |

The near-fold prefactor at the illustrative physical scaling
is approximately
\[
\Delta E \simeq 1.663\ \mathrm{pJ}\,
\left[\frac{V-V_{\mathrm{fold}}}{1\,\mathrm{nL}}\right]^{3/2}.
\]

**Important physical distinction:** a high activation barrier
far from the fold is entirely compatible with a *vanishing*
barrier at the fold. A driven meniscus need not jump from the
distant 7 pJ saddle. Volume evolution can remove its local
metastable state, after which there is no activation barrier
left at that bifurcation.

### Reliability

The result was independently recomputed from:
- Full analytic gradients of the unconstrained four-variable
  surface and wall-wetting energy.
- A separate high-resolution cubic-spline antiderivative of
  the wall length, with agreement to the existing Round 10
  energy curves at approximately \(10^{-9}\) model-energy units.
- Full stationary KKT roots and three projected Hessian
  eigenvalues, not merely a one-dimensional energy graph.
- Two independent stationary solutions at each fixed volume
  above the second fold. One has strictly positive constrained
  eigenvalues, the other exactly one negative eigenvalue.
- Independent saddle-node minimization of the volume along
  the equal-pressure branch, yielding a vanishing lowest
  Hessian eigenvalue.

The repository's reduced regression test repeats the
pressure-balanced branch and barrier-exponent checks.
This is robust model verification, **not** experimental evidence.

## Why it matters for the original question

The strongest practical prediction is a **two-stage response**
from a geometrical Fourier-phase design:
1. A continuous symmetry-breaking branch emerges as liquid
   volume passes \(V_*\approx0.3159444\,\mu\mathrm{L}\).
2. After decreasing by another 2.595 nL, that local
   asymmetric equilibrium terminates at the fold.
3. The barrier toward that fold scales as the universal
   \(3/2\) law, giving a concrete signature for volume-resolved
   microscopy and for possible controlled forcing experiments.

The paired negative-phase geometry instead has
a **subcritical local onset** at
\(V_*\approx0.3110371\,\mu\mathrm{L}\).
The two channels have **equal roughness power spectra**,
but different relative Fourier phases.

Do **not** claim experimental jump lengths, rates or universal
outcomes. The true threshold could be pre-empted by a
nonneutral 3D mode, chemical pinning, corner films, or
another saddle and driving protocol. Global minimum-energy
paths and real surface topographies have not been validated.

## Remaining gates

1. Mesh-converged, non-neutral plate Young–Laplace 3D
   continuation using the *actual* 2/4-harmonic geometry,
   with a bordered quartic solver and verified competing modes.
2. Independent measured sidewall and plate Young angles,
   contact-angle hysteresis, 3D sidewall roughness and
   corner rounding.
3. Synchronized tracks of both contact lines and liquid
   volume: experimentally resolve a 2.6 nL window
   (0.8% of the plug volume), accounting for errors.
4. Blind prediction of event volumes and local trajectory
   types, using published Wells/Bisquert results as baseline.
5. Separate *global basin and minimum-energy-path* checks
   from local saddle classification; test whether observed
   dynamics follow the predicted branches.

**Professor assessment:** A meaningful, falsifiable
planar-model mechanism. Internal mathematical/numerical
confidence has increased; underlying fold/barrier exponents
are well-known and not independently novel.
Novelty vs existing wetting literature remains uncertain.
Experimental evidence: absent.
Journal-readiness remains **4/10**, potentially
*Langmuir*, *Physics of Fluids*, or *Physical Review Fluids*
after a strong 3D/external experimental test.

Run the existing core test, followed by the extension:

```bash
python cah-statics/round10_full_equilibria.py
python cah-statics/round10_fold_barrier_scaling.py
```

Outputs include `round10_secondary_fold_barrier.csv` and
`round10_secondary_fold_barrier.json`. The new calculation
is added to the CAH-scoped PR workflow only; no recurring
GitHub Actions schedule is introduced.
