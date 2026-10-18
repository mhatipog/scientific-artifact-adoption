# A pressure–volume transfer identity for two identical capillary fronts

**Exact theorem inside an additive, reflection-symmetric equilibrium model; not a claim of a new physical mechanism or continuum 3D stability.**

Let H(x) be a smooth strictly increasing single-front volume coordinate, q = H(x), h = H'(s) > 0. Let the single-front pressure law be P(q)=p(H^{-1}(q)) and let its reduced energy obey E'(q)=-P(q), under a consistent liquid-minus-gas pressure sign convention.

For two *independent identical* menisci, at fixed total front-associated volume 2q0, define q1=q0+delta and q2=q0-delta. Additivity gives

    F(delta)=E(q0+delta)+E(q0-delta).
    F(delta)-F(0) = -P'(q0) delta^2
                     -P'''(q0) delta^4/12
                     -P'''''(q0) delta^6/360 + O(delta^8).

Only **odd derivatives of the single-front pressure–volume law** enter the even constrained energy because the two arms have exchanged volumes symmetrically. Thus the first antisymmetric quadratic loss of stability occurs at P'(q0)=0, i.e. p'(s)=0 when H'(s)>0. At this stationary pressure:

- P'''(q0)<0 produces a *positive* local quartic coefficient (locally supercritical when higher constraints permit);
- P'''(q0)>0 produces a *negative* local quartic coefficient (subcritical-capable; global jump still unproved);
- P'''(q0)=0 is a possible higher-order degeneracy. If also P'(q0)=0 and P'''''(q0) is nonzero, the leading local even term is sixth order with coefficient -P'''''/360 in delta.

For the original center shift c and adjusted half-length s(c) under fixed total volume, delta=H'(s)c+O(c^3). Applying the inverse-function chain rule at p'(s)=0 yields

    P'''(q0) = p'''(s)/[H'(s)]^3
                - 3 p''(s) H''(s)/[H'(s)]^4,

and the exact quartic coefficient in c is

    B4 = -(H')^4 P'''/12
       = [3 H'' p'' - H' p''']/12
       = (H'' p'' - H' p'''/3)/4.

This reproduces the pre-existing two-meniscus normal-form coefficient and shows it is a **coordinate-invariant, third-derivative criterion** when expressed in individual-front liquid volumes. Geometry enters by shaping P(q), not through an independent universal wall-slope metric. Numerically, the accompanying test integrates the analytic polynomial energy exactly, enforces the full volume constraint via Newton iteration, verifies convergence of F(c)/c^4 for both quartic signs, and checks the sixth-order degeneracy.

A possible practical experiment: measure the *single-front* P(q) curve with fully characterized contact angles/geometry. Its stationary-point derivatives would then predict the **local** two-front antisymmetric onset class **without fitting two-front motion**. This is a falsifiable transfer only if the two interfaces really have independent additive energies, exchange volumes at fixed total volume, and their other constrained modes remain stable.

**Objections and limits:**

1. Additivity need not survive interacting 3D interfaces, film coupling, trapped gas, dynamic pressure losses, surfactants, pinning heterogeneity or shared contact-line constraints.
2. Pressure P(q) third derivatives are noise-sensitive, especially where the pore fillet is only 2.25 camera pixels wide. Regularization would need preregistration and derivative-error propagation.
3. This is a standard-looking Landau symmetry reduction; claims of mathematical novelty require a dedicated prior-art audit. Wells et al. (2018) and Bisquert et al. (2024) already show fold/pitchfork competition in capillarity.
4. The quartic sign is local and does not by itself establish global minimum, continuum Hessian coercivity or an observable finite jump.

References: [Hatipogullari et al. (2019)](https://doi.org/10.1103/PhysRevFluids.4.044008), [Wells et al. (2018)](https://doi.org/10.1038/s41467-018-03840-6), [Bisquert et al. (2024)](https://doi.org/10.1021/acs.langmuir.4c02908).
