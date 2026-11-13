/** @file gsMaterialMatrixTFTSplit.hpp

    @brief Tension-field material decorator based on the Zhang-Kiendl split (linearised or exact wrinkled stress)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

#pragma once

#include <limits>

#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsCore/gsFunction.h>

namespace gismo
{

template <short_t dim, class T, bool linear >
gsMaterialMatrixTFTSplit<dim,T,linear>::gsMaterialMatrixTFTSplit(const material_ptr & materialMatrix)
:
m_materialMat(materialMatrix),
m_nThetaFallbacks(0),
m_nTFTFailures(0),
m_nExactNonconverged(0),
m_nExactSaddle(0),
m_nExactWrongBranch(0),
m_nExactStagnated(0),
m_nExactReclassifiedTaut(0),
m_nExactReclassifiedSlack(0),
m_nExactSkipped(0),
m_nExactUsed2x2(0),
m_nExactMinimised(0)
{
    if (m_materialMat->hasUndeformed())
        this->setUndeformed(m_materialMat->getUndeformed());
    if (m_materialMat->hasDeformed())
        this->setDeformed(m_materialMat->getDeformed());

    this->setThickness(m_materialMat->getThickness());

    // Root finder for the fallback; it does not own the wrapped material, which m_materialMat keeps alive
    m_tft = memory::make_shared(new gsMaterialMatrixTFT<dim,T,linear>(m_materialMat.get()));

    m_options.addReal("Eta","Regularisation: fraction of the unmodified stress/tangent kept in wrinkled and slack states; must lie in [0,1] (Zhang 2026 Eq. 40)",0);
    m_options.addReal("SlackMultiplier","Slack state: S_eff = SlackMultiplier*S if > 0, else Eta*S; must be >= 0",0);
    m_options.addReal("S2crit","Display-only threshold: a wrinkled point is reported wrinkled in eval3D_tensionfield iff Smin < S2crit",0);
    m_options.addSwitch("Midsurface","Evaluate the state on the mid-surface only (MatIntegration::Constant); membrane analysis only",false);
    m_options.addSwitch("Exact","Exact (non-truncated) wrinkled stress S(E + gamma n x n) by a local Newton solve (linear=false); no effect for linear=true, where both coincide",false);
    m_options.addReal("ExactTol","Exact solve: residual tolerance relative to the relaxed stress at the iterate (with a round-off floor)",1e-12);
    m_options.addReal("ExactStagTol","Exact solve: relative residual tolerance up to which a stalled iteration is accepted",1e-9);
    m_options.addInt("ExactMaxIt","Exact solve: budget of provider calls (S/C evaluations) per pass",50);
    m_options.addInt("ExactLocalSolve","Exact solve: local solver; 0: scalar pass + 2x2 Newton (legacy); 1: scalar pass + globalised minimisation over P = v (x) v where the 2x2 pass ran and at failed directions (default); 2: minimisation only. Modes 1 and 2 need the wrapped law's energy (eval3D_energy_C); without a W provider the legacy path runs",1);
}

template <short_t dim, class T, bool linear >
void gsMaterialMatrixTFTSplit<dim,T,linear>::setOptions(gsOptionList opt)
{
    m_options.update(opt,gsOptionList::addIfUnknown);
    if (linear && m_options.getSwitch("Exact"))
        gsInfo<<"gsMaterialMatrixTFTSplit: Exact has no effect for linear=true (exact and linearised coincide)\n";
}

template <short_t dim, class T, bool linear >
detail_tftsplit::WrinkleOptions<T> gsMaterialMatrixTFTSplit<dim,T,linear>::wrinkleOptions() const
{
    detail_tftsplit::WrinkleOptions<T> opt;
    opt.tensionField    = m_options.getInt("TensionField");
    opt.eta             = m_options.getReal("Eta");
    opt.slackMultiplier = m_options.getReal("SlackMultiplier");
    opt.exact           = !linear && m_options.getSwitch("Exact");
    opt.exactTol        = m_options.getReal("ExactTol");
    opt.exactStagTol    = m_options.getReal("ExactStagTol");
    opt.exactMaxIt      = m_options.getInt("ExactMaxIt");
    opt.exactLocalSolve = m_options.getInt("ExactLocalSolve");
    GISMO_ENSURE(opt.tensionField==0 || opt.tensionField==1 || opt.tensionField==2,
                 "gsMaterialMatrixTFTSplit: TensionField must be 0, 1 or 2, got "<<opt.tensionField);
    GISMO_ENSURE(opt.exactLocalSolve==0 || opt.exactLocalSolve==1 || opt.exactLocalSolve==2,
                 "gsMaterialMatrixTFTSplit: ExactLocalSolve must be 0, 1 or 2, got "<<opt.exactLocalSolve);
    GISMO_ENSURE(opt.eta>=0 && opt.eta<=1,
                 "gsMaterialMatrixTFTSplit: Option Eta must lie in [0,1], got "<<opt.eta);
    GISMO_ENSURE(opt.slackMultiplier>=0,
                 "gsMaterialMatrixTFTSplit: Option SlackMultiplier must be >= 0, got "<<opt.slackMultiplier);
    return opt;
}

template <short_t dim, class T, bool linear >
std::ostream & gsMaterialMatrixTFTSplit<dim,T,linear>::print(std::ostream &os) const
{
    os  <<"---------------------------------------------------------------------\n"
        <<"----------------Tension-Field Theory Material (split)----------------\n"
        <<"---------------------------------------------------------------------\n\n";

    m_materialMat->print(os);

    os  <<"---------------------------------------------------------------------\n\n";
    return os;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_pstretch(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    return m_materialMat->eval3D_pstretch( patch,u,z );
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_pstretchDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const
{
    return m_materialMat->eval3D_pstretchDir( patch,u,z );
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_detF(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return m_materialMat->eval3D_detF( patch,u,z,out );
}

/*
    Pointwise core.

    Notation: n = (c,s) = (cos th, sin th) is a covector, m = (-s,c), and in engineering Voigt form
        n1 = (c^2, s^2, 2cs)          (n x n),
        n2 = (-sc, sc, c^2 - s^2)     (sym(n x m)),
        n3 = dn2/dth = (s^2 - c^2, c^2 - s^2, -4cs),   dn1/dth = 2 n2.
    Then n1.S = S^{ab} n_a n_b and n2.S = S^{ab} n_a m_b.

    Wrinkle direction. The wrinkled state has S'^{ab} n_a m_b = 0 and, to first order in the
    wrinkle strain gamma n1 (S' = S + gamma C n1, n1.S' = 0), th is the root of
        f(th) = n2.S + gamma n2.C n1,   gamma(th) = -n1.S / (n1.C n1).
    For a material that is isotropic with respect to the reference metric G, S and C share their
    principal directions and the root is the eigendirection of the minor principal stress,
        S^{ab} n_b = lambda_min G^{ab} n_b   <=>   (G S_con) n = lambda_min n,
    which is evaluated in closed form. It is accepted only if |f| <= 1e-10 |S|_2 + 1e-8 (Smax-Smin)
    and n1.C n1 > 0; otherwise the root finder of gsMaterialMatrixTFT::eval_theta is used, and a
    root of that finder is accepted under the same test.

    Tangent. With H.col(d) = (dC/dE_d) n1 (H = 0 for a material linear in E), the derivatives at fixed
    S(E), C(E) read
        dgamma/dE |_th = -(C n1 + gamma H^T n1) / b,                       b = n1.C n1
        dgamma/dth |_E = -(2/b) (n2.S + 2 gamma n2.C n1)
        df/dE |_th     = C n2 + (n2.C n1) dgamma/dE|_th + gamma H^T n2
        df/dth |_E     = n3.S + (dgamma/dth)(n2.C n1) + gamma (n3.C n1 + 2 n2.C n2)
        dth/dE         = -(df/dE) / (df/dth)
        dgamma/dE      = dgamma/dE|_th + (dgamma/dth) dth/dE.
    With D_kd = d(gamma n1)_k / dE_d = n1_k dgamma/dE_d + 2 gamma n2_k dth/dE_d, the consistent tangent
    of S' = S + gamma C n1 is (A_ij = dS'_i/dE_j)
        J = C + C D + gamma H,
    which is not symmetric in general. See Verhelst et al., CMAME 441 (2025) 117955, Sect. 4, and
    Zhang & Kiendl (2026), Eq. 34-38. The regularised result is S_eff = (1-eta) S' + eta S and
    C_eff = (1-eta) J + eta C (Eq. 35, 40).

    Complexity: O(1) per point, plus one eval_theta call (O(50) residual evaluations) on guard failure.

    Accuracy of the linearised mode. With H as above,
        S(E + gamma n1) = S + gamma C n1 + 1/2 gamma^2 H n1 + O(gamma^3),
    so the truncation error of S' is 1/2 gamma^2 H n1: O(gamma^2) absolute, O(gamma) relative to
    gamma C n1; it vanishes for a material linear in E.

    Exact mode (opt.exact, linear=false). Instead of the Taylor-truncated S' = S + gamma C n1, the
    elastic strain E'(gamma,th) = E + gamma n1(th) is used and (gamma,th) solve
        R(gamma,th) = (n1.S(E'), n2.S(E')) = 0,
    the stationarity of W(gamma,th) = W(E + gamma n1(th)) on the rank-one family; under the class-level
    hypotheses, a solution with gamma > 0 and n4.S' > 0 is the rank-one point of the KKT partition of
    Pipkin's W_r = min_{P PSD} W(E + P). S(E') and
    C(E') are evaluated by the wrapped material at the prescribed metric
        Cmat' = (G11 + 2E'_0, G22 + 2E'_1, G12 + E'_2)  (= Cmat + 2 gamma (c^2, s^2, cs)),
    with the tensor g12 (E'_2 is the engineering shear). Newton, with dE'/dgamma = n1, dE'/dth = 2 gamma n2:
        J11 = n1.C'n1,                      J12 = 2 n2.S' + 2 gamma n1.C'n2,
        J21 = n2.C'n1,                      J22 = n3.S' + 2 gamma n2.C'n2.
    With W(gamma,th) = W(E + gamma n1(th)): dW/dgamma = n1.S(E'), dW/dth = 2 gamma n2.S(E'), so that
    at R = 0 the Hessian of W(gamma,th) is diag(1, 2 gamma) Jac0.
    For a material isotropic w.r.t. G, S and C are coaxial, R2 vanishes along gamma at the closed-form th
    and a scalar Newton step -R1/J11 in gamma suffices; the 2x2 system is used if the guard rejected the
    closed form or if |R2| exceeds the tolerance at the scalar solution (ExactLocalSolve = 0; under
    ExactLocalSolve 1 the minimisation of _exactMinimise takes over at those points, and under 2 it replaces
    both passes).

    Tolerance. At an iterate E' = E + gamma n1 with S' = S(E'), C' = C(E'), a residual is judged against
        floorR = 64 eps ||C'||_F max(1,||E'||_2),   tol(k) = max(k ||S'||_2, floorR),
    with tolC = tol(ExactTol) for convergence and tolS = tol(ExactStagTol) for stall acceptance. The
    reference is the relaxed stress S' at the iterate itself (not the unrelaxed S at E, which is large
    where the wrinkled stress is small), so the test is relative to the quantity being solved for. floorR is
    the round-off floor: the round-off of E' maps to a residual of about eps ||C'|| max(1,||E'||), which
    takes over where S' -> 0 (max(1,||E'||) stands for the metric
    scale, assuming ||G|| = O(1)). Stall acceptance at tolS is meant for providers whose S' and C' carry noise
    above tolC (e.g. an inner iteration of the wrapped law stopped at a loose tolerance). Each pass has a budget
    of ExactMaxIt provider calls; the options satisfy ExactStagTol >= ExactTol > 0, ExactMaxIt >= 1.

    Scalar pass. g(gamma) = n1.S(E + gamma n1) at the closed-form th is solved by a bracketed
    Newton-bisection (rtsafe, Numerical Recipes) started at the linearised gamma. If g(gamma_lin) > 0 the
    g(0) is evaluated: g(0) >= 0 means the point is taut at gamma = 0 (no wrinkle), otherwise [0, gamma_lin]
    brackets the root. If g(gamma_lin) < 0 the upper end is open and the iterate is expanded by Newton
    steps limited to a doubling of gamma. A linearised gamma <= 0 (strain criterion only) is decided by g(0) alone. A non-finite S' or C' at a
    trial means the trial is too far (upper limit) and the next trial is the bracket midpoint. The scalar
    pass ends on |g| <= tolC (followed by one polishing Newton step, kept if |g| does not grow), or, after a
    Newton step that fails to reduce |g| below the three previous evaluations, or when the bracket has
    collapsed to round-off, on the best iterate if its |g| is within its own tolS (stall acceptance). If
    n2.S' is within the tolerance at which the scalar iterate was accepted the solve is done, else the 2x2
    pass continues from (gamma, th) (ExactLocalSolve = 0; under ExactLocalSolve 1 the scalar pass hands off
    to the minimisation instead, see _exactSolve and _exactMinimise).

    2x2 pass. Plain Newton in (gamma,th), converged on ||R||_inf <= tolC (followed by one polishing step); a
    residual that does not decrease over two iterations after the third ends it with the best iterate if
    that is within its tolS. Non-finite S' or C', or a singular Jacobian, end it unconverged.

    Acceptance. A stationary point with gamma > 0 is tested for a minimum first (M11 > 0, det Jac0 > 0, see
    below); one that is not a minimum is a saddle, unless its wrinkle tension sw is not positive and the
    point is slack within the band below. A minimum needs sw > 0. A point found taut at gamma = 0, or at
    gamma <= 0 within the taut band, is returned as taut (state 1, unmodified S, E, C); a point with sw <= 0
    within the slack band as slack. Bands: with band = 1e-8 and the principal values of the unmodified
    S (or E, for the strain criterion), taut means Smin >= -band max(|Smin|,|Smax|) (Emin for TensionField 0),
    slack means Emax <= band max(|Emin|,|Emax|) (Smax for TensionField 1). A stationary point outside the
    band on a branch inconsistent with the classifier (gamma <= 0, or a minimum with sw <= 0) is a wrong
    branch. Non-convergence, a saddle and a wrong branch all return the linearised result.

    Tangent. By the envelope theorem dW_r/dE = S(E'), hence the tangent of the relaxed energy is symmetric
    whenever C' is. Differentiating R(gamma(E),th(E);E) = 0 gives d(gamma,th)/dE = -Jac^{-1} K C' with
    K = [n1^T; n2^T], and dS'/dE = C'(I + [n1, 2 gamma n2] d(gamma,th)/dE). At R = 0, n2.S' = 0 and
    n3.S' = n4.S' =: sw, n4 = n3 + n1 = (s^2,c^2,-2cs) (m x m), the tension along the wrinkles. With
        M = K C' K^T,  Jac0 = [[M11, 2 gamma M12], [M21, sw + 2 gamma M22]],  X = diag(1,2 gamma) Jac0^{-1},
        J_ex = C' - [C'n1, C'n2] X [n1^T C'; n2^T C'].
    For gamma > 0 this equals C' - C'K^T (K C' K^T + diag(0, sw/(2 gamma)))^{-1} K C'; the unscaled form is
    finite for gamma -> 0. For gamma > 0 the stationary point is a strict local minimum of W(gamma,th) on the rank-one
    family iff M11 > 0 and det Jac0 > 0.
    This is guaranteed if sw > 0 and C' is positive definite on span{n1, n2}.
    In the major principal direction stationarity gives gamma < 0 (to first order -Smax/(n1.C n1)), which
    is not an admissible wrinkle; hence the start from the minor direction. The linearised gamma is the
    first scalar Newton step in gamma from gamma = 0 at fixed th (not so for the 2x2 iteration at a
    non-coaxial th, where R2 != 0 at gamma = 0). In every mode the tangent jumps at the taut/wrinkled
    boundary (W_r is C^1, not C^2).

    Energy. With W the energy of the wrapped material (S = dW/dE) the effective energy is
        taut:   W_eff = W(E),                       slack:  W_eff = sigma W(E),
        exact wrinkled:  W_eff = (1-Eta) W(E') + Eta W(E),
    sigma = SlackMultiplier if > 0, else Eta. At an accepted exact state R = 0, so by the envelope theorem
    d/dE W(E') = S' (the wrinkle variables drop out), and dW_eff/dE = (1-Eta) S' + Eta S = S_eff. For linear=true
    W is quadratic, W(E) = 1/2 E.S, the linearised state is the exact minimiser along n1 and
    W(E') = W(E) + gamma n1.S + 1/2 gamma^2 b, with b = n1.C n1 (no cancellation of the E'-quadratic form).
    The remaining wrinkled points (linearised mode of linear=false, exact fallbacks, and under
    ExactLocalSolve = 0 indefinite b or a failed direction) have a non-symmetric tangent, hence no potential:
    W_eff is NaN there, as for non-finite data. Under ExactLocalSolve 1 or 2 such a point enters the
    minimisation and, if accepted, has W_eff like any accepted exact state.

    Complexity: two wrapped _C evaluations per provider call (each with an inner thickness Newton for
    compressible laws); typically 2-5 calls from the linearised start plus one polishing call per accepted
    solve, at most ExactMaxIt calls per pass. Under ExactLocalSolve 1 and 2 the minimisation adds one wrapped
    eval3D_energy_C evaluation per W call and has its own budget of ExactMaxIt S/C calls (see _exactMinimise).
*/
template <short_t dim, class T, bool linear >
detail_tftsplit::WrinkleResult<T> gsMaterialMatrixTFTSplit<dim,T,linear>::wrinkle(
                                                const gsMatrix<T> & S,
                                                const gsMatrix<T> & C,
                                                const gsMatrix<T> & E,
                                                const gsMatrix<T> & Gcov,
                                                const detail_tftsplit::WrinkleOptions<T> & opt,
                                                bool needTangent,
                                                const dC_provider & dCdE,
                                                const SC_provider & SCatE,
                                                const W_provider & WatE,
                                                bool needEnergy) const
{
    using namespace detail_tftsplit;
    GISMO_ASSERT(S.rows()*S.cols()==3,"S must have 3 entries");
    GISMO_ASSERT(E.rows()*E.cols()==3,"E must have 3 entries");
    GISMO_ASSERT(Gcov.rows()==2 && Gcov.cols()==2,"Gcov must be 2x2");
    GISMO_ASSERT((C.rows()==3 && C.cols()==3) || (C.rows()*C.cols()==9),"C must be 3x3 or 9x1");
    if (needEnergy && !linear && !WatE)
        GISMO_ERROR("gsMaterialMatrixTFTSplit::wrinkle: the energy was requested without a W provider");

    const T nan = std::numeric_limits<T>::quiet_NaN();
    const T pi  = static_cast<T>(EIGEN_PI);

    gsVector<T,3> S3, E3;
    gsMatrix<T,3,3> C3;
    for (index_t i = 0; i!=3; ++i)
    {
        S3(i) = S.data()[i];
        E3(i) = E.data()[i];
    }
    if (C.rows()==3 && C.cols()==3)
        for (index_t l = 0; l!=3; ++l)
            for (index_t i = 0; i!=3; ++i)
                C3(i,l) = C(i,l);
    else
        for (index_t l = 0; l!=3; ++l)
            for (index_t i = 0; i!=3; ++i)
                C3(i,l) = C.data()[i+3*l];

    WrinkleResult<T> res;
    res.state = 1;
    res.theta = res.gamma = res.gammaHat = res.residual = nan;
    res.W = nan;
    res.guardFallback = res.tftFailure = false;
    res.exactFallback = res.exactReclassified = false;
    res.exactSaddle = res.exactWrongBranch = res.exactStagnated = res.exactSkipped = res.exactUsed2x2 = false;
    res.exactMinimised = false;
    res.exactIterations = 0;
    res.Gcov = Gcov;

    // Principal values with respect to the reference metric
    const Mat2<T> Gcon = invert2x2(Gcov);
    const Mat2<T> MS   = stressPencil(Gcov,S3(0),S3(1),S3(2));
    const Mat2<T> ME   = strainPencil(Gcon,E3(0),E3(1),E3(2));
    principalValues(MS,res.Smin,res.Smax);
    principalValues(ME,res.Emin,res.Emax);

    // State; the slack boundary is strict
    const T tauS = T(1e-12) * math::max(math::abs(res.Smin),math::abs(res.Smax));
    const T tauE = T(1e-12) * math::max(math::abs(res.Emin),math::abs(res.Emax));
    bool taut, slack;
    if (opt.tensionField==0)
    {
        taut  = res.Emin > 0 || math::abs(res.Emin) <= tauE;
        slack = !taut && res.Emax < -tauE;
    }
    else if (opt.tensionField==1)
    {
        taut  = res.Smin > 0 || math::abs(res.Smin) <= tauS;
        slack = !taut && res.Smax < -tauS;
    }
    else
    {
        taut  = res.Smin > 0 || math::abs(res.Smin) <= tauS;
        slack = !taut && res.Emax < -tauE;
    }

    // Non-finite data (e.g. an inverted element) is passed through unmodified, so that it stays
    // visible downstream instead of being fed to the root finder.
    bool finiteData = true;
    if (!(math::isfinite(res.Smin) && math::isfinite(res.Smax) && math::isfinite(res.Emin) && math::isfinite(res.Emax)))
    {
        finiteData = false;
        taut  = true;
        slack = false;
    }

    // W of the wrapped material at the elastic strain Ee; only called if needEnergy
    auto energyAt = [&](const gsVector<T,3> & Ee) -> T
    {   // linear: S3 = C3 E3, so W = 1/2 E.S is the wrapped law's potential
        if (linear) return T(0.5) * Ee.dot(S3);   // only ever called with Ee = E3 for linear
        gsMatrix<T> Em(3,1); Em = Ee; return WatE(Em);
    };

    if (taut)
    {
        res.state = 1;
        res.S = S3;
        res.E = E3;
        if (needTangent)
            res.C = C3;
        if (needEnergy && finiteData)
            res.W = energyAt(E3);
        return res;
    }

    if (slack)
    {
        const T sigma = opt.slackMultiplier > 0 ? opt.slackMultiplier : opt.eta;
        res.state = -1;
        res.S = sigma * S3;
        res.E = sigma * E3;
        if (needTangent)
            res.C = sigma * C3;
        if (needEnergy)
            res.W = sigma * energyAt(E3);
        return res;
    }

    res.state = 0;

    // b = n1.C n1, gamma = -n1.S/b, f = n2.S + gamma n2.C n1 at the direction theta
    auto evalAt = [&](const T th, T & b, T & gam, T & f)
    {
        const T c = math::cos(th), s = math::sin(th);
        const VoigtDirs<T> n(c,s);
        const gsVector<T,3> Cn1 = C3 * n.n1;
        b   = n.n1.dot(Cn1);
        gam = -n.n1.dot(S3) / b;
        f   = n.n2.dot(S3) + gam * n.n2.dot(Cn1);
    };

    // Closed-form direction: right eigenvector of MS for the minor stress
    const T lam = res.Smin;
    const T r0a = MS.a - lam, r0b = MS.b;
    const T r1a = MS.c,       r1b = MS.d - lam;
    const T r0n = math::sqrt(r0a*r0a + r0b*r0b);
    const T r1n = math::sqrt(r1a*r1a + r1b*r1b);
    T nx, ny;
    if (r0n >= r1n)
    {
        nx = MS.b;
        ny = lam - MS.a;
    }
    else
    {
        nx = lam - MS.d;
        ny = MS.c;
    }
    const T nlen = math::sqrt(nx*nx + ny*ny);
    // The negated comparisons also reject NaN
    const bool degenerate = !((res.Smax - res.Smin) > T(1e-12) * math::max(math::abs(res.Smin),math::abs(res.Smax)))
                            || !(math::max(r0n,r1n) > 0) || !(nlen > 0);

    auto wrapAngle = [&](T th)
    {
        if (th < 0)   th += pi;
        if (th >= pi) th -= pi;
        return th;
    };

    T theta = nan, b = nan, gam = nan, f = nan;
    const T tol = T(1e-10) * S3.norm() + T(1e-8) * (res.Smax - res.Smin);
    bool cfOK = false;
    T thetaCF = nan, bCF = nan, gamCF = nan, fCF = nan;
    if (!degenerate)
    {
        thetaCF = wrapAngle(math::atan2(ny/nlen,nx/nlen));
        evalAt(thetaCF,bCF,gamCF,fCF);
        cfOK = bCF > 0 && math::abs(fCF) <= tol;
    }

    bool tftFailed = false;
    if (cfOK)
    {
        theta = thetaCF; b = bCF; gam = gamCF; f = fCF;
    }
    else
    {
        // The closed form is not the root of f: take the root finder of TFT
        res.guardFallback = true;
#       pragma omp atomic
        ++m_nThetaFallbacks;

        bool warn = false;
#       pragma omp critical (gsMaterialMatrixTFTSplit_guardFallback)
        {
            static bool warned = false;
            warn = !warned;
            warned = true;
        }
        if (warn)
            gsWarn<<"gsMaterialMatrixTFTSplit: the closed-form wrinkle direction was rejected by the isotropy guard; "
                    "falling back to the root finder of gsMaterialMatrixTFT (reported once, see numThetaFallbacks())\n";

        gsMatrix<T> C9(9,1), Sm(3,1), Em(3,1);
        for (index_t l = 0; l!=3; ++l)
            for (index_t i = 0; i!=3; ++i)
                C9(i+3*l,0) = C3(i,l);
        for (index_t i = 0; i!=3; ++i)
        {
            Sm(i,0) = S3(i);
            Em(i,0) = E3(i);
        }
        const T thetaT = wrapAngle(m_tft->eval_theta(C9,Sm,Em)(0,0));
        T bT, gamT, fT;
        evalAt(thetaT,bT,gamT,fT);

        if (bT > 0 && math::abs(fT) <= tol)
        {
            theta = thetaT; b = bT; gam = gamT; f = fT;
        }
        else
        {
            // Neither root passes: keep the closed form when it is usable
            tftFailed = true;
            if (!degenerate && bCF > 0)
            {
                theta = thetaCF; b = bCF; gam = gamCF; f = fCF;
            }
            else
            {
                theta = thetaT; b = bT; gam = gamT; f = fT;
            }
        }
    }

    // A wrapped tangent that is not positive definite along n1 has no wrinkle correction
    const bool indefinite = !(b > 0);
    if (indefinite)
        tftFailed = true;

    if (tftFailed)
    {
        res.tftFailure = true;
#       pragma omp atomic
        ++m_nTFTFailures;

        bool warn = false;
#       pragma omp critical (gsMaterialMatrixTFTSplit_tftFailure)
        {
            static bool warned = false;
            warn = !warned;
            warned = true;
        }
        if (warn)
            gsWarn<<"gsMaterialMatrixTFTSplit: no admissible wrinkle direction was found (root finder residual too large "
                    "or n1^T C n1 <= 0) (reported once, see numTFTFailures())\n";
    }

    res.theta    = theta;
    res.residual = math::abs(f);

    // An empty W provider runs the legacy local solver in every mode
    const index_t localSolve = WatE ? opt.exactLocalSolve : 0;
    const bool directionOK = !indefinite && !tftFailed;

    if (!linear && opt.exact && !directionOK && localSolve == 0)
    {
        res.exactSkipped = true;
#       pragma omp atomic
        ++m_nExactSkipped;
    }

    if (!linear && opt.exact && (directionOK || localSolve != 0))
    {
        if (!SCatE)
            GISMO_ERROR("gsMaterialMatrixTFTSplit::wrinkle: exact mode requires an S/C provider");
        if (!(opt.exactStagTol >= opt.exactTol && opt.exactTol > 0 && opt.exactMaxIt >= 1))
            GISMO_ERROR("gsMaterialMatrixTFTSplit::wrinkle: the options ExactTol, ExactStagTol and ExactMaxIt must satisfy "
                        "ExactStagTol >= ExactTol > 0 and ExactMaxIt >= 1, got ExactTol = "<<opt.exactTol
                        <<", ExactStagTol = "<<opt.exactStagTol<<", ExactMaxIt = "<<opt.exactMaxIt);
        if (!(opt.exactLocalSolve >= 0 && opt.exactLocalSolve <= 2))
            GISMO_ERROR("gsMaterialMatrixTFTSplit::wrinkle: ExactLocalSolve must be 0, 1 or 2, got "<<opt.exactLocalSolve);

        T gE = gam, thE = theta, resNorm = nan;
        gsVector<T,3> Spx;
        gsMatrix<T,3,3> Cpx;
        index_t nEval = 0;
        bool tautAtZero = false, stagnated = false, used2x2 = false;

        // Taut/slack bands: a stationary point of the wrong kind is still accepted as taut/slack if the
        // classifier's own quantity is within band of the criterion boundary (relative to the larger principal value)
        const T band = T(1e-8);
        const T Sref = math::max(math::abs(res.Smin),math::abs(res.Smax));
        const T Eref = math::max(math::abs(res.Emin),math::abs(res.Emax));
        const bool tautBand  = (opt.tensionField==0) ? (res.Emin >= -band * Eref) : (res.Smin >= -band * Sref);
        const bool slackBand = (opt.tensionField==1) ? (res.Smax <=  band * Sref) : (res.Emax <=  band * Eref);

        bool minimised = false, solved = false;
        if (localSolve == 0)
            solved = this->template _exactSolve<linear>(E3,opt,SCatE,cfOK,gE,thE,Spx,Cpx,resNorm,nEval,
                                                         tautAtZero,stagnated,used2x2);
        else
        {
            index_t handOff = 0;
            if (localSolve == 1 && directionOK && cfOK)
                solved = this->template _exactSolve<linear>(E3,opt,SCatE,cfOK,gE,thE,Spx,Cpx,resNorm,nEval,
                                                             tautAtZero,stagnated,used2x2,&handOff);
            if (localSolve == 2 || !directionOK || !cfOK || handOff != 0)
            {
                // Start of the minimisation: the first matching rule wins
                bool haveStart = false, admissible = false;
                T theta0 = nan, gamma0 = nan;
                if (handOff == 1 || handOff == 2)
                {
                    theta0 = thE; gamma0 = gE; admissible = (handOff == 1);
                    haveStart = true;
                }
                else if (!indefinite && gam > 0)
                {
                    theta0 = theta; gamma0 = gam; admissible = false;
                    haveStart = true;
                }
                else if (res.Smin >= -band * Sref)
                {
                    tautAtZero = true;
                    solved = true;
                }
                else
                {
                    const T mid = T(0.5) * (S3(0) + S3(1));
                    const T dif = T(0.5) * (S3(0) - S3(1));
                    const T q   = S3(2);
                    const T rad = math::sqrt(dif*dif + q*q);
                    const T lamS = mid - rad;
                    if (lamS < 0)
                    {
                        // Eigenvector of Shat(E) = [[S0,S2],[S2,S1]] for lamS
                        const T ph = T(0.5) * math::atan2(2*q,S3(0) - S3(1));
                        theta0 = math::atan2(math::cos(ph),-math::sin(ph));
                        const VoigtDirs<T> n0(math::cos(theta0),math::sin(theta0));
                        const T b0 = n0.n1.dot(C3 * n0.n1);
                        if (b0 > 0 && math::isfinite(b0))
                            gamma0 = -lamS / b0;
                        else
                            gamma0 = (math::abs(res.Emin) > 0) ? math::abs(res.Emin) : Eref;
                        admissible = false;
                        haveStart = true;
                    }
                }
                if (haveStart)
                {
                    minimised = true;
                    solved = this->template _exactMinimise<linear>(E3,opt,SCatE,WatE,theta0,gamma0,admissible,gE,thE,Spx,Cpx,
                                                                   resNorm,nEval,stagnated);
                }
            }
        }
        res.exactIterations = nEval;
        res.exactStagnated  = stagnated;
        res.exactUsed2x2    = used2x2;
        res.exactMinimised  = minimised;
        if (minimised)
        {
#           pragma omp atomic
            ++m_nExactMinimised;
        }
        if (stagnated)
        {
#           pragma omp atomic
            ++m_nExactStagnated;
        }
        if (used2x2)
        {
#           pragma omp atomic
            ++m_nExactUsed2x2;
        }

        // Taut: the unmodified point; slack: scaled by SlackMultiplier (or Eta)
        auto reclassify = [&](const bool toTaut)
        {
            const T sigma = opt.slackMultiplier > 0 ? opt.slackMultiplier : opt.eta;
            const T fac   = toTaut ? T(1) : sigma;
            res.state = toTaut ? 1 : -1;
            res.theta = res.gamma = res.gammaHat = res.residual = nan;
            res.exactReclassified = true;
            res.S = fac * S3;
            res.E = fac * E3;
            if (needTangent)
                res.C = fac * C3;
            if (needEnergy)
                res.W = fac * energyAt(E3);
        };

        // 0: accepted, 1: not converged, 2: saddle, 3: wrong branch, 4: slack within the slack band
        int failure = 0;
        if (!solved)
            failure = 1;
        else if (tautAtZero)
        {
            reclassify(true);
#           pragma omp atomic
            ++m_nExactReclassifiedTaut;
            return res;
        }
        else if (!(gE > 0))
        {
            if (tautBand)
            {
                reclassify(true);
#               pragma omp atomic
                ++m_nExactReclassifiedTaut;
                return res;
            }
            failure = 3;
        }
        else
        {
            thE = wrapAngle(thE);
            const T ce = math::cos(thE), se = math::sin(thE);
            const VoigtDirs<T> ne(ce,se);
            const gsVector<T,3> ne4 = ne.n3 + ne.n1;
            const T sw = ne4.dot(Spx);

            const gsVector<T,3> Cn1 = Cpx * ne.n1;
            const gsVector<T,3> Cn2 = Cpx * ne.n2;
            const T M11 = ne.n1.dot(Cn1), M12 = ne.n1.dot(Cn2);
            const T M21 = ne.n2.dot(Cn1), M22 = ne.n2.dot(Cn2);
            const T a = M11, bb = 2 * gE * M12, cc = M21, d = sw + 2 * gE * M22;
            const T det0 = a*d - bb*cc;

            const bool minimum = math::isfinite(M11) && math::isfinite(det0) && M11 > 0 && det0 > 0;
            if (!minimum)
                failure = (!(sw > 0) && slackBand) ? 4 : 2;
            else if (!(sw > 0))
                failure = slackBand ? 4 : 3;

            if (failure == 4)
            {
                reclassify(false);
#               pragma omp atomic
                ++m_nExactReclassifiedSlack;
                return res;
            }

            if (failure == 0)
            {
                res.theta    = thE;
                res.gamma    = gE;
                res.gammaHat = gE * (ce*(Gcon.a*ce + Gcon.b*se) + se*(Gcon.c*ce + Gcon.d*se));
                res.residual = resNorm;

                res.S = (1 - opt.eta) * Spx + opt.eta * S3;
                res.E = E3 + gE * ne.n1;
                if (needEnergy)
                    res.W = (1 - opt.eta) * WatE(res.E) + opt.eta * energyAt(E3);
                if (needTangent)
                {
                    // X = diag(1,2 gamma) Jac0^{-1}
                    gsMatrix<T,2,2> X;
                    X << d, -bb, 2*gE*(-cc), 2*gE*a;
                    X /= det0;
                    gsMatrix<T,3,2> L;
                    L.col(0) = Cn1;
                    L.col(1) = Cn2;
                    gsMatrix<T,2,3> Rt;
                    Rt.row(0) = (Cpx.transpose() * ne.n1).transpose();
                    Rt.row(1) = (Cpx.transpose() * ne.n2).transpose();
                    const gsMatrix<T,3,3> Jex = Cpx - L * X * Rt;
                    res.C = (1 - opt.eta) * Jex + opt.eta * C3;
                }
                return res;
            }
        }

        // Not converged, saddle or wrong branch: the linearised result below is returned
        res.exactFallback = true;
        if (failure == 1)
        {
#           pragma omp atomic
            ++m_nExactNonconverged;
        }
        else if (failure == 2)
        {
            res.exactSaddle = true;
#           pragma omp atomic
            ++m_nExactSaddle;
        }
        else
        {
            res.exactWrongBranch = true;
#           pragma omp atomic
            ++m_nExactWrongBranch;
        }

        bool warn = false;
#       pragma omp critical (gsMaterialMatrixTFTSplit_exactFallback)
        {
            static bool warned = false;
            warn = !warned;
            warned = true;
        }
        if (warn)
            gsWarn<<"gsMaterialMatrixTFTSplit: the exact wrinkling solve did not converge, reached a saddle or a wrong branch; "
                    "the linearised result is used (reported once, see numExactFallbacks())\n";
    }

    gsVector<T,3> Sp = S3, Ep = E3;
    gsMatrix<T,3,3> J = C3;
    if (!indefinite)
    {
        const T c = math::cos(theta), s = math::sin(theta);
        const VoigtDirs<T> n(c,s);
        const gsVector<T,3> & n1 = n.n1;
        const gsVector<T,3> & n2 = n.n2;
        const gsVector<T,3> & n3 = n.n3;
        const gsVector<T,3> Cn1 = C3 * n1;

        // S' = S + gamma C n1 and the elastic strain E' = E + gamma n (x) n
        Sp += gam * Cn1;
        Ep += gam * n1;

        res.gamma    = gam;
        res.gammaHat = gam * (c*(Gcon.a*c + Gcon.b*s) + s*(Gcon.c*c + Gcon.d*s));

        if (linear && needEnergy && !tftFailed)
        {
            const T WE = T(0.5) * E3.dot(S3);
            res.W = (1 - opt.eta) * (WE + gam * n1.dot(S3) + T(0.5) * gam * gam * b) + opt.eta * WE;
        }

        if (needTangent)
        {
            gsMatrix<T,3,3> H = gsMatrix<T,3,3>::Zero();
            if (dCdE)
            {
                std::vector<gsMatrix<T>> dC;
                dCdE(dC);
                GISMO_ENSURE(dC.size()==3,"gsMaterialMatrixTFTSplit: the dC provider must return 3 matrices");
                for (index_t d = 0; d!=3; ++d)
                    H.col(d) = dC[d] * n1;
            }

            const gsVector<T,3> Cn2   = C3 * n2;
            const gsVector<T,3> HTn1  = H.transpose() * n1;
            const gsVector<T,3> HTn2  = H.transpose() * n2;
            const T n2Cn1 = n2.dot(Cn1);
            const T n3Cn1 = n3.dot(Cn1);
            const T n2Cn2 = n2.dot(Cn2);

            const gsVector<T,3> dgdE  = -(Cn1 + gam * HTn1) / b;
            const T dgdth             = -(2/b) * (n2.dot(S3) + 2 * gam * n2Cn1);
            const gsVector<T,3> dfdE  = Cn2 + n2Cn1 * dgdE + gam * HTn2;
            const T dfdth             = n3.dot(S3) + dgdth * n2Cn1 + gam * (n3Cn1 + 2 * n2Cn2);
            const gsVector<T,3> dthdE = -dfdE / dfdth;
            const gsVector<T,3> dgamdE = dgdE + dgdth * dthdE;

            const gsMatrix<T,3,3> D = n1 * dgamdE.transpose() + (2 * gam) * (n2 * dthdE.transpose());
            J = C3 + C3 * D + gam * H;
        }
    }

    // Regularisation (Zhang 2026 Eq. 35/40)
    res.S = (1 - opt.eta) * Sp + opt.eta * S3;
    res.E = Ep;
    if (needTangent)
        res.C = (1 - opt.eta) * J + opt.eta * C3;
    return res;
}

/*
    Exact wrinkled state; see the theory comment of wrinkle(). Two passes, each with a budget of
    opt.exactMaxIt provider calls (all calls also count in nEval).

    Scalar pass (only if \a scalar0, i.e. the closed-form theta passed the isotropy guard): a bracketed
    Newton-bisection (rtsafe) for g(gamma) = n1.S(E + gamma n1) at the fixed theta, started at the
    incoming linearised gamma. It ends with a taut-at-zero verdict (tautAtZero, g(0) >= 0), a converged
    iterate (|g| <= tolC, followed by one polishing Newton step), or the best iterate accepted at the
    stall tolerance (stagnated). If n2.S' at the result is within the tolerance the iterate was accepted
    at, the solve is done; otherwise the 2x2 pass continues from (gamma, theta). With \a handOff, the routine
    instead returns false at that point, with *handOff = 1 (the scalar pass converged) or 2 (it was accepted
    at the stall tolerance) and gamma set to the kept scalar iterate; the 2x2 pass does not run and used2x2
    stays false. *handOff is 0 on every other exit; \a handOff requires \a scalar0.

    2x2 pass (Newton in (gamma,theta), not globalised): converged when ||R||_inf <= tolC, followed by one
    polishing step; a stalled residual ends it with the best iterate if that is within the stall tolerance.

    Returns true if converged, accepted at the stall tolerance or taut at gamma = 0; false otherwise (also on
    non-finite S' or C'). gamma and theta are not modified unless true is returned, and not at all for
    tautAtZero. On true, Sp, Cp are S', C' at the returned iterate and resNorm = max(|n1.S'|, |n2.S'|) there.
    used2x2 is set whenever the 2x2 pass ran.
*/
template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<!_linear, bool>::type
gsMaterialMatrixTFTSplit<dim,T,linear>::_exactSolve(const gsVector<T,3> & E3, const detail_tftsplit::WrinkleOptions<T> & opt,
                                                    const SC_provider & SCatE, const bool scalar0,
                                                    T & gamma, T & theta, gsVector<T,3> & Sp, gsMatrix<T,3,3> & Cp, T & resNorm,
                                                    index_t & nEval, bool & tautAtZero, bool & stagnated, bool & used2x2, index_t * handOff) const
{
    using namespace detail_tftsplit;
    const T nan = std::numeric_limits<T>::quiet_NaN();
    const T eps = std::numeric_limits<T>::epsilon();
    const index_t maxIt = opt.exactMaxIt;

    GISMO_ASSERT(!handOff || scalar0,"gsMaterialMatrixTFTSplit::_exactSolve: handOff requires the scalar pass");
    tautAtZero = stagnated = used2x2 = false;
    if (handOff)
        *handOff = 0;

    T g = gamma, th = theta;
    gsMatrix<T> Eel(3,1), Se, Ce;
    index_t nPass = 0;

    // One provider call at E' = E + gg n1; false if S' or C' is not finite. tolC/tolS are the tolerances at
    // this iterate; floorR is the round-off of E' mapped to a stress, eps ||C'|| max(1,||E'||) (scaled by 64)
    auto provide = [&](const T gg, const VoigtDirs<T> & n, T & tolC, T & tolS) -> bool
    {
        Eel = E3 + gg * n.n1;
        SCatE(Eel,Se,Ce);
        ++nEval;
        ++nPass;

        for (index_t i = 0; i != 3; ++i)
            Sp(i) = Se.data()[i];
        if (Ce.rows() == 3 && Ce.cols() == 3)
            for (index_t l = 0; l != 3; ++l)
                for (index_t i = 0; i != 3; ++i)
                    Cp(i,l) = Ce(i,l);
        else
            for (index_t l = 0; l != 3; ++l)
                for (index_t i = 0; i != 3; ++i)
                    Cp(i,l) = Ce.data()[i+3*l];

        if (!(Sp.allFinite() && Cp.allFinite()))
            return false;

        const T Sn     = Sp.norm();
        const T floorR = 64 * eps * Cp.norm() * math::max(T(1),Eel.norm());
        tolC = math::max(opt.exactTol * Sn,floorR);
        tolS = math::max(opt.exactStagTol * Sn,floorR);
        return true;
    };

    // An evaluated iterate; gv = n1.S' (scalar pass) or ||R||_inf (2x2 pass), a = n1.C'n1
    struct Iterate
    {
        T g, th, gv, a, tolC, tolS, tolAcc;
        gsVector<T,3> S;
        gsMatrix<T,3,3> C;
    };
    Iterate cur{}, best{}, kept{};

    if (scalar0)
    {
        const VoigtDirs<T> n(math::cos(th),math::sin(th));
        T tolC, tolS;

        // gamma_lin <= 0 (strain-based criterion only): the only admissible outcome is taut at gamma = 0
        if (!(g > 0))
        {
            if (!provide(T(0),n,tolC,tolS))
                return false;
            if (n.n1.dot(Sp) >= 0)
            {
                tautAtZero = true;
                return true;
            }
            return false;
        }

        enum Bracket { noBracket, openUpper, closedBracket } bs = noBracket;
        T lo = 0, hi = 0, dx = 0, dxold = 0;
        T hist[3] = {nan,nan,nan};  // |g| of the three most recent finite evaluations at gamma > 0
        index_t nh = 0;
        bool first = true, byNewton = false, converged = false;
        T gt = g;

        auto closeBracket = [&]()
        {
            if (bs != closedBracket)
            {
                bs = closedBracket;
                dx = dxold = hi - lo;
            }
        };
        auto bisect = [&]()
        {
            dxold = dx;
            dx    = T(0.5) * (hi - lo);
            gt    = lo + dx;
            byNewton = false;
        };

        for (;;)
        {
            if (nPass >= maxIt)
                return false;

            if (!provide(gt,n,tolC,tolS))
            {
                // A non-finite trial is too far
                if (first)
                    return false;
                hi = gt;
                closeBracket();
                bisect();
                continue;
            }

            cur.g = gt;
            cur.gv = n.n1.dot(Sp);
            cur.a = n.n1.dot(Cp * n.n1);
            cur.tolC = tolC;
            cur.tolS = tolS;
            cur.S = Sp;
            cur.C = Cp;
            if (first || math::abs(cur.gv) < math::abs(best.gv))
                best = cur;

            if (math::abs(cur.gv) <= tolC)
            {
                converged = true;
                break;
            }

            // Stall: a Newton step that did not reduce |g| below the recent evaluations
            if (byNewton && nh > 0)
            {
                T mn = hist[0];
                for (index_t i = 1; i < nh; ++i)
                    mn = math::min(mn,hist[i]);
                if (math::abs(cur.gv) >= mn && math::abs(best.gv) <= best.tolS)
                {
                    kept = best;
                    kept.tolAcc = best.tolS;
                    stagnated = true;
                    break;
                }
            }
            hist[2] = hist[1]; hist[1] = hist[0]; hist[0] = math::abs(cur.gv);
            nh = (nh < 3) ? nh + 1 : 3;

            if (first)
            {
                first = false;
                if (cur.gv > 0)
                {
                    // Taut check at gamma = 0; the iterate stays gamma_lin
                    T tC, tS;
                    if (nPass >= maxIt || !provide(T(0),n,tC,tS))
                        return false;
                    if (n.n1.dot(Sp) >= 0)
                    {
                        tautAtZero = true;
                        return true;
                    }
                    lo = 0;
                    hi = cur.g;
                    closeBracket();
                }
                else
                {
                    lo = cur.g;
                    bs = openUpper;
                }
            }
            else if (cur.gv < 0)
                lo = cur.g;
            else
            {
                hi = cur.g;
                closeBracket();
            }

            // Bracket collapsed to round-off
            if (bs == closedBracket && (hi - lo) <= 4 * eps * math::max(lo,hi) && math::abs(best.gv) <= best.tolS)
            {
                kept = best;
                kept.tolAcc = best.tolS;
                stagnated = true;
                break;
            }

            const bool newtonOK = cur.a > 0 && math::isfinite(cur.a);
            const T gn = newtonOK ? cur.g - cur.gv / cur.a : nan;
            if (bs == openUpper)
            {
                // Newton, at most doubling; a clipped step counts as a doubling
                if (newtonOK && gn <= 2 * cur.g)
                {
                    gt = gn;
                    byNewton = true;
                }
                else
                {
                    gt = 2 * cur.g;
                    byNewton = false;
                }
            }
            else
            {
                // rtsafe: bisect unless the Newton point is inside the bracket and the step before last halved
                if (newtonOK && gn > lo && gn < hi && !(math::abs(2 * cur.gv) > math::abs(dxold * cur.a)))
                {
                    dxold = dx;
                    dx    = gn - cur.g;
                    gt    = gn;
                    byNewton = true;
                }
                else
                    bisect();
            }
        }

        if (converged)
        {
            // One polishing Newton step (the point is kept only if |g| does not grow)
            kept = cur;
            kept.tolAcc = cur.tolC;
            if (cur.a > 0 && math::isfinite(cur.a) && nPass < maxIt)
            {
                const T gn = cur.g - cur.gv / cur.a;
                const bool inside = gn > 0 && (bs != closedBracket || (gn > lo && gn < hi));
                if (inside && provide(gn,n,tolC,tolS) && math::abs(n.n1.dot(Sp)) <= math::abs(cur.gv))
                {
                    kept.g = gn;
                    kept.gv = n.n1.dot(Sp);
                    kept.tolAcc = tolC;
                    kept.S = Sp;
                    kept.C = Cp;
                }
            }
        }

        const T R2 = n.n2.dot(kept.S);
        if (math::abs(R2) <= kept.tolAcc)
        {
            gamma   = kept.g;
            Sp      = kept.S;
            Cp      = kept.C;
            resNorm = math::max(math::abs(kept.gv),math::abs(R2));
            return true;
        }

        if (handOff)
        {
            // 1: the scalar pass converged at tolC; 2: it was accepted at the stall tolerance
            *handOff = stagnated ? 2 : 1;
            gamma = kept.g;
            stagnated = false;
            return false;
        }

        g = kept.g;
        stagnated = false;
    }

    // 2x2 pass
    used2x2 = true;
    nPass = 0;
    T rn1 = nan, rn2 = nan;  // ||R_{k-1}||, ||R_{k-2}||
    bool haveBest = false;

    for (;;)
    {
        if (nPass >= maxIt)
            return false;

        const VoigtDirs<T> n(math::cos(th),math::sin(th));
        T tolC, tolS;
        if (!provide(g,n,tolC,tolS))
            return false;

        const T R1 = n.n1.dot(Sp);
        const T R2 = n.n2.dot(Sp);
        cur.g = g;
        cur.th = th;
        cur.gv = math::max(math::abs(R1),math::abs(R2));
        cur.tolC = tolC;
        cur.tolS = tolS;
        cur.S = Sp;
        cur.C = Cp;
        if (!haveBest || cur.gv < best.gv)
        {
            best = cur;
            haveBest = true;
        }

        // Newton step on R = (n1.S', n2.S')
        const gsVector<T,3> Cn1 = Cp * n.n1;
        const gsVector<T,3> Cn2 = Cp * n.n2;
        const T J11 = n.n1.dot(Cn1);
        const T J12 = 2 * R2 + 2 * g * n.n1.dot(Cn2);
        const T J21 = n.n2.dot(Cn1);
        const T J22 = n.n3.dot(Sp) + 2 * g * n.n2.dot(Cn2);
        const T det = J11*J22 - J12*J21;
        const bool stepOK = math::isfinite(det) && math::abs(det) > T(1e-14) * (math::abs(J11*J22) + math::abs(J12*J21));
        T dg = nan, dth = nan;
        if (stepOK)
        {
            dg  = -(J22*R1 - J12*R2) / det;
            dth = -(J11*R2 - J21*R1) / det;
        }

        if (cur.gv <= tolC)
        {
            // One polishing step (the point is kept only if ||R||_inf does not grow)
            kept = cur;
            const T gn = g + dg, thn = th + dth;
            if (stepOK && gn > 0 && math::isfinite(thn) && nPass < maxIt)
            {
                const VoigtDirs<T> nn(math::cos(thn),math::sin(thn));
                if (provide(gn,nn,tolC,tolS))
                {
                    const T rnew = math::max(math::abs(nn.n1.dot(Sp)),math::abs(nn.n2.dot(Sp)));
                    if (rnew <= cur.gv)
                    {
                        kept.g = gn;
                        kept.th = thn;
                        kept.gv = rnew;
                        kept.S = Sp;
                        kept.C = Cp;
                    }
                }
            }
            gamma   = kept.g;
            theta   = kept.th;
            Sp      = kept.S;
            Cp      = kept.C;
            resNorm = kept.gv;
            return true;
        }

        if (!stepOK)
            return false;

        if (nPass >= 4 && cur.gv >= rn2)
        {
            if (best.gv > best.tolS)
                return false;
            stagnated = true;
            gamma   = best.g;
            theta   = best.th;
            Sp      = best.S;
            Cp      = best.C;
            resNorm = best.gv;
            return true;
        }

        g  += dg;
        th += dth;
        if (!(math::isfinite(g) && math::isfinite(th)))
            return false;

        rn2 = rn1;
        rn1 = cur.gv;
    }
}

template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if< _linear, bool>::type
gsMaterialMatrixTFTSplit<dim,T,linear>::_exactSolve(const gsVector<T,3> &, const detail_tftsplit::WrinkleOptions<T> &, const SC_provider &, const bool,
                                                    T &, T &, gsVector<T,3> &, gsMatrix<T,3,3> &, T &, index_t &,
                                                    bool &, bool &, bool &, index_t *) const
{
    GISMO_ERROR("gsMaterialMatrixTFTSplit: exact solve called for linear=true");
}

/*
    Exact wrinkled state by minimisation of Phi(v) = W(E + P(v)), P(v) = (v1^2, v2^2, 2 v1 v2) = v (x) v in
    engineering Voigt notation, over the covector v = sqrt(gamma) n in R^2 (gamma = |v|^2, n = (cos theta, sin theta)),
    with W the potential of the wrapped law. A stationary point with gamma > 0 is the wrinkled state S(E + gamma n1);
    unlike the (gamma,theta) parametrisation, v is regular at gamma = 0. The iteration is regular there;
    theta = atan2(v2,v1) itself has a condition number ~ 1/sqrt(gamma).

    Derivatives. With E' = E + P(v), S' = S(E'), C' = C(E'), Shat' = [[S'_0, S'_2], [S'_2, S'_1]] and
    B = dP/dv = [[2 v1, 0], [0, 2 v2], [2 v2, 2 v1]] (3x2):
        grad Phi = B^T S' = 2 Shat' v,
        Hes      = B^T C' B + 2 Shat'   (the second derivative of P contracted with S' is 2 Shat').

    Relations. With n = v/|v| = (c,s), m = (-s,c), Q = [n m] and R = (R1,R2) = (n1.S', n2.S'):
        n.Shat' n = R1, m.Shat' n = R2, hence grad Phi = 2 sqrt(gamma) Q R, so that
        |grad Phi| = 2 sqrt(gamma) |R| and, for v != 0, grad Phi = 0 iff R = 0;
        B n = 2 sqrt(gamma) n1, B m = 2 sqrt(gamma) n2, hence
        Q^T Hes Q = 4 gamma M + 2 [[R1, R2], [R2, sw]],
        M = [[n1.C'n1, n1.C'n2], [n2.C'n1, n2.C'n2]], sw = m.Shat' m = n4.S', n4 = n3 + n1.

    Lemma. At R = 0, gamma > 0 and symmetric C': Q^T Hes Q = Jac0 diag(4 gamma, 2), with
    Jac0 = [[M11, 2 gamma M12], [M21, sw + 2 gamma M22]] the matrix of the acceptance test in wrinkle().
    Proof: at R = 0 the second summand is 2 diag(0, sw), and Q^T Hes Q = 4 gamma M + 2 diag(0, sw) has the
    columns 4 gamma (M11, M21) and (4 gamma M12, 4 gamma M22 + 2 sw) = 2 (2 gamma M12, sw + 2 gamma M22),
    which are the columns of Jac0 scaled by 4 gamma and 2. Hence det Hes = 8 gamma det Jac0
    and (Q^T Hes Q)_11 = 4 gamma n1.C'n1: Hes is positive definite iff M11 > 0 and det Jac0 > 0, which is the
    minimum test of wrinkle(). The test does not imply n4.S' >= 0: a rank-one minimum with n4.S' < 0 is reported
    as a wrong branch by wrinkle().
    At v = 0, Hes = 2 Shat(E): a direction with n.Shat(E) n < 0 is a direction of negative curvature of Phi
    at v = 0, and S(E) positive definite makes v = 0 a strict local minimum.

    Start rule (in wrinkle(); the first match wins). (1) ExactLocalSolve 1, scalar pass converged at tolC but
    n2.S' not within it (hand-off 1): its (gamma, theta), admissible. (2) Scalar pass accepted at the stall
    tolerance (hand-off 2): its (gamma, theta), with the sign test. (3) n1.C n1 > 0 and the linearised
    gamma > 0, i.e. n1.S(E) < 0, a direction of negative curvature at v = 0 (also the direction kept after a
    failed direction search): the linearised (gamma, theta), with the sign test. (4) The pencil
    S_min >= -band max(|S_min|,|S_max|) (band = 1e-8), for every TensionField: taut at gamma = 0, the minimisation does not
    run; this includes S(E) positive definite, where v = 0 is a strict local minimum. (5) Otherwise theta0 is
    the Euclidean eigenvector of Shat(E) for its smaller eigenvalue lamS; the congruence with G preserves
    inertia, so S_min < 0 in the pencil iff lamS < 0. gamma0 = -lamS/(n1.C n1) if n1.C n1 at theta0 is positive and
    finite, else |E_min|, or max(|E_min|,|E_max|) if E_min = 0; with the sign test. If round-off gives lamS >= 0
    there, the solve is counted as not converged.

    Start (gamma0, theta0). If \a admissible, the start is used as is (it lies on the ray minimiser, which
    has Phi < W(E)); otherwise the sign test n1.S(E + gamma0 n1) < 0 must hold, and gamma0 is halved
    (at most ten times) until it does.

    Iteration. At the iterate, ||R||_inf <= tolC ends the iteration (polish below). Otherwise the Newton step
    is taken with eigenvalue flipping: with Hes = lam1 q1 q1^T + lam2 q2 q2^T (closed form, lam1 <= lam2),
    d = -sum_i q_i (q_i.grad Phi)/max(|lam_i|, 1e-8 ||Hes||), ||Hes|| = max|lam_i|; the step is a descent direction
    even where Hes is indefinite. Where the gradient is at round-off level and lam1 < 0, d = |v| q1 (signed to
    descend) leaves the saddle. With Hes = 0 and grad Phi != 0 the step is the steepest-descent step
    d = -|v| grad Phi/|grad Phi| (length |v|); with Hes = 0 and grad Phi = 0 the iterate is stationary and the
    pass ends as a stall (see Stall). The step is capped at |d| <= |v| (trust radius), which keeps v from
    crossing the origin, where theta is undefined. Where grad Phi is above round-off, the line search is Armijo
    backtracking on Phi, Phi(v + a d) <= Phi(v) + 1e-4 a grad Phi.d, a = 2^-j, j <= 10.

    Floor. Phi is only known to floorW = 64 eps (||C'||_F max(1,||E'||)^2 + |Phi|). Once the predicted decrease
    |grad Phi.d| is below floorW and Hes is positive definite, Phi cannot discriminate between trial points:
    the unit step is judged on the residual (accepted iff ||R||_inf does not grow, without an energy
    evaluation) and Armijo is used only if it is rejected.

    Polish. After convergence one Newton step is tried only if budget is left (fewer than ExactMaxIt S/C calls
    made) and Hes at the converged iterate is positive definite with a regular step (not the zero-Hessian
    case); it is kept iff ||R||_inf does not grow. Otherwise the converged iterate is returned as is.

    Stall. Three consecutive iterates that do not reduce the best ||R||_inf, a failed line search, or a
    stationary iterate with zero Hessian and zero gradient, end the pass: the best iterate is accepted
    (\a stagnated) iff its residual is within the stall tolerance.
    All tolerances follow the pattern of _exactSolve: tolC = max(ExactTol ||S'||, floorR), tolS likewise with
    ExactStagTol, floorR = 64 eps ||C'||_F max(1,||E'||).

    Budget. At most ExactMaxIt trial points (provider calls S/C) per call. A budget exit (a further S/C call is
    needed when none is left) before convergence is non-convergence; at convergence the budget only skips the polish.
    A stall exit on the last allowed trial is judged by the stall rule above.
    Energy evaluations are not counted. The budget is separate from that of the scalar pass, so a handed-off
    point can make up to 2 ExactMaxIt S/C calls in total.

    Returns true if converged or accepted at the stall tolerance. On true, gamma = |v|^2, theta = atan2(v2,v1)
    (not wrapped), Sp, Cp are S', C' and resNorm = ||R||_inf there; on false they are not modified. nEval
    counts every S/C call.

    Complexity: O(1) 2x2 algebra per trial. An iteration whose unit step is accepted costs one S/C call, plus one
    W call above the floor; one with a full backtrack up to 11 S/C and 12 W calls. Every W call is at a point that
    already had its S/C call, so a pass makes at most ExactMaxIt S/C calls and at most as many W calls. Each S/C call
    is two wrapped _C evaluations and each W call one, each with its own thickness Newton for compressible laws.

    The energy is the wrapped law's eval3D_energy_C through \a WatE; a law without one (any law that keeps the
    gsMaterialMatrixBase default) raises GISMO_NO_IMPLEMENTATION at the first evalEnergy call and needs
    ExactLocalSolve = 0.
*/
template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<!_linear, bool>::type
gsMaterialMatrixTFTSplit<dim,T,linear>::_exactMinimise(const gsVector<T,3> & E3, const detail_tftsplit::WrinkleOptions<T> & opt,
                                                       const SC_provider & SCatE, const W_provider & WatE,
                                                       const T theta0, const T gamma0, const bool admissible,
                                                       T & gamma, T & theta, gsVector<T,3> & Sp, gsMatrix<T,3,3> & Cp, T & resNorm,
                                                       index_t & nEval, bool & stagnated) const
{
    using namespace detail_tftsplit;
    const T nan = std::numeric_limits<T>::quiet_NaN();
    const T eps = std::numeric_limits<T>::epsilon();
    const index_t maxIt = opt.exactMaxIt;
    const T delta = T(1e-8), c1 = T(1e-4);

    index_t nPass = 0;

    // An evaluated trial point; R1 = n1.S', Rn = ||R||_inf, Phi valid iff haveW
    struct Trial
    {
        gsVector<T,2> v;
        gsVector<T,3> E, S;
        gsMatrix<T,3,3> C;
        T R1, Rn, tolC, tolS, Phi;
        bool haveW;
    };

    gsMatrix<T> Eel(3,1), Se, Ce;
    gsVector<T,3> Pv;

    // The only place that calls SCatE; false if S', C' or the residual is not finite
    auto evalTrial = [&](const gsVector<T,2> & v, Trial & tr) -> bool
    {
        Pv << v(0)*v(0), v(1)*v(1), 2*v(0)*v(1);   // P(v) = v (x) v, engineering Voigt
        Eel = E3 + Pv;
        SCatE(Eel,Se,Ce); ++nEval; ++nPass;

        tr.v = v;
        tr.E = Eel;
        tr.haveW = false;
        tr.Phi = nan;
        for (index_t i = 0; i != 3; ++i)
            tr.S(i) = Se.data()[i];
        if (Ce.rows() == 3 && Ce.cols() == 3)
            for (index_t l = 0; l != 3; ++l)
                for (index_t i = 0; i != 3; ++i)
                    tr.C(i,l) = Ce(i,l);
        else
            for (index_t l = 0; l != 3; ++l)
                for (index_t i = 0; i != 3; ++i)
                    tr.C(i,l) = Ce.data()[i+3*l];

        tr.R1 = tr.Rn = tr.tolC = tr.tolS = nan;
        if (!(tr.S.allFinite() && tr.C.allFinite()))
            return false;

        const T r = v.norm();
        const VoigtDirs<T> nd(v(0)/r,v(1)/r);
        tr.R1 = nd.n1.dot(tr.S);
        tr.Rn = math::max(math::abs(tr.R1),math::abs(nd.n2.dot(tr.S)));

        const T Sn     = tr.S.norm();
        const T floorR = 64 * eps * tr.C.norm() * math::max(T(1),Eel.norm());
        tr.tolC = math::max(opt.exactTol * Sn,floorR);
        tr.tolS = math::max(opt.exactStagTol * Sn,floorR);
        return math::isfinite(tr.Rn);
    };

    // Phi at a trial; false if not finite
    auto evalEnergy = [&](Trial & tr) -> bool
    {
        gsMatrix<T> Em(3,1);
        Em = tr.E;
        tr.Phi = WatE(Em);
        tr.haveW = math::isfinite(tr.Phi);
        return tr.haveW;
    };

    // Step at a trial: 0 = ok, 1 = non-finite Hessian (abort), 2 = stationary with a zero Hessian (stall)
    gsVector<T,2> gradPhi;
    bool HesPD = false;
    auto stepAt = [&](const Trial & tr, gsVector<T,2> & d) -> int
    {
        const gsVector<T,2> & v = tr.v;
        gsMatrix<T,2,2> Shat, Hes;
        gsMatrix<T,3,2> B;
        Shat << tr.S(0), tr.S(2), tr.S(2), tr.S(1);
        gradPhi = 2 * Shat * v;
        B << 2*v(0), 0, 0, 2*v(1), 2*v(1), 2*v(0);
        Hes = B.transpose() * tr.C * B + 2 * Shat;
        Hes = T(0.5) * (Hes + Hes.transpose());

        // Eigen-decomposition of [[a,q],[q,d]]: lam1 <= lam2, q1 and q2 the unit eigenvectors
        const T a = Hes(0,0), dd = Hes(1,1), q = Hes(0,1);
        const T mid = T(0.5) * (a + dd);
        const T rad = math::sqrt(T(0.25) * (a - dd) * (a - dd) + q * q);
        const T lam1 = mid - rad, lam2 = mid + rad;
        const T phi = T(0.5) * math::atan2(2*q,a - dd);
        gsVector<T,2> q2, q1;
        q2 << math::cos(phi), math::sin(phi);
        q1 << -math::sin(phi), math::cos(phi);
        const T nHes = math::max(math::abs(lam1),math::abs(lam2));
        HesPD = lam1 > 0;
        if (!math::isfinite(nHes))
            return 1;

        const T nv = v.norm();
        const T ng = gradPhi.norm();
        if (nHes == 0)
        {
            if (ng == 0)
                return 2;
            d = -(nv / ng) * gradPhi;
        }
        else if (lam1 < 0 && ng <= 64 * eps * nHes * nv)
        {
            d = nv * q1;
            if (gradPhi.dot(d) > 0)
                d = -d;
        }
        else
        {
            const T lt1 = math::max(math::abs(lam1),delta * nHes);
            const T lt2 = math::max(math::abs(lam2),delta * nHes);
            d = -(q1 * (q1.dot(gradPhi)) / lt1 + q2 * (q2.dot(gradPhi)) / lt2);
        }

        const T nd = d.norm();
        if (nd > nv)
            d *= nv / nd;
        return 0;
    };

    if (!(gamma0 > 0 && math::isfinite(gamma0)))
        return false;

    // Start: the sign test n1.S(E + gamma0 n1) < 0 unless admissible
    Trial cur, best, next;
    {
        gsVector<T,2> v;
        v << math::sqrt(gamma0) * math::cos(theta0), math::sqrt(gamma0) * math::sin(theta0);
        for (index_t k = 0; ; ++k)
        {
            if (nPass >= maxIt)
                return false;
            const bool fin = evalTrial(v,cur);
            if (fin && (admissible || cur.R1 < 0))
                break;
            if (admissible || k == 10)
                return false;
            v *= math::sqrt(T(0.5));
        }
    }
    if (!evalEnergy(cur))
        return false;
    best = cur;
    index_t noRed = 0;
    T PhiRef = cur.Phi;
    bool polish = false;

    for (;;)
    {
        if (cur.Rn <= cur.tolC)
        {
            polish = true;
            break;
        }
        if (noRed >= 3)
            break;

        gsVector<T,2> d;
        const int sc = stepAt(cur,d);
        if (sc == 1)
            return false;
        if (sc == 2)
            break;

        const T pred = gradPhi.dot(d);
        const T enorm  = math::max(T(1),cur.E.norm());
        const T floorW = 64 * eps * (cur.C.norm() * enorm * enorm + math::abs(PhiRef));

        bool accepted = false, haveReuse = false, reuseFinite = false;
        Trial reuse;
        if (math::abs(pred) <= floorW && HesPD)
        {
            if (nPass >= maxIt)
                return false;
            reuseFinite = evalTrial(cur.v + d,reuse);
            haveReuse = true;
            if (reuseFinite && reuse.Rn <= cur.Rn)
            {
                next = reuse;
                accepted = true;
            }
        }

        if (!accepted)
        {
            if (!cur.haveW)
            {
                if (!evalEnergy(cur))
                    return false;
                PhiRef = cur.Phi;
            }
            T a = 1;
            for (index_t j = 0; j != 11 && !accepted; ++j, a *= T(0.5))
            {
                Trial tt;
                bool fin;
                if (j == 0 && haveReuse)
                {
                    tt = reuse;
                    fin = reuseFinite;
                }
                else
                {
                    if (nPass >= maxIt)
                        return false;
                    fin = evalTrial(cur.v + a * d,tt);
                }
                if (fin && evalEnergy(tt) && tt.Phi <= cur.Phi + c1 * a * pred)
                {
                    next = tt;
                    accepted = true;
                }
            }
            if (!accepted)
                break;
        }

        if (next.Rn < best.Rn)
        {
            noRed = 0;
            best = next;
        }
        else
            ++noRed;
        if (next.haveW)
            PhiRef = next.Phi;
        cur = next;
    }

    const Trial * kept = &best;
    Trial polished;
    if (polish)
    {
        polished = cur;
        gsVector<T,2> d;
        if (nPass < maxIt && stepAt(cur,d) == 0 && HesPD)
        {
            Trial tt;
            if (evalTrial(cur.v + d,tt) && tt.Rn <= cur.Rn)
                polished = tt;
        }
        kept = &polished;
    }
    else
    {
        if (!(best.Rn <= best.tolS))
            return false;
        stagnated = true;
    }

    gamma   = kept->v.squaredNorm();
    theta   = math::atan2(kept->v(1),kept->v(0));
    Sp      = kept->S;
    Cp      = kept->C;
    resNorm = kept->Rn;
    return true;
}

template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if< _linear, bool>::type
gsMaterialMatrixTFTSplit<dim,T,linear>::_exactMinimise(const gsVector<T,3> &, const detail_tftsplit::WrinkleOptions<T> &,
                                                       const SC_provider &, const W_provider &,
                                                       const T, const T, const bool,
                                                       T &, T &, gsVector<T,3> &, gsMatrix<T,3,3> &, T &,
                                                       index_t &, bool &) const
{
    GISMO_ERROR("gsMaterialMatrixTFTSplit: exact minimisation called for linear=true");
}

template <short_t dim, class T, bool linear >
std::vector<detail_tftsplit::WrinkleResult<T>> gsMaterialMatrixTFTSplit<dim,T,linear>::_evalCore(
                                                const index_t patch,
                                                const gsMatrix<T> & u,
                                                const gsMatrix<T> & z,
                                                bool needTangent,
                                                bool needEnergy) const
{
    const detail_tftsplit::WrinkleOptions<T> opt = this->wrinkleOptions();

    this->_computePoints(patch,u);

    std::vector<detail_tftsplit::WrinkleResult<T>> res(u.cols()*z.rows());
    gsMatrix<T> Gori(2,2), Gdef(2,2);
    gsMatrix<T> Gb(3,2), gb(3,2);
    gsMatrix<T> Ev(3,1);
    for (index_t k = 0; k!=u.cols(); ++k)
    {
        const gsVector<T> uk = u.col(k);
        for (index_t j = 0; j!=z.rows(); ++j)
        {
            const index_t colIdx = j*u.cols()+k;
            const T zz = z(j,k);

            this->_getMetric(k, zz * m_data.mine().m_Tmat(0,k));
            if (linear)
            {
                Gori = m_data.mine().m_Acov_ori;
                Gdef = m_data.mine().m_Acov_def;
            }
            else
            {
                Gori = m_data.mine().m_Gcov_ori.block(0,0,2,2);
                Gdef = m_data.mine().m_Gcov_def.block(0,0,2,2);
            }
            Ev(0,0) = T(0.5) * (Gdef(0,0) - Gori(0,0));
            Ev(1,0) = T(0.5) * (Gdef(1,1) - Gori(1,1));
            Ev(2,0) =          (Gdef(0,1) - Gori(0,1));

            // Cartesian covariant basis vectors (columns): mid-surface for linear, fibre otherwise
            if (linear && dim==3)
            {
                Gb = m_data.mine().m_acov_ori;
                gb = m_data.mine().m_acov_def;
            }
            else
            {
                Gb = m_data.mine().m_gcov_ori.leftCols(2);
                gb = m_data.mine().m_gcov_def.leftCols(2);
            }

            const gsMatrix<T> S = m_materialMat->eval3D_vector(patch,uk,zz,MaterialOutput::VectorN);
            const gsMatrix<T> C = m_materialMat->eval3D_matrix(patch,uk,zz,MaterialOutput::MatrixA);

            // dC/dE_d from the wrapped material: eval3D_dmatrix has entry d + 3*(i + 3*l) = dC_il/dC_d
            // and dC/dE_d = 2 dC/dC_d.
            dC_provider provider;
            if (!linear)
                provider = [&](std::vector<gsMatrix<T>> & dC)
                {
                    const gsMatrix<T> v = m_materialMat->eval3D_dmatrix(patch,uk,zz,MaterialOutput::MatrixA);
                    dC.assign(3,gsMatrix<T>(3,3));
                    for (index_t d = 0; d!=3; ++d)
                        for (index_t l = 0; l!=3; ++l)
                            for (index_t i = 0; i!=3; ++i)
                                dC[d](i,l) = 2 * v((i+3*l)*3 + d);
                };

            SC_provider scProvider;
            if (!linear && opt.exact)
                scProvider = [&](const gsMatrix<T> & Eel, gsMatrix<T> & Se, gsMatrix<T> & Ce)
                {
                    gsMatrix<T> Cm(3,1);
                    Cm(0,0) = Gori(0,0) + 2 * Eel(0,0);
                    Cm(1,0) = Gori(1,1) + 2 * Eel(1,0);
                    Cm(2,0) = Gori(0,1) +     Eel(2,0);
                    Se = m_materialMat->eval3D_vector_C(Cm,patch,uk,zz,MaterialOutput::VectorN);
                    Ce = m_materialMat->eval3D_matrix_C(Cm,patch,uk,zz,MaterialOutput::MatrixA);
                };

            W_provider wProvider;
            if (!linear && (needEnergy || opt.exact))
                wProvider = [&](const gsMatrix<T> & Eel) -> T
                {
                    gsMatrix<T> Cm(3,1);
                    Cm(0,0) = Gori(0,0) + 2 * Eel(0,0);
                    Cm(1,0) = Gori(1,1) + 2 * Eel(1,0);
                    Cm(2,0) = Gori(0,1) +     Eel(2,0);
                    return m_materialMat->eval3D_energy_C(Cm,patch,uk,zz,MaterialOutput::Generic)(0,0);
                };

            res[colIdx] = this->wrinkle(S,C,Ev,Gori,opt,needTangent,provider,scProvider,wProvider,needEnergy);
            res[colIdx].gdef = Gdef;
            res[colIdx].Gbasis = Gb;
            res[colIdx].gbasis = gb;
        }
    }
    return res;
}

/*
    Prescribed-metric path. The wrapped material provides S(Cmat) and C(Cmat), Cmat = (g11,g22,g12)
    with the tensor g12. The derivative of C with respect to the engineering strain is taken by central
    differences of eval3D_matrix_C in x = (g11, g22, 2 g12), with step h on x. A step h on x_2 is a step
    h/2 on Cmat(2). Since E_1 = (g11 - G11)/2 and E_3 = 2E_12 = g12 - G12 (so dE = dx/2 for each
    component), dC/dE_r = 2 dC/dx_r = (C(x + h e_r) - C(x - h e_r)) / h.
    Error: O(h^2 |d^4 W|) + (Newton tolerance 1e-10 of the wrapped C33 solve)/h, about 1e-5 relative.
*/
template <short_t dim, class T, bool linear >
detail_tftsplit::WrinkleResult<T> gsMaterialMatrixTFTSplit<dim,T,linear>::_evalCore_C(
                                                const gsMatrix<T> & Cmat,
                                                const index_t patch,
                                                const gsVector<T> & u,
                                                const T z,
                                                bool needTangent,
                                                bool needEnergy) const
{
    const detail_tftsplit::WrinkleOptions<T> opt = this->wrinkleOptions();

    const gsMatrix<T> um = u;
    this->_computePoints(patch,um);
    this->_getMetricUndeformed(0, z * m_data.mine().m_Tmat(0,0));

    gsMatrix<T> Gori(2,2);
    if (linear)
        Gori = m_data.mine().m_Acov_ori;
    else
        Gori = m_data.mine().m_Gcov_ori.block(0,0,2,2);

    gsMatrix<T> Ev(3,1);
    Ev(0,0) = T(0.5) * (Cmat(0,0) - Gori(0,0));
    Ev(1,0) = T(0.5) * (Cmat(1,0) - Gori(1,1));
    Ev(2,0) =          (Cmat(2,0) - Gori(0,1));

    const gsMatrix<T> S = m_materialMat->eval3D_vector_C(Cmat,patch,u,z,MaterialOutput::VectorN);
    const gsMatrix<T> C = m_materialMat->eval3D_matrix_C(Cmat,patch,u,z,MaterialOutput::MatrixA);

    dC_provider provider;
    if (!linear)
        provider = [&](std::vector<gsMatrix<T>> & dC)
        {
            const T h = 1e-5;
            dC.assign(3,gsMatrix<T>(3,3));
            for (index_t r = 0; r!=3; ++r)
            {
                const T step = (r < 2) ? h : h/2;
                gsMatrix<T> Cp = Cmat, Cm = Cmat;
                Cp(r,0) += step;
                Cm(r,0) -= step;
                const gsMatrix<T> Ap = m_materialMat->eval3D_matrix_C(Cp,patch,u,z,MaterialOutput::MatrixA);
                const gsMatrix<T> Am = m_materialMat->eval3D_matrix_C(Cm,patch,u,z,MaterialOutput::MatrixA);
                for (index_t l = 0; l!=3; ++l)
                    for (index_t i = 0; i!=3; ++i)
                        dC[r](i,l) = ( Ap(i+3*l,0) - Am(i+3*l,0) ) / h;
            }
        };

    SC_provider scProvider;
    if (!linear && opt.exact)
        scProvider = [&](const gsMatrix<T> & Eel, gsMatrix<T> & Se, gsMatrix<T> & Ce)
        {
            gsMatrix<T> Cm(3,1);
            Cm(0,0) = Gori(0,0) + 2 * Eel(0,0);
            Cm(1,0) = Gori(1,1) + 2 * Eel(1,0);
            Cm(2,0) = Gori(0,1) +     Eel(2,0);
            Se = m_materialMat->eval3D_vector_C(Cm,patch,u,z,MaterialOutput::VectorN);
            Ce = m_materialMat->eval3D_matrix_C(Cm,patch,u,z,MaterialOutput::MatrixA);
        };

    W_provider wProvider;
    if (!linear && (needEnergy || opt.exact))
        wProvider = [&](const gsMatrix<T> & Eel) -> T
        {
            gsMatrix<T> Cm(3,1);
            Cm(0,0) = Gori(0,0) + 2 * Eel(0,0);
            Cm(1,0) = Gori(1,1) + 2 * Eel(1,0);
            Cm(2,0) = Gori(0,1) +     Eel(2,0);
            return m_materialMat->eval3D_energy_C(Cm,patch,u,z,MaterialOutput::Generic)(0,0);
        };

    return this->wrinkle(S,C,Ev,Gori,opt,needTangent,provider,scProvider,wProvider,needEnergy);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::_collectStress(const std::vector<detail_tftsplit::WrinkleResult<T>> & res) const
{
    gsMatrix<T> result(3,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result.col(c) = res[c].S;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::_collectCauchy(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result = this->_collectStress(res);
    const gsMatrix<T> detF = m_materialMat->eval3D_detF(patch,u,z,out);
    for (size_t c = 0; c!=res.size(); ++c)
        result.col(c) /= detF(0,c);
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,true);
    gsMatrix<T> result(9,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        for (index_t l = 0; l!=3; ++l)
            for (index_t i = 0; i!=3; ++i)
                result(i+3*l,c) = res[c].C(i,l);
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    if (m_options.getSwitch("Midsurface") && out==MaterialOutput::VectorM)
        GISMO_ERROR("gsMaterialMatrixTFTSplit: Midsurface=true supports membrane analysis only (VectorM undefined)");
    return this->_collectStress(this->_evalCore(patch,u,z,false));
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput /*out*/) const
{
    const detail_tftsplit::WrinkleResult<T> r = this->_evalCore_C(Cmat,patch,u,z,true);
    gsMatrix<T> result(9,1);
    for (index_t l = 0; l!=3; ++l)
        for (index_t i = 0; i!=3; ++i)
            result(i+3*l,0) = r.C(i,l);
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_vector_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput /*out*/) const
{
    const detail_tftsplit::WrinkleResult<T> r = this->_evalCore_C(Cmat,patch,u,z,false);
    return r.S;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    return this->_collectStress(this->_evalCore(patch,u,z,false));
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_strain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(3,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result.col(c) = res[c].E;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(2,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
    {
        const detail_tftsplit::Mat2<T> M = detail_tftsplit::stressPencil(res[c].Gcov,res[c].S(0,0),res[c].S(1,0),res[c].S(2,0));
        detail_tftsplit::principalValues(M,result(0,c),result(1,c));
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_pstrain(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(3,res.size());
    result.setZero();
    for (size_t c = 0; c!=res.size(); ++c)
    {
        const detail_tftsplit::Mat2<T> Gcon = detail_tftsplit::invert2x2(res[c].Gcov);
        const detail_tftsplit::Mat2<T> M = detail_tftsplit::strainPencil(Gcon,res[c].E(0,0),res[c].E(1,0),res[c].E(2,0));
        detail_tftsplit::principalValues(M,result(0,c),result(1,c));
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return this->_collectCauchy(patch,u,z,out);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_CauchyVector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    return this->_collectCauchy(patch,u,z,out);
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    const gsMatrix<T> detF = m_materialMat->eval3D_detF(patch,u,z,out);
    gsMatrix<T> result(2,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
    {
        detail_tftsplit::Mat2<T> M = detail_tftsplit::stressPencil(res[c].gdef,res[c].S(0,0),res[c].S(1,0),res[c].S(2,0));
        const T Jinv = 1 / detF(0,c);
        M.a *= Jinv; M.b *= Jinv; M.c *= Jinv; M.d *= Jinv;
        detail_tftsplit::principalValues(M,result(0,c),result(1,c));
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_tensionfield(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const T S2crit = m_options.getReal("S2crit");
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(1,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
    {
        if (res[c].state==0)
            result(0,c) = (res[c].Smin < S2crit) ? 0 : 1;
        else
            result(0,c) = res[c].state;
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_theta(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(1,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result(0,c) = res[c].theta;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_gamma(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(1,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result(0,c) = res[c].gamma;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_gammaHat(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(1,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result(0,c) = res[c].gammaHat;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_energy(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false,true);
    gsMatrix<T> result(1,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result(0,c) = res[c].W;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_energy_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput /*out*/) const
{
    gsMatrix<T> result(1,1);
    result(0,0) = this->_evalCore_C(Cmat,patch,u,z,false,true).W;
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_wrinkleIntensity(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(1,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
        result(0,c) = (res[c].state==0) ? detail_tftsplit::wrinkleIntensity(res[c].Emin) : T(0);
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_wrinkleDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(3,res.size());
    for (size_t c = 0; c!=res.size(); ++c)
    {
        if (res[c].state==0 && math::isfinite(res[c].theta))
            result.col(c) = detail_tftsplit::wrinkleDirection(res[c].theta,res[c].gbasis);
        else
            result.col(c).setZero();
    }
    return result;
}

template <short_t dim, class T, bool linear >
gsMatrix<T> gsMaterialMatrixTFTSplit<dim,T,linear>::eval3D_pstressDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    const std::vector<detail_tftsplit::WrinkleResult<T>> res = this->_evalCore(patch,u,z,false);
    gsMatrix<T> result(9,res.size());
    gsMatrix<T> F;
    for (size_t c = 0; c!=res.size(); ++c)
    {
        detail_tftsplit::principalStressFrame(res[c].S,res[c].Gbasis,F);
        result.col(c) = F.reshape(9,1);
    }
    return result;
}

} // namespace gismo
