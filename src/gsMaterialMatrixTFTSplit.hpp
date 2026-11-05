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
m_nExactFallbacks(0),
m_nExactReclassified(0)
{
    if (m_materialMat->hasUndeformed())
        this->setUndeformed(m_materialMat->getUndeformed());
    if (m_materialMat->hasDeformed())
        this->setDeformed(m_materialMat->getDeformed());

    this->setThickness(m_materialMat->getThickness());

    // Root finder for the fallback; it does not own the wrapped material, which m_materialMat keeps alive
    m_tft = memory::make_shared(new gsMaterialMatrixTFT<dim,T,linear>(m_materialMat.get()));

    m_options.addReal("Eta","Regularisation: fraction of the unmodified stress/tangent kept in wrinkled and slack states (Zhang 2026 Eq. 40)",0);
    m_options.addReal("SlackMultiplier","Slack state: S_eff = SlackMultiplier*S if > 0, else Eta*S",0);
    m_options.addReal("S2crit","Display-only threshold: a wrinkled point is reported wrinkled in eval3D_tensionfield iff Smin < S2crit",0);
    m_options.addSwitch("Midsurface","Evaluate the state on the mid-surface only (MatIntegration::Constant); membrane analysis only",false);
    m_options.addSwitch("Exact","Exact (non-truncated) wrinkled stress S(E + gamma n x n) by a local Newton solve (linear=false); no effect for linear=true, where both coincide",false);
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
    GISMO_ENSURE(opt.tensionField==0 || opt.tensionField==1 || opt.tensionField==2,
                 "gsMaterialMatrixTFTSplit: TensionField must be 0, 1 or 2, got "<<opt.tensionField);
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
    the stationarity of the relaxed energy W_r(E) = min_{gamma,th} W(E + gamma n x n) (Pipkin). S(E') and
    C(E') are evaluated by the wrapped material at the prescribed metric
        Cmat' = (G11 + 2E'_0, G22 + 2E'_1, G12 + E'_2)  (= Cmat + 2 gamma (c^2, s^2, cs)),
    with the tensor g12 (E'_2 is the engineering shear). Newton, with dE'/dgamma = n1, dE'/dth = 2 gamma n2:
        J11 = n1.C'n1,                      J12 = 2 n2.S' + 2 gamma n1.C'n2,
        J21 = n2.C'n1,                      J22 = n3.S' + 2 gamma n2.C'n2.
    With W(gamma,th) = W(E + gamma n1(th)): dW/dgamma = n1.S(E'), dW/dth = 2 gamma n2.S(E'), so that
    at R = 0 the Hessian of W(gamma,th) is diag(1, 2 gamma) Jac0.
    For a material isotropic w.r.t. G, S and C are coaxial, R2 vanishes along gamma at the closed-form th
    and a scalar Newton step -R1/J11 in gamma suffices; the 2x2 system is used if the guard rejected the
    closed form or if |R2| > tolR at the scalar solution.

    Convergence: tolR = 1e-8 |S|_2 (S at E), at most 20 iterations per mode. A point is accepted when
    ||R||_inf <= tolR and either the step is below 1e-10 (relative in gamma) or the residual no longer
    contracts by half (stagnation at the noise floor: a compressible wrapped law solves its thickness
    stretch to 1e-10, so S' and C' carry that noise). Non-finite S' or C', a non-positive J11 (scalar), a
    singular J (2x2), or a residual that does not decrease over two iterations after the third ends the
    solve unconverged.

    Tangent. By the envelope theorem dW_r/dE = S(E'), hence the tangent of the relaxed energy is symmetric
    whenever C' is. Differentiating R(gamma(E),th(E);E) = 0 gives d(gamma,th)/dE = -Jac^{-1} K C' with
    K = [n1^T; n2^T], and dS'/dE = C'(I + [n1, 2 gamma n2] d(gamma,th)/dE). At R = 0, n2.S' = 0 and
    n3.S' = n4.S' =: sw, n4 = n3 + n1 = (s^2,c^2,-2cs) (m x m), the tension along the wrinkles. With
        M = K C' K^T,  Jac0 = [[M11, 2 gamma M12], [M21, sw + 2 gamma M22]],  X = diag(1,2 gamma) Jac0^{-1},
        J_ex = C' - [C'n1, C'n2] X [n1^T C'; n2^T C'].
    For gamma > 0 this equals C' - C'K^T (K C' K^T + diag(0, sw/(2 gamma)))^{-1} K C'; the unscaled form is
    finite for gamma -> 0. The stationary point is a minimum iff M11 > 0 and det Jac0 > 0.
    This is guaranteed if sw > 0 and C' is positive definite on span{n1, n2}.
    In the major principal direction stationarity gives gamma < 0 (to first order -Smax/(n1.C n1)), which
    is not an admissible wrinkle; hence the start from the minor direction. The linearised gamma is the
    first scalar Newton step in gamma from gamma = 0 at fixed th (not so for the 2x2 iteration at a
    non-coaxial th, where R2 != 0 at gamma = 0). In every mode the tangent jumps at the taut/wrinkled
    boundary (W_r is C^1, not C^2).

    Complexity: two wrapped _C evaluations per Newton iteration (each with an inner thickness Newton for
    compressible laws), typically 2-5 iterations from the linearised start.
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
                                                const SC_provider & SCatE) const
{
    using namespace detail_tftsplit;
    GISMO_ASSERT(S.rows()*S.cols()==3,"S must have 3 entries");
    GISMO_ASSERT(E.rows()*E.cols()==3,"E must have 3 entries");
    GISMO_ASSERT(Gcov.rows()==2 && Gcov.cols()==2,"Gcov must be 2x2");
    GISMO_ASSERT((C.rows()==3 && C.cols()==3) || (C.rows()*C.cols()==9),"C must be 3x3 or 9x1");

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
    res.guardFallback = res.tftFailure = false;
    res.exactFallback = res.exactReclassified = false;
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
    if (!(math::isfinite(res.Smin) && math::isfinite(res.Smax) && math::isfinite(res.Emin) && math::isfinite(res.Emax)))
    {
        taut  = true;
        slack = false;
    }

    if (taut)
    {
        res.state = 1;
        res.S = S3;
        res.E = E3;
        if (needTangent)
            res.C = C3;
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

    if (!linear && opt.exact && !indefinite && !tftFailed)
    {
        if (!SCatE)
            GISMO_ERROR("gsMaterialMatrixTFTSplit::wrinkle: exact mode requires an S/C provider");

        T gE = gam, thE = theta, resNorm = nan;
        gsVector<T,3> Spx;
        gsMatrix<T,3,3> Cpx;
        index_t nEval = 0;
        const T tolR = T(1e-8) * S3.norm();
        bool exactFailed = !this->template _exactSolve<linear>(E3,tolR,SCatE,cfOK,gE,thE,Spx,Cpx,resNorm,nEval);
        res.exactIterations = nEval;

        if (!exactFailed)
        {
            thE = wrapAngle(thE);
            const T ce = math::cos(thE), se = math::sin(thE);
            const VoigtDirs<T> ne(ce,se);
            const gsVector<T,3> ne4 = ne.n3 + ne.n1;
            const T sw = ne4.dot(Spx);

            if (gE <= 0 || !(sw > 0))
            {
                const T sigma = opt.slackMultiplier > 0 ? opt.slackMultiplier : opt.eta;
                const T fac   = (gE <= 0) ? T(1) : sigma;
                res.state = (gE <= 0) ? 1 : -1;
                res.theta = res.gamma = res.gammaHat = res.residual = nan;
                res.exactReclassified = true;
#               pragma omp atomic
                ++m_nExactReclassified;
                res.S = fac * S3;
                res.E = fac * E3;
                if (needTangent)
                    res.C = fac * C3;
                return res;
            }

            const gsVector<T,3> Cn1 = Cpx * ne.n1;
            const gsVector<T,3> Cn2 = Cpx * ne.n2;
            const T M11 = ne.n1.dot(Cn1), M12 = ne.n1.dot(Cn2);
            const T M21 = ne.n2.dot(Cn1), M22 = ne.n2.dot(Cn2);
            const T a = M11, bb = 2 * gE * M12, cc = M21, d = sw + 2 * gE * M22;
            const T det0 = a*d - bb*cc;

            if (math::isfinite(M11) && math::isfinite(det0) && M11 > 0 && det0 > 0)
            {
                res.theta    = thE;
                res.gamma    = gE;
                res.gammaHat = gE * (ce*(Gcon.a*ce + Gcon.b*se) + se*(Gcon.c*ce + Gcon.d*se));
                res.residual = resNorm;

                res.S = (1 - opt.eta) * Spx + opt.eta * S3;
                res.E = E3 + gE * ne.n1;
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
            exactFailed = true;
        }

        if (exactFailed)
        {
            res.exactFallback = true;
#           pragma omp atomic
            ++m_nExactFallbacks;

            bool warn = false;
#           pragma omp critical (gsMaterialMatrixTFTSplit_exactFallback)
            {
                static bool warned = false;
                warn = !warned;
                warned = true;
            }
            if (warn)
                gsWarn<<"gsMaterialMatrixTFTSplit: the exact wrinkling solve did not converge (or reached a saddle); "
                        "the linearised result is used (reported once, see numExactFallbacks())\n";
        }
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
    Newton iteration of the exact wrinkled state; see the theory comment of wrinkle(). The scalar mode
    updates gamma only (theta fixed); the 2x2 mode updates (gamma,theta). If the scalar solution has
    |R2| > tolR the iteration continues in 2x2 mode from that point. gamma and theta are not modified
    unless the solve converges. S', C' and ||R||_inf are those at the returned iterate.
*/
template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if<!_linear, bool>::type
gsMaterialMatrixTFTSplit<dim,T,linear>::_exactSolve(const gsVector<T,3> & E3, const T tolR, const SC_provider & SCatE, const bool scalar0,
                                                    T & gamma, T & theta, gsVector<T,3> & Sp, gsMatrix<T,3,3> & Cp, T & resNorm, index_t & nEval) const
{
    using namespace detail_tftsplit;
    const index_t maxIt = 20;
    const T epsG = T(1e-12);
    const T nan  = std::numeric_limits<T>::quiet_NaN();

    T g = gamma, th = theta;
    gsMatrix<T> Eel(3,1), Se, Ce;

    for (int pass = 0; pass != 2; ++pass)
    {
        const bool scalar = (pass == 0) && scalar0;
        bool converged = false;
        T rn1 = nan, rn2 = nan;  // ||R_{k-1}||, ||R_{k-2}||
        T R2 = nan;

        for (index_t k = 0; k != maxIt; ++k)
        {
            const T c = math::cos(th), s = math::sin(th);
            const VoigtDirs<T> n(c,s);
            Eel = E3 + g * n.n1;
            SCatE(Eel,Se,Ce);
            ++nEval;

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

            const T R1 = n.n1.dot(Sp);
            R2 = n.n2.dot(Sp);
            const T rn = scalar ? math::abs(R1) : math::max(math::abs(R1),math::abs(R2));

            const gsVector<T,3> Cn1 = Cp * n.n1;
            const gsVector<T,3> Cn2 = Cp * n.n2;
            const T J11 = n.n1.dot(Cn1);

            if (rn <= tolR && k >= 1 && rn > T(0.5) * rn1)
            {
                resNorm = math::max(math::abs(R1),math::abs(R2));
                converged = true;
                break;
            }

            T dg, dth = 0;
            if (scalar)
            {
                if (!math::isfinite(J11) || !(J11 > 0))
                    return false;
                dg = -R1 / J11;
            }
            else
            {
                const T J12 = 2 * n.n2.dot(Sp) + 2 * g * n.n1.dot(Cn2);
                const T J21 = n.n2.dot(Cn1);
                const T J22 = n.n3.dot(Sp) + 2 * g * n.n2.dot(Cn2);
                const T det = J11*J22 - J12*J21;
                if (!math::isfinite(det) || !(math::abs(det) > T(1e-14) * (math::abs(J11*J22) + math::abs(J12*J21))))
                    return false;
                dg  = -(J22*R1 - J12*R2) / det;
                dth = -(J11*R2 - J21*R1) / det;
            }

            if (rn <= tolR && math::abs(dg) <= T(1e-10) * math::max(math::abs(g),epsG) && math::abs(dth) <= T(1e-10))
            {
                resNorm = math::max(math::abs(R1),math::abs(R2));
                converged = true;
                break;
            }

            if (rn > tolR && k >= 3 && rn >= rn2)
                return false;

            g  += dg;
            th += dth;
            if (!(math::isfinite(g) && math::isfinite(th)))
                return false;

            rn2 = rn1;
            rn1 = rn;
        }

        if (!converged)
            return false;
        if (!scalar || math::abs(R2) <= tolR)
        {
            gamma = g;
            theta = th;
            return true;
        }
    }
    return false;
}

template <short_t dim, class T, bool linear >
template <bool _linear>
typename std::enable_if< _linear, bool>::type
gsMaterialMatrixTFTSplit<dim,T,linear>::_exactSolve(const gsVector<T,3> &, const T, const SC_provider &, const bool,
                                                    T &, T &, gsVector<T,3> &, gsMatrix<T,3,3> &, T &, index_t &) const
{
    GISMO_ERROR("gsMaterialMatrixTFTSplit: exact solve called for linear=true");
}

template <short_t dim, class T, bool linear >
std::vector<detail_tftsplit::WrinkleResult<T>> gsMaterialMatrixTFTSplit<dim,T,linear>::_evalCore(
                                                const index_t patch,
                                                const gsMatrix<T> & u,
                                                const gsMatrix<T> & z,
                                                bool needTangent) const
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

            res[colIdx] = this->wrinkle(S,C,Ev,Gori,opt,needTangent,provider,scProvider);
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
                                                bool needTangent) const
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

    return this->wrinkle(S,C,Ev,Gori,opt,needTangent,provider,scProvider);
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
