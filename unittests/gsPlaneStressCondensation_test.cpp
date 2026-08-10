/** @file gsPlaneStressCondensation_test.cpp

    @brief Unit tests for gsPlaneStressCondensation (task 12).

    All references here are ANALYTIC (closed-form isotropic elasticity / manufactured
    strain states / central finite differences), never the code's own output re-pasted
    as truth. The suite is geometry-free: it de-risks the constitutive core consumed by
    the frame-transform work in tasks 14/15.

    == Analytic references (isotropic, Lamé  λ = Eν/((1+ν)(1−2ν)),  μ = E/(2(1+ν))) ==
      - Linear / St.Venant–Kirchhoff plane stress, S33 = 0:
            S33 = λ(E11+E22+E33) + 2μ E33 = 0  ⇒  E33 = −λ(E11+E22)/(λ+2μ).
      - Condensed SvK tangent, ENGINEERING Voigt [11,22,12] (shear input = 2E12):
            C2D = E/(1−ν²) · [[1, ν, 0],
                              [ν, 1, 0],
                              [0, 0, (1−ν)/2]]   (slot (2,2) = μ, so S12 = μ·2E12 = 2μE12).
      - S2D (linear law) = C2D · e2D.
      - Any compressible hyperelastic linearised at F = I reduces to isotropic linear
        elasticity ⇒ the NH condensed tangent at zero strain equals the SvK matrix above.

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsPlaneStressCondensation.h>
#include <gsPhaseFieldFracture/materials/gsLinearMaterial.h>
#include <gsPhaseFieldFracture/materials/gsNeoHookeMaterial.h>

SUITE(gsPlaneStressCondensation)
{

// ---- Small local helpers (analytic / FD, independent of the class internals) ----

// Constant (E,nu) parameter rows, one row of length N per parameter [E, nu].
std::vector<gsMatrix<real_t>> makeParams(index_t N, real_t E, real_t nu)
{
    std::vector<gsMatrix<real_t>> p(2);
    p[0].setConstant(1, N, E);
    p[1].setConstant(1, N, nu);
    return p;
}

// Analytic plane-stress SvK tangent (engineering Voigt [11,22,12]).
gsMatrix<real_t> svkC2D(real_t E, real_t nu)
{
    gsMatrix<real_t> C(3,3);
    const real_t f = E / (1.0 - nu*nu);
    C << f,      f*nu,   0.0,
         f*nu,   f,      0.0,
         0.0,    0.0,    f*(1.0-nu)/2.0;
    return C;
}

// Analytic out-of-plane strain E33 for the linear/SvK law.
real_t svkE33(real_t E, real_t nu, real_t E11, real_t E22)
{
    const real_t lambda = E*nu / ((1.0+nu)*(1.0-2.0*nu));
    const real_t mu     = E / (2.0*(1.0+nu));
    return -lambda*(E11+E22) / (lambda + 2.0*mu);
}

// Deterministic in-plane strain states (loop-generated, NO rand()). Engineering shear.
gsMatrix<real_t> makeStates(real_t shearMax = 0.03)
{
    const real_t vals[4] = {-0.05, 0.0, 0.02, 0.1};
    const real_t shear[2] = {0.0, shearMax};
    std::vector<gsVector<real_t,3> > states;
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
        for (int s=0; s<2; ++s)
        {
            gsVector<real_t,3> e; e << vals[i], vals[j], shear[s];
            states.push_back(e);
        }
    gsMatrix<real_t> e2D(3, states.size());
    for (size_t k=0; k!=states.size(); ++k) e2D.col(k) = states[k];
    return e2D;
}

// Build hand-filled material data (strain + F from C=2E+I Cholesky) for a direct law
// query at prescribed in-plane strain e2D and out-of-plane E33 — mirrors _fillPoint.
gsMaterialData<real_t> makeData(const gsMatrix<real_t> & e2D,
                                const gsMatrix<real_t> & E33,
                                const std::vector<gsMatrix<real_t>> & params)
{
    const index_t N = e2D.cols();
    gsMaterialData<real_t> data;
    data.dim = 3; data.size = N; data.patch = 0;
    data.parameters = params;
    data.strain.resize(9,N);
    data.deformationGradient.resize(9,N);
    for (index_t k=0; k!=N; ++k)
    {
        // Native Eigen fixed-size types: the G+Smo gsMatrix<T,3,3> wrapper has no
        // working Eigen evaluator for the TriangularView (matrixL) assignment path.
        gsEigen::Matrix<real_t,3,3> E; E.setZero();
        E(0,0) = e2D(0,k);
        E(1,1) = e2D(1,k);
        E(2,2) = E33(0,k);
        E(0,1) = E(1,0) = 0.5*e2D(2,k);
        gsEigen::Matrix<real_t,3,3> C = 2.0*E;
        C(0,0) += 1.0; C(1,1) += 1.0; C(2,2) += 1.0;
        gsEigen::LLT<gsEigen::Matrix<real_t,3,3> > llt(C);
        data.strain.reshapeCol(k,3,3) = E;
        data.deformationGradient.reshapeCol(k,3,3) = llt.matrixL().transpose();
    }
    return data;
}

// Relative Frobenius error, robust for near-zero references.
real_t relFro(const gsMatrix<real_t> & a, const gsMatrix<real_t> & ref)
{
    return (a - ref).norm() / (1.0 + ref.norm());
}

/////////////////////////////////////////////////////////////////////////////////////

// 1. E33, condensed C2D and S2D vs the closed-form isotropic plane-stress references.
TEST(SvK_analytic_E33_and_C)
{
    const real_t E = 200.0, nu = 0.3;   // O(1) magnitudes; tolerances are relative.
    const gsMatrix<real_t> e2D = makeStates();
    const index_t N = e2D.cols();
    const gsMatrix<real_t> Cref = svkC2D(E,nu);

    gsLinearMaterial<real_t> law(E, nu, 3);
    gsPlaneStressCondensation<real_t> cond(&law);

    gsMatrix<real_t> S2D, C2D, E33;
    cond.condense(e2D, makeParams(N,E,nu), S2D, C2D, E33);

    real_t maxE33 = 0, maxC = 0, maxS = 0;
    for (index_t k=0; k!=N; ++k)
    {
        // (a) E33 vs analytic.
        const real_t e33ref = svkE33(E, nu, e2D(0,k), e2D(1,k));
        const real_t eE33 = math::abs(E33(0,k)-e33ref) / (1.0 + math::abs(e33ref));
        maxE33 = math::max(maxE33, eE33);
        CHECK_CLOSE(e33ref, E33(0,k), 1e-12*(1.0+math::abs(e33ref)));

        // (b) condensed C2D vs analytic (per point; the linear tangent is constant).
        gsMatrix<real_t> Ck = C2D.reshapeCol(k,3,3);
        const real_t eC = relFro(Ck, Cref);
        maxC = math::max(maxC, eC);
        CHECK(eC < 1e-12);

        // (c) S2D vs C2D_analytic · e2D.
        gsMatrix<real_t> Sref = Cref * e2D.col(k);
        gsMatrix<real_t> Sk = S2D.col(k);
        const real_t eS = relFro(Sk, Sref);
        maxS = math::max(maxS, eS);
        CHECK(eS < 1e-12);
    }
    gsInfo<<"[SvK_analytic] max rel err  E33="<<maxE33<<"  C2D="<<maxC<<"  S2D="<<maxS<<"\n";
}

// 2. A linear (affine-in-E33) law converges in exactly ONE Newton update (task-12 property).
TEST(SvK_one_step)
{
    const real_t E = 200.0, nu = 0.3;
    const gsMatrix<real_t> e2D = makeStates();
    const index_t N = e2D.cols();

    gsLinearMaterial<real_t> law(E, nu, 3);
    gsPlaneStressCondensation<real_t> cond(&law);

    gsMatrix<real_t> S2D, C2D, E33;
    index_t iters = -1;
    cond.condense(e2D, makeParams(N,E,nu), S2D, C2D, E33, nullptr, &iters);
    gsInfo<<"[SvK_one_step] itersOut="<<iters<<"\n";
    CHECK_EQUAL(1, iters);
}

// 3. Neo-Hooke: S33≈0 at the converged state; condensed tangent and dE33dE match central FD.
TEST(NH_S33_zero_and_FD_tangent)
{
    const real_t E = 200.0, nu = 0.3;
    const gsMatrix<real_t> e2D = makeStates();   // states up to |E11|=0.1 (moderate)
    const index_t N = e2D.cols();

    gsNeoHookeMaterial<real_t> law(E, nu, 3);
    // Tight Newton tol (on |dE33|): the plane-stress residual is bounded by
    // |S33| <= C3333 * tol, so to verify S33 ~ 0 at 1e-9 the constraint must be
    // solved to ~1e-12 (the default 1e-10 targets only the ~1e-8 parity in task 15).
    gsPlaneStressCondensation<real_t> cond(&law, 1e-12);

    gsMatrix<real_t> S2D, C2D, E33, dE33dE;
    cond.condense(e2D, makeParams(N,E,nu), S2D, C2D, E33, &dE33dE);

    // (a) Independently re-evaluate the FULL 3D stress at the converged (e2D,E33) and
    //     check the out-of-plane component S33 vanishes.
    gsMaterialData<real_t> data = makeData(e2D, E33, makeParams(N,E,nu));
    gsMatrix<real_t> Sfull;
    law.compute_stress_into(data, Sfull);
    real_t maxS33 = 0;
    for (index_t k=0; k!=N; ++k)
    {
        const real_t S33 = Sfull(8,k);          // (2,2) of 3x3, col-major flat 8
        const real_t S11 = Sfull(0,k);
        maxS33 = math::max(maxS33, math::abs(S33) / (1.0 + math::abs(S11)));
        CHECK(math::abs(S33) <= 1e-9*(1.0+math::abs(S11)));
    }
    gsInfo<<"[NH] max |S33|/(1+|S11|) = "<<maxS33<<"\n";

    // (b) central FD of the CONDENSED stress: (S2D(e+h δ_a) − S2D(e−h δ_a))/(2h) vs C2D·δ_a.
    //     This is convention-free: it perturbs the real engineering input e2D(2)=2E12 and
    //     measures real S12, so it independently pins the Voigt map and the shear factor.
    const real_t h = 1e-6;
    real_t maxTang = 0;
    for (index_t k=0; k!=N; ++k)
    {
        gsMatrix<real_t> Ck = C2D.reshapeCol(k,3,3);
        for (index_t a=0; a!=3; ++a)
        {
            gsMatrix<real_t> ep = e2D.col(k), em = e2D.col(k);
            ep(a,0) += h;  em(a,0) -= h;
            gsMatrix<real_t> Sp, Cp, E33p, Sm, Cm, E33m;
            cond.condense(ep, makeParams(1,E,nu), Sp, Cp, E33p);
            cond.condense(em, makeParams(1,E,nu), Sm, Cm, E33m);
            gsMatrix<real_t> fd = (Sp - Sm) / (2.0*h);   // column a of dS2D/de2D
            gsMatrix<real_t> ref = Ck.col(a);
            maxTang = math::max(maxTang, relFro(fd, ref));
        }
    }
    gsInfo<<"[NH] max rel FD-vs-C2D tangent err = "<<maxTang<<"\n";
    CHECK(maxTang < 1e-5);

    // (c) central FD of E33 w.r.t. e2D vs the reported dE33dE.
    real_t maxdE33 = 0;
    for (index_t k=0; k!=N; ++k)
    {
        for (index_t a=0; a!=3; ++a)
        {
            gsMatrix<real_t> ep = e2D.col(k), em = e2D.col(k);
            ep(a,0) += h;  em(a,0) -= h;
            gsMatrix<real_t> Sp, Cp, E33p, Sm, Cm, E33m;
            cond.condense(ep, makeParams(1,E,nu), Sp, Cp, E33p);
            cond.condense(em, makeParams(1,E,nu), Sm, Cm, E33m);
            const real_t fd  = (E33p(0,0) - E33m(0,0)) / (2.0*h);
            const real_t ref = dE33dE(a,k);
            maxdE33 = math::max(maxdE33, math::abs(fd-ref)/(1.0+math::abs(ref)));
        }
    }
    gsInfo<<"[NH] max rel FD-vs-dE33dE err = "<<maxdE33<<"\n";
    CHECK(maxdE33 < 1e-6);
}

// 4. Zero-strain limits: E33 = 0 for both laws; the NH tangent at F = I equals the
//    isotropic SvK plane-stress matrix (hyperelastic linearises to linear elasticity).
TEST(zero_strain_limits)
{
    const real_t E = 200.0, nu = 0.3;
    const index_t N = 1;
    gsMatrix<real_t> e2D = gsMatrix<real_t>::Zero(3,N);
    const gsMatrix<real_t> Cref = svkC2D(E,nu);

    // Linear law.
    {
        gsLinearMaterial<real_t> law(E, nu, 3);
        gsPlaneStressCondensation<real_t> cond(&law);
        gsMatrix<real_t> S2D, C2D, E33;
        cond.condense(e2D, makeParams(N,E,nu), S2D, C2D, E33);
        CHECK(math::abs(E33(0,0)) <= 1e-14);
        gsMatrix<real_t> C = C2D.reshapeCol(0,3,3);
        gsInfo<<"[zero linear] |E33|="<<math::abs(E33(0,0))<<"  relC="<<relFro(C,Cref)<<"\n";
        CHECK(relFro(C, Cref) < 1e-12);
    }
    // Neo-Hooke law: linearises to the SAME isotropic elasticity tensor at the identity.
    {
        gsNeoHookeMaterial<real_t> law(E, nu, 3);
        gsPlaneStressCondensation<real_t> cond(&law);
        gsMatrix<real_t> S2D, C2D, E33;
        cond.condense(e2D, makeParams(N,E,nu), S2D, C2D, E33);
        CHECK(math::abs(E33(0,0)) <= 1e-14);
        gsMatrix<real_t> C = C2D.reshapeCol(0,3,3);
        gsInfo<<"[zero NH] |E33|="<<math::abs(E33(0,0))<<"  relC="<<relFro(C,Cref)<<"\n";
        CHECK(relFro(C, Cref) < 1e-10);
    }
}

// 5. Negative control: a reference perturbed by (1+1e-6) must FAIL the 1e-12 comparison,
//    proving the tolerance actually bites (pattern: gsMaterialProvider_test negative control).
TEST(negative_control)
{
    const real_t E = 200.0, nu = 0.3;
    const index_t N = 1;
    gsMatrix<real_t> e2D(3,N); e2D << 0.02, -0.01, 0.03;

    gsLinearMaterial<real_t> law(E, nu, 3);
    gsPlaneStressCondensation<real_t> cond(&law);
    gsMatrix<real_t> S2D, C2D, E33;
    cond.condense(e2D, makeParams(N,E,nu), S2D, C2D, E33);

    gsMatrix<real_t> C = C2D.reshapeCol(0,3,3);
    gsMatrix<real_t> Cbad = svkC2D(E,nu) * (1.0 + 1e-6);   // deliberately wrong reference
    const real_t rel = relFro(C, Cbad);
    gsInfo<<"[negative_control] rel err vs perturbed ref = "<<rel<<" (must be > 1e-8)\n";
    CHECK(rel > 1e-8);
}

} // SUITE

#endif // gsPhaseFieldFracture_ENABLED
