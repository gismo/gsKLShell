/** @file gsMaterialMatrix3D_test.cpp

    @brief Pointwise adapter-vs-legacy parity tests for gsMaterialMatrix3D (task 15).

    This suite is THE frame-transform gate for the task-14 adapter
    (gsMaterialMatrix3D), which wraps a 3D gsPhaseFieldFracture (PFF) material law
    behind the classic gsMaterialMatrixBase interface. Covariant/contravariant
    convention mismatches are the classic silent bug of this migration, so the
    references here are machine-precision ORACLES:

      - The legacy gsMaterialMatrixLinear<3,real_t> St.Venant-Kirchhoff response is
        used as an exact pointwise oracle for the adapter running gsLinearMaterial.
        On a FLAT undeformed plate the legacy z-constant tangent equals the exact
        metric tangent at every z; on a doubly-curved patch it is exact ONLY at
        z=0 (Gcov=Acov there) — that z=0 curved test cannot be passed by a
        frame-flipped adapter and is the true gate.
      - The legacy compressible Neo-Hooke (gsMaterialMatrixNonlinear, Analytical,
        Compressibility=true) is the oracle for the model-matched PFF NH law. The
        match is established numerically in TEST(NH_model_match) below.

    Conventions (both sides): contravariant curvilinear components, shell Voigt-3
    order [11,22,12], column layout colIdx = j*u.cols()+k (eval3D_matrix 9 rows,
    eval3D_vector 3 rows). All material PARAMETERS are CONSTANT (E, nu): the adapter
    samples parameters at the parametric points u while the legacy material samples
    them at the physical points (task-14 report hand-off) — the two coincide only
    for constant parameters, so constant E/nu is used everywhere by construction.

    Fixed-size Eigen math uses gsEigen::Matrix<T,3,3> (the gsMatrix<T,3,3> wrapper
    has no Eigen evaluator for the decomposition paths — task-12 addendum).

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsMaterialMatrix3D.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/getMaterialMatrix.h>          // legacy NH via getMaterialMatrix<3,real_t>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>

#include <gsPhaseFieldFracture/materials/gsLinearMaterial.h>
#include <gsPhaseFieldFracture/materials/gsNeoHookeMaterial.h>
#include <gsPhaseFieldFracture/materials/gsNeoHookeQuadMaterial.h>   // QuadraticVolumetric alias
#include <gsPhaseFieldFracture/materials/gsNeoHookeLogMaterial.h>    // LogarithmicVolumetric alias

SUITE(gsMaterialMatrix3D)
{

// =====================================================================
// Shared constants and small helpers (analytic / geometric, no rand())
// =====================================================================

// O(1) material constants; all tolerances below are relative.
const real_t E_MOD = 200.0;
const real_t NU    = 0.3;

// Relative Frobenius error, robust for near-zero references (sibling-test style).
real_t relFro(const gsMatrix<real_t> & a, const gsMatrix<real_t> & ref)
{
    return (a - ref).norm() / (1.0 + ref.norm());
}

// Flat unit-square shell, 3D-embedded, degree-2, one refinement. Undeformed
// curvature B_ori = 0 ⇒ the legacy z-constant Linear tangent is exact at ALL z.
gsMultiPatch<real_t> makeFlat()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.embed(3);
    mp.patch(0).degreeElevate(1);
    mp.patch(0).uniformRefine();
    return mp;
}

// Doubly-curved patch: a smooth, NON-rational, well-conditioned paraboloid-style
// bump z = a X^2 + b Y^2 + c X Y over the unit square, degree-3, one refinement.
// This has genuine double curvature (Gcov ≠ Acov for z ≠ 0, non-zero B_ori) so a
// covariant/contravariant frame flip in the adapter produces O(1) errors — while
// avoiding the conditioning pitfalls of the rational, radius-10 eighth_sphere
// fixture (a near-degenerate pole corner floors pointwise parity at ~1e-6). |g|~1
// keeps the metric well-scaled so the SvK oracle is a genuine machine-precision
// gate. The undeformed bump is fixed; deformSmooth() supplies the deformation.
gsMultiPatch<real_t> makeCurved()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.embed(3);
    mp.patch(0).degreeElevate(2);          // bicubic in-plane
    mp.patch(0).uniformRefine();
    // Lift the control net onto a smooth doubly-curved bump.
    gsMatrix<real_t> & c = mp.patch(0).coefs();
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0), Y = c(i,1);
        c(i,2) = 0.30*X*X + 0.20*Y*Y + 0.15*X*Y;   // z(X,Y): both principal curvatures
    }
    GISMO_ENSURE(mp.nPatches()==1 && mp.targetDim()==3, "curved patch build failed.");
    return mp;
}

// Deterministic, smooth, NON-homogeneous displacement of the control net (no
// rand()): each control point is moved by a low-order polynomial of its own
// position, so the strain field varies across the patch while the deformed map
// stays a valid spline. The polynomial arguments are NORMALIZED by the geometry
// scale (max |coord|) so the resulting strains are O(s) REGARDLESS of the patch
// size — crucial here because the eighth_sphere fixture has radius 10, and an
// un-normalized X*X term would produce an unphysical (~30x) distortion. The
// scale @a s tunes the strain magnitude (O(s) relative strain).
gsMultiPatch<real_t> deformSmooth(const gsMultiPatch<real_t> & mp, real_t s)
{
    gsMultiPatch<real_t> def = mp;
    gsMatrix<real_t> & c = def.patch(0).coefs();
    const real_t L = math::max(c.cwiseAbs().maxCoeff(), (real_t)1.0); // geometry scale
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0)/L, Y = c(i,1)/L, Z = c(i,2)/L; // normalized in [-1,1]
        c(i,0) += s * L * ( 0.20*X*X + 0.10*X*Y + 0.05*Y );
        c(i,1) += s * L * (-0.15*Y*Y + 0.08*X   + 0.03*Z );
        c(i,2) += s * L * ( 0.12*X*Y - 0.06*Y   + 0.04*X );
    }
    return def;
}

// Homogeneous diagonal stretch (λx, λy) on a flat plate: curvilinear == Cartesian,
// used to establish the NH model match at a single point.
gsMultiPatch<real_t> stretchDiag(const gsMultiPatch<real_t> & mp, real_t lx, real_t ly)
{
    gsMultiPatch<real_t> def = mp;
    gsMatrix<real_t> & c = def.patch(0).coefs();
    for (index_t i=0; i!=c.rows(); ++i)
    {
        c(i,0) *= lx;
        c(i,1) *= ly;
    }
    return def;
}

// A handful of interior parametric points (in (0,1)^2), deterministic.
gsMatrix<real_t> interiorPoints()
{
    gsMatrix<real_t> u(2,6);
    u.col(0) << 0.20, 0.30;
    u.col(1) << 0.50, 0.40;
    u.col(2) << 0.70, 0.60;
    u.col(3) << 0.35, 0.75;
    u.col(4) << 0.60, 0.25;
    u.col(5) << 0.80, 0.85;
    return u;
}

// Build the (nz x N) z-matrix the integrator passes: per point a column of the
// SAME dimensionless z nodes.
gsMatrix<real_t> makeZ(const std::vector<real_t> & zs, index_t N)
{
    gsMatrix<real_t> z((index_t)zs.size(), N);
    for (size_t j=0; j!=zs.size(); ++j)
        z.row((index_t)j).setConstant(zs[j]);
    return z;
}

// Constant thickness / parameter FunctionExprs on the PHYSICAL domain (dim 3):
// the legacy metric engine evaluates thickness and legacy parameters at the
// physical map points (gsMaterialMatrixBaseDim.hpp:447-455).
gsFunctionExpr<real_t> constFun3(real_t v)
{ return gsFunctionExpr<real_t>(util::to_string(v), 3); }

// =====================================================================
// TEST 1 : SvK flat, pointwise, all z rows (matrix) + z=0 (vector)
// =====================================================================
TEST(SvK_flat_pointwise)
{
    gsMultiPatch<real_t> mp     = makeFlat();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.30);

    const real_t thick = 0.01;
    gsFunctionExpr<real_t> t  = constFun3(thick);
    gsFunctionExpr<real_t> Ef = constFun3(E_MOD);
    gsFunctionExpr<real_t> nf = constFun3(NU);

    // Legacy SvK oracle.
    gsMaterialMatrixLinear<3,real_t> legacy(mp, t, Ef, nf);
    legacy.setDeformed(&mp_def);

    // Adapter over the PFF linear law. The law's parameter functions live on the
    // PARAMETRIC domain (dim 2): the adapter evaluates them at the parametric u.
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
    adapter.setDeformed(&mp_def);

    const gsMatrix<real_t> u = interiorPoints();
    const index_t N = u.cols();

    // Matrix: on a FLAT plate the undeformed metric is z-independent, so the
    // legacy z-constant C is exact at every z row and the adapter (undeformed-
    // metric transform + strain-independent linear tangent) must reproduce it.
    const std::vector<real_t> zs = {-0.5, 0.0, 0.37};
    const gsMatrix<real_t> z = makeZ(zs, N);

    gsMatrix<real_t> Cleg = legacy.eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    gsMatrix<real_t> Cad  = adapter.eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    CHECK_EQUAL(Cleg.rows(), Cad.rows());
    CHECK_EQUAL(Cleg.cols(), Cad.cols());

    real_t maxC = 0;
    for (index_t col=0; col!=Cleg.cols(); ++col)
    {
        gsMatrix<real_t> a = Cad .reshapeCol(col,3,3);
        gsMatrix<real_t> r = Cleg.reshapeCol(col,3,3);
        maxC = math::max(maxC, relFro(a,r));
    }
    gsInfo << "[SvK_flat] max rel matrix err (all z) = " << maxC << "\n";
    CHECK(maxC < 1e-10);

    // Vector at z=0: the adapter's Gcov strain equals the legacy VectorN membrane
    // strain there, so the condensed S^ij must match the legacy N-response.
    const gsMatrix<real_t> z0 = makeZ({0.0}, N);
    gsMatrix<real_t> Sleg = legacy.eval3D_vector(0, u, z0, MaterialOutput::VectorN);
    gsMatrix<real_t> Sad  = adapter.eval3D_vector(0, u, z0, MaterialOutput::VectorN);
    real_t maxS = 0;
    for (index_t col=0; col!=Sleg.cols(); ++col)
        maxS = math::max(maxS, relFro(Sad.col(col), Sleg.col(col)));
    gsInfo << "[SvK_flat] max rel vector err (z=0) = " << maxS << "\n";
    CHECK(maxS < 1e-10);
}

// =====================================================================
// TEST 2 : SvK curved, z=0 — THE frame gate
// =====================================================================
TEST(SvK_curved_z0_pointwise)
{
    // On a DOUBLY-CURVED patch Gcov=Acov holds only at z=0, so the legacy Linear
    // A-metric tangent is exact there. A covariant/contravariant frame flip in the
    // adapter's two _transformation orientations CANNOT pass this test: it would
    // rotate/scale S^ij and C^ijkl by the (non-identity) curvilinear metric,
    // producing O(1) errors while leaving the invariants intact.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    const real_t thick = 0.01;
    gsFunctionExpr<real_t> t  = constFun3(thick);
    gsFunctionExpr<real_t> Ef = constFun3(E_MOD);
    gsFunctionExpr<real_t> nf = constFun3(NU);

    gsMaterialMatrixLinear<3,real_t> legacy(mp, t, Ef, nf);
    legacy.setDeformed(&mp_def);

    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
    adapter.setDeformed(&mp_def);

    const gsMatrix<real_t> u = interiorPoints();
    const index_t N = u.cols();
    const gsMatrix<real_t> z0 = makeZ({0.0}, N);

    gsMatrix<real_t> Cleg = legacy.eval3D_matrix(0, u, z0, MaterialOutput::MatrixA);
    gsMatrix<real_t> Cad  = adapter.eval3D_matrix(0, u, z0, MaterialOutput::MatrixA);
    real_t maxC = 0;
    for (index_t col=0; col!=Cleg.cols(); ++col)
    {
        gsMatrix<real_t> a = Cad .reshapeCol(col,3,3);
        gsMatrix<real_t> r = Cleg.reshapeCol(col,3,3);
        maxC = math::max(maxC, relFro(a,r));
    }
    gsInfo << "[SvK_curved] max rel matrix err (z=0) = " << maxC << "\n";
    CHECK(maxC < 1e-10);

    gsMatrix<real_t> Sleg = legacy.eval3D_vector(0, u, z0, MaterialOutput::VectorN);
    gsMatrix<real_t> Sad  = adapter.eval3D_vector(0, u, z0, MaterialOutput::VectorN);
    real_t maxS = 0;
    for (index_t col=0; col!=Sleg.cols(); ++col)
        maxS = math::max(maxS, relFro(Sad.col(col), Sleg.col(col)));
    gsInfo << "[SvK_curved] max rel vector err (z=0) = " << maxS << "\n";
    CHECK(maxS < 1e-10);
}

// =====================================================================
// TEST 3 : SvK integrated moments (honest model-difference + convergence)
// =====================================================================

// Rel Frobenius diff between adapter and legacy integrated MatrixA at thickness t.
static real_t matrixA_reldiff(const gsMultiPatch<real_t> & mp,
                              const gsMultiPatch<real_t> & mp_def,
                              real_t thick, const gsMatrix<real_t> & u)
{
    gsFunctionExpr<real_t> t  = constFun3(thick);
    gsFunctionExpr<real_t> Ef = constFun3(E_MOD);
    gsFunctionExpr<real_t> nf = constFun3(NU);

    gsMaterialMatrixLinear<3,real_t> legacy(mp, t, Ef, nf);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);

    // The integrator ctor's deformed argument triggers setDeformed on the material.
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> iL(&legacy , &mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> iA(&adapter, &mp_def);

    gsMatrix<real_t> AL, AA;
    iL.piece(0).eval_into(u, AL);
    iA.piece(0).eval_into(u, AA);
    return relFro(AA, AL);
}

TEST(SvK_integrated_moments)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    const gsMatrix<real_t> u = interiorPoints();

    // Legacy integrates a z-CONSTANT C analytically (Constant mode); the adapter
    // Gauss-integrates the EXACT z-dependent metric tangent. The gap is model error
    // O(t^2/R^2), NOT roundoff — small but nonzero at t=0.01 on this unit-curvature
    // fixture, and it must shrink ~quadratically as t is halved.
    const real_t d1 = matrixA_reldiff(mp, mp_def, 0.01,  u);
    const real_t d2 = matrixA_reldiff(mp, mp_def, 0.005, u);
    gsInfo << "[SvK_integrated] MatrixA rel diff: t=0.01 -> " << d1
           << ",  t=0.005 -> " << d2 << ",  ratio = " << (d1/d2) << "\n";

    // (a) bounded model difference at t=0.01.
    CHECK(d1 < 5e-4);
    // Sanity: it is a genuine (nonzero) model difference, not machine noise.
    CHECK(d1 > 1e-9);
    // (b) quadratic trend: halving t shrinks the difference by >= 3x (expect ~4x).
    CHECK(d1 / d2 > 3.0);

    // Also drive VectorN through the integrator to exercise the moment/vector path.
    {
        gsFunctionExpr<real_t> t  = constFun3(0.01);
        gsFunctionExpr<real_t> Ef = constFun3(E_MOD);
        gsFunctionExpr<real_t> nf = constFun3(NU);
        gsMaterialMatrixLinear<3,real_t> legacy(mp, t, Ef, nf);
        gsLinearMaterial<real_t> law(E_MOD, NU, 2);
        gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
        gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> nL(&legacy , &mp_def);
        gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> nA(&adapter, &mp_def);
        gsMatrix<real_t> NL, NA;
        nL.piece(0).eval_into(u, NL);
        nA.piece(0).eval_into(u, NA);
        const real_t dn = relFro(NA, NL);
        gsInfo << "[SvK_integrated] VectorN rel diff (t=0.01) = " << dn << "\n";
        CHECK(dn < 5e-4);
    }
}

// =====================================================================
// TEST 4a : NH model matching (do this FIRST — it decides test 4b)
// =====================================================================
TEST(NH_model_match)
{
    // Single point on a FLAT plate with a homogeneous stretch (λx=1.1, λy=0.95):
    // curvilinear == Cartesian, so the legacy compressible NH pointwise stress at
    // z=0 is directly comparable to the PFF NH candidates through the adapter.
    //
    // Analytic expectation (verified in code below): the legacy compressible NH
    // (gsMaterialMatrixNonlinear, Analytical, comp=true) uses an isochoric part
    // mu*J^{-2/3}(I - I1/3 C^{-1}) and a volumetric part with bulk modulus
    // K = E/(3(1-2nu)) and dPsi_vol/dC = K/4 (J^2-1) C^{-1}, i.e.
    // U(J) = (K/4)(J^2 - 1 - 2 ln J) — the Simo-Taylor-Pister QUADRATIC policy.
    // The PFF gsNeoHooke has the identical isochoric split and kappa = lambda+2mu/3
    // = E/(3(1-2nu)) = K, so:
    //   - gsNeoHookeQuadMaterial (QuadraticVolumetric) MATCHES,
    //   - gsNeoHookeLogMaterial  (U = kappa/2 ln^2 J)   is REJECTED (O(1) error).
    gsMultiPatch<real_t> mp     = makeFlat();
    gsMultiPatch<real_t> mp_def = stretchDiag(mp, 1.1, 0.95);

    const real_t thick = 0.01;
    gsFunctionExpr<real_t> t   = constFun3(thick);
    gsFunctionExpr<real_t> Ef  = constFun3(E_MOD);
    gsFunctionExpr<real_t> nf  = constFun3(NU);
    gsFunctionExpr<real_t> rho = constFun3(1.0);

    // Legacy compressible NH via getMaterialMatrix (Material=NH, comp=true, Analytical).
    std::vector<gsFunctionSet<real_t>*> pars(2);
    pars[0] = &Ef; pars[1] = &nf;
    gsOptionList opts;
    opts.addInt   ("Material","",(index_t)Material::NH);
    opts.addSwitch("Compressibility","",true);
    opts.addInt   ("Implementation","",(index_t)Implementation::Analytical);
    gsMaterialMatrixBase<real_t>::uPtr legacyNH =
        getMaterialMatrix<3,real_t>(mp, t, pars, rho, opts);
    legacyNH->setDeformed(&mp_def);

    const gsMatrix<real_t> u = interiorPoints();
    const index_t N = u.cols();
    const gsMatrix<real_t> z0 = makeZ({0.0}, N);
    gsMatrix<real_t> Sref = legacyNH->eval3D_vector(0, u, z0, MaterialOutput::VectorN);

    // Candidate 1: QuadraticVolumetric (expected match).
    real_t errQuad = 0;
    {
        gsNeoHookeQuadMaterial<real_t> law(E_MOD, NU, 2);
        gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
        adapter.setDeformed(&mp_def);
        gsMatrix<real_t> S = adapter.eval3D_vector(0, u, z0, MaterialOutput::VectorN);
        for (index_t col=0; col!=N; ++col)
            errQuad = math::max(errQuad, relFro(S.col(col), Sref.col(col)));
    }
    // Candidate 2: LogarithmicVolumetric (expected rejection).
    real_t errLog = 0;
    {
        gsNeoHookeLogMaterial<real_t> law(E_MOD, NU, 2);
        gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
        adapter.setDeformed(&mp_def);
        gsMatrix<real_t> S = adapter.eval3D_vector(0, u, z0, MaterialOutput::VectorN);
        for (index_t col=0; col!=N; ++col)
            errLog = math::max(errLog, relFro(S.col(col), Sref.col(col)));
    }

    gsInfo << "[NH_model_match] candidate rel errors vs legacy comp-NH:\n"
           << "    gsNeoHookeQuadMaterial (QuadraticVolumetric) = " << errQuad << "  <- MATCHED\n"
           << "    gsNeoHookeLogMaterial  (LogarithmicVolumetric) = " << errLog << "  <- rejected\n";

    // The quadratic policy is the matched law (hard-coded choice used by test 4b).
    CHECK(errQuad < 1e-8);
    // The log policy must be clearly distinguishable (negative control on the match).
    CHECK(errLog  > 1e-3);
}

// =====================================================================
// TEST 4b : NH curved pointwise (matched law = gsNeoHookeQuadMaterial)
// =====================================================================
TEST(NH_curved_pointwise)
{
    // Matched law: gsNeoHookeQuadMaterial (established in NH_model_match). Curved
    // patch, moderate strains, z rows {-0.5, 0.25}. Both sides solve an INDEPENDENT
    // C33 Newton at tol 1e-10 (the condenser's default and the legacy Analytical
    // C33 iteration), so the parity floor is ~1e-8, not 1e-12: two independent
    // ~1e-10 root-finds contribute their tolerances to the compared response.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.12);   // moderate strains

    const real_t thick = 0.01;
    gsFunctionExpr<real_t> t   = constFun3(thick);
    gsFunctionExpr<real_t> Ef  = constFun3(E_MOD);
    gsFunctionExpr<real_t> nf  = constFun3(NU);
    gsFunctionExpr<real_t> rho = constFun3(1.0);

    std::vector<gsFunctionSet<real_t>*> pars(2);
    pars[0] = &Ef; pars[1] = &nf;
    gsOptionList opts;
    opts.addInt   ("Material","",(index_t)Material::NH);
    opts.addSwitch("Compressibility","",true);
    opts.addInt   ("Implementation","",(index_t)Implementation::Analytical);
    gsMaterialMatrixBase<real_t>::uPtr legacyNH =
        getMaterialMatrix<3,real_t>(mp, t, pars, rho, opts);
    legacyNH->setDeformed(&mp_def);

    gsNeoHookeQuadMaterial<real_t> law(E_MOD, NU, 2);
    gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
    adapter.setDeformed(&mp_def);

    const gsMatrix<real_t> u = interiorPoints();
    const index_t N = u.cols();
    const gsMatrix<real_t> z = makeZ({-0.5, 0.25}, N);

    gsMatrix<real_t> Cleg = legacyNH->eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    gsMatrix<real_t> Cad  = adapter .eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    real_t maxC = 0;
    for (index_t col=0; col!=Cleg.cols(); ++col)
    {
        gsMatrix<real_t> a = Cad .reshapeCol(col,3,3);
        gsMatrix<real_t> r = Cleg.reshapeCol(col,3,3);
        maxC = math::max(maxC, relFro(a,r));
    }
    gsInfo << "[NH_curved] max rel matrix err (z=-0.5,0.25) = " << maxC << "\n";
    CHECK(maxC < 1e-8);

    gsMatrix<real_t> Sleg = legacyNH->eval3D_vector(0, u, z, MaterialOutput::VectorN);
    gsMatrix<real_t> Sad  = adapter .eval3D_vector(0, u, z, MaterialOutput::VectorN);
    real_t maxS = 0;
    for (index_t col=0; col!=Sleg.cols(); ++col)
        maxS = math::max(maxS, relFro(Sad.col(col), Sleg.col(col)));
    gsInfo << "[NH_curved] max rel vector err (z=-0.5,0.25) = " << maxS << "\n";
    CHECK(maxS < 1e-8);
}

// =====================================================================
// TEST 5 : cross-output cache — six integrators collapse to ONE sweep
// =====================================================================
TEST(cache_six_to_one)
{
    // Pinned to the DIRECTLY-MEASURED task-14 pattern (tasks/14-report.md, "MEASURED
    // 6->1"): all six legacy integrators route eval3D through the SAME NumGauss
    // z-grid, so per element one constitutive SWEEP + five cache HITS. The ctors each
    // bump m_configRev via setDeformed, so counters are reset AFTER construction.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.12);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> iA(&adapter, &mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> iB(&adapter, &mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> iC(&adapter, &mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> iD(&adapter, &mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> iN(&adapter, &mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> iM(&adapter, &mp_def);

    const gsMatrix<real_t> u = interiorPoints();

    // Reset AFTER construction: each ctor's setDeformed bump would otherwise be
    // mis-read as extra sweeps.
    gsMaterialMatrix3DResetSweeps();
    gsMaterialMatrix3DResetHits();

    gsMatrix<real_t> tmp;
    iA.piece(0).eval_into(u, tmp);
    iB.piece(0).eval_into(u, tmp);
    iC.piece(0).eval_into(u, tmp);
    iD.piece(0).eval_into(u, tmp);
    iN.piece(0).eval_into(u, tmp);
    iM.piece(0).eval_into(u, tmp);

    const size_t sweeps = gsMaterialMatrix3DSweeps();
    const size_t hits   = gsMaterialMatrix3DHits();
    gsInfo << "[cache_six_to_one] sweeps=" << sweeps << "  hits=" << hits
           << " (task-14 measured: sweeps=1, hits=5)\n";
    CHECK_EQUAL((size_t)1, sweeps);
    CHECK_EQUAL((size_t)5, hits);
}

// =====================================================================
// TEST 6 : setDeformed invalidates the adapter cache (stale-cache regression)
// =====================================================================
TEST(setDeformed_invalidates)
{
    // The task-11-style regression, now on the adapter cross-output cache: mutating
    // the deformed coefs IN PLACE (stable pointer) and re-calling setDeformed must
    // bump m_configRev, force a NEW sweep, and change the response.
    //
    // A NONLINEAR law (NH) is used deliberately: the MatrixA TANGENT of a linear
    // SvK law is strain-INDEPENDENT (built from the undeformed metric only), so its
    // value would NOT change when only the deformed config changes — that would test
    // nothing about the response. The NH tangent depends on the deformed state, so a
    // genuine value change confirms the cache actually recomputed (not a stale hit).
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.12);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsNeoHookeQuadMaterial<real_t> law(E_MOD, NU, 2);
    gsMaterialMatrix3D<3,real_t> adapter(mp, t, law);
    adapter.setDeformed(&mp_def);

    const gsMatrix<real_t> u = interiorPoints();
    const gsMatrix<real_t> z = makeZ({-0.3, 0.2}, u.cols());

    gsMaterialMatrix3DResetSweeps();
    gsMaterialMatrix3DResetHits();

    gsMatrix<real_t> C0 = adapter.eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    CHECK_EQUAL((size_t)1, gsMaterialMatrix3DSweeps());      // first eval = one sweep
    gsMatrix<real_t> C0b = adapter.eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    CHECK_EQUAL((size_t)1, gsMaterialMatrix3DHits());        // repeat = one hit, no new sweep
    CHECK_EQUAL((size_t)1, gsMaterialMatrix3DSweeps());

    // Mutate the deformed geometry in place (same pointer), then re-set it.
    {
        gsMatrix<real_t> & c = mp_def.patch(0).coefs();
        for (index_t i=0; i!=c.rows(); ++i)
            c(i,0) += 0.05 * c(i,1);   // deterministic shear-like perturbation
    }
    adapter.setDeformed(&mp_def);

    gsMaterialMatrix3DResetSweeps();
    gsMaterialMatrix3DResetHits();
    gsMatrix<real_t> C1 = adapter.eval3D_matrix(0, u, z, MaterialOutput::MatrixA);
    CHECK_EQUAL((size_t)1, gsMaterialMatrix3DSweeps());      // config changed => new sweep
    CHECK_EQUAL((size_t)0, gsMaterialMatrix3DHits());

    const real_t diff = (C1 - C0).norm();
    gsInfo << "[setDeformed_invalidates] |C1 - C0| = " << diff << "\n";
    CHECK(diff > 0.0);                                       // the response actually changed
}

} // SUITE

#endif // gsPhaseFieldFracture_ENABLED
