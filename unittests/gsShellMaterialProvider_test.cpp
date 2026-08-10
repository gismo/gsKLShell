/** @file gsShellMaterialProvider_test.cpp

    @brief THE numerical gate for gsShellMaterialProvider (task 23).

    Until this suite existed, the correctness of gsShellMaterialProvider — the
    Phase-4 keystone that replaces SIX gsMaterialMatrixIntegrate coefficients by
    ONE batched per-element sweep — rested entirely on code review.

    ### The oracle, and honestly what it does and does not cover
    The reference is the LEGACY integration path
    (gsMaterialMatrixIntegrate<MaterialOutput::X> over the task-14 adapter
    gsMaterialMatrix3D), driven with the SAME PFF law, the SAME thickness, the
    SAME NumGauss z-grid and the SAME points. Both sides therefore share the
    constitutive law and gsPlaneStressCondensation by construction: this suite
    does NOT re-validate those (tasks 12/16 and gsPlaneStressCondensation_test
    do). What it DOES gate is everything the provider re-implements:

      - the metric engine swap (gsShellKinematics fed by INJECTED gsMapData
        instead of gsMaterialMatrixBaseDim's owned geometry),
      - the moment weights and the dimensionless-z convention
        w(j)*(z_j*t)^m * val * t, for m = 0,1,1,2 (A,B,C,D) and 0,1 (N,M),
      - the col = j*N+k grid ordering shared with the condenser,
      - the thickness evaluated at the UNDEFORMED PHYSICAL map values (the
        TEST(NonConstantThickness) expression depends on the out-of-plane
        coordinate on purpose: on this fixture x,y ARE the parametric
        coordinates to machine precision, so only a z-dependence separates a
        physical from a parametric evaluation),
      - the local-Cartesian frame transform R and its back-transform R^T,
      - the 42-row block layout (A,B,C,D column-major 3x3, then N,M 3x1),
      - the output-request mask semantics, including the wantC2D==false
        condensation branch that NO test executed before this one,
      - the per-element fill/cache protocol.

    ### Fixture
    makeCurved() — the NON-orthogonal doubly-curved bicubic bump
    z = 0.30X^2 + 0.20Y^2 + 0.15XY. This is the fixture class that exposed the
    task-15 frame bug: the near-orthogonal sphere masks a strain/stress
    transform swap to ~1e-6, this one does not. Fixture non-degeneracy
    (curvature and non-orthogonality) is inherited from
    gsMaterialMatrix3D_test.cpp / gsShellKinematics_test.cpp, where it is
    asserted; here the reference-block norms are asserted nonzero instead, so a
    vacuous "0 == 0" comparison cannot pass silently.

    ### Tolerance policy
    Both routes preserve the legacy multiplication grouping ((w*pow)*val)*t and
    the metric transplant was proven bit-exact in task 21, so agreement is
    EXPECTED to be exact. The gate is nonetheless the relative form
    (A-B).norm() <= 1e-12*(1+|B|)  — CHECK_CLOSE is ABSOLUTE
    (optional/gsUnitTest/Checks.h:39) and would be meaningless across blocks
    whose magnitudes differ by t^2 ~ 1e-4. Every block prints its worst ABSOLUTE
    and RELATIVE deviation and whether it is exactly 0.0; tolerances are never
    widened to make a test pass.

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"

#ifdef gsPhaseFieldFracture_ENABLED

// gsShellMaterialExpr.h includes only <gsExpressions/gsExpressions.h>, where
// gsExprHelper is merely FORWARD-DECLARED (gsExpressions.h:23) -- so a TU that
// includes nothing else fails inside the views' parse() on evList.add(...).
// Pulling the assembler in first is the framework convention (identical to the
// solids twin material_expr.h) and is what makes the view header usable here.
#include <gsAssembler/gsExprAssembler.h>

#include <gsKLShell/src/gsShellMaterialProvider.h>
#include <gsKLShell/src/gsShellMaterialExpr.h>
#include <gsKLShell/src/gsMaterialMatrix3D.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>

#include <gsPhaseFieldFracture/materials/gsLinearMaterial.h>
#include <gsPhaseFieldFracture/materials/gsNeoHookeQuadMaterial.h>

#ifdef _OPENMP
#include <omp.h>
#endif

SUITE(gsShellMaterialProvider)
{

// =====================================================================
// Shared constants
// =====================================================================

const real_t E_MOD = 200.0;
const real_t NU    = 0.3;

/// Relative gate. The expectation is EXACT agreement; this is the honest
/// upper bound, not a tuned number (see the file header).
const real_t REL_TOL = 1e-12;

/// Through-thickness Gauss nodes. Must equal the legacy "NumGauss" option,
/// which is ASSERTED in TEST(SixMoments_SvK_curved) rather than assumed.
const index_t NZ = 4;

typedef gsShellMaterialProvider<3,real_t> Provider;

// =====================================================================
// Fixtures (patterns of gsMaterialMatrix3D_test.cpp:85,111,141,166)
// =====================================================================

// Doubly-curved, NON-orthogonal bicubic bump over the unit square. Adapted
// (not included) from gsMaterialMatrix3D_test.cpp:85 — see the file header for
// why this and not a near-orthogonal sphere.
gsMultiPatch<real_t> makeCurved()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.embed(3);
    mp.patch(0).degreeElevate(2);          // bicubic in-plane
    mp.patch(0).uniformRefine();
    gsMatrix<real_t> & c = mp.patch(0).coefs();
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0), Y = c(i,1);
        c(i,2) = 0.30*X*X + 0.20*Y*Y + 0.15*X*Y;   // both principal curvatures
    }
    GISMO_ENSURE(mp.nPatches()==1 && mp.targetDim()==3, "curved patch build failed.");
    return mp;
}

// Deterministic, smooth, NON-homogeneous control-net displacement (no rand()):
// the strain field varies across the patch, so no moment block is constant.
// gsMaterialMatrix3D_test.cpp:111.
gsMultiPatch<real_t> deformSmooth(const gsMultiPatch<real_t> & mp, real_t s)
{
    gsMultiPatch<real_t> def = mp;
    gsMatrix<real_t> & c = def.patch(0).coefs();
    const real_t L = math::max(c.cwiseAbs().maxCoeff(), (real_t)1.0); // geometry scale
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0)/L, Y = c(i,1)/L, Z = c(i,2)/L;
        c(i,0) += s * L * ( 0.20*X*X + 0.10*X*Y + 0.05*Y );
        c(i,1) += s * L * (-0.15*Y*Y + 0.08*X   + 0.03*Z );
        c(i,2) += s * L * ( 0.12*X*Y - 0.06*Y   + 0.04*X );
    }
    return def;
}

// A handful of interior parametric points, deterministic (the "element".)
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

// A DIFFERENT element: different points AND a different point count, so the
// per-element resize path is exercised too (TEST(FillCountPerElement)).
gsMatrix<real_t> otherPoints()
{
    gsMatrix<real_t> u(2,4);
    u.col(0) << 0.11, 0.13;
    u.col(1) << 0.42, 0.91;
    u.col(2) << 0.88, 0.17;
    u.col(3) << 0.55, 0.55;
    return u;
}

// Thickness / parameter expressions on the PHYSICAL domain (dim 3): the legacy
// engine evaluates the thickness at the physical map points
// (gsMaterialMatrixBaseDim.hpp:447) and gsShellMaterialProvider ENFORCES
// domainDim()==dim for exactly that reason (task-22 report, deviation 7).
gsFunctionExpr<real_t> constFun3(real_t v)
{ return gsFunctionExpr<real_t>(util::to_string(v), 3); }

// =====================================================================
// Map injection
// =====================================================================

// BaseDim's flag word, copied VERBATIM from gsMaterialMatrixBaseDim.hpp:471
// (== :591 for the undeformed twin). Paraphrasing it is the one way this test
// could silently lose bit-exactness with the legacy side.
// NEED_VALUE is MANDATORY for the provider: the thickness is evaluated at
// mapOri.values[0] and gsShellKinematics derives its point count from
// values[0].cols() (map.points is unusable in the real pipeline -- task-22
// repair round 1). NB NEED_JACOBIAN == NEED_DERIV; computeMap only ORs flags in.
unsigned baseDimFlags()
{ return NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2; }

void computeMapLike(const gsGeometry<real_t> & g, const gsMatrix<real_t> & u,
                    index_t patch, gsMapData<real_t> & md)
{
    md.flags  = baseDimFlags();
    md.points = u;
    static_cast<const gsFunction<real_t>&>(g).computeMap(md);
    // _fill compares patchId against the element's patch. Setting it here is
    // what makes the injected data honest: these maps are built by hand, so
    // nothing else would ever fill patchId.
    //
    // NB the protocol guards ARE live and ARE testable in the DEV build — the
    // earlier claim here, that they are compiled out, was wrong twice over:
    //  - CMAKE_CXX_FLAGS_RELWITHDEBINFO is plain "-O2 -g" with NO -DNDEBUG, and
    //    gsDebug.h:105-113 gates GISMO_ASSERT on #ifndef NDEBUG, so this very
    //    patchId check throws std::logic_error there (in a genuine Release build,
    //    -O3 -DNDEBUG, it does not exist -- hence the #ifndef NDEBUG in
    //    TEST(ProtocolGuards));
    //  - most of the provider's guards (the NEED_VALUE checks at
    //    gsShellMaterialProvider.h:664/671, the ctor invariants at :467-474 and,
    //    since task 47, the caller-sequencing checks at :411/:423 and the
    //    request-mask checks at :359/:360) are
    //    GISMO_ENSURE, which is live in EVERY build and throws
    //    std::runtime_error (gsDebug.h:120-124) — two DIFFERENT exception types.
    // TEST(ProtocolGuards) at the end of this file exercises both kinds.
    md.patchId = patch;
}

// Drives ONE element through the provider exactly as the real pipeline does:
// bind the fresh map addresses, signal the element boundary through the PIECE's
// eval_into (which is what calls beginElement), then fill.
// @a mdOri / @a mdDef must outlive the fill: the provider stores their addresses.
void fillElement(const Provider & prov,
                 const gsMultiPatch<real_t> & mp, const gsMultiPatch<real_t> & mp_def,
                 const gsMatrix<real_t> & u,
                 gsMapData<real_t> & mdOri, gsMapData<real_t> & mdDef)
{
    computeMapLike(mp    .patch(0), u, 0, mdOri);
    computeMapLike(mp_def.patch(0), u, 0, mdDef);
    prov.bindKinematics(&mdOri,&mdDef);
    gsMatrix<real_t> tmp;
    prov.piece(0).eval_into(u,tmp);      // element-boundary signal == beginElement
    CHECK_EQUAL((index_t)Provider::TargetDim, tmp.rows());
    CHECK_EQUAL(u.cols(), tmp.cols());
    prov.ensureFilled();
}

// One complete run on a fresh provider with a given request mask. The returned
// matrix is a COPY of the 42 x N cache (the provider and its maps die here).
gsMatrix<real_t> providerRun(const gsMultiPatch<real_t> & mp,
                             const gsMultiPatch<real_t> & mp_def,
                             const gsFunctionSet<real_t> & thick,
                             const gsMaterialBase<real_t> * law,
                             const gsMatrix<real_t> & u,
                             unsigned mask       = Provider::ShellReq_All,
                             size_t * fills      = nullptr,
                             size_t * matFills   = nullptr,
                             size_t * strFills   = nullptr)
{
    Provider prov(law, mp, thick, NZ);
    prov.resetCounters();
    prov.setRequested(mask);
    gsMapData<real_t> mdOri, mdDef;
    fillElement(prov, mp, mp_def, u, mdOri, mdDef);
    if (fills)    *fills    = prov.fillCount();
    if (matFills) *matFills = prov.matrixMomentFills();
    if (strFills) *strFills = prov.stressMomentFills();
    return prov.cachedMoments();
}

// =====================================================================
// The legacy oracle: six independent gsMaterialMatrixIntegrate coefficients
// =====================================================================

struct LegacySix
{
    gsMatrix<real_t> A,B,C,D,N,M;   ///< 9 x N (A..D), 3 x N (N,M)
};

LegacySix legacySix(const gsMultiPatch<real_t> & mp,
                    const gsMultiPatch<real_t> & mp_def,
                    const gsFunctionSet<real_t> & thick,
                    const gsMaterialBase<real_t> & law,
                    const gsMatrix<real_t> & u)
{
    gsMaterialMatrix3D<3,real_t> adapter(mp, thick, law);
    // Each integrator ctor calls setDeformed on the adapter (bump-on-set).
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> iA(&adapter,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> iB(&adapter,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> iC(&adapter,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> iD(&adapter,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> iN(&adapter,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> iM(&adapter,&mp_def);

    LegacySix out;
    iA.piece(0).eval_into(u,out.A);
    iB.piece(0).eval_into(u,out.B);
    iC.piece(0).eval_into(u,out.C);
    iD.piece(0).eval_into(u,out.D);
    iN.piece(0).eval_into(u,out.N);
    iM.piece(0).eval_into(u,out.M);
    return out;
}

// =====================================================================
// Comparison helpers
// =====================================================================

const char * const BLOCK_NAME[6] = {"A","B","C","D","N","M"};
const index_t      BLOCK_OFF [6] = {Provider::OffsetA, Provider::OffsetB, Provider::OffsetC,
                                    Provider::OffsetD, Provider::OffsetN, Provider::OffsetM};
const index_t      BLOCK_ROWS[6] = {9,9,9,9,3,3};

/// Block @a b (0..5) of a 42 x N provider cache, as its own matrix.
gsMatrix<real_t> blockOf(const gsMatrix<real_t> & M42, int b)
{
    GISMO_ENSURE(M42.rows()==(index_t)Provider::TargetDim,
                 "provider cache has "<<M42.rows()<<" rows, expected 42.");
    return M42.middleRows(BLOCK_OFF[b], BLOCK_ROWS[b]);
}

/// Compares all six blocks of the provider cache against the legacy six.
/// @a refMustBeNonzero guards against a vacuous 0-vs-0 comparison; it is
/// switched off only where the reference is legitimately zero (zero strain).
void checkSixBlocks(const char * tag, const gsMatrix<real_t> & M42, const LegacySix & L,
                    bool refMustBeNonzero = true)
{
    const gsMatrix<real_t> * ref[6] = {&L.A,&L.B,&L.C,&L.D,&L.N,&L.M};
    CHECK_EQUAL((index_t)Provider::TargetDim, M42.rows());
    for (int b=0; b!=6; ++b)
    {
        const gsMatrix<real_t> got = blockOf(M42,b);
        const gsMatrix<real_t> & r = *ref[b];
        // Shapes FIRST: (A-B).norm() of two 0x0 matrices is 0, i.e. an unsized
        // reference would make this whole comparison green and vacuous.
        CHECK_EQUAL(BLOCK_ROWS[b], r.rows());
        CHECK_EQUAL(M42.cols(),    r.cols());
        if (r.rows()!=BLOCK_ROWS[b] || r.cols()!=M42.cols()) continue;

        const real_t absDev = (got - r).norm();
        const real_t relDev = absDev / ((real_t)1.0 + r.norm());
        gsInfo << "[" << tag << "] block " << BLOCK_NAME[b]
               << " : |ref| = " << r.norm()
               << " , abs dev = " << absDev
               << " , rel dev = " << relDev
               << (absDev == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
        CHECK(relDev <= REL_TOL);
        if (refMustBeNonzero) CHECK(r.norm() > (real_t)0.0);
    }
}

/// Frobenius norm of block @a b of a provider cache.
real_t blockNorm(const gsMatrix<real_t> & M42, int b)
{ return blockOf(M42,b).norm(); }

// =====================================================================
// TEST 1 : the six moments, SvK on the curved patch — THE parity gate
// =====================================================================
TEST(SixMoments_SvK_curved)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);   // PFF laws: parametric parameter domain
    const gsMatrix<real_t> u = interiorPoints();

    // The two sides must integrate the SAME z-grid. The legacy grid size comes
    // from the material's "NumGauss" option (gsMaterialMatrixIntegrate.hpp:286);
    // gsMaterialMatrix3D registers no such option, so it falls back to 4. Pin the
    // coupling instead of assuming it -- a drift here would look exactly like a
    // provider bug in every moment block below.
    {
        gsMaterialMatrix3D<3,real_t> probe(mp, t, law);
        Provider provProbe(&law, mp, t, NZ);
        gsInfo << "[SvK_curved] NumGauss: legacy = " << probe.options().askInt("NumGauss",4)
               << " , provider = " << provProbe.numGauss() << "\n";
        CHECK_EQUAL(NZ, probe.options().askInt("NumGauss",4));
        CHECK_EQUAL(probe.options().askInt("NumGauss",4), provProbe.numGauss());
    }

    const LegacySix L = legacySix(mp, mp_def, t, law, u);
    const gsMatrix<real_t> M42 = providerRun(mp, mp_def, t, &law, u);

    CHECK_EQUAL(u.cols(), M42.cols());
    checkSixBlocks("SvK_curved", M42, L);
}

// =====================================================================
// TEST 2 : the six moments, Neo-Hooke (C33 Newton inside the condensation)
// =====================================================================
TEST(SixMoments_NH_curved)
{
    // gsNeoHookeQuadMaterial is the law pinned as the legacy compressible-NH
    // twin in gsMaterialMatrix3D_test.cpp TEST(NH_model_match). Here it is used
    // on BOTH sides -- the point is not the constitutive model but that the
    // provider's sweep reproduces the legacy moment integration for a law whose
    // condensation runs a genuine per-point C33 Newton on a deformed state.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.12);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsNeoHookeQuadMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();

    const LegacySix L = legacySix(mp, mp_def, t, law, u);
    const gsMatrix<real_t> M42 = providerRun(mp, mp_def, t, &law, u);

    CHECK_EQUAL(u.cols(), M42.cols());
    checkSixBlocks("NH_curved", M42, L);
}

// =====================================================================
// TEST 3 : non-constant thickness (physical-point evaluation + (z t)^m t)
// =====================================================================
TEST(NonConstantThickness)
{
    // The thickness is a function of the PHYSICAL coordinates (domain dim 3).
    // This test catches (a) a thickness evaluated somewhere other than the
    // UNDEFORMED PHYSICAL map values and (b) the (z*t)^m * t weight
    // mis-assembled when t varies per point -- both invisible for a constant
    // thickness.
    //
    // The z TERM IS LOAD-BEARING, do not drop it. On this fixture the in-plane
    // physical coordinates are the parametric ones to machine precision
    // (BSplineSquare(1) is the identity in x,y and only the z-coefficients are
    // lifted; measured max|x-u| = 2.2e-16, max|y-v| = 1.1e-16). A thickness of
    // x and y ALONE would therefore return the very same numbers whether it was
    // evaluated at physical or at parametric points, and the physical-point
    // convention -- the whole point of this test -- would go ungated. Only the
    // out-of-plane coordinate distinguishes the two, and no evaluation on the
    // 2D parameter domain can produce it.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    gsFunctionExpr<real_t> t("0.006 + 0.004*x + 0.003*y + 0.020*z", 3);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();

    // Fixture non-vacuity: the thickness must genuinely VARY over the element
    // (else this degenerates into a copy of TEST(SixMoments_SvK_curved)) AND the
    // out-of-plane contribution must be a substantial share of it (else the
    // physical-vs-parametric distinction is numerically negligible).
    {
        gsMapData<real_t> md;
        computeMapLike(mp.patch(0), u, 0, md);
        gsMatrix<real_t> Tv;
        t.eval_into(md.values[0], Tv);
        const real_t tmin = Tv.minCoeff(), tmax = Tv.maxCoeff();
        const real_t zmin = md.values[0].row(2).minCoeff();
        const real_t zmax = md.values[0].row(2).maxCoeff();
        gsInfo << "[NonConstThick] thickness over the element: min = " << tmin
               << " , max = " << tmax << " , ratio = " << (tmax/tmin) << "\n";
        gsInfo << "[NonConstThick] out-of-plane coordinate z in [" << zmin << ", " << zmax
               << "] -> z-share of t between " << (0.020*zmin/tmax)
               << " and " << (0.020*zmax/tmax) << "\n";
        CHECK(tmin > (real_t)0.0);
        CHECK(tmax/tmin > (real_t)1.8);
        CHECK(0.020*zmax/tmax > (real_t)0.20);   // the z term really carries the thickness
    }

    const LegacySix L = legacySix(mp, mp_def, t, law, u);
    const gsMatrix<real_t> M42 = providerRun(mp, mp_def, t, &law, u);

    checkSixBlocks("NonConstThick", M42, L);
}

// =====================================================================
// TEST 4 : B and C are bit-identical (the copy shortcut is valid)
// =====================================================================
TEST(BEqualsC)
{
    // B and C are moment 1 of the SAME symmetric tangent
    // (gsMaterialMatrixIntegrate.h:246-260, where MatrixC is annotated
    // "must be 1"), so the provider accumulates the moment-1 block ONCE into
    // rows [9,18) and copies it into [18,27). If this ever fails, the copy
    // shortcut is invalid and B/C must be integrated separately.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsNeoHookeQuadMaterial<real_t> law(E_MOD, NU, 2);  // nonlinear: C is not trivially symmetric
    const gsMatrix<real_t> u = interiorPoints();

    const gsMatrix<real_t> M42 = providerRun(mp, mp_def, t, &law, u);
    const gsMatrix<real_t> B = blockOf(M42,1);
    const gsMatrix<real_t> C = blockOf(M42,2);

    const real_t d = (B - C).norm();
    gsInfo << "[BEqualsC] |B| = " << B.norm() << " , |B - C| = " << d
           << (d == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
    CHECK(B.norm() > (real_t)0.0);      // not a vacuous 0 == 0
    CHECK(d == (real_t)0.0);            // element-wise identical, not merely close
}

// =====================================================================
// TEST 5 : mask N|M — FIRST NUMERICAL VALIDATION OF wantC2D == false
// =====================================================================
TEST(MaskGating_VectorOnly)
{
    // A residual-only request (ShellReq_N|ShellReq_M) is the ONLY caller that
    // makes gsPlaneStressCondensation::condense run with wantC2D == false: the
    // final tangent evaluation and the per-point Schur condensation are skipped
    // entirely. The task-22 review established that this branch EXECUTES but
    // that nothing asserted its output was right. That is what this test does:
    // the N and M blocks are compared against the LEGACY integrators (an oracle
    // outside the provider) as well as against the provider's own all-mask run.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.12);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsNeoHookeQuadMaterial<real_t> law(E_MOD, NU, 2);   // the Newton loop must still run
    const gsMatrix<real_t> u = interiorPoints();

    const LegacySix L = legacySix(mp, mp_def, t, law, u);

    size_t fills=0, matFills=0, strFills=0;
    const gsMatrix<real_t> Mall = providerRun(mp, mp_def, t, &law, u,
                                              Provider::ShellReq_All);
    const gsMatrix<real_t> Mnm  = providerRun(mp, mp_def, t, &law, u,
                                              Provider::ShellReq_N | Provider::ShellReq_M,
                                              &fills, &matFills, &strFills);

    // (a) the condensed tangent was never formed
    gsInfo << "[MaskGating_NM] fillCount = " << fills
           << " , matrixMomentFills = " << matFills
           << " , stressMomentFills = " << strFills << "\n";
    CHECK_EQUAL((size_t)1, fills);
    CHECK_EQUAL((size_t)0, matFills);      // wantC2D == false was taken
    CHECK_EQUAL((size_t)1, strFills);

    // (b) A/B/C/D rows are EXACTLY zero (not stale, not partially filled)
    for (int b=0; b!=4; ++b)
    {
        const real_t nb = blockNorm(Mnm,b);
        gsInfo << "[MaskGating_NM] |" << BLOCK_NAME[b] << "| = " << nb << " (must be 0)\n";
        CHECK(nb == (real_t)0.0);
    }

    // (c) N and M are RIGHT, against the EXTERNAL legacy oracle
    for (int b=4; b!=6; ++b)
    {
        const gsMatrix<real_t> got = blockOf(Mnm,b);
        const gsMatrix<real_t> & r = (b==4 ? L.N : L.M);
        CHECK_EQUAL(BLOCK_ROWS[b], r.rows());
        CHECK_EQUAL(u.cols(),      r.cols());
        const real_t absDev = (got - r).norm();
        const real_t relDev = absDev / ((real_t)1.0 + r.norm());
        gsInfo << "[MaskGating_NM] block " << BLOCK_NAME[b] << " vs LEGACY : |ref| = "
               << r.norm() << " , abs dev = " << absDev << " , rel dev = " << relDev
               << (absDev == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
        CHECK(r.norm() > (real_t)0.0);
        CHECK(relDev <= REL_TOL);
    }

    // (d) and identical to the all-mask (wantC2D == true) run: skipping the
    //     tangent must not perturb the stress the Newton loop converged to.
    for (int b=4; b!=6; ++b)
    {
        const real_t d = (blockOf(Mnm,b) - blockOf(Mall,b)).norm();
        gsInfo << "[MaskGating_NM] block " << BLOCK_NAME[b] << " vs ALL-MASK : abs dev = " << d
               << (d == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
        CHECK(d <= REL_TOL * ((real_t)1.0 + blockNorm(Mall,b)));
    }
}

// =====================================================================
// TEST 6 : single-block masks, including the C-only regression
// =====================================================================
TEST(MaskGating_SingleBlocks)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();

    const gsMatrix<real_t> Mall = providerRun(mp, mp_def, t, &law, u, Provider::ShellReq_All);

    const unsigned single[4] = {Provider::ShellReq_A, Provider::ShellReq_B,
                                Provider::ShellReq_C, Provider::ShellReq_D};
    for (int b=0; b!=4; ++b)
    {
        // The C-only case is the regression: repair round 1 of task 22 fixed a
        // real bug where a C-only mask ALSO left the shared moment-1 accumulator
        // sitting in the B rows, contradicting "non-requested rows are ZERO".
        size_t fills=0, matFills=0, strFills=0;
        const gsMatrix<real_t> M = providerRun(mp, mp_def, t, &law, u, single[b],
                                               &fills, &matFills, &strFills);
        CHECK_EQUAL((size_t)1, fills);
        CHECK_EQUAL((size_t)1, matFills);
        CHECK_EQUAL((size_t)0, strFills);   // no stress moment requested

        const real_t absDev = (blockOf(M,b) - blockOf(Mall,b)).norm();
        const real_t relDev = absDev / ((real_t)1.0 + blockNorm(Mall,b));
        gsInfo << "[MaskGating_" << BLOCK_NAME[b] << "] requested block: |ref| = "
               << blockNorm(Mall,b) << " , abs dev = " << absDev
               << " , rel dev = " << relDev
               << (absDev == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
        CHECK(blockNorm(Mall,b) > (real_t)0.0);
        CHECK(relDev <= REL_TOL);

        // every OTHER block (matrix AND stress) must be exactly zero
        for (int o=0; o!=6; ++o)
        {
            if (o==b) continue;
            const real_t nb = blockNorm(M,o);
            if (nb != (real_t)0.0)
                gsInfo << "[MaskGating_" << BLOCK_NAME[b] << "] LEAK into block "
                       << BLOCK_NAME[o] << " : |" << BLOCK_NAME[o] << "| = " << nb << "\n";
            CHECK(nb == (real_t)0.0);
        }
        gsInfo << "[MaskGating_" << BLOCK_NAME[b] << "] all other blocks exactly zero\n";
    }
}

// =====================================================================
// TEST 7 : a narrower request leaves NO stale matrix data behind
// =====================================================================
TEST(MaskGating_NoStaleData)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();

    Provider prov(&law, mp, t, NZ);
    gsMapData<real_t> mdOri, mdDef;

    // (1) full-mask fill on this element
    prov.setRequested(Provider::ShellReq_All);
    fillElement(prov, mp, mp_def, u, mdOri, mdDef);
    const gsMatrix<real_t> Mall = prov.cachedMoments();
    real_t matNorm = 0;
    for (int b=0; b!=4; ++b) matNorm = math::max(matNorm, blockNorm(Mall,b));
    gsInfo << "[NoStale] after full fill: max |A..D| = " << matNorm << " (must be > 0)\n";
    CHECK(matNorm > (real_t)0.0);

    // (2) narrower N|M fill on the SAME provider and the SAME element
    prov.setRequested(Provider::ShellReq_N | Provider::ShellReq_M);
    fillElement(prov, mp, mp_def, u, mdOri, mdDef);
    const gsMatrix<real_t> & Mnm = prov.cachedMoments();

    const real_t stale = Mnm.topRows((index_t)Provider::OffsetN).norm();
    gsInfo << "[NoStale] after N|M fill: |rows [0,36)| = " << stale
           << (stale == (real_t)0.0 ? "   (EXACTLY 0)" : "   (STALE DATA!)") << "\n";
    CHECK(stale == (real_t)0.0);
    // ... while the stress moments are genuinely there
    CHECK(blockNorm(Mnm,4) > (real_t)0.0);
    CHECK(blockNorm(Mnm,5) > (real_t)0.0);
}

// =====================================================================
// TEST 8 : one fill per element, and the cache actually caches
// =====================================================================
TEST(FillCountPerElement)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u1 = interiorPoints();   // 6 points
    const gsMatrix<real_t> u2 = otherPoints();      // 4 points: also a resize

    Provider prov(&law, mp, t, NZ);
    prov.resetCounters();
    gsMapData<real_t> mdOri, mdDef;

    // element 1
    fillElement(prov, mp, mp_def, u1, mdOri, mdDef);
    CHECK_EQUAL((size_t)1, prov.fillCount());
    CHECK_EQUAL(u1.cols(), prov.cachedPoints());
    const gsMatrix<real_t> M1 = prov.cachedMoments();

    // a second ensureFilled() WITHOUT a new beginElement is a cache HIT
    prov.ensureFilled();
    CHECK_EQUAL((size_t)1, prov.fillCount());
    CHECK((prov.cachedMoments() - M1).norm() == (real_t)0.0);

    // element 2 (both maps recomputed: values[0].cols() must follow the element)
    fillElement(prov, mp, mp_def, u2, mdOri, mdDef);
    CHECK_EQUAL((size_t)2, prov.fillCount());
    CHECK_EQUAL(u2.cols(), prov.cachedPoints());
    const gsMatrix<real_t> M2 = prov.cachedMoments();
    prov.ensureFilled();
    CHECK_EQUAL((size_t)2, prov.fillCount());

    gsInfo << "[FillCount] two elements -> fillCount = " << prov.fillCount()
           << " , matrixMomentFills = " << prov.matrixMomentFills()
           << " , stressMomentFills = " << prov.stressMomentFills()
           << " , cachedPoints = " << prov.cachedPoints() << "\n";
    CHECK_EQUAL((size_t)2, prov.matrixMomentFills());
    CHECK_EQUAL((size_t)2, prov.stressMomentFills());

    // element 2 must be a genuinely DIFFERENT element, else the cache-hit check
    // above proves nothing about invalidation.
    CHECK(M2.cols() != M1.cols());

    // and it is the right answer for element 2, not element 1's leftovers
    const LegacySix L2 = legacySix(mp, mp_def, t, law, u2);
    checkSixBlocks("FillCount_elem2", M2, L2);
}

// =====================================================================
// TEST 9 : ori == def aliasing (what a linear assembly does)
// =====================================================================
TEST(OriEqualsDefAliasing)
{
    // gsExprHelper::getMap DEDUPLICATES, so a linear (undeformed) assembly binds
    // the very same gsMapData address as both ori and def. The two
    // gsShellKinematics compute calls then fill their disjoint member sets from
    // identical data: the strain is identically zero, so N and M must vanish
    // while A..D remain the undeformed tangent moments.
    gsMultiPatch<real_t> mp = makeCurved();

    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();

    // Legacy oracle on an undeformed == deformed configuration.
    const LegacySix L = legacySix(mp, mp, t, law, u);

    Provider prov(&law, mp, t, NZ);
    prov.resetCounters();
    gsMapData<real_t> md;
    computeMapLike(mp.patch(0), u, 0, md);
    prov.bindKinematics(&md,&md);            // THE SAME ADDRESS TWICE
    CHECK(prov.mapOri() == prov.mapDef());
    gsMatrix<real_t> tmp;
    prov.piece(0).eval_into(u,tmp);
    prov.ensureFilled();                     // must not throw
    CHECK_EQUAL((size_t)1, prov.fillCount());

    const gsMatrix<real_t> M42 = prov.cachedMoments();

    // The zero-strain reference for N/M is legitimately ~0, so the non-vacuity
    // guard is switched off here and replaced by the explicit checks below.
    checkSixBlocks("Aliased", M42, L, false);

    const real_t nA = blockNorm(M42,0), nN = blockNorm(M42,4), nM = blockNorm(M42,5);
    gsInfo << "[Aliased] |A| = " << nA << " (must be > 0) , |N| = " << nN
           << " , |M| = " << nM << " (both must be ~ 0)\n";
    CHECK(nA > (real_t)0.0);                        // the tangent moments survive
    CHECK(nN <= REL_TOL * ((real_t)1.0 + nA));      // zero strain => zero stress
    CHECK(nM <= REL_TOL * ((real_t)1.0 + nA));
}

// =====================================================================
// TEST 10 : the protocol guards actually FIRE
// =====================================================================

/// Which exception kind a guard produced. G+Smo uses TWO of them, and they are
/// siblings (neither catches the other): GISMO_ASSERT throws std::logic_error
/// while GISMO_ENSURE throws std::runtime_error. Getting the type wrong would
/// make CHECK_THROW report "expected exception not thrown" for a guard that did
/// fire, so each case below pins the exact kind as well.
///
/// The two kinds ALSO differ in which builds they exist in, and that decides
/// what may be asserted here:
///  - GISMO_ASSERT (gsDebug.h:105-113) is wrapped in `#ifndef NDEBUG`. It is
///    live in the dev build (RelWithDebInfo, from which G+Smo strips -DNDEBUG:
///    cmake/gsConfig.cmake:12-16, CMakeLists.txt:113-115) and ABSENT in a
///    Release build (-O3 -DNDEBUG). Absent means the checked condition is not
///    tested at all, so the offending call proceeds into code the guard was
///    protecting -- an ASSERT case is therefore compiled out here TOGETHER with
///    its guard, body included, not just its CHECK_THROW.
///  - GISMO_ENSURE (gsDebug.h:120-124) is ungated and fires in EVERY build,
///    Release included, so its case runs unconditionally below.
enum ThrowKind { Threw_none = 0, Threw_logic = 1, Threw_runtime = 2, Threw_other = 3 };

const char * throwKindName(int k)
{
    switch (k)
    {
    case Threw_logic  : return "std::logic_error   (GISMO_ASSERT)";
    case Threw_runtime: return "std::runtime_error (GISMO_ENSURE)";
    case Threw_other  : return "some OTHER exception";
    default           : return "NOTHING -- the guard did not fire";
    }
}

template<class Fn>
int throwKindOf(Fn f)
{
    try                                { f(); }
    catch (const std::logic_error   &) { return Threw_logic;   }
    catch (const std::runtime_error &) { return Threw_runtime; }
    catch (...)                        { return Threw_other;   }
    return Threw_none;
}

// Each case must violate the FIRST guard it reaches, otherwise CHECK_THROW goes
// green while proving something else. The order inside the provider is
// (E = GISMO_ENSURE, live in every build; A = GISMO_ASSERT, dev build only):
//   ensureFilled : :411 material!=nullptr (E) -> :415 cacheValid
//                  -> :423 both maps bound (E)
//   _fill        : :664 NEED_VALUE(ori) (E) -> :671 NEED_VALUE(def) (E)
//                  -> :675/:678 point count(ori/def) (A)
//                  -> :685/:688 patchId(ori/def)     (A)
// so every case below perturbs exactly ONE of them and leaves the rest valid.
TEST(ProtocolGuards)
{
    gsMultiPatch<real_t> mp = makeCurved();
    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u  = interiorPoints();   // 6 points

#ifndef NDEBUG
    gsInfo << "[ProtocolGuards] the five cases below trip provider guards ON "
              "PURPOSE. Each fires its guard TWICE (once through throwKindOf, "
              "once through CHECK_THROW), so TEN Assert/Ensure dumps on stderr "
              "are EXPECTED in a PASSING run.\n";
#else
    gsInfo << "[ProtocolGuards] NDEBUG build: only the two GISMO_ASSERT cases "
              "(d),(e) are compiled out together with the guard they test; the "
              "three GISMO_ENSURE cases (a),(b),(c) and the clean control (f) "
              "run here. Each of the three fires its guard TWICE (throwKindOf + "
              "CHECK_THROW), so SIX Ensure dumps on stderr are EXPECTED.\n";
#endif

    // Cases (a) and (b) check CALLER-SEQUENCING contracts and are GISMO_ENSURE
    // (gsShellMaterialProvider.h:411 and :423), hence live -- and tested -- in
    // EVERY build, Release included. They were GISMO_ASSERT until task 47, and
    // that was not merely a coverage gap: under -DNDEBUG the macro expands to
    // NOTHING (gsDebug.h:105-113), so the guard was absent and ensureFilled()
    // ran on into _fill with a null material (a) or dereferenced null map
    // pointers (b, `_fill(c,*mapOri,*mapDef)` at gsShellMaterialProvider.h:428).
    // The Release run of this suite then died inside case (a) with "Unhandled
    // exception: test crashed", UnitTest++'s report for a translated SIGSEGV
    // (Posix/SignalTranslator.h, ExecuteTest.h:48-51). Promoting the two guards
    // is what turned that crash into the clean throw checked below, so these two
    // cases are no longer wrapped in #ifndef NDEBUG -- and they expect
    // std::runtime_error, the GISMO_ENSURE type (gsDebug.h:120-124).

    // (a) a fill that never saw beginElement() -- gsShellMaterialProvider.h:411
    {
        Provider prov(&law, mp, t, NZ);
        gsMapData<real_t> md;
        computeMapLike(mp.patch(0), u, 0, md);
        prov.bindKinematics(&md,&md);               // bound, but no element yet
        const int kind = throwKindOf([&](){ prov.ensureFilled(); });
        gsInfo << "[ProtocolGuards] (a) no beginElement  -> " << throwKindName(kind) << "\n";
        CHECK_EQUAL((int)Threw_runtime, kind);
        CHECK_THROW(prov.ensureFilled(), std::runtime_error);
    }

    // (b) an UNBOUND provider: beginElement ran, bindKinematics did not -- :423.
    //     This is the stale/missing-bind scenario the views must prevent by
    //     re-binding in EVERY parse.
    {
        Provider prov(&law, mp, t, NZ);
        gsMatrix<real_t> tmp;
        prov.piece(0).eval_into(u,tmp);             // beginElement only
        CHECK(prov.mapOri() == nullptr);
        const int kind = throwKindOf([&](){ prov.ensureFilled(); });
        gsInfo << "[ProtocolGuards] (b) unbound maps    -> " << throwKindName(kind) << "\n";
        CHECK_EQUAL((int)Threw_runtime, kind);
        CHECK_THROW(prov.ensureFilled(), std::runtime_error);
    }

    // (c) a map computed WITHOUT NEED_VALUE -- :664. GISMO_ENSURE, so this one
    //     fires in a true release build as well: it is deliberately NOT wrapped
    //     in #ifndef NDEBUG, and it is the case that keeps this test honest in
    //     the Release build.
    {
        Provider prov(&law, mp, t, NZ);
        gsMapData<real_t> md;
        md.flags  = baseDimFlags() & ~(unsigned)NEED_VALUE;
        md.points = u;
        static_cast<const gsFunction<real_t>&>(mp.patch(0)).computeMap(md);
        md.patchId = 0;
        // computeMap ORs flags IN (gsFunction.hpp:819-822), so check that it did
        // not silently re-add NEED_VALUE -- otherwise this case is vacuous.
        CHECK((md.flags & NEED_VALUE) == 0);
        prov.bindKinematics(&md,&md);
        gsMatrix<real_t> tmp;
        prov.piece(0).eval_into(u,tmp);
        const int kind = throwKindOf([&](){ prov.ensureFilled(); });
        gsInfo << "[ProtocolGuards] (c) no NEED_VALUE   -> " << throwKindName(kind) << "\n";
        CHECK_EQUAL((int)Threw_runtime, kind);
        CHECK_THROW(prov.ensureFilled(), std::runtime_error);
    }

    // (d) and (e) STAY compiled out under NDEBUG, deliberately. Their detectors
    // are still GISMO_ASSERT -- gsShellMaterialProvider.h:675/:678 (point counts)
    // for (d) and :685/:688 (patchId) for (e) -- and they were kept that way on
    // purpose: unlike (a)/(b) these are INTERNAL consistency checks between the
    // bound maps and the cached element, largely the library's own doing, so a
    // permanent per-fill runtime cost was judged not worth it. The bodies must go
    // with the guards: without them _fill would read past the shorter map (d) or
    // silently mix patch data (e), so the CHECK_THROW alone could not be kept.
    // Promoting :675/:678 and :685/:688 to GISMO_ENSURE is the ONLY change that
    // would let this block be unguarded.
#ifndef NDEBUG
    // (d) ori and def computed on DIFFERENT point counts -- :678. Both carry
    //     NEED_VALUE and the right patchId, so the count is the only violation.
    {
        // declared here, not at the top of the test: under NDEBUG this whole
        // block is gone and a file-scope u2 would be an unused variable.
        const gsMatrix<real_t> u2 = otherPoints();  // 4 points
        Provider prov(&law, mp, t, NZ);
        gsMapData<real_t> mdOri, mdDef;
        computeMapLike(mp.patch(0), u , 0, mdOri);   // 6 points
        computeMapLike(mp.patch(0), u2, 0, mdDef);   // 4 points
        CHECK(mdOri.values[0].cols() != mdDef.values[0].cols());
        prov.bindKinematics(&mdOri,&mdDef);
        gsMatrix<real_t> tmp;
        prov.piece(0).eval_into(u,tmp);             // the element cached 6 points
        const int kind = throwKindOf([&](){ prov.ensureFilled(); });
        gsInfo << "[ProtocolGuards] (d) 6 vs 4 points  -> " << throwKindName(kind) << "\n";
        CHECK_EQUAL((int)Threw_logic, kind);
        CHECK_THROW(prov.ensureFilled(), std::logic_error);
    }

    // (e) a bound map whose patchId contradicts the element's patch -- :685.
    //     Same points, same flags: only patchId differs.
    {
        Provider prov(&law, mp, t, NZ);
        gsMapData<real_t> mdOri, mdDef;
        computeMapLike(mp.patch(0), u, 1, mdOri);   // patchId 1 ...
        computeMapLike(mp.patch(0), u, 0, mdDef);
        prov.bindKinematics(&mdOri,&mdDef);
        gsMatrix<real_t> tmp;
        prov.piece(0).eval_into(u,tmp);             // ... but the element is patch 0
        const int kind = throwKindOf([&](){ prov.ensureFilled(); });
        gsInfo << "[ProtocolGuards] (e) patchId 1 vs 0 -> " << throwKindName(kind) << "\n";
        CHECK_EQUAL((int)Threw_logic, kind);
        CHECK_THROW(prov.ensureFilled(), std::logic_error);
    }
#endif // NDEBUG -- (d) and (e) test the GISMO_ASSERTs at :675/:678 and :685/:688

    // Sanity: the very same driver WITHOUT any perturbation must not throw, so
    // the cases above are attributable to their perturbation and not to the
    // harness itself. Runs in every build.
    {
        Provider prov(&law, mp, t, NZ);
        gsMapData<real_t> mdOri, mdDef;
        const int kind = throwKindOf([&](){ fillElement(prov, mp, mp, u, mdOri, mdDef); });
        gsInfo << "[ProtocolGuards] (f) clean control  -> " << throwKindName(kind) << "\n";
        CHECK_EQUAL((int)Threw_none, kind);
        CHECK_EQUAL((size_t)1, prov.fillCount());
    }
}

// =====================================================================
// TESTS 11-16 : the VIEW expressions (task 25)
//
// Everything above drives the PROVIDER directly (hand-built gsMapData,
// bindKinematics/beginElement/ensureFilled). The six tests below drive the six
// shellMaterialView expressions through a REAL gsExprHelper parse/precompute
// cycle -- the same cycle gsExprAssembler runs -- so what is gated here is the
// VIEW layer: the block table, the zero-copy mapping, the shared sweep, the
// request-mask assert, and the copy semantics that make the views usable under
// OpenMP.
// =====================================================================

typedef MaterialOutput MO;

typedef expr::shellMaterialView_expr<MO::MatrixA,3,real_t> ViewA;
typedef expr::shellMaterialView_expr<MO::MatrixB,3,real_t> ViewB;
typedef expr::shellMaterialView_expr<MO::MatrixC,3,real_t> ViewC;
typedef expr::shellMaterialView_expr<MO::MatrixD,3,real_t> ViewD;
typedef expr::shellMaterialView_expr<MO::VectorN,3,real_t> ViewN;
typedef expr::shellMaterialView_expr<MO::VectorM,3,real_t> ViewM;

/// One gsExprHelper plus the three symbols the views are built from. Held
/// together because the maps must outlive every parse that binds them.
struct ViewRig
{
    gsExprHelper<real_t>::uPtr  ev;
    expr::gsGeometryMap<real_t> G;    ///< undeformed midsurface
    expr::gsGeometryMap<real_t> Gd;   ///< deformed midsurface
    expr::gsFeVariable<real_t>  pv;   ///< the provider as a PLAIN variable

    ViewRig(const gsMultiPatch<real_t> & mp, const gsMultiPatch<real_t> & mpDef,
            const Provider & prov)
    :
    ev(gsExprHelper<real_t>::make()),
    G (ev->getMap(mp)),
    Gd(ev->getMap(mpDef)),
    pv(ev->getVar(prov))
    { }
};

/// The six views on one rig. The implicitly generated copy constructor invokes
/// each view's own copy constructor -- exactly what gsExprAssembler::assemble
/// does per thread with std::make_tuple(args...) (gsExprAssembler.h:1119).
struct SixViews
{
    ViewA A; ViewB B; ViewC C; ViewD D; ViewN N; ViewM M;

    SixViews(const ViewRig & r, const Provider * p)
    :
    A(expr::shellMaterialView<MO::MatrixA>(r.pv,r.G,r.Gd,p)),
    B(expr::shellMaterialView<MO::MatrixB>(r.pv,r.G,r.Gd,p)),
    C(expr::shellMaterialView<MO::MatrixC>(r.pv,r.G,r.Gd,p)),
    D(expr::shellMaterialView<MO::MatrixD>(r.pv,r.G,r.Gd,p)),
    N(expr::shellMaterialView<MO::VectorN>(r.pv,r.G,r.Gd,p)),
    M(expr::shellMaterialView<MO::VectorM>(r.pv,r.G,r.Gd,p))
    { }
};

/// Registers the six views AND the rig's own maps. The second half is NOT
/// decoration: parse() binds the view's PRIVATE deep copies of the maps, so a
/// gsGeometryMap that the caller passes in but never parses elsewhere stays
/// unbound and G.data() trips gsGeometryMap.h:63. In a real assembly the
/// integrand parses them anyway; here it is done explicitly.
void parseRig(ViewRig & r, SixViews & v)
{ r.ev->parse(v.A,v.B,v.C,v.D,v.N,v.M,r.G,r.Gd); }

/// Worst RELATIVE deviation of view @a v from the legacy coefficient matrix
/// @a ref (rows*cols x npts), over all quadrature points. No CHECK inside: this
/// is also called from inside an OpenMP region, where UnitTest++'s shared
/// TestResults must not be touched.
template<class E>
real_t viewDevVsLegacy(const E & v, const gsMatrix<real_t> & ref, index_t npts)
{
    const index_t r = v.rows(), c = v.cols();
    real_t w = 0;
    for (index_t k=0; k!=npts; ++k)
    {
        const gsMatrix<real_t> got = v.eval(k);
        const gsMatrix<real_t> rk  = ref.reshapeCol(k,r,c);
        w = math::max(w, (got-rk).norm() / ((real_t)1.0 + rk.norm()));
    }
    return w;
}

/// Same, but also prints the worst ABSOLUTE deviation and whether it is exactly
/// 0 (the format of checkSixBlocks, so the log reads uniformly), and gates on
/// the file's REL_TOL. Returns the worst relative deviation.
template<class E>
real_t compareViewToLegacy(const char * tag, const char * blockName,
                           const E & v, const gsMatrix<real_t> & ref, index_t npts)
{
    const index_t r = v.rows(), c = v.cols();
    // Shapes FIRST: a 0x0 reference would make every norm below 0 and the whole
    // comparison vacuously green.
    CHECK_EQUAL(r*c,  ref.rows());
    CHECK_EQUAL(npts, ref.cols());
    if (ref.rows()!=r*c || ref.cols()!=npts) return (real_t)1.0;

    real_t worstAbs = 0, worstRel = 0;
    for (index_t k=0; k!=npts; ++k)
    {
        const gsMatrix<real_t> got = v.eval(k);
        const gsMatrix<real_t> rk  = ref.reshapeCol(k,r,c);
        const real_t a = (got-rk).norm();
        worstAbs = math::max(worstAbs, a);
        worstRel = math::max(worstRel, a / ((real_t)1.0 + rk.norm()));
    }
    gsInfo << "[" << tag << "] view " << blockName << " vs LEGACY : |ref| = " << ref.norm()
           << " , worst abs dev = " << worstAbs << " , worst rel dev = " << worstRel
           << (worstAbs == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
    CHECK(ref.norm() > (real_t)0.0);      // not a vacuous 0-vs-0 comparison
    CHECK(worstRel <= REL_TOL);
    return worstRel;
}

/// All six views of @a v against the six legacy integrator coefficients.
void compareSixViews(const char * tag, SixViews & v, const LegacySix & L, index_t npts)
{
    compareViewToLegacy(tag,"A",v.A,L.A,npts);
    compareViewToLegacy(tag,"B",v.B,L.B,npts);
    compareViewToLegacy(tag,"C",v.C,L.C,npts);
    compareViewToLegacy(tag,"D",v.D,L.D,npts);
    compareViewToLegacy(tag,"N",v.N,L.N,npts);
    compareViewToLegacy(tag,"M",v.M,L.M,npts);
}

// =====================================================================
// TEST 11 : the six views vs the six legacy integrators, SvK and NH
// =====================================================================
TEST(ViewsMatchLegacy)
{
    // The oracle is the LEGACY gsMaterialMatrixIntegrate path -- code OUTSIDE
    // the provider/view stack, i.e. not the code's own output re-pasted as
    // truth. The gate is the file's REL_TOL; the ACTUAL worst deviation of each
    // block is printed above, and task 23/24 measured exact zeros throughout.
    // A regression to ~1e-12 would still pass this gate but would show up in the
    // log as "(NONZERO)" -- that is deliberate, it is information, and the
    // tolerance is never widened to make a test pass.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    gsFunctionExpr<real_t> t = constFun3(0.01);
    const gsMatrix<real_t> u = interiorPoints();

    gsLinearMaterial<real_t>       svk(E_MOD, NU, 2);
    gsNeoHookeQuadMaterial<real_t> nh (E_MOD, NU, 2);   // C33 Newton inside the condensation
    const gsMaterialBase<real_t> * laws [2] = {&svk, &nh};
    const char *                   names[2] = {"Views_SvK", "Views_NH"};

    for (int m=0; m!=2; ++m)
    {
        const LegacySix L = legacySix(mp, mp_def, t, *laws[m], u);

        Provider prov(laws[m], mp, t, NZ);
        prov.resetCounters();
        ViewRig  rig(mp, mp_def, prov);
        SixViews V(rig, &prov);

        parseRig(rig, V);
        // parse() handed the provider THIS thread's map-data addresses, and the
        // view's private copies resolve to the same per-function-set slot as the
        // rig's own maps.
        CHECK(prov.mapOri() == &rig.G .data());
        CHECK(prov.mapDef() == &rig.Gd.data());
        CHECK((rig.G .data().flags & NEED_VALUE ) != 0);
        CHECK((rig.G .data().flags & NEED_NORMAL) != 0);
        CHECK((rig.G .data().flags & NEED_DERIV2) != 0);
        CHECK((rig.Gd.data().flags & NEED_VALUE ) != 0);

        rig.ev->points() = u;
        rig.ev->precompute(0);
        // The element-boundary signal ran (beginElement), but the sweep is
        // DEFERRED to the first view evaluation.
        CHECK_EQUAL((size_t)0, prov.fillCount());

        compareSixViews(names[m], V, L, u.cols());
        CHECK_EQUAL((size_t)1, prov.fillCount());
    }
}

// =====================================================================
// TEST 12 : the views map the cache IN PLACE (pointer identity)
// =====================================================================
TEST(ViewsAreZeroCopy)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();          // 6 points

    Provider prov(&law, mp, t, NZ);
    prov.resetCounters();
    ViewRig  rig(mp, mp_def, prov);
    SixViews V(rig, &prov);
    parseRig(rig, V);
    rig.ev->points() = u;
    rig.ev->precompute(0);

    // The FIRST evaluation runs the sweep, which RESIZES the 42 x N cache; the
    // base pointer is therefore only stable afterwards.
    const gsAsConstMatrix<real_t> a0 = V.A.eval(0);
    const real_t * base = prov.cachedMoments().data();
    CHECK_EQUAL((size_t)1, prov.fillCount());
    CHECK_EQUAL((index_t)42, prov.cachedMoments().rows());
    CHECK_EQUAL((index_t)6,  prov.cachedMoments().cols());

    // The offsets below are LITERALS ON PURPOSE. Spelling them as
    // Provider::OffsetM / Provider::TargetDim -- which is what
    // shellMomentBlock<> itself does -- would move BOTH sides of the comparison
    // together, so a poisoned block table would sail through. The literals are
    // the independent statement of the layout documented at
    // gsShellMaterialProvider.h:206-217 and in the task spec: stride 42, and
    // A=0, B=9, C=18, D=27, N=36, M=39.
    gsInfo << "[ZeroCopy] cache base = " << (const void*)base
           << " , A.eval(0) = "  << (const void*)a0.data()
           << " , M.eval(3) = "  << (const void*)V.M.eval(3).data() << "\n";
    CHECK(a0.data()            == base + 0*42 +  0);   // A at k = 0
    CHECK(V.M.eval(3).data()   == base + 3*42 + 39);   // M at k = 3
    // Two more, at other columns and other blocks: a UNIFORMLY shifted table
    // would still have to survive all four.
    CHECK(V.D.eval(2).data()   == base + 2*42 + 27);
    CHECK(V.N.eval(1).data()   == base + 1*42 + 36);
    CHECK(V.B.eval(5).data()   == base + 5*42 +  9);
    CHECK(V.C.eval(4).data()   == base + 4*42 + 18);

    // Shapes, and non-vacuity: the mapped memory must actually carry data.
    CHECK_EQUAL((index_t)3, V.A.rows());  CHECK_EQUAL((index_t)3, V.A.cols());
    CHECK_EQUAL((index_t)3, V.N.rows());  CHECK_EQUAL((index_t)1, V.N.cols());
    CHECK(V.A.eval(0).norm() > (real_t)0.0);
    CHECK(V.M.eval(3).norm() > (real_t)0.0);

    // Zero copy also means WRITE-THROUGH visibility: a change of the underlying
    // cache is seen by the view without re-evaluating anything. Re-precomputing
    // a DIFFERENT element and re-reading through the same view object proves the
    // view is a window, not a snapshot.
    const real_t before = V.A.eval(0).norm();
    rig.ev->points() = otherPoints();                 // 4 points
    rig.ev->precompute(0);
    const real_t after = V.A.eval(0).norm();
    gsInfo << "[ZeroCopy] |A(0)| element 1 = " << before << " , element 2 = " << after << "\n";
    CHECK_EQUAL((size_t)2, prov.fillCount());
    CHECK_EQUAL((index_t)4, prov.cachedMoments().cols());
    CHECK(V.A.eval(0).data() == prov.cachedMoments().data());
}

// =====================================================================
// TEST 13 : six coefficients, ONE sweep -- the point of the whole phase
// =====================================================================
TEST(ViewsShareOneSweep)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u  = interiorPoints();     // 6 points
    const gsMatrix<real_t> u2 = otherPoints();        // 4 points

    Provider prov(&law, mp, t, NZ);
    prov.resetCounters();
    ViewRig  rig(mp, mp_def, prov);
    SixViews V(rig, &prov);
    parseRig(rig, V);
    rig.ev->points() = u;
    rig.ev->precompute(0);
    CHECK_EQUAL((size_t)0, prov.fillCount());

    // Six views x six quadrature points = 36 evaluations. The LEGACY route pays
    // one full geometric precomputation per COEFFICIENT (six of them); here all
    // 36 evaluations must share ONE per-element sweep.
    index_t nEval = 0;
    real_t  acc   = 0;
    for (index_t k=0; k!=u.cols(); ++k)
    {
        acc += V.A.eval(k).norm() + V.B.eval(k).norm() + V.C.eval(k).norm()
             + V.D.eval(k).norm() + V.N.eval(k).norm() + V.M.eval(k).norm();
        nEval += 6;
    }
    gsInfo << "[OneSweep] " << nEval << " view evaluations -> fillCount = " << prov.fillCount()
           << " , matrixMomentFills = " << prov.matrixMomentFills()
           << " , stressMomentFills = " << prov.stressMomentFills() << "\n";
    CHECK_EQUAL((index_t)36, nEval);            // the loop really ran 36 times
    CHECK(acc > (real_t)0.0);                   // and returned data, not zeros
    CHECK_EQUAL((size_t)1, prov.fillCount());
    CHECK_EQUAL((size_t)1, prov.matrixMomentFills());
    CHECK_EQUAL((size_t)1, prov.stressMomentFills());

    // A NEW element costs exactly ONE more sweep, no matter how many views read
    // it -- the cache is invalidated per element, not per view.
    rig.ev->points() = u2;
    rig.ev->precompute(0);
    CHECK_EQUAL((size_t)1, prov.fillCount());   // precompute alone still does not fill
    for (index_t k=0; k!=u2.cols(); ++k)
    {
        V.A.eval(k); V.B.eval(k); V.C.eval(k);
        V.D.eval(k); V.N.eval(k); V.M.eval(k);
    }
    gsInfo << "[OneSweep] a second element (24 more evaluations) -> fillCount = "
           << prov.fillCount() << "\n";
    CHECK_EQUAL((size_t)2, prov.fillCount());
    CHECK_EQUAL(u2.cols(), prov.cachedPoints());

    // And the second element is genuinely a different one (else the fill count
    // above would prove nothing about invalidation).
    CHECK(u2.cols() != u.cols());
    const LegacySix L2 = legacySix(mp, mp_def, t, law, u2);
    compareSixViews("OneSweep_elem2", V, L2, u2.cols());
    CHECK_EQUAL((size_t)2, prov.fillCount());   // ... still ONE sweep for element 2
}

// =====================================================================
// TEST 14 : a view whose block is not requested ASSERTS
// =====================================================================
TEST(ViewsRequestMaskAssert)
{
    // Non-requested rows are ZERO, not stale, so a mask/view mismatch would
    // otherwise surface as a silently vanishing stiffness contribution. The
    // GISMO_ASSERT in shellMaterialView_expr::eval (gsShellMaterialExpr.h:263)
    // is the ONLY detector -- it throws std::logic_error, NOT std::runtime_error.
    // Being a GISMO_ASSERT it is live in the dev build (RelWithDebInfo carries no
    // -DNDEBUG) and compiled out under -DNDEBUG (gsDebug.h:105-113), so the
    // negative case is asserted conditionally below. Unlike the provider guards
    // in TEST(ProtocolGuards), the unguarded call is harmless here (the view just
    // returns its -- zeroed -- rows), so only the expectations are switched.
    //
    // This detector STAYS a GISMO_ASSERT even though the provider's own
    // caller-sequencing guards were promoted in task 47: it sits in the per-point
    // hot path (shellMaterialView_expr::eval(k), one call per point of every
    // view), where a permanent branch is exactly what the assembly work of task
    // 42 exists to remove. NB it does NOT depend on the request-mask checks at
    // gsShellMaterialProvider.h:359/:360 (those validate setRequested's argument,
    // not a view/mask mismatch), so promoting those changes nothing here.
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();

#ifndef NDEBUG
    gsInfo << "[RequestMask] the case below trips the view's request-mask assert ON "
              "PURPOSE: one Assert dump on stderr is EXPECTED in a PASSING run.\n";
#else
    gsInfo << "[RequestMask] NDEBUG build: the request-mask GISMO_ASSERT does not "
              "exist here, so NO Assert dump is expected and the mismatch below "
              "goes undetected -- see the #else branch of the check.\n";
#endif

    Provider prov(&law, mp, t, NZ);
    prov.resetCounters();
    prov.setRequested(Provider::ShellReq_N | Provider::ShellReq_M);
    ViewRig  rig(mp, mp_def, prov);
    SixViews V(rig, &prov);
    parseRig(rig, V);
    rig.ev->points() = u;
    rig.ev->precompute(0);

    // POSITIVE CONTROL FIRST: a view whose bit IS in the mask evaluates fine.
    // Without it, the negative case below could be produced by any breakage of
    // the pipeline rather than by the mask.
    const int okKind = throwKindOf([&](){ V.M.eval(0); });
    gsInfo << "[RequestMask] M requested   -> " << throwKindName(okKind) << "\n";
    CHECK_EQUAL((int)Threw_none, okKind);
    CHECK(V.M.eval(0).norm() > (real_t)0.0);
    CHECK_EQUAL((size_t)1, prov.fillCount());

    // A is NOT in the mask.
    CHECK(!prov.requested((unsigned)Provider::ShellReq_A));
    const int badKind = throwKindOf([&](){ V.A.eval(0); });
    gsInfo << "[RequestMask] A NOT requested -> " << throwKindName(badKind) << "\n";
#ifndef NDEBUG // depends on the GISMO_ASSERT at gsShellMaterialExpr.h:263
    CHECK_EQUAL((int)Threw_logic, badKind);
    CHECK_THROW(V.A.eval(0), std::logic_error);
#else
    // -DNDEBUG: the detector does not exist (gsDebug.h:105-113), so the only
    // truthful statement about this build is that the evaluation returns
    // silently. Pinning it keeps the case non-vacuous and documents exactly what
    // Release loses: the mask/view mismatch goes UNDETECTED here.
    CHECK_EQUAL((int)Threw_none, badKind);
#endif

    // Widen the mask and drive a NEW element: the cache is only re-filled at an
    // element boundary, so re-precomputing is what makes the wider mask take
    // effect. The same view object must now evaluate, and be RIGHT.
    prov.setRequested(Provider::ShellReq_All);
    rig.ev->points() = u;
    rig.ev->precompute(0);
    const int wideKind = throwKindOf([&](){ V.A.eval(0); });
    gsInfo << "[RequestMask] A requested   -> " << throwKindName(wideKind) << "\n";
    CHECK_EQUAL((int)Threw_none, wideKind);
    const LegacySix L = legacySix(mp, mp_def, t, law, u);
    compareViewToLegacy("RequestMask","A",V.A,L.A,u.cols());
}

// =====================================================================
// TEST 15 : the COPY CONSTRUCTOR -- what gsExprAssembler runs per thread
// =====================================================================

// gsExprAssembler::assemble builds `std::make_tuple(args...)` INSIDE its
// `#pragma omp parallel` block (gsExprAssembler.h:1117-1120), i.e. every thread
// copy-constructs the expressions and then parses ITS OWN copies. That is the
// only reason per-thread map resolution works at all: gsExprHelper::add stores
// `&m_mdata[fs].mine()`, the address of the CALLING THREAD's slot
// (gsExprHelper.h:353-360 through util::gsThreaded::operator C&). An expression
// that held a map BY REFERENCE would have all threads rebind ONE shared symbol.
//
// This test has two parts, and they prove DIFFERENT things -- see the report:
//   (A) DETERMINISTIC: the view's maps are its OWN objects, not the caller's.
//   (B) STATISTICAL  : under threads, each copy resolves to a DISTINCT slot.
TEST(ViewCopySemanticsPerThread)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();
    const LegacySix L = legacySix(mp, mp_def, t, law, u);

    Provider prov(&law, mp, t, NZ);
    prov.resetCounters();

    // ---- (A) deterministic: distinct map OBJECTS, both resolving to data ----
    ViewRig r1(mp, mp_def, prov);
    ViewA   V (expr::shellMaterialView<MO::MatrixA>(r1.pv, r1.G, r1.Gd, &prov));
    ViewA   Vc(V);                       // <-- THE COPY CONSTRUCTOR under test

    r1.ev->parse(V, r1.G, r1.Gd);
    const gsMapData<real_t> * a1 = prov.mapOri();
    const gsMapData<real_t> * b1 = prov.mapDef();
    CHECK(a1 == &r1.G .data());
    CHECK(b1 == &r1.Gd.data());
    CHECK(a1 != b1);                     // two different function sets, two slots

    // A SECOND helper, and the COPY parsed into it.
    ViewRig r2(mp, mp_def, prov);
    r2.ev->parse(Vc, r2.G, r2.Gd);
    const gsMapData<real_t> * a2 = prov.mapOri();
    const gsMapData<real_t> * b2 = prov.mapDef();
    gsInfo << "[CopySemantics] helper 1 slots = " << (const void*)a1 << " / " << (const void*)b1
           << " , helper 2 slots = "              << (const void*)a2 << " / " << (const void*)b2 << "\n";
    CHECK(a2 == &r2.G .data());
    CHECK(b2 == &r2.Gd.data());
    CHECK(a1 != a2);                     // the two helpers own separate slots
    CHECK(b1 != b2);

    // THE DISCRIMINATOR of part (A): parsing the COPY into helper 2 did NOT drag
    // the original caller's symbols along. Had the view stored `const
    // gsGeometryMap<T>&` members (the design the class doc rules out), Vc would
    // reference r1.G itself and the line below would now read helper 2's slot.
    CHECK(&r1.G .data() == a1);
    CHECK(&r1.Gd.data() == b1);

    // ... and BOTH must resolve to valid data. Helper 2 is the one currently
    // bound, so drive it first.
    r2.ev->points() = u;
    r2.ev->precompute(0);
    compareViewToLegacy("CopySem_h2","A",Vc,L.A,u.cols());

    r1.ev->parse(V, r1.G, r1.Gd);        // re-parse re-binds helper 1's addresses
    r1.ev->points() = u;
    r1.ev->precompute(0);
    compareViewToLegacy("CopySem_h1","A",V,L.A,u.cols());

    // ---- (B) statistical: one copy per thread -> one slot per thread --------
#ifdef _OPENMP
    // gsThreaded sizes its array with omp_get_max_threads() AT CONSTRUCTION, and
    // gsExprHelper::add indexes it with omp_get_thread_num(): asking for more
    // threads than the helper was built for is an out-of-bounds WRITE, not an
    // assertion. Hence the clamp.
    const int nt = math::min(4, omp_get_max_threads());
    gsInfo << "[CopySemantics] omp_get_max_threads() = " << omp_get_max_threads()
           << " , using " << nt << " threads\n";
    if (nt < 2)
        gsInfo << "[CopySemantics] only one thread available: part (B) is vacuous here "
                  "and is SKIPPED (part (A) above still ran).\n";
    else
    {
        const int R = 32;                // rounds: cleanUp() re-allocates the
                                         // slots every round, so each round is a
                                         // fresh draw of the addresses
        ViewRig r3(mp, mp_def, prov);
        ViewA   V3(expr::shellMaterialView<MO::MatrixA>(r3.pv, r3.G, r3.Gd, &prov));

        std::vector<const gsMapData<real_t>*> ori((size_t)R*nt, nullptr);
        std::vector<const gsMapData<real_t>*> def((size_t)R*nt, nullptr);
        std::vector<int> threw((size_t)nt, 0);

        // NO CHECK/CHECK_EQUAL inside this region: UnitTest++ accumulates into a
        // shared TestResults, and an exception escaping an OpenMP structured
        // block is undefined behaviour. Results are collected per thread and
        // asserted after the join.
#       pragma omp parallel num_threads(nt)
        {
            const int tid = omp_get_thread_num();
            for (int r=0; r!=R; ++r)
            {
                try
                {
                    ViewA Vt(V3);                 // the per-thread copy
                    r3.ev->parse(Vt);             // ... parsed by THIS thread only
                    ori[(size_t)r*nt+tid] = prov.mapOri();
                    def[(size_t)r*nt+tid] = prov.mapDef();
                }
                catch (...) { threw[(size_t)tid] = 1; }
                // Mandatory: gsExprHelper::parse begins with cleanUp(), whose
                // `omp single` does NOT stop the executing thread from clearing
                // m_mdata while the others are still reading this round.
#               pragma omp barrier
            }
        }

        int nThrew = 0;
        for (int i=0; i!=nt; ++i) nThrew += threw[(size_t)i];
        CHECK_EQUAL(0, nThrew);

        index_t collisions = 0, nulls = 0, oriDefClash = 0;
        for (int r=0; r!=R; ++r)
            for (int i=0; i!=nt; ++i)
            {
                if (nullptr==ori[(size_t)r*nt+i] || nullptr==def[(size_t)r*nt+i]) ++nulls;
                if (ori[(size_t)r*nt+i] == def[(size_t)r*nt+i]) ++oriDefClash;
                for (int j=i+1; j!=nt; ++j)
                    if (ori[(size_t)r*nt+i] == ori[(size_t)r*nt+j] ||
                        def[(size_t)r*nt+i] == def[(size_t)r*nt+j]) ++collisions;
            }
        gsInfo << "[CopySemantics] " << R << " rounds x " << nt
               << " threads: address collisions = " << collisions
               << " , null bindings = " << nulls
               << " , ori==def clashes = " << oriDefClash << "\n";
        CHECK_EQUAL((index_t)0, collisions);   // every thread bound its OWN slot
        CHECK_EQUAL((index_t)0, nulls);
        CHECK_EQUAL((index_t)0, oriDefClash);
    }
#else
    gsInfo << "[CopySemantics] built without OpenMP: part (B) is not meaningful "
              "(util::gsThreaded degenerates to a single slot) and is SKIPPED.\n";
#endif
}

// =====================================================================
// TEST 16 : the six views under OMP_NUM_THREADS = 4
// =====================================================================
TEST(ViewsUnderOpenMP)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.15);
    gsFunctionExpr<real_t> t = constFun3(0.01);
    gsLinearMaterial<real_t> law(E_MOD, NU, 2);
    const gsMatrix<real_t> u = interiorPoints();
    const LegacySix L = legacySix(mp, mp_def, t, law, u);

    Provider prov(&law, mp, t, NZ);

    // ---- the SINGLE-THREADED reference, on this very provider ---------------
    gsMatrix<real_t> Msingle;
    real_t worstSingle = 0;
    {
        prov.resetCounters();
        ViewRig  rig(mp, mp_def, prov);
        SixViews V(rig, &prov);
        parseRig(rig, V);
        rig.ev->points() = u;
        rig.ev->precompute(0);
        worstSingle = math::max(math::max(viewDevVsLegacy(V.A,L.A,u.cols()),
                                          viewDevVsLegacy(V.B,L.B,u.cols())),
                     math::max(math::max(viewDevVsLegacy(V.C,L.C,u.cols()),
                                          viewDevVsLegacy(V.D,L.D,u.cols())),
                               math::max(viewDevVsLegacy(V.N,L.N,u.cols()),
                                          viewDevVsLegacy(V.M,L.M,u.cols()))));
        Msingle = prov.cachedMoments();          // a COPY of the 42 x N cache
        CHECK_EQUAL((size_t)1, prov.fillCount());
        CHECK(Msingle.norm() > (real_t)0.0);
    }
    gsInfo << "[UnderOMP] single-threaded worst rel dev vs LEGACY = " << worstSingle << "\n";
    CHECK(worstSingle <= REL_TOL);

#ifdef _OPENMP
    const int nt = math::min(4, omp_get_max_threads());
    gsInfo << "[UnderOMP] omp_get_max_threads() = " << omp_get_max_threads()
           << " , using " << nt << " threads\n";
    if (nt < 2)
    {
        gsInfo << "[UnderOMP] only one thread available: the threaded part is vacuous "
                  "here and is SKIPPED.\n";
        return;
    }

    ViewRig  rig(mp, mp_def, prov);
    SixViews V0(rig, &prov);

    std::vector<real_t> worstLegacy((size_t)nt, (real_t)-1.0);
    std::vector<real_t> devVsSingle((size_t)nt, (real_t)-1.0);
    std::vector<size_t> fills      ((size_t)nt, (size_t)0);
    std::vector<int>    threw      ((size_t)nt, 0);

    // Same structure as gsExprAssembler::assemble: copy the expressions per
    // thread, parse the COPIES, then run the element. The rig's own G/Gd are
    // deliberately NOT parsed here -- they are ONE shared symbol, and having
    // every thread rebind it is precisely the race the by-value design exists to
    // avoid. The views carry their own copies, so the maps are computed anyway.
    // No CHECK inside the region (see TEST(ViewCopySemanticsPerThread)).
#   pragma omp parallel num_threads(nt)
    {
        const int tid = omp_get_thread_num();
        SixViews Vt(V0);                        // per-thread copy constructors
        try { rig.ev->parse(Vt.A,Vt.B,Vt.C,Vt.D,Vt.N,Vt.M); }
        catch (...) { threw[(size_t)tid] = 1; }

        // Every thread reaches this: no thread may still be inserting into the
        // shared m_mdata/m_fdata maps while another iterates them in precompute.
#       pragma omp barrier

        if (0 == threw[(size_t)tid])
        {
            try
            {
                rig.ev->points() = u;            // gsThreaded: this thread's points
                rig.ev->precompute(0);
                real_t w = 0;
                w = math::max(w, viewDevVsLegacy(Vt.A,L.A,u.cols()));
                w = math::max(w, viewDevVsLegacy(Vt.B,L.B,u.cols()));
                w = math::max(w, viewDevVsLegacy(Vt.C,L.C,u.cols()));
                w = math::max(w, viewDevVsLegacy(Vt.D,L.D,u.cols()));
                w = math::max(w, viewDevVsLegacy(Vt.N,L.N,u.cols()));
                w = math::max(w, viewDevVsLegacy(Vt.M,L.M,u.cols()));
                worstLegacy[(size_t)tid] = w;
                // The whole per-thread cache against the single-threaded one.
                const gsMatrix<real_t> & Mt = prov.cachedMoments();
                devVsSingle[(size_t)tid] = (Mt.rows()==Msingle.rows() &&
                                            Mt.cols()==Msingle.cols())
                                         ? (Mt - Msingle).norm() : (real_t)1.0;
                fills[(size_t)tid] = 1;          // this thread completed its element
            }
            catch (...) { threw[(size_t)tid] = 2; }
        }
    }

    int nThrew = 0; size_t nRan = 0;
    real_t wLeg = 0, wSing = 0;
    for (int i=0; i!=nt; ++i)
    {
        nThrew += threw[(size_t)i];
        nRan   += fills[(size_t)i];
        gsInfo << "[UnderOMP] thread " << i << " : threw = " << threw[(size_t)i]
               << " , worst rel dev vs LEGACY = " << worstLegacy[(size_t)i]
               << " , |cache - single-threaded| = " << devVsSingle[(size_t)i]
               << (devVsSingle[(size_t)i] == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)")
               << "\n";
        wLeg  = math::max(wLeg,  worstLegacy[(size_t)i]);
        wSing = math::max(wSing, devVsSingle[(size_t)i]);
    }
    CHECK_EQUAL(0, nThrew);                       // no crash, no exception
    CHECK_EQUAL((size_t)nt, nRan);                // every thread really ran
    CHECK(wLeg  <= REL_TOL);                      // every thread got the RIGHT answer
    // Deliberately NOT a bit-identity assertion across thread counts: where
    // gsExprAssembler accumulates a form it does so under `#pragma omp atomic`
    // (gsExprAssembler.h:739/746/759), which is atomic but NOT ordered. Here
    // nothing is accumulated -- each thread runs the same deterministic
    // per-element sweep -- so exact agreement is EXPECTED and reported, but the
    // gate stays relative.
    CHECK(wSing <= REL_TOL * ((real_t)1.0 + Msingle.norm()));
    gsInfo << "[UnderOMP] worst over all threads: rel dev vs LEGACY = " << wLeg
           << " , |cache - single-threaded| = " << wSing
           << (wSing == (real_t)0.0 ? "   (EXACTLY 0)" : "   (NONZERO)") << "\n";
    // fillCount is a plain counter with benign races under OpenMP
    // (gsShellMaterialProvider.h:568), so it is REPORTED, never asserted here.
    gsInfo << "[UnderOMP] fillCount after the threaded run = " << prov.fillCount()
           << " (racy counter, reported only)\n";
#else
    gsInfo << "[UnderOMP] built without OpenMP: the threaded part is SKIPPED.\n";
#endif
}

} // SUITE

#endif // gsPhaseFieldFracture_ENABLED
