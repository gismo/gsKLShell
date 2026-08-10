/** @file gsThinShellAssembler2_test.cpp

    @brief Linear-stiffness / rhs / mass parity gate for gsThinShellAssembler2
           (task 27), and the FIRST execution of the dim==2 code paths.

    gsThinShellAssembler2 (task 26) reproduces gsThinShellAssembler's LINEAR
    assembly verbatim except for ONE structural change: the six
    gsMaterialMatrixIntegrate coefficients are replaced by ONE
    gsShellMaterialProvider plus zero-copy moment views. This suite is the
    numerical gate for that swap.

    ### The two comparison rings (they are NOT the same claim)

    1. SAME-MATERIAL (ring 1), the real gate. gsThinShellAssembler2 driven by a
       PFF law, against the LEGACY gsThinShellAssembler driven by
       gsMaterialMatrix3D<d,real_t> (the task-14 adapter) wrapping the SAME law.
       Both sides then run the same constitutive law, the same
       gsPlaneStressCondensation and the same NumGauss=4 z-Gauss rule, so the
       only thing that differs is the assembler's material plumbing. Agreement
       is EXPECTED to be exact; the gate is 1e-12 relative.

    2. VS PURE LEGACY (ring 2), model-level only. gsThinShellAssembler2 with
       gsLinearMaterial against legacy gsMaterialMatrixLinear. These are
       DIFFERENT MODELS: gsMaterialMatrixLinear reports
       MatIntegration::Constant, so the integrator uses the analytic moments
       (t, 0, t^3/12) of a z-CONSTANT tangent, while the new path
       Gauss-integrates the exact z-dependent metric. The gap is O(t^2/R^2) and
       is PRINTED, never gated tightly (TEST(Linear3D_vs_legacy_model) asserts
       only a loose bound plus the t-halving signature).

    ### Tolerance policy
    CHECK_CLOSE is ABSOLUTE (optional/gsUnitTest/Checks.h), and the quantities
    here range over ~9 orders of magnitude (|K| ~ 1e9, |M| ~ 3e1), so every
    gate is written by hand in the relative form
    (A-B).norm() <= tol*(1+|B|). TOL_EXACT is expressed in units of
    machine epsilon so that a float / multiprecision real_t build scales with
    the arithmetic instead of failing or over-tightening; in double it is
    exactly 1.0e-12. No tolerance was ever widened to make a test pass; the
    measured deviations are printed by every test.

    ### Hazards designed around (all measured in task 26, all documented in
    ###  gsThinShellAssembler2.h, section "Known inherited hazards", :131-204
    ###  -- and all INHERITED, legacy has them too)
      - m_rhs is never cleared, so rhs() can return a stale vector. Every
        scenario below uses a FRESH assembler instance.
      - gsExprAssembler accumulates under #pragma omp atomic, so bit-identity
        across runs is not a thread-safe claim; the one test that needs an EXACT
        integer (materialFills) pins the team size with SingleThreadScope.

    ### Two of those hazards were FIXED by task 44, in BOTH assemblers
      - assembleMass() used to leave the trial space HOMOGENIZED, so that
        assemble() -> assembleMass() -> assemble() silently returned the
        homogenized rhs. It now restores the l2Projection setup before returning.
        TEST(MassRestoresDirichlet) is the direct regression gate;
        TEST(NeumannAndPointLoads) and TEST(Modal_vs_analytic) still apply the old
        remedy (updateBCs) explicitly, which is now redundant but harmless -- and
        keeps asserting that the remedy itself has not broken.
      - assembleMass(lumped=true) did NOT lump. The ledger's "bit-identical to the
        consistent matrix" was a TRUE observation of a FALSE mechanism: the
        rowSum() expression is vector-valued (rowsum_expr.h:32 keeps
        Space = E::Space) and gsExprAssembler dispatched it into the RHS at
        compile time, after which matrix() returned an unmanaged CACHE -- empty,
        sized-but-zero, or STALE, depending on the caller's call history. The
        routine now assembles the consistent matrix unconditionally and row-sums
        it (the gsProjection.hpp:59-69 pattern). TEST(MassLumping_legacy) and
        TEST(MassLumping_new) gate all three of those cache states plus mass
        conservation, and (task 62) the PER-ENTRY identity M_lump(i,i) == the
        i-th row sum of the consistent matrix -- without which "all the mass on
        dof 0" would pass every other gate; TEST(Mass_vs_legacy_varying_t_rho)
        still compares the consistent mass only.

    ### One hazard was GUARDED by task 62, in the legacy assembler
      - a point mass at a parameter point OUTSIDE the domain used to NaN the whole
        mass matrix while assembleMass() reported Success. _applyMass now rejects
        it with a GISMO_ENSURE and the caller sees AssemblyError.
        TEST(MassPointMassOutOfDomain_legacy) gates all four halves of that: the
        rejection, that a point ON the boundary is still accepted, that a point
        inside the guard's tolerance band is CLAMPED onto the domain rather than
        silently losing its mass (task 63), and that a point outside the band is
        still rejected.

    ### And a second one was REPAIRED by task 63, in the legacy assembler
      - _applyMass numbered its actives in the INTEGRATION basis m_basis while
        resolving them through the mapper of the SPACE basis *m_spaceBasis. Under
        setSpaceBasis() the point mass then landed on the wrong dofs, silently,
        under Success. It now dispatches over *m_spaceBasis with the same ladder
        _applyLoads uses. TEST(MassPointMassSpaceBasis_legacy) is the gate, in both
        divergence directions, against a control assembler with no divergence.

    ### Domain conventions of the input functions (three different domains!)
      - the PFF law's own parameters : 2D PARAMETRIC domain (dim 2, always)
      - thickness                    : PHYSICAL, domainDim()==d
      - density                      : PHYSICAL, domainDim()==d

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"

#ifdef gsPhaseFieldFracture_ENABLED

// gsShellMaterialExpr.h (pulled in by the assembler) only forward-declares
// gsExprHelper; the assembler header must come after a real gsExprAssembler
// declaration. Same convention as gsShellMaterialProvider_test.cpp.
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

#include <gsKLShell/src/gsThinShellAssembler2.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrix3D.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/getMaterialMatrix.h>

#include <gsPhaseFieldFracture/materials/gsLinearMaterial.h>
#include <gsPhaseFieldFracture/materials/gsNeoHookeQuadMaterial.h>

#include <limits>
#include <functional>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

SUITE(gsThinShellAssembler2)
{

// =====================================================================
// Shared constants
// =====================================================================

/// Scordelis-Lo roof standard modulus. nu is taken 0.3 (not the textbook 0.0)
/// so that the plane-stress condensation of the 3D law is genuinely active on
/// both sides of ring 1 -- with nu = 0 the C33 condensation is a no-op and the
/// gate would be weaker than it looks.
const real_t E_MOD  = 4.32e8;
const real_t NU     = 0.3;
const real_t THICK  = 0.25;    // ~ 0.01 * L, L ~ 50 (the benchmark's default)
const index_t NREF  = 3;       // 8x8 = 64 elements, 260 dofs (task-26 fixture)

/// Machine-precision gate. Expressed in units of eps so that a float or
/// multiprecision real_t scales with the arithmetic; == 1.0e-12 for double.
const real_t TOL_EXACT = 4.5e3 * std::numeric_limits<real_t>::epsilon();

// =====================================================================
// Small helpers
// =====================================================================

/// Relative Frobenius deviation, robust for near-zero references. The +1 in
/// the denominator is what makes this usable across |K|~1e9 and |M|~3e1.
real_t relDiff(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
{ return (a-b).norm() / (1.0 + b.norm()); }

real_t relDiff(const gsSparseMatrix<real_t> & a, const gsSparseMatrix<real_t> & b)
{
    gsSparseMatrix<real_t> D = a - b;
    return D.norm() / (1.0 + b.norm());
}

/// Absolute Frobenius deviation (so the report can state "exactly 0").
real_t absDiff(const gsSparseMatrix<real_t> & a, const gsSparseMatrix<real_t> & b)
{
    gsSparseMatrix<real_t> D = a - b;
    return D.norm();
}

/// Largest STORED off-diagonal magnitude. A row-summed (lumped) matrix stores
/// diagonal entries ONLY, so this is exactly 0 for it -- structurally, not
/// numerically. Explicit zeros stored off the diagonal would also read 0 here,
/// which is fine: the claim being gated is "no off-diagonal MASS".
real_t maxOffDiag(const gsSparseMatrix<real_t> & M)
{
    real_t m = 0.0;
    for (index_t k = 0; k != M.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(M,k); it; ++it)
            if (it.row() != it.col())
                m = math::max(m, math::abs(it.value()));
    return m;
}

/// Grand sum sum_i sum_j M_ij == the TOTAL MASS carried by the operator. This is
/// the invariant that row-sum lumping conserves by construction, and the one the
/// cheaper rhs() lumping route would NOT conserve on a constrained problem.
real_t grandSum(const gsSparseMatrix<real_t> & M)
{
    real_t s = 0.0;
    for (index_t k = 0; k != M.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(M,k); it; ++it)
            s += it.value();
    return s;
}

/// Row sums rs(i) = sum_j M_ij of a sparse matrix. THE definition of row-sum
/// lumping is the PER-ENTRY identity M_lump(i,i) == rs(i): the weaker gates
/// "diagonal" and "grand sum conserved" do NOT pin WHICH diagonal -- a bug that
/// wrote the whole mass onto dof 0, or that permuted the row sums, satisfies both
/// of them. (Task 44's reviewer, note N2.)
gsMatrix<real_t> consistentRowSums(const gsSparseMatrix<real_t> & M)
{
    gsMatrix<real_t> rs = gsMatrix<real_t>::Zero(M.rows(),1);
    for (index_t k = 0; k != M.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(M,k); it; ++it)
            rs(it.row(),0) += it.value();
    return rs;
}

/// The gsExprAssembler accumulates under #pragma omp atomic: atomic but
/// UNORDERED, so the last bits of two runs of the same form may differ. Only
/// the fill-count invariant needs a pinned team size (the provider's counter is
/// a plain non-atomic member); every other test runs with the ambient team.
struct SingleThreadScope
{
#ifdef _OPENMP
    int m_saved;
    SingleThreadScope() : m_saved(omp_get_max_threads()) { omp_set_num_threads(1); }
    ~SingleThreadScope() { omp_set_num_threads(m_saved); }
#endif
};

// ---------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------

/// The Scordelis-Lo roof from the main filedata: a RATIONAL (NURBS), singly
/// curved, degree-2 surface patch. Curvature is essential here: it is what
/// makes the O(t^2/R^2) model difference of ring 2 measurable at all, and what
/// makes a covariant/contravariant frame flip in the new path visible.
gsMultiPatch<real_t> makeRoof(index_t nref = NREF)
{
    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>("surfaces/scordelis_lo_roof.xml", mp);
    GISMO_ENSURE(mp.nPatches()==1,"Failed to read 'surfaces/scordelis_lo_roof.xml'");
    for (index_t r=0; r!=nref; ++r)
        mp.uniformRefine();
    return mp;
}

/// Scordelis-Lo diaphragm boundary conditions (the fixture of
/// optional/gsKLShell/examples/shell_material_benchmark.cpp:160-170).
/// @a gD is the (scalar) Dirichlet datum; pass nullptr for the homogeneous case.
void roofBCs(const gsMultiPatch<real_t> & mp,
             gsBoundaryConditions<real_t> & bc,
             gsFunctionSet<real_t> * gD = nullptr)
{
    bc.setGeoMap(mp);
    bc.addCondition(boundary::west, condition_type::dirichlet, gD, 0, false, 1);
    bc.addCondition(boundary::west, condition_type::dirichlet, gD, 0, false, 2);
    bc.addCondition(boundary::east, condition_type::dirichlet, gD, 0, false, 1);
    bc.addCondition(boundary::east, condition_type::dirichlet, gD, 0, false, 2);
}

/// A WARPED planar unit square, targetDim 2. A bare BSplineSquare has an
/// IDENTITY Jacobian, which makes every frame convention agree trivially (the
/// task-21 lesson), so the control net is pushed through a smooth non-affine,
/// NON-ORTHOGONAL map. det J stays ~1.05..1.25 over the patch (checked
/// analytically for the underlying warp), so the geometry is non-degenerate.
gsMultiPatch<real_t> makeWarpedSquare2D()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.patch(0).degreeElevate(1);          // biquadratic
    mp.patch(0).uniformRefine();
    mp.patch(0).uniformRefine();           // 4x4 = 16 elements

    gsMatrix<real_t> & c = mp.patch(0).coefs();
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0), Y = c(i,1);
        c(i,0) =        X + 0.25*Y + 0.10*X*Y;
        c(i,1) = 0.20*X + 1.10*Y   - 0.15*X*X;
    }
    GISMO_ENSURE(mp.targetDim()==2 && mp.domainDim()==2,"planar fixture build failed");
    return mp;
}

/// Non-orthogonality / non-identity check on the 2D fixture, so that
/// TEST(Membrane2D_vs_legacy) cannot silently degenerate into the trivial
/// identity-Jacobian case that proves nothing about frame conventions.
void assertWarped(const gsMultiPatch<real_t> & mp)
{
    gsMatrix<real_t> u(2,3);
    u.col(0) << 0.25, 0.25;
    u.col(1) << 0.50, 0.60;
    u.col(2) << 0.80, 0.35;
    gsMatrix<real_t> J;
    real_t worstOffDiag = 0.0, worstDev = 0.0, minDet = 1e30;
    for (index_t k=0; k!=u.cols(); ++k)
    {
        mp.patch(0).jacobian_into(u.col(k), J);
        const real_t g12 = J.col(0).dot(J.col(1));
        worstOffDiag = math::max(worstOffDiag, math::abs(g12));
        worstDev     = math::max(worstDev, (J - gsMatrix<real_t>::Identity(2,2)).norm());
        minDet       = math::min(minDet, J.determinant());
    }
    gsInfo << "  [2D fixture] max |g12| = "<<worstOffDiag
           << " , max |J - I| = "<<worstDev
           << " , min det J = "<<minDet<<"\n";
    CHECK(worstOffDiag > 1e-2);   // genuinely non-orthogonal parametrisation
    CHECK(worstDev     > 1e-1);   // genuinely not the identity map
    CHECK(minDet       > 0.0);    // non-degenerate
}

// =====================================================================
// TEST 1 : ring 1 -- SvK on the roof, stiffness and rhs, machine tight
// =====================================================================
//
// The REAL gate of this task. Both sides run gsLinearMaterial through
// gsPlaneStressCondensation on a NumGauss=4 z-grid; only the assembler-side
// plumbing (one provider + four views vs four gsMaterialMatrixIntegrate
// coefficients) differs.
TEST(Linear3D_vs_adapter)
{
    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);   // PHYSICAL, dim 3

    // The PFF law: its own parameters live on the 2D PARAMETRIC domain.
    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    // NEW path
    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,law);
    // LEGACY path, driven by the task-14 adapter over the SAME law
    gsMaterialMatrix3D<3,real_t> adapter(mp,t,law);
    gsThinShellAssembler<3,real_t,true> AL(mp,dbasis,bc,force,&adapter);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    // Copy OUT: matrix()/rhs() return references into state that a later
    // assembly overwrites (task 26 was bitten by exactly this).
    gsSparseMatrix<real_t> K2 = A2.matrix();
    gsMatrix<real_t>       r2 = A2.rhs();

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsSparseMatrix<real_t> KL = AL.matrix();
    gsMatrix<real_t>       rL = AL.rhs();

    const real_t aK = absDiff(K2,KL), rK = relDiff(K2,KL);
    const real_t aR = (r2-rL).norm(), rR = relDiff(r2,rL);

    gsInfo << "[Linear3D_vs_adapter] #elem = "<<dbasis.totalElements()
           << " , #dofs = "<<A2.numDofs()<<"\n";
    gsInfo << "  |KL| = "<<KL.norm()<<"  abs dev = "<<aK<<"  REL dev = "<<rK
           << (aK==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    gsInfo << "  |rL| = "<<rL.norm()<<"  abs dev = "<<aR<<"  REL dev = "<<rR
           << (aR==0.0 ? "   (EXACTLY 0)" : "")<<"\n";

    // Non-vacuity: a comparison of two zero matrices would pass any tolerance.
    CHECK(KL.norm() > 1.0);
    CHECK(rL.norm() > 1.0);
    CHECK_EQUAL(K2.rows(),KL.rows());
    CHECK_EQUAL(K2.cols(),KL.cols());

    CHECK(rK <= TOL_EXACT);
    CHECK(rR <= TOL_EXACT);
}

// =====================================================================
// TEST 2 : ring 2 -- the O(t^2/R^2) MODEL difference, printed not gated
// =====================================================================
//
// gsMaterialMatrixLinear is MatIntegration::Constant: the integrator applies
// the analytic moments of a z-CONSTANT tangent. gsThinShellAssembler2
// Gauss-integrates the exact z-dependent metric. The residue is a MODEL
// difference of O(t^2/R^2), so the assertions are (a) a loose upper bound,
// (b) a non-triviality floor and (c) the signature that identifies it as a
// model difference rather than a bug: it SHRINKS when the thickness is halved.
real_t modelDiffAtThickness(real_t thick)
{
    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(thick),3);

    // Legacy pure-SvK material: its parameters are PHYSICAL (dim 3).
    gsFunctionExpr<real_t> Ef(util::to_string(E_MOD),3);
    gsFunctionExpr<real_t> nf(util::to_string(NU),3);
    gsFunctionExpr<real_t> rf("1.0",3);
    gsMaterialMatrixLinear<3,real_t> legacyMat(mp,t,Ef,nf,rf);

    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,law);
    gsThinShellAssembler<3,real_t,true>  AL(mp,dbasis,bc,force,&legacyMat);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsSparseMatrix<real_t> K2 = A2.matrix();
    gsMatrix<real_t>       r2 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsSparseMatrix<real_t> KL = AL.matrix();
    gsMatrix<real_t>       rL = AL.rhs();

    // The rhs carries NO material at all, so it must agree EXACTLY even across
    // the two models. This is what tells a genuine model difference apart from
    // a broken assembly.
    const real_t rR = relDiff(r2,rL);
    gsInfo << "  t = "<<thick<<" : rhs REL dev = "<<rR
           << ((r2-rL).norm()==0.0 ? "  (EXACTLY 0)" : "")<<"\n";
    CHECK(rR <= TOL_EXACT);

    return relDiff(K2,KL);
}

TEST(Linear3D_vs_legacy_model)
{
    gsInfo << "[Linear3D_vs_legacy_model] gsThinShellAssembler2(gsLinearMaterial)"
           << " vs legacy(gsMaterialMatrixLinear): MODEL difference, not parity\n";

    const real_t d0 = modelDiffAtThickness(THICK);        // t = 0.25
    const real_t d1 = modelDiffAtThickness(0.5*THICK);    // t = 0.125
    const real_t ratio = d1>0.0 ? d0/d1 : 0.0;

    gsInfo << "  ||K2-KL||/(1+||KL||) at t = "<<THICK      <<" : "<<d0<<"\n";
    gsInfo << "  ||K2-KL||/(1+||KL||) at t = "<<0.5*THICK  <<" : "<<d1<<"\n";
    gsInfo << "  halving ratio d(t)/d(t/2) = "<<ratio<<"   (O(t^2) model signature)\n";

    // (a) LOOSE sanity bound: a model difference, not a broken assembly.
    CHECK(d0 < 1e-2);
    // (b) non-triviality floor: the two models really are different, so this
    //     test can never be mistaken for a parity test that silently passes.
    CHECK(d0 > 1e-6);
    // (c) the model signature: halving the thickness must SHRINK the gap.
    //     O(t^2) predicts ~4; 2.0 is a safe floor that still excludes both a
    //     constant offset (ratio 1) and a growing one.
    CHECK(d1 < d0);
    CHECK(ratio > 2.0);
}

// =====================================================================
// TEST 3 : ring 1 with the pinned NH twin
// =====================================================================
//
// HONEST SCOPE. assemble_impl hardcodes defpatches = m_patches, so the LINEAR
// assembly is always evaluated at ZERO strain (F = I). At F = I the
// quadratic-volumetric Neo-Hooke tangent IS the Lame tangent, so this test does
// NOT discriminate the material model -- what it gates is that the NH law flows
// through the provider/view pipeline (metrics, triad, batched condensation
// including its C33 Newton) identically to the legacy integrator route. The
// nonzero-strain discrimination belongs to assembleVector/assembleMatrix
// (task 28); there is no nonzero-strain mode in this class to test.
// The measured K_NH vs K_SvK gap is printed to document exactly that.
TEST(Linear3D_NH_twin)
{
    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);

    // The task-15 pinned twin of the legacy compressible Neo-Hooke.
    gsNeoHookeQuadMaterial<real_t> law(E_MOD,NU,2);

    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,law);
    gsMaterialMatrix3D<3,real_t> adapter(mp,t,law);
    gsThinShellAssembler<3,real_t,true> AL(mp,dbasis,bc,force,&adapter);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsSparseMatrix<real_t> K2 = A2.matrix();
    gsMatrix<real_t>       r2 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsSparseMatrix<real_t> KL = AL.matrix();
    gsMatrix<real_t>       rL = AL.rhs();

    const real_t aK = absDiff(K2,KL), rK = relDiff(K2,KL);
    const real_t rR = relDiff(r2,rL);

    gsInfo << "[Linear3D_NH_twin] ring 1, gsNeoHookeQuadMaterial on both sides\n";
    gsInfo << "  |KL| = "<<KL.norm()<<"  abs dev = "<<aK<<"  REL dev = "<<rK
           << (aK==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    gsInfo << "  rhs REL dev = "<<rR<<"\n";

    CHECK(KL.norm() > 1.0);
    CHECK(rK <= TOL_EXACT);
    CHECK(rR <= TOL_EXACT);

    // ---- ring 2 (printed): the pinned twin against the LEGACY compressible NH.
    // Both are MatIntegration::NotIntegrated, so this compares the z-integration
    // path of two independent implementations (two independent C33 Newtons).
    gsFunctionExpr<real_t> Ef(util::to_string(E_MOD),3);
    gsFunctionExpr<real_t> nf(util::to_string(NU),3);
    gsFunctionExpr<real_t> rf("1.0",3);
    std::vector<gsFunctionSet<real_t>*> pars(2);
    pars[0] = &Ef; pars[1] = &nf;
    gsOptionList opts;
    opts.addInt   ("Material","Material model",(index_t)Material::NH);
    opts.addSwitch("Compressibility","Compressible",true);
    opts.addInt   ("Implementation","Implementation",1);      // Analytical
    typename gsMaterialMatrixBase<real_t>::uPtr legacyNH =
        getMaterialMatrix<3,real_t>(mp,t,pars,rf,opts);
    GISMO_ENSURE(legacyNH!=nullptr,"getMaterialMatrix returned null");

    gsThinShellAssembler<3,real_t,true> ANH(mp,dbasis,bc,force,legacyNH.get());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ANH.assemble());
    gsSparseMatrix<real_t> KNH = ANH.matrix();
    const real_t rModel = relDiff(K2,KNH);
    gsInfo << "  [ring 2] vs legacy compressible NH : REL dev = "<<rModel<<"\n";
    CHECK(rModel < 1e-4);          // LOOSE: two independent C33 Newtons

    // And the honest scope statement, measured: at zero strain the NH tangent
    // must coincide with the SvK one, so this suite's NH test cannot claim
    // model discrimination.
    gsLinearMaterial<real_t> svk(E_MOD,NU,2);
    gsThinShellAssembler2<3,real_t,true> ASVK(mp,dbasis,bc,force,t,svk);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ASVK.assemble());
    gsSparseMatrix<real_t> KSVK = ASVK.matrix();
    gsInfo << "  [scope] ||K_NH - K_SvK||/(1+||K_SvK||) = "<<relDiff(K2,KSVK)
           << "   (F = I: the linear assembly cannot separate the two models)\n";
}

// =====================================================================
// TEST 4 : mass matrix with NON-CONSTANT thickness AND density
// =====================================================================
//
// The mass involves NO constitutive law: the legacy Density output is exactly
// thickness(x)*rho(x) at PHYSICAL points, and gsThinShellAssembler2 replaces the
// gsMaterialMatrixIntegrate<Density> coefficient by two plain COMPOSITIONS with
// the undeformed map. Constant t and rho would prove nothing about that
// substitution -- a parametric-vs-physical mix-up cancels for constants. On the
// roof x in [0,50] and z in [0,13.5] while the parametric coordinates are in
// [0,1], so an x/z dependence separates the two evaluations by O(1).
TEST(Mass_vs_legacy_varying_t_rho)
{
    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    // PHYSICAL (dim 3), strictly positive, O(1) variation over the patch.
    gsFunctionExpr<real_t> tVar  ("0.25*(1 + 0.4*sin(x/8) + 0.3*cos(z/6))",3);
    gsFunctionExpr<real_t> rhoVar("1.0 + 0.5*cos(x/9) + 0.25*sin(z/5)",3);
    gsFunctionExpr<real_t> tCst  ("0.25",3);
    gsFunctionExpr<real_t> rhoCst("1.0",3);

    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,tVar,law,rhoVar);
    gsMaterialMatrix3D<3,real_t> adapter(mp,tVar,rhoVar,law);
    gsThinShellAssembler<3,real_t,true> AL(mp,dbasis,bc,force,&adapter);

    // The NEW assembler's status IS deterministic (task-26 divergence 2).
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleMass());
    gsSparseMatrix<real_t> M2 = A2.massMatrix();     // copy OUT

    // The LEGACY status used to be an UNINITIALISED READ here: legacy's
    // assembleMass never assigned m_status on success (that assignment sat inside
    // the commented-out block at gsThinShellAssembler.hpp:1468-1484) and no legacy
    // constructor initialises the member either. Observed live before task 44:
    // this test passed in isolation (the value happened to be 0 == Success) and
    // FAILED inside the full suite with "Expected 0 but was -1330503344". Task 44
    // fixed the success path, so the status is asserted again -- but note that
    // status() BEFORE any assembly is still an uninitialised read on legacy
    // (gsThinShellAssembler.h:704 has no default member initialiser), which is why
    // TEST(MassStatus_legacy) drives it from a DETERMINISTIC AssemblyError instead
    // of from a virgin instance.
    const int legacyMassStatus = (int)AL.assembleMass();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,legacyMassStatus);
    gsSparseMatrix<real_t> ML = AL.massMatrix();

    const real_t aM = absDiff(M2,ML), rM = relDiff(M2,ML);
    gsInfo << "[Mass_vs_legacy_varying_t_rho] varying t(x,z) and rho(x,z)\n";
    gsInfo << "  |ML| = "<<ML.norm()<<"  abs dev = "<<aM<<"  REL dev = "<<rM
           << (aM==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    gsInfo << "  [legacy] assembleMass() returned status "<<legacyMassStatus
           << "  (must be Success = "<<(int)ThinShellAssemblerStatus::Success<<")\n";

    CHECK(ML.norm() > 1e-3);
    CHECK(rM <= TOL_EXACT);

    // CONTROL: the same assembly with CONSTANT t and rho. If the variation did
    // not actually reach the integrand, this difference would be ~0 and the
    // test above would be vacuous.
    gsThinShellAssembler2<3,real_t,true> AC(mp,dbasis,bc,force,tCst,law,rhoCst);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AC.assembleMass());
    gsSparseMatrix<real_t> MC = AC.massMatrix();
    const real_t ctrl = relDiff(M2,MC);
    gsInfo << "  [control] varying vs constant t/rho : REL dev = "<<ctrl
           << "   (must be O(1): the variation reaches the integrand)\n";
    CHECK(ctrl > 0.1);
}

// =====================================================================
// TEST 5 : THE FIRST dim==2 EXECUTION IN THIS PROJECT
// =====================================================================
//
// gsShellMaterialProvider<2,real_t>, the dim-2 view layout and the dim-2 branch
// of shellMaterialView_expr::_mapFlags() have only ever been COMPILE-checked
// (task-26 report, coverage note). Ring 1 is available here because BOTH
// gsMaterialMatrix3D<2,real_t> and gsThinShellAssembler<2,real_t,false> are
// explicitly instantiated, so this is a machine-tight gate, not a model
// comparison.
TEST(Membrane2D_vs_legacy)
{
    gsMultiPatch<real_t> mp = makeWarpedSquare2D();
    assertWarped(mp);

    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1);

    // A Neumann side, so that the dim==2 branch of _assembleNeumann_impl does
    // REAL work. Without it assembleBdr() iterates an empty container and the
    // 2D boundary integrand (getBdrFunction(ori) * meas(ori) on a boundary of a
    // planar patch) is entered but never evaluated -- which is exactly where a
    // dim-2 measure/normal convention error would live.
    gsFunctionExpr<real_t> gN2("1.0e2","-2.0e2",2);
    bc.addCondition(boundary::north, condition_type::neumann, &gN2);

    // targetDim must be d == 2. NOTE that for d == 2 the force is ALWAYS taken
    // as parametric: m_parametricForce = (domainDim()==2 && (d==2||d==3)), and
    // for a planar shell the PHYSICAL space is 2D too, so no force function can
    // ever select the m_physforce branch. Legacy has the identical line
    // (gsThinShellAssembler.hpp:50/91), so parity is unaffected.
    gsFunctionExpr<real_t> force("0.0","-1.0e3",2);
    gsFunctionExpr<real_t> t("0.05",2);        // PHYSICAL for d == 2: domainDim 2

    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    gsThinShellAssembler2<2,real_t,false> A2(mp,dbasis,bc,force,t,law);
    gsMaterialMatrix3D<2,real_t> adapter(mp,t,law);
    gsThinShellAssembler<2,real_t,false> AL(mp,dbasis,bc,force,&adapter);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsSparseMatrix<real_t> K2 = A2.matrix();
    gsMatrix<real_t>       r2 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsSparseMatrix<real_t> KL = AL.matrix();
    gsMatrix<real_t>       rL = AL.rhs();

    const real_t aK = absDiff(K2,KL), rK = relDiff(K2,KL);
    const real_t aR = (r2-rL).norm(),  rR = relDiff(r2,rL);

    gsInfo << "[Membrane2D_vs_legacy] FIRST dim==2 execution; #elem = "
           << dbasis.totalElements()<<" , #dofs = "<<A2.numDofs()<<"\n";
    gsInfo << "  |KL| = "<<KL.norm()<<"  abs dev = "<<aK<<"  REL dev = "<<rK
           << (aK==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    gsInfo << "  |rL| = "<<rL.norm()<<"  abs dev = "<<aR<<"  REL dev = "<<rR
           << (aR==0.0 ? "   (EXACTLY 0)" : "")<<"\n";

    CHECK(KL.norm() > 1.0);
    CHECK(rL.norm() > 1e-8);
    CHECK_EQUAL(K2.rows(),KL.rows());
    CHECK(rK <= TOL_EXACT);
    CHECK(rR <= TOL_EXACT);

    // The Neumann side must actually contribute, otherwise the dim==2 boundary
    // integrand above is "covered" only in the sense of being entered.
    gsBoundaryConditions<real_t> bcNoN;
    bcNoN.setGeoMap(mp);
    bcNoN.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0);
    bcNoN.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1);
    gsThinShellAssembler2<2,real_t,false> A2n(mp,dbasis,bcNoN,force,t,law);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2n.assemble());
    gsMatrix<real_t> r2n = A2n.rhs();
    gsInfo << "  [Neumann control] |r without Neumann| = "<<r2n.norm()
           << " ; REL separation = "<<relDiff(r2,r2n)
           << "   (the dim==2 boundary integrand really evaluated)\n";
    CHECK(relDiff(r2,r2n) > 1e-2);
}

// =====================================================================
// TEST 6 : inhomogeneous Dirichlet + Neumann + point load, rhs parity
// =====================================================================
//
// This is where hazard 1 (stale m_rhs) lives, and where the FORMER hazard 2
// (assembleMass homogenizes the space, fixed by task 44) used to. The test is
// structured so that a regression in EITHER makes it FAIL:
//
//  (i)  HOMOGENIZED-SPACE regression. rHom is an explicit control assembled with
//       ZERO Dirichlet data. The test asserts both that r2 == rL (ring 1) and
//       that ||r2 - rHom|| is O(1) relative -- so an rhs that silently drops the
//       Dirichlet lifting cannot pass. After assembleMass() the OLD remedy
//       updateBCs(bc) is applied -- redundant since task 44 made assembleMass
//       restore the space itself, but kept so that the remedy keeps being
//       exercised -- and the rhs is asserted to be BIT-equal to the pre-mass one
//       AND still far from rHom. TEST(MassRestoresDirichlet) is the gate for the
//       fix WITHOUT the remedy.
//  (ii) STALE / ACCUMULATING m_rhs regression. assemble() is called twice on the
//       SAME loaded instance and the two rhs vectors must be IDENTICAL: m_rhs is
//       assigned, never accumulated. A '+=' regression fails here immediately.
//       Additionally a fresh UNLOADED instance must produce a DIFFERENT rhs than
//       the loaded one (so a silently dropped point load cannot pass), and must
//       itself match its legacy twin.
TEST(NeumannAndPointLoads)
{
    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);
    gsFunctionExpr<real_t> rho("1.0",3);
    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    // NON-ZERO Dirichlet datum (scalar, PHYSICAL) + a Neumann traction.
    gsFunctionExpr<real_t> gD("0.02",3);
    gsFunctionExpr<real_t> gN("1.0e2","-2.0e2","3.0e2",3);

    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc,&gD);
    bc.addCondition(boundary::north, condition_type::neumann, &gN);

    // The homogeneous control: identical structure, ZERO Dirichlet data.
    gsBoundaryConditions<real_t> bcH;
    roofBCs(mp,bcH,nullptr);
    bcH.addCondition(boundary::north, condition_type::neumann, &gN);

    gsPointLoads<real_t> pLoads;
    gsVector<real_t> pt(2); pt << 0.5, 0.5;
    gsVector<real_t> pv(3); pv << 0.0, 0.0, -1.0e6;
    pLoads.addLoad(pt,pv,0,true);

    gsMaterialMatrix3D<3,real_t> adapter(mp,t,rho,law);

    // ---- ring 1, LOADED --------------------------------------------------
    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,law,rho);
    A2.setPointLoads(pLoads);
    gsThinShellAssembler<3,real_t,true>  AL(mp,dbasis,bc,force,&adapter);
    AL.setPointLoads(pLoads);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsMatrix<real_t> r2 = A2.rhs();                       // copy OUT
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsMatrix<real_t> rL = AL.rhs();

    const real_t rR = relDiff(r2,rL);
    gsInfo << "[NeumannAndPointLoads] |rL| = "<<rL.norm()
           << "  abs dev = "<<(r2-rL).norm()<<"  REL dev = "<<rR
           << ((r2-rL).norm()==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    CHECK(rL.norm() > 1.0);
    CHECK(rR <= TOL_EXACT);

    // ---- (ii) m_rhs is ASSIGNED, not accumulated -------------------------
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsMatrix<real_t> r2b = A2.rhs();
    gsInfo << "  [reassign] ||rhs(2nd assemble) - rhs(1st)|| = "<<(r2b-r2).norm()<<"\n";
    CHECK(relDiff(r2b,r2) <= TOL_EXACT);

    // ---- (ii) the point load really contributes --------------------------
    gsThinShellAssembler2<3,real_t,true> A2n(mp,dbasis,bc,force,t,law,rho);  // FRESH, no loads
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2n.assemble());
    gsMatrix<real_t> r2n = A2n.rhs();
    gsThinShellAssembler<3,real_t,true> ALn(mp,dbasis,bc,force,&adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALn.assemble());
    gsMatrix<real_t> rLn = ALn.rhs();
    gsInfo << "  [no loads] |rLn| = "<<rLn.norm()<<"  REL dev vs new = "<<relDiff(r2n,rLn)
           << " ; ||r_loaded - r_unloaded||/(1+|r_unloaded|) = "<<relDiff(r2,r2n)<<"\n";
    CHECK(relDiff(r2n,rLn) <= TOL_EXACT);
    CHECK(relDiff(r2,r2n)  > 1e-2);       // the load is genuinely in the rhs

    // ---- (i) the rhs is the NON-homogenized one --------------------------
    gsThinShellAssembler2<3,real_t,true> A2h(mp,dbasis,bcH,force,t,law,rho); // FRESH
    A2h.setPointLoads(pLoads);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2h.assemble());
    gsMatrix<real_t> rHom = A2h.rhs();
    gsInfo << "  [homogeneous control] |rHom| = "<<rHom.norm()
           << " vs |r2| = "<<r2.norm()<<" ; REL separation = "<<relDiff(r2,rHom)<<"\n";
    CHECK(relDiff(r2,rHom) > 1e-2);       // the Dirichlet lifting is genuinely in the rhs

    // ---- (i) assembleMass() + the documented remedy ----------------------
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleMass());
    A2.updateBCs(bc);                     // the verified remedy (task-26 review)
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsMatrix<real_t> r2c = A2.rhs();
    gsInfo << "  [after mass + updateBCs] |r2c| = "<<r2c.norm()
           << " ; REL dev vs r2 = "<<relDiff(r2c,r2)
           << " ; REL separation from rHom = "<<relDiff(r2c,rHom)<<"\n";
    CHECK(relDiff(r2c,r2)   <= TOL_EXACT); // NOT the homogenized rhs
    CHECK(relDiff(r2c,rHom) > 1e-2);
}

// =====================================================================
// TEST 7 : the central claim of the phase, at assembler level
// =====================================================================
//
// ONE per-element constitutive sweep per element, instead of the legacy six
// gsMaterialMatrixIntegrate coefficients each redoing the full geometric
// precomputation. The provider's fill counter is a plain non-atomic member, so
// the EXACT-integer claim is a single-threaded claim and the team size is
// pinned here.
TEST(MaterialFillsOncePerElement)
{
    SingleThreadScope pin;
    GISMO_UNUSED(pin);   // the struct is EMPTY without _OPENMP

    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);
    gsFunctionExpr<real_t> rho("1.0",3);
    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    const size_t nElements = (size_t)dbasis.totalElements();
    CHECK(nElements > 1);

    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,law,rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsInfo << "[MaterialFillsOncePerElement] #elem = "<<nElements
           << " , materialFills() after 1 assemble = "<<A2.materialFills()<<"\n";
    CHECK_EQUAL(nElements, A2.materialFills());

    // The counter ACCUMULATES over assemblies (it is the benchmark observable),
    // and assembleMass() uses no material at all, so it must add nothing.
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleMass());
    gsInfo << "  after assembleMass() : "<<A2.materialFills()<<" (mass uses no material)\n";
    CHECK_EQUAL(nElements, A2.materialFills());

    A2.updateBCs(bc);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsInfo << "  after a 2nd assemble()  : "<<A2.materialFills()<<"\n";
    CHECK_EQUAL(2*nElements, A2.materialFills());

    A2.resetMaterialFills();
    CHECK_EQUAL((size_t)0, A2.materialFills());

    // The MEMBRANE branch asks the provider for ShellReq_A only; the fill count
    // per element must be the same ONE sweep.
    gsThinShellAssembler2<3,real_t,false> AM(mp,dbasis,bc,force,t,law,rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AM.assemble());
    gsInfo << "  membrane <3,real_t,false> : "<<AM.materialFills()<<"\n";
    CHECK_EQUAL(nElements, AM.materialFills());
}

// #####################################################################
// ##  TASK 29 : the NONLINEAR pipeline                               ##
// #####################################################################
//
// Everything below gates the DEFORMED-configuration routines
// (assembleMatrix / assembleVector / the follower pressure) and the Newton
// loops built on them. Two kinds of claim, never mixed:
//
//  PARITY (tests 8-10, 15): ring 1 exactly as above -- gsThinShellAssembler2
//      against the LEGACY gsThinShellAssembler driven by gsMaterialMatrix3D
//      over the SAME law, on the SAME deformed configuration. Gate: TOL_EXACT.
//      Each parity test carries a DISCRIMINATOR assertion measured on the
//      LEGACY side alone (the two laws must genuinely disagree there, or the
//      pressure must genuinely move the operator), so a degenerate comparison
//      cannot pass silently.
//
//  ANALYTIC (tests 11-13): a full Newton loop, or a modal analysis, against a
//      CLOSED-FORM solution. No legacy object is involved; the reference is
//      solid mechanics. These are what make the pipeline right rather than
//      merely equal to something else.
//
// ### HARD SEQUENCING CONSTRAINT (measured in task 28, INHERITED from legacy)
// assembleVector must NOT be the first assembly on a fresh instance:
// initVector(1) does not size the system matrix, gsExprAssembler::assemble
// asserts m_fmatrix.cols()==numDofs() (gsExprAssembler.h:1107), and the
// routine's catch turns that into ThinShellAssemblerStatus::AssemblyError. A
// fresh LEGACY assembler behaves identically, so this is not a divergence and
// it is deliberately NOT asserted here (asserting it would pin an inherited
// wart). Every residual assembly below is preceded by assemble() or
// assembleMatrix() on the SAME instance; thereafter {assembleMatrix,
// assembleVector} may be called in any order -- cleanUp() touches neither
// m_fmatrix nor m_matrix, so the tangent survives a residual assembly.
//
// ### NOT a global claim
// matrixMomentFills()/stressMomentFills() are fed ONLY by assembleMatrix and
// assembleVector; the linear assemble() contributes nothing to them (it was
// left byte-identical for the parity suite above). Test 14 therefore resets
// the counters after its priming assemble().

/// The deformed configuration used by every nonlinear parity test: the control
/// net scaled by @a f, exactly as optional/gsKLShell/examples/
/// shell_material_benchmark.cpp:157 does. On the roof this is a ~1% inflation,
/// which is far outside the linear regime of the constitutive comparison: it
/// separates gsLinearMaterial from gsNeoHookeQuadMaterial by ~6% of ||K||
/// (the LEGACY-side discriminator asserted below), whereas at F = I the two
/// coincide to machine precision (TEST(Linear3D_NH_twin), finding 2 of task 27).
gsMultiPatch<real_t> scaleCoefs(const gsMultiPatch<real_t> & mp, const real_t f)
{
    gsMultiPatch<real_t> def = mp;
    for (size_t p=0; p!=def.nPatches(); ++p)
        def.patch(p).coefs() *= f;
    return def;
}

/// Prints a parity line in the format used throughout this file and returns the
/// relative deviation, so a test reads as "measure, print, gate".
real_t reportParity(const std::string & tag,
                    const gsSparseMatrix<real_t> & a, const gsSparseMatrix<real_t> & b)
{
    const real_t ad = absDiff(a,b), rd = relDiff(a,b);
    gsInfo << "  "<<tag<<" : |ref| = "<<b.norm()<<"  abs dev = "<<ad
           << "  REL dev = "<<rd << (ad==0.0 ? "   (EXACTLY 0)" : "") <<"\n";
    return rd;
}

real_t reportParity(const std::string & tag,
                    const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
{
    const real_t ad = (a-b).norm(), rd = relDiff(a,b);
    gsInfo << "  "<<tag<<" : |ref| = "<<b.norm()<<"  abs dev = "<<ad
           << "  REL dev = "<<rd << (ad==0.0 ? "   (EXACTLY 0)" : "") <<"\n";
    return rd;
}

// =====================================================================
// TEST 8 : ring 1 for assembleMatrix on a DEFORMED configuration
// =====================================================================
//
// The task-27 suite could only ever evaluate the material at F = I (the linear
// assemble hardcodes defpatches = m_patches), where the NH and the SvK tangents
// are the SAME matrix -- so it gated the plumbing, not the law. Here the
// configuration is deformed, the two laws genuinely disagree, and the parity
// claim becomes a claim about the constitutive path as well.
TEST(NonlinearMatrix_vs_adapter)
{
    gsMultiPatch<real_t> mp     = makeRoof();
    gsMultiPatch<real_t> mp_def = scaleCoefs(mp,1.01);
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);

    gsLinearMaterial<real_t>       svk(E_MOD,NU,2);
    gsNeoHookeQuadMaterial<real_t> nh (E_MOD,NU,2);

    gsInfo << "[NonlinearMatrix_vs_adapter] assembleMatrix at 1.01x coefs, #elem = "
           << dbasis.totalElements()<<"\n";

    // ---- SvK -------------------------------------------------------------
    gsThinShellAssembler2<3,real_t,true> A2s(mp,dbasis,bc,force,t,svk);
    gsMaterialMatrix3D<3,real_t> adS(mp,t,svk);
    gsThinShellAssembler<3,real_t,true>  ALs(mp,dbasis,bc,force,&adS);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2s.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> K2s = A2s.matrix();          // copy OUT
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALs.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> KLs = ALs.matrix();

    // ---- NH (the task-15 pinned twin of the legacy compressible Neo-Hooke) --
    gsThinShellAssembler2<3,real_t,true> A2n(mp,dbasis,bc,force,t,nh);
    gsMaterialMatrix3D<3,real_t> adN(mp,t,nh);
    gsThinShellAssembler<3,real_t,true>  ALn(mp,dbasis,bc,force,&adN);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2n.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> K2n = A2n.matrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALn.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> KLn = ALn.matrix();

    const real_t rS = reportParity("SvK matrix",K2s,KLs);
    const real_t rN = reportParity("NH  matrix",K2n,KLn);

    CHECK(KLs.norm() > 1.0);
    CHECK(KLn.norm() > 1.0);
    CHECK_EQUAL(K2s.rows(),KLs.rows());
    CHECK(rS <= TOL_EXACT);
    CHECK(rN <= TOL_EXACT);

    // ---- DISCRIMINATOR, measured on the LEGACY side ALONE -----------------
    // Without this the two parity gates could both be satisfied by a path that
    // silently ignores the law. The legacy assembler itself must separate the
    // two materials by a healthy margin at this configuration.
    const real_t disc = relDiff(KLn,KLs);
    gsInfo << "  [discriminator] legacy ||K_NH - K_SvK||/(1+||K_SvK||) = "<<disc
           << "   (|K_SvK| = "<<KLs.norm()<<" , |K_NH| = "<<KLn.norm()<<")\n";
    CHECK(disc > 1e-2);          // measured 0.0591 (task-28 review, [P5])
}

// =====================================================================
// TEST 9 : ring 1 for assembleVector on a DEFORMED configuration
// =====================================================================
TEST(NonlinearVector_vs_adapter)
{
    gsMultiPatch<real_t> mp     = makeRoof();
    gsMultiPatch<real_t> mp_def = scaleCoefs(mp,1.01);
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);

    gsLinearMaterial<real_t>       svk(E_MOD,NU,2);
    gsNeoHookeQuadMaterial<real_t> nh (E_MOD,NU,2);

    gsInfo << "[NonlinearVector_vs_adapter] assembleVector at 1.01x coefs\n";

    gsThinShellAssembler2<3,real_t,true> A2s(mp,dbasis,bc,force,t,svk);
    gsMaterialMatrix3D<3,real_t> adS(mp,t,svk);
    gsThinShellAssembler<3,real_t,true>  ALs(mp,dbasis,bc,force,&adS);

    // SEQUENCING: size the system first (see the block comment above). Both
    // sides need it -- the precondition is inherited, not introduced.
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2s.assemble());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALs.assemble());

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2s.assembleVector(mp_def));
    gsMatrix<real_t> r2s = A2s.rhs();                   // copy OUT
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALs.assembleVector(mp_def));
    gsMatrix<real_t> rLs = ALs.rhs();

    gsThinShellAssembler2<3,real_t,true> A2n(mp,dbasis,bc,force,t,nh);
    gsMaterialMatrix3D<3,real_t> adN(mp,t,nh);
    gsThinShellAssembler<3,real_t,true>  ALn(mp,dbasis,bc,force,&adN);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2n.assemble());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALn.assemble());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2n.assembleVector(mp_def));
    gsMatrix<real_t> r2n = A2n.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALn.assembleVector(mp_def));
    gsMatrix<real_t> rLn = ALn.rhs();

    const real_t rS = reportParity("SvK residual",r2s,rLs);
    const real_t rN = reportParity("NH  residual",r2n,rLn);

    CHECK(rLs.norm() > 1.0);
    CHECK(rLn.norm() > 1.0);
    CHECK(rS <= TOL_EXACT);
    CHECK(rN <= TOL_EXACT);

    // ---- DISCRIMINATOR on the LEGACY side ---------------------------------
    const real_t disc = relDiff(rLn,rLs);
    gsInfo << "  [discriminator] legacy ||R_NH - R_SvK||/(1+||R_SvK||) = "<<disc<<"\n";
    CHECK(disc > 1e-2);          // measured 0.0290 for the same quantity in task 28

    // ---- CONTROL: the DEFORMED map really reaches the integrand -----------
    // The internal force vanishes identically at F = I, so a residual assembled
    // on the UNDEFORMED configuration must be the pure external force, and it
    // must be far from the deformed one. This separates "the law is wrong" from
    // "the deformed geometry never arrived".
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2s.assembleVector(mp));
    gsMatrix<real_t> r0 = A2s.rhs();
    gsInfo << "  [control] |R(F=I)| = "<<r0.norm()<<" vs |R(1.01x)| = "<<r2s.norm()
           << " ; REL separation = "<<relDiff(r2s,r0)<<"\n";
    CHECK(relDiff(r2s,r0) > 1.0);
}

// =====================================================================
// TEST 10 : the FOLLOWER PRESSURE, matrix and vector
// =====================================================================
//
// p * space * sn(m_def).normalized() * meas(m_ori): the pressure follows the
// DEFORMED normal but is integrated over the UNDEFORMED measure, so the true
// pressure it represents is p * A_0/A_def (this is exactly what
// TEST(Balloon_newton_analytic) has to undo). The term is material-independent,
// so the discriminator here cannot be an NH-vs-SvK one: it is that switching
// the pressure on genuinely MOVES both operators.
TEST(FollowerPressure_vs_adapter)
{
    gsMultiPatch<real_t> mp     = makeRoof();
    gsMultiPatch<real_t> mp_def = scaleCoefs(mp,1.01);
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);
    gsConstantFunction<real_t> pressFun(1.0e4,3);
    gsLinearMaterial<real_t> svk(E_MOD,NU,2);

    gsInfo << "[FollowerPressure_vs_adapter] setPressure(1e4) on both assemblers\n";

    gsMaterialMatrix3D<3,real_t> adapter(mp,t,svk);

    // ---- WITH pressure ----------------------------------------------------
    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,svk);
    gsThinShellAssembler<3,real_t,true>  AL(mp,dbasis,bc,force,&adapter);
    A2.setPressure(pressFun);
    AL.setPressure(pressFun);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> K2 = A2.matrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> KL = AL.matrix();

    // the matrix assembly above has sized the system, so the residual may follow
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleVector(mp_def));
    gsMatrix<real_t> r2 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleVector(mp_def));
    gsMatrix<real_t> rL = AL.rhs();

    const real_t rK = reportParity("pressure matrix  ",K2,KL);
    const real_t rR = reportParity("pressure residual",r2,rL);
    CHECK(KL.norm() > 1.0);
    CHECK(rL.norm() > 1.0);
    CHECK(rK <= TOL_EXACT);
    CHECK(rR <= TOL_EXACT);

    // ---- DISCRIMINATOR: the pressure term is LIVE on both sides ------------
    // FRESH instances without setPressure (m_pressInd has no setter back to
    // false, and m_rhs is never cleared -- hazard 1).
    gsThinShellAssembler2<3,real_t,true> A20(mp,dbasis,bc,force,t,svk);
    gsThinShellAssembler<3,real_t,true>  AL0(mp,dbasis,bc,force,&adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A20.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> K20 = A20.matrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A20.assembleVector(mp_def));
    gsMatrix<real_t> r20 = A20.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL0.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> KL0 = AL0.matrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL0.assembleVector(mp_def));
    gsMatrix<real_t> rL0 = AL0.rhs();

    const real_t dK2 = relDiff(K2,K20), dKL = relDiff(KL,KL0);
    const real_t dr2 = relDiff(r2,r20), drL = relDiff(rL,rL0);
    gsInfo << "  [discriminator] new    : ||K(p)-K(0)||/(1+|K|) = "<<dK2
           << " ; ||R(p)-R(0)||/(1+|R|) = "<<dr2<<"\n";
    gsInfo << "  [discriminator] legacy : ||K(p)-K(0)||/(1+|K|) = "<<dKL
           << " ; ||R(p)-R(0)||/(1+|R|) = "<<drL<<"\n";
    // The stiffness contribution of the pressure is a SMALL fraction of ||K||
    // (measured 1.4e-4 on this fixture -- it is a geometric term only), the
    // residual contribution is large (measured 4.1e-2). The floors are set per
    // leg from those measurements; a shared 1e-2 floor would be wrong for the
    // matrix leg.
    CHECK(dK2 > 1e-5);
    CHECK(dKL > 1e-5);
    CHECK(dr2 > 1e-2);
    CHECK(drL > 1e-2);
    // and the two sides must move by the SAME amount
    CHECK(relDiff(K20,KL0) <= TOL_EXACT);
    CHECK(relDiff(r20,rL0) <= TOL_EXACT);
}

// =====================================================================
// TEST 14 : the residual path does NO matrix-moment work
// =====================================================================
//
// The efficiency claim of the phase, made permanent. assembleVector requests
// only ShellReq_N|ShellReq_M (bending) resp. ShellReq_N (membrane), so the
// provider must skip the four tangent moments A/B/C/D entirely -- one
// constitutive sweep per element either way, but a strictly cheaper one.
// The counters are plain non-atomic members, so this is a single-threaded claim.
TEST(VectorPathSkipsMatrixMoments)
{
    SingleThreadScope pin;
    GISMO_UNUSED(pin);

    gsMultiPatch<real_t> mp     = makeRoof();
    gsMultiPatch<real_t> mp_def = scaleCoefs(mp,1.01);
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);
    gsLinearMaterial<real_t> svk(E_MOD,NU,2);

    const size_t nElements = (size_t)dbasis.totalElements();
    CHECK(nElements > 1);

    // ---- the RESIDUAL path, bending <3,true> ------------------------------
    gsThinShellAssembler2<3,real_t,true> AV(mp,dbasis,bc,force,t,svk);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AV.assemble()); // sizes the system
    AV.resetMaterialFills();                 // the linear assemble feeds only materialFills
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AV.assembleVector(mp_def));
    gsInfo << "[VectorPathSkipsMatrixMoments] #elem = "<<nElements
           << " ; after assembleVector : fills = "<<AV.materialFills()
           << " , matrixMoments = "<<AV.matrixMomentFills()
           << " , stressMoments = "<<AV.stressMomentFills()<<"\n";
    CHECK(AV.rhs().norm() > 1.0);            // the routine did real work
    CHECK_EQUAL(nElements, AV.materialFills());
    CHECK_EQUAL((size_t)0,  AV.matrixMomentFills());   // <-- THE CLAIM
    CHECK_EQUAL(nElements, AV.stressMomentFills());

    // ---- POSITIVE CONTROL: the same instance, the MATRIX path --------------
    // Without this the claim above would also hold for a provider that never
    // integrates anything at all.
    AV.resetMaterialFills();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AV.assembleMatrix(mp_def));
    gsInfo << "  [positive control] after assembleMatrix : fills = "<<AV.materialFills()
           << " , matrixMoments = "<<AV.matrixMomentFills()
           << " , stressMoments = "<<AV.stressMomentFills()<<"\n";
    CHECK_EQUAL(nElements, AV.materialFills());
    CHECK_EQUAL(nElements, AV.matrixMomentFills());
    CHECK_EQUAL(nElements, AV.stressMomentFills());

    // ---- the MEMBRANE branch <3,false> asks for ShellReq_N only ------------
    gsThinShellAssembler2<3,real_t,false> AM(mp,dbasis,bc,force,t,svk);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AM.assemble());
    AM.resetMaterialFills();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AM.assembleVector(mp_def));
    gsInfo << "  [membrane <3,false>] fills = "<<AM.materialFills()
           << " , matrixMoments = "<<AM.matrixMomentFills()
           << " , stressMoments = "<<AM.stressMomentFills()<<"\n";
    CHECK_EQUAL(nElements, AM.materialFills());
    CHECK_EQUAL((size_t)0,  AM.matrixMomentFills());
    CHECK_EQUAL(nElements, AM.stressMomentFills());
}

// =====================================================================
// TEST 15 : dim == 2 with a VARYING thickness, linear AND nonlinear
// =====================================================================
//
// Two holes closed at once (both found by the task-27/28 reviews):
//
//  (a) The PHYSICAL-vs-PARAMETRIC thickness substitution was never exercised at
//      dim 2 -- TEST(Membrane2D_vs_legacy) uses the constant "0.05", and a
//      constant cancels the distinction exactly. On the warped square the
//      physical and parametric coordinates differ by O(1) (parametric
//      (0.5,0.6) maps to physical ~(0.68,0.72)), so a mix-up between the two
//      breaks ring-1 parity by O(1). gsShellMaterialProvider evaluates the
//      thickness at mapOri.values[0], i.e. PHYSICAL points, and hard-asserts
//      thickness->domainDim()==d.
//
//  (b) <2,real_t,false> had no NONLINEAR execution at all (task-28 report, §8:
//      compile-only). The deformed legs below are the first ones.
TEST(Membrane2D_varying_thickness)
{
    gsMultiPatch<real_t> mp     = makeWarpedSquare2D();
    gsMultiPatch<real_t> mp_def = scaleCoefs(mp,1.01);
    assertWarped(mp);

    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1);

    gsFunctionExpr<real_t> force("0.0","-1.0e3",2);
    // PHYSICAL (domainDim == d == 2), strictly positive: 0.05*[0.5 .. 1.5].
    gsFunctionExpr<real_t> tVar("0.05*(1 + 0.3*sin(2*x) + 0.2*cos(3*y))",2);
    gsFunctionExpr<real_t> tCst("0.05",2);
    gsLinearMaterial<real_t> law(E_MOD,NU,2);

    gsInfo << "[Membrane2D_varying_thickness] warped planar square, #elem = "
           << dbasis.totalElements()<<"\n";

    gsThinShellAssembler2<2,real_t,false> A2(mp,dbasis,bc,force,tVar,law);
    gsMaterialMatrix3D<2,real_t> adapter(mp,tVar,law);
    gsThinShellAssembler<2,real_t,false>  AL(mp,dbasis,bc,force,&adapter);

    // ---- LINEAR leg -------------------------------------------------------
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsSparseMatrix<real_t> K2 = A2.matrix();
    gsMatrix<real_t>       r2 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsSparseMatrix<real_t> KL = AL.matrix();
    gsMatrix<real_t>       rL = AL.rhs();

    const real_t rK = reportParity("linear K   ",K2,KL);
    const real_t rR = reportParity("linear rhs ",r2,rL);
    CHECK(KL.norm() > 1.0);
    CHECK(rK <= TOL_EXACT);
    CHECK(rR <= TOL_EXACT);

    // ---- CONTROL: the thickness VARIATION reaches the integrand ------------
    // Without this the parity above would be "covered" only in the sense of
    // being entered -- exactly the dim-2 Neumann hole task 27 found.
    gsThinShellAssembler2<2,real_t,false> AC(mp,dbasis,bc,force,tCst,law);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AC.assemble());
    gsSparseMatrix<real_t> KC = AC.matrix();
    const real_t ctrl = relDiff(K2,KC);
    gsInfo << "  [control] varying vs constant t : REL dev = "<<ctrl
           << "   (must be O(1): the variation reaches the integrand)\n";
    CHECK(ctrl > 0.1);

    // ---- NONLINEAR legs: the FIRST <2,real_t,false> deformed execution -----
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> K2d = A2.matrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleMatrix(mp_def));
    gsSparseMatrix<real_t> KLd = AL.matrix();
    // the matrix assembly sized the system on both, so the residuals may follow
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleVector(mp_def));
    gsMatrix<real_t> r2d = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleVector(mp_def));
    gsMatrix<real_t> rLd = AL.rhs();

    const real_t rKd = reportParity("deformed K  ",K2d,KLd);
    const real_t rRd = reportParity("deformed rhs",r2d,rLd);
    CHECK(rKd <= TOL_EXACT);
    CHECK(rRd <= TOL_EXACT);
    // and the deformed configuration is genuinely a different operator
    gsInfo << "  [control] ||K(1.01x)-K(F=I)||/(1+|K|) = "<<relDiff(K2d,K2)<<"\n";
    CHECK(relDiff(K2d,K2) > 1e-3);
}

// =====================================================================
// ANALYTIC MACHINERY for the three closed-form tests below
// =====================================================================
//
// gsNeoHookeQuadMaterial<T> is gsNeoHookeMaterial<T,QuadraticVolumetric>:
//
//   Psi = mu/2 (J^{-2/3} I1 - 3) + kappa/4 (J^2 - 1 - 2 ln J)
//   kappa = lambda_lame + 2 mu/3 = E/(3(1-2nu)) ,  I1 = tr(C) ,  J = det F
//
// always evaluated in 3D: gsPlaneStressCondensation.h:164 sets data.dim = 3
// regardless of the shell's embedding dimension d, so ONE closed form serves
// both the dim-2 uniaxial test and the dim-3 balloon.
//
// This is EXACTLY the model behind the legacy compressible Neo-Hooke's analytic
// solution (gsThinShellAssembler_test.cpp:695-701): that code uses
// K = 2 mu (1+nu)/(3-6 nu) == kappa and the volumetric derivative
// 0.25*K*(2 J^2/l - 2/l) == d/dl [kappa/4 (J^2-1-2 ln J)]. The identification is
// ASSERTED, not merely asserted-by-comment: TEST(UAT_newton_analytic) re-derives
// the legacy hard-coded root J = 1.105598565 from the Psi above and CHECKs the
// two against each other. That root is an EXTERNAL reference -- it predates this
// project and was computed for a different implementation.

/// dPsi/dlambda_i in principal stretches (i = 0,1,2).
real_t nhDPsi(const index_t i, const real_t l1, const real_t l2, const real_t l3,
              const real_t mu, const real_t kappa)
{
    const real_t J  = l1*l2*l3;
    const real_t I1 = l1*l1 + l2*l2 + l3*l3;
    const real_t li = (0==i ? l1 : (1==i ? l2 : l3));
    const real_t a  = math::pow(J,(real_t)(-2.0/3.0));         // J^{-2/d}, d = 3
    return (real_t)0.5*mu*a*(2*li - 2*I1/(3*li))               // isochoric part
         + (real_t)0.5*kappa*(J - 1/J) * J/li;                 // U'(J) * dJ/dlambda_i
}

/// Bisection root of a scalar residual bracketed on [lo,hi]. 200 halvings take
/// any double bracket to its last bit, and unlike a Newton on the J^{-2/3}
/// nonlinearity it cannot wander out of the physical range lambda > 0.
real_t bisect(const std::function<real_t(real_t)> & f, real_t lo, real_t hi)
{
    real_t flo = f(lo);
    GISMO_ENSURE(flo*f(hi) < 0.0,"bisect: root not bracketed on ["<<lo<<","<<hi<<"]");
    for (index_t k=0; k!=200; ++k)
    {
        const real_t mid = (real_t)0.5*(lo+hi), fm = f(mid);
        if ((fm<0.0) == (flo<0.0)) { lo = mid; flo = fm; }
        else                         hi = mid;
    }
    return (real_t)0.5*(lo+hi);
}

/// Surface area of @a mp, needed because gsThinShellAssembler2 (deliberately)
/// has no getArea: the follower pressure is integrated over the UNDEFORMED
/// measure, so the TRUE pressure it represents is p * A_0/A_def.
real_t surfaceArea(const gsMultiPatch<real_t> & mp, const gsMultiBasis<real_t> & dbasis)
{
    gsExprEvaluator<real_t> ev;
    ev.setIntegrationElements(dbasis);
    auto G = ev.getMap(mp);
    return ev.integral(meas(G));
}

/// CHECKs an assembly status AND reports whether to go on. UnitTest++'s
/// CHECK/CHECK_EQUAL record a failure and CONTINUE, so a status that is not
/// acted on would leave the Newton loop grinding on whatever the aborted
/// routine left behind -- and the real diagnostic (a sequencing regression:
/// status 1 = AssemblyError) would be buried under cascading failures of the
/// analytic gates. Every status below is therefore both CHECKed and obeyed.
bool okStatus(const ThinShellAssemblerStatus st, const char * routine)
{
    const bool ok = (ThinShellAssemblerStatus::Success == st);
    if (!ok)
        gsWarn << "newtonSolve: "<<routine<<" returned status "<<(int)st
               << " (1 = AssemblyError) -- aborting the Newton loop\n";
    CHECK(ok);
    return ok;
}

/**
 * @brief One Newton loop on a gsThinShellAssembler2.
 *
 * SEQUENCING (the task-28 constraint): assemble() runs FIRST and sizes the
 * system; only then may assembleVector be called. Inside the loop
 * assembleMatrix and assembleVector alternate freely -- cleanUp() touches
 * neither m_fmatrix nor m_matrix, so the tangent survives the residual assembly
 * and vice versa.
 *
 * @a homogenize is left at its default true: passing false gives a
 * converged-but-wrong answer that only an analytic oracle would catch.
 *
 * @param[out] its     Newton iterations actually taken; set to @a maxIt on any
 *                     failed assembly, so the caller's "its < N" gate also
 *                     catches a stall and an aborted loop
 * @param[out] resRel  final ||residual|| relative to (1 + the FIRST one)
 * @param[out] res0    the first residual norm itself, so that @a resRel is
 *                     unambiguous (if res0 were small, 1+res0 ~ 1 and resRel
 *                     would silently be an ABSOLUTE norm)
 */
template<short_t d, bool bending>
gsMatrix<real_t> newtonSolve(gsThinShellAssembler2<d,real_t,bending> & A,
                             gsMultiPatch<real_t>                    & mp_def,
                             index_t & its, real_t & resRel, real_t & res0,
                             const index_t maxIt = 40,
                             const real_t  tol   = 1e-13)
{
    its = maxIt; resRel = -1.0; res0 = 0.0;
    gsMatrix<real_t> sol;

    // Linear predictor. This is also what SIZES the system.
    if (!okStatus(A.assemble(),"assemble")) return sol;
    gsSparseSolver<real_t>::LU solver;
    solver.compute( A.matrix() );
    sol = solver.solve( A.rhs() );

    A.constructSolution(sol,mp_def);
    if (!okStatus(A.assembleVector(mp_def),"assembleVector (predictor)")) return sol;
    gsMatrix<real_t> res = A.rhs();          // copy OUT: rhs() is a reference
    res0 = res.norm();

    its = 0;
    for (index_t it=0; it!=maxIt; ++it)
    {
        if (!okStatus(A.assembleMatrix(mp_def),"assembleMatrix")) { its = maxIt; break; }
        solver.compute( A.matrix() );
        gsMatrix<real_t> du = solver.solve(res);
        sol += du;

        A.constructSolution(sol,mp_def);
        if (!okStatus(A.assembleVector(mp_def),"assembleVector")) { its = maxIt; break; }
        res = A.rhs();
        ++its;
        if (du.norm() <= tol*(1.0+sol.norm())) break;
    }
    // A loop that never met the update criterion exits with its == maxIt, which
    // every caller gates on.
    resRel = res.norm()/(1.0+res0);
    return sol;
}

// =====================================================================
// TEST 11 : uniaxial tension, full Newton against the ANALYTIC solution
// =====================================================================
//
// The fixture of gsThinShellAssembler_test.cpp:468-753 (UAT_numerical /
// UAT_analytical), with ONE declared substitution: the legacy suite runs the
// INCOMPRESSIBLE Neo-Hooke, which has no PFF counterpart, so this uses the
// COMPRESSIBLE branch (nu = 0.45) -- the branch whose legacy analytic root
// J = 1.105598565 we can re-derive and cross-check.
//
// Why the oracle is not tautological: the continuous solution is the AFFINE map
// (u,v) -> (lambda u, s v) with s from a scalar plane-stress condition. It is
// representable exactly by the discrete space and satisfies equilibrium
// exactly, so the discrete solution must BE it. The test therefore checks the
// whole displacement field pointwise, not just a derived scalar -- and s is
// computed here by bisection on a closed-form Psi, sharing no code with the
// assembler's own gsPlaneStressCondensation Newton.
//
// This is also a dim-2 NONLINEAR execution of gsThinShellAssembler2<2,.,false>,
// which task 28 could only compile-check.
TEST(UAT_newton_analytic)
{
    const real_t mu        = 1.5e6;
    const real_t nu        = 0.45;      // COMPRESSIBLE (see the note above)
    const real_t E         = 2*mu*(1+nu);
    const real_t kappa     = E/(3*(1-2*nu));
    const real_t thickness = 0.001;
    const real_t LAM       = 2.0;       // imposed stretch in x

    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.degreeElevate(1);
    mp.uniformRefine();
    gsMultiPatch<real_t> mp_def = mp;
    gsMultiBasis<real_t> dbasis(mp);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    gsConstantFunction<real_t> displx(LAM-1.0,2);
    bc.addCondition(boundary::west,  condition_type::dirichlet, 0,        0, false, 0);
    bc.addCondition(boundary::east,  condition_type::dirichlet, &displx,  0, false, 0);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0,        0, false, 1);

    gsVector<real_t> zero(2); zero.setZero();
    gsConstantFunction<real_t> force(zero,2);
    gsFunctionExpr<real_t> t(util::to_string(thickness),2);   // PHYSICAL, dim d = 2
    gsNeoHookeQuadMaterial<real_t> law(E,nu,2);               // parameters: PARAMETRIC dim 2

    gsThinShellAssembler2<2,real_t,false> A(mp,dbasis,bc,force,t,law);

    index_t its; real_t resRel, res0;
    gsMatrix<real_t> sol = newtonSolve(A,mp_def,its,resRel,res0);

    // ---- ANALYTIC: lambda2 = lambda3 = s from dPsi/dlambda2 = 0 -----------
    // Both transverse directions are traction free (the lateral edges carry no
    // load and the plane-stress condensation kills S33), and the material is
    // isotropic, so the two transverse stretches coincide.
    const real_t s  = bisect([&](real_t x){ return nhDPsi(1,LAM,x,x,mu,kappa); },
                             (real_t)1e-3, LAM);
    const real_t J  = LAM*s*s;

    gsInfo << "[UAT_newton_analytic] Newton its = "<<its
           << " , ||res_0|| = "<<res0
           << " , final ||res||/(1+||res_0||) = "<<resRel
           << " , #dofs = "<<A.numDofs()<<"\n";
    gsInfo << "  analytic  lambda2 = lambda3 = "<<s<<"  (J = "<<J<<")\n";
    gsInfo << "  legacy hard-coded J (gsThinShellAssembler_test.cpp:698) = 1.105598565"
           << " ; deviation = "<<math::abs(J-1.105598565)<<"\n";

    CHECK(its  < 20);
    CHECK(resRel < 1e-10);
    // The EXTERNAL cross-check: our closed form reproduces the legacy root.
    CHECK(math::abs(J-1.105598565) < 1e-8);

    // ---- the whole deformed FIELD against the affine analytic map ----------
    // BSplineSquare(1) is the identity map on [0,1]^2, so the analytic deformed
    // position of the parametric point (u,v) is (LAM*u, s*v).
    gsMatrix<real_t> uu(2,6);
    uu.col(0)<<0.0,0.0;  uu.col(1)<<1.0,0.0;  uu.col(2)<<0.0,1.0;
    uu.col(3)<<1.0,1.0;  uu.col(4)<<0.35,0.8; uu.col(5)<<0.7,0.25;
    gsMatrix<real_t> X = mp_def.patch(0).eval(uu);
    gsMatrix<real_t> Xa(2,uu.cols());
    for (index_t k=0; k!=uu.cols(); ++k)
    {
        Xa(0,k) = LAM*uu(0,k);
        Xa(1,k) = s  *uu(1,k);
    }
    const real_t fieldErr = (X-Xa).cwiseAbs().maxCoeff();
    gsInfo << "  max |x_h - x_analytic| over 6 points = "<<fieldErr
           << "   (|x| ~ "<<X.cwiseAbs().maxCoeff()<<")\n";
    gsInfo << "  numerical lambda2 = "<<X(1,3)<<" vs analytic "<<s
           << " , REL = "<<math::abs(X(1,3)-s)/s<<"\n";

    // Non-vacuity: the bar really stretched (a zero solution would trivially
    // satisfy nothing, but a WRONG s of, say, the incompressible 1/sqrt(2)
    // = 0.7071 would be 5e-2 away and is excluded by the gate below).
    CHECK(math::abs(X(0,3)-LAM) < 1e-9);
    CHECK(math::abs(X(1,3)-s)/s < 1e-9);
    CHECK(fieldErr < 1e-9);
    CHECK(sol.norm() > 1e-3);
}

// =====================================================================
// TEST 12 : inflated balloon -- Newton WITH the follower pressure
// =====================================================================
//
// The fixture of gsThinShellAssembler_test.cpp:201-466 (balloon_numerical /
// balloon_analytical): eighth sphere R = 10, thickness 0.1, mu = 4.225e5.
// TWO declared substitutions:
//   (a) the legacy balloon runs INCOMPRESSIBLE materials, which have no PFF
//       counterpart, so this uses the compressible Neo-Hooke twin at nu = 0.45
//       and a matching closed-form inflation curve is derived here;
//   (b) the pressure is 5e3 rather than 10e3, which puts the balloon at
//       lambda ~ 1.15 instead of ~1.46: still a genuine finite deformation
//       (the thickness drops by 15%) but comfortably below the limit point of
//       the TRUE-pressure branch, so the Newton needs no continuation.
//
// The oracle: for an equibiaxially stretched thin sphere,
//   p_true = 2 sigma_1 t_cur / r_cur = 2 sigma_1 (t_0 l3) / (R lambda)
// with sigma_1 = lambda dPsi/dlambda_1 / J the Cauchy stress and l3 the
// plane-stress thickness stretch (dPsi/dlambda_3 = 0). The assembler integrates
// p over the UNDEFORMED measure (p * space * sn(m_def).normalized() * meas(ori)),
// so the input p represents p_true * A_def/A_0 -- exactly the conversion legacy
// makes with getArea(mp)/getArea(mp_def). Both the area identity and the
// uniformity of the inflation are ASSERTED, not assumed.
TEST(Balloon_newton_analytic)
{
    const real_t mu       = 4.225e5;
    const real_t nu       = 0.45;      // COMPRESSIBLE (see the note above)
    const real_t E        = 2*mu*(1+nu);
    const real_t kappa    = E/(3*(1-2*nu));
    const real_t t0       = 0.1;
    const real_t pressure = 5.0e3;

    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>("surfaces/eighth_sphere.xml", mp);
    GISMO_ENSURE(mp.nPatches()==1,"Failed to read 'surfaces/eighth_sphere.xml'");
    mp.patch(0).degreeElevate();
    mp.patch(0).uniformRefine();
    gsMultiPatch<real_t> mp_def = mp;
    gsMultiBasis<real_t> dbasis(mp);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2);
    // symmetry in x
    bc.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 0);
    bc.addCondition(boundary::east,  condition_type::clamped,   0, 0, false, 1);
    bc.addCondition(boundary::east,  condition_type::clamped,   0, 0, false, 2);
    // symmetry in y
    bc.addCondition(boundary::west,  condition_type::clamped,   0, 0, false, 0);
    bc.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, 1);
    bc.addCondition(boundary::west,  condition_type::clamped,   0, 0, false, 2);

    gsVector<real_t> zero(3); zero.setZero();
    gsConstantFunction<real_t> force(zero,3);
    gsConstantFunction<real_t> pressFun(pressure,3);
    gsFunctionExpr<real_t> t(util::to_string(t0),3);
    gsNeoHookeQuadMaterial<real_t> law(E,nu,2);

    // ---- SANITY LEG: this fixture (its clamped symmetry conditions in
    // particular) must reproduce the legacy assembler EXACTLY before any
    // analytic claim is made about it. If the analytic leg below ever fails,
    // this line says whether the cause is the fixture or the derivation.
    {
        gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,law);
        gsMaterialMatrix3D<3,real_t> adapter(mp,t,law);
        gsThinShellAssembler<3,real_t,true>  AL(mp,dbasis,bc,force,&adapter);
        CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
        gsSparseMatrix<real_t> K2 = A2.matrix();
        CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
        gsSparseMatrix<real_t> KL = AL.matrix();
        gsInfo << "[Balloon_newton_analytic] fixture parity (clamped symmetry BCs):\n";
        const real_t rK = reportParity("linear K",K2,KL);
        CHECK(KL.norm() > 1.0);
        CHECK_EQUAL(A2.numDofs(),AL.numDofs());
        CHECK(rK <= TOL_EXACT);
    }

    gsThinShellAssembler2<3,real_t,true> A(mp,dbasis,bc,force,t,law);
    A.setPressure(pressFun);

    index_t its; real_t resRel, res0;
    gsMatrix<real_t> sol = newtonSolve(A,mp_def,its,resRel,res0);

    // ---- measured stretch: UNIFORM inflation, verified point by point ------
    gsMatrix<real_t> uu(2,5);
    uu.col(0)<<0.2,0.2;  uu.col(1)<<0.5,0.5;  uu.col(2)<<0.8,0.4;
    uu.col(3)<<0.4,0.9;  uu.col(4)<<1.0,1.0;
    gsMatrix<real_t> X0 = mp    .patch(0).eval(uu);
    gsMatrix<real_t> X1 = mp_def.patch(0).eval(uu);
    real_t lamMin = 1e30, lamMax = -1e30, lamSum = 0.0, Rsum = 0.0;
    real_t Rmin = 1e30, Rmax = -1e30;
    for (index_t k=0; k!=uu.cols(); ++k)
    {
        const real_t Rk = X0.col(k).norm();
        const real_t lk = X1.col(k).norm()/Rk;
        lamMin = math::min(lamMin,lk);  lamMax = math::max(lamMax,lk);
        Rmin   = math::min(Rmin,Rk);    Rmax   = math::max(Rmax,Rk);
        lamSum += lk;  Rsum += Rk;
    }
    const real_t lamNum = lamSum/uu.cols();
    const real_t R      = Rsum  /uu.cols();
    const real_t A0     = surfaceArea(mp,dbasis);
    const real_t Adef   = surfaceArea(mp_def,dbasis);

    gsInfo << "  Newton its = "<<its<<" , ||res_0|| = "<<res0
           << " , final ||res||/(1+||res_0||) = "<<resRel
           << " , #dofs = "<<A.numDofs()<<"\n";
    gsInfo << "  R = "<<R<<" (fixture spread "<<(Rmax-Rmin)/R
           << ") , lambda = "<<lamNum
           << " (spread over 5 points = "<<(lamMax-lamMin)/lamNum<<")\n";
    gsInfo << "  A_def/A_0 = "<<Adef/A0<<" vs lambda^2 = "<<lamNum*lamNum
           << " , REL = "<<math::abs(Adef/A0-lamNum*lamNum)/(lamNum*lamNum)<<"\n";
    CHECK(its < 25);
    CHECK(resRel < 1e-10);
    CHECK(lamNum > 1.05);                                   // it really inflated
    // The inflation must be UNIFORM -- that is what makes "one lambda" and the
    // area identity A_def/A_0 = lambda^2 meaningful. It is uniform to O(1e-4),
    // not to machine precision, and deliberately not gated tighter: the fixture
    // is only an APPROXIMATE sphere (eighth_sphere.xml perturbs the degenerate
    // pole by 0.01 to keep the patch regular, and the exactly-spherical octant
    // is anyway distorted by the h-refinement of a rational patch), so a
    // spread of this size is fixture geometry, not solver error. Both measured
    // values are printed above; 1e-3 keeps ~8x margin over them while still
    // excluding any qualitatively non-uniform deformation.
    CHECK((lamMax-lamMin)/lamNum < 1e-3);
    CHECK(math::abs(Adef/A0-lamNum*lamNum)/(lamNum*lamNum) < 1e-3);

    // ---- analytic inflation curve -----------------------------------------
    // pApplied(l) is the pressure the assembler must be GIVEN to reach stretch
    // l, i.e. the true membrane pressure times A_def/A_0 = l^2.
    std::function<real_t(real_t)> pApplied = [&](real_t l)
    {
        const real_t l3 = bisect([&](real_t x){ return nhDPsi(2,l,l,x,mu,kappa); },
                                 (real_t)1e-3, l);
        const real_t Jl = l*l*l3;
        const real_t sig1 = l*nhDPsi(0,l,l,l3,mu,kappa)/Jl;   // Cauchy stress
        return 2*sig1*(t0*l3)/(R*l) * l*l;
    };

    const real_t lamAna = bisect([&](real_t l){ return pApplied(l)-pressure; },
                                 (real_t)1.0001, (real_t)1.6);
    const real_t l3Ana  = bisect([&](real_t x){ return nhDPsi(2,lamAna,lamAna,x,mu,kappa); },
                                 (real_t)1e-3, lamAna);
    const real_t Ptrue    = pressure*A0/Adef;                 // measured, legacy's form
    const real_t PtrueAna = pApplied(lamNum)/(lamNum*lamNum);  // analytic at the MEASURED lambda

    gsInfo << "  analytic lambda = "<<lamAna<<" (l3 = "<<l3Ana<<") vs numerical "<<lamNum
           << " , REL = "<<math::abs(lamNum-lamAna)/lamAna<<"\n";
    gsInfo << "  true pressure: measured p*A_0/A_def = "<<Ptrue
           << " vs analytic "<<PtrueAna
           << " , REL = "<<math::abs(Ptrue-PtrueAna)/PtrueAna<<"\n";

    // The residual model error is the KIRCHHOFF-LOVE shell against a pure
    // MEMBRANE oracle: O((t/R)^2) = 1e-4 here. The stretch form is the primary
    // gate because d(ln p)/d(ln lambda) ~ 6 on this branch, so the pressure form
    // AMPLIFIES the same discrepancy sixfold.
    CHECK(math::abs(lamNum-lamAna)/lamAna < 1e-3);
    CHECK(math::abs(Ptrue-PtrueAna)/PtrueAna < 1e-2);

    // ---- CONTROL: without the pressure nothing moves -----------------------
    // Guards against an oracle that would be satisfied by any inflated sphere.
    gsMultiPatch<real_t> mp_def0 = mp;
    gsThinShellAssembler2<3,real_t,true> A0a(mp,dbasis,bc,force,t,law);
    index_t its0; real_t resRel0, res00;
    gsMatrix<real_t> sol0 = newtonSolve(A0a,mp_def0,its0,resRel0,res00);
    gsInfo << "  [control] no pressure: |sol| = "<<sol0.norm()
           << " vs |sol(p)| = "<<sol.norm()<<"\n";
    CHECK(sol0.norm() < 1e-8*(1.0+sol.norm()));
}

// =====================================================================
// TEST 13 : modal analysis -- assembleMass against ANALYTIC frequencies
// =====================================================================
//
// The fixture of gsThinShellAssembler_test.cpp:755-924 (Modal_numerical /
// Modal_analytical, isotropic branch): a simply-supported unit square plate,
// t = 0.01, E = 1e5, rho = 1, nu = 0.3, and
//     omega_mn = pi^2 (m^2 + n^2) sqrt(D/(rho t)) ,  D = E t^3/(12(1-nu^2)) .
//
// The plate is FLAT ON PURPOSE. The Kirchhoff-Love plate constant D above is a
// z-CONSTANT-tangent moment, which is the very thing TEST(Linear3D_vs_legacy_model)
// shows the new path does NOT use -- it Gauss-integrates the exact z-dependent
// metric, and the two differ by O(t^2/R^2). At R = infinity that difference is
// identically zero, which is what makes the classical oracle legitimate here.
// Do not "improve" this fixture by curving it.
//
// Constraint 2 (task 28) USED TO BE: assembleMass() leaves the trial space
// HOMOGENIZED and no routine re-runs _assembleDirichlet. Task 44 fixed that --
// assembleMass now restores the l2Projection setup itself -- but the mass-first
// ordering is still asserted below on the STIFFNESS that feeds the eigenproblem:
// it must reproduce the safe-order K bit for bit. The updateBCs(bc) call is now
// redundant and is kept only so that the OLD remedy stays exercised. (The rhs
// half is gated by TEST(NeumannAndPointLoads) above, which carries inhomogeneous
// Dirichlet data, and without the remedy by TEST(MassRestoresDirichlet).)
TEST(Modal_vs_analytic)
{
    const real_t thickness = 0.01;
    const real_t E         = 1e5;
    const real_t density   = 1.0;
    const real_t nu        = 0.3;
    const index_t nElev    = 2, nRef = 4;      // degree 3, 16x16 elements

    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.addAutoBoundaries();
    mp.embed(3);
    for (index_t i=0; i!=nElev; ++i) mp.patch(0).degreeElevate();
    for (index_t i=0; i!=nRef;  ++i) mp.patch(0).uniformRefine();
    gsMultiBasis<real_t> dbasis(mp);

    gsBoundaryConditions<real_t> bc;
    for (index_t c=0; c!=3; ++c)
    {
        bc.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, c);
        bc.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, c);
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, c);
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, c);
    }
    bc.setGeoMap(mp);

    gsVector<real_t> zero(3); zero.setZero();
    gsConstantFunction<real_t> force(zero,3);
    gsFunctionExpr<real_t> t(util::to_string(thickness),3);
    gsConstantFunction<real_t> rho(density,3);
    gsLinearMaterial<real_t> law(E,nu,2);

    gsThinShellAssembler2<3,real_t,true> A(mp,dbasis,bc,force,t,law,rho);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A.assemble());
    gsSparseMatrix<real_t> K = A.matrix();       // copy OUT before the mass
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A.assembleMass());
    gsSparseMatrix<real_t> M = A.massMatrix();

    gsEigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base > eig;
    eig.compute(K,M);
    gsMatrix<real_t> omega = eig.eigenvalues().cwiseSqrt();

    // analytic, sorted
    const real_t D = E*math::pow(thickness,3)/(12*(1-nu*nu));
    const real_t pi = 4*math::atan((real_t)1.0);
    std::vector<real_t> ana;
    for (index_t m=1; m!=8; ++m)
        for (index_t n=1; n!=8; ++n)
            ana.push_back((m*m+n*n)*pi*pi*math::sqrt(D/(density*thickness)));
    std::sort(ana.begin(),ana.end());

    const index_t nModes = 10;
    real_t worst = 0.0;
    gsInfo << "[Modal_vs_analytic] #dofs = "<<A.numDofs()
           << " , #elem = "<<dbasis.totalElements()<<"\n";
    for (index_t k=0; k!=nModes; ++k)
    {
        const real_t rel = math::abs(omega(k,0)-ana[k])/ana[k];
        worst = math::max(worst,rel);
        gsInfo << "   mode "<<k<<" : numerical "<<omega(k,0)
               << " , analytic "<<ana[k]<<" , REL = "<<rel<<"\n";
    }
    gsInfo << "  worst relative frequency error over "<<nModes<<" modes = "<<worst<<"\n";
    CHECK(M.norm() > 1e-8);
    CHECK(ana[0] > 0.0);
    CHECK(worst < 1e-3);

    // ---- CONTROL: the MASS is what sets the scale --------------------------
    // sqrt of a generalized eigenvalue scales as 1/sqrt(rho): assembling with a
    // 4x density must divide every frequency by exactly 2. A mass matrix that
    // ignored the density (or the thickness) would fail this.
    gsConstantFunction<real_t> rho4(4*density,3);
    gsThinShellAssembler2<3,real_t,true> A4(mp,dbasis,bc,force,t,law,rho4);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A4.assembleMass());
    gsSparseMatrix<real_t> M4 = A4.massMatrix();
    gsEigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base > eig4;
    eig4.compute(K,M4);
    gsMatrix<real_t> omega4 = eig4.eigenvalues().cwiseSqrt();
    const real_t scale = omega(0,0)/omega4(0,0);
    gsInfo << "  [control] omega(rho)/omega(4 rho) = "<<scale<<"   (must be 2)\n";
    CHECK(math::abs(scale-2.0) < 1e-8);

    // ---- constraint 2: assembleMass first, then the REMEDY ------------------
    gsThinShellAssembler2<3,real_t,true> AB(mp,dbasis,bc,force,t,law,rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AB.assembleMass());
    AB.updateBCs(bc);                            // the verified remedy
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AB.assemble());
    const real_t rK = relDiff(AB.matrix(),K);
    gsInfo << "  [mass-first + updateBCs] ||K_b - K||/(1+||K||) = "<<rK
           << ((AB.matrix()-K).norm()==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    CHECK(rK <= TOL_EXACT);
}

// =====================================================================
// TEST 16 : the LINEAR assemble() must apply the follower pressure
// =====================================================================
//
// REGRESSION GATE for the task-29 finding: gsThinShellAssembler2::assemble()
// silently assembled a ZERO pressure load, with no error and no warning.
// Legacy's linear routine applies the pressure through the
// UNDEFORMED-geometry _assemblePressure overload (called for BOTH _matrix
// values at gsThinShellAssembler.hpp:1793-1797 and :1871-1875):
//     load   : p * space * usn(m_def) * meas(m_ori)      (.hpp:389-407)
//     matrix : EMPTY -- "No matrix contribution for the linear case" (.hpp:381)
// In the linear routine defpatches == m_patches, so m_def and m_ori are the
// SAME geometry: the load is the pressure on the UNDEFORMED normal.
//
// Three legs, and they are three different claims:
//  (1) PARITY of the loaded rhs against LEGACY + gsMaterialMatrix3D over the
//      same law -- the usual ring-1 gate;
//  (2) the REGRESSION control: switching the pressure on must MOVE the linear
//      rhs. This is the leg that FAILS against the pre-fix code, where
//      r(p) - r(0) is exactly the zero vector;
//  (3) the linear stiffness must be UNTOUCHED by the pressure -- the only
//      assertion anywhere on the empty _matrix body, and what keeps someone
//      from "fixing" this by wiring in the DEFORMED follower stiffness.
//
// Single-threaded: legs (1) and (3) are exactness claims and gsExprAssembler
// accumulates under an unordered #pragma omp atomic.
TEST(LinearFollowerPressure)
{
    SingleThreadScope pin;
    GISMO_UNUSED(pin);

    gsMultiPatch<real_t> mp = makeRoof();
    gsMultiBasis<real_t> dbasis(mp);
    gsBoundaryConditions<real_t> bc;
    roofBCs(mp,bc);

    gsVector<real_t> tmp(3); tmp << 0, 0, -90;
    gsConstantFunction<real_t> force(tmp,3);
    gsFunctionExpr<real_t> t(util::to_string(THICK),3);
    gsConstantFunction<real_t> pressFun(1.0e4,3);      // as in TEST 10
    gsLinearMaterial<real_t> svk(E_MOD,NU,2);

    gsMaterialMatrix3D<3,real_t> adapter(mp,t,svk);

    gsInfo << "[LinearFollowerPressure] setPressure(1e4) + assemble(), #elem = "
           << dbasis.totalElements()<<"\n";

    // ---- WITH pressure: new vs legacy -------------------------------------
    gsThinShellAssembler2<3,real_t,true> A2(mp,dbasis,bc,force,t,svk);
    gsThinShellAssembler<3,real_t,true>  AL(mp,dbasis,bc,force,&adapter);
    A2.setPressure(pressFun);
    AL.setPressure(pressFun);

    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsSparseMatrix<real_t> K2 = A2.matrix();          // copy OUT
    gsMatrix<real_t>       r2 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsSparseMatrix<real_t> KL = AL.matrix();
    gsMatrix<real_t>       rL = AL.rhs();

    const real_t rK = reportParity("linear K   (p)",K2,KL);
    const real_t rR = reportParity("linear rhs (p)",r2,rL);
    CHECK(KL.norm() > 1.0);
    CHECK(rL.norm() > 1.0);
    CHECK_EQUAL(A2.numDofs(),AL.numDofs());
    CHECK(rK <= TOL_EXACT);
    CHECK(rR <= TOL_EXACT);                                        // leg (1)

    // ---- WITHOUT pressure: FRESH instances --------------------------------
    // m_pressInd has no way back to false, and m_rhs is never cleared
    // (hazard 1), so the unpressurised reference needs its own assemblers.
    gsThinShellAssembler2<3,real_t,true> A20(mp,dbasis,bc,force,t,svk);
    gsThinShellAssembler<3,real_t,true>  AL0(mp,dbasis,bc,force,&adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A20.assemble());
    gsSparseMatrix<real_t> K20 = A20.matrix();
    gsMatrix<real_t>       r20 = A20.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL0.assemble());
    gsMatrix<real_t>       rL0 = AL0.rhs();

    const real_t d2 = relDiff(r2,r20), dL = relDiff(rL,rL0);
    gsInfo << "  [regression] new    : ||r(p)-r(0)||/(1+|r(0)|) = "<<d2
           << "   (|r(0)| = "<<r20.norm()<<" , |r(p)| = "<<r2.norm()<<")\n";
    gsInfo << "  [regression] legacy : ||r(p)-r(0)||/(1+|r(0)|) = "<<dL<<"\n";
    // leg (2). p = 1e4 against a surface force of -90 dominates the rhs by two
    // orders of magnitude, so the floor is set an order of magnitude below the
    // measured margin (110.7 on this fixture: |r| goes from 17274.8 to
    // 1.93e6) and remains vastly above any conceivable round-off. Pre-fix, the
    // new side measures EXACTLY 0 here while legacy still measures 110.7.
    CHECK(d2 > 10.0);
    CHECK(dL > 10.0);
    // ... and the two UNPRESSURISED rhs must still agree to machine precision,
    // so the margin above is a pressure effect and not a drifting baseline.
    CHECK(relDiff(r20,rL0) <= TOL_EXACT);

    // ---- leg (3): the linear stiffness carries NO pressure contribution ----
    // Bit-identity across two instances is a legitimate claim here, and only
    // here: the team size is pinned to 1 (so gsExprAssembler's unordered atomic
    // accumulation is out of play), and the extra registrations the pressure
    // body makes in A2 (one coefficient, both maps, the trial space) never
    // enter the PARSED stiffness expression, so they cannot move a quadrature
    // point or an accumulation order. If this ever reads 1e-16 instead of 0,
    // that premise broke and is worth investigating rather than relaxing.
    const real_t dK = absDiff(K2,K20);
    gsInfo << "  [matrix no-op] ||K(p)-K(0)|| = "<<dK
           << (dK==0.0 ? "   (EXACTLY 0)" : "")<<"\n";
    CHECK(dK == 0.0);
}

// =====================================================================
// TASK 44 : the four assembleMass() defects
// =====================================================================
//
// All four lived in ONE routine (gsThinShellAssembler::assembleMass, and its
// gsThinShellAssembler2 twin for two of them) and none of them had a single
// assertion anywhere: a green suite proved nothing about any of them.
//
//   D1  m_status was returned UNINITIALISED on the success path
//   D2  lumped = true did not lump
//   D3  the trial space was left HOMOGENIZED
//   D4  point masses were silently dropped (legacy only -- gsThinShellAssembler2
//       has no point-mass API at all)
//
// Five of the seven tests below are the regression gates for those four defects,
// each written so that it FAILS against the pre-fix code for a DETERMINISTIC
// reason -- see the individual comments, and in particular the note on D1
// (Success == 0, so a virgin-instance status check is a coin flip, not a test).
// The sixth, TEST(MassFailurePath_legacy), gates the FAILURE path of the fixed
// routine, which has no pre-fix counterpart to poison.
//
// TASK 62 added the seventh, TEST(MassPointMassOutOfDomain_legacy), plus the
// per-entry "WHICH diagonal" assertion inside the two lumping tests and
// TEST(MassPointMasses_legacy): restoring the _applyMass() call made a previously
// dead body reachable, and an out-of-domain point mass then NaN-ed the whole mass
// matrix under a Success status.
//
// TASK 63 added the eighth, TEST(MassPointMassSpaceBasis_legacy), and two more
// halves to the seventh: the same dead body also numbered its actives in the wrong
// basis whenever setSpaceBasis() had been used, and task 62's own tolerance band
// admitted points whose mass was then silently dropped.

/// The Scordelis-Lo self-weight vector of the fixtures above, as a value.
gsVector<real_t> roofForce()
{
    gsVector<real_t> v(3);
    v << 0, 0, -90;
    return v;
}

/// Shared fixture pieces for the task-44 tests. The members are declared in
/// dependency order: the adapter and every assembler below hold REFERENCES into
/// this object, so it must outlive them (it does -- it is the first local of each
/// test).
struct MassFixture
{
    gsMultiPatch<real_t>          mp;
    gsMultiBasis<real_t>          dbasis;
    gsConstantFunction<real_t>    force;
    gsFunctionExpr<real_t>        t, rho;
    gsLinearMaterial<real_t>      law;
    gsMaterialMatrix3D<3,real_t>  adapter;

    MassFixture()
    :
    mp(makeRoof()), dbasis(mp), force(roofForce(),3),
    t(util::to_string(THICK),3), rho("1.0",3),
    law(E_MOD,NU,2), adapter(mp,t,rho,law)
    { }
};

/// ONE parametric point mass of magnitude @a value at (@a u, @a v) on patch 0.
/// _applyMass asserts the value is a scalar, hence the size-1 vector.
gsPointLoads<real_t> onePointMass(real_t u, real_t v, real_t value)
{
    gsPointLoads<real_t> pl;
    gsVector<real_t> pt(2);  pt << u, v;
    gsVector<real_t> pv(1);  pv << value;
    pl.addLoad(pt,pv,0,true);
    return pl;
}

/// Total (grand-sum) mass of a FRESH roof assembler carrying @a pm (nullptr for
/// none). The status is handed back through @a status and the grand sum is taken
/// of whatever m_mass holds afterwards, so a rejected point mass is observable
/// both ways -- as the status AND as the matrix left behind.
real_t roofMassTotal(MassFixture & f, const gsBoundaryConditions<real_t> & bc,
                     const gsPointLoads<real_t> * pm, int & status)
{
    gsThinShellAssembler<3,real_t,true> A(f.mp,f.dbasis,bc,f.force,&f.adapter);
    if (pm!=nullptr) A.setPointMass(*pm);
    status = (int)A.assembleMass();
    return grandSum(A.massMatrix());
}

// ---------------------------------------------------------------------
// D1 : the status of a successful assembleMass() must be Success
// ---------------------------------------------------------------------
//
// NON-VACUITY, and why this test is shaped the way it is: Success == 0, so
// reading an uninitialised m_status on a virgin instance has a real chance of
// yielding 0 and PASSING against the broken code. The test therefore first drives
// m_status to a DETERMINISTIC non-Success value and then demands that the mass
// assembly overwrite it. The lever is hazard 4 of the class doc: assembleVector
// cannot be the first assembly (initVector(1) does not size m_fmatrix, the
// gsExprAssembler assert fires and is swallowed by the routine's own catch), so
// a first-call assembleVector returns AssemblyError every time. That single
// sequence covers BOTH halves of the criterion -- the failure path and the
// success path.
TEST(MassStatus_legacy)
{
    MassFixture f;
    gsBoundaryConditions<real_t> bc;
    roofBCs(f.mp,bc);

    gsThinShellAssembler<3,real_t,true> AL(f.mp,f.dbasis,bc,f.force,&f.adapter);

    const int st0 = (int)AL.assembleVector(f.mp);
    gsInfo << "[MassStatus_legacy] first-call assembleVector -> status "<<st0
           << "   (deterministic AssemblyError = "
           << (int)ThinShellAssemblerStatus::AssemblyError<<")\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,st0);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,(int)AL.status());

    // ... and now the SUCCESS path must overwrite it. Pre-fix this returned the
    // AssemblyError left above, because assembleMass assigned m_status nowhere.
    const int st1 = (int)AL.assembleMass();
    gsInfo << "  consistent assembleMass() -> status "<<st1
           << "   (must be Success = "<<(int)ThinShellAssemblerStatus::Success<<")\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,st1);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.status());
    CHECK(AL.massMatrix().nonZeros() > 0);

    // the lumped path reports through the same member
    const int st2 = (int)AL.assembleMass(true);
    gsInfo << "  lumped assembleMass(true) -> status "<<st2<<"\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,st2);
}

// ---------------------------------------------------------------------
// D2 : lumped = true must return a DIAGONAL matrix conserving the total mass
// ---------------------------------------------------------------------
//
// The pre-fix symptom was PROBE-DEPENDENT, so a single call sequence is not a
// baseline. The lumped expression was vector-valued and landed in the expression
// assembler's RHS; m_mass was then read from m_assembler.matrix(), a cache that
// no reset path clears. Three outcomes, by call history:
//
//   (a) assembleMass(true) FIRST on a fresh instance      -> m_mass 0x0 EMPTY
//   (b) assemble(); assembleMass(true)                    -> sized, nonZeros == 0
//   (c) assemble(); (void)matrix(); assembleMass(true)    -> the STALE stiffness
//       (matrix() sets m_modified = false, so the next matrix() hands back the
//       same cached object -- this is what the ledger measured as
//       "bit-identical to the consistent matrix")
//
// All three are gated. (a) and (b) fail pre-fix on the shape/nnz assertions,
// (c) fails on the off-diagonal assertion AND on the "lumped != stiffness" one.
//
// The fixture has REAL Dirichlet elimination (roofBCs), which is the case where
// the cheaper "read m_assembler.rhs()" lumping would NOT conserve the total mass:
// push<false> has no column loop, so it also deposits the eliminated columns'
// contributions. Row-summing the assembled matrix conserves it by construction.
TEST(MassLumping_legacy)
{
    MassFixture f;
    gsBoundaryConditions<real_t> bc;
    roofBCs(f.mp,bc);

    // ---- (a) lumped mass as the FIRST assembly on a FRESH instance ----------
    gsThinShellAssembler<3,real_t,true> Aa(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Aa.assembleMass(true));
    gsSparseMatrix<real_t> Ma = Aa.massMatrix();      // copy OUT
    gsInfo << "[MassLumping_legacy] (a) fresh, lumped-first : "<<Ma.rows()<<"x"<<Ma.cols()
           << " , nnz = "<<Ma.nonZeros()<<" , #dofs = "<<Aa.numDofs()
           << "   (pre-fix: 0x0)\n";
    CHECK_EQUAL(Aa.numDofs(),Ma.rows());
    CHECK_EQUAL(Aa.numDofs(),Ma.cols());
    CHECK(Ma.nonZeros() > 0);
    CHECK(maxOffDiag(Ma) == 0.0);

    // The fixture must actually eliminate dofs, or the conservation claim below
    // is the trivial one.
    CHECK(Aa.numDofs() < 3*(index_t)f.dbasis.totalSize());

    // ---- (b) assemble() first, but WITHOUT touching matrix() ---------------
    gsThinShellAssembler<3,real_t,true> Ab(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ab.assemble());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ab.assembleMass(true));
    gsSparseMatrix<real_t> Mb = Ab.massMatrix();
    gsInfo << "  (b) assemble(); lumped : nnz = "<<Mb.nonZeros()
           << "   (pre-fix: 0)\n";
    CHECK(Mb.nonZeros() > 0);
    CHECK(maxOffDiag(Mb) == 0.0);

    // ---- (c) assemble(); matrix(); lumped  -- the STALE-CACHE sequence -----
    gsThinShellAssembler<3,real_t,true> Ac(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ac.assemble());
    gsSparseMatrix<real_t> K = Ac.matrix();           // flips m_modified to false
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ac.assembleMass(true));
    gsSparseMatrix<real_t> Mc = Ac.massMatrix();
    gsInfo << "  (c) assemble(); matrix(); lumped : nnz = "<<Mc.nonZeros()
           << " , max|offdiag| = "<<maxOffDiag(Mc)
           << " , ||M_l - K|| = "<<absDiff(Mc,K)<<"   (pre-fix: M_l WAS K)\n";
    CHECK(Mc.nonZeros() > 0);
    CHECK(maxOffDiag(Mc) == 0.0);
    CHECK(absDiff(Mc,K) > 1.0);                       // not the stale stiffness

    // the three sequences must all deliver the SAME lumped matrix
    CHECK(relDiff(Ma,Mb) <= TOL_EXACT);
    CHECK(relDiff(Ma,Mc) <= TOL_EXACT);

    // ---- mass conservation, on a fixture WITH eliminated dofs ---------------
    gsThinShellAssembler<3,real_t,true> Ad(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ad.assembleMass());
    gsSparseMatrix<real_t> Mcons = Ad.massMatrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ad.assembleMass(true));
    gsSparseMatrix<real_t> Mlump = Ad.massMatrix();

    const real_t tc = grandSum(Mcons), tl = grandSum(Mlump);
    gsInfo << "  total mass: consistent "<<tc<<" vs lumped "<<tl
           << " , REL = "<<math::abs(tc-tl)/math::abs(tc)<<"\n";
    CHECK(math::abs(tc) > 1e-3);
    CHECK(math::abs(tc-tl) <= 1e-12*math::abs(tc));

    // ... and the lumped matrix is genuinely DIFFERENT from the consistent one,
    // so the diagonality gate above is not passing on a degenerate input.
    CHECK(maxOffDiag(Mcons) > 1e-6);
    CHECK(relDiff(Mlump,Mcons) > 1e-2);

    // ---- WHICH diagonal (task 62) ------------------------------------------
    // Everything above pins that the result is diagonal, conserves the grand sum
    // and differs from the consistent matrix -- and NONE of it pins which
    // diagonal: an "all the mass on dof 0" bug passes all three. This is the
    // per-entry definition of row-sum lumping, and it is the gate that such a bug
    // fails. Measured deviation is exactly 0 at OMP=1; TOL_EXACT absorbs the few
    // ulp that the two assembleMass() calls may differ by under the assembler's
    // unordered atomic accumulation.
    const gsMatrix<real_t> rs = consistentRowSums(Mcons);
    real_t maxDev = 0.0;
    for (index_t i = 0; i != Mcons.rows(); ++i)
    {
        maxDev = math::max(maxDev, math::abs(Mlump.coeff(i,i)-rs(i,0))/(1.0+math::abs(rs(i,0))));
        CHECK(math::abs(Mlump.coeff(i,i) - rs(i,0)) <= TOL_EXACT*(1.0 + math::abs(rs(i,0))));
    }
    gsInfo << "  per-entry: max_i |M_l(i,i) - rowsum_i(M_c)| (rel) = "<<maxDev
           << " , the row sums spread over ["<<rs.minCoeff()<<", "<<rs.maxCoeff()
           << "]\n";
}

/// The same four claims on gsThinShellAssembler2, which carried the identical
/// defect (gsThinShellAssembler2.hpp:533-539, Space chain 0+0+1+0 = 1).
TEST(MassLumping_new)
{
    MassFixture f;
    gsBoundaryConditions<real_t> bc;
    roofBCs(f.mp,bc);

    gsThinShellAssembler2<3,real_t,true> Aa(f.mp,f.dbasis,bc,f.force,f.t,f.law,f.rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Aa.assembleMass(true));
    gsSparseMatrix<real_t> Ma = Aa.massMatrix();
    gsInfo << "[MassLumping_new] (a) fresh, lumped-first : "<<Ma.rows()<<"x"<<Ma.cols()
           << " , nnz = "<<Ma.nonZeros()<<"   (pre-fix: 0x0)\n";
    CHECK_EQUAL(Aa.numDofs(),Ma.rows());
    CHECK(Ma.nonZeros() > 0);
    CHECK(maxOffDiag(Ma) == 0.0);
    CHECK(Aa.numDofs() < 3*(index_t)f.dbasis.totalSize());

    gsThinShellAssembler2<3,real_t,true> Ab(f.mp,f.dbasis,bc,f.force,f.t,f.law,f.rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ab.assemble());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ab.assembleMass(true));
    gsSparseMatrix<real_t> Mb = Ab.massMatrix();
    gsInfo << "  (b) assemble(); lumped : nnz = "<<Mb.nonZeros()<<"   (pre-fix: 0)\n";
    CHECK(Mb.nonZeros() > 0);
    CHECK(maxOffDiag(Mb) == 0.0);

    gsThinShellAssembler2<3,real_t,true> Ac(f.mp,f.dbasis,bc,f.force,f.t,f.law,f.rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ac.assemble());
    gsSparseMatrix<real_t> K = Ac.matrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ac.assembleMass(true));
    gsSparseMatrix<real_t> Mc = Ac.massMatrix();
    gsInfo << "  (c) assemble(); matrix(); lumped : max|offdiag| = "<<maxOffDiag(Mc)
           << " , ||M_l - K|| = "<<absDiff(Mc,K)<<"   (pre-fix: M_l WAS K)\n";
    CHECK(maxOffDiag(Mc) == 0.0);
    CHECK(absDiff(Mc,K) > 1.0);

    CHECK(relDiff(Ma,Mb) <= TOL_EXACT);
    CHECK(relDiff(Ma,Mc) <= TOL_EXACT);

    gsThinShellAssembler2<3,real_t,true> Ad(f.mp,f.dbasis,bc,f.force,f.t,f.law,f.rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ad.assembleMass());
    gsSparseMatrix<real_t> Mcons = Ad.massMatrix();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Ad.assembleMass(true));
    gsSparseMatrix<real_t> Mlump = Ad.massMatrix();

    const real_t tc = grandSum(Mcons), tl = grandSum(Mlump);
    gsInfo << "  total mass: consistent "<<tc<<" vs lumped "<<tl
           << " , REL = "<<math::abs(tc-tl)/math::abs(tc)<<"\n";
    CHECK(math::abs(tc) > 1e-3);
    CHECK(math::abs(tc-tl) <= 1e-12*math::abs(tc));
    CHECK(maxOffDiag(Mcons) > 1e-6);
    CHECK(relDiff(Mlump,Mcons) > 1e-2);

    // WHICH diagonal (task 62) -- see TEST(MassLumping_legacy) for the reasoning.
    const gsMatrix<real_t> rs = consistentRowSums(Mcons);
    real_t maxDev = 0.0;
    for (index_t i = 0; i != Mcons.rows(); ++i)
    {
        maxDev = math::max(maxDev, math::abs(Mlump.coeff(i,i)-rs(i,0))/(1.0+math::abs(rs(i,0))));
        CHECK(math::abs(Mlump.coeff(i,i) - rs(i,0)) <= TOL_EXACT*(1.0 + math::abs(rs(i,0))));
    }
    gsInfo << "  per-entry: max_i |M_l(i,i) - rowsum_i(M_c)| (rel) = "<<maxDev
           << " , the row sums spread over ["<<rs.minCoeff()<<", "<<rs.maxCoeff()
           << "]\n";

    // ring 1: the two assemblers must agree on the LUMPED matrix too, not only
    // on the consistent one (TEST(Mass_vs_legacy_varying_t_rho) covers that).
    gsThinShellAssembler<3,real_t,true> AL(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleMass(true));
    gsSparseMatrix<real_t> MlumpL = AL.massMatrix();
    gsInfo << "  [ring 1] ||M_l(new) - M_l(legacy)||/(1+|M_l|) = "
           << relDiff(Mlump,MlumpL)<<"\n";
    CHECK(relDiff(Mlump,MlumpL) <= TOL_EXACT);
}

// ---------------------------------------------------------------------
// D3 : assembleMass() must not leave the trial space HOMOGENIZED
// ---------------------------------------------------------------------
//
// assembleMass sets the space up with dirichlet::homogeneous (correctly -- a mass
// operator has no Dirichlet lifting) but used to leave it that way, and neither
// assemble() nor assembleMatrix()/assembleVector(false) re-runs _assembleDirichlet.
// With NON-ZERO Dirichlet data the sequence assemble() -> assembleMass() ->
// assemble() therefore returned the homogenized rhs, four orders of magnitude off.
//
// NON-VACUITY: the homogeneous CONTROL (identical structure, zero Dirichlet data)
// is asserted to be O(1) away from the real rhs, so the "same rhs" assertion
// cannot pass by both sides being homogenized. Pre-fix, r1 equals the CONTROL,
// not r0.
TEST(MassRestoresDirichlet)
{
    MassFixture f;

    gsFunctionExpr<real_t> gD("0.02",3);          // NON-ZERO Dirichlet datum
    gsBoundaryConditions<real_t> bc;   roofBCs(f.mp,bc,&gD);
    gsBoundaryConditions<real_t> bcH;  roofBCs(f.mp,bcH,nullptr);

    // ---- the homogeneous control (a separate, FRESH instance) --------------
    gsThinShellAssembler<3,real_t,true> ALh(f.mp,f.dbasis,bcH,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALh.assemble());
    gsMatrix<real_t> rHom = ALh.rhs();            // copy OUT

    // ---- legacy ------------------------------------------------------------
    gsThinShellAssembler<3,real_t,true> AL(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsMatrix<real_t> r0 = AL.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleMass());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsMatrix<real_t> r1 = AL.rhs();
    gsInfo << "[MassRestoresDirichlet] legacy : |r0| = "<<r0.norm()
           << " , |r(after mass)| = "<<r1.norm()<<" , |r(homogeneous)| = "<<rHom.norm()
           << "\n   REL dev = "<<relDiff(r1,r0)
           << " , REL separation from the control = "<<relDiff(r0,rHom)<<"\n";
    CHECK(relDiff(r0,rHom) > 1e-2);               // the control is genuinely apart
    CHECK(relDiff(r1,r0)  <= TOL_EXACT);          // FAILS pre-fix

    // the LUMPED path takes the same restore
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assembleMass(true));
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    CHECK(relDiff(AL.rhs(),r0) <= TOL_EXACT);

    // ---- gsThinShellAssembler2, identical defect ---------------------------
    gsThinShellAssembler2<3,real_t,true> A2(f.mp,f.dbasis,bc,f.force,f.t,f.law,f.rho);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsMatrix<real_t> s0 = A2.rhs();
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assembleMass());
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A2.assemble());
    gsMatrix<real_t> s1 = A2.rhs();
    gsInfo << "  new    : |s0| = "<<s0.norm()<<" , |s(after mass)| = "<<s1.norm()
           << " , REL dev = "<<relDiff(s1,s0)<<"\n";
    CHECK(relDiff(s0,rHom) > 1e-2);
    CHECK(relDiff(s1,s0)  <= TOL_EXACT);
}

// ---------------------------------------------------------------------
// D4 : point masses must reach the mass matrix
// ---------------------------------------------------------------------
//
// setPointMass() is public API on the legacy assembler, but the _applyMass() call
// that consumes it lived ONLY inside the commented-out block, so the point masses
// were silently dropped. This is a REAL defect, not a deliberate simplification:
// the setter, the member and the whole _applyMass() body were all still there and
// reachable, only never called. Restored.
//
// The ORACLE is exact. _applyMass deposits bVals(k)*bVals(l)*value into each of
// the d component blocks, and sum_k sum_l N_k N_l = (sum_k N_k)^2 = 1 by partition
// of unity at an interior point where every active is free, so the grand sum of
// the mass matrix must grow by EXACTLY d*value.
//
// gsThinShellAssembler2 has NO point-mass API at all (no setPointMass, no m_pMass,
// no _applyMass), so D4 does not exist there and there is nothing to gate.
TEST(MassPointMasses_legacy)
{
    MassFixture f;
    gsBoundaryConditions<real_t> bc;
    roofBCs(f.mp,bc);

    const real_t pmass = 7.5;

    // ---- without point masses ----------------------------------------------
    gsThinShellAssembler<3,real_t,true> A0(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A0.assembleMass());
    const real_t total0 = grandSum(A0.massMatrix());

    // ---- with ONE interior point mass --------------------------------------
    gsPointLoads<real_t> pMass;
    gsVector<real_t> pt(2);  pt << 0.5, 0.5;      // interior: all actives are free
    gsVector<real_t> pv(1);  pv << pmass;         // _applyMass asserts size 1
    pMass.addLoad(pt,pv,0,true);

    gsThinShellAssembler<3,real_t,true> A1(f.mp,f.dbasis,bc,f.force,&f.adapter);
    A1.setPointMass(pMass);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A1.assembleMass());
    gsSparseMatrix<real_t> M1 = A1.massMatrix();  // copy OUT
    const real_t total1 = grandSum(M1);

    gsInfo << "[MassPointMasses_legacy] total mass "<<total0<<" -> "<<total1
           << " , delta = "<<(total1-total0)<<"   (must be d*value = "<<3*pmass
           << " ; pre-fix: exactly 0)\n";
    CHECK(math::abs(total0) > 1e-3);
    CHECK(math::abs((total1-total0) - 3*pmass) <= 1e-9*(1.0+math::abs(total0)));

    // ---- and the LUMPED path carries them too ------------------------------
    // That is precisely what the ordering "_applyMass() BEFORE the row-sum" buys:
    // the point-mass block is a CONSISTENT rank-one contribution, so lumping after
    // it folds the point masses into the diagonal instead of dropping them (the
    // "To do: add point masses in lumped case" of the old dead block).
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)A1.assembleMass(true));
    gsSparseMatrix<real_t> ML1 = A1.massMatrix();
    gsInfo << "  lumped with point masses: max|offdiag| = "<<maxOffDiag(ML1)
           << " , total = "<<grandSum(ML1)<<"\n";
    CHECK(maxOffDiag(ML1) == 0.0);
    CHECK(math::abs(grandSum(ML1)-total1) <= 1e-12*math::abs(total1));

    // ---- WHICH diagonal, WITH the point mass present (task 62) --------------
    // The grand-sum gate above pins only that d*value ARRIVED, not WHERE. This
    // pins the destination per entry: the point-mass block is deposited before
    // the lumping, so each lumped diagonal entry must still equal its own row sum
    // of the consistent matrix INCLUDING that block. Note that a point mass is
    // smeared over the actives in proportion to N_i -- the B-spline basis is not
    // interpolatory, so it never lands entirely on one dof, even at a "node".
    const gsMatrix<real_t> rs = consistentRowSums(M1);
    real_t maxDev = 0.0;
    for (index_t i = 0; i != M1.rows(); ++i)
    {
        maxDev = math::max(maxDev, math::abs(ML1.coeff(i,i)-rs(i,0))/(1.0+math::abs(rs(i,0))));
        CHECK(math::abs(ML1.coeff(i,i) - rs(i,0)) <= TOL_EXACT*(1.0 + math::abs(rs(i,0))));
    }
    gsInfo << "  per-entry: max_i |M_l(i,i) - rowsum_i(M_c)| (rel) = "<<maxDev
           << " , the row sums spread over ["<<rs.minCoeff()<<", "<<rs.maxCoeff()
           << "]\n";
}

// ---------------------------------------------------------------------
// The OUT-OF-DOMAIN point mass (task 62)
// ---------------------------------------------------------------------
//
// Task 44 made _applyMass live again, and its reviewer then measured that a point
// mass at a parametric point OUTSIDE the domain silently NaN-ed the ENTIRE mass
// matrix while assembleMass() still returned Success -- silent corruption reported
// as success, and every downstream consumer branches on that status. _applyMass now
// rejects such a point with a GISMO_ENSURE (a caller-supplied value on a public API,
// evaluated once per point mass, so it is NOT NDEBUG-gated: unlike
// TEST(MassStatus_legacy) and TEST(MassFailurePath_legacy) this test keeps its lever
// in a Release build). The std::runtime_error is thrown INSIDE assembleMass's try,
// so the caller sees AssemblyError rather than the exception -- the same shape as
// TEST(MassFailurePath_legacy).
//
// TASK 63 added the CLAMP: the guard admits points within
// tol = 1e3*eps*extent (~2.2e-13 in double) of the domain, and inside that band the
// point mass was NOT harmless -- at u = 1+1e-14 the tensor basis returns a SHIFTED
// active set on which every N vanishes, so on a polynomial basis the k x l loop
// deposits nothing (mass silently lost) and on a RATIONAL one -- which this roof
// fixture is -- gsRationalBasis additionally divides by sum_i w_i N_i = 0 and the
// WHOLE mass matrix goes NaN, still under Success. Measured in the task-63 poison
// round: the two in-band deltas below read -nan, not 0. So task 62's tolerance band
// re-opened its own defect over a 2.2e-13-wide sliver. An accepted point is now
// clamped onto the domain before it is evaluated.
//
// NON-VACUITY has four halves and ALL of them are asserted:
//   (1) the out-of-domain point is REJECTED -- pre-guard this read Success;
//   (2) a point ON the boundary of the domain is still ACCEPTED and yields a finite
//       mass. A guard that rejected valid input would be worse than the NaN it
//       replaces, so the closure is gated explicitly -- and gated on the DELTA, not
//       merely on finiteness, so that a silent drop cannot pass it;
//   (3) a point INSIDE the tolerance band deposits the SAME mass as the boundary
//       point itself -- pre-clamp its delta read -nan here, under Success;
//   (4) a point OUTSIDE the band is still rejected, i.e. the clamp did not turn the
//       guard into a no-op.
// The domain is READ FROM THE BASIS and printed, not assumed to be [0,1]^2.
//
// CONTRACT WIDENING recorded here (task 62's reviewer, note N4): a point mass whose
// VALUE is zero at an out-of-domain point used to be a harmless no-op, because the
// "value != 0" test inside _applyMass gated the write and the NaN basis values were
// never consumed. It now raises AssemblyError like any other out-of-domain point.
// So the guard is NOT confined to the path that previously produced a NaN matrix.
TEST(MassPointMassOutOfDomain_legacy)
{
    MassFixture f;
    gsBoundaryConditions<real_t> bc;
    roofBCs(f.mp,bc);

    const gsMatrix<real_t> supp = f.dbasis.basis(0).support();
    gsInfo << "[MassPointMassOutOfDomain_legacy] parameter domain of patch 0 : ["
           << supp(0,0)<<", "<<supp(0,1)<<"] x ["<<supp(1,0)<<", "<<supp(1,1)<<"]\n";

    gsVector<real_t> pv(1);  pv << 7.5;           // _applyMass asserts size 1

    // ---- (1) a point far OUTSIDE the domain is rejected ---------------------
    gsPointLoads<real_t> outside;
    gsVector<real_t> ptOut(2);  ptOut << supp(0,1)+2.5, supp(1,0)-2.0;
    outside.addLoad(ptOut,pv,0,true);

    gsThinShellAssembler<3,real_t,true> Aout(f.mp,f.dbasis,bc,f.force,&f.adapter);
    Aout.setPointMass(outside);
    const int st = (int)Aout.assembleMass();
    gsInfo << "  point ("<<ptOut.transpose()<<") -> status "<<st
           << "   (must be AssemblyError = "
           << (int)ThinShellAssemblerStatus::AssemblyError
           << " ; pre-guard: Success, with a NaN mass matrix)\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,st);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,(int)Aout.status());

    // ... and the CORRUPTION half of the defect: pre-guard the whole matrix read
    // NaN. Inspecting m_mass after a caught failure is legitimate HERE, unlike in
    // TEST(MassFailurePath_legacy): the guard fires before _applyMass writes
    // anything, so what is left behind is the clean consistent matrix -- which is
    // exactly the claim, "the bad input reached nothing".
    const real_t tOut = grandSum(Aout.massMatrix());
    gsInfo << "  the mass matrix left behind: total = "<<tOut
           << "   (pre-guard: NaN throughout)\n";
    CHECK(0 == math::isnan(tOut));

    // ---- (2) ... and a point ON the boundary is NOT -------------------------
    gsPointLoads<real_t> onEdge;
    gsVector<real_t> ptEdge(2);  ptEdge << supp(0,1), 0.5*(supp(1,0)+supp(1,1));
    onEdge.addLoad(ptEdge,pv,0,true);

    gsThinShellAssembler<3,real_t,true> Aedge(f.mp,f.dbasis,bc,f.force,&f.adapter);
    Aedge.setPointMass(onEdge);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)Aedge.assembleMass());
    const real_t tEdge = grandSum(Aedge.massMatrix());
    gsInfo << "  boundary point ("<<ptEdge.transpose()<<") -> Success , total mass = "
           << tEdge<<"\n";
    CHECK(0 == math::isnan(tEdge));
    CHECK(math::abs(tEdge) > 1e-3);

    // The baseline the two deltas below are measured against.
    int stBase = -1;
    const real_t tBase = roofMassTotal(f,bc,nullptr,stBase);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stBase);

    // WHAT CORRECT CODE RETURNS AT THE SAME STATE, before any number below is read
    // as a signature. On the EAST edge roofBCs eliminates components 1 and 2 and a
    // clamped B-spline is interpolatory at u = 1, so only the free component of the
    // single non-vanishing active contributes: the delta is value = 7.5, NOT
    // d*value. On the NORTH edge nothing is eliminated, every active is free and
    // partition of unity gives sum_k sum_l N_k N_l = 1 per component, i.e. the full
    // d*value = 22.5 -- the same oracle TEST(MassPointMasses_legacy) uses.
    const real_t dEdge = tEdge - tBase;
    gsInfo << "  east-edge delta = "<<dEdge<<"   (value = 7.5: east eliminates"
           << " components 1 and 2)\n";
    CHECK(math::abs(dEdge - 7.5) <= TOL_EXACT*(1.0+math::abs(tBase)));

    int stN = -1;
    const gsPointLoads<real_t> north = onePointMass(0.5,supp(1,1),7.5);
    const real_t dNorth = roofMassTotal(f,bc,&north,stN) - tBase;
    gsInfo << "  north-edge delta = "<<dNorth<<"   (must be d*value = 22.5)\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stN);
    CHECK(math::abs(dNorth - 3*7.5) <= TOL_EXACT*(1.0+math::abs(tBase)));

    // ---- (3) a point INSIDE the tolerance band is CLAMPED, not dropped -------
    // tol = 1e3*eps*extent ~ 2.2e-13 here, so 1e-14 past the boundary is admitted.
    // Pre-clamp BOTH deltas below read -nan with status Success (poison round B).
    // The two absolute gates above (dEdge == 7.5, dNorth == 22.5) are what makes
    // these two relative ones non-vacuous: a poison that rejects the boundary point
    // collapses BOTH sides of the comparison to 0 and it would pass on nothing.
    const real_t inBand = 1e-14;
    int stE2 = -1, stN2 = -1;
    const gsPointLoads<real_t> edgeBand  = onePointMass(supp(0,1)+inBand,ptEdge(1),7.5);
    const gsPointLoads<real_t> northBand = onePointMass(0.5,supp(1,1)+inBand,7.5);
    const real_t dEdgeBand  = roofMassTotal(f,bc,&edgeBand ,stE2) - tBase;
    const real_t dNorthBand = roofMassTotal(f,bc,&northBand,stN2) - tBase;
    gsInfo << "  in-band (+"<<inBand<<") : east delta = "<<dEdgeBand
           << " (must be "<<dEdge<<") , north delta = "<<dNorthBand
           << " (must be "<<dNorth<<")   [pre-clamp: -nan and -nan, under Success]\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stE2);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stN2);
    // Relative form, not CHECK_EQUAL: these are separate assembleMass() calls and
    // the consistent part accumulates under an unordered #pragma omp atomic.
    CHECK(math::abs(dEdgeBand  - dEdge ) <= TOL_EXACT*(1.0+math::abs(tBase)));
    CHECK(math::abs(dNorthBand - dNorth) <= TOL_EXACT*(1.0+math::abs(tBase)));

    // ---- (4) ... but a point OUTSIDE the band still is rejected -------------
    // 1e-9 is ~4000x the band and still 8 orders below the element size 1/8, so
    // this pins that the clamp did not silently widen the guard.
    int stFar = -1;
    const gsPointLoads<real_t> justOut = onePointMass(supp(0,1)+1e-9,ptEdge(1),7.5);
    const real_t tFar = roofMassTotal(f,bc,&justOut,stFar);
    gsInfo << "  out-of-band (+1e-9) -> status "<<stFar<<"   (must be AssemblyError = "
           << (int)ThinShellAssemblerStatus::AssemblyError<<") , total = "<<tFar<<"\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,stFar);
    CHECK(0 == math::isnan(tFar));
}

// ---------------------------------------------------------------------
// The SPACE basis, not the integration basis (task 63)
// ---------------------------------------------------------------------
//
// _applyMass resolves its actives through m_mapper, which comes from m_space, which
// _initialize() builds from *m_spaceBasis (gsThinShellAssembler.hpp:259). m_basis is
// only the INTEGRATION basis (:248). The two are the same object until the public
// setSpaceBasis() is called, and _applyMass used to compute its actives from
// m_basis.front().basis(patch) -- numbering them in one basis and resolving them in
// the other. _applyLoads has carried a dynamic_cast ladder over *m_spaceBasis for
// exactly this since forever; _applyMass had none.
//
// MEASURED before the repair, with m_basis = the roof basis and
// setSpaceBasis(a degree-elevated copy): the point mass landed on a completely
// different dof set and only 15.9375 of its 22.5 units arrived (the rest resolved
// onto ELIMINATED dofs and was discarded by is_free_index) -- silently, under
// Success. The reverse split lost even more (7.5 of 22.5) and read the mapper's
// per-patch dof array past its own length. Task 66 re-measured both numbers on
// the degree-elevated fixture used below, by reverting _applyMass to its
// pre-task-63 m_basis.front().basis(patch) form and rebuilding: they are
// IDENTICAL (4.24471 / 15.9375 and 3.75 / 7.5) to the ones the original
// uniformRefine fixture produced, at OMP_NUM_THREADS = 1, 4 AND 8.
//
// THE ORACLE IS A CONTROL ASSEMBLER, not a hand-derived number: the same SPACE with
// no divergence at all. Its mapper is identical (same basis, same BCs, same
// continuity), so the point-mass INCREMENT M(with) - M(without) must agree
// ENTRY-FOR-ENTRY. The increment is used rather than the full matrix because it
// ISOLATES the point-mass contribution: it removes the O(1e3) background mass, so
// the gate is scaled to |D_ref| ~ 4.3 instead of ||M|| ~ 1e2, and a mis-placed 7.5
// cannot hide inside it. (Under the ORIGINAL uniformRefine fixture the increment
// was also the only comparable quantity, because the two configurations then
// integrated on different meshes. That is no longer so -- with the degreeElevate
// fixture below the two FULL matrices are bit-identical, measured; see the note
// on task 66. The increment is kept for the scaling argument alone.)
// NOTE that the grand sum alone is blind to this defect whenever the mis-resolved
// dofs all happen to be free, so the entry-wise comparison is the load-bearing gate
// and the grand sum is only a second opinion.
//
// HOW THE DIVERGENT BASIS IS BUILT, AND WHY IT IS *NOT* uniformRefine (task 66)
// -----------------------------------------------------------------------------
// The two bases must differ in their DOF NUMBERING (that is what the defect is
// about) while sharing the SAME MESH. degreeElevate() does exactly that here:
// 100 -> 324 functions on the same 64 elements. uniformRefine(), which this test
// originally used, also gives 324 -- but on 256 elements, and that configuration
// is NOT one gsExprAssembler supports:
//
//   gsExprAssembler::assemble() sets SAME_ELEMENT (gsExprAssembler.h:1227), under
//   which gsFunctionSet::compute() takes the active set ONCE, at the FIRST
//   quadrature point (gsFunctionSet.hpp:186-189), for the whole element. With an
//   integration mesh COARSER than the space basis' mesh that assumption is false:
//   MEASURED on this fixture, in ALL 64 integration elements the space basis'
//   actives at the first and last quadrature point differ. Two consequences, both
//   measured at OMP_NUM_THREADS=1 (i.e. not a threading artefact):
//     (i)  the assembled mass matrix is silently WRONG --
//          ||M_split - M_true||_F / ||M_true||_F = 0.628 , total mass 1245.66
//          against the true 1258.09, under Success;
//     (ii) the sparsity pattern MISSED 4160 of the (i,j) entries the assembly
//          loop actually touches, so _eval::push had to STRUCTURALLY INSERT into
//          the gsFiberMatrix columns. That insert is protected by nothing -- the
//          "#pragma omp atomic" at gsExprAssembler.h:739 covers the "+=" only,
//          while _pattern::push takes a per-column omp lock
//          (gsExprAssembler.h:854-864) precisely because insertion is not thread
//          safe.
//          HISTORY, since the citations above are to a MOVING file: at the time
//          this was measured, _pattern::push sampled the space basis' actives at
//          ONE point of the element (the element centre). Task 67 changed that:
//          _computePattern now maps the SAME quadrature rule the assembly loop
//          uses (gsExprAssembler.h:1053-1056) and _pattern::push UNIONS the
//          actives over all of its points (the "for k" at :837-868), so the
//          volume pattern is no longer incomplete for a fine space basis. The
//          measured numbers here therefore describe the PRE-67 behaviour and are
//          kept as the record of the defect; the mesh-coincidence requirement
//          below is what this test still enforces.
//   From ~6 threads the concurrent inserts collide and trip Eigen's
//   "you cannot insert an element that already exists" assert inside assembleMass,
//   which surfaces as AssemblyError and an EMPTY massMatrix(). That is exactly why
//   this test used to pass at 4 threads and fail at 8 (0/10, 3/10, 5/10 and 10/10
//   runs failed at 4, 6, 7 and 8 threads).
// With degreeElevate the meshes coincide, the SAME_ELEMENT contract holds, the
// pattern is complete (0 missing entries, measured), M_split is BIT-IDENTICAL to
// M_true, and the test is thread-count independent. The dof-mapping coverage is
// unchanged: the actives are still numbered in one basis and resolved in the
// other's mapper, in both directions. It is unchanged STRUCTURALLY, not just
// empirically: the refined and the elevated basis are both 18x18 tensor bases of
// size 324 with IDENTICAL west and east boundary index sets (measured), so the two
// mappers are the same index set with the same free/eliminated partition -- which
// is why reverting _applyMass to its pre-task-63 form reproduces the elevated
// fixture's failure numbers BIT-IDENTICALLY to the refined one's.
// The equal-mesh CHECK below is load bearing -- it stops this fixture from being
// "simplified" back into the unsupported configuration.
TEST(MassPointMassSpaceBasis_legacy)
{
    MassFixture f;
    gsBoundaryConditions<real_t> bc;
    roofBCs(f.mp,bc);

    const real_t pmass = 7.5;
    const gsPointLoads<real_t> pm = onePointMass(0.5,0.5,pmass);   // interior

    gsMultiBasis<real_t> elev(f.dbasis);
    elev.degreeElevate();
    gsInfo << "[MassPointMassSpaceBasis_legacy] basis sizes: coarse "
           << f.dbasis.basis(0).size()<<" , elevated "<<elev.basis(0).size()
           << "   (elements "<<f.dbasis.basis(0).numElements()<<" vs "
           << elev.basis(0).numElements()<<")\n";
    CHECK(elev.basis(0).size() > f.dbasis.basis(0).size());   // the split is real
    // ... and the integration meshes MUST agree, see the note above.
    CHECK_EQUAL(f.dbasis.basis(0).numElements(),elev.basis(0).numElements());

    // The point-mass increment of a fresh pair of assemblers. @a intBasis is the
    // integration basis (m_basis); @a spaceBasis is passed to setSpaceBasis unless
    // it is null, in which case the two coincide and there is nothing to diverge.
    struct Local
    {
        static gsSparseMatrix<real_t> increment(MassFixture & f,
                                                const gsBoundaryConditions<real_t> & bc,
                                                const gsMultiBasis<real_t> & intBasis,
                                                const gsMultiBasis<real_t> * spaceBasis,
                                                const gsPointLoads<real_t> & pm,
                                                int & status)
        {
            gsThinShellAssembler<3,real_t,true> A0(f.mp,intBasis,bc,f.force,&f.adapter);
            if (spaceBasis!=nullptr) A0.setSpaceBasis(*spaceBasis);
            status = (int)A0.assembleMass();
            gsSparseMatrix<real_t> M0 = A0.massMatrix();      // copy OUT

            gsThinShellAssembler<3,real_t,true> A1(f.mp,intBasis,bc,f.force,&f.adapter);
            if (spaceBasis!=nullptr) A1.setSpaceBasis(*spaceBasis);
            A1.setPointMass(pm);
            const int st1 = (int)A1.assembleMass();
            if (st1 != (int)ThinShellAssemblerStatus::Success) status = st1;
            return gsSparseMatrix<real_t>(A1.massMatrix() - M0);
        }
    };

    // ---- direction 1 : integration COARSE, space ELEVATED -------------------
    // The actives are numbered in the coarse basis (0..99) and resolved in the
    // elevated mapper: in range, but a different set of functions.
    int stR = -1, stS = -1;
    const gsSparseMatrix<real_t> Dref1   = Local::increment(f,bc,elev    ,nullptr,pm,stR);
    const gsSparseMatrix<real_t> Dsplit1 = Local::increment(f,bc,f.dbasis,&elev  ,pm,stS);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stR);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stS);
    gsInfo << "  coarse integration / elevated space: |D_split - D_oracle| = "
           << absDiff(Dsplit1,Dref1)<<"   grand sums "<<grandSum(Dsplit1)<<" vs "
           << grandSum(Dref1)<<"   (pre-fix: 4.24471 , 15.9375 vs 22.5)\n";
    CHECK(math::abs(grandSum(Dref1)   - 3*pmass) <= TOL_EXACT*(1.0+3*pmass));
    CHECK(math::abs(grandSum(Dsplit1) - 3*pmass) <= TOL_EXACT*(1.0+3*pmass));
    CHECK(absDiff(Dsplit1,Dref1) <= TOL_EXACT*(1.0+Dref1.norm()));

    // ---- direction 2 : integration ELEVATED, space COARSE -------------------
    // This is the dangerous one: the actives then exceed the space's per-patch dof
    // count, so the mapper is indexed past the end of its own array.
    int stR2 = -1, stS2 = -1;
    const gsSparseMatrix<real_t> Dref2   = Local::increment(f,bc,f.dbasis,nullptr  ,pm,stR2);
    const gsSparseMatrix<real_t> Dsplit2 = Local::increment(f,bc,elev    ,&f.dbasis,pm,stS2);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stR2);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stS2);
    gsInfo << "  elevated integration / coarse space: |D_split - D_oracle| = "
           << absDiff(Dsplit2,Dref2)<<"   grand sums "<<grandSum(Dsplit2)<<" vs "
           << grandSum(Dref2)<<"   (pre-fix: 3.75 , 7.5 vs 22.5)\n";
    CHECK(math::abs(grandSum(Dref2)   - 3*pmass) <= TOL_EXACT*(1.0+3*pmass));
    CHECK(math::abs(grandSum(Dsplit2) - 3*pmass) <= TOL_EXACT*(1.0+3*pmass));
    CHECK(absDiff(Dsplit2,Dref2) <= TOL_EXACT*(1.0+Dref2.norm()));
}

// ---------------------------------------------------------------------
// The FAILURE path of assembleMass itself
// ---------------------------------------------------------------------
//
// TEST(MassStatus_legacy) breaks assembleVector, not assembleMass, so the
// routine's own catch branch -- and, more importantly, the Dirichlet RESTORE that
// task 44 put AFTER the try/catch -- would otherwise never be executed by any
// test. The code comment there claims the restore "runs on the error path too";
// this is the assertion behind that claim.
//
// The lever is _applyMass's own guard: a point mass must carry a SCALAR value
// (GISMO_ASSERT at gsThinShellAssembler.hpp:1491), so a 2-vector throws INSIDE
// the try. NOTE that this build carries no -DNDEBUG; on a build that strips
// assertions the guard is absent and this test would read Success instead, which
// is a loud failure rather than a silent pass.
//
// Three claims, and the middle one is the real unknown:
//   (1) the failure is REPORTED as AssemblyError, not thrown past the routine;
//   (2) _assembleDirichlet() on a just-cleanUp()ed expression assembler does not
//       throw -- simply reaching the line after assembleMass() proves it;
//   (3) the space is restored DESPITE the failure, so the next assemble() still
//       returns the NON-homogenized rhs.
TEST(MassFailurePath_legacy)
{
    MassFixture f;

    gsFunctionExpr<real_t> gD("0.02",3);
    gsBoundaryConditions<real_t> bc;   roofBCs(f.mp,bc,&gD);
    gsBoundaryConditions<real_t> bcH;  roofBCs(f.mp,bcH,nullptr);

    gsThinShellAssembler<3,real_t,true> ALh(f.mp,f.dbasis,bcH,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)ALh.assemble());
    gsMatrix<real_t> rHom = ALh.rhs();            // copy OUT

    gsThinShellAssembler<3,real_t,true> AL(f.mp,f.dbasis,bc,f.force,&f.adapter);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsMatrix<real_t> r0 = AL.rhs();

    // an ILLEGAL point mass: _applyMass wants a scalar value
    gsPointLoads<real_t> bad;
    gsVector<real_t> pt(2);  pt << 0.5, 0.5;
    gsVector<real_t> bv(2);  bv << 1.0, 2.0;
    bad.addLoad(pt,bv,0,true);
    AL.setPointMass(bad);

    const int st = (int)AL.assembleMass();        // (1)
    gsInfo << "[MassFailurePath_legacy] broken assembleMass -> status "<<st
           << "   (must be AssemblyError = "
           << (int)ThinShellAssemblerStatus::AssemblyError<<")\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,st);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,(int)AL.status());

    // (2) reaching here at all is the claim; (3) is asserted below. m_mass is
    // deliberately NOT inspected -- after a caught failure it legitimately holds a
    // half-built object, which is pre-existing catch semantics everywhere in this
    // class.
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,(int)AL.assemble());
    gsMatrix<real_t> r1 = AL.rhs();
    gsInfo << "  after the FAILED mass: |r| = "<<r1.norm()<<" vs |r0| = "<<r0.norm()
           << " , |r(homogeneous)| = "<<rHom.norm()
           << " ; REL dev = "<<relDiff(r1,r0)<<"\n";
    CHECK(relDiff(r0,rHom) > 1e-2);               // the control is genuinely apart
    CHECK(relDiff(r1,r0)  <= TOL_EXACT);          // (3)
}

// =====================================================================
// TASK 61 / D5 : _applyMass on a MULTIPATCH
// =====================================================================
//
// Two pieces of the hardened _applyMass rested on SCRATCH PROBES only -- task 63's
// report measured both, but nothing in the suite gated either:
//   (1) a point mass on patch > 0 WORKS. Before task 63 the routine resolved its
//       basis as m_basis.front().basis(patch), i.e. it asked PATCH 0's basis for a
//       piece it does not have, and gsBasis::piece()'s GISMO_ENSURE(0==k)
//       (src/gsCore/gsBasis.h:105-109) threw for EVERY patch != 0. Dispatching over
//       *m_spaceBasis (a gsMultiBasis, which indexes patches properly) removed that.
//       Task 63 probed a DISCONNECTED two-patch fixture; the fixture here is GLUED,
//       so the interface dofs are genuinely shared and the mapper path the fix
//       touched is the one exercised.
//   (2) an OUT-OF-RANGE patch index is rejected by the new GISMO_ENSURE at
//       gsThinShellAssembler.hpp:1500.
//
// ASSERT ON THE STATUS, NOT WITH CHECK_THROW: the GISMO_ENSURE throws INSIDE
// assembleMass's try, whose catch(...) (:1575 region) converts it into
// ThinShellAssemblerStatus::AssemblyError. A CHECK_THROW here would fail against
// correct code.

/// Two unit squares glued along x = 1, embedded in 3D, degree 2, 4x4 elements each.
/// FLAT on purpose: the claim under test is about dof BOOKKEEPING across patches,
/// and a flat mid-surface keeps the mass oracle a pure partition-of-unity statement.
gsMultiPatch<real_t> makeGluedStrip()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1,0,0) );   // [0,1] x [0,1]
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1,1,0) );   // [1,2] x [0,1]
    mp.embed(3);
    mp.degreeElevate();       // degree 2: the active set at a point is wider than 1
    mp.uniformRefine();
    mp.uniformRefine();
    mp.computeTopology();     // <- this is what GLUES them
    return mp;
}

/// Clamp the two OUTER edges (patch 0 west, patch 1 east) in all three components.
/// The interface edge and the north/south edges stay FREE, which is what makes the
/// partition-of-unity oracle exact at the points used below.
void stripBCs(const gsMultiPatch<real_t> & mp, gsBoundaryConditions<real_t> & bc)
{
    bc.setGeoMap(mp);
    for (int c = 0; c != 3; ++c)
    {
        bc.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0,false,c);
        bc.addCondition(1,boundary::east,condition_type::dirichlet,nullptr,0,false,c);
    }
}

struct StripMassFixture
{
    gsMultiPatch<real_t>          mp;
    gsMultiBasis<real_t>          dbasis;
    gsConstantFunction<real_t>    force;
    gsFunctionExpr<real_t>        t, rho;
    gsLinearMaterial<real_t>      law;
    gsMaterialMatrix3D<3,real_t>  adapter;

    StripMassFixture()
    :
    mp(makeGluedStrip()), dbasis(mp), force(roofForce(),3),
    t(util::to_string(THICK),3), rho("1.0",3),
    law(E_MOD,NU,2), adapter(mp,t,rho,law)
    { }
};

/// ONE parametric point mass of magnitude @a value at (@a u, @a v) on @a patch.
gsPointLoads<real_t> onePointMassOn(index_t patch, real_t u, real_t v, real_t value)
{
    gsPointLoads<real_t> pl;
    gsVector<real_t> pt(2);  pt << u, v;
    gsVector<real_t> pv(1);  pv << value;
    pl.addLoad(pt,pv,patch,true);
    return pl;
}

/// Set the assembler's interface continuity. MEASURED, and NOT what the name of a
/// glued fixture suggests: gsFeSpace::setup (src/gsExpressions/gsFeSpace.h:155-162)
/// calls matchInterface ONLY when interfaceCont() == 0, and gsThinShellAssembler
/// defaults "Continuity" to -1 (gsThinShellAssembler.hpp:209). So a multipatch whose
/// TOPOLOGY has an interface still gets a mapper with INDEPENDENT dofs on either
/// side of it unless this is set. setOptions() is the route because it is what
/// re-runs _initialize() when the value changes (:229-240).
void setContinuity(gsThinShellAssemblerBase<real_t> & A, index_t cont)
{
    gsOptionList o;
    o.addInt("Continuity","Set the continuity for the space",cont);
    A.setOptions(o);
}

/// M(with @a pm) - M(without), from a FRESH pair of assemblers on the strip, both at
/// interface continuity @a cont. The increment is used rather than the matrix itself
/// so that the consistent part cancels exactly and the point-mass block is left.
gsSparseMatrix<real_t> stripMassIncrement(StripMassFixture & f,
                                          const gsBoundaryConditions<real_t> & bc,
                                          const gsPointLoads<real_t> & pm,
                                          int & status,
                                          index_t cont = 0)
{
    gsThinShellAssembler<3,real_t,true> A0(f.mp,f.dbasis,bc,f.force,&f.adapter);
    setContinuity(A0,cont);
    status = (int)A0.assembleMass();
    gsSparseMatrix<real_t> M0 = A0.massMatrix();          // copy OUT

    gsThinShellAssembler<3,real_t,true> A1(f.mp,f.dbasis,bc,f.force,&f.adapter);
    setContinuity(A1,cont);
    A1.setPointMass(pm);
    const int st1 = (int)A1.assembleMass();
    if (st1 != (int)ThinShellAssemblerStatus::Success) status = st1;
    return gsSparseMatrix<real_t>(A1.massMatrix() - M0);
}

/// Which rows carry any of the deposited mass. This is the PLACEMENT observable --
/// the grand sum is blind to it (task 63's reviewer, and the reason that note asked
/// for an independent placement assertion here).
std::vector<bool> touchedRows(const gsSparseMatrix<real_t> & M)
{
    // Threshold on the SCALE of the increment in units of eps (the TOL_EXACT idiom
    // of :142), NOT a fixed 1e-14: the increment is a difference of two separately
    // assembled consistent parts, whose cancellation residue is ~eps*|M| -- ~1e-18
    // in double but ~1e-9 at real_t = float, i.e. ABOVE any fixed 1e-14, which would
    // report spurious touched rows there. The smallest genuine entry deposited here
    // is 7.5*0.25*0.25 ~ 0.47, so this keeps ~3 orders of headroom even at float.
    const real_t thresh = TOL_EXACT*(1.0 + M.norm());
    std::vector<bool> r((size_t)M.rows(),false);
    for (index_t k = 0; k != M.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(M,k); it; ++it)
            if (math::abs(it.value()) > thresh) r[(size_t)it.row()] = true;
    return r;
}

index_t countTrue(const std::vector<bool> & v)
{
    index_t n = 0;
    for (size_t i = 0; i != v.size(); ++i) if (v[i]) ++n;
    return n;
}

// ---------------------------------------------------------------------
// (1) a point mass on patch > 0, on a GLUED multipatch
// ---------------------------------------------------------------------
//
// THREE claims, each with its own oracle, and none of them is the code's own output:
//
//   (a) ANALYTIC TOTAL. At a point where every active dof is free, the deposited
//       block is value * N_k N_l per component and sum_k sum_l N_k N_l =
//       (sum_k N_k)^2 = 1 by partition of unity, so the grand sum must grow by
//       EXACTLY d*value = 22.5 -- on patch 1 just as on patch 0. Pre-fix the patch-1
//       call did not produce a wrong number, it produced AssemblyError.
//
//   (b) PLACEMENT, independently of any total. The patch-0 and patch-1 interior
//       increments must touch DISJOINT, non-empty row sets: with degree 2 and four
//       elements the actives at u = 0.5 are strictly interior to their own patch and
//       never reach the interface functions. A routine that ignored the patch index
//       and always deposited on patch 0 satisfies (a) exactly and fails this.
//
//   (c) THE GLUE ITSELF. The SAME physical point on the shared edge, addressed as
//       (patch 0, u = 1) and as (patch 1, u = 0), must give ENTRY-FOR-ENTRY the same
//       increment: a clamped B-spline is interpolatory at its end knot, so only the
//       interface functions carry a nonzero value there, and at Continuity = 0 the
//       mapper has identified those dofs across the interface. This is the assertion
//       that actually needs the patches to be GLUED.
//
//       WHAT CORRECT CODE DOES AT THE OTHER SETTING, checked before this number was
//       read as a signature: at the DEFAULT Continuity = -1 the mapper does NOT
//       match the interface at all (gsFeSpace.h:155), so the two sides land on
//       DISJOINT dofs and |DL - DR| is 9.18559 = sqrt(2)*|DL| -- correct behaviour
//       for that setting, not a defect. Both settings are asserted below, and the
//       contrast between them is what makes the equality claim non-trivial.
TEST(MassPointMassMultipatchGlued_legacy)
{
    StripMassFixture f;
    gsBoundaryConditions<real_t> bc;
    stripBCs(f.mp,bc);

    // The fixture is what the test claims it is.
    CHECK_EQUAL(2,(int)f.mp.nPatches());
    CHECK_EQUAL(1,(int)f.mp.nInterfaces());
    gsInfo << "[MassPointMassMultipatchGlued_legacy] "<<f.mp.nPatches()<<" patches, "
           << f.mp.nInterfaces()<<" interface(s), basis sizes "
           << f.dbasis.basis(0).size()<<" + "<<f.dbasis.basis(1).size()<<"\n";

    const real_t pmass = 7.5;
    const real_t oracle = 3*pmass;      // d * value

    // ---- (a) + (b) : interior of patch 1, and of patch 0 --------------------
    int st0 = -1, st1 = -1;
    const gsPointLoads<real_t> pm0 = onePointMassOn(0,0.5,0.5,pmass);
    const gsPointLoads<real_t> pm1 = onePointMassOn(1,0.5,0.5,pmass);
    const gsSparseMatrix<real_t> D0 = stripMassIncrement(f,bc,pm0,st0);
    const gsSparseMatrix<real_t> D1 = stripMassIncrement(f,bc,pm1,st1);

    gsInfo << "  interior patch 0 -> status "<<st0<<" , grand sum "<<grandSum(D0)<<"\n"
           << "  interior patch 1 -> status "<<st1<<" , grand sum "<<grandSum(D1)
           << "   (both must be d*value = "<<oracle
           << " ; pre-task-63 the patch-1 call was AssemblyError)\n";
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,st0);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,st1);
    CHECK(math::abs(grandSum(D0) - oracle) <= TOL_EXACT*(1.0+oracle));
    CHECK(math::abs(grandSum(D1) - oracle) <= TOL_EXACT*(1.0+oracle));

    // (b) placement: disjoint, non-empty row sets
    const std::vector<bool> r0 = touchedRows(D0), r1 = touchedRows(D1);
    index_t shared = 0;
    for (size_t i = 0; i != r0.size(); ++i) if (r0[i] && r1[i]) ++shared;
    gsInfo << "  rows touched: patch 0 -> "<<countTrue(r0)<<" , patch 1 -> "
           << countTrue(r1)<<" , SHARED -> "<<shared<<"   (must be 0: the two"
           << " interior points are on different patches)\n";
    CHECK(countTrue(r0) > 0);
    CHECK(countTrue(r1) > 0);
    CHECK_EQUAL(0,(int)shared);

    // ---- (c) : the shared edge, addressed from either side ------------------
    const gsPointLoads<real_t> pmL = onePointMassOn(0,1.0,0.5,pmass);  // east of patch 0
    const gsPointLoads<real_t> pmR = onePointMassOn(1,0.0,0.5,pmass);  // west of patch 1

    // (c1) GLUED mapper: the same physical point, so the same entries.
    int stL = -1, stR = -1;
    const gsSparseMatrix<real_t> DL = stripMassIncrement(f,bc,pmL,stL,0);
    const gsSparseMatrix<real_t> DR = stripMassIncrement(f,bc,pmR,stR,0);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stL);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stR);

    gsInfo << "  interface, Continuity = 0 (glued) : grand sums "<<grandSum(DL)
           << " and "<<grandSum(DR)<<" , |DL - DR| = "<<absDiff(DL,DR)
           << "   (the same physical point on shared dofs, so this must be 0)\n";
    // Non-vacuity: two empty increments would satisfy the comparison on nothing.
    CHECK(DL.norm() > 1.0);
    CHECK(math::abs(grandSum(DL) - oracle) <= TOL_EXACT*(1.0+oracle));
    CHECK(math::abs(grandSum(DR) - oracle) <= TOL_EXACT*(1.0+oracle));
    CHECK(absDiff(DL,DR) <= TOL_EXACT*(1.0+DL.norm()));

    // (c2) the CONTRAST at the default Continuity = -1, where the mapper leaves the
    // interface unmatched: the very same two point masses must then land on disjoint
    // dofs, i.e. |DL - DR| = sqrt(2)*|DL|. This is what tells the equality above
    // apart from "the comparison cannot distinguish anything".
    int stLu = -1, stRu = -1;
    const gsSparseMatrix<real_t> DLu = stripMassIncrement(f,bc,pmL,stLu,-1);
    const gsSparseMatrix<real_t> DRu = stripMassIncrement(f,bc,pmR,stRu,-1);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stLu);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stRu);
    gsInfo << "  interface, Continuity = -1 (unmatched) : |DLu - DRu| = "
           << absDiff(DLu,DRu)<<"   (must be sqrt(2)*|DLu| = "
           << math::sqrt(2.0)*DLu.norm()<<" : disjoint dofs)\n";
    CHECK(math::abs(absDiff(DLu,DRu) - math::sqrt(2.0)*DLu.norm())
          <= TOL_EXACT*(1.0+DLu.norm()));
}

// ---------------------------------------------------------------------
// (2) an OUT-OF-RANGE patch index
// ---------------------------------------------------------------------
//
// Making patch > 0 work removed the accidental net that gsBasis::piece()'s own
// GISMO_ENSURE used to provide, so _applyMass carries an explicit range test
// (gsThinShellAssembler.hpp:1500). It is a GISMO_ENSURE and not a GISMO_ASSERT
// deliberately: under -DNDEBUG an ASSERT would leave patch >= nPatches() falling
// through into gsMultiBasis::basis(patch) -> *m_bases[patch], an out-of-range
// dereference rather than a throw.
//
// WHAT THIS TEST DISCRIMINATES, stated honestly: it pins AssemblyError, i.e. that a
// bad patch index is REPORTED and neither succeeds nor escapes. It does NOT prove
// the ENSURE at :1500 is what produced it -- this build keeps assertions, so
// gsMultiBasis::basis()'s own GISMO_ASSERT (gsMultiBasis.h:271) is a second net that
// would give the same status if :1500 were deleted. The measured mutant it DOES
// catch is a guard that CLAMPS the index instead of rejecting it, which returns
// Success (see the report of task 61). The Release-build half of the claim -- where
// only the ENSURE survives -- is a named gap, not something this build can exercise.
TEST(MassPointMassPatchOutOfRange_legacy)
{
    StripMassFixture f;
    gsBoundaryConditions<real_t> bc;
    stripBCs(f.mp,bc);

    const real_t pmass = 7.5;

    // The CONTROL first: the last VALID patch index must succeed, otherwise
    // "AssemblyError for everything" would pass the two checks below.
    int stOk = -1;
    const gsPointLoads<real_t> good = onePointMassOn(1,0.5,0.5,pmass);
    const gsSparseMatrix<real_t> Dok = stripMassIncrement(f,bc,good,stOk);
    CHECK_EQUAL((int)ThinShellAssemblerStatus::Success,stOk);
    CHECK(math::abs(grandSum(Dok) - 3*pmass) <= TOL_EXACT*(1.0+3*pmass));

    // ... and now one PAST the end, and far past it.
    const index_t bad[2] = {2,5};                 // nPatches() == 2
    for (int i = 0; i != 2; ++i)
    {
        gsThinShellAssembler<3,real_t,true> A(f.mp,f.dbasis,bc,f.force,&f.adapter);
        setContinuity(A,0);      // the SAME configuration as the control leg above
        A.setPointMass(onePointMassOn(bad[i],0.5,0.5,pmass));
        const int st = (int)A.assembleMass();
        const real_t tot = grandSum(A.massMatrix());
        gsInfo << "  point mass on patch "<<bad[i]<<" of "<<f.mp.nPatches()
               << " -> status "<<st<<"   (must be AssemblyError = "
               << (int)ThinShellAssemblerStatus::AssemblyError
               << ") , mass left behind = "<<tot<<"\n";
        CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,st);
        CHECK_EQUAL((int)ThinShellAssemblerStatus::AssemblyError,(int)A.status());
        CHECK(0 == math::isnan(tot));
    }
}

} // SUITE

#endif // gsPhaseFieldFracture_ENABLED
