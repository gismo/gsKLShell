/** @file gsMaterialMatrixBaseDim_test.cpp

    @brief Analytic oracle for the contravariant members published by
           gsMaterialMatrixBaseDim::_getMetric (task 45).

    ### Why this suite exists
    _getMetricDeformed_impl / _getMetricUndeformed_impl do not COMPUTE the
    metric: the heavy lifting happens once per element in
    _computeMetricUndeformed / _computeMetricDeformed, which store Acov/Acon and
    acov/acon into the batched m_*_mat arrays (gsMaterialMatrixBaseDim.hpp:522,
    568, 640, 685 -- there Acon = Acov.inverse() and a^i = A^{ij} a_j, both
    correct). _getMetric* then merely COPIES column k of those arrays into the
    per-point scalar members. Four of those copies fed a CONTRAVARIANT member
    from its COVARIANT sibling:

      | dim | member      | was fed  | site (pre-fix line) |
      |-----|-------------|----------|---------------------|
      |  3  | m_acon_def  | acov_def | :1174               |
      |  2  | m_Acon_def  | Acov_def | :1203               |
      |  3  | m_acon_ori  | acov_ori | :1305               |
      |  2  | m_Acon_ori  | Acov_ori | :1332               |

    Three of the four are latent (m_acon_ori/m_acon_def/m_Acon_def are written
    but never read in-tree). m_Acon_ori is NOT latent: it is the SvK material
    tensor's only geometric input in gsMaterialMatrixLinear.hpp:554 and
    gsMaterialMatrixNonlinear.hpp:1304, and gsMaterialMatrixLinear<2,real_t> /
    gsMaterialMatrixNonlinear<2,real_t,...> are both instantiated. So the dim-2
    row of this table was a live wrong-stiffness bug, not a dormant one.

    ### The oracle
    Not a golden file and not a cross-implementation comparison: the two
    DEFINING identities of a contravariant basis / an inverse metric, which hold
    at every point of every non-degenerate patch by construction.

      A_{ab} A^{bc} = delta      ->  Acov * Acon == I_2
      a^i . a_j     = delta      ->  acon^T * acov == I_2

    Both are 2x2 in-plane statements; for dim==3 the basis vectors are 3x2
    (surface vectors in space) and the product is still 2x2.

    ### Why the fixtures must be curved / warped -- and why that makes the
    ### suite automatically non-vacuous
    On an un-warped BSplineSquare the Jacobian is the identity, so
    Acov == Acon == I and acov == acon: EVERY mis-assignment in the table above
    would satisfy both identities exactly. The fixtures are therefore the
    non-orthogonal doubly-curved bump (3D) and the smoothly-warped square (2D)
    of gsShellKinematics_test.cpp:188/210, and each test ASSERTS that its metric
    is neither the identity nor diagonal.

    That assertion is exactly the non-vacuity proof, because under the bug the
    two checks degenerate into precisely those quantities:
      acon := acov  =>  acon^T * acov == Acov , deviation ||Acov - I||
      Acon := Acov  =>  Acov  * Acon  == Acov^2, deviation ||Acov^2 - I||
    so "the fixture metric is not the identity" and "the poisoned check fails"
    are the same statement. A green run on a degenerate fixture is impossible.
    NB this must hold per CONFIGURATION: a poison at a _def site deviates by
    ||Acov_def - I||, which "Acov_ori != I" does NOT bound below. Both tests
    therefore assert the ori AND the def metric to be away from the identity.

    ### Coverage note
    Each dimension is checked for BOTH identities even though only one of them
    was broken per dimension (dim-3: the a-vectors; dim-2: the A-tensors). The
    already-correct half is kept as a regression guard on the sibling lines.

    ### The SvK 2D-vs-3D differential oracle (task 51)
    The two identity tests above certify the metric MEMBERS. They do not
    witness that a material quantity BUILT from them comes out right: the
    suite-wide poison experiment of task 50 broke :1332 and
    gsMaterialMatrixTFT_test::MM_SvK did not move, because every assertion in
    that suite is a self-consistency / finite-difference check of the material
    against itself -- a wrong C^{ijkl} shifts both sides and cancels.

    The last three tests close that gap with TWO oracles per quantity:

    (a) DIFFERENTIAL. _Cijkl is one dim-agnostic template body
        (gsMaterialMatrixLinear.hpp:552-554) whose ONLY geometric input is
        m_Acon_ori, and its Cconstant = 2*lambda*mu/(lambda+2*mu) is the same
        plane-stress-condensed constant for dim 2 and dim 3 (no condensation
        step, gsPlaneStressCondensation is not involved). So a dim-2 material
        on the warped planar patch and a dim-3 material on the SAME patch
        embedded in 3D (embed(3) zero-fills z, hence FLAT, hence B == 0 and
        its surface first fundamental form is the planar metric identically)
        must produce the SAME C^{ijkl} and the same stress, to machine
        precision. They travel through different dim dispatches to get there.

    (b) ANALYTIC. The shell-Voigt C and the stress are also recomputed inside
        the test from the patch Jacobian (gsFunction::jacobian -> deriv_into,
        a different evaluation path from the gsMapData::computeMap used by
        _computeMetric*) and the textbook plane-stress SvK formula. This half
        is needed because (a) is structurally BLIND to any defect in the
        shared _Cijkl body: both dimensions instantiate the same template, so
        a wrong Cconstant moves both sides equally. Scope of the independence
        claim: independent of the metric plumbing and of the dim dispatch --
        it never reads m_Acon_ori, it inverts a 2x2 metric of its own -- but
        NOT independent of a transcription error in the constitutive formula
        itself, which is copied from the model definition
        (gsMaterialMatrixLinear.h:265 doxygen).

    Non-vacuity is again structural: pre-fix, dim 2 wrote Acon_ori := Acov_ori
    while dim 3 wrote the true contravariant tensor, and on this warped patch
    Acov != Acon, so both oracles are violated. The warp is LOAD-BEARING; an
    un-warped BSplineSquare has Acov == Acon == I and every one of these
    checks would pass for any frame convention. That is what
    SvK_2D_vs_3D_fixture_is_warped guards.

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"

#include <gsKLShell/src/gsMaterialMatrixLinear.h>

SUITE(gsMaterialMatrixBaseDim)
{

// =====================================================================
// Oracle access
// =====================================================================

// gsMaterialMatrixBaseDim exposes _computePoints / _getMetric publicly, but the
// RESULTS live in the protected per-thread gsMaterialMatrixBaseDimData m_data
// (gsMaterialMatrixBaseDim.h:482). A minimal derived probe publishes them.
// gsMaterialMatrixLinear is used only because it is the cheapest CONCRETE
// BaseDim derivative with a (mp,thickness,E,nu) constructor -- no material
// response is evaluated here, so the parameter VALUES are irrelevant.
// Two non-template probes (rather than one templated over dim) keep m_data a
// non-dependent name, so it needs no this-> qualification.
// (Pattern shared with gsShellKinematics_test.cpp:64/80 -- adapted, not
// included, following this project's one-fixture-per-test-file convention.)
class BaseDimProbe3 : public gsMaterialMatrixLinear<3,real_t>
{
public:
    BaseDimProbe3(const gsFunctionSet<real_t> & mp,  const gsFunctionSet<real_t> & thick,
                  const gsFunctionSet<real_t> & E,   const gsFunctionSet<real_t> & nu)
    : gsMaterialMatrixLinear<3,real_t>(mp,thick,E,nu) { }

    // _getMetric is public in gsMaterialMatrixBaseDim (.h:263) and is INHERITED
    // publicly here: gsMaterialMatrixLinear used to re-declare it in a protected
    // block (`using Base::_getMetric;`), which lowered its access, but task 46
    // removed that re-declaration -- see the note now standing at
    // gsMaterialMatrixLinear.h:290-295. This using-declaration is therefore a
    // no-op today; it is kept so the probe does not depend on which of the two
    // arrangements gsMaterialMatrixLinear happens to carry.
    using gsMaterialMatrixBaseDim<3,real_t>::_getMetric;

    const gsMaterialMatrixBaseDimData<3,real_t> & data() const { return m_data.mine(); }
};

class BaseDimProbe2 : public gsMaterialMatrixLinear<2,real_t>
{
public:
    BaseDimProbe2(const gsFunctionSet<real_t> & mp,  const gsFunctionSet<real_t> & thick,
                  const gsFunctionSet<real_t> & E,   const gsFunctionSet<real_t> & nu)
    : gsMaterialMatrixLinear<2,real_t>(mp,thick,E,nu) { }

    using gsMaterialMatrixBaseDim<2,real_t>::_getMetric;   // see BaseDimProbe3

    const gsMaterialMatrixBaseDimData<2,real_t> & data() const { return m_data.mine(); }
};

// =====================================================================
// Shared constants and helpers
// =====================================================================

// The identities are exact in exact arithmetic; the only error is the rounding
// of one 2x2 inverse and one 2-term linear combination on an O(1)
// well-conditioned metric, i.e. a few ulp. 1e-12 is a generous ceiling, never a
// tuned number: the observed deviation is expected at the 1e-16 level and is
// PRINTED by every test so that silent drift toward the ceiling is visible.
const real_t TOL = 1e-12;

// Material constants (never used numerically -- see BaseDimProbe*).
const real_t E_MOD = 200.0;
const real_t NU    = 0.3;
const real_t THICK = 0.01;

// Relative Frobenius deviation, robust for near-zero references.
real_t relFro(const gsMatrix<real_t> & a, const gsMatrix<real_t> & ref)
{
    return (a - ref).norm() / ((real_t)1.0 + ref.norm());
}

// Shape gate. A product of wrongly-shaped operands asserts inside Eigen, and a
// 0x0 field would make the identity check compare nothing at all, so every
// operand is validated BEFORE it is multiplied.
bool shapeIs(const char * name, const gsMatrix<real_t> & M, index_t r, index_t c)
{
    CHECK_EQUAL(r, M.rows());
    CHECK_EQUAL(c, M.cols());
    if (M.rows()!=r || M.cols()!=c)
    {
        gsInfo<<"[gsMaterialMatrixBaseDim] SHAPE failure in "<<name<<" : got "
              <<M.rows()<<"x"<<M.cols()<<", expected "<<r<<"x"<<c<<"\n";
        return false;
    }
    return true;
}

// The actual assertion: @a M must be the 2x2 identity. @a worst accumulates the
// deviation and @a nCmp counts the checks actually performed, so a run that
// silently skipped everything is visible in the printed summary.
void checkIsIdentity2(const char * name, index_t k, real_t z,
                      const gsMatrix<real_t> & M, real_t & worst, index_t & nCmp)
{
    if (!shapeIs(name, M, 2, 2))
        return;                       // already reported; a value check would be noise

    gsMatrix<real_t> Id = gsMatrix<real_t>::Identity(2,2);
    const real_t d = (M - Id).norm();
    worst = math::max(worst, d);
    ++nCmp;

    if (d > TOL)
        gsInfo<<"[gsMaterialMatrixBaseDim] IDENTITY failure in "<<name
              <<" at k="<<k<<", z="<<z<<" : |M - I|_F = "<<d<<"\n  M = "<<M<<"\n";

    CHECK(d <= TOL);
}

// =====================================================================
// Fixtures (adapted from gsShellKinematics_test.cpp:188/210/231/249)
// =====================================================================

// Doubly-curved, NON-orthogonal, non-rational bicubic patch:
// z = 0.30 X^2 + 0.20 Y^2 + 0.15 XY over the unit square. a_1.a_2 != 0, so
// Acov is neither the identity nor diagonal and a covariant/contravariant flip
// cannot hide.
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
        c(i,2) = 0.30*X*X + 0.20*Y*Y + 0.15*X*Y;
    }
    GISMO_ENSURE(mp.nPatches()==1 && mp.targetDim()==3, "curved patch build failed.");
    return mp;
}

// Planar (2D) counterpart: the unit square, degree-2, refined, with a mild
// SMOOTH in-plane warp of the control net. The warp is what makes the test a
// real gate: an un-warped BSplineSquare has the identity Jacobian, so
// Acov == Acon == I and every frame convention would agree trivially. The
// amplitude is small enough to keep det J > 0 everywhere.
gsMultiPatch<real_t> makePlanar()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
    mp.patch(0).degreeElevate(1);
    mp.patch(0).uniformRefine();
    gsMatrix<real_t> & c = mp.patch(0).coefs();
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0), Y = c(i,1);
        c(i,0) += 0.10*X*Y + 0.05*Y;       // shear + stretch: non-orthogonal, non-constant J
        c(i,1) += 0.12*Y*Y - 0.04*X;
    }
    GISMO_ENSURE(mp.nPatches()==1 && mp.targetDim()==2, "planar patch build failed.");
    return mp;
}

// Deterministic, smooth, NON-homogeneous control-net displacement (no rand()):
// each control point moves by a low-order polynomial of its own NORMALIZED
// position, so the DEFORMED metric differs from the undeformed one everywhere
// (otherwise the _def half of every check would duplicate the _ori half).
// Works for targetDim 2 and 3.
gsMultiPatch<real_t> deformSmooth(const gsMultiPatch<real_t> & mp, real_t s)
{
    gsMultiPatch<real_t> def = mp;
    gsMatrix<real_t> & c = def.patch(0).coefs();
    const index_t d = c.cols();
    const real_t  L = math::max(c.cwiseAbs().maxCoeff(), (real_t)1.0);  // geometry scale
    for (index_t i=0; i!=c.rows(); ++i)
    {
        const real_t X = c(i,0)/L, Y = c(i,1)/L, Z = (d>2 ? c(i,2)/L : (real_t)0.0);
        c(i,0) += s * L * ( 0.20*X*X + 0.10*X*Y + 0.05*Y );
        c(i,1) += s * L * (-0.15*Y*Y + 0.08*X   + 0.03*Z );
        if (d>2)
            c(i,2) += s * L * ( 0.12*X*Y - 0.06*Y + 0.04*X );
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

// Constant FunctionExpr on a @a domainDim-dimensional domain: BaseDim's
// _computePoints evaluates thickness and parameters at the PHYSICAL map points
// (gsMaterialMatrixBaseDim.hpp:447-455), whose dimension is dim.
gsFunctionExpr<real_t> constFun(real_t v, short_t domainDim)
{ return gsFunctionExpr<real_t>(util::to_string(v), domainDim); }

// =====================================================================
// TEST 1 : dim == 3, curved surface -- gates :1174 (m_acon_def) and
//          :1305 (m_acon_ori) via the dual-basis identity
// =====================================================================
TEST(ContravariantIdentities_curved3D)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.10);

    gsFunctionExpr<real_t> t  = constFun(THICK, 3);
    gsFunctionExpr<real_t> Ef = constFun(E_MOD, 3);
    gsFunctionExpr<real_t> nf = constFun(NU,    3);

    BaseDimProbe3 probe(mp, t, Ef, nf);
    probe.setDeformed(&mp_def);

    const gsMatrix<real_t> uv = interiorPoints();
    probe._computePoints(0, uv);

    // The a/A members are z-INDEPENDENT (they are midsurface quantities copied
    // straight out of the batched arrays). Sweeping z anyway proves the copies
    // are not contaminated by the through-thickness terms that share the body.
    const std::vector<real_t> zs = {-0.4, 0.0, 0.35};

    real_t  worst = 0.0;
    index_t nCmp  = 0;

    for (index_t k=0; k!=uv.cols(); ++k)
        for (size_t j=0; j!=zs.size(); ++j)
        {
            const real_t z = zs[j];
            probe._getMetric(k,z);                  // fills probe.data().m_*
            const gsMaterialMatrixBaseDimData<3,real_t> & d = probe.data();

            // --- dual bases: a^i . a_j == delta  (THE gate for :1174/:1305) ---
            if (shapeIs("acon_ori", d.m_acon_ori, 3,2) && shapeIs("acov_ori", d.m_acov_ori, 3,2))
            {
                gsMatrix<real_t> P = d.m_acon_ori.transpose() * d.m_acov_ori;
                checkIsIdentity2("acon_ori^T * acov_ori", k, z, P, worst, nCmp);
            }
            if (shapeIs("acon_def", d.m_acon_def, 3,2) && shapeIs("acov_def", d.m_acov_def, 3,2))
            {
                gsMatrix<real_t> P = d.m_acon_def.transpose() * d.m_acov_def;
                checkIsIdentity2("acon_def^T * acov_def", k, z, P, worst, nCmp);
            }

            // --- inverse metric: A_ab A^bc == delta (regression guard: the
            //     dim-3 A-copies at :1152/:1283 were already correct) ---
            if (shapeIs("Acov_ori", d.m_Acov_ori, 2,2) && shapeIs("Acon_ori", d.m_Acon_ori, 2,2))
            {
                gsMatrix<real_t> P = d.m_Acov_ori * d.m_Acon_ori;
                checkIsIdentity2("Acov_ori * Acon_ori", k, z, P, worst, nCmp);
            }
            if (shapeIs("Acov_def", d.m_Acov_def, 2,2) && shapeIs("Acon_def", d.m_Acon_def, 2,2))
            {
                gsMatrix<real_t> P = d.m_Acov_def * d.m_Acon_def;
                checkIsIdentity2("Acov_def * Acon_def", k, z, P, worst, nCmp);
            }
        }

    CHECK_EQUAL(72, nCmp);      // 6 points x 3 heights x 4 identities
    gsInfo<<"[gsMaterialMatrixBaseDim] ContravariantIdentities_curved3D : "<<nCmp
          <<" identity checks, worst |M - I|_F = "<<worst<<"\n";
    CHECK(worst <= TOL);

    // --- Non-vacuity: see the file header. Under the pre-fix code the two
    //     _ori checks reduce to ||Acov_ori - I|| and the fixture makes that
    //     O(0.1), i.e. 11 orders of magnitude above TOL. ---
    probe._getMetric(0, 0.0);
    const gsMaterialMatrixBaseDimData<3,real_t> & d0 = probe.data();
    gsMatrix<real_t> Id2 = gsMatrix<real_t>::Identity(2,2);
    CHECK(relFro(d0.m_Acov_ori, Id2) > 1e-3);                 // metric != identity
    CHECK(math::abs(d0.m_Acov_ori(0,1)) > 1e-3);              // non-orthogonal: a_1.a_2 != 0
    CHECK(relFro(d0.m_Acov_def, d0.m_Acov_ori) > 1e-3);       // deformed != undeformed
    // The _def half needs its OWN lower bound: a poison at a _def site deviates
    // by ||Acov_def - I||, which neither of the two asserts above bounds below
    // (Acov_def could in principle sit closer to I than Acov_ori does).
    CHECK(relFro(d0.m_Acov_def, Id2) > 1e-3);                 // deformed metric != identity
    gsInfo<<"[gsMaterialMatrixBaseDim]   fixture: ||Acov_ori - I||rel = "
          <<relFro(d0.m_Acov_ori, Id2)<<", Acov_ori(0,1) = "<<d0.m_Acov_ori(0,1)<<"\n";
}

// =====================================================================
// TEST 2 : dim == 2, warped planar patch -- gates :1203 (m_Acon_def) and
//          :1332 (m_Acon_ori) via the inverse-metric identity.
//          m_Acon_ori is the SvK material tensor's only geometric input
//          (gsMaterialMatrixLinear.hpp:554), so THIS is the live one.
// =====================================================================
TEST(ContravariantIdentities_planar2D)
{
    gsMultiPatch<real_t> mp     = makePlanar();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.10);

    gsFunctionExpr<real_t> t  = constFun(THICK, 2);
    gsFunctionExpr<real_t> Ef = constFun(E_MOD, 2);
    gsFunctionExpr<real_t> nf = constFun(NU,    2);

    BaseDimProbe2 probe(mp, t, Ef, nf);
    probe.setDeformed(&mp_def);

    const gsMatrix<real_t> uv = interiorPoints();
    probe._computePoints(0, uv);

    // The dim==2 metric is z-independent in the implementation, so one height
    // suffices (matches gsShellKinematics_test.cpp:373).
    const real_t z = 0.0;

    real_t  worst = 0.0;
    index_t nCmp  = 0;

    for (index_t k=0; k!=uv.cols(); ++k)
    {
        probe._getMetric(k,z);
        const gsMaterialMatrixBaseDimData<2,real_t> & d = probe.data();

        // --- inverse metric (THE gate for :1203/:1332) ---
        if (shapeIs("Acov_ori", d.m_Acov_ori, 2,2) && shapeIs("Acon_ori", d.m_Acon_ori, 2,2))
        {
            gsMatrix<real_t> P = d.m_Acov_ori * d.m_Acon_ori;
            checkIsIdentity2("Acov_ori * Acon_ori", k, z, P, worst, nCmp);
        }
        if (shapeIs("Acov_def", d.m_Acov_def, 2,2) && shapeIs("Acon_def", d.m_Acon_def, 2,2))
        {
            gsMatrix<real_t> P = d.m_Acov_def * d.m_Acon_def;
            checkIsIdentity2("Acov_def * Acon_def", k, z, P, worst, nCmp);
        }

        // --- dual bases (regression guard: the dim-2 a-copies at :1235/:1364
        //     were already correct). Here the basis vectors are 2x2. ---
        if (shapeIs("acon_ori", d.m_acon_ori, 2,2) && shapeIs("acov_ori", d.m_acov_ori, 2,2))
        {
            gsMatrix<real_t> P = d.m_acon_ori.transpose() * d.m_acov_ori;
            checkIsIdentity2("acon_ori^T * acov_ori", k, z, P, worst, nCmp);
        }
        if (shapeIs("acon_def", d.m_acon_def, 2,2) && shapeIs("acov_def", d.m_acov_def, 2,2))
        {
            gsMatrix<real_t> P = d.m_acon_def.transpose() * d.m_acov_def;
            checkIsIdentity2("acon_def^T * acov_def", k, z, P, worst, nCmp);
        }
    }

    CHECK_EQUAL(24, nCmp);      // 6 points x 4 identities
    gsInfo<<"[gsMaterialMatrixBaseDim] ContravariantIdentities_planar2D : "<<nCmp
          <<" identity checks, worst |M - I|_F = "<<worst<<"\n";
    CHECK(worst <= TOL);

    // --- Non-vacuity: under the pre-fix code the A-checks reduce to
    //     ||Acov^2 - I||, and the warp guarantees Acov is far from I. ---
    probe._getMetric(0, z);
    const gsMaterialMatrixBaseDimData<2,real_t> & d0 = probe.data();
    gsMatrix<real_t> Id2 = gsMatrix<real_t>::Identity(2,2);
    CHECK(relFro(d0.m_Acov_ori, Id2) > 1e-3);                 // metric != identity
    CHECK(math::abs(d0.m_Acov_ori(0,1)) > 1e-3);              // non-orthogonal warp
    CHECK(relFro(d0.m_Acov_def, d0.m_Acov_ori) > 1e-3);       // deformed != undeformed
    CHECK(relFro(d0.m_Acov_def, Id2) > 1e-3);                 // deformed metric != identity (see TEST 1)
    gsInfo<<"[gsMaterialMatrixBaseDim]   fixture: ||Acov_ori - I||rel = "
          <<relFro(d0.m_Acov_ori, Id2)<<", Acov_ori(0,1) = "<<d0.m_Acov_ori(0,1)<<"\n";
}

// =====================================================================
// Task 51 : the SvK 2D-vs-3D differential oracle + its analytic twin.
//           See the file header for the design; below are only the tools.
// =====================================================================

// Machine-precision ceilings for the two comparisons. The differential one is
// expected to be EXACT (both paths execute the same arithmetic on the same
// numbers), the analytic one costs one extra 2x2 inverse and a handful of
// products on an O(1) well-conditioned metric. Both are PRINTED so drift
// toward the ceiling is visible.
// The values are the task's 1e-13 / 1e-12, floored at a precision-aware level
// so that a float or multiprecision real_t build stays meaningful instead of
// impossible. With real_t == double the floor is INACTIVE (100*eps = 2.2e-14,
// below 1e-13), i.e. these are literally 1e-13 and 1e-12 here.
const real_t TOL_MAT = math::max( (real_t)1e-13, (real_t)(100*math::limits::epsilon()) );
const real_t TOL_VEC = math::max( (real_t)1e-12, (real_t)(100*math::limits::epsilon()) );

// Covariant metric A_ab = a_a . a_b straight from the patch Jacobian. This is
// deliberately NOT the route _computeMetric* takes (that one goes through
// gsMapData::computeMap); gsFunction::jacobian dispatches to deriv_into.
// Returns 2x2 for a 2- and for a 3-column control net alike (the third row of
// J is identically zero for the flat twin, contributing an exact 0.0).
gsMatrix<real_t> covMetric(const gsMultiPatch<real_t> & mp, const gsMatrix<real_t> & pt)
{
    gsMatrix<real_t> J = mp.patch(0).jacobian(pt);      // targetDim x 2
    gsMatrix<real_t> A = J.transpose() * J;
    return A;
}

// Plane-stress Saint-Venant-Kirchhoff shell tensor in the 3x3 shell-Voigt
// layout (11,22,12) used by eval3D_matrix:
//      C = [C1111 C1122 C1112 ; . C2222 C2212 ; . . C1212]
// from the MODEL definition C^{ijkl} = 2 lambda mu/(lambda+2mu) a^ij a^kl
// + mu (a^ik a^jl + a^il a^jk), with a^ij the CONTRAVARIANT metric.
gsMatrix<real_t> svkVoigt(const gsMatrix<real_t> & Acov, real_t Emod, real_t nu)
{
    const gsMatrix<real_t> A = Acov.inverse();          // a^ij : the whole point
    const real_t mu     = Emod / (2.*(1.+nu));
    const real_t lambda = Emod*nu / ((1.+nu)*(1.-2.*nu));
    const real_t Cc     = 2.*lambda*mu/(lambda+2.*mu);  // already condensed

    // deliberately a closure over the formula, so the six Voigt slots below
    // read exactly like eval3D_matrix's own six lines
    auto Cijkl = [&A,&mu,&Cc](index_t i, index_t j, index_t k, index_t l)
                 { return Cc*A(i,j)*A(k,l) + mu*(A(i,k)*A(j,l) + A(i,l)*A(j,k)); };

    gsMatrix<real_t> C(3,3);
    C(0,0)          = Cijkl(0,0,0,0);
    C(1,1)          = Cijkl(1,1,1,1);
    C(2,2)          = Cijkl(0,1,0,1);
    C(1,0) = C(0,1) = Cijkl(0,0,1,1);
    C(2,0) = C(0,2) = Cijkl(0,0,0,1);
    C(2,1) = C(1,2) = Cijkl(1,1,0,1);
    return C;
}

// Membrane 2nd Piola-Kirchhoff stress S^ab = C^{abkl} E_kl with the
// Green-Lagrange membrane strain E = 1/2 (A_def - A_ori) -- i.e. exactly what
// _E(z,VectorN) builds (gsMaterialMatrixLinear.hpp:606). Layout (11,22,12).
gsVector<real_t> svkStress(const gsMatrix<real_t> & Acov_ori, const gsMatrix<real_t> & Acov_def,
                           real_t Emod, real_t nu)
{
    const gsMatrix<real_t> A = Acov_ori.inverse();
    const real_t mu     = Emod / (2.*(1.+nu));
    const real_t lambda = Emod*nu / ((1.+nu)*(1.-2.*nu));
    const real_t Cc     = 2.*lambda*mu/(lambda+2.*mu);
    const gsMatrix<real_t> E = 0.5*(Acov_def - Acov_ori);

    auto Cijkl = [&A,&mu,&Cc](index_t i, index_t j, index_t k, index_t l)
                 { return Cc*A(i,j)*A(k,l) + mu*(A(i,k)*A(j,l) + A(i,l)*A(j,k)); };
    auto Sij   = [&Cijkl,&E](index_t i, index_t j)
                 {
                     real_t s = 0.0;
                     for (index_t k=0; k!=2; ++k)
                         for (index_t l=0; l!=2; ++l)
                             s += Cijkl(i,j,k,l)*E(k,l);
                     return s;
                 };

    gsVector<real_t> S(3);
    S << Sij(0,0), Sij(1,1), Sij(0,1);
    return S;
}

// The 2D fixture and its FLAT 3D twin, built so that the twin is the same
// control net with an exactly-zero third coordinate. NB the deformation must
// be applied in 2D and embedded afterwards: deformSmooth() on an already-3D
// net adds an out-of-plane term (its `if (d>2)` branch), which would lift the
// twin off the plane and make the comparison undefined.
void makeTwins(gsMultiPatch<real_t> & mp2, gsMultiPatch<real_t> & def2,
               gsMultiPatch<real_t> & mp3, gsMultiPatch<real_t> & def3, real_t s)
{
    mp2  = makePlanar();
    def2 = deformSmooth(mp2, s);
    mp3  = mp2;   mp3.embed(3);
    def3 = def2;  def3.embed(3);
}

// Premise guard for the twin construction: flatness (exact) and in-plane
// identity of the control nets (exact). Everything below is meaningless
// without these two.
void checkTwinIsFlat(const gsMultiPatch<real_t> & mp2, const gsMultiPatch<real_t> & mp3)
{
    const gsMatrix<real_t> & c2 = mp2.patch(0).coefs();
    const gsMatrix<real_t> & c3 = mp3.patch(0).coefs();
    CHECK_EQUAL(2, c2.cols());
    CHECK_EQUAL(3, c3.cols());
    CHECK_EQUAL(c2.rows(), c3.rows());
    if (c3.cols()==3 && c2.cols()==2 && c2.rows()==c3.rows())
    {
        CHECK( c3.col(2).cwiseAbs().maxCoeff() == 0.0 );                     // z == 0 exactly
        CHECK( (c3.leftCols(2)-c2).cwiseAbs().maxCoeff() == 0.0 );           // same in-plane net
    }
}

// =====================================================================
// TEST 3 : premise guard. Without this, TESTs 4-5 are silently vacuous
//          the day someone flattens the warp in makePlanar().
// =====================================================================
TEST(SvK_2D_vs_3D_fixture_is_warped)
{
    gsMultiPatch<real_t> mp2, def2, mp3, def3;
    makeTwins(mp2, def2, mp3, def3, 0.10);
    checkTwinIsFlat(mp2, mp3);
    checkTwinIsFlat(def2, def3);

    gsFunctionExpr<real_t> t2 = constFun(THICK,2), E2 = constFun(E_MOD,2), n2 = constFun(NU,2);
    gsFunctionExpr<real_t> t3 = constFun(THICK,3), E3 = constFun(E_MOD,3), n3 = constFun(NU,3);
    BaseDimProbe2 p2(mp2,t2,E2,n2);   p2.setDeformed(&def2);
    BaseDimProbe3 p3(mp3,t3,E3,n3);   p3.setDeformed(&def3);

    const gsMatrix<real_t> uv = interiorPoints();
    p2._computePoints(0, uv);
    p3._computePoints(0, uv);

    const gsMatrix<real_t> Id2 = gsMatrix<real_t>::Identity(2,2);
    real_t  worstFlat = 0.0, minAwayFromId = math::limits::max(), minShear = math::limits::max(), worstTwin = 0.0;
    index_t nCmp = 0;

    for (index_t k=0; k!=uv.cols(); ++k)
    {
        p2._getMetric(k, 0.0);
        p3._getMetric(k, 0.0);
        const gsMaterialMatrixBaseDimData<2,real_t> & d2 = p2.data();
        const gsMaterialMatrixBaseDimData<3,real_t> & d3 = p3.data();

        // (i) the fixture is genuinely warped: Acov_ori is neither I nor diagonal
        const real_t awayId = relFro(d2.m_Acov_ori, Id2);
        const real_t shear  = math::abs(d2.m_Acov_ori(0,1));
        minAwayFromId = math::min(minAwayFromId, awayId);
        minShear      = math::min(minShear, shear);
        CHECK(awayId > 1e-3);
        CHECK(shear  > 1e-3);

        // (ii) the DEFORMED metric is warped too, and differs from the
        //      undeformed one -- the stress test needs both
        CHECK(relFro(d2.m_Acov_def, Id2) > 1e-3);
        CHECK(relFro(d2.m_Acov_def, d2.m_Acov_ori) > 1e-3);

        // (iii) the twin really is the same surface: its first fundamental
        //       form must equal the planar metric (this is what makes TESTs
        //       4-5 a comparison of MATERIAL code rather than of geometries)
        worstTwin = math::max(worstTwin, relFro(d3.m_Acov_ori, d2.m_Acov_ori));
        worstTwin = math::max(worstTwin, relFro(d3.m_Acov_def, d2.m_Acov_def));
        CHECK(relFro(d3.m_Acov_ori, d2.m_Acov_ori) <= TOL_MAT);
        CHECK(relFro(d3.m_Acov_def, d2.m_Acov_def) <= TOL_MAT);

        // (iv) flatness of the twin, in the metric members: B == 0 exactly,
        //      so Gcov_ori == Acov_ori at ANY height and the dim-3 path has
        //      no bending contamination to contribute
        worstFlat = math::max(worstFlat, d3.m_Bcov_ori.cwiseAbs().maxCoeff());
        worstFlat = math::max(worstFlat, d3.m_Bcov_def.cwiseAbs().maxCoeff());
        CHECK(worstFlat == 0.0);
        ++nCmp;
    }

    CHECK_EQUAL(6, nCmp);
    gsInfo<<"[gsMaterialMatrixBaseDim] SvK_2D_vs_3D_fixture_is_warped : "<<nCmp
          <<" points; min ||Acov_ori - I||rel = "<<minAwayFromId
          <<", min |Acov_ori(0,1)| = "<<minShear
          <<"; worst twin metric deviation = "<<worstTwin
          <<", worst |Bcov(3D)| = "<<worstFlat<<"\n";
}

// =====================================================================
// TEST 4 : the material tensor. Differential (dim 2 vs dim 3) AND analytic.
//          THE gate for :1332 on an evaluated material quantity.
// =====================================================================
TEST(SvK_2D_vs_3D_material_tensor)
{
    gsMultiPatch<real_t> mp2, def2, mp3, def3;
    makeTwins(mp2, def2, mp3, def3, 0.10);
    checkTwinIsFlat(mp2, mp3);
    checkTwinIsFlat(def2, def3);

    gsFunctionExpr<real_t> t2 = constFun(THICK,2), E2 = constFun(E_MOD,2), n2 = constFun(NU,2);
    gsFunctionExpr<real_t> t3 = constFun(THICK,3), E3 = constFun(E_MOD,3), n3 = constFun(NU,3);
    gsMaterialMatrixLinear<2,real_t> mm2(mp2,t2,E2,n2);   mm2.setDeformed(&def2);
    gsMaterialMatrixLinear<3,real_t> mm3(mp3,t3,E3,n3);   mm3.setDeformed(&def3);

    const gsMatrix<real_t> uv = interiorPoints();

    // z is (nz x npts) -- eval3D_matrix loops j < z.rows() and reads z(j,k).
    // The heights are MULTIPLIED by the thickness inside (THICK=0.01), so
    // z=0.3 is the effective height 0.003. Two heights are used to
    // demonstrate that C is z-independent on both paths (C^{ijkl} sees only
    // the midsurface Acon_ori; the flat twin additionally has B == 0).
    gsMatrix<real_t> z(2, uv.cols());
    z.row(0).setZero();
    z.row(1).setConstant(0.3);

    gsMatrix<real_t> C2 = mm2.eval3D_matrix(0,uv,z,MaterialOutput::Generic);
    gsMatrix<real_t> C3 = mm3.eval3D_matrix(0,uv,z,MaterialOutput::Generic);
    CHECK_EQUAL(9, C2.rows());   CHECK_EQUAL(uv.cols()*z.rows(), C2.cols());
    CHECK_EQUAL(9, C3.rows());   CHECK_EQUAL(uv.cols()*z.rows(), C3.cols());

    real_t  worstDiff = 0.0, worstRef2 = 0.0, worstRef3 = 0.0, minNorm = math::limits::max();
    index_t nCmp = 0;

    for (index_t k=0; k!=uv.cols(); ++k)
    {
        gsMatrix<real_t> pt = uv.col(k);
        // analytic reference, from the 2D patch Jacobian only
        gsMatrix<real_t> Cref = svkVoigt(covMetric(mp2,pt), E_MOD, NU);
        minNorm = math::min(minNorm, Cref.norm());

        for (index_t j=0; j!=z.rows(); ++j)
        {
            gsMatrix<real_t> Ck2 = C2.reshapeCol(j*uv.cols()+k,3,3);
            gsMatrix<real_t> Ck3 = C3.reshapeCol(j*uv.cols()+k,3,3);
            if (!shapeIs("C(2D)",Ck2,3,3) || !shapeIs("C(3D)",Ck3,3,3))
                continue;

            // symmetry of the shell-Voigt block (cheap sanity on the layout)
            CHECK(relFro(Ck2, Ck2.transpose()) <= TOL_MAT);

            const real_t dDiff = relFro(Ck2, Ck3);      // (a) differential
            const real_t dRef2 = relFro(Ck2, Cref);     // (b) analytic, dim 2
            const real_t dRef3 = relFro(Ck3, Cref);     // (b) analytic, dim 3
            worstDiff = math::max(worstDiff, dDiff);
            worstRef2 = math::max(worstRef2, dRef2);
            worstRef3 = math::max(worstRef3, dRef3);
            ++nCmp;

            if (dDiff > TOL_MAT || dRef2 > TOL_MAT || dRef3 > TOL_MAT)
                gsInfo<<"[gsMaterialMatrixBaseDim] SvK C failure at k="<<k<<", z="<<z(j,k)
                      <<" : rel(2D,3D) = "<<dDiff<<", rel(2D,ref) = "<<dRef2
                      <<", rel(3D,ref) = "<<dRef3
                      <<"\n  C(2D) = "<<Ck2<<"\n  C(3D) = "<<Ck3<<"\n  Cref  = "<<Cref<<"\n";

            CHECK(dDiff <= TOL_MAT);
            CHECK(dRef2 <= TOL_MAT);
            CHECK(dRef3 <= TOL_MAT);
        }
    }

    CHECK_EQUAL(12, nCmp);      // 6 points x 2 heights
    CHECK(minNorm > 1.0);       // the tensor is not the zero matrix
    gsInfo<<"[gsMaterialMatrixBaseDim] SvK_2D_vs_3D_material_tensor : "<<nCmp
          <<" comparisons, worst rel dev: 2D-vs-3D = "<<worstDiff
          <<", 2D-vs-analytic = "<<worstRef2
          <<", 3D-vs-analytic = "<<worstRef3
          <<" (tol "<<TOL_MAT<<", ||Cref|| >= "<<minNorm<<")\n";
}

// =====================================================================
// TEST 5 : the end-to-end half -- stress, i.e. C^{ijkl} contracted with the
//          strain built by _E. Same two oracles.
// =====================================================================
TEST(SvK_2D_vs_3D_stress)
{
    gsMultiPatch<real_t> mp2, def2, mp3, def3;
    makeTwins(mp2, def2, mp3, def3, 0.10);
    checkTwinIsFlat(mp2, mp3);
    checkTwinIsFlat(def2, def3);

    gsFunctionExpr<real_t> t2 = constFun(THICK,2), E2 = constFun(E_MOD,2), n2 = constFun(NU,2);
    gsFunctionExpr<real_t> t3 = constFun(THICK,3), E3 = constFun(E_MOD,3), n3 = constFun(NU,3);
    gsMaterialMatrixLinear<2,real_t> mm2(mp2,t2,E2,n2);   mm2.setDeformed(&def2);
    gsMaterialMatrixLinear<3,real_t> mm3(mp3,t3,E3,n3);   mm3.setDeformed(&def3);

    const gsMatrix<real_t> uv = interiorPoints();
    gsMatrix<real_t> z(1, uv.cols());
    z.setZero();

    // MEMBRANE branch: _E(z,VectorN) = 1/2 (Acov_def - Acov_ori). On this flat
    // fixture the Generic branch (which goes through Gcov) is algebraically the
    // same thing, and that equality is asserted below rather than assumed.
    gsMatrix<real_t> S2 = mm2.eval3D_vector(0,uv,z,MaterialOutput::VectorN);
    gsMatrix<real_t> S3 = mm3.eval3D_vector(0,uv,z,MaterialOutput::VectorN);
    gsMatrix<real_t> G2 = mm2.eval3D_vector(0,uv,z,MaterialOutput::Generic);
    gsMatrix<real_t> G3 = mm3.eval3D_vector(0,uv,z,MaterialOutput::Generic);
    CHECK_EQUAL(3, S2.rows());   CHECK_EQUAL(uv.cols(), S2.cols());
    CHECK_EQUAL(3, S3.rows());   CHECK_EQUAL(uv.cols(), S3.cols());

    real_t  worstDiff = 0.0, worstRef2 = 0.0, worstRef3 = 0.0, worstGen = 0.0;
    real_t  minShearFrac = math::limits::max(), minNorm = math::limits::max();
    index_t nCmp = 0;

    for (index_t k=0; k!=uv.cols(); ++k)
    {
        gsMatrix<real_t> pt   = uv.col(k);
        gsMatrix<real_t> Aori = covMetric(mp2,pt);
        gsMatrix<real_t> Adef = covMetric(def2,pt);
        gsVector<real_t> Sref = svkStress(Aori, Adef, E_MOD, NU);
        gsMatrix<real_t> Eref = 0.5*(Adef - Aori);

        // the deformation must actually shear: E_12 is a real fraction of the
        // strain, otherwise the C1112 / C2212 slots ride along untested
        const real_t shearFrac = math::abs(Eref(0,1)) / Eref.norm();
        minShearFrac = math::min(minShearFrac, shearFrac);
        CHECK(shearFrac > 1e-2);

        gsMatrix<real_t> Sk2 = S2.col(k), Sk3 = S3.col(k), Sk = Sref;
        minNorm = math::min(minNorm, Sref.norm());

        const real_t dDiff = relFro(Sk2, Sk3);
        const real_t dRef2 = relFro(Sk2, Sk);
        const real_t dRef3 = relFro(Sk3, Sk);
        worstDiff = math::max(worstDiff, dDiff);
        worstRef2 = math::max(worstRef2, dRef2);
        worstRef3 = math::max(worstRef3, dRef3);
        worstGen  = math::max(worstGen , relFro(G2.col(k), Sk2));
        worstGen  = math::max(worstGen , relFro(G3.col(k), Sk3));
        ++nCmp;

        if (dDiff > TOL_VEC || dRef2 > TOL_VEC || dRef3 > TOL_VEC)
            gsInfo<<"[gsMaterialMatrixBaseDim] SvK S failure at k="<<k
                  <<" : rel(2D,3D) = "<<dDiff<<", rel(2D,ref) = "<<dRef2
                  <<", rel(3D,ref) = "<<dRef3
                  <<"\n  S(2D) = "<<Sk2.transpose()<<"\n  S(3D) = "<<Sk3.transpose()
                  <<"\n  Sref  = "<<Sk.transpose()<<"\n";

        CHECK(dDiff <= TOL_VEC);
        CHECK(dRef2 <= TOL_VEC);
        CHECK(dRef3 <= TOL_VEC);
        // Generic == VectorN on a flat fixture (documented, not assumed)
        CHECK(relFro(G2.col(k), Sk2) <= TOL_VEC);
        CHECK(relFro(G3.col(k), Sk3) <= TOL_VEC);
    }

    CHECK_EQUAL(6, nCmp);
    CHECK(minNorm > 1e-3);      // a non-trivial stress state
    gsInfo<<"[gsMaterialMatrixBaseDim] SvK_2D_vs_3D_stress : "<<nCmp
          <<" comparisons, worst rel dev: 2D-vs-3D = "<<worstDiff
          <<", 2D-vs-analytic = "<<worstRef2
          <<", 3D-vs-analytic = "<<worstRef3
          <<", Generic-vs-VectorN = "<<worstGen
          <<" (tol "<<TOL_VEC<<"); min |E_12|/||E|| = "<<minShearFrac
          <<", min ||Sref|| = "<<minNorm<<"\n";
}

} // SUITE
