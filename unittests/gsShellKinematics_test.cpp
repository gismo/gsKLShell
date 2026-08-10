/** @file gsShellKinematics_test.cpp

    @brief Parity gate: gsShellKinematics<dim,T> vs the legacy gsMaterialMatrixBaseDim
           metric engine (task 21).

    gsShellKinematics (task 20) is a stateless re-implementation of the
    fundamental-form / through-thickness metric machinery of
    gsMaterialMatrixBaseDim, driven by an INJECTED gsMapData instead of by owned
    geometry. Its entire value is FIDELITY: every consumer of the new material
    provider inherits whatever this class computes. So the oracle here is the
    legacy engine itself, and the reference is NOT a tolerance:

      *** The headline assertion is BIT-FOR-BIT equality, (A-B).norm() == 0. ***

    That is legitimate (not a lucky-rounding accident) because:
      - the transplanted bodies are statement-for-statement identical expression
        trees over the same dynamic gsMatrix<T> types (never a fixed-size
        gsMatrix<T,3,3>, so Eigen dispatches the same inverse/determinant kernels);
      - both sides are fed a gsMapData computed on the same geometry with the SAME
        flag word as gsMaterialMatrixBaseDim.hpp:471/591 uses, so the inputs are
        the identical doubles;
      - gsKLShell is compiled with -O2 -std=c++14, no -ffast-math and no -march
        override, so there is no FMA contraction or excess-precision freedom that
        could perturb the last bit between the two translation units.
    A relative-Frobenius check against REL_TOL is kept as a clearly-labelled
    SECONDARY fallback, and the worst relative deviation of each test is printed;
    any nonzero deviation is reported as a genuine finding rather than tolerated.

    Fixtures. The undeformed patch is the NON-ORTHOGONAL doubly-curved bump
    z = 0.30X^2 + 0.20Y^2 + 0.15XY of gsMaterialMatrix3D_test.cpp:85 (adapted, not
    included across test files). That is the fixture class that exposed the task-15
    frame bug: it has B_ori != 0 (so the z-sweep is live) and a_1.a_2 != 0 (so a
    covariant/contravariant frame flip cannot hide). The fixture properties are
    ASSERTED in the tests, so an exact comparison of degenerate data cannot pass
    silently.

    Compared fields (the consumer contract of PointMetrics): Gcov_ori, Gcov_def,
    gcov_ori, gcon_ori and the scalar J0_sq. gcov_def/gcon_def are deliberately not
    computed by gsShellKinematics (outside the contract) and are not tested.

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"

#include <gsKLShell/src/gsShellKinematics.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>

SUITE(gsShellKinematics)
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

// Secondary (fallback) tolerance only; the primary assertion is exact equality,
// so this stays at the spec's machine-precision level for every real_t: if the
// transplant is faithful the deviation is identically 0 and the check is met
// regardless of the precision of real_t.
const real_t REL_TOL = 1e-14;

// Material constants (never used numerically -- see BaseDimProbe*).
const real_t E_MOD = 200.0;
const real_t NU    = 0.3;
const real_t THICK = 0.01;

// Relative Frobenius deviation, robust for near-zero references.
real_t relFro(const gsMatrix<real_t> & a, const gsMatrix<real_t> & ref)
{
    return (a - ref).norm() / ((real_t)1.0 + ref.norm());
}

// BaseDim's flag word, copied VERBATIM from gsMaterialMatrixBaseDim.hpp:471
// (== :591 for the undeformed twin). Paraphrasing it is the one way this test
// could silently lose bit-exactness, so it lives in exactly one place.
// NB NEED_JACOBIAN == NEED_DERIV; computeMap only ORs flags in.
unsigned baseDimFlags()
{ return NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2; }

// Hand-build the injected map exactly as BaseDim builds its internal one.
void computeMapLike(const gsGeometry<real_t> & g, const gsMatrix<real_t> & u,
                    gsMapData<real_t> & md)
{
    md.flags  = baseDimFlags();
    md.points = u;
    static_cast<const gsFunction<real_t>&>(g).computeMap(md);
}

// Exact-equality parity check for one 3x3 field, with the SHAPE checked FIRST:
// (A-B).norm() of two 0x0 matrices is 0, i.e. an unsized field would produce a
// green test that compared nothing, and mismatched shapes would assert inside
// Eigen. @a worst accumulates the relative deviation, @a nCmp counts the
// comparisons actually performed (printed, so a vacuous run is visible).
void checkFieldExact(const char * name, index_t k, real_t z,
                     const gsMatrix<real_t> & got, const gsMatrix<real_t> & ref,
                     real_t & worst, index_t & nCmp)
{
    CHECK_EQUAL(3, ref.rows());  CHECK_EQUAL(3, ref.cols());
    CHECK_EQUAL(3, got.rows());  CHECK_EQUAL(3, got.cols());
    if (ref.rows()!=3 || ref.cols()!=3 || got.rows()!=3 || got.cols()!=3)
    {
        gsInfo<<"[gsShellKinematics] SHAPE failure in "<<name<<" at k="<<k<<", z="<<z
              <<" : got "<<got.rows()<<"x"<<got.cols()
              <<", ref "<<ref.rows()<<"x"<<ref.cols()<<"\n";
        return; // shape already reported; a difference would be meaningless
    }

    const real_t d   = (got - ref).norm();
    const real_t rel = d / ((real_t)1.0 + ref.norm());
    worst = math::max(worst, rel);
    ++nCmp;

    if (d != (real_t)0)
        gsInfo<<"[gsShellKinematics] DEVIATION in "<<name<<" at k="<<k<<", z="<<z
              <<" : |diff|_F = "<<d<<" (relative "<<rel<<")\n"
              <<"  got = "<<got<<"\n  ref = "<<ref<<"\n";

    CHECK(d == (real_t)0);      // PRIMARY  : bit-for-bit
    CHECK(rel <= REL_TOL);      // SECONDARY: labelled fallback
}

// Same, for the scalar J0_sq.
void checkScalarExact(const char * name, index_t k, real_t z,
                      real_t got, real_t ref, real_t & worst, index_t & nCmp)
{
    const real_t d   = math::abs(got - ref);
    const real_t rel = d / ((real_t)1.0 + math::abs(ref));
    worst = math::max(worst, rel);
    ++nCmp;

    if (d != (real_t)0)
        gsInfo<<"[gsShellKinematics] DEVIATION in "<<name<<" at k="<<k<<", z="<<z
              <<" : got "<<got<<", ref "<<ref<<" (relative "<<rel<<")\n";

    CHECK(d == (real_t)0);      // PRIMARY  : bit-for-bit
    CHECK(rel <= REL_TOL);      // SECONDARY: labelled fallback
}

// =====================================================================
// Fixtures (adapted from gsMaterialMatrix3D_test.cpp:85/111/141)
// =====================================================================

// Doubly-curved, NON-orthogonal, non-rational bicubic patch:
// z = 0.30 X^2 + 0.20 Y^2 + 0.15 XY over the unit square. Genuine double
// curvature (B_ori != 0, so the z-sweep is live) and a_1.a_2 != 0 (so a frame
// flip cannot hide), while |g|~1 keeps the metric well-conditioned. Both
// properties are asserted in TEST(Curved3D_parity).
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
// position, so the strain field varies across the patch while the deformed map
// stays a valid spline. Works for targetDim 2 and 3.
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
// TEST 1 : curved, non-orthogonal surface in 3D -- the real gate
// =====================================================================
TEST(Curved3D_parity)
{
    gsMultiPatch<real_t> mp     = makeCurved();
    gsMultiPatch<real_t> mp_def = deformSmooth(mp, 0.10);

    gsFunctionExpr<real_t> t  = constFun(THICK, 3);
    gsFunctionExpr<real_t> Ef = constFun(E_MOD, 3);
    gsFunctionExpr<real_t> nf = constFun(NU,    3);

    // --- Oracle: the legacy engine on its own geometry ---
    BaseDimProbe3 probe(mp, t, Ef, nf);
    probe.setDeformed(&mp_def);

    const gsMatrix<real_t> uv = interiorPoints();
    probe._computePoints(0, uv);

    // --- Class under test: same geometry, injected maps ---
    gsMapData<real_t> mdOri, mdDef;
    computeMapLike(mp.patch(0),     uv, mdOri);
    computeMapLike(mp_def.patch(0), uv, mdDef);

    gsShellKinematics<3,real_t> kin;
    kin.computeUndeformed(mdOri);
    kin.computeDeformed  (mdDef);
    CHECK_EQUAL(uv.cols(), kin.nPoints());

    // PHYSICAL through-thickness heights, small against the geometry scale.
    const std::vector<real_t> zs = {-0.4, 0.0, 0.35};

    real_t  worst = 0.0;
    index_t nCmp  = 0;
    typename gsShellKinematics<3,real_t>::PointMetrics pm;

    for (index_t k=0; k!=uv.cols(); ++k)
        for (size_t j=0; j!=zs.size(); ++j)
        {
            const real_t z = zs[j];

            probe._getMetric(k,z);        // fills probe.data().m_*
            kin.getMetric(k,z,pm);

            // Compared in diagnostic order: Gcov_ori at z=0 is J^T J plus
            // (2,2)=1 -- no z-term, no inverse, no normal. If THAT one deviates,
            // the injected map (flags/points/geometry) is the cause, not the
            // transplanted formulas.
            checkFieldExact ("Gcov_ori", k, z, pm.Gcov_ori, probe.data().m_Gcov_ori, worst, nCmp);
            checkFieldExact ("Gcov_def", k, z, pm.Gcov_def, probe.data().m_Gcov_def, worst, nCmp);
            checkFieldExact ("gcov_ori", k, z, pm.gcov_ori, probe.data().m_gcov_ori, worst, nCmp);
            checkFieldExact ("gcon_ori", k, z, pm.gcon_ori, probe.data().m_gcon_ori, worst, nCmp);
            checkScalarExact("J0_sq"   , k, z, pm.J0_sq   , probe.data().m_J0_sq   , worst, nCmp);
        }

    // 6 points x 3 heights x 5 fields
    CHECK_EQUAL(90, nCmp);
    gsInfo<<"[gsShellKinematics] Curved3D_parity : "<<nCmp
          <<" comparisons, worst relative deviation = "<<worst<<"\n";
    CHECK(worst <= REL_TOL);

    // --- The fixture must exercise what it claims (else the exact comparisons
    //     above could be a green comparison of degenerate data) ---
    typename gsShellKinematics<3,real_t>::PointMetrics p0, pz;
    kin.getMetric(0, 0.0 , p0);
    kin.getMetric(0, 0.35, pz);
    // (a) genuine double curvature: G_ori really depends on z (B_ori != 0)
    CHECK(relFro(pz.Gcov_ori, p0.Gcov_ori) > 1e-3);
    // (b) NON-orthogonal parameterisation: a_1 . a_2 != 0
    CHECK(math::abs(p0.Gcov_ori(0,1)) > 1e-3);
    // (c) genuine deformation: the deformed metric is not the undeformed one
    CHECK(relFro(p0.Gcov_def, p0.Gcov_ori) > 1e-3);
    // (d) J0_sq is a live, finite, positive quantity
    CHECK(p0.J0_sq > 0.0);
    CHECK(math::abs(p0.J0_sq - 1.0) > 1e-3);
}

// =====================================================================
// TEST 2 : planar 2D patch (the dim==2 branch: B=0, n=e3)
// =====================================================================
TEST(Planar2D_parity)
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

    gsMapData<real_t> mdOri, mdDef;
    computeMapLike(mp.patch(0),     uv, mdOri);
    computeMapLike(mp_def.patch(0), uv, mdDef);

    gsShellKinematics<2,real_t> kin;
    kin.computeUndeformed(mdOri);
    kin.computeDeformed  (mdDef);
    CHECK_EQUAL(uv.cols(), kin.nPoints());

    // The dim==2 metric is z-independent in BOTH implementations, so a single
    // height suffices (z=0).
    const real_t z = 0.0;

    real_t  worst = 0.0;
    index_t nCmp  = 0;
    typename gsShellKinematics<2,real_t>::PointMetrics pm;

    for (index_t k=0; k!=uv.cols(); ++k)
    {
        probe._getMetric(k,z);
        kin.getMetric(k,z,pm);

        checkFieldExact ("Gcov_ori", k, z, pm.Gcov_ori, probe.data().m_Gcov_ori, worst, nCmp);
        checkFieldExact ("Gcov_def", k, z, pm.Gcov_def, probe.data().m_Gcov_def, worst, nCmp);
        checkFieldExact ("gcov_ori", k, z, pm.gcov_ori, probe.data().m_gcov_ori, worst, nCmp);
        checkFieldExact ("gcon_ori", k, z, pm.gcon_ori, probe.data().m_gcon_ori, worst, nCmp);
        checkScalarExact("J0_sq"   , k, z, pm.J0_sq   , probe.data().m_J0_sq   , worst, nCmp);
    }

    CHECK_EQUAL(30, nCmp);      // 6 points x 5 fields
    gsInfo<<"[gsShellKinematics] Planar2D_parity : "<<nCmp
          <<" comparisons, worst relative deviation = "<<worst<<"\n";
    CHECK(worst <= REL_TOL);

    // Fixture sanity: the warp really is non-trivial (otherwise Acov == I and
    // every frame convention would agree by accident), and the deformation is real.
    typename gsShellKinematics<2,real_t>::PointMetrics p0;
    kin.getMetric(0, z, p0);
    CHECK(math::abs(p0.Gcov_ori(0,1)) > 1e-3);              // non-orthogonal
    CHECK(relFro(p0.Gcov_ori, gsMatrix<real_t>::Identity(3,3)) > 1e-3); // not identity
    CHECK(relFro(p0.Gcov_def, p0.Gcov_ori) > 1e-3);         // deformed != undeformed
    CHECK(p0.J0_sq > 0.0);
}

// =====================================================================
// TEST 3 : ori == def aliasing -- an internal consistency check
// =====================================================================
// NOTE: this test exercises gsShellKinematics ALONE (no oracle), so it is not
// parity evidence; it pins the aliasing invariant that the deformed and
// undeformed branches implement the SAME formula, which the provider relies on
// when it evaluates an undeformed configuration.
TEST(UndeformedOnly)
{
    gsMultiPatch<real_t> mp = makeCurved();

    const gsMatrix<real_t> uv = interiorPoints();
    gsMapData<real_t> md;
    computeMapLike(mp.patch(0), uv, md);

    gsShellKinematics<3,real_t> kin;
    kin.computeUndeformed(md);
    kin.computeDeformed  (md);          // SAME map for both configurations

    const std::vector<real_t> zs = {-0.4, 0.0, 0.35};

    real_t  worst = 0.0;
    index_t nCmp  = 0;
    typename gsShellKinematics<3,real_t>::PointMetrics pm;

    for (index_t k=0; k!=uv.cols(); ++k)
        for (size_t j=0; j!=zs.size(); ++j)
        {
            const real_t z = zs[j];
            kin.getMetric(k,z,pm);

            // Gcov_def and Gcov_ori are built by two SEPARATE _impl bodies from
            // two separate batched member sets: with identical input they must
            // agree bit-for-bit, and hence J0_sq must be exactly 1.
            checkFieldExact ("Gcov_def(alias)", k, z, pm.Gcov_def, pm.Gcov_ori, worst, nCmp);
            checkScalarExact("J0_sq(alias)"   , k, z, pm.J0_sq   , (real_t)1.0 , worst, nCmp);
        }

    CHECK_EQUAL(36, nCmp);      // 6 points x 3 heights x 2 quantities
    gsInfo<<"[gsShellKinematics] UndeformedOnly  : "<<nCmp
          <<" comparisons, worst relative deviation = "<<worst<<"\n";
    CHECK(worst <= REL_TOL);

    // Non-vacuity: the metric really varies with z and over the patch.
    typename gsShellKinematics<3,real_t>::PointMetrics p0, pz;
    kin.getMetric(0, 0.0 , p0);
    kin.getMetric(0, 0.35, pz);
    CHECK(relFro(pz.Gcov_ori, p0.Gcov_ori) > 1e-3);
}

} // SUITE
