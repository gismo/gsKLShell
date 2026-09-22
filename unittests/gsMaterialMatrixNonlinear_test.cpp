/** @file gsMaterialMatrixNonlinear_test.cpp

    @brief Provides unittests for the compressible plane-stress through-thickness
           Newton solve in gsMaterialMatrixNonlinear (_eval3D_Compressible_C33).

    The routine iterates a scalar Newton solve for the through-thickness squared
    stretch C33 such that S33(C33)=0. Its two overloads share the same loop body:
    the geometry overload (reached through eval3D_detF) and the Cmat overload
    (reached through eval3D_matrix_C, which takes the deformation tensor
    directly). On exhaustion of itmax=100 iterations the routine no longer
    throws: it writes a quiet NaN and returns normally, because it is evaluated
    inside gsExprAssembler's OpenMP region, where an exception cannot be caught
    above the parallel construct. Arms A and B exercise that failure path
    directly at material level; arm C confirms it is visible through
    gsThinShellAssembler::assembleMatrix as ThinShellAssemblerStatus::AssemblyError.

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_MATRIX_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Two further arms gate the compressible plane-stress through-thickness
    solve's per-station indexing and its degenerate-metric handling:
      - A1 is a true regression gate: it pins eval3D_pstretch's per-station
        packing (colIdx = j*u.cols()+k) against a one-station-at-a-time
        reference, and is observed to FAIL when that packed read is
        mis-indexed.
      - A3 is a contract/behaviour-preservation test, NOT a regression gate:
        it pins the observable outcome (NaN, no throw) of a degenerate
        deformed metric (m_J0_sq==0); the pre-hardening code already reaches
        the same NaN by a slower path, so this test passes identically with
        or without that hardening and must never be read as evidence it is
        exercised.
      - A search for a reachable inadmissible converged through-thickness
        iterate (accepting a non-positive C33) under NH_ext/Analytical, swept
        over eval3D_matrix_C across m_J0_sq spanning 1e-100 to 1e100, found no
        such state under the current backtracking guard, with or without the
        convergence break's positivity condition. No test for that case is
        written here, since a test cannot be observed to fail.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

using namespace gismo;

namespace {

// Compressible Neo-Hookean, analytical implementation, plane stress -- the
// configuration measured to reproduce the through-thickness Newton failure
// of the snapping-thesis run (E=78e6, nu=0.4, t=1e-3;
// optional/gsStructuralAnalysis/examples/snapping_example_shell.cpp:343-350).
// The parameter functions are members so their lifetime matches
// materialMatrix's: gsMaterialMatrixNonlinear takes its gsFunctionSet<T>
// arguments by const-reference and keeps that reference rather than copying.
struct CompressibleNH
{
    gsMultiPatch<real_t> mp;
    gsFunctionExpr<real_t> t;
    gsFunctionExpr<real_t> E;
    gsConstantFunction<real_t> nu;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    explicit CompressibleNH(real_t nuValue = 0.4)
    :
    t("1e-3",2),
    E("78e6",2),
    nu(nuValue,2)
    {
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) ); // unit square, degree 1
        mp.addAutoBoundaries();

        std::vector<gsFunctionSet<real_t>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsOptionList options;
        options.addInt("Material","Material model: (1): NH",1);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        options.addSwitch("Compressibility","Compressibility: (true): Compressible",true);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
    }
};

// Compressible Neo-Hookean surface (dim==3), analytical implementation --
// the fixture the A1 gate needs. Under dim==2 every through-thickness
// station shares the same in-plane metric (gsMaterialMatrixBaseDim.hpp
// _getGcov_def_impl<_dim==2> carries no z-term at all), so C33 is identical
// station to station and a packing defect that swaps which station's C33 is
// read is invisible: dim==2 cannot discriminate A1. Only dim==3 carries the
// -2z*Bcov + z^2*n(x)n curvature term that makes C33 vary through the
// thickness. t=0.5 is thick relative to the unit-square footprint so that,
// combined with a curved deformed configuration, the through-thickness
// spread in C33 clears the comparison tolerance (see C33_MultiStation_pstretch).
struct BentShellNH
{
    gsMultiPatch<real_t> mp;
    gsFunctionExpr<real_t> t;
    gsFunctionExpr<real_t> E;
    gsConstantFunction<real_t> nu;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    explicit BentShellNH(real_t nuValue = 0.4)
    :
    t("0.5",3),
    E("78e6",3),
    nu(nuValue,3)
    {
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) ); // unit square, degree 1
        mp.addAutoBoundaries();
        mp.embed(3);

        std::vector<gsFunctionSet<real_t>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsOptionList options;
        options.addInt("Material","Material model: (1): NH",1);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        options.addSwitch("Compressibility","Compressibility: (true): Compressible",true);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);
    }
};

} // namespace

SUITE(gsMaterialMatrixNonlinear_test)
{

// Undeformed state: m_J0_sq=1, the Newton initial guess c(2,2)=1 already
// satisfies S33=0, so the loop breaks at it==0 and detF = sqrt(1*1) = 1
// exactly (up to solver round-off). Exercises the geometry overload
// (_eval3D_Compressible_C33:1865) through eval3D_detF, but not its loop body.
TEST(C33_Undeformed_Control)
{
    CompressibleNH fix;
    gsMultiPatch<real_t> mp_def = fix.mp;
    fix.materialMatrix->setDeformed(&mp_def);

    gsVector<real_t> pt(2); pt.setConstant(0.25);
    gsMatrix<real_t> z(1,1); z.setZero();

    gsMatrix<real_t> detF = fix.materialMatrix->eval3D_detF(0,pt,z,MaterialOutput::Generic);
    CHECK_CLOSE(1.0, detF(0,0), 1e-10);
}

// A moderate biaxial stretch (uniform coefficient scaling by s=1.2, which
// keeps det(Gcov_def)=s^4>0 by construction) that actually iterates the
// Newton loop. Measured: detF=1.10747 (s=1.2, thickness t=1e-3 does not
// enter the plane-stress C33 solve), so
// c33 = detF^2/m_J0_sq = 1.10747^2/1.2^4 = 0.5916, consistent with the
// physical expectation that the sheet thins under in-plane stretch
// (c33 < 1). Asserted as a loose bound, not the calibrated ladder value.
TEST(C33_Stretched_Control)
{
    CompressibleNH fix;
    gsMultiPatch<real_t> mp_def = fix.mp;
    mp_def.patch(0).coefs() *= 1.2;
    fix.materialMatrix->setDeformed(&mp_def);

    gsVector<real_t> pt(2); pt.setConstant(0.25);
    gsMatrix<real_t> z(1,1); z.setZero();

    gsMatrix<real_t> detF = fix.materialMatrix->eval3D_detF(0,pt,z,MaterialOutput::Generic);
    CHECK(math::isfinite(detF(0,0)));
    CHECK(detF(0,0) > 0.0);
    CHECK(detF(0,0) < 1.2*1.2); // strictly less than the in-plane stretch factor squared
}

// Calibrated non-convergent state: isotropic collapse of the deformed
// coefficients by s=1e-1 (mildest rung of the search ladder that fails --
// s in {1e-1,...,1e-6} were tried, s=1e-1 already diverges). Uniform scaling
// keeps det(Gcov_def)=s^4>0, so this never touches the negative-determinant
// guard in gsMaterialMatrixBaseDim::_getMetric; the failure is purely the
// Newton solve's own exhaustion (m_J0_sq=1e-4 blows up the initial guess
// c(2,2)=m_J0_sq^-1=1e4, and the iteration never recovers within itmax=100).
// The measured warning trail (backtracking, then non-convergence) is quoted
// in the task report.
TEST(C33_NonConvergent_Geometry)
{
    const real_t s = 1e-1;

    CompressibleNH fix;
    gsMultiPatch<real_t> mp_def = fix.mp;
    mp_def.patch(0).coefs() *= s;
    fix.materialMatrix->setDeformed(&mp_def);

    gsVector<real_t> pt(2); pt.setConstant(0.25);
    gsMatrix<real_t> z(1,1); z.setZero();

    bool threw = false;
    gsMatrix<real_t> detF;
    try
    {
        detF = fix.materialMatrix->eval3D_detF(0,pt,z,MaterialOutput::Generic);
    }
    catch (...)
    {
        threw = true;
    }
    CHECK(!threw);
    // Guarded: on a (regression) throw, detF is still the default-constructed
    // 0x0 matrix, and detF(0,0) would read through a null gsMatrix::data().
    if (!threw)
        CHECK(!math::isfinite(detF(0,0)));
}

// Healthy control for the Cmat overload (_eval3D_Compressible_C33:1973,
// reached through eval3D_matrix_C). Cmat is the undeformed metric of the
// patch at pt, computed from the geometry map rather than hard-coded.
// Cross-checked against eval3D_matrix on the undeformed configuration --
// the two APIs must agree on the healthy path.
TEST(C33_Cmat_Healthy_Control)
{
    CompressibleNH fix;
    fix.materialMatrix->setDeformed(&fix.mp);

    gsVector<real_t> pt(2); pt.setConstant(0.25);
    gsMatrix<real_t> z(1,1); z.setZero();

    gsMapData<real_t> map;
    map.flags = NEED_JACOBIAN;
    map.points = pt;
    static_cast<const gsFunction<real_t>&>(fix.mp.patch(0)).computeMap(map);
    gsMatrix<real_t> g = map.jacobian(0);
    g.resize(2,2);
    gsMatrix<real_t> G = g.transpose() * g;

    gsMatrix<real_t> Cmat(3,1);
    Cmat(0,0) = G(0,0);
    Cmat(1,0) = G(1,1);
    Cmat(2,0) = G(0,1);

    gsMatrix<real_t> C_Cmat = fix.materialMatrix->eval3D_matrix_C(Cmat,0,pt,z(0,0),MaterialOutput::Generic);
    CHECK(C_Cmat.allFinite());

    gsMatrix<real_t> C_ref = fix.materialMatrix->eval3D_matrix(0,pt,z,MaterialOutput::Generic);
    gsMatrix<real_t> C_ref_3x3 = C_ref.reshape(3,3);
    gsMatrix<real_t> C_Cmat_3x3 = C_Cmat.reshape(3,3);
    CHECK_MATRIX_CLOSE(C_ref_3x3, C_Cmat_3x3, 1e-8);
}

// Calibrated non-convergent Cmat: near-degenerate shear, Cmat=(1,1,1-eps)
// with eps=1e-2 (mildest rung of {1e-2,1e-4,1e-6,1e-8} that fails --
// det(Gcov_def)=1-(1-eps)^2=0.0199 -> 0 while both diagonal entries stay 1).
// All nine entries of the measured 3x3 result were observed non-finite at
// this state (the warning trail and the NaN matrix are quoted in the task
// report), so the blanket assertion below reflects what was actually
// measured, not an assumption that _Cijkl must touch every slot.
TEST(C33_NonConvergent_Cmat)
{
    const real_t eps = 1e-2;

    CompressibleNH fix;
    fix.materialMatrix->setDeformed(&fix.mp);

    gsVector<real_t> pt(2); pt.setConstant(0.25);
    gsMatrix<real_t> z(1,1); z.setZero();

    gsMatrix<real_t> Cmat(3,1);
    Cmat<<1,1,1-eps;

    bool threw = false;
    gsMatrix<real_t> C;
    try
    {
        C = fix.materialMatrix->eval3D_matrix_C(Cmat,0,pt,z(0,0),MaterialOutput::Generic);
    }
    catch (...)
    {
        threw = true;
    }
    CHECK(!threw);
    // Guarded: on a (regression) throw, C is still the default-constructed
    // 0x0 matrix, and both C.allFinite() and C.reshape(3,3) (which
    // gsMatrix::reshape documents as assuming pre-allocated storage) would
    // read through a null gsMatrix::data().
    if (!threw)
    {
        CHECK(!C.allFinite());
        gsMatrix<real_t> C33 = C.reshape(3,3);
        for (index_t i = 0; i!=3; i++)
            for (index_t j = 0; j!=3; j++)
                CHECK(!math::isfinite(C33(i,j)));
    }
}

// Assembler-level status: the same isotropic collapse (s=1e-1) that fails
// at material level in C33_NonConvergent_Geometry, now observed through
// gsThinShellAssembler::assembleMatrix. Crosses the OpenMP boundary that
// arms A and B cannot -- the only part of this suite exercising the actual
// failure mode the run exists for.
TEST(C33_Assembler_Status)
{
    CompressibleNH fix;
    fix.mp.degreeElevate(1);
    fix.mp.uniformRefine();
    gsMultiBasis<real_t> dbasis(fix.mp);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(fix.mp);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

    gsVector<real_t> tmp(2); tmp.setZero();
    gsConstantFunction<real_t> force(tmp,2);

    gsThinShellAssembler<2, real_t, false> assembler(fix.mp,dbasis,bc,force,fix.materialMatrix);
    assembler.assemble();

    gsMultiPatch<real_t> mp_def_ctrl = fix.mp;
    ThinShellAssemblerStatus statusCtrl = assembler.assembleMatrix(mp_def_ctrl);
    CHECK(statusCtrl==ThinShellAssemblerStatus::Success);

    gsMultiPatch<real_t> mp_def_fail = fix.mp;
    mp_def_fail.patch(0).coefs() *= 1e-1;
    ThinShellAssemblerStatus statusFail = assembler.assembleMatrix(mp_def_fail);
    CHECK(statusFail==ThinShellAssemblerStatus::AssemblyError);
}

// Regression gate for the eval3D_pstretch packing (eval3D_pstretch_impl,
// compressible branch): result.col(colIdx) with colIdx=j*u.cols()+k must
// carry through-thickness station j's own C33, not another station's.
// Row 2 of the result is the sharpest instrument for this: inside
// _evalStretch (gsMaterialMatrixBaseDim.hpp:1370-1419) the third slot is
// assigned directly, stretches.at(2) = C(2,2), and then square-rooted --
// so result(2,colIdx) == sqrt(C33 as read for that packed column), exactly,
// with no eigenvector sort or near-identity branch (:1391) able to hide a
// mis-indexed read. Rows 0-1 come from a SelfAdjointEigenSolver on the
// in-plane block, which is refreshed per station regardless of the C33
// packing, so they cannot fail under this defect; they are kept only as a
// guard against a wider accidental change.
//
// The reference is built one through-thickness station at a time: with
// z.rows()==1 the packed index collapses to colIdx=k, so the buggy read
// C33s(0,k) and the correct read C33s(0,colIdx) are the same element for a
// single-station call -- the reference is valid whether or not the packing
// defect is present, which is what makes it usable as ground truth.
TEST(C33_MultiStation_pstretch)
{
    BentShellNH fix;

    // Deformed configuration: stretched in-plane (so C33 departs from 1) and
    // curved (one corner lifted out of plane, so the through-thickness
    // metric picks up the z-dependent Bcov/curvature terms that dim==3
    // carries and dim==2 cannot). Calibrated so all six eval3D_detF entries
    // below stay finite and strictly positive.
    gsMultiPatch<real_t> mp_def = fix.mp;
    mp_def.patch(0).coefs() *= 1.2;      // in-plane stretch, so C33 != 1
    mp_def.patch(0).coefs()(3,2) = 0.6;  // lift one corner out of plane -> curvature
    fix.materialMatrix->setDeformed(&mp_def);

    // Two in-plane points, u.cols()==2, exercise the k-index of the packing.
    gsMatrix<real_t> u(2,2);
    u.col(0) << 0.25, 0.25;
    u.col(1) << 0.75, 0.75;

    // Three through-thickness stations, and the two in-plane points carry
    // DIFFERENT station ladders: a transposed read (z(k,j) for z(j,k)) or a
    // mis-packed colIdx (k*z.rows()+j for j*u.cols()+k) would both be
    // invisible if both columns shared one ladder.
    const index_t nz = 3;
    gsMatrix<real_t> z(nz,2);
    z.col(0) << -0.5, 0.0, 0.5;
    z.col(1) << -0.25, 0.1, 0.4;

    // State sanity: an independent per-station probe untouched by the A1
    // defect (eval3D_detF reads C33s at colIdx already, gsMaterialMatrixNonlinear.hpp:996).
    gsMatrix<real_t> detF = fix.materialMatrix->eval3D_detF(0,u,z,MaterialOutput::Generic);
    CHECK_EQUAL(1, detF.rows());
    CHECK_EQUAL(u.cols()*z.rows(), detF.cols());
    for (index_t c = 0; c!=detF.cols(); c++)
    {
        CHECK(math::isfinite(detF(0,c)));
        CHECK(detF(0,c) > 0.0);
    }

    // Reference, one station at a time.
    std::vector<gsMatrix<real_t>> ref(nz);
    for (index_t j = 0; j!=nz; j++)
    {
        gsMatrix<real_t> zj(1,u.cols());
        for (index_t k = 0; k!=u.cols(); k++)
            zj(0,k) = z(j,k);
        ref[j] = fix.materialMatrix->eval3D_pstretch(0,u,zj);
        CHECK(ref[j].allFinite());
    }

    // Anti-vacuity: row 2 (=sqrt(C33)) must genuinely differ between the
    // first and last station at every in-plane point -- otherwise a packing
    // defect that swaps stations would be numerically invisible and the
    // gate below would pass vacuously.
    const real_t tol = 1e-10;
    for (index_t k = 0; k!=u.cols(); k++)
        CHECK( math::abs( ref[0](2,k) - ref[nz-1](2,k) ) > 1e-6 );

    // The comparison: one call with the full multi-station z.
    gsMatrix<real_t> all = fix.materialMatrix->eval3D_pstretch(0,u,z);
    CHECK_EQUAL(3, all.rows());
    CHECK_EQUAL(u.cols()*z.rows(), all.cols());
    CHECK(all.allFinite());

    for (index_t j = 0; j!=nz; j++)
    {
        for (index_t k = 0; k!=u.cols(); k++)
        {
            const index_t colIdx = j*u.cols()+k;
            // Primary: row 2 is the value the A1 defect mis-indexes.
            CHECK_CLOSE( ref[j](2,k), all(2,colIdx), tol );
            // Secondary guard only -- not the gate itself (see file header).
            CHECK_CLOSE( ref[j](0,k), all(0,colIdx), tol );
            CHECK_CLOSE( ref[j](1,k), all(1,colIdx), tol );
        }
    }
}

// Contract test, NOT a regression gate: pins the observable behaviour of a
// degenerate deformed metric (m_J0_sq==0) -- returns without throwing and
// yields a non-finite value. Collapsing the deformed square's y-coordinates
// onto a line makes Acov_def rank-deficient, so det(Gcov_def)==0 while
// det(Gcov_ori)==1 (undeformed unit square); in _getMetric
// (gsMaterialMatrixBaseDim.hpp) only det_def is zero, so the
// (det_ori==0 && det_def==0) branch is not taken, ratio=0 exactly, and
// "!(ratio>=0)" lets 0 through unchanged -- m_J0_sq==0, precisely.
// The pre-hardening code already returns this same NaN, the slow way: the
// Newton initial guess c(2,2)=pow(0,-1)=+inf is non-finite, so both S33 and
// C3333 are NaN from the first evaluation, math::lessthan(NaN,tol) is always
// false, the loop exhausts itmax, and the itmax branch writes NaN anyway.
// This test therefore passes identically with or without that hardening; it
// exists to pin the contract (no throw, NaN out), not to catch a regression.
TEST(C33_DegenerateDeformedMetric_YieldsNaN)
{
    CompressibleNH fix; // dim==2 is fine here: the defect is in _getMetric, not eval3D_pstretch's packing
    gsMultiPatch<real_t> mp_def = fix.mp;
    mp_def.patch(0).coefs().col(1).setZero(); // collapse the square onto a segment
    fix.materialMatrix->setDeformed(&mp_def);

    gsVector<real_t> pt(2); pt.setConstant(0.25);
    gsMatrix<real_t> z(1,1); z.setZero();

    bool threw = false;
    gsMatrix<real_t> detF;
    try
    {
        detF = fix.materialMatrix->eval3D_detF(0,pt,z,MaterialOutput::Generic);
    }
    catch (...)
    {
        threw = true;
    }
    CHECK(!threw);
    // Guarded as in C33_NonConvergent_Geometry: on a (regression) throw,
    // detF is still the default-constructed 0x0 matrix.
    if (!threw)
        CHECK(!math::isfinite(detF(0,0)));
}

}
