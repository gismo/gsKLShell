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

    Four further tests (dCdC_*) gate eval3D_dmatrix, the derivative dC/dC of
    the material matrix with respect to the deformed metric, at a state with a
    nonzero deformed shear metric g12 (pure shear x <- x + 0.4 y, so g12=0.4).
    For all materials but incompressible Neo-Hookean the derivative is a finite
    difference of the Voigt-input functor, whose third entry is 2*g12; the
    tests pin that the derivative is taken at the true metric. The reference is
    a central finite difference of the public eval3D_matrix_C (layout 3*c+j,
    j=0,1 the diagonal metric entries, j=2 the derivative w.r.t. 2*g12).
      - dCdC_Anchor_IncompNH_Analytic is an anchor, NOT a regression gate: the
        analytic incompressible Neo-Hookean overload does not use the functor,
        so it validates the reference's convention (including the shear-row
        scaling) against analytic code and passes with or without the defect.
      - dCdC_ShearScaling_IncompMR_vs_NH compares the functor path
        (incompressible Mooney-Rivlin, Ratio=1e8, which tends to Neo-Hookean)
        against the analytic Neo-Hookean dC/dC, all entries and the shear row.
      - dCdC_EvalPoint_IncompMR and dCdC_EvalPoint_CompNH are regression gates
        for the incompressible and the compressible branch of eval3D_matrix_C;
        each also checks that the state discriminates (the reference evaluated
        at g12/2 instead of g12 differs from the true one by far more than the
        tolerance).

    TEST(CauchyStress_Compressible_IsStressOverDetF) gates the compressible
    Cauchy stress: eval3D_CauchyStress must equal the second Piola-Kirchhoff
    components eval3D_stress divided column-wise by eval3D_detF, for
    compressible Neo-Hookean and Mooney-Rivlin on a twisted shell (unit
    square stretched by 1.2 with one corner lifted out of plane), where detF
    differs from 1 and varies over the in-plane points and through the
    thickness (z = 0 and z = 0.2). The incompressible Neo-Hookean material
    (J=1) on an isotropically stretched flat shell anchors the convention:
    its Cauchy stress equals its second Piola-Kirchhoff stress.

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

// Incompressible (plane-stress) material on the unit square, analytical
// implementation. Material 1 (NH) takes parameters {E,nu} and selects the
// analytic dC/dC overload; Material 3 (MR) takes {E,nu,Ratio} and selects the
// finite-difference one. Parameter functions are members for the lifetime
// reason given for CompressibleNH.
struct IncompressibleShell
{
    gsMultiPatch<real_t> mp;
    gsFunctionExpr<real_t> t;
    gsFunctionExpr<real_t> E;
    gsConstantFunction<real_t> nu;
    gsConstantFunction<real_t> ratio;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    IncompressibleShell(index_t material, bool withRatio, real_t ratioValue)
    :
    t("1e-3",2),
    E("1.0",2),
    nu(0.5,2),
    ratio(ratioValue,2)
    {
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) ); // unit square, degree 1
        mp.addAutoBoundaries();

        std::vector<gsFunctionSet<real_t>*> parameters;
        parameters.push_back(&E);
        parameters.push_back(&nu);
        if (withRatio)
            parameters.push_back(&ratio);

        gsOptionList options;
        options.addInt("Material","Material model: (1): NH, (3): MR",material);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        options.addSwitch("Compressibility","Compressibility: (true): Compressible",false);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
    }
};

struct IncompressibleNH : IncompressibleShell
{
    IncompressibleNH() : IncompressibleShell(1,false,0.0) {}
};

struct IncompressibleMR : IncompressibleShell
{
    // Ratio=c1/c2 -> infinity recovers Neo-Hookean (c2=mu/(Ratio+1)).
    explicit IncompressibleMR(real_t ratioValue = 1e8) : IncompressibleShell(3,true,ratioValue) {}
};

// Pure-shear deformed configuration x <- x + gamma*y of the unit square, and
// the true deformed metric G = g^T g at pt. Here a1=(1,0), a2=(gamma,1), so
// G = [1 gamma; gamma 1+gamma^2] and det G = 1.
struct ShearedState
{
    static constexpr real_t gamma = 0.4;

    gsMultiPatch<real_t> mp_def;
    gsVector<real_t> pt;
    gsMatrix<real_t> z;
    gsMatrix<real_t> G;

    ShearedState(const gsMultiPatch<real_t> & mp, gsMaterialMatrixBase<real_t> & mm)
    :
    mp_def(mp),
    pt(2),
    z(1,1)
    {
        mp_def.patch(0).coefs().col(0) += gamma * mp_def.patch(0).coefs().col(1);
        mm.setDeformed(&mp_def);

        pt.setConstant(0.25);
        z.setZero();

        gsMapData<real_t> map;
        map.flags = NEED_JACOBIAN;
        map.points = pt;
        static_cast<const gsFunction<real_t>&>(mp_def.patch(0)).computeMap(map);
        gsMatrix<real_t> g = map.jacobian(0);
        g.resize(2,2);
        G = g.transpose() * g;
    }

    // Unhalved Voigt metric [g11, g22, g12] as taken by eval3D_matrix_C
    gsMatrix<real_t> voigt() const
    {
        gsMatrix<real_t> Cmat(3,1);
        Cmat << G(0,0), G(1,1), G(0,1);
        return Cmat;
    }
};
constexpr real_t ShearedState::gamma;

// Central finite difference of the public eval3D_matrix_C, layout 3*c+j for
// output component c (column-major 3x3 material matrix) and input variable j.
// Variable j=2 is the derivative w.r.t. 2*g12, i.e. half the derivative w.r.t.
// the unhalved g12 entry. Cost: 6 evaluations of eval3D_matrix_C.
gsMatrix<real_t> dCdC_central_fd(const gsMaterialMatrixBase<real_t> & mm,
                                 const gsMatrix<real_t> & Cmat,
                                 const gsVector<real_t> & pt,
                                 const real_t delta)
{
    gsMatrix<real_t> result(27,1);
    for (index_t j = 0; j != 3; ++j)
    {
        gsMatrix<real_t> Cp = Cmat, Cm = Cmat;
        Cp(j,0) += delta;
        Cm(j,0) -= delta;
        const real_t h = (j==2 ? 4.0 : 2.0) * delta;
        gsMatrix<real_t> Cplus  = mm.eval3D_matrix_C(Cp,0,pt,0.0,MaterialOutput::Generic);
        gsMatrix<real_t> Cminus = mm.eval3D_matrix_C(Cm,0,pt,0.0,MaterialOutput::Generic);
        for (index_t c = 0; c != 9; ++c)
            result(3*c+j,0) = (Cplus(c,0) - Cminus(c,0)) / h;
    }
    return result;
}

// ||A-B||_max / ||B||_max
real_t relDiff(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
{
    return (A-B).cwiseAbs().maxCoeff() / B.cwiseAbs().maxCoeff();
}

// Shear row (variable j=2) of a 27x1 dC/dC vector, as a 1x9 row
gsMatrix<real_t> shearRow(const gsMatrix<real_t> & D)
{
    gsMatrix<real_t> D39 = D.reshape(3,9);
    return D39.row(2);
}

// Reference at the metric [g11, g22, g12/2], i.e. at the wrong evaluation point
gsMatrix<real_t> dCdC_central_fd_halfShear(const gsMaterialMatrixBase<real_t> & mm,
                                           const ShearedState & st,
                                           const real_t delta)
{
    gsMatrix<real_t> Cwrong = st.voigt();
    Cwrong(2,0) *= 0.5;
    return dCdC_central_fd(mm,Cwrong,st.pt,delta);
}

// Unit-square shell, analytical implementation, with parameters {E,nu} (NH)
// or {E,nu,Ratio} (MR), in-plane stretched by 'lambda' about the origin.
struct StretchedShell
{
    gsMultiPatch<real_t> mp;
    gsMultiPatch<real_t> mp_def;
    gsFunctionExpr<real_t> t;
    gsFunctionExpr<real_t> E;
    gsConstantFunction<real_t> nu;
    gsConstantFunction<real_t> ratio;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    StretchedShell(index_t material, bool compressible, real_t lambda)
    :
    t("1e-1",2),
    E("1.0",2),
    nu(compressible ? 0.3 : 0.5,2),
    ratio(1.5,2)
    {
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        mp.addAutoBoundaries();

        std::vector<gsFunctionSet<real_t>*> parameters;
        parameters.push_back(&E);
        parameters.push_back(&nu);
        if (material==3)
            parameters.push_back(&ratio);

        gsOptionList options;
        options.addInt("Material","Material model: (1): NH, (3): MR",material);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        options.addSwitch("Compressibility","Compressibility: (true): Compressible",compressible);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);

        mp_def = mp;
        mp_def.patch(0).coefs() *= lambda;
        materialMatrix->setDeformed(&mp_def);
    }
};

// Compressible shell surface (dim==3) with parameters {E,nu} (NH) or
// {E,nu,Ratio} (MR), stretched by 'lambda' with one corner lifted out of
// plane: x = (lambda*u, lambda*v, lift*u*v). The deformed Jacobian and the
// curvature (b12 != 0) vary over the patch, so detF depends on both the
// in-plane point and the through-thickness coordinate.
struct TwistedShell
{
    gsMultiPatch<real_t> mp;
    gsMultiPatch<real_t> mp_def;
    gsFunctionExpr<real_t> t;
    gsFunctionExpr<real_t> E;
    gsConstantFunction<real_t> nu;
    gsConstantFunction<real_t> ratio;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    TwistedShell(index_t material, real_t lambda, real_t lift)
    :
    t("1e-1",3),
    E("1.0",3),
    nu(0.3,3),
    ratio(1.5,3)
    {
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        mp.addAutoBoundaries();
        mp.embed(3);

        std::vector<gsFunctionSet<real_t>*> parameters;
        parameters.push_back(&E);
        parameters.push_back(&nu);
        if (material==3)
            parameters.push_back(&ratio);

        gsOptionList options;
        options.addInt("Material","Material model: (1): NH, (3): MR",material);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        options.addSwitch("Compressibility","Compressibility: (true): Compressible",true);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,options);

        mp_def = mp;
        mp_def.patch(0).coefs() *= lambda;
        mp_def.patch(0).coefs()(3,2) += lift;
        materialMatrix->setDeformed(&mp_def);
    }
};

// Checks eval3D_CauchyStress == eval3D_stress / eval3D_detF column-wise on a
// state where detF differs from 1 at every column and varies over points and
// through the thickness.
void checkCauchyIsStressOverDetF(const index_t material)
{
    TwistedShell fix(material,1.2,1.0);
    gsMatrix<real_t> u(2,3);
    u << 0.2, 0.5, 0.8,
         0.3, 0.6, 0.1;
    gsMatrix<real_t> z(2,u.cols());
    z.row(0).setConstant(0.0);
    z.row(1).setConstant(0.2);

    gsMatrix<real_t> S    = fix.materialMatrix->eval3D_stress(0,u,z,MaterialOutput::Generic);
    gsMatrix<real_t> detF = fix.materialMatrix->eval3D_detF(0,u,z,MaterialOutput::Generic);
    gsMatrix<real_t> sig  = fix.materialMatrix->eval3D_CauchyStress(0,u,z,MaterialOutput::Generic);
    CHECK(S.allFinite());
    CHECK(detF.allFinite());
    CHECK(sig.allFinite());
    CHECK_EQUAL(u.cols()*z.rows(),sig.cols());
    CHECK_EQUAL(S.cols(),sig.cols());
    CHECK_EQUAL(sig.cols(),detF.cols());

    CHECK(S.cwiseAbs().maxCoeff() > 0.0);
    CHECK((detF.array()-1.0).abs().minCoeff() >= 0.05);
    CHECK(detF.maxCoeff() - detF.minCoeff() > 1e-3);
    // detF at z=0 and z=0.2 differs at some point (curvature term)
    CHECK((detF.leftCols(u.cols()) - detF.rightCols(u.cols())).cwiseAbs().maxCoeff() > 1e-4);

    real_t err = 0;
    for (index_t c = 0; c != sig.cols(); ++c)
    {
        const real_t e = (sig.col(c) - S.col(c)/detF(0,c)).cwiseAbs().maxCoeff()
                         / S.col(c).cwiseAbs().maxCoeff();
        if (!(e <= err))
            err = e;
    }
    gsInfo<<"CauchyStress_Compressible material "<<material<<": err="<<err
          <<", detF in ["<<detF.minCoeff()<<","<<detF.maxCoeff()<<"]\n"<<detF<<"\n";
    CHECK(err < 1e-10);
}

// Flat reference patch (skewed, so G_ab varies with u; tilted by a fixed rotation for d = 3)
// and the affine deformation F = I + (lam-1) a a^T + (mu-1) b b^T + shear a b^T of its plane.
// The Cartesian stress is then uniform over the patch.
template<short_t d>
struct PSDirNL_Scene
{
    gsMultiPatch<real_t> ref, def;
    gsMatrix<real_t> R3;                ///< Plane orientation (identity for d = 2)
    gsVector<real_t,3> a, b;            ///< In-plane unit directions, padded to 3D

    PSDirNL_Scene()
    {
        const real_t deg = 3.14159265358979323846 / 180.0;
        R3 = gsMatrix<real_t>::Identity(3,3);
        if (d==3)
        {
            gsMatrix<real_t> Rz(3,3), Rx(3,3);
            Rz << std::cos(20*deg), -std::sin(20*deg), 0,
                  std::sin(20*deg),  std::cos(20*deg), 0,
                  0, 0, 1;
            Rx << 1, 0, 0,
                  0, std::cos(35*deg), -std::sin(35*deg),
                  0, std::sin(35*deg),  std::cos(35*deg);
            R3 = Rz * Rx;
        }

        ref.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        ref.patch(0).coefs().col(0) *= 2.0;
        ref.patch(0).coefs().col(1) *= 0.5;
        ref.patch(0).coefs()(3,0) += 0.7;
        if (d==3)
        {
            ref.embed(3);
            ref.patch(0).coefs() = ref.patch(0).coefs() * R3.transpose();
        }
        ref.addAutoBoundaries();
        def = ref;

        a = R3 * (gsVector<real_t,3>() << std::cos(25*deg),  std::sin(25*deg), 0).finished();
        b = R3 * (gsVector<real_t,3>() << -std::sin(25*deg), std::cos(25*deg), 0).finished();
    }

    void setState(const real_t lam, const real_t mu, const real_t shear)
    {
        const gsVector<real_t> ad = a.head(d), bd = b.head(d);
        gsMatrix<real_t> Fm = gsMatrix<real_t>::Identity(d,d);
        Fm += (lam-1)*ad*ad.transpose() + (mu-1)*bd*bd.transpose() + shear*ad*bd.transpose();
        def.patch(0).coefs() = ref.patch(0).coefs() * Fm.transpose();
    }

    // Reference Jacobian at pt, padded to 3x2
    gsMatrix<real_t> jac3(const gsMatrix<real_t> & pt) const
    {
        gsMapData<real_t> map;
        map.flags = NEED_JACOBIAN | NEED_VALUE;
        map.points = pt;
        static_cast<const gsFunction<real_t>&>(ref.patch(0)).computeMap(map);
        gsMatrix<real_t> J = map.jacobian(0);
        J.resize(d,2);
        gsMatrix<real_t> Jc = gsMatrix<real_t>::Zero(3,2);
        Jc.topRows(d) = J;
        return Jc;
    }
};

// Cartesian in-plane stress T = S^{kl} G_k G_l^T at one point and fibre height
template<short_t d>
gsMatrix<real_t> PSDirNL_T(const gsMaterialMatrixBase<real_t> & mm, const PSDirNL_Scene<d> & sc,
                           const gsMatrix<real_t> & pt, const real_t z)
{
    gsMatrix<real_t> zz(1,1);
    zz(0,0) = z;
    const gsMatrix<real_t> S1 = mm.eval3D_stress(0,pt,zz,MaterialOutput::Generic);
    const gsMatrix<real_t> Jc = sc.jac3(pt);
    const gsVector<real_t> G1 = Jc.col(0), G2 = Jc.col(1);
    gsMatrix<real_t> T = S1(0,0)*G1*G1.transpose() + S1(1,0)*G2*G2.transpose()
                       + S1(2,0)*(G1*G2.transpose() + G2*G1.transpose());
    return T;
}

// Ascending in-plane eigenvalues of T (tangent basis a, b)
template<short_t d>
gsVector<real_t> PSDirNL_inplane(const gsMatrix<real_t> & T, const PSDirNL_Scene<d> & sc)
{
    gsMatrix<real_t> Q(3,2);
    Q.col(0) = sc.a;
    Q.col(1) = sc.b;
    const gsMatrix<real_t> T2 = Q.transpose() * T * Q;
    gsEigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base > es;
    es.compute(T2);
    gsVector<real_t> ev = es.eigenvalues();
    return ev;
}

// Lateral stretch mu for which S_b = 0 at axial stretch lam (bisection on mu in [0.5,1])
template<short_t d>
real_t PSDirNL_findMu(gsMaterialMatrixBase<real_t> & mm, PSDirNL_Scene<d> & sc, const real_t lam)
{
    gsMatrix<real_t> pt0(2,1);
    pt0 << 0.5, 0.5;
    real_t lo = 0.5, hi = 1.0, mid = 0.5*(lo+hi);
    for (index_t it = 0; it != 100; ++it)
    {
        mid = 0.5*(lo+hi);
        sc.setState(lam,mid,0.0);
        mm.setDeformed(&sc.def);
        const gsMatrix<real_t> T = PSDirNL_T<d>(mm,sc,pt0,0.0);
        const real_t f = sc.b.dot(T*sc.b);
        if (std::abs(f) <= 1e-13*std::abs(sc.a.dot(T*sc.a)))
            break;
        if (f < 0)
            lo = mid;
        else
            hi = mid;
    }
    return mid;
}

template<short_t d>
void PSDirNL_check(const gsMaterialMatrixBase<real_t> & mm, const PSDirNL_Scene<d> & sc,
                   const bool uniaxial, const char * label)
{
    gsMatrix<real_t> u(2,3), z(2,3);
    u << 0.2, 0.5, 0.8,
         0.3, 0.6, 0.1;
    z.row(0).setConstant(0.0);
    z.row(1).setConstant(0.25);

    const gsMatrix<real_t> D = mm.eval3D_pstressDir(0,u,z,MaterialOutput::Generic);
    CHECK_EQUAL(9, D.rows());
    CHECK_EQUAL(u.cols()*z.rows(), D.cols());
    if (D.rows()!=9 || D.cols()!=u.cols()*z.rows())
        return;
    CHECK(D.allFinite());

    real_t wOrth = 0, wDet = 0, wN = 0, wSingle = 0, wEig = 0, wUnif = 0, wUni = 0;
    gsVector<real_t,3> N0, dmin0;
    for (index_t c = 0; c != D.cols(); ++c)
    {
        const index_t k = c % u.cols(), j = c / u.cols();
        const gsMatrix<real_t> uk = u.col(k);
        gsMatrix<real_t> zz(1,1);
        zz(0,0) = z(j,k);

        gsMatrix<real_t> F = D.col(c);
        F.resize(3,3);
        const gsVector<real_t,3> dmin = F.col(0), dmax = F.col(1), N = F.col(2);

        const gsMatrix<real_t> Jc = sc.jac3(uk);
        gsVector<real_t,3> n = gsVector<real_t,3>(Jc.col(0)).cross(gsVector<real_t,3>(Jc.col(1)));
        n.normalize();

        wOrth = std::max(wOrth, (F.transpose()*F - gsMatrix<real_t>::Identity(3,3)).cwiseAbs().maxCoeff());
        wDet  = std::max(wDet, 1.0 - F.determinant());
        wN    = std::max(wN, (N - n).cwiseAbs().maxCoeff());

        const gsMatrix<real_t> D1 = mm.eval3D_pstressDir(0,uk,zz,MaterialOutput::Generic);
        const gsMatrix<real_t> P1 = mm.eval3D_pstress(0,uk,zz,MaterialOutput::Generic);
        CHECK(D1.rows()==9 && D1.cols()==1);
        CHECK(P1.rows()==2 && P1.cols()==1);
        if (D1.rows()!=9 || D1.cols()!=1 || P1.rows()!=2 || P1.cols()!=1)
            return;
        wSingle = std::max(wSingle, (D.col(c) - D1.col(0)).cwiseAbs().maxCoeff());

        const gsMatrix<real_t> T = PSDirNL_T<d>(mm,sc,uk,z(j,k));
        const real_t scale = std::max(std::max(std::abs(P1(0,0)),std::abs(P1(1,0))), 1e-300);
        CHECK(P1(0,0) <= P1(1,0));
        wEig = std::max(wEig, (T*dmin - P1(0,0)*dmin).norm() / scale);
        wEig = std::max(wEig, (T*dmax - P1(1,0)*dmax).norm() / scale);

        if (c==0)
        {
            N0 = N;
            dmin0 = dmin;
        }
        wUnif = std::max(wUnif, (N - N0).cwiseAbs().maxCoeff());
        if (P1(1,0) - P1(0,0) > 1e-8*scale)
            wUnif = std::max(wUnif, 1.0 - std::abs(dmin.dot(dmin0)));

        if (uniaxial)
        {
            wUni = std::max(wUni, 1.0 - std::abs(dmin.dot(sc.b)));
            wUni = std::max(wUni, 1.0 - std::abs(dmax.dot(sc.a)));
            wUni = std::max(wUni, std::abs(dmin.dot(N)));
        }
    }
    gsInfo<<"PStressDir "<<label<<": |FtF-I|="<<wOrth<<", 1-det="<<wDet<<", |N-n|="<<wN
          <<", multi-vs-single="<<wSingle<<", eig residual="<<wEig<<", uniform="<<wUnif
          <<", uniaxial="<<wUni<<"\n";
    CHECK(wOrth   <= 1e-12);
    CHECK(wDet    <= 1e-12);
    CHECK(wN      <= 1e-12);
    CHECK(wSingle <= 1e-10);
    CHECK(wEig    <= 1e-10);
    CHECK(wUnif   <= 1e-10);
    CHECK(wUni    <= 1e-10);
}

template<short_t d>
void PSDirNL_run(const bool compressible)
{
    PSDirNL_Scene<d> sc;

    gsFunctionExpr<real_t> t("1e-1",d);
    gsFunctionExpr<real_t> E("1.0",d);
    gsConstantFunction<real_t> nu(compressible ? 0.3 : 0.5,d);
    std::vector<gsFunctionSet<real_t>*> parameters;
    parameters.push_back(&E);
    parameters.push_back(&nu);
    gsOptionList options;
    options.addInt("Material","Material model: (1): NH",1);
    options.addInt("Implementation","Implementation: (1): Analytical",1);
    options.addSwitch("Compressibility","Compressibility: (true): Compressible",compressible);
    gsMaterialMatrixBase<real_t>::uPtr mm = getMaterialMatrix<d,real_t>(sc.ref,t,parameters,options);
    const std::string name = compressible ? "NH compressible" : "NH incompressible";

    gsMatrix<real_t> pt0(2,1);
    pt0 << 0.5, 0.5;

    // generic state
    sc.setState(1.10, 0.95, 0.12);
    mm->setDeformed(&sc.def);
    gsVector<real_t> ev = PSDirNL_inplane<d>(PSDirNL_T<d>(*mm,sc,pt0,0.0),sc);
    CHECK(ev(1) - ev(0) >= 1e-2*std::max(std::abs(ev(0)),std::abs(ev(1))));
    PSDirNL_check<d>(*mm,sc,false,(name + " generic").c_str());

    // uniaxial state: lateral stretch such that S_b = 0
    const real_t mu = PSDirNL_findMu<d>(*mm,sc,1.2);
    sc.setState(1.2, mu, 0.0);
    mm->setDeformed(&sc.def);
    ev = PSDirNL_inplane<d>(PSDirNL_T<d>(*mm,sc,pt0,0.0),sc);
    CHECK(ev(1) > 0 && std::abs(ev(0)) <= 1e-8*ev(1));
    PSDirNL_check<d>(*mm,sc,true,(name + " uniaxial").c_str());

    // undeformed state
    sc.setState(1.0, 1.0, 0.0);
    mm->setDeformed(&sc.def);
    PSDirNL_check<d>(*mm,sc,false,(name + " undeformed").c_str());
}

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

// Anchor, not a regression gate: the analytic incompressible Neo-Hookean
// dC/dC does not go through the Voigt-input functor, so it is correct with or
// without the evaluation-point defect. It validates the finite-difference
// reference (layout 3*c+j and the shear-row scaling) against analytic code at
// a state with g12 != 0.
TEST(dCdC_Anchor_IncompNH_Analytic)
{
    const real_t tol = 1e-6;
    IncompressibleNH fix;
    ShearedState st(fix.mp,*fix.materialMatrix);
    CHECK_CLOSE(ShearedState::gamma, st.G(0,1), 1e-12);

    gsMatrix<real_t> D = fix.materialMatrix->eval3D_dmatrix(0,st.pt,st.z,MaterialOutput::Generic);
    gsMatrix<real_t> Dref = dCdC_central_fd(*fix.materialMatrix,st.voigt(),st.pt,1e-6);
    CHECK(D.allFinite());
    CHECK(Dref.allFinite());
    CHECK_EQUAL(27,D.rows());

    CHECK(shearRow(Dref).cwiseAbs().maxCoeff() > 0.0);
    const real_t errAll   = relDiff(D,Dref);
    const real_t errShear = relDiff(shearRow(D),shearRow(Dref));
    gsInfo<<"dCdC_Anchor_IncompNH_Analytic: err(all)="<<errAll<<", err(shear row)="<<errShear<<"\n";
    CHECK(errAll   < tol);
    CHECK(errShear < tol);
}

// Shear-row scaling against analytic incompressible Neo-Hookean: Mooney-Rivlin
// with Ratio=1e8 (c2/mu ~ 1e-8) takes the finite-difference path and must
// reproduce the analytic Neo-Hookean dC/dC, including the row of the
// derivative w.r.t. g12 (a halved or doubled shear row would differ by a
// factor 2, far above the tolerance).
TEST(dCdC_ShearScaling_IncompMR_vs_NH)
{
    const real_t tol = 1e-6;
    IncompressibleNH nh;
    IncompressibleMR mr;
    ShearedState stNH(nh.mp,*nh.materialMatrix);
    ShearedState stMR(mr.mp,*mr.materialMatrix);
    CHECK_CLOSE(ShearedState::gamma, stMR.G(0,1), 1e-12);

    gsMatrix<real_t> Dnh = nh.materialMatrix->eval3D_dmatrix(0,stNH.pt,stNH.z,MaterialOutput::Generic);
    gsMatrix<real_t> Dmr = mr.materialMatrix->eval3D_dmatrix(0,stMR.pt,stMR.z,MaterialOutput::Generic);
    CHECK(Dnh.allFinite());
    CHECK(Dmr.allFinite());

    CHECK(shearRow(Dnh).cwiseAbs().maxCoeff() > 0.0);
    const real_t errAll   = relDiff(Dmr,Dnh);
    const real_t errShear = relDiff(shearRow(Dmr),shearRow(Dnh));
    gsInfo<<"dCdC_ShearScaling_IncompMR_vs_NH: err(all)="<<errAll<<", err(shear row)="<<errShear<<"\n";
    CHECK(errAll   < tol);
    CHECK(errShear < tol);
}

// Regression gate, incompressible branch of eval3D_matrix_C: dC/dC from the
// functor path must equal the finite difference of eval3D_matrix_C at the
// true deformed metric, and that state must discriminate the evaluation
// point (g12/2 gives a derivative far from the true one).
TEST(dCdC_EvalPoint_IncompMR)
{
    const real_t tol = 1e-6;
    IncompressibleMR fix;
    ShearedState st(fix.mp,*fix.materialMatrix);
    CHECK_CLOSE(ShearedState::gamma, st.G(0,1), 1e-12);

    gsMatrix<real_t> D    = fix.materialMatrix->eval3D_dmatrix(0,st.pt,st.z,MaterialOutput::Generic);
    gsMatrix<real_t> Dref = dCdC_central_fd(*fix.materialMatrix,st.voigt(),st.pt,1e-6);
    gsMatrix<real_t> Dwrong = dCdC_central_fd_halfShear(*fix.materialMatrix,st,1e-6);
    CHECK(D.allFinite());
    CHECK(Dref.allFinite());

    CHECK(shearRow(Dref).cwiseAbs().maxCoeff() > 0.0);
    const real_t err      = relDiff(D,Dref);
    const real_t errWrong = relDiff(Dwrong,Dref);
    gsInfo<<"dCdC_EvalPoint_IncompMR: err(true point)="<<err<<", discrimination (g12/2 vs g12)="<<errWrong<<"\n";
    CHECK(err < tol);
    CHECK(errWrong > 100*tol);
}

// Regression gate, compressible branch of eval3D_matrix_C (through-thickness
// Newton solve inside every evaluation; the finite-difference step 1e-5 keeps
// the solve tolerance noise ~1e-10/delta at about 1e-5 relative).
TEST(dCdC_EvalPoint_CompNH)
{
    const real_t tol = 1e-4;
    CompressibleNH fix;
    ShearedState st(fix.mp,*fix.materialMatrix);
    CHECK_CLOSE(ShearedState::gamma, st.G(0,1), 1e-12);

    gsMatrix<real_t> D    = fix.materialMatrix->eval3D_dmatrix(0,st.pt,st.z,MaterialOutput::Generic);
    gsMatrix<real_t> Dref = dCdC_central_fd(*fix.materialMatrix,st.voigt(),st.pt,1e-5);
    gsMatrix<real_t> Dwrong = dCdC_central_fd_halfShear(*fix.materialMatrix,st,1e-5);
    CHECK(D.allFinite());
    CHECK(Dref.allFinite());
    CHECK(Dwrong.allFinite());

    CHECK(shearRow(Dref).cwiseAbs().maxCoeff() > 0.0);
    const real_t err      = relDiff(D,Dref);
    const real_t errWrong = relDiff(Dwrong,Dref);
    gsInfo<<"dCdC_EvalPoint_CompNH: err(true point)="<<err<<", discrimination (g12/2 vs g12)="<<errWrong<<"\n";
    CHECK(err < tol);
    CHECK(errWrong > 100*tol);
}

// The compressible Cauchy stress is S/detF.
TEST(CauchyStress_Compressible_IsStressOverDetF)
{
    checkCauchyIsStressOverDetF(1); // Neo-Hookean
    checkCauchyIsStressOverDetF(3); // Mooney-Rivlin
}

// Anchor: for the incompressible material (J=1) the Cauchy stress equals the
// second Piola-Kirchhoff stress, the convention the compressible test uses.
TEST(CauchyStress_Incompressible_EqualsStress)
{
    StretchedShell fix(1,false,1.2);
    gsMatrix<real_t> u(2,2);
    u << 0.2, 0.7,
         0.3, 0.4;
    gsMatrix<real_t> z(2,u.cols());
    z.row(0).setConstant(0.0);
    z.row(1).setConstant(0.03);
    gsMatrix<real_t> S   = fix.materialMatrix->eval3D_stress(0,u,z,MaterialOutput::Generic);
    gsMatrix<real_t> sig = fix.materialMatrix->eval3D_CauchyStress(0,u,z,MaterialOutput::Generic);
    CHECK(S.allFinite());
    CHECK(sig.allFinite());
    CHECK(sig.cwiseAbs().maxCoeff() > 0.0);
    CHECK_EQUAL(S.cols(),sig.cols());
    CHECK((sig-S).cwiseAbs().maxCoeff() <= 1e-10*S.cwiseAbs().maxCoeff());
}

TEST(PStressDir_NH_Compressible)
{
    PSDirNL_run<2>(true);
    PSDirNL_run<3>(true);
}

TEST(PStressDir_NH_Incompressible)
{
    PSDirNL_run<2>(false);
    PSDirNL_run<3>(false);
}

}
