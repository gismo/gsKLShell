/** @file gsMaterialMatrixLinear_test.cpp

    @brief Gates the 2D Saint Venant-Kirchhoff material (gsMaterialMatrixLinear<2>)
    against the same flat geometry evaluated as a 3D surface.

    On a flat surface in the z = 0 plane the 3D metric code reduces to the 2D one, so
    MatrixA and VectorN must coincide. The reference parametrisation is skewed and
    scaled so that the covariant metric A_{ab} and the contravariant metric A^{ab}
    differ by O(1); this catches a 2D code path that loads the wrong one of the two.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

using namespace gismo;

namespace {

// Saint Venant-Kirchhoff shell on the unit square (E = 1, nu = 0.3), with an
// orthonormal reference metric G = I or a sheared one, G = [1 0.4; 0.4 1.16]
// (coefficient column 0 += 0.4 * column 1). The parameter functions are members
// because the material keeps references to them.
struct EnergyLinear
{
    gsMultiPatch<real_t> mp, def;
    gsFunctionExpr<real_t> t, E, nu;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    explicit EnergyLinear(bool skewed)
    :
    t("1e-1",2), E("1.0",2), nu("0.3",2)
    {
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        if (skewed)
            mp.patch(0).coefs().col(0) += 0.4 * mp.patch(0).coefs().col(1);
        mp.addAutoBoundaries();

        std::vector<gsFunctionSet<real_t>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsOptionList options;
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);
        materialMatrix->setDeformed(&mp);
    }

    // Deformed configuration x <- x + 0.3*y, x <- 1.1*x, and the corner (1,1)
    // moved by (0.15,-0.1), so that the strain varies over the patch.
    void deform()
    {
        def = mp;
        def.patch(0).coefs().col(0) += 0.3 * def.patch(0).coefs().col(1);
        def.patch(0).coefs() *= 1.1;
        def.patch(0).coefs()(3,0) += 0.15;
        def.patch(0).coefs()(3,1) -= 0.1;
        materialMatrix->setDeformed(&def);
    }
};

// Metric G = g^T g of patch 0 of mp at pt
gsMatrix<real_t> energyLinMetric(const gsMultiPatch<real_t> & mp, const gsVector<real_t> & pt)
{
    gsMapData<real_t> map;
    map.flags = NEED_JACOBIAN;
    map.points = pt;
    static_cast<const gsFunction<real_t>&>(mp.patch(0)).computeMap(map);
    gsMatrix<real_t> g = map.jacobian(0);
    g.resize(2,2);
    return g.transpose() * g;
}

// Deformed Voigt metric [G11+2E1, G22+2E2, G12+E3] of the engineering strain E = [E11, E22, 2 E12]
gsMatrix<real_t> energyLinCm(const gsMatrix<real_t> & G, const gsVector<real_t> & E)
{
    gsMatrix<real_t> Cm(3,1);
    Cm(0,0) = G(0,0) + 2 * E(0);
    Cm(1,0) = G(1,1) + 2 * E(1);
    Cm(2,0) = G(0,1) +     E(2);
    return Cm;
}

// Test strains: every component is at least 1e-2 in magnitude
gsVector<real_t> energyLinStrain(index_t i)
{
    gsVector<real_t> E(3);
    if      (i == 0) E <<  0.05, -0.02,  0.03;
    else if (i == 1) E <<  0.12,  0.08, -0.06;
    else             E << -0.03, -0.02,  0.015;
    return E;
}

real_t energyLinAt(const gsMaterialMatrixBase<real_t> & mm, const gsMatrix<real_t> & Cm, const gsVector<real_t> & pt)
{
    return mm.eval3D_energy_C(Cm,0,pt,0.0,MaterialOutput::Generic)(0,0);
}

real_t energyLinMaxC(const gsMaterialMatrixBase<real_t> & mm, const gsMatrix<real_t> & Cm, const gsVector<real_t> & pt)
{
    return mm.eval3D_matrix_C(Cm,0,pt,0.0,MaterialOutput::Generic).cwiseAbs().maxCoeff();
}

} // namespace

SUITE(gsMaterialMatrixLinear_test)                 // The suite should have the same name as the file
{
    /// Reference and deformed flat patches with a clearly non-orthonormal reference metric.
    static void LinSvK_geometries(gsMultiPatch<> & ref, gsMultiPatch<> & def)
    {
        ref.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        ref.patch(0).coefs().col(0) *= 2.0;
        ref.patch(0).coefs().col(1) *= 0.5;
        ref.patch(0).coefs()(3,0) += 0.7; // skew
        ref.addAutoBoundaries();

        def = ref;
        def.patch(0).coefs().col(0) *= 1.05;
        def.patch(0).coefs()(1,1)   += 0.04;
        def.patch(0).coefs()(3,0)   += 0.08;
        def.patch(0).coefs()(3,1)   += 0.03;
    }

    /// Relative max-abs discrepancy
    static real_t LinSvK_relDiff(const gsMatrix<> & a, const gsMatrix<> & b)
    {
        return (a-b).cwiseAbs().maxCoeff() / a.cwiseAbs().maxCoeff();
    }

    template<short_t d>
    static void LinSvK_evaluate(const gsMultiPatch<> & ref, const gsMultiPatch<> & def,
                                const gsVector<> & pt, const real_t z,
                                gsMatrix<> & A, gsMatrix<> & N)
    {
        gsFunctionExpr<> t("0.01",2);
        gsFunctionExpr<> E("2.0e5",2);
        gsFunctionExpr<> nu("0.3",2);
        std::vector<gsFunctionSet<>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;

        gsOptionList options;
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);

        gsMaterialMatrixBase<real_t>::uPtr mm = getMaterialMatrix<d,real_t>(ref,t,parameters,options);
        mm->setDeformed(&def);

        gsMatrix<> zmat(1,1);
        zmat(0,0) = z;
        A = mm->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        N = mm->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
    }

    TEST(SvK_2D_equals_3D_embedding)
    {
        gsMultiPatch<> ref2, def2;
        LinSvK_geometries(ref2,def2);

        gsMultiPatch<> ref3 = ref2, def3 = def2;
        ref3.embed(3);
        def3.embed(3);

        gsVector<> pt(2);
        pt << 0.3, 0.6;

        // Non-vacuity: A_{ab} and A^{ab} differ by O(1) at the test point
        gsMapData<> map;
        map.flags = NEED_JACOBIAN | NEED_VALUE;
        map.points = pt;
        static_cast<const gsFunction<>&>(ref2.patch(0)).computeMap(map);
        gsMatrix<> J = map.jacobian(0);
        J.resize(2,2);
        gsMatrix<> Acov = J.transpose() * J;
        gsMatrix<> Acon = Acov.inverse();
        CHECK(Acov.cwiseAbs().maxCoeff() > 0);
        CHECK((Acov - Acon).cwiseAbs().maxCoeff() >= 0.1);

        gsMatrix<> A2, N2, A3, N3;
        LinSvK_evaluate<2>(ref2,def2,pt,0.0,A2,N2);
        LinSvK_evaluate<3>(ref3,def3,pt,0.0,A3,N3);

        const real_t dA = LinSvK_relDiff(A3,A2);
        const real_t dN = LinSvK_relDiff(N3,N2);
        gsInfo<<"SvK 2D vs 3D: rel diff MatrixA = "<<dA<<", VectorN = "<<dN<<"\n";
        CHECK(N3.cwiseAbs().maxCoeff() > 0); // deformed state is non-trivial
        CHECK(dA <= 1e-10);
        CHECK(dN <= 1e-10);
    }

    /// Flat reference patch (skewed, so G_ab varies with u; tilted by a fixed rotation for d = 3)
    /// and the affine deformation F = I + (lam-1) a a^T + (mu-1) b b^T + shear a b^T of its plane.
    /// The Cartesian stress is then uniform over the patch.
    template<short_t d>
    struct PSDirLin_Scene
    {
        gsMultiPatch<> ref, def;
        gsMatrix<> R3;                  ///< Plane orientation (identity for d = 2)
        gsVector<real_t,3> a, b;        ///< In-plane unit directions, padded to 3D

        explicit PSDirLin_Scene(const bool unitSquare = false)
        {
            const real_t deg = 3.14159265358979323846 / 180.0;
            R3 = gsMatrix<>::Identity(3,3);
            if (d==3)
            {
                gsMatrix<> Rz(3,3), Rx(3,3);
                Rz << std::cos(20*deg), -std::sin(20*deg), 0,
                      std::sin(20*deg),  std::cos(20*deg), 0,
                      0, 0, 1;
                Rx << 1, 0, 0,
                      0, std::cos(35*deg), -std::sin(35*deg),
                      0, std::sin(35*deg),  std::cos(35*deg);
                R3 = Rz * Rx;
            }

            ref.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
            if (!unitSquare)
            {
                ref.patch(0).coefs().col(0) *= 2.0;
                ref.patch(0).coefs().col(1) *= 0.5;
                ref.patch(0).coefs()(3,0) += 0.7;
            }
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
            const gsVector<> ad = a.head(d), bd = b.head(d);
            gsMatrix<> Fm = gsMatrix<>::Identity(d,d);
            Fm += (lam-1)*ad*ad.transpose() + (mu-1)*bd*bd.transpose() + shear*ad*bd.transpose();
            def.patch(0).coefs() = ref.patch(0).coefs() * Fm.transpose();
        }

        /// Reference Jacobian at pt, padded to 3x2
        gsMatrix<> jac3(const gsMatrix<> & pt) const
        {
            gsMapData<> map;
            map.flags = NEED_JACOBIAN | NEED_VALUE;
            map.points = pt;
            static_cast<const gsFunction<>&>(ref.patch(0)).computeMap(map);
            gsMatrix<> J = map.jacobian(0);
            J.resize(d,2);
            gsMatrix<> Jc = gsMatrix<>::Zero(3,2);
            Jc.topRows(d) = J;
            return Jc;
        }
    };

    /// Cartesian in-plane stress T = S^{kl} G_k G_l^T at one point and fibre height
    template<short_t d>
    static gsMatrix<> PSDirLin_T(const gsMaterialMatrixBase<real_t> & mm, const PSDirLin_Scene<d> & sc,
                                 const gsMatrix<> & pt, const real_t z)
    {
        gsMatrix<> zz(1,1);
        zz(0,0) = z;
        const gsMatrix<> S1 = mm.eval3D_stress(0,pt,zz,MaterialOutput::Generic);
        const gsMatrix<> Jc = sc.jac3(pt);
        const gsVector<> G1 = Jc.col(0), G2 = Jc.col(1);
        gsMatrix<> T = S1(0,0)*G1*G1.transpose() + S1(1,0)*G2*G2.transpose()
                     + S1(2,0)*(G1*G2.transpose() + G2*G1.transpose());
        return T;
    }

    /// Ascending in-plane eigenvalues of T (tangent basis a, b)
    template<short_t d>
    static gsVector<> PSDirLin_inplane(const gsMatrix<> & T, const PSDirLin_Scene<d> & sc)
    {
        gsMatrix<> Q(3,2);
        Q.col(0) = sc.a;
        Q.col(1) = sc.b;
        const gsMatrix<> T2 = Q.transpose() * T * Q;
        gsEigen::SelfAdjointEigenSolver< gsMatrix<>::Base > es;
        es.compute(T2);
        gsVector<> ev = es.eigenvalues();
        return ev;
    }

    template<short_t d>
    static void PSDirLin_check(const gsMaterialMatrixBase<real_t> & mm, const PSDirLin_Scene<d> & sc,
                               const bool uniaxial, const char * label)
    {
        gsMatrix<> u(2,3), z(2,3);
        u << 0.2, 0.5, 0.8,
             0.3, 0.6, 0.1;
        z.row(0).setConstant(0.0);
        z.row(1).setConstant(0.25);

        const gsMatrix<> D = mm.eval3D_pstressDir(0,u,z,MaterialOutput::Generic);
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
            const gsMatrix<> uk = u.col(k);
            gsMatrix<> zz(1,1);
            zz(0,0) = z(j,k);

            gsMatrix<> F = D.col(c);
            F.resize(3,3);
            const gsVector<real_t,3> dmin = F.col(0), dmax = F.col(1), N = F.col(2);

            const gsMatrix<> Jc = sc.jac3(uk);
            gsVector<real_t,3> n = gsVector<real_t,3>(Jc.col(0)).cross(gsVector<real_t,3>(Jc.col(1)));
            n.normalize();

            wOrth = std::max(wOrth, (F.transpose()*F - gsMatrix<>::Identity(3,3)).cwiseAbs().maxCoeff());
            wDet  = std::max(wDet, 1.0 - F.determinant());
            wN    = std::max(wN, (N - n).cwiseAbs().maxCoeff());

            const gsMatrix<> D1 = mm.eval3D_pstressDir(0,uk,zz,MaterialOutput::Generic);
            const gsMatrix<> P1 = mm.eval3D_pstress(0,uk,zz,MaterialOutput::Generic);
            CHECK(D1.rows()==9 && D1.cols()==1);
            CHECK(P1.rows()==2 && P1.cols()==1);
            if (D1.rows()!=9 || D1.cols()!=1 || P1.rows()!=2 || P1.cols()!=1)
                return;
            wSingle = std::max(wSingle, (D.col(c) - D1.col(0)).cwiseAbs().maxCoeff());

            const gsMatrix<> T = PSDirLin_T<d>(mm,sc,uk,z(j,k));
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
    static void PSDirLin_run()
    {
        PSDirLin_Scene<d> sc;

        gsFunctionExpr<> t("0.01",2);
        gsFunctionExpr<> E("2.0e5",2);
        gsFunctionExpr<> nu("0.3",2);
        std::vector<gsFunctionSet<real_t>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;
        gsOptionList options;
        options.addInt("Material","Material model: (0): SvK",0);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        gsMaterialMatrixBase<real_t>::uPtr mm = getMaterialMatrix<d,real_t>(sc.ref,t,parameters,options);

        gsMatrix<> pt0(2,1);
        pt0 << 0.5, 0.5;

        // generic state
        sc.setState(1.10, 0.95, 0.12);
        mm->setDeformed(&sc.def);
        gsVector<> ev = PSDirLin_inplane<d>(PSDirLin_T<d>(*mm,sc,pt0,0.0),sc);
        CHECK(ev(1) - ev(0) >= 1e-2*std::max(std::abs(ev(0)),std::abs(ev(1))));
        PSDirLin_check<d>(*mm,sc,false,"SvK generic");

        // uniaxial state: S_b = 0 for E_bb = -nu E_aa
        sc.setState(1.2, std::sqrt(1.0 - 0.3*(1.2*1.2 - 1.0)), 0.0);
        mm->setDeformed(&sc.def);
        ev = PSDirLin_inplane<d>(PSDirLin_T<d>(*mm,sc,pt0,0.0),sc);
        CHECK(ev(1) > 0 && std::abs(ev(0)) <= 1e-8*ev(1));
        PSDirLin_check<d>(*mm,sc,true,"SvK uniaxial");

        // undeformed state
        sc.setState(1.0, 1.0, 0.0);
        mm->setDeformed(&sc.def);
        PSDirLin_check<d>(*mm,sc,false,"SvK undeformed");
    }

    /// Calls _evalPStress directly on S = diag(1, 1e-15): the eigenvalues of S^{kl} G_k G_l^T
    /// are (0, 1e-15, 1), so a 3x3 solve cannot tell the in-plane near-zero principal stress
    /// from the out-of-plane zero.
    template<short_t d>
    static void PSDirLin_zeroPrincipal()
    {
        PSDirLin_Scene<d> sc(true);

        gsFunctionExpr<> t("0.01",2);
        gsFunctionExpr<> E("2.0e5",2);
        gsFunctionExpr<> nu("0.3",2);
        std::vector<gsFunctionSet<real_t>*> parameters(2);
        parameters[0] = &E;
        parameters[1] = &nu;
        gsOptionList options;
        options.addInt("Material","Material model: (0): SvK",0);
        options.addInt("Implementation","Implementation: (1): Analytical",1);
        gsMaterialMatrixBase<real_t>::uPtr mm = getMaterialMatrix<d,real_t>(sc.ref,t,parameters,options);
        mm->setDeformed(&sc.def);

        gsMatrix<> pt(2,1), zz(1,1);
        pt << 0.5, 0.5;
        zz(0,0) = 0.0;
        mm->eval3D_stress(0,pt,zz,MaterialOutput::Generic); // loads the metric of the point

        gsMaterialMatrixBaseDim<d,real_t> * bd = dynamic_cast<gsMaterialMatrixBaseDim<d,real_t>*>(mm.get());
        CHECK(bd!=nullptr);
        if (bd==nullptr)
            return;

        gsMatrix<> S = gsMatrix<>::Zero(3,3);
        S(0,0) = 1.0;
        S(1,1) = 1e-15;
        const std::pair<gsVector<real_t>,gsMatrix<real_t> > res = bd->_evalPStress(S);
        CHECK(res.second.rows()==3 && res.second.cols()==3);
        if (res.second.rows()!=3 || res.second.cols()!=3)
            return;

        const gsVector<real_t,3> e2 = sc.R3.col(1), e3 = sc.R3.col(2);
        const gsVector<real_t,3> dmin = res.second.col(0), N = res.second.col(2);
        gsInfo<<"PStressDir _evalPStress d="<<d<<": |dmin.e2|="<<std::abs(dmin.dot(e2))
              <<", |N-e3|="<<(N-e3).cwiseAbs().maxCoeff()<<", values "<<res.first.transpose()<<"\n";
        CHECK(std::abs(dmin.dot(e2)) >= 1.0 - 1e-12);
        CHECK((N - e3).cwiseAbs().maxCoeff() <= 1e-12);
        CHECK(std::abs(res.first(0)) <= 1e-14);
        CHECK(std::abs(res.first(1) - 1.0) <= 1e-12);
    }

    TEST(PStressDir_SvK_2D)
    {
        PSDirLin_run<2>();
    }

    TEST(PStressDir_SvK_3D)
    {
        PSDirLin_run<3>();
    }

    TEST(PStressDir_EvalPStress_ZeroPrincipal)
    {
        PSDirLin_zeroPrincipal<2>();
        PSDirLin_zeroPrincipal<3>();
    }

    typedef std::pair<gsVector<real_t>,gsMatrix<real_t> > PGuard_Pair;

    /// SvK material on a copy of a reference patch, with direct access to the BaseDim interface.
    template<short_t d>
    struct PGuard_Setup
    {
        gsFunctionExpr<> t, E, nu;
        gsMultiPatch<> ref;
        gsMaterialMatrixBase<real_t>::uPtr mm;
        gsMaterialMatrixBaseDim<d,real_t> * bd;

        explicit PGuard_Setup(const gsMultiPatch<> & r)
        : t("0.01",2), E("2.0e5",2), nu("0.3",2), ref(r), bd(nullptr)
        {
            std::vector<gsFunctionSet<real_t>*> parameters(2);
            parameters[0] = &E;
            parameters[1] = &nu;
            gsOptionList options;
            options.addInt("Material","Material model: (0): SvK",0);
            options.addInt("Implementation","Implementation: (1): Analytical",1);
            mm = getMaterialMatrix<d,real_t>(ref,t,parameters,options);
            bd = dynamic_cast<gsMaterialMatrixBaseDim<d,real_t>*>(mm.get());
        }

        PGuard_Setup(const PGuard_Setup &) = delete;
        PGuard_Setup & operator=(const PGuard_Setup &) = delete;

        /// Loads the reference metric of the parametric point pt (2x1) at z = 0
        void load(const gsMatrix<> & pt) const
        {
            bd->_computePoints(0,pt);
            bd->_getMetricUndeformed(0,0.0);
        }
    };

    /// Plane orientation used by PSDirLin_Scene<3>
    static gsMatrix<> PGuard_R3()
    {
        PSDirLin_Scene<3> sc(true);
        return sc.R3;
    }

    /// Embeds and tilts a 2D patch for d = 3, and sets up its boundaries
    template<short_t d>
    static void PGuard_orient(gsMultiPatch<> & mp)
    {
        if (d==3)
        {
            const gsMatrix<> R3 = PGuard_R3();
            mp.embed(3);
            mp.patch(0).coefs() = mp.patch(0).coefs() * R3.transpose();
        }
        mp.addAutoBoundaries();
    }

    /// Unit square whose edge v = 1 is collapsed to the point (0,1); for the near collapse the
    /// collapsed edge is opened by delta in x
    template<short_t d>
    static gsMultiPatch<> PGuard_collapsed(const bool exact, const real_t delta)
    {
        gsMultiPatch<> mp;
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp.patch(0).coefs().row(3) = mp.patch(0).coefs().row(2);
        if (!exact)
            mp.patch(0).coefs()(3,0) += delta;
        PGuard_orient<d>(mp);
        return mp;
    }

    /// Unit square stretched by 1e-6 in x (regular but strongly anisotropic frame)
    template<short_t d>
    static gsMultiPatch<> PGuard_anisotropic()
    {
        gsMultiPatch<> mp;
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        mp.patch(0).coefs().col(0) *= 1e-6;
        PGuard_orient<d>(mp);
        return mp;
    }

    /// Reference Jacobian of patch 0 at pt, padded to 3x2
    static gsMatrix<> PGuard_jac3(const gsMultiPatch<> & mp, const gsMatrix<> & pt)
    {
        const index_t dd = mp.patch(0).geoDim();
        gsMapData<> map;
        map.flags = NEED_JACOBIAN | NEED_VALUE;
        map.points = pt;
        static_cast<const gsFunction<>&>(mp.patch(0)).computeMap(map);
        gsMatrix<> J = map.jacobian(0);
        J.resize(dd,2);
        gsMatrix<> Jc = gsMatrix<>::Zero(3,2);
        Jc.topRows(dd) = J;
        return Jc;
    }

    /// Non-vacuity of a collapsed frame: exactly zero G1 (exact), or a nonzero cross product
    /// below the relative threshold used by _evalPStress/_evalPStrain (near)
    static void PGuard_checkCollapse(const gsMatrix<> & Jc, const bool exact, const char * label)
    {
        const gsVector<real_t,3> G1 = Jc.col(0), G2 = Jc.col(1);
        const real_t c = G1.cross(G2).norm();
        const real_t s = G1.squaredNorm() + G2.squaredNorm();
        const real_t tol = 10 * std::numeric_limits<real_t>::epsilon();
        gsInfo<<"  "<<label<<" frame: |G1|="<<G1.norm()<<", |G1xG2|/(|G1|^2+|G2|^2)="<<c/s
              <<", threshold "<<tol<<"\n";
        if (exact)
        {
            CHECK(G1.norm() == 0);
        }
        else
        {
            CHECK(c > 0);
            CHECK(c <= tol*s);
        }
    }

    template<class M>
    static bool PGuard_allNaN(const M & m)
    {
        for (index_t j = 0; j != m.cols(); ++j)
            for (index_t i = 0; i != m.rows(); ++i)
                if (!std::isnan(m(i,j)))
                    return false;
        return true;
    }

    static bool PGuard_shape(const PGuard_Pair & r, const index_t nvals)
    {
        const bool ok = (r.first.size()==nvals && r.second.rows()==3 && r.second.cols()==3);
        CHECK(ok);
        return ok;
    }

    /// Loads the point and evaluates _evalPStress; Jc returns the reference Jacobian of the point
    template<short_t d>
    static PGuard_Pair PGuard_evalStress(const gsMultiPatch<> & mp, const gsMatrix<> & pt,
                                         const gsMatrix<> & S, gsMatrix<> & Jc)
    {
        PGuard_Setup<d> su(mp);
        CHECK(su.bd != nullptr);
        if (su.bd == nullptr)
            return PGuard_Pair();
        su.load(pt);
        Jc = PGuard_jac3(mp,pt);
        return su.bd->_evalPStress(S);
    }

    /// Loads the point and evaluates _evalPStrain; Jc returns the reference Jacobian of the point
    template<short_t d>
    static PGuard_Pair PGuard_evalStrain(const gsMultiPatch<> & mp, const gsMatrix<> & pt,
                                         const gsMatrix<> & E, gsMatrix<> & Jc)
    {
        PGuard_Setup<d> su(mp);
        CHECK(su.bd != nullptr);
        if (su.bd == nullptr)
            return PGuard_Pair();
        su.load(pt);
        Jc = PGuard_jac3(mp,pt);
        return su.bd->_evalPStrain(E);
    }

    /// Cartesian tensor T = sum_{kl} S^{kl} G_k G_l^T
    static gsMatrix<> PGuard_T(const gsMatrix<> & S, const gsMatrix<> & Jc)
    {
        gsMatrix<> T = gsMatrix<>::Zero(3,3);
        for (index_t k = 0; k != 2; ++k)
            for (index_t l = 0; l != 2; ++l)
                T += S(k,l) * Jc.col(k) * Jc.col(l).transpose();
        return T;
    }

    /// Principal frame of T checked against the ascending eigenvalues of Q^T T Q (Q: 3x2 tangent basis)
    static void PGuard_checkFrame(const char * label, const PGuard_Pair & res, const gsMatrix<> & T,
                                  const gsMatrix<> & Q, const real_t tol, const gsVector<real_t,3> * N)
    {
        if (!PGuard_shape(res,2))
            return;
        const gsMatrix<> T2 = Q.transpose() * T * Q;
        gsEigen::SelfAdjointEigenSolver< gsMatrix<>::Base > es;
        es.compute(T2);
        const gsVector<> lam = es.eigenvalues();
        const real_t scale = std::max(lam.cwiseAbs().maxCoeff(), 1e-300);

        const gsMatrix<> F = res.second;
        const bool finite = res.first.allFinite() && F.allFinite();
        real_t wEig = 0;
        for (index_t c = 0; c != 2; ++c)
            wEig = std::max(wEig, (T*F.col(c) - res.first(c)*F.col(c)).norm() / scale);
        const real_t wVal  = (res.first - lam).cwiseAbs().maxCoeff() / scale;
        const real_t wOrth = (F.transpose()*F - gsMatrix<>::Identity(3,3)).cwiseAbs().maxCoeff();
        gsInfo<<"  "<<label<<": values "<<res.first.transpose()<<" (reference "<<lam.transpose()
              <<"), |FtF-I|="<<wOrth<<", value err="<<wVal<<", eig residual="<<wEig<<"\n";
        CHECK(finite);
        CHECK(res.first(0) <= res.first(1));
        CHECK(wOrth <= 1e-12);
        CHECK(wVal  <= tol);
        CHECK(wEig  <= tol);
        if (N != nullptr)
            CHECK((F.col(2) - *N).cwiseAbs().maxCoeff() == 0);
    }

    static void PGuard_checkAllNaN(const char * label, const PGuard_Pair & res)
    {
        if (!PGuard_shape(res,2))
            return;
        gsInfo<<"  "<<label<<": values "<<res.first.transpose()<<"\n";
        CHECK(PGuard_allNaN(res.first));
        CHECK(PGuard_allNaN(res.second));
    }

    /// E with a zero in-plane principal strain (1e-15): a 3x3 solve cannot tell it from the out-of-plane zero
    template<short_t d>
    static void PGuard_zeroPrincipal()
    {
        PSDirLin_Scene<d> sc(true);
        PGuard_Setup<d> su(sc.ref);
        CHECK(su.bd != nullptr);
        if (su.bd == nullptr)
            return;
        gsMatrix<> pt(2,1);
        pt << 0.5, 0.5;
        su.load(pt);

        gsMatrix<> E = gsMatrix<>::Zero(3,3);
        E(0,0) = 1.0;
        E(1,1) = 1e-15;
        E(2,2) = 0.123;
        const PGuard_Pair res = su.bd->_evalPStrain(E);
        if (!PGuard_shape(res,3))
            return;

        const gsVector<real_t,3> v0 = res.second.col(0), v1 = res.second.col(1), nn = res.second.col(2);
        const gsVector<real_t,3> e1 = sc.R3.col(0), e2 = sc.R3.col(1), e3 = sc.R3.col(2);
        gsInfo<<"PStrain_ZeroPrincipal d="<<d<<": values "<<res.first.transpose()
              <<", |col0.e2|="<<std::abs(v0.dot(e2))<<", |col1.e1|="<<std::abs(v1.dot(e1))
              <<", |col2-e3|="<<(nn-e3).cwiseAbs().maxCoeff()<<"\n";
        CHECK(std::abs(v0.dot(e2)) >= 1.0 - 1e-12);
        CHECK(std::abs(v1.dot(e1)) >= 1.0 - 1e-12);
        CHECK(std::abs(res.first(0) - 1e-15) <= 1e-14);
        CHECK(std::abs(res.first(1) - 1.0) <= 1e-12);
        CHECK(res.first(2) == 0.123);
        CHECK((nn - e3).cwiseAbs().maxCoeff() <= 1e-12);
    }

    /// Non-finite entries of E give NaN principal values and directions; E(2,2) and the normal pass through
    template<short_t d>
    static void PGuard_nonFinite()
    {
        PSDirLin_Scene<d> sc(true);
        PGuard_Setup<d> su(sc.ref);
        CHECK(su.bd != nullptr);
        if (su.bd == nullptr)
            return;
        gsMatrix<> pt(2,1);
        pt << 0.5, 0.5;
        su.load(pt);

        const real_t bad[2] = { std::numeric_limits<real_t>::quiet_NaN(),
                                std::numeric_limits<real_t>::infinity() };
        const char * name[2] = { "NaN", "inf" };
        for (index_t c = 0; c != 2; ++c)
        {
            gsMatrix<> E = gsMatrix<>::Zero(3,3);
            E(0,0) = bad[c];
            E(1,1) = 0.5;
            E(0,1) = E(1,0) = 0.1;
            E(2,2) = 0.123;
            const PGuard_Pair res = su.bd->_evalPStrain(E);
            if (!PGuard_shape(res,3))
                return;
            const gsVector<real_t,3> nn = res.second.col(2), e3 = sc.R3.col(2);
            gsInfo<<"PStrain_NonFinite d="<<d<<" E(0,0)="<<name[c]<<": values "<<res.first.transpose()
                  <<", |col2-e3|="<<(nn-e3).cwiseAbs().maxCoeff()<<"\n";
            CHECK(res.first(2) == 0.123);
            CHECK((nn - e3).cwiseAbs().maxCoeff() <= 1e-12);
            CHECK(std::isnan(res.first(0)));
            CHECK(std::isnan(res.first(1)));
            CHECK(PGuard_allNaN(res.second.leftCols(2)));
        }
    }

    /// Finite input on a nearly degenerate parallelogram frame (angle theta between G1 and G2) with the
    /// physical strain b b^T: probes whether the out-of-plane entry survives
    static void PGuard_finiteNoZero()
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        const real_t deg = 3.14159265358979323846 / 180.0;
        const gsMatrix<> R3 = PGuard_R3();
        const gsVector<real_t,3> b = R3 * (gsVector<real_t,3>() << -std::sin(25*deg), std::cos(25*deg), 0).finished();
        const real_t thetas[4] = { 1e-2, 1e-3, 1e-4, 1e-5 };
        gsMatrix<> pt(2,1);
        pt << 0.5, 0.5;

        for (index_t i = 0; i != 4; ++i)
        {
            const real_t th = thetas[i];
            gsMultiPatch<> mp;
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
            mp.patch(0).coefs().row(2) << std::cos(th), std::sin(th);
            mp.patch(0).coefs().row(3) << 1.0 + std::cos(th), std::sin(th);
            PGuard_orient<3>(mp);

            const gsMatrix<> Jc = PGuard_jac3(mp,pt);
            const gsVector<real_t,3> G1 = Jc.col(0), G2 = Jc.col(1);
            gsVector<real_t,3> n = G1.cross(G2);
            n.normalize();
            const real_t g1b = G1.dot(b), g2b = G2.dot(b);
            gsMatrix<> E = gsMatrix<>::Zero(3,3);
            E(0,0) = g1b*g1b;
            E(1,1) = g2b*g2b;
            E(0,1) = E(1,0) = g1b*g2b;
            E(2,2) = 0.123;

            gsMatrix<> Jc2;
            const PGuard_Pair res = PGuard_evalStrain<3>(mp,pt,E,Jc2);
            if (!PGuard_shape(res,3))
                return;
            const gsVector<real_t,3> v0 = res.second.col(0), v1 = res.second.col(1), nn = res.second.col(2);
            const real_t tolDir = 1e3 * eps / std::sin(th);
            const real_t tolVal = 1e3 * eps / (std::sin(th)*std::sin(th));
            gsInfo<<"PStrain_FiniteNoZero theta="<<th<<": values "<<res.first.transpose()
                  <<", first(2)-0.123="<<res.first(2)-0.123<<", |col2-n|="<<(nn-n).cwiseAbs().maxCoeff()
                  <<", |col0.n|="<<std::abs(v0.dot(n))<<", |col1.n|="<<std::abs(v1.dot(n))
                  <<" (tolDir "<<tolDir<<", tolVal "<<tolVal<<")\n";
            CHECK(res.first(2) == 0.123);
            CHECK((nn - n).cwiseAbs().maxCoeff() <= tolDir);
            CHECK(std::abs(v0.dot(n)) <= tolDir);
            CHECK(std::abs(v1.dot(n)) <= tolDir);
            CHECK(std::abs(res.first(0)) <= tolVal);
            CHECK(std::abs(res.first(1) - 1.0) <= tolVal);
        }
    }

    /// Preservation check of the regular path: eval3D_pstrain against the principal strains of
    /// eps = (F^T F - I)/2 for an affine deformation
    template<short_t d>
    static void PGuard_values()
    {
        PSDirLin_Scene<d> sc;
        PGuard_Setup<d> su(sc.ref);

        gsMatrix<> u(2,3);
        u << 0.2, 0.5, 0.8,
             0.3, 0.6, 0.1;
        const real_t states[2][3] = { {1.10, 0.95, 0.12}, {1.2, 1.0, 0.0} };
        const char * names[2] = { "generic", "uniaxial" };

        for (index_t s = 0; s != 2; ++s)
        {
            const real_t lam = states[s][0], mu = states[s][1], sh = states[s][2];
            sc.setState(lam,mu,sh);
            su.mm->setDeformed(&sc.def);

            gsMatrix<> C(2,2);
            C << lam*lam, lam*sh,
                 lam*sh,  mu*mu + sh*sh;
            const gsMatrix<> eps2 = 0.5 * (C - gsMatrix<>::Identity(2,2));
            gsEigen::SelfAdjointEigenSolver< gsMatrix<>::Base > es;
            es.compute(eps2);
            const gsVector<> ev = es.eigenvalues();

            real_t wVal = 0, wOut = 0, wUni = 0;
            for (index_t k = 0; k != u.cols(); ++k)
            {
                const gsMatrix<> uk = u.col(k);
                gsMatrix<> zz(1,1);
                zz(0,0) = 0.0;
                const gsMatrix<> P = su.mm->eval3D_pstrain(0,uk,zz);
                CHECK(P.rows()==3 && P.cols()==1);
                if (P.rows()!=3 || P.cols()!=1)
                    return;
                CHECK(P.allFinite());
                const real_t scale = ev.cwiseAbs().maxCoeff();
                wVal = std::max(wVal, std::abs(P(0,0) - ev(0)) / scale);
                wVal = std::max(wVal, std::abs(P(1,0) - ev(1)) / scale);
                wOut = std::max(wOut, std::abs(P(2,0)));
                if (s==1)
                    wUni = std::max(wUni, std::abs(P(0,0)) / P(1,0));
            }
            gsInfo<<"PStrain_Values d="<<d<<" "<<names[s]<<": reference "<<ev.transpose()
                  <<", max rel value err="<<wVal<<", max |P(2)|="<<wOut<<", uniaxial |P(0)|/P(1)="<<wUni<<"\n";
            CHECK(wVal <= 1e-10);
            CHECK(wOut == 0);
            if (s==1)
                CHECK(wUni <= 1e-12);
        }
    }

    TEST(PStrain_ZeroPrincipal)
    {
        PGuard_zeroPrincipal<2>();
        PGuard_zeroPrincipal<3>();
    }

    TEST(PStrain_NonFinite)
    {
        PGuard_nonFinite<2>();
        PGuard_nonFinite<3>();
    }

    TEST(PStrain_FiniteNoZero)
    {
        PGuard_finiteNoZero();
    }

    TEST(PStrain_Values)
    {
        PGuard_values<2>();
        PGuard_values<3>();
    }

    TEST(PStress_DegenerateFrame)
    {
        gsMatrix<> S = gsMatrix<>::Zero(3,3);
        S(0,0) = 2.0;
        S(1,1) = -1.0;
        S(0,1) = S(1,0) = 0.3;
        gsMatrix<> pt(2,1), ptReg(2,1);
        pt    << 0.5, 1.0;
        ptReg << 0.5, 0.5;
        const real_t delta2 = std::ldexp(1.0,-60);
        const real_t delta3 = std::ldexp(1.0,-50);
        gsMatrix<> Q2 = gsMatrix<>::Zero(3,2);
        Q2(0,0) = 1.0;
        Q2(1,1) = 1.0;
        gsVector<real_t,3> e3;
        e3 << 0.0, 0.0, 1.0;
        gsMatrix<> Jc;

        gsInfo<<"PStress_DegenerateFrame\n";
        {   // (i) 2D exact collapse: finite fallback with N = e3
            const gsMultiPatch<> mp = PGuard_collapsed<2>(true,0.0);
            const PGuard_Pair res = PGuard_evalStress<2>(mp,pt,S,Jc);
            PGuard_checkCollapse(Jc,true,"(i) 2D exact");
            PGuard_checkFrame("(i) 2D exact",res,PGuard_T(S,Jc),Q2,1e-12,&e3);
        }
        {   // (ii) 2D near collapse: the finite fallback equals the exact principal decomposition
            const gsMultiPatch<> mp = PGuard_collapsed<2>(false,delta2);
            const PGuard_Pair res = PGuard_evalStress<2>(mp,pt,S,Jc);
            PGuard_checkCollapse(Jc,false,"(ii) 2D near");
            PGuard_checkFrame("(ii) 2D near",res,PGuard_T(S,Jc),Q2,1e-12,&e3);
        }
        {   // (iii) 3D exact collapse: NaN (contract pin: the unfixed library gives incidental NaN)
            const gsMultiPatch<> mp = PGuard_collapsed<3>(true,0.0);
            const PGuard_Pair res = PGuard_evalStress<3>(mp,pt,S,Jc);
            PGuard_checkCollapse(Jc,true,"(iii) 3D exact");
            PGuard_checkAllNaN("(iii) 3D exact",res);
        }
        {   // (iv) 3D near collapse: NaN
            const gsMultiPatch<> mp = PGuard_collapsed<3>(false,delta3);
            const PGuard_Pair res = PGuard_evalStress<3>(mp,pt,S,Jc);
            PGuard_checkCollapse(Jc,false,"(iv) 3D near");
            PGuard_checkAllNaN("(iv) 3D near",res);
        }
        {   // (v) 3D regular anisotropic control (preservation check of the regular path)
            const gsMultiPatch<> mp = PGuard_anisotropic<3>();
            const PGuard_Pair res = PGuard_evalStress<3>(mp,ptReg,S,Jc);
            const gsVector<real_t,3> G1 = Jc.col(0), G2 = Jc.col(1);
            gsVector<real_t,3> n = G1.cross(G2);
            n.normalize();
            const gsVector<real_t,3> e1 = G1.normalized(), e2 = n.cross(e1);
            gsMatrix<> Q(3,2);
            Q.col(0) = e1;
            Q.col(1) = e2;
            PGuard_checkFrame("(v) 3D regular anisotropic",res,PGuard_T(S,Jc),Q,1e-10,nullptr);
        }
        {   // (vi) non-finite input on the regular frame of (v): NaN
            const gsMultiPatch<> mp = PGuard_anisotropic<3>();
            gsMatrix<> Sbad = S;
            Sbad(0,0) = std::numeric_limits<real_t>::quiet_NaN();
            const PGuard_Pair res = PGuard_evalStress<3>(mp,ptReg,Sbad,Jc);
            PGuard_checkAllNaN("(vi) 3D regular, S(0,0)=NaN",res);
        }
    }

    TEST(PStrain_DegenerateFrame)
    {
        gsMatrix<> E = gsMatrix<>::Zero(3,3);
        E(0,0) = 0.2;
        E(1,1) = -0.1;
        E(0,1) = E(1,0) = 0.05;
        E(2,2) = 0.123;
        gsMatrix<> pt(2,1), ptReg(2,1);
        pt    << 0.5, 1.0;
        ptReg << 0.5, 0.5;
        const real_t delta2 = std::ldexp(1.0,-60);
        const real_t delta3 = std::ldexp(1.0,-50);
        gsMatrix<> Jc;

        gsInfo<<"PStrain_DegenerateFrame\n";
        for (index_t c = 0; c != 4; ++c)
        {
            const bool is2D = (c < 2), exact = (c % 2 == 0);
            const char * label = (c==0 ? "(i) 2D exact" : c==1 ? "(ii) 2D near"
                                : c==2 ? "(iii) 3D exact" : "(iv) 3D near");
            gsMultiPatch<> mp = is2D ? PGuard_collapsed<2>(exact,delta2) : PGuard_collapsed<3>(exact,delta3);
            const PGuard_Pair res = is2D ? PGuard_evalStrain<2>(mp,pt,E,Jc) : PGuard_evalStrain<3>(mp,pt,E,Jc);
            PGuard_checkCollapse(Jc,exact,label);
            if (!PGuard_shape(res,3))
                continue;
            gsInfo<<"  "<<label<<": values "<<res.first.transpose()<<"\n";
            CHECK(res.first(2) == 0.123);
            CHECK(std::isnan(res.first(0)));
            CHECK(std::isnan(res.first(1)));
            CHECK(PGuard_allNaN(res.second.leftCols(2)));
            if (is2D)
                CHECK(res.second(0,2)==0 && res.second(1,2)==0 && res.second(2,2)==1);
        }
        {   // (v) 3D regular anisotropic control
            const gsMultiPatch<> mp = PGuard_anisotropic<3>();
            const PGuard_Pair res = PGuard_evalStrain<3>(mp,ptReg,E,Jc);
            if (!PGuard_shape(res,3))
                return;
            const gsVector<real_t,3> G1 = Jc.col(0), G2 = Jc.col(1);
            gsVector<real_t,3> n = G1.cross(G2);
            n.normalize();
            const gsVector<real_t,3> v0 = res.second.col(0), v1 = res.second.col(1);
            gsInfo<<"  (v) 3D regular anisotropic: values "<<res.first.transpose()<<", |col0.n|="
                  <<std::abs(v0.dot(n))<<", |col1.n|="<<std::abs(v1.dot(n))<<"\n";
            CHECK(res.first(2) == 0.123);
            CHECK(res.first.allFinite());
            CHECK(res.second.allFinite());
            CHECK(std::abs(v0.dot(n)) <= 1e-12);
            CHECK(std::abs(v1.dot(n)) <= 1e-12);
        }
    }

    // The energy of the quadratic law is W = 1/2 E.S in engineering Voigt form
    // (E3 S3 = 2 E12 S12), W(2E) = 4 W(E) and W(0) = 0.
    TEST(Energy_Linear_HalfES)
    {
        gsVector<real_t> pt(2);
        pt << 0.25, 0.25;

        real_t worstHalf = 0, worstQuad = 0, worstZero = 0;
        for (index_t skew = 0; skew != 2; ++skew)
        {
            EnergyLinear law(skew==1);
            const gsMatrix<real_t> G = energyLinMetric(law.mp,pt);
            CHECK_CLOSE(skew==1 ? 0.4 : 0.0, G(0,1), 1e-12);
            for (index_t s = 0; s != 3; ++s)
            {
                const gsVector<real_t> E = energyLinStrain(s);
                const gsMatrix<real_t> Cm = energyLinCm(G,E);
                const gsMatrix<real_t> S = law.materialMatrix->eval3D_vector_C(Cm,0,pt,0.0,MaterialOutput::Generic);
                CHECK_EQUAL(3,S.rows());
                if (S.rows()!=3)
                    continue;
                CHECK(S.cwiseAbs().maxCoeff() > 0);

                const real_t W = energyLinAt(*law.materialMatrix,Cm,pt);
                const real_t WES = 0.5 * ( E(0)*S(0,0) + E(1)*S(1,0) + E(2)*S(2,0) );
                const real_t relHalf = math::abs(W-WES) / math::abs(WES);
                CHECK(WES > 0);
                CHECK(relHalf <= 1e-12);

                const real_t W2 = energyLinAt(*law.materialMatrix,energyLinCm(G,2.0*E),pt);
                const real_t relQuad = math::abs(W2-4*W) / math::abs(4*W);
                CHECK(relQuad <= 1e-12);

                worstHalf = std::max(worstHalf,relHalf);
                worstQuad = std::max(worstQuad,relQuad);
            }
            const gsMatrix<real_t> Cm0 = energyLinCm(G,gsVector<real_t>::Zero(3));
            const real_t W0 = energyLinAt(*law.materialMatrix,Cm0,pt);
            const real_t cmax = energyLinMaxC(*law.materialMatrix,Cm0,pt);
            CHECK(cmax > 0);
            CHECK(math::abs(W0) <= 1e-14*cmax);
            worstZero = std::max(worstZero,math::abs(W0)/cmax);
        }
        gsInfo<<"Energy_Linear_HalfES: worst |W - 1/2 E.S|/|1/2 E.S| = "<<worstHalf
              <<", worst |W(2E) - 4W(E)|/|4W(E)| = "<<worstQuad
              <<", worst |W(0)|/max|C| = "<<worstZero<<"\n";
    }

    // S = dW/dE by central finite differences, with the criterion
    // |S - FD| <= 1e-6 max|S| + 1e-10 max|C|. The shear entry of S must be at
    // least 1e-3 max|S|, i.e. a wrong shear factor would miss the criterion by
    // three orders of magnitude (on the sheared reference the shear stress of
    // state 0 is only 2e-3 max|S| because the contravariant metric couples it
    // to the normal stresses).
    TEST(Energy_Linear_FD)
    {
        gsVector<real_t> pt(2);
        pt << 0.25, 0.25;
        const real_t h = 1e-6;

        real_t worstRel = 0, worstRatio = 0, worstShear = 1;
        for (index_t skew = 0; skew != 2; ++skew)
        {
            EnergyLinear law(skew==1);
            const gsMatrix<real_t> G = energyLinMetric(law.mp,pt);
            for (index_t s = 0; s != 3; ++s)
            {
                const gsVector<real_t> E = energyLinStrain(s);
                const gsMatrix<real_t> Cm = energyLinCm(G,E);
                CHECK(Cm(0,0) > 0 && Cm(0,0)*Cm(1,0) - Cm(2,0)*Cm(2,0) > 0);
                const gsMatrix<real_t> S = law.materialMatrix->eval3D_vector_C(Cm,0,pt,0.0,MaterialOutput::Generic);
                CHECK_EQUAL(3,S.rows());
                if (S.rows()!=3)
                    continue;
                real_t err = 0;
                for (index_t j = 0; j != 3; ++j)
                {
                    gsVector<real_t> Ep = E, Em = E;
                    Ep(j) += h;
                    Em(j) -= h;
                    const real_t fd = ( energyLinAt(*law.materialMatrix,energyLinCm(G,Ep),pt)
                                      - energyLinAt(*law.materialMatrix,energyLinCm(G,Em),pt) ) / (2*h);
                    err = std::max(err,math::abs(S(j,0)-fd));
                }
                const real_t Smax = S.cwiseAbs().maxCoeff();
                const real_t tol = 1e-6*Smax + 1e-10*energyLinMaxC(*law.materialMatrix,Cm,pt);
                if (!(err <= tol) || !(math::abs(S(2,0)) >= 1e-3*Smax))
                    gsInfo<<"  failed at "<<(skew ? "skewed" : "orthonormal")<<" state "<<s<<": |S-FD|="<<err<<", tol="<<tol<<"\n";
                CHECK(err <= tol);
                CHECK(math::abs(S(2,0)) >= 1e-3*Smax);
                worstRel   = std::max(worstRel,err/Smax);
                worstRatio = std::max(worstRatio,err/tol);
                worstShear = std::min(worstShear,math::abs(S(2,0))/Smax);
            }
        }
        gsInfo<<"Energy_Linear_FD: worst |S-FD|/max|S| = "<<worstRel<<", worst |S-FD|/tol = "<<worstRatio
              <<", smallest |S_2|/max|S| = "<<worstShear<<"\n";
    }

    // The geometry-based eval3D_energy packs column j*u.cols()+k for point k and
    // fibre row j, and agrees with eval3D_energy_C at the deformed metric.
    TEST(Energy_Linear_GeometryPath)
    {
        gsMatrix<real_t> u(2,3);
        u << 0.2, 0.55, 0.8,
             0.35, 0.1, 0.7;
        gsMatrix<real_t> z(2,3);
        z.row(0).setZero();
        z.row(1).setConstant(0.03);

        real_t worst = 0;
        for (index_t skew = 0; skew != 2; ++skew)
        {
            EnergyLinear law(skew==1);
            law.deform();

            const gsMatrix<real_t> W = law.materialMatrix->eval3D_energy(0,u,z);
            CHECK_EQUAL(1,W.rows());
            CHECK_EQUAL(u.cols()*z.rows(),W.cols());
            if (W.rows()!=1 || W.cols()!=u.cols()*z.rows())
                continue;
            CHECK(W.allFinite());
            // the energy varies over the points, so a packing error cannot go unnoticed
            CHECK(W.leftCols(u.cols()).maxCoeff() - W.leftCols(u.cols()).minCoeff() >= 1e-2 * W.leftCols(u.cols()).maxCoeff());

            for (index_t k = 0; k != u.cols(); ++k)
            {
                const gsVector<real_t> pt = u.col(k);
                const gsMatrix<real_t> Cm = energyLinCm(energyLinMetric(law.def,pt),gsVector<real_t>::Zero(3));
                const real_t Wc = energyLinAt(*law.materialMatrix,Cm,pt);
                CHECK(Wc > 0);
                const real_t rel = math::abs(W(0,k)-Wc) / math::abs(Wc);
                if (!(rel <= 1e-12))
                    gsInfo<<"  failed at "<<(skew ? "skewed" : "orthonormal")<<", point "<<k<<": W="<<W(0,k)<<", W_C="<<Wc<<"\n";
                CHECK(rel <= 1e-12);
                worst = std::max(worst,rel);
            }
        }
        gsInfo<<"Energy_Linear_GeometryPath: worst relative difference to eval3D_energy_C = "<<worst<<"\n";
    }
}
