/** @file gsMaterialMatrixTFTSplit_test.cpp

    @brief Unit tests of the linearised and exact gsMaterialMatrixTFTSplit.

    The oracles are independent of the decorator's own output:
    - central finite differences of the decorator's prescribed-metric stress (the tangent is
      compared with the transposed finite-difference Jacobian, A_ij = dS_i/dE_j),
    - the wrapped material evaluated on its own (taut and slack identities, regularisation blend),
    - the root finder of gsMaterialMatrixTFT for the wrinkle direction,
    - the Voigt direction vectors n1 = (c^2,s^2,2cs), n2 = (-sc,sc,c^2-s^2) for the wrinkled-stress
      conditions n1.S' = n2.S' = 0, and exact binary strains for the criterion boundaries,
    - the wrapped material evaluated at the returned elastic strain E' = E + gamma n1 (exact mode:
      n1.S_w(E') = n2.S_w(E') = 0),
    - synthetic stress/tangent providers with a known solution, which drive every branch of the exact
      solve through the pointwise core,
    - generalized symmetric eigenvalue solves for the principal-value outputs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>
#include <gsParallel/gsOpenMP.h>

#include <cmath>
#include <iomanip>
#include <limits>

using namespace gismo;

SUITE(gsMaterialMatrixTFTSplit_test)                 // The suite should have the same name as the file
{
    /// Largest entry-wise difference of \a A and \a B relative to the largest entry of the reference \a B
    static real_t SplitRelErr(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
    {
        return (A-B).cwiseAbs().maxCoeff() / B.cwiseAbs().maxCoeff();
    }

    /// Relative comparison of two matrices (reference \a B) with the achieved error printed
#   define SPLIT_CHECK_REL(label,A,B,tol) \
    do { \
        gsInfo << "  " << (label) << ": rel.err = " << SplitRelErr((A),(B)) << " (tol " << (tol) << ")\n"; \
        CHECK((A).allFinite()); \
        CHECK((B).allFinite()); \
        CHECK(((A)-(B)).cwiseAbs().maxCoeff() <= (tol) * (B).cwiseAbs().maxCoeff()); \
    } while(0)

    /// Bitwise equality of two finite matrices of the same size
#   define SPLIT_CHECK_SAME(A,B) \
    do { \
        CHECK_EQUAL((A).rows(),(B).rows()); \
        CHECK_EQUAL((A).cols(),(B).cols()); \
        if ((A).rows()==(B).rows() && (A).cols()==(B).cols()) \
        { \
            CHECK((A).allFinite()); \
            CHECK((B).allFinite()); \
            CHECK(((A)-(B)).cwiseAbs().maxCoeff() == 0); \
        } \
    } while(0)

    /// Absolute difference of two angles modulo pi, folded into [0,pi/2]
    static real_t SplitThetaDiff(const real_t a, const real_t b)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        real_t d = std::fmod(a - b, pi);
        if (d >  pi/2) d -= pi;
        if (d < -pi/2) d += pi;
        return math::abs(d);
    }

    /// 1x1 through-thickness coordinate matrix
    static gsMatrix<real_t> SplitZ(const real_t z)
    {
        gsMatrix<real_t> zmat(1,1);
        zmat(0,0) = z;
        return zmat;
    }

    enum SplitMaterial { SplitMatSvK, SplitMatNHComp, SplitMatNHIncomp, SplitMatMRIncomp };

    static real_t SplitPoisson(const SplitMaterial m)
    {
        return (m==SplitMatSvK || m==SplitMatNHComp) ? 0.45 : 0.5;
    }

    /// Young's modulus from the shear modulus mu = 4.225e5
    static real_t SplitYoung(const SplitMaterial m)
    {
        return 2 * 4.225e5 * (1 + SplitPoisson(m));
    }

    /// Principal stretches of the deformed state and the criterion state they must give (the rotation is passed separately)
    struct SplitStateSpec
    {
        real_t lam1, lam2;
        index_t state;
    };

    static SplitStateSpec SplitWrinkledSpec() { SplitStateSpec s = { 1.2, 0.8 / math::sqrt(1.2), 0 }; return s; }
    static SplitStateSpec SplitTautSpec()     { SplitStateSpec s = { 1.10, 1.05, 1 }; return s; }
    static SplitStateSpec SplitSlackSpec()    { SplitStateSpec s = { 0.95, 0.90, -1 }; return s; }

    /// Bilinear unit-square patch mapped by the 2x2 Jacobian \a J (x = J X)
    static gsMultiPatch<real_t> SplitAffinePatch(const gsMatrix<real_t> & J)
    {
        gsMultiPatch<real_t> mp;
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        mp.patch(0).coefs() = mp.patch(0).coefs() * J.transpose();
        mp.addAutoBoundaries();
        return mp;
    }

    /// Deformed Jacobian R diag(lam1,lam2) R^T Jref: principal stretches lam1, lam2 along the axes rotated by \a alpha
    static gsMatrix<real_t> SplitJdef(const gsMatrix<real_t> & Jref, const real_t lam1, const real_t lam2, const real_t alpha)
    {
        gsMatrix<real_t> R(2,2), D(2,2);
        R << math::cos(alpha), -math::sin(alpha), math::sin(alpha), math::cos(alpha);
        D << lam1, 0.0, 0.0, lam2;
        return R * D * R.transpose() * Jref;
    }

    /**
     * Owns a wrapped material, its parameter functions, the patches, the decorator under test and
     * the TFT reference, all acting on the same wrapped object. The parameter functions are members
     * because the nonlinear materials keep references to them.
     */
    template<short_t dim, bool linear>
    struct SplitFixture
    {
        typedef gsMaterialMatrixTFTSplit<dim,real_t,linear> SplitType;
        typedef gsMaterialMatrixTFT<dim,real_t,linear>      TFTType;

        gsMultiPatch<real_t> mp, mp_def;
        gsConstantFunction<real_t> t, E, nu, ratio;
        gsMaterialMatrixBase<real_t>::uPtr wrapped;
        memory::unique_ptr<SplitType> split;
        memory::unique_ptr<TFTType>   tft;

        /// Metrics of an affine pair (empty otherwise): reference 2x2, deformed 2x2 and reference (G11,G22,G12)
        gsMatrix<real_t> Gref, Gdef, gori;

        SplitFixture(const SplitMaterial mat, const gsMultiPatch<real_t> & ref, const gsMultiPatch<real_t> & def, const real_t thickness)
        :
        mp(ref), mp_def(def),
        t(thickness,dim), E(SplitYoung(mat),dim), nu(SplitPoisson(mat),dim), ratio(7.0,dim)
        {
            std::vector<gsFunctionSet<real_t>*> parameters(2);
            parameters[0] = &E;
            parameters[1] = &nu;
            if (mat==SplitMatMRIncomp)
                parameters.push_back(&ratio);

            gsOptionList options;
            options.addInt("Material","Material model: (0): SvK | (1): NH | (3): MR",
                           mat==SplitMatSvK ? 0 : (mat==SplitMatMRIncomp ? 3 : 1));
            options.addInt("Implementation","Implementation: (1): Analytical",1);
            if (mat!=SplitMatSvK)
                options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",mat==SplitMatNHComp);
            wrapped = getMaterialMatrix<dim,real_t>(mp,t,parameters,options);

            split = memory::make_unique(new SplitType(wrapped.get()));
            tft   = memory::make_unique(new TFTType(wrapped.get()));
            // The decorators take the geometries explicitly: the constructors copy them from the wrapped
            // material only when the latter already has a deformed geometry
            split->setUndeformed(&mp);
            tft->setUndeformed(&mp);
            split->setDeformed(&mp_def);
            tft->setDeformed(&mp_def);
        }

        /// Affine 2D pair: Jref and the deformed Jacobian R diag(lam1,lam2) R^T Jref
        SplitFixture(const SplitMaterial mat, const gsMatrix<real_t> & Jref, const real_t lam1, const real_t lam2, const real_t alpha, const real_t thickness = 0.01)
        :
        SplitFixture(mat,SplitAffinePatch(Jref),SplitAffinePatch(SplitJdef(Jref,lam1,lam2,alpha)),thickness)
        {
            const gsMatrix<real_t> Jdef = SplitJdef(Jref,lam1,lam2,alpha);
            Gref = Jref.transpose() * Jref;
            Gdef = Jdef.transpose() * Jdef;
            gori.resize(3,1);
            gori << Gref(0,0), Gref(1,1), Gref(0,1);
        }

        /// Prescribed metric (g11,g22,g12) of the deformed state (tensor g12)
        gsMatrix<real_t> CmatDef() const
        {
            gsMatrix<real_t> c(3,1);
            c << Gdef(0,0), Gdef(1,1), Gdef(0,1);
            return c;
        }

        /// Engineering Voigt strain of the deformed state
        gsMatrix<real_t> EengDef() const
        {
            gsMatrix<real_t> e(3,1);
            e << 0.5*(Gdef(0,0)-Gref(0,0)), 0.5*(Gdef(1,1)-Gref(1,1)), Gdef(0,1)-Gref(0,1);
            return e;
        }

        /// Prescribed metric Cmat = 2 E + g_ori for an engineering Voigt strain (the shear entry is halved first)
        gsMatrix<real_t> CmatFromE(const gsMatrix<real_t> & Eeng) const
        {
            gsMatrix<real_t> Et = Eeng;
            Et(2,0) *= 0.5;
            return 2*Et + gori;
        }

        /// Criterion state (1 taut, 0 wrinkled, -1 slack) that eval3D_vector_C classifies at \a Cmat,
        /// reproduced from the wrapped material's own stress and tangent (before any exact-mode reclassification)
        index_t stateAt(const gsMatrix<real_t> & Cmat, const gsVector<real_t> & pt, const real_t z) const
        {
            gsMatrix<real_t> Eeng(3,1);
            Eeng << 0.5*(Cmat(0,0)-Gref(0,0)), 0.5*(Cmat(1,0)-Gref(1,1)), Cmat(2,0)-Gref(0,1);
            const gsMatrix<real_t> Sw = wrapped->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::VectorN);
            const gsMatrix<real_t> Cw = wrapped->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::MatrixA);
            detail_tftsplit::WrinkleOptions<real_t> opt = split->wrinkleOptions();
            opt.exact = false;
            return split->wrinkle(Sw,Cw,Eeng,Gref,opt,false).state;
        }

        /// Displayed state of the geometry path at a single point
        index_t displayedState(const gsVector<real_t> & pt, const real_t z) const
        {
            return static_cast<index_t>(split->eval3D_tensionfield(0,pt,SplitZ(z),MaterialOutput::Generic)(0,0));
        }
    };

    /// Prescribed-metric stress S_eff of the decorator as a function of the engineering Voigt strain.
    /// For every evaluated column the criterion state is classified from the wrapped material's own
    /// S and C, so that the finite-difference stencil is known to stay in the intended state.
    template<short_t dim, bool linear>
    class SplitStressStencil : public gsFunction<real_t>
    {
    public:
        SplitStressStencil(const SplitFixture<dim,linear> & fx, const gsVector<real_t> & pt, const real_t z, const index_t expected)
        :
        m_fx(fx), m_pt(pt), m_z(z), m_expected(expected), m_columns(0), m_allExpected(true)
        {
        }

        void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const
        {
            result.resize(3,u.cols());
            for (index_t k=0; k!=u.cols(); ++k)
            {
                const gsMatrix<real_t> Eeng = u.col(k);
                const gsMatrix<real_t> Cmat = m_fx.CmatFromE(Eeng);
                result.col(k) = m_fx.split->eval3D_vector_C(Cmat,0,m_pt,m_z,MaterialOutput::Generic);
                m_allExpected = m_allExpected && (m_fx.stateAt(Cmat,m_pt,m_z) == m_expected);
                ++m_columns;
            }
        }

        short_t domainDim() const {return 3;}
        short_t targetDim() const {return 3;}

        index_t columns() const { return m_columns; }
        bool allExpected() const { return m_allExpected; }

    private:
        const SplitFixture<dim,linear> & m_fx;
        const gsVector<real_t> m_pt;
        const real_t m_z;
        const index_t m_expected;
        mutable index_t m_columns;
        mutable bool m_allExpected;
    };

    /// Transposed finite-difference Jacobian J_ij = dS_i/dE_j of the stencil function at \a Eeng
    /// (deriv_into stores dS_c/dE_j at j+3c)
    template<short_t dim, bool linear>
    static gsMatrix<real_t> SplitFDJacobian(const SplitFixture<dim,linear> & fx, const gsVector<real_t> & pt, const real_t z,
                                            const index_t expected, const gsMatrix<real_t> & Eeng,
                                            index_t & columns, bool & allExpected)
    {
        SplitStressStencil<dim,linear> stencil(fx,pt,z,expected);
        gsMatrix<real_t> res;
        stencil.deriv_into(Eeng,res);
        columns = stencil.columns();
        allExpected = stencil.allExpected();
        return res.reshape(3,3).transpose();
    }

    /// Isotropic plane-stress SvK tangent (E = 1, Poisson ratio nu) in engineering Voigt form for a unit reference metric
    static gsMatrix<real_t> SplitIsoC(const real_t nu = 0.3)
    {
        gsMatrix<real_t> C(3,3);
        C << 1, nu, 0,
             nu, 1, 0,
             0,  0, (1-nu)/2;
        return C / (1-nu*nu);
    }

    /// Inputs of the pointwise core for which the isotropy guard rejects the closed-form direction:
    /// C = T^T C0 T is a rotated orthotropic tangent, S = C E.
    static void SplitAnisoInput(gsMatrix<real_t> & S, gsMatrix<real_t> & C, gsMatrix<real_t> & E)
    {
        const real_t beta = 0.4;
        const real_t c = math::cos(beta), s = math::sin(beta);
        gsMatrix<real_t> C0(3,3), T(3,3);
        C0 << 10, 1, 0,
              1,  2, 0,
              0,  0, 1.5;
        T << c*c,    s*s,    c*s,
             s*s,    c*c,   -c*s,
            -2*c*s,  2*c*s,  c*c-s*s;
        C = T.transpose() * C0 * T;
        E.resize(3,1);
        E << 0.02, -0.03, 0.01;
        S = C * E;
    }

    /// Reference Jacobian of the skewed affine pair (g12 != 0 in the reference metric)
    static gsMatrix<real_t> SplitSkewedJref()
    {
        gsMatrix<real_t> Jref(2,2);
        Jref << 1.0, 0.3, 0.0, 0.8;
        return Jref;
    }

    static gsMatrix<real_t> SplitIdentity2()
    {
        gsMatrix<real_t> I(2,2);
        I.setIdentity();
        return I;
    }

    /// Single-point evaluation point of the affine tests
    static gsVector<real_t> SplitPoint()
    {
        gsVector<real_t> pt(2);
        pt.setConstant(0.25);
        return pt;
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// FD of the decorator's prescribed-metric stress against its tangent, over the three criterion
    /// states and eta in {0,1e-3}. The slack state at eta = 0 is identically zero and is checked as such.
    template<bool linear>
    static void SplitFDCheck(const SplitMaterial mat, const real_t tol)
    {
        const SplitStateSpec specs[3] = { SplitWrinkledSpec(), SplitTautSpec(), SplitSlackSpec() };
        const real_t etas[2] = { 0.0, 1e-3 };
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const char * names[3] = { "wrinkled", "taut", "slack" };

        for (index_t s = 0; s!=3; ++s)
            for (index_t e = 0; e!=2; ++e)
            {
                gsInfo << "state " << names[s] << ", eta " << etas[e] << "\n";
                SplitFixture<2,linear> fx(mat,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,static_cast<real_t>(EIGEN_PI)/6);
                fx.split->options().setReal("Eta",etas[e]);
                fx.split->resetCounters();

                const gsMatrix<real_t> Cmat = fx.CmatDef();
                const gsMatrix<real_t> Eeng = fx.EengDef();
                CHECK_EQUAL(specs[s].state, fx.stateAt(Cmat,pt,z));

                const gsMatrix<real_t> Seff = fx.split->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::Generic);
                const gsMatrix<real_t> Ceff = fx.split->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::Generic).reshape(3,3);
                if (specs[s].state==-1 && etas[e]==0)
                {
                    const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::VectorN);
                    CHECK(Sw.allFinite());
                    CHECK(Seff.allFinite());
                    CHECK(Ceff.allFinite());
                    CHECK(Sw.cwiseAbs().maxCoeff() > 0);
                    CHECK(Seff.cwiseAbs().maxCoeff() == 0);
                    CHECK(Ceff.cwiseAbs().maxCoeff() == 0);
                }
                else
                {
                    index_t columns = 0;
                    bool allExpected = false;
                    const gsMatrix<real_t> J_FD = SplitFDJacobian(fx,pt,z,specs[s].state,Eeng,columns,allExpected);
                    CHECK_EQUAL(12,columns);
                    CHECK(allExpected);
                    SPLIT_CHECK_REL("J_C vs FD",Ceff,J_FD,tol);
                    if (specs[s].state==-1)
                        CHECK(Ceff.cwiseAbs().maxCoeff() > 0);
                }
                CHECK_EQUAL(0,fx.split->numThetaFallbacks());
                CHECK_EQUAL(0,fx.split->numTFTFailures());
            }
    }

    TEST(FD_C_SvK)
    {
        SplitFDCheck<true>(SplitMatSvK,1e-6);
    }

    TEST(FD_C_NH_Comp)
    {
        SplitFDCheck<false>(SplitMatNHComp,1e-4);
    }

    TEST(FD_C_NH_Incomp)
    {
        SplitFDCheck<false>(SplitMatNHIncomp,1e-4);
    }

    TEST(FD_C_MR_Incomp)
    {
        SplitFDCheck<false>(SplitMatMRIncomp,1e-4);

        // The Mooney-Rivlin stress differs from the Neo-Hookean one at the same metric, so the material id took effect
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,false> mr(SplitMatMRIncomp,SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6);
        SplitFixture<2,false> nh(SplitMatNHIncomp,SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6);
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> Smr = mr.wrapped->eval3D_vector_C(mr.CmatDef(),0,pt,0,MaterialOutput::VectorN);
        const gsMatrix<real_t> Snh = nh.wrapped->eval3D_vector_C(nh.CmatDef(),0,pt,0,MaterialOutput::VectorN);
        gsInfo << "  MR vs NH stress: rel.diff = " << SplitRelErr(Smr,Snh) << "\n";
        CHECK(SplitRelErr(Smr,Snh) > 1e-3);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Orientation of the tangent: the finite-difference Jacobian of the decorator's own stress matches
    /// A_ij = dS_i/dE_j only after transposing deriv_into's layout, and the tangent is not symmetric.
    TEST(Orientation_NH_Incomp)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),w.lam1,w.lam2,pi/6);
        fx.split->resetCounters();
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const gsMatrix<real_t> Cmat = fx.CmatDef();
        const gsMatrix<real_t> Eeng = fx.EengDef();
        CHECK_EQUAL(0,fx.stateAt(Cmat,pt,z));

        SplitStressStencil<2,false> stencil(fx,pt,z,0);
        gsMatrix<real_t> res;
        stencil.deriv_into(Eeng,res);
        CHECK_EQUAL(12,stencil.columns());
        CHECK(stencil.allExpected());

        const gsMatrix<real_t> J_FD = res.reshape(3,3).transpose();
        const gsMatrix<real_t> J_untransposed = res.reshape(3,3);
        const gsMatrix<real_t> J = fx.split->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::Generic).reshape(3,3);

        SPLIT_CHECK_REL("J_C vs transposed FD",J,J_FD,1e-4);
        CHECK((J_FD - J_FD.transpose()).cwiseAbs().maxCoeff() >= 100 * 1e-4 * J_FD.cwiseAbs().maxCoeff());
        gsInfo << "  untransposed FD vs J_C: rel.err = " << SplitRelErr(J_untransposed,J) << "\n";
        CHECK(SplitRelErr(J_untransposed,J) > 1e-2);
        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Geometry path against prescribed-metric path at the same point, fibre and metric (eta = 0)
    template<bool linear>
    static void SplitPathsCheck(const SplitMaterial mat, const real_t tolMatC, const real_t tolFD)
    {
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,linear> fx(mat,SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6);
        fx.split->resetCounters();
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const gsMatrix<real_t> zmat = SplitZ(z);
        const gsMatrix<real_t> Cmat = fx.CmatDef();
        const gsMatrix<real_t> Eeng = fx.EengDef();
        CHECK_EQUAL(0,fx.stateAt(Cmat,pt,z));

        const gsMatrix<real_t> Sgeo = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> SC   = fx.split->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::Generic);
        SPLIT_CHECK_REL("S geometry vs S_C",Sgeo,SC,1e-10);

        const gsMatrix<real_t> Ageo = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        const gsMatrix<real_t> AC   = fx.split->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::Generic);
        SPLIT_CHECK_REL("A geometry vs A_C",Ageo,AC,tolMatC);

        index_t columns = 0;
        bool allExpected = false;
        const gsMatrix<real_t> J_FD = SplitFDJacobian(fx,pt,z,0,Eeng,columns,allExpected);
        CHECK_EQUAL(12,columns);
        CHECK(allExpected);
        const gsMatrix<real_t> Ageo33 = Ageo.reshape(3,3);
        SPLIT_CHECK_REL("A geometry vs FD of S_C",Ageo33,J_FD,tolFD);

        CHECK_EQUAL(0,fx.displayedState(pt,z));
        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
    }

    TEST(Paths_SvK)
    {
        SplitPathsCheck<true>(SplitMatSvK,1e-10,1e-5);
    }

    TEST(Paths_NH_Incomp)
    {
        SplitPathsCheck<false>(SplitMatNHIncomp,1e-4,1e-6);
    }

    TEST(Paths_NH_Comp)
    {
        SplitPathsCheck<false>(SplitMatNHComp,1e-4,1e-5);
    }

    TEST(Paths_MR_Incomp)
    {
        SplitPathsCheck<false>(SplitMatMRIncomp,1e-4,1e-6);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Equivalence with gsMaterialMatrixTFT (eta = 0, SlackMultiplier = 0.5 on both) on the same wrapped
    /// material: closed-form wrinkle direction against TFT's root, then S and C at wrinkled, taut and slack points.
    template<bool linear>
    static void SplitEquivCheck(const SplitMaterial mat, const real_t tolC)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const gsMatrix<real_t> zmat = SplitZ(z);
        const SplitStateSpec specs[3] = { SplitWrinkledSpec(), SplitTautSpec(), SplitSlackSpec() };

        for (index_t s = 0; s!=3; ++s)
        {
            SplitFixture<2,linear> fx(mat,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,pi/6);
            fx.split->options().setReal("Eta",0.0);
            fx.split->options().setReal("SlackMultiplier",0.5);
            fx.tft->options().setReal("SlackMultiplier",0.5);
            fx.split->resetCounters();
            CHECK_EQUAL(specs[s].state,fx.stateAt(fx.CmatDef(),pt,z));

            const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const gsMatrix<real_t> Ew = fx.wrapped->eval3D_strain(0,pt,zmat);

            const gsMatrix<real_t> Ssplit = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Asplit = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const gsMatrix<real_t> Stft   = fx.tft->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Atft   = fx.tft->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);

            if (specs[s].state==0)
            {
                const real_t thetaSplit = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
                const real_t thetaTFT   = fx.tft->eval_theta(Cw,Sw,Ew)(0,0);
                gsInfo << "  theta split = " << thetaSplit << ", theta TFT = " << thetaTFT
                       << ", diff mod pi = " << SplitThetaDiff(thetaSplit,thetaTFT) << "\n";
                CHECK(SplitThetaDiff(thetaSplit,thetaTFT) <= 1e-8);
                CHECK(thetaSplit >= 0 && thetaSplit < pi);
                if (SplitThetaDiff(thetaSplit,thetaTFT) <= 1e-8)
                {
                    SPLIT_CHECK_REL("S split vs TFT",Ssplit,Stft,1e-10);
                    SPLIT_CHECK_REL("C split vs TFT",Asplit,Atft,tolC);
                }
            }
            else
            {
                SPLIT_CHECK_REL("S split vs TFT",Ssplit,Stft,1e-13);
                SPLIT_CHECK_REL("C split vs TFT",Asplit,Atft,1e-13);
                if (specs[s].state==1)
                {
                    SPLIT_CHECK_REL("S split vs wrapped",Ssplit,Sw,1e-13);
                    SPLIT_CHECK_REL("C split vs wrapped",Asplit,Cw,1e-13);
                }
                else
                {
                    const gsMatrix<real_t> halfS = 0.5 * Sw;
                    SPLIT_CHECK_REL("S split vs 0.5 S wrapped",Ssplit,halfS,1e-13);
                }
            }
            CHECK_EQUAL(0,fx.split->numThetaFallbacks());
            CHECK_EQUAL(0,fx.split->numTFTFailures());
        }
    }

    TEST(EquivTFT_SvK_skewed)
    {
        SplitEquivCheck<true>(SplitMatSvK,1e-10);
    }

    TEST(EquivTFT_NH_Incomp)
    {
        SplitEquivCheck<false>(SplitMatNHIncomp,1e-7);
    }

    TEST(EquivTFT_NH_Comp)
    {
        SplitEquivCheck<false>(SplitMatNHComp,1e-7);
    }

    TEST(EquivTFT_MR_Incomp)
    {
        SplitEquivCheck<false>(SplitMatMRIncomp,1e-7);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// The wrinkled stress has one vanishing principal value and satisfies n1.S' = n2.S' = 0 for the
    /// Voigt direction vectors of the reported angle (eta = 0)
    template<bool linear>
    static void SplitWrinkledStressCheck(const SplitMaterial mat)
    {
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,linear> fx(mat,SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6);
        fx.split->resetCounters();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        CHECK_EQUAL(0,fx.displayedState(pt,0));

        const gsMatrix<real_t> p = fx.split->eval3D_pstress(0,pt,zmat,MaterialOutput::Generic);
        CHECK_EQUAL(2,p.rows());
        gsInfo << "  principal stresses of S_eff: " << p(0,0) << ", " << p(1,0) << "\n";
        const real_t pmax = p.cwiseAbs().maxCoeff();
        const index_t zero  = math::abs(p(0,0)) < math::abs(p(1,0)) ? 0 : 1;
        const index_t other = 1 - zero;
        CHECK(math::abs(p(zero,0)) <= 1e-8 * pmax);
        CHECK(p(other,0) > 0);

        const real_t theta = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
        const gsMatrix<real_t> S = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const real_t c = math::cos(theta), s = math::sin(theta);
        gsVector<real_t,3> n1, n2;
        n1 << c*c, s*s, 2*c*s;
        n2 << -s*c, s*c, c*c-s*s;
        const real_t n1S = n1.dot(S.col(0)), n2S = n2.dot(S.col(0));
        gsInfo << "  n1.S = " << n1S << ", n2.S = " << n2S << ", |S| = " << S.norm() << "\n";
        CHECK(math::abs(n1S) <= 1e-8 * S.norm());
        CHECK(math::abs(n2S) <= 1e-8 * S.norm());
        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
    }

    TEST(Wrinkled_Smin_zero)
    {
        SplitWrinkledStressCheck<true>(SplitMatSvK);
        SplitWrinkledStressCheck<false>(SplitMatNHIncomp);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Eta and SlackMultiplier limits on the geometry path against the wrapped material evaluated on its own
    template<bool linear>
    static void SplitEtaCheck(const SplitMaterial mat)
    {
        const real_t eta = 1e-3;
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);

        // taut: exact copy of the wrapped S and C
        {
            const SplitStateSpec sp = SplitTautSpec();
            SplitFixture<2,linear> fx(mat,SplitSkewedJref(),sp.lam1,sp.lam2,pi/6);
            fx.split->options().setReal("Eta",eta);
            const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const gsMatrix<real_t> S = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> C = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            CHECK_EQUAL(1,fx.displayedState(pt,0));
            SPLIT_CHECK_REL("taut S",S,Sw,1e-13);
            SPLIT_CHECK_REL("taut C",C,Cw,1e-13);
        }

        // slack: eta S and eta C, and SlackMultiplier S and SlackMultiplier C when it is positive
        {
            const SplitStateSpec sp = SplitSlackSpec();
            SplitFixture<2,linear> fx(mat,SplitSkewedJref(),sp.lam1,sp.lam2,pi/6);
            fx.split->options().setReal("Eta",eta);
            fx.split->options().setReal("SlackMultiplier",0.0);
            const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            CHECK_EQUAL(-1,fx.displayedState(pt,0));

            gsMatrix<real_t> S = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            gsMatrix<real_t> C = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            gsMatrix<real_t> Sexp = eta * Sw, Cexp = eta * Cw;
            SPLIT_CHECK_REL("slack S, SlackMultiplier 0",S,Sexp,1e-13);
            SPLIT_CHECK_REL("slack C, SlackMultiplier 0",C,Cexp,1e-13);

            fx.split->options().setReal("SlackMultiplier",0.25);
            S = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            C = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            Sexp = 0.25 * Sw;
            Cexp = 0.25 * Cw;
            SPLIT_CHECK_REL("slack S, SlackMultiplier 0.25",S,Sexp,1e-13);
            SPLIT_CHECK_REL("slack C, SlackMultiplier 0.25",C,Cexp,1e-13);
        }

        // wrinkled: (1-eta) X(eta=0) + eta X_wrapped
        {
            const SplitStateSpec sp = SplitWrinkledSpec();
            SplitFixture<2,linear> fx(mat,SplitSkewedJref(),sp.lam1,sp.lam2,pi/6);
            const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            fx.split->options().setReal("Eta",0.0);
            CHECK_EQUAL(0,fx.displayedState(pt,0));
            const gsMatrix<real_t> S0 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> C0 = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            fx.split->options().setReal("Eta",eta);
            const gsMatrix<real_t> S = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> C = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const gsMatrix<real_t> Sexp = (1-eta) * S0 + eta * Sw;
            const gsMatrix<real_t> Cexp = (1-eta) * C0 + eta * Cw;
            SPLIT_CHECK_REL("wrinkled S blend",S,Sexp,1e-13);
            SPLIT_CHECK_REL("wrinkled C blend",C,Cexp,1e-13);
            // the blend is not the unregularised result
            CHECK(SplitRelErr(S,S0) > 1e-6);
        }
    }

    TEST(EtaLimits)
    {
        SplitEtaCheck<false>(SplitMatNHComp);
        SplitEtaCheck<true>(SplitMatSvK);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// The tension-field criterion switches with the TensionField option. A stress-taut point whose
    /// strain has a negative principal value is taut for TensionField 1 and 2 and wrinkled for 0.
    TEST(TensionField_switch)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);

        // Pointwise core, unit metric, isotropic C
        {
            SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),1.0,1.0,0.0);
            gsMatrix<real_t> S(3,1), E(3,1);
            S << 2.0, 1.0, 0.0;
            E << 0.01, -0.001, 0.0;
            const gsMatrix<real_t> C = SplitIsoC();
            const gsMatrix<real_t> I2 = SplitIdentity2();
            detail_tftsplit::WrinkleOptions<real_t> opt;
            opt.tensionField = 2; CHECK_EQUAL(1,fx.split->wrinkle(S,C,E,I2,opt,false).state);
            opt.tensionField = 0; CHECK_EQUAL(0,fx.split->wrinkle(S,C,E,I2,opt,false).state);
            opt.tensionField = 1; CHECK_EQUAL(1,fx.split->wrinkle(S,C,E,I2,opt,false).state);
        }

        // Through the class: SvK on a unit-square reference with Jdef = diag(sqrt(1.02), sqrt(0.998))
        {
            SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),math::sqrt(1.02),math::sqrt(0.998),0.0);
            fx.split->resetCounters();
            const gsVector<real_t> pt = SplitPoint();
            const gsMatrix<real_t> zmat = SplitZ(0);
            const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const detail_tftsplit::WrinkleResult<real_t> r =
                fx.split->wrinkle(Sw,Cw,fx.EengDef(),SplitIdentity2(),fx.split->wrinkleOptions(),false);
            gsInfo << "  E principal values (" << r.Emin << ", " << r.Emax << "), S principal values ("
                   << r.Smin << ", " << r.Smax << ")\n";
            CHECK(r.Emin < 0 && r.Emax > 0 && r.Smin > 0);

            // TensionField = 2 (default): taut
            CHECK(std::isnan(fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0)));
            const gsMatrix<real_t> S2 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            SPLIT_CHECK_REL("TF=2 S vs wrapped S",S2,Sw,1e-13);

            // TensionField = 0: wrinkled, with a negative wrinkle amplitude
            gsOptionList o;
            o.addInt("TensionField","",0);
            fx.split->setOptions(o);
            const real_t theta = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
            CHECK(std::isfinite(theta) && theta >= 0 && theta < pi);
            CHECK(fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0) < 0);
            const gsMatrix<real_t> S0 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            gsInfo << "  TF=0 vs TF=2 stress: rel.diff = " << SplitRelErr(S0,S2) << "\n";
            CHECK(SplitRelErr(S0,S2) > 1e-3);

            // Displayed state, with the S2crit threshold removed
            fx.split->options().setReal("S2crit",std::numeric_limits<real_t>::max());
            CHECK_EQUAL(0,fx.displayedState(pt,0));
            fx.split->options().setInt("TensionField",2);
            CHECK_EQUAL(1,fx.displayedState(pt,0));

            CHECK_EQUAL(0,fx.split->numThetaFallbacks());
            CHECK_EQUAL(0,fx.split->numTFTFailures());

            // An unknown TensionField value is rejected (setOptions itself does not validate)
            fx.split->options().setInt("TensionField",3);
            CHECK_THROW(fx.split->wrinkleOptions(), std::runtime_error);
        }
    }

    /// The slack boundary is strict: E_max < 0 with a relative tolerance 1e-12 max|E_i|. Power-of-two strains
    /// make the pencil arithmetic exact.
    TEST(SlackBoundary_strict)
    {
        SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),1.0,1.0,0.0);
        const gsMatrix<real_t> I2 = SplitIdentity2();
        const gsMatrix<real_t> C = SplitIsoC();
        gsMatrix<real_t> S(3,1), E(3,1);
        S << 0.5, -1.0, 0.0;
        detail_tftsplit::WrinkleOptions<real_t> opt;
        opt.tensionField = 2;

        E << 0.0, -std::ldexp(1.0,-6), 0.0;
        detail_tftsplit::WrinkleResult<real_t> r = fx.split->wrinkle(S,C,E,I2,opt,false);
        gsInfo << "  Emax = " << r.Emax << ", Smin = " << r.Smin << "\n";
        CHECK(r.Smin < 0);
        CHECK(r.Emax == 0);
        CHECK_EQUAL(0,r.state);

        E << -std::ldexp(1.0,-20), -std::ldexp(1.0,-6), 0.0;
        r = fx.split->wrinkle(S,C,E,I2,opt,false);
        gsInfo << "  Emax = " << r.Emax << "\n";
        CHECK(r.Emax < 0);
        CHECK_EQUAL(-1,r.state);

        E << -1e-16, -std::ldexp(1.0,-6), 0.0;
        r = fx.split->wrinkle(S,C,E,I2,opt,false);
        gsInfo << "  Emax = " << r.Emax << ", Emin = " << r.Emin << "\n";
        CHECK(r.Emax < 0 && math::abs(r.Emax) <= 1e-12 * math::abs(r.Emin));
        CHECK_EQUAL(0,r.state);
    }

    /// The displayed state of a wrinkled point is 0 iff the unmodified S_min < S2crit (strict), and S2crit
    /// leaves S and C untouched.
    TEST(S2crit_strict)
    {
        SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),1.2,0.8,0.0);
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const gsMatrix<real_t> zmat = SplitZ(z);

        // The core's S_min, reproduced from the same single-point wrapped calls
        const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,z,MaterialOutput::VectorN);
        const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,z,MaterialOutput::MatrixA);
        const real_t Smin = fx.split->wrinkle(Sw,Cw,fx.EengDef(),SplitIdentity2(),fx.split->wrinkleOptions(),false).Smin;
        gsInfo << "  Smin = " << std::setprecision(17) << Smin << std::setprecision(6) << "\n";
        CHECK(Smin < 0);

        const real_t inf = std::numeric_limits<real_t>::infinity();
        CHECK_EQUAL(0,fx.displayedState(pt,z));                                   // default S2crit = 0
        fx.split->options().setReal("S2crit",Smin);
        CHECK_EQUAL(1,fx.displayedState(pt,z));                                   // Smin < Smin fails
        fx.split->options().setReal("S2crit",std::nextafter(Smin,inf));
        CHECK_EQUAL(0,fx.displayedState(pt,z));
        fx.split->options().setReal("S2crit",std::nextafter(Smin,-inf));
        CHECK_EQUAL(1,fx.displayedState(pt,z));
        fx.split->options().setReal("S2crit",Smin + 1e-6*math::abs(Smin));
        CHECK_EQUAL(0,fx.displayedState(pt,z));

        // S2crit is display-only
        const gsMatrix<real_t> S1 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> A1 = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        fx.split->options().setReal("S2crit",0.0);
        const gsMatrix<real_t> S0 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> A0 = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        CHECK(S1.allFinite());
        CHECK(S0.allFinite());
        CHECK(A1.allFinite());
        CHECK(A0.allFinite());
        CHECK((S1-S0).cwiseAbs().maxCoeff() == 0);
        CHECK((A1-A0).cwiseAbs().maxCoeff() == 0);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// An anisotropic tangent makes the closed-form (isotropic) wrinkle direction violate n2.S' = 0, so
    /// the isotropy guard falls back to the root finder of gsMaterialMatrixTFT.
    TEST(Guard_anisotropic_fallback)
    {
        SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),1.0,1.0,0.0);
        fx.split->resetCounters();
        gsMatrix<real_t> S, C, E;
        SplitAnisoInput(S,C,E);
        const gsMatrix<real_t> I2 = SplitIdentity2();

        // Closed-form direction: minor eigenvector of [[S0,S2],[S2,S1]] for the unit metric
        const real_t a = S(0,0), b = S(2,0), d = S(1,0);
        const real_t lmin = 0.5*(a+d) - math::sqrt(0.25*(a-d)*(a-d) + b*b);
        const real_t lmax = 0.5*(a+d) + math::sqrt(0.25*(a-d)*(a-d) + b*b);
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        real_t thetaCF = std::atan2(lmin - a, b);
        if (thetaCF < 0) thetaCF += pi;
        const real_t c = math::cos(thetaCF), s = math::sin(thetaCF);
        gsVector<real_t,3> n1, n2;
        n1 << c*c, s*s, 2*c*s;
        n2 << -s*c, s*c, c*c-s*s;
        const gsVector<real_t,3> Cn1 = C * n1;
        const real_t gam = -n1.dot(S.col(0)) / n1.dot(Cn1);
        const real_t fCF = n2.dot(S.col(0)) + gam * n2.dot(Cn1);
        const real_t tol = 1e-10 * S.norm() + 1e-8 * (lmax - lmin);
        gsInfo << "  |f(theta_cf)| = " << math::abs(fCF) << ", guard tolerance = " << tol << "\n";
        CHECK(math::abs(fCF) > 1e3 * tol);

        detail_tftsplit::WrinkleOptions<real_t> opt;
        opt.tensionField = 2;
        const detail_tftsplit::WrinkleResult<real_t> r = fx.split->wrinkle(S,C,E,I2,opt,false);
        CHECK_EQUAL(0,r.state);
        CHECK(r.guardFallback);
        CHECK_EQUAL(1,fx.split->numThetaFallbacks());
        CHECK(!r.tftFailure);
        CHECK_EQUAL(0,fx.split->numTFTFailures());

        const gsMatrix<real_t> C9 = C.reshape(9,1);
        const real_t thetaTFT = fx.tft->eval_theta(C9,S,E)(0,0);
        gsInfo << "  theta = " << r.theta << ", theta TFT = " << thetaTFT << ", closed form = " << thetaCF
               << ", residual = " << r.residual << "\n";
        CHECK(SplitThetaDiff(r.theta,thetaTFT) <= 1e-10);
        CHECK(r.residual <= tol);
    }

    /// A tangent that is not positive definite along the wrinkle direction has no admissible correction:
    /// the stress is returned unmodified and the failure is counted.
    TEST(Guard_indefinite_C)
    {
        SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),1.0,1.0,0.0);
        fx.split->resetCounters();
        gsMatrix<real_t> S(3,1), E(3,1);
        S << 1.0, -0.5, 0.0;
        E << 0.01, -0.01, 0.0;
        const gsMatrix<real_t> C = -SplitIsoC();
        const gsMatrix<real_t> I2 = SplitIdentity2();
        detail_tftsplit::WrinkleOptions<real_t> opt;
        opt.tensionField = 2;
        opt.eta = 0;
        const detail_tftsplit::WrinkleResult<real_t> r = fx.split->wrinkle(S,C,E,I2,opt,false);
        CHECK_EQUAL(0,r.state);
        CHECK(r.tftFailure);
        CHECK_EQUAL(1,fx.split->numTFTFailures());
        CHECK(r.S.allFinite());
        CHECK(S.allFinite());
        CHECK((r.S - S).cwiseAbs().maxCoeff() == 0);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Output sizes and the layout colIdx = j*numPoints + k (point k, fibre j) of all multi-point outputs,
    /// against single-point evaluations. \a z is (fibres x points). With \a mixed, some columns must not be wrinkled,
    /// so that the NaN entries of the wrinkle outputs are exercised too. The non-wrinkled columns are counted by
    /// displayed state in \a numTaut and \a numSlack.
    template<short_t dim>
    static void SplitLayoutCheck(const SplitFixture<dim,false> & fx, const gsMatrix<real_t> & u, const gsMatrix<real_t> & z, const bool mixed,
                                 index_t & numTaut, index_t & numSlack)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const index_t np = u.cols(), nf = z.rows(), n = np*nf;
        CHECK_EQUAL(np,z.cols());

        const gsMatrix<real_t> S  = fx.split->eval3D_vector(0,u,z,MaterialOutput::VectorN);
        const gsMatrix<real_t> A  = fx.split->eval3D_matrix(0,u,z,MaterialOutput::MatrixA);
        const gsMatrix<real_t> th = fx.split->eval3D_theta(0,u,z,MaterialOutput::Generic);
        const gsMatrix<real_t> ga = fx.split->eval3D_gamma(0,u,z,MaterialOutput::Generic);
        const gsMatrix<real_t> tf = fx.split->eval3D_tensionfield(0,u,z,MaterialOutput::Generic);
        const gsMatrix<real_t> ps = fx.split->eval3D_pstress(0,u,z,MaterialOutput::Generic);
        const gsMatrix<real_t> st = fx.split->eval3D_strain(0,u,z);
        const gsMatrix<real_t> cv = fx.split->eval3D_CauchyVector(0,u,z,MaterialOutput::Generic);

        CHECK_EQUAL(3,S.rows());  CHECK_EQUAL(n,S.cols());
        CHECK_EQUAL(9,A.rows());  CHECK_EQUAL(n,A.cols());
        CHECK_EQUAL(1,th.rows()); CHECK_EQUAL(n,th.cols());
        CHECK_EQUAL(1,ga.rows()); CHECK_EQUAL(n,ga.cols());
        CHECK_EQUAL(1,tf.rows()); CHECK_EQUAL(n,tf.cols());
        CHECK_EQUAL(2,ps.rows()); CHECK_EQUAL(n,ps.cols());
        CHECK_EQUAL(3,st.rows()); CHECK_EQUAL(n,st.cols());
        CHECK_EQUAL(3,cv.rows()); CHECK_EQUAL(n,cv.cols());
        if (S.cols()!=n || A.cols()!=n || th.cols()!=n || ga.cols()!=n || tf.cols()!=n)
            return;

        // The points are distinct states (first fibre)
        for (index_t a = 0; a!=np; ++a)
            for (index_t b = a+1; b!=np; ++b)
                CHECK((S.col(a) - S.col(b)).norm() > 1e-6 * S.col(a).norm());

        index_t numWrinkled = 0;
        numTaut = numSlack = 0;
        for (index_t j = 0; j!=nf; ++j)
            for (index_t k = 0; k!=np; ++k)
            {
                const index_t col = j*np + k;
                const gsVector<real_t> uk = u.col(k);
                const gsMatrix<real_t> zk = SplitZ(z(j,k));
                const gsMatrix<real_t> S1  = fx.split->eval3D_vector(0,uk,zk,MaterialOutput::VectorN);
                const gsMatrix<real_t> A1  = fx.split->eval3D_matrix(0,uk,zk,MaterialOutput::MatrixA);
                const real_t th1 = fx.split->eval3D_theta(0,uk,zk,MaterialOutput::Generic)(0,0);
                const real_t ga1 = fx.split->eval3D_gamma(0,uk,zk,MaterialOutput::Generic)(0,0);
                const real_t tf1 = fx.split->eval3D_tensionfield(0,uk,zk,MaterialOutput::Generic)(0,0);

                const gsMatrix<real_t> Scol = S.col(col);
                const gsMatrix<real_t> Acol = A.col(col);
                CHECK(Scol.allFinite());
                CHECK(S1.allFinite());
                CHECK(Acol.allFinite());
                CHECK(A1.allFinite());
                CHECK((Scol-S1).cwiseAbs().maxCoeff() <= 1e-14 * S1.cwiseAbs().maxCoeff());
                CHECK((Acol-A1).cwiseAbs().maxCoeff() <= 1e-14 * A1.cwiseAbs().maxCoeff());
                CHECK(tf(0,col) == tf1);
                if (std::isnan(th1))
                {
                    CHECK(std::isnan(th(0,col)) && std::isnan(ga(0,col)) && std::isnan(ga1));
                }
                else
                {
                    CHECK(math::abs(th(0,col)-th1) <= 1e-14 * pi);
                    CHECK(math::abs(ga(0,col)-ga1) <= 1e-14 * math::abs(ga1));
                    CHECK(th(0,col) >= 0 && th(0,col) < pi);
                    ++numWrinkled;
                }
                // The wrinkle angle is defined exactly at the displayed-wrinkled points
                CHECK((tf(0,col) == 0) == !std::isnan(th(0,col)));
                if (tf(0,col) == 1)
                    ++numTaut;
                else if (tf(0,col) == -1)
                    ++numSlack;
            }
        gsInfo << "  wrinkled columns: " << numWrinkled << " of " << n << ", taut: " << numTaut << ", slack: " << numSlack << "\n";
        CHECK(numWrinkled >= 1);
        if (mixed)
            CHECK(numWrinkled < n);
    }

    /// Non-affine pair: the skewed affine deformation with one deformed control point moved, so that the
    /// state varies from wrinkled to taut over the points. The through-thickness coordinate differs per point.
    TEST(Layout_multipoint_2D)
    {
        const SplitStateSpec w = SplitWrinkledSpec();
        gsMultiPatch<real_t> ref = SplitAffinePatch(SplitSkewedJref());
        gsMultiPatch<real_t> def = SplitAffinePatch(SplitJdef(SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6));
        def.patch(0).coefs()(3,1) += 0.8;
        SplitFixture<2,false> fx(SplitMatNHIncomp,ref,def,0.01);

        gsMatrix<real_t> u(2,4), z(2,4);
        u << 0.2, 0.7, 0.4, 0.85,
             0.3, 0.2, 0.8, 0.65;
        for (index_t k = 0; k!=4; ++k)
        {
            z(0,k) = -0.5 + 0.1 * k;
            z(1,k) =  0.5 - 0.1 * k;
        }
        index_t numTaut, numSlack;
        SplitLayoutCheck<2>(fx,u,z,true,numTaut,numSlack);
        CHECK_EQUAL(0,numSlack);
        CHECK(numTaut >= 1);
    }

    /// Curved 3D pair: unit square stretched in x by 1.2, contracted in y by 0.75, one corner lifted. The
    /// metric varies through the thickness, so a fibre-index error shows.
    TEST(Layout_multipoint_3D)
    {
        gsMultiPatch<real_t> ref;
        ref.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        ref.addAutoBoundaries();
        ref.embed(3);
        gsMultiPatch<real_t> def = ref;
        def.patch(0).coefs().col(0) *= 1.2;
        def.patch(0).coefs().col(1) *= 0.75;
        def.patch(0).coefs()(3,2) += 0.2;
        SplitFixture<3,false> fx(SplitMatNHComp,ref,def,0.1);

        gsMatrix<real_t> u(2,4), z(2,4);
        u << 0.2, 0.7, 0.4, 0.85,
             0.3, 0.2, 0.8, 0.65;
        for (index_t k = 0; k!=4; ++k)
        {
            z(0,k) = -0.5 + 0.1 * k;
            z(1,k) =  0.5 - 0.1 * k;
        }
        index_t numTaut, numSlack;
        SplitLayoutCheck<3>(fx,u,z,false,numTaut,numSlack);

        // The metric varies through the thickness: the stress at a point depends on its own z, so that
        // using the z of another column would be caught by the column comparisons
        real_t maxZEffect = 0;
        for (index_t j = 0; j!=2; ++j)
            for (index_t k = 1; k!=4; ++k)
            {
                const gsVector<real_t> uk = u.col(k);
                const gsMatrix<real_t> Sk = fx.split->eval3D_vector(0,uk,SplitZ(z(j,k)),MaterialOutput::VectorN);
                const gsMatrix<real_t> S0 = fx.split->eval3D_vector(0,uk,SplitZ(z(j,0)),MaterialOutput::VectorN);
                maxZEffect = math::max(maxZEffect,(Sk-S0).norm() / S0.norm());
            }
        gsInfo << "  largest relative stress change with the fibre coordinate of another column: " << maxZEffect << "\n";
        CHECK(maxZEffect > 1e-10);

        const gsMatrix<real_t> S = fx.split->eval3D_vector(0,u,z,MaterialOutput::VectorN);
        real_t maxFibreDiff = 0;
        for (index_t k = 0; k!=4; ++k)
            maxFibreDiff = math::max(maxFibreDiff,(S.col(k) - S.col(4+k)).norm() / S.col(k).norm());
        gsInfo << "  largest relative fibre difference: " << maxFibreDiff << "\n";
        CHECK(maxFibreDiff > 1e-10);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Midsurface evaluates once per point (Constant integration) and has no bending output
    TEST(Midsurface_option)
    {
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6);
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);

        CHECK(fx.split->isMatIntegrated() == MatIntegration::Linear);
        CHECK(fx.split->isVecIntegrated() == MatIntegration::Linear);
        const gsMatrix<real_t> S0 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);

        gsOptionList o;
        o.addSwitch("Midsurface","",true);
        fx.split->setOptions(o);
        CHECK(fx.split->isMatIntegrated() == MatIntegration::Constant);
        CHECK(fx.split->isVecIntegrated() == MatIntegration::Constant);
        CHECK_THROW(fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorM), std::runtime_error);

        const gsMatrix<real_t> S1 = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        CHECK(S1.allFinite());
        CHECK(S0.allFinite());
        CHECK((S1-S0).cwiseAbs().maxCoeff() == 0);
        CHECK(S0.cwiseAbs().maxCoeff() > 0);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Single-point evaluations under an OpenMP parallel loop write the same values as the serial loop
    /// and as one multi-point call.
    TEST(OpenMP_deterministic)
    {
        const SplitStateSpec w = SplitWrinkledSpec();
        gsMultiPatch<real_t> ref = SplitAffinePatch(SplitSkewedJref());
        gsMultiPatch<real_t> def = SplitAffinePatch(SplitJdef(SplitSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6));
        def.patch(0).coefs()(3,1) += 0.8;
        SplitFixture<2,false> fx(SplitMatNHIncomp,ref,def,0.01);

        const index_t np = 64;
        gsMatrix<real_t> u(2,np), z(1,np);
        for (index_t k = 0; k!=np; ++k)
        {
            u(0,k) = 0.05 + 0.1 * (k % 8);
            u(1,k) = 0.05 + 0.1 * (k / 8);
        }
        z.setZero();
        gsInfo << "  omp_get_max_threads() = " << omp_get_max_threads() << "\n";
        if (omp_get_max_threads()==1)
            gsInfo << "  single thread: the parallel comparison is vacuous\n";

        gsMatrix<real_t> Sser(3,np), Aser(9,np), Spar(3,np), Apar(9,np);
        const gsMatrix<real_t> zmat = SplitZ(0);
        for (index_t k = 0; k!=np; ++k)
        {
            const gsVector<real_t> uk = u.col(k);
            Sser.col(k) = fx.split->eval3D_vector(0,uk,zmat,MaterialOutput::VectorN);
            Aser.col(k) = fx.split->eval3D_matrix(0,uk,zmat,MaterialOutput::MatrixA);
        }

#       pragma omp parallel for
        for (index_t k = 0; k<np; ++k)
        {
            const gsVector<real_t> uk = u.col(k);
            const gsMatrix<real_t> zm = SplitZ(0);
            Spar.col(k) = fx.split->eval3D_vector(0,uk,zm,MaterialOutput::VectorN);
            Apar.col(k) = fx.split->eval3D_matrix(0,uk,zm,MaterialOutput::MatrixA);
        }

        CHECK(Sser.allFinite());
        CHECK(Spar.allFinite());
        CHECK(Aser.allFinite());
        CHECK(Apar.allFinite());
        CHECK(Sser.cwiseAbs().maxCoeff() > 0);
        CHECK((Sser-Spar).cwiseAbs().maxCoeff() == 0);
        CHECK((Aser-Apar).cwiseAbs().maxCoeff() == 0);

        const gsMatrix<real_t> Smulti = fx.split->eval3D_vector(0,u,z,MaterialOutput::VectorN);
        const gsMatrix<real_t> Amulti = fx.split->eval3D_matrix(0,u,z,MaterialOutput::MatrixA);
        gsInfo << "  multi-point vs serial: S " << SplitRelErr(Smulti,Sser) << ", A " << SplitRelErr(Amulti,Aser) << "\n";
        CHECK(Smulti.allFinite());
        CHECK(Amulti.allFinite());
        CHECK((Smulti-Sser).cwiseAbs().maxCoeff() <= 1e-14 * Sser.cwiseAbs().maxCoeff());
        CHECK((Amulti-Aser).cwiseAbs().maxCoeff() <= 1e-14 * Aser.cwiseAbs().maxCoeff());
    }

    /// The guard-fallback counter is incremented atomically. The count is necessary but not sufficient
    /// evidence of atomicity: every call runs a root find and enters a named critical section, so the
    /// contention on the counter is low. The contended check is OpenMP_exact_counters_atomic.
    TEST(OpenMP_counter_atomic)
    {
        SplitFixture<2,true> fx(SplitMatSvK,SplitIdentity2(),1.0,1.0,0.0);
        gsMatrix<real_t> S, C, E;
        SplitAnisoInput(S,C,E);
        const gsMatrix<real_t> I2 = SplitIdentity2();
        detail_tftsplit::WrinkleOptions<real_t> opt;
        opt.tensionField = 2;
        gsInfo << "  omp_get_max_threads() = " << omp_get_max_threads() << "\n";

        fx.split->resetCounters();
        const index_t N = 2000;
#       pragma omp parallel for
        for (index_t i = 0; i<N; ++i)
            fx.split->wrinkle(S,C,E,I2,opt,false);
        CHECK_EQUAL(N,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Voigt direction n1 = (c^2,s^2,2cs) of the wrinkle covector (cos theta, sin theta)
    static gsVector<real_t,3> SplitN1(const real_t theta)
    {
        const real_t c = math::cos(theta), s = math::sin(theta);
        gsVector<real_t,3> n;
        n << c*c, s*s, 2*c*s;
        return n;
    }

    /// Voigt direction n2 = (-sc,sc,c^2-s^2) of the wrinkle covector (cos theta, sin theta)
    static gsVector<real_t,3> SplitN2(const real_t theta)
    {
        const real_t c = math::cos(theta), s = math::sin(theta);
        gsVector<real_t,3> n;
        n << -s*c, s*c, c*c-s*s;
        return n;
    }

    /// Symmetric 2x2 matrix [[a,b],[b,d]]
    static gsMatrix<real_t> SplitSym2(const real_t a, const real_t b, const real_t d)
    {
        gsMatrix<real_t> M(2,2);
        M << a, b, b, d;
        return M;
    }

    /// Ascending eigenvalues (2x1) of A x = lambda B x, B symmetric positive definite. This is the
    /// reference for all principal-value outputs, independent of the closed-form formulas of the class.
    static gsMatrix<real_t> SplitPencilEigenvalues(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
    {
        gsMatrix<real_t>::GenSelfAdjEigenSolver ges(A,B);
        const gsMatrix<real_t> ev = ges.eigenvalues();
        return ev;
    }

    /// Ascending principal values \a P (first two rows, one column) against the eigenvalues \a lam,
    /// relative to the largest eigenvalue magnitude
    static void SplitCheckPrincipal(const char * label, const gsMatrix<real_t> & P, const gsMatrix<real_t> & lam, const real_t tol)
    {
        CHECK(P.allFinite());
        CHECK(lam.allFinite());
        const real_t scale = lam.cwiseAbs().maxCoeff();
        const real_t err = math::max(math::abs(P(0,0)-lam(0,0)),math::abs(P(1,0)-lam(1,0))) / scale;
        gsInfo << "  " << label << ": rel.err = " << err << " (tol " << tol << ")\n";
        CHECK(err <= tol);
    }

    /// Largest antisymmetric part of \a C relative to its largest entry
    static real_t SplitAsymmetry(const gsMatrix<real_t> & C)
    {
        return (C - C.transpose()).cwiseAbs().maxCoeff() / C.cwiseAbs().maxCoeff();
    }

    /// Exact mode on a wrinkled point: the wrapped material, evaluated on its own at the returned elastic strain
    /// E' = E + gamma n1(theta), has n1.S_w(E') = n2.S_w(E') = 0. The linearised E'_lin = E + gamma_lin n1(theta_lin)
    /// of a twin fixture does not, which shows that the residual check can fail. The stress scale is ||S_w(E)||.
    static void SplitExactResidualCheck(const SplitMaterial mat, const real_t tolS)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec w = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        SplitFixture<2,false> fx(mat,SplitSkewedJref(),w.lam1,w.lam2,pi/6);
        SplitFixture<2,false> lin(mat,SplitSkewedJref(),w.lam1,w.lam2,pi/6);
        fx.split->options().setSwitch("Exact",true);
        fx.split->options().setInt("ExactLocalSolve",0);
        fx.split->resetCounters();

        const gsMatrix<real_t> Ep = fx.split->eval3D_strain(0,pt,zmat);
        const real_t theta = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
        const real_t gamma = fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);
        CHECK_EQUAL(0,fx.displayedState(pt,0));
        CHECK(gamma > 0);
        CHECK(theta >= 0 && theta < pi);

        const gsVector<real_t,3> n1 = SplitN1(theta), n2 = SplitN2(theta);
        const gsMatrix<real_t> gn1 = gamma * n1;
        const gsMatrix<real_t> dE = Ep - fx.EengDef();
        CHECK(dE.allFinite());
        CHECK((dE - gn1).cwiseAbs().maxCoeff() <= 1e-12 * gn1.cwiseAbs().maxCoeff());

        const gsMatrix<real_t> Sw  = fx.wrapped->eval3D_vector_C(fx.CmatDef(),0,pt,0,MaterialOutput::VectorN);
        const gsMatrix<real_t> Swp = fx.wrapped->eval3D_vector_C(fx.CmatFromE(Ep),0,pt,0,MaterialOutput::VectorN);
        CHECK(Sw.allFinite());
        CHECK(Swp.allFinite());
        const real_t scale = Sw.norm();
        const real_t r1 = math::abs(n1.dot(Swp.col(0))), r2 = math::abs(n2.dot(Swp.col(0)));
        gsInfo << "  gamma = " << gamma << ", theta = " << theta << ", |n1.S_w'|/|S_w| = " << r1/scale
               << ", |n2.S_w'|/|S_w| = " << r2/scale << " (tol 1e-8)\n";
        CHECK(r1 <= 1e-8 * scale);
        CHECK(r2 <= 1e-8 * scale);

        const gsMatrix<real_t> Seff = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        SPLIT_CHECK_REL("S_eff vs S_w(E')",Seff,Swp,tolS);

        const gsMatrix<real_t> p = fx.split->eval3D_pstress(0,pt,zmat,MaterialOutput::Generic);
        CHECK(p.allFinite());
        const index_t zero = math::abs(p(0,0)) < math::abs(p(1,0)) ? 0 : 1;
        CHECK(math::abs(p(zero,0)) <= 1e-8 * p.cwiseAbs().maxCoeff());
        CHECK(p(1-zero,0) > 0);

        const real_t thetaL = lin.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
        const real_t gammaL = lin.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);
        const gsVector<real_t,3> n1L = SplitN1(thetaL), n2L = SplitN2(thetaL);
        const gsMatrix<real_t> EpL = lin.EengDef() + gammaL * n1L;
        const gsMatrix<real_t> SwL = lin.wrapped->eval3D_vector_C(lin.CmatFromE(EpL),0,pt,0,MaterialOutput::VectorN);
        CHECK(SwL.allFinite());
        const real_t rL = math::max(math::abs(n1L.dot(SwL.col(0))),math::abs(n2L.dot(SwL.col(0))));
        gsInfo << "  linearised contrast: max|n.S_w(E'_lin)|/|S_w| = " << rL/scale << " (must exceed 1e-5)\n";
        CHECK(rL > 1e-5 * scale);

        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
        CHECK_EQUAL(0,fx.split->numExactFallbacks());
        CHECK_EQUAL(0,fx.split->numExactReclassified());
    }

    /// Wrapped-material residual at E' for the three nonlinear materials (see SplitExactResidualCheck)
    TEST(Exact_residual_wrapped)
    {
        SplitExactResidualCheck(SplitMatNHComp,1e-8);
        SplitExactResidualCheck(SplitMatNHIncomp,1e-10);
        SplitExactResidualCheck(SplitMatMRIncomp,1e-10);
    }

    /// Exact mode: the geometry path agrees with the prescribed-metric path at the same point, and the
    /// exact tangent is symmetric (the Hessian of Pipkin's relaxed energy) for eta in {0,1e-3}. With \a contrast,
    /// the linearised tangent at the same point must be visibly unsymmetric, so that the symmetry check can fail.
    static void SplitExactPathsCheck(const SplitMaterial mat, const gsMatrix<real_t> & Jref, const real_t tolPath, const bool contrast)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,false> fx(mat,Jref,w.lam1,w.lam2,pi/6);
        fx.split->options().setSwitch("Exact",true);
        fx.split->options().setInt("ExactLocalSolve",0);
        fx.split->resetCounters();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        const gsMatrix<real_t> Cmat = fx.CmatDef();
        CHECK_EQUAL(0,fx.stateAt(Cmat,pt,0));

        const real_t etas[2] = { 0.0, 1e-3 };
        for (index_t e = 0; e!=2; ++e)
        {
            gsInfo << "eta " << etas[e] << "\n";
            fx.split->options().setReal("Eta",etas[e]);
            const gsMatrix<real_t> Sgeo = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> SC   = fx.split->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::Generic);
            SPLIT_CHECK_REL("S geometry vs S_C",Sgeo,SC,tolPath);

            const gsMatrix<real_t> Ageo = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const gsMatrix<real_t> AC   = fx.split->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::Generic);
            SPLIT_CHECK_REL("A geometry vs A_C",Ageo,AC,tolPath);

            const gsMatrix<real_t> Cgeo = Ageo.reshape(3,3);
            const gsMatrix<real_t> Cpre = AC.reshape(3,3);
            gsInfo << "  asymmetry: geometry " << SplitAsymmetry(Cgeo) << ", prescribed metric " << SplitAsymmetry(Cpre) << " (tol 1e-6)\n";
            CHECK(Cgeo.allFinite());
            CHECK(Cpre.allFinite());
            CHECK(SplitAsymmetry(Cgeo) <= 1e-6);
            CHECK(SplitAsymmetry(Cpre) <= 1e-6);
        }

        if (contrast)
        {
            SplitFixture<2,false> lin(mat,Jref,w.lam1,w.lam2,pi/6);
            const gsMatrix<real_t> Alin = lin.split->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::Generic);
            const gsMatrix<real_t> Clin = Alin.reshape(3,3);
            gsInfo << "  linearised asymmetry: " << SplitAsymmetry(Clin) << " (must be at least 1e-2)\n";
            CHECK(Clin.allFinite());
            CHECK(SplitAsymmetry(Clin) >= 1e-2);
        }

        CHECK_EQUAL(0,fx.split->numExactFallbacks());
        CHECK_EQUAL(0,fx.split->numExactReclassified());
    }

    /// Path agreement and tangent symmetry in exact mode; the point with identity reference metric and rotated principal axes
    /// also shows the unsymmetric linearised tangent
    TEST(Exact_paths_symmetry)
    {
        SplitExactPathsCheck(SplitMatNHComp,SplitSkewedJref(),1e-7,false);
        SplitExactPathsCheck(SplitMatNHIncomp,SplitSkewedJref(),1e-10,false);
        SplitExactPathsCheck(SplitMatMRIncomp,SplitSkewedJref(),1e-10,false);
        SplitExactPathsCheck(SplitMatNHIncomp,SplitIdentity2(),1e-10,true);
    }

    /// Exact-mode analogue of SplitFDCheck: finite differences of the decorator's exact prescribed-metric stress
    /// against its exact tangent over wrinkled, taut and slack states and eta in {0,1e-3}. Every stencil column
    /// must stay in the intended state, and no exact solve may fail or reclassify a point. At the wrinkled eta = 0
    /// point the linearised tangent of a twin fixture must differ visibly, so that the exact tangent is tested
    /// where it matters.
    static void SplitExactFDCheck(const SplitMaterial mat, const real_t tol = 1e-4)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec specs[3] = { SplitWrinkledSpec(), SplitTautSpec(), SplitSlackSpec() };
        const real_t etas[2] = { 0.0, 1e-3 };
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const char * names[3] = { "wrinkled", "taut", "slack" };

        for (index_t s = 0; s!=3; ++s)
            for (index_t e = 0; e!=2; ++e)
            {
                gsInfo << "state " << names[s] << ", eta " << etas[e] << "\n";
                SplitFixture<2,false> fx(mat,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,pi/6);
                fx.split->options().setSwitch("Exact",true);
                fx.split->options().setInt("ExactLocalSolve",0);
                fx.split->options().setReal("Eta",etas[e]);
                fx.split->resetCounters();

                const gsMatrix<real_t> Cmat = fx.CmatDef();
                const gsMatrix<real_t> Eeng = fx.EengDef();
                CHECK_EQUAL(specs[s].state, fx.stateAt(Cmat,pt,z));

                const gsMatrix<real_t> Seff = fx.split->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::Generic);
                const gsMatrix<real_t> Ceff = fx.split->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::Generic).reshape(3,3);
                if (specs[s].state==-1 && etas[e]==0)
                {
                    const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::VectorN);
                    CHECK(Sw.allFinite());
                    CHECK(Seff.allFinite());
                    CHECK(Ceff.allFinite());
                    CHECK(Sw.cwiseAbs().maxCoeff() > 0);
                    CHECK(Seff.cwiseAbs().maxCoeff() == 0);
                    CHECK(Ceff.cwiseAbs().maxCoeff() == 0);
                }
                else
                {
                    index_t columns = 0;
                    bool allExpected = false;
                    const gsMatrix<real_t> J_FD = SplitFDJacobian(fx,pt,z,specs[s].state,Eeng,columns,allExpected);
                    CHECK_EQUAL(12,columns);
                    CHECK(allExpected);
                    SPLIT_CHECK_REL("exact J_C vs FD",Ceff,J_FD,tol);
                    if (specs[s].state==-1)
                        CHECK(Ceff.cwiseAbs().maxCoeff() > 0);
                }

                if (specs[s].state==0 && etas[e]==0)
                {
                    SplitFixture<2,false> lin(mat,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,pi/6);
                    const gsMatrix<real_t> Clin = lin.split->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::Generic).reshape(3,3);
                    gsInfo << "  exact vs linearised tangent: rel.diff = " << SplitRelErr(Clin,Ceff) << " (must exceed 1e-3)\n";
                    CHECK(Clin.allFinite());
                    CHECK(SplitRelErr(Clin,Ceff) > 1e-3);
                }
                CHECK_EQUAL(0,fx.split->numThetaFallbacks());
                CHECK_EQUAL(0,fx.split->numTFTFailures());
                CHECK_EQUAL(0,fx.split->numExactFallbacks());
                CHECK_EQUAL(0,fx.split->numExactReclassified());
            }
    }

    /// Compressible Neo-Hookean, exact tangent against finite differences
    TEST(Exact_FD_C_NH_Comp)
    {
        SplitExactFDCheck(SplitMatNHComp);
    }

    /// Incompressible Neo-Hookean, exact tangent against finite differences
    TEST(Exact_FD_C_NH_Incomp)
    {
        SplitExactFDCheck(SplitMatNHIncomp);
    }

    /// Incompressible Mooney-Rivlin, exact tangent against finite differences
    TEST(Exact_FD_C_MR_Incomp)
    {
        SplitExactFDCheck(SplitMatMRIncomp);
    }

    /// For a wrapped material that is linear in E the exact and the linearised wrinkled state coincide. No
    /// nonlinear St. Venant-Kirchhoff material is constructible in 2D, so the wrapped material is the linear SvK
    /// (S = C E, constant C) inside the nonlinear decorator, where a test-side provider drives the exact branch.
    TEST(Exact_SvK_equals_linearised)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,false> ex(SplitMatSvK,SplitSkewedJref(),w.lam1,w.lam2,pi/6);
        SplitFixture<2,false> li(SplitMatSvK,SplitSkewedJref(),w.lam1,w.lam2,pi/6);
        ex.split->options().setSwitch("Exact",true);
        ex.split->options().setInt("ExactLocalSolve",0);
        ex.split->resetCounters();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        const gsMatrix<real_t> Cmat = ex.CmatDef();
        CHECK_EQUAL(0,ex.displayedState(pt,0));

        const gsMatrix<real_t> Sge = ex.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> Sgl = li.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> Age = ex.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        const gsMatrix<real_t> Agl = li.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        SPLIT_CHECK_REL("S geometry, exact vs linearised",Sge,Sgl,1e-8);
        SPLIT_CHECK_REL("C geometry, exact vs linearised",Age,Agl,1e-8);

        const gsMatrix<real_t> Sce = ex.split->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::Generic);
        const gsMatrix<real_t> Scl = li.split->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::Generic);
        const gsMatrix<real_t> Ace = ex.split->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::Generic);
        const gsMatrix<real_t> Acl = li.split->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::Generic);
        SPLIT_CHECK_REL("S prescribed metric, exact vs linearised",Sce,Scl,1e-8);
        SPLIT_CHECK_REL("C prescribed metric, exact vs linearised",Ace,Acl,1e-8);

        // The exact branch itself: a provider that evaluates the wrapped material at the strain E'
        const gsMatrix<real_t> S3 = ex.wrapped->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::VectorN);
        const gsMatrix<real_t> C3 = ex.wrapped->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::MatrixA);
        const gsMatrix<real_t> E3 = ex.EengDef();
        SplitFixture<2,false>::SplitType::SC_provider provider =
            [&](const gsMatrix<real_t> & Eel, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            const gsMatrix<real_t> Cm = ex.CmatFromE(Eel);
            Se = ex.wrapped->eval3D_vector_C(Cm,0,pt,0,MaterialOutput::VectorN);
            Ce = ex.wrapped->eval3D_matrix_C(Cm,0,pt,0,MaterialOutput::MatrixA);
        };
        detail_tftsplit::WrinkleOptions<real_t> opt = ex.split->wrinkleOptions();
        CHECK(opt.exact);
        const detail_tftsplit::WrinkleResult<real_t> re =
            ex.split->wrinkle(S3,C3,E3,ex.Gref,opt,true,SplitFixture<2,false>::SplitType::dC_provider(),provider);
        opt.exact = false;
        const detail_tftsplit::WrinkleResult<real_t> rl = ex.split->wrinkle(S3,C3,E3,ex.Gref,opt,true);
        gsInfo << "  provider calls of the exact solve: " << re.exactIterations << "\n";
        CHECK_EQUAL(0,re.state);
        CHECK(re.exactIterations >= 1);
        CHECK(!re.exactFallback);
        CHECK(!re.exactReclassified);
        SPLIT_CHECK_REL("S wrinkle(), exact vs linearised",re.S,rl.S,1e-8);
        SPLIT_CHECK_REL("C wrinkle(), exact vs linearised",re.C,rl.C,1e-8);

        CHECK_EQUAL(0,ex.split->numExactFallbacks());
        CHECK_EQUAL(0,ex.split->numExactReclassified());
    }

    /// Incompressible Neo-Hookean over the family lambda_2(delta) = (1-delta) lambda_1^(-1/2): delta > 0 is wrinkled,
    /// delta = 0 is the taut/wrinkled boundary. The exact and the linearised stress differ at second order in the
    /// wrinkle amplitude, S_ex - S_lin = O(gamma^2); the slope of log d versus log gamma_lin is fitted by least squares.
    TEST(Exact_differs_Ogamma2)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t lam1 = SplitWrinkledSpec().lam1;
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);

        const real_t deltas[4] = { 0.01, 0.02, 0.04, 0.2 };
        real_t gl[4], gx[4], d[4];
        gsInfo << "  delta, gamma_lin, gamma_ex, d\n";
        for (index_t i = 0; i!=4; ++i)
        {
            SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),lam1,(1-deltas[i])/math::sqrt(lam1),pi/6);
            fx.split->resetCounters();
            const gsMatrix<real_t> Sl = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Al = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            gl[i] = fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);

            fx.split->options().setSwitch("Exact",true);
            fx.split->options().setInt("ExactLocalSolve",0);
            CHECK_EQUAL(0,fx.displayedState(pt,0));
            const gsMatrix<real_t> Sx = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Ax = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            gx[i] = fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);
            CHECK(Sl.allFinite());
            CHECK(Sx.allFinite());
            CHECK(Al.allFinite());
            CHECK(Ax.allFinite());
            CHECK(gl[i] > 0);
            CHECK(gx[i] > 0);

            d[i] = (Sx-Sl).norm() / Sl.norm();
            gsInfo << "  " << deltas[i] << ", " << gl[i] << ", " << gx[i] << ", " << d[i] << "\n";
            CHECK_EQUAL(0,fx.split->numExactFallbacks());
            CHECK_EQUAL(0,fx.split->numExactReclassified());

            if (i==3)
            {
                gsInfo << "  20 % over-contraction: d = " << d[i] << ", tangent rel.diff = " << SplitRelErr(Ax,Al) << "\n";
                CHECK(d[i] > 1e-4);
                CHECK(SplitRelErr(Ax,Al) > 1e-3);
            }
        }

        // Least-squares slope of log d against log gamma_lin over the small-amplitude values
        const index_t n = 3;
        real_t sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (index_t i = 0; i!=n; ++i)
        {
            const real_t x = math::log(gl[i]), y = math::log(d[i]);
            sx += x; sy += y; sxx += x*x; sxy += x*y;
        }
        const real_t slope = (n*sxy - sx*sy) / (n*sxx - sx*sx);
        gsInfo << "  gamma_lin spans a factor " << gl[n-1]/gl[0] << ", slope of log d vs log gamma_lin = " << slope << "\n";
        CHECK(gl[n-1] >= 2 * gl[0]);
        CHECK(slope >= 1.5 && slope <= 2.5);
    }

    /// Incompressible Neo-Hookean as the wrinkle amplitude gamma tends to 0 on the wrinkled side. The difference
    /// between the exact and the linearised stress and tangent does not grow as delta decreases, and both tangents tend
    /// to the reduced tangent C_w - C_w n1 n1^T C_w / (n1^T C_w n1) of the wrapped material. The tangent jump to the taut
    /// wrapped tangent at delta = 0 is intrinsic and is not compared.
    TEST(Exact_continuity_gamma0)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t lam1 = SplitWrinkledSpec().lam1;
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);

        const real_t deltas[4] = { 1e-2, 1e-4, 1e-6, 1e-8 };
        real_t dSprev = 0, dCprev = 0;
        gsInfo << "  delta, gamma_ex, d_S, d_C\n";
        for (index_t i = 0; i!=4; ++i)
        {
            SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),lam1,(1-deltas[i])/math::sqrt(lam1),pi/6);
            fx.split->resetCounters();
            const gsMatrix<real_t> Sl = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Al = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);

            fx.split->options().setSwitch("Exact",true);
            fx.split->options().setInt("ExactLocalSolve",0);
            CHECK_EQUAL(0,fx.displayedState(pt,0));
            const gsMatrix<real_t> Sx = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Ax = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
            const real_t gamma = fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);
            const real_t theta = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
            CHECK(Sl.allFinite());
            CHECK(Sx.allFinite());
            CHECK(Al.allFinite());
            CHECK(Ax.allFinite());
            CHECK(gamma > 0);

            const real_t dS = (Sx-Sl).norm() / Sl.norm();
            const real_t dC = SplitRelErr(Ax,Al);
            gsInfo << "  " << deltas[i] << ", " << gamma << ", " << dS << ", " << dC << "\n";
            if (i > 0)
            {
                CHECK(dS <= math::max(dSprev,real_t(1e-14)));
                CHECK(dC <= math::max(dCprev,real_t(1e-14)));
            }
            dSprev = dS;
            dCprev = dC;

            if (i==3)
            {
                const gsMatrix<real_t> Cw = fx.wrapped->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA).reshape(3,3);
                const gsVector<real_t,3> n1 = SplitN1(theta);
                const gsVector<real_t,3> Cn1 = Cw * n1;
                const gsVector<real_t,3> CTn1 = Cw.transpose() * n1;
                const gsMatrix<real_t> Clim = Cw - Cn1 * CTn1.transpose() / n1.dot(Cn1);
                const gsMatrix<real_t> Cx = Ax.reshape(3,3);
                const gsMatrix<real_t> Cl = Al.reshape(3,3);
                SPLIT_CHECK_REL("C_ex vs reduced tangent",Cx,Clim,1e-5);
                SPLIT_CHECK_REL("C_lin vs reduced tangent",Cl,Clim,1e-5);
            }
            CHECK_EQUAL(0,fx.split->numExactFallbacks());
            CHECK_EQUAL(0,fx.split->numExactReclassified());
        }

        // Taut side: nothing is modified by the exact mode
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),lam1,(1+1e-3)/math::sqrt(lam1),pi/6);
        const gsMatrix<real_t> Sl = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> Al = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        fx.split->options().setSwitch("Exact",true);
        CHECK_EQUAL(1,fx.displayedState(pt,0));
        const gsMatrix<real_t> Sx = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> Ax = fx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA);
        SPLIT_CHECK_SAME(Sx,Sl);
        SPLIT_CHECK_SAME(Ax,Al);
    }

    /// A wrapped material that is linear in E (linear = true) ignores the Exact option: no error, the option
    /// does not reach the core, all outputs are bitwise those of the linearised mode and no counter moves.
    TEST(Exact_linear_true_ignored)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec w = SplitWrinkledSpec();
        SplitFixture<2,true> fx(SplitMatSvK,SplitSkewedJref(),w.lam1,w.lam2,pi/6);
        fx.split->resetCounters();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        const gsMatrix<real_t> Cmat = fx.CmatDef();
        CHECK_EQUAL(0,fx.displayedState(pt,0));

        const SplitFixture<2,true> & cfx = fx;
        auto evalAll = [&]()
        {
            std::vector<gsMatrix<real_t> > v;
            v.push_back(cfx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN));
            v.push_back(cfx.split->eval3D_matrix(0,pt,zmat,MaterialOutput::MatrixA));
            v.push_back(cfx.split->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::Generic));
            v.push_back(cfx.split->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::Generic));
            v.push_back(cfx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic));
            v.push_back(cfx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic));
            v.push_back(cfx.split->eval3D_strain(0,pt,zmat));
            return v;
        };
        const std::vector<gsMatrix<real_t> > off = evalAll();

        gsOptionList o;
        o.addSwitch("Exact","",true);
        bool thrown = false;
        try { fx.split->setOptions(o); }
        catch (...) { thrown = true; }
        CHECK(!thrown);
        CHECK(fx.split->options().getSwitch("Exact"));
        CHECK(!fx.split->wrinkleOptions().exact);

        const std::vector<gsMatrix<real_t> > on = evalAll();
        CHECK_EQUAL(off.size(),on.size());
        for (size_t i = 0; i!=off.size() && i!=on.size(); ++i)
            SPLIT_CHECK_SAME(on[i],off[i]);

        // The pointwise core needs no provider for linear = true even when asked for the exact mode
        const gsMatrix<real_t> S3 = fx.wrapped->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::VectorN);
        const gsMatrix<real_t> C3 = fx.wrapped->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::MatrixA);
        const gsMatrix<real_t> E3 = fx.EengDef();
        detail_tftsplit::WrinkleOptions<real_t> opt = fx.split->wrinkleOptions();
        opt.exact = false;
        const detail_tftsplit::WrinkleResult<real_t> rl = fx.split->wrinkle(S3,C3,E3,fx.Gref,opt,true);
        opt.exact = true;
        detail_tftsplit::WrinkleResult<real_t> re;
        bool thrownCore = false;
        try { re = fx.split->wrinkle(S3,C3,E3,fx.Gref,opt,true); }
        catch (...) { thrownCore = true; }
        CHECK(!thrownCore);
        if (!thrownCore)
        {
            CHECK_EQUAL(0,re.state);
            CHECK(!re.exactFallback);
            CHECK(!re.exactReclassified);
            CHECK_EQUAL(0,re.exactIterations);
            SPLIT_CHECK_SAME(re.S,rl.S);
            SPLIT_CHECK_SAME(re.C,rl.C);
        }

        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
        CHECK_EQUAL(0,fx.split->numExactFallbacks());
        CHECK_EQUAL(0,fx.split->numExactReclassified());
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Input of the pointwise core (Gcov = I) with a linearised reference call, and synthetic providers of the
    /// wrapped stress and tangent whose exact-mode solution is known in closed form.
    ///
    /// With theta_lin, gamma_lin of the linearised call, c = cos theta_lin, s = sin theta_lin, the stress-like vectors
    ///   T_m = (s^2, c^2, -sc):    n1.T_m = n2.T_m = 0,  n4.T_m = 1,
    ///   w   = (-2sc, 2sc, c^2-s^2): n1.w = 0,           n2.w = 1
    /// (n4 = n3 + n1 = (s^2,c^2,-2cs) is the stress along the wrinkles) shift the residual components or the
    /// wrinkle tension of an affine provider S(E') = S + C (E' - E) independently of each other.
    struct SplitSynthetic
    {
        typedef SplitFixture<2,false>::SplitType         SplitType;
        typedef SplitType::SC_provider                   Provider;
        typedef detail_tftsplit::WrinkleResult<real_t>   Result;
        typedef detail_tftsplit::WrinkleOptions<real_t>  Options;

        /// affine S + C (E'-E); taut: -S + C (E'-E); slack(alpha): affine - alpha T_m; NaN; negJ11: tangent -C;
        /// saddle(beta): tangent C - beta w w^T; switch(kappa): affine + kappa ||S|| w
        enum Kind { Affine, Taut, Slack, NaN, NegJ11, Saddle, Switch };

        SplitFixture<2,false> fx;
        gsMatrix<real_t> I2, S, C, E, Tm, w;
        bool aniso;
        Options optLin, optEx;
        Result lin, linNoC;
        real_t thetaLin, gammaLin, swLin, sNorm;

        /// Isotropic input (the closed-form direction passes the guard) or anisotropic one (it does not)
        explicit SplitSynthetic(const bool anisotropic)
        :
        fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0), I2(SplitIdentity2()), aniso(anisotropic)
        {
            if (aniso)
                SplitAnisoInput(S,C,E);
            else
            {
                C = SplitIsoC(0.3);
                E.resize(3,1);
                E << 0.02, -0.01, 0.006;
                S = C * E;
            }
            sNorm = S.norm();

            optLin.tensionField = 2;
            optLin.eta = 0;
            optLin.slackMultiplier = 0;
            optLin.exact = false;
            optEx = optLin;
            optEx.exact = true;

            lin    = fx.split->wrinkle(S,C,E,I2,optLin,true);
            linNoC = fx.split->wrinkle(S,C,E,I2,optLin,false);
            CHECK_EQUAL(0,lin.state);
            CHECK_EQUAL(aniso,lin.guardFallback);
            thetaLin = lin.theta;
            gammaLin = lin.gamma;

            const real_t c = math::cos(thetaLin), s = math::sin(thetaLin);
            Tm.resize(3,1);
            Tm << s*s, c*c, -s*c;
            w.resize(3,1);
            w << -2*s*c, 2*s*c, c*c-s*s;
            gsVector<real_t,3> n4;
            n4 << s*s, c*c, -2*c*s;
            const gsVector<real_t,3> n1 = SplitN1(thetaLin), n2 = SplitN2(thetaLin);
            CHECK(math::abs(n1.dot(Tm.col(0))) <= 1e-15);
            CHECK(math::abs(n2.dot(Tm.col(0))) <= 1e-15);
            CHECK(math::abs(n4.dot(Tm.col(0)) - 1) <= 1e-15);
            CHECK(math::abs(n1.dot(w.col(0))) <= 1e-15);
            CHECK(math::abs(n2.dot(w.col(0)) - 1) <= 1e-15);

            swLin = n4.dot(lin.S.col(0));
            CHECK(gammaLin > 0);
            CHECK(swLin > 0);
        }

        /// Provider of the given kind; \a par is alpha, beta or kappa
        Provider provider(const Kind kind, const real_t par = 0) const
        {
            const gsMatrix<real_t> S0 = S, C0 = C, E0 = E, Tm0 = Tm, w0 = w;
            const real_t sn = sNorm;
            const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
            return [=](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
            {
                switch (kind)
                {
                case Taut:
                    Se = -S0 + C0 * (Ep - E0);
                    Ce = C0;
                    break;
                case Slack:
                    Se = S0 + C0 * (Ep - E0) - par * Tm0;
                    Ce = C0;
                    break;
                case NaN:
                    Se.setConstant(3,1,nan);
                    Ce.setConstant(3,3,nan);
                    break;
                case NegJ11:
                    Se = S0 + C0 * (Ep - E0);
                    Ce = -C0;
                    break;
                case Saddle:
                    Se = S0 + C0 * (Ep - E0);
                    Ce = C0 - par * (w0 * w0.transpose());
                    break;
                case Switch:
                    Se = S0 + C0 * (Ep - E0) + (par * sn) * w0;
                    Ce = C0;
                    break;
                default:
                    Se = S0 + C0 * (Ep - E0);
                    Ce = C0;
                    break;
                }
            };
        }

        /// The pointwise core on the stored input
        Result run(const Options & opt, const bool tangent, const Provider & prov = Provider()) const
        {
            return fx.split->wrinkle(S,C,E,I2,opt,tangent,SplitType::dC_provider(),prov);
        }

        /// n2^T C n2 at theta_lin
        real_t m22() const
        {
            const gsVector<real_t,3> n2 = SplitN2(thetaLin);
            const gsVector<real_t,3> Cn2 = C * n2;
            return n2.dot(Cn2);
        }
    };

    /// Without a provider, the exact branch of a wrinkled point is an error; taut and slack points return before it
    TEST(Exact_requires_provider)
    {
        SplitSynthetic sy(false);
        CHECK_THROW(sy.run(sy.optEx,false), std::runtime_error);
        CHECK_THROW(sy.run(sy.optEx,true), std::runtime_error);

        gsMatrix<real_t> E(3,1), S;
        for (index_t t = 0; t!=2; ++t)
        {
            E << 0.02, 0.01, 0.0;
            S = sy.C * E;
            CHECK_EQUAL(1,sy.fx.split->wrinkle(S,sy.C,E,sy.I2,sy.optEx,t==0).state);
            E << -0.02, -0.01, 0.0;
            S = sy.C * E;
            CHECK_EQUAL(-1,sy.fx.split->wrinkle(S,sy.C,E,sy.I2,sy.optEx,t==0).state);
        }
    }

    /// An affine provider reproduces the linearised result: the exact solve of an affine problem is exact and its
    /// tangent reduces to the linearised one. The isotropic input runs the scalar mode, the anisotropic one the 2x2 mode.
    TEST(Exact_synthetic_affine)
    {
        for (index_t a = 0; a!=2; ++a)
        {
            gsInfo << (a==0 ? "isotropic" : "anisotropic") << " input\n";
            SplitSynthetic sy(a==1);
            sy.fx.split->resetCounters();
            const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Affine);

            const SplitSynthetic::Result r = sy.run(sy.optEx,true,prov);
            CHECK_EQUAL(0,r.state);
            CHECK(!r.exactFallback);
            CHECK(!r.exactReclassified);
            CHECK(r.exactIterations >= 1);
            CHECK_EQUAL(a==1,r.guardFallback);
            SPLIT_CHECK_REL("S exact vs linearised",r.S,sy.lin.S,1e-12);
            SPLIT_CHECK_REL("C exact vs linearised",r.C,sy.lin.C,1e-12);
            gsInfo << "  asymmetry of the exact tangent: " << SplitAsymmetry(r.C) << " (tol 1e-12)\n";
            CHECK(SplitAsymmetry(r.C) <= 1e-12);

            const SplitSynthetic::Result r0 = sy.run(sy.optEx,false,prov);
            CHECK_EQUAL(0,r0.C.size());
            CHECK_EQUAL(0,r0.state);
            SPLIT_CHECK_REL("S exact vs linearised, no tangent",r0.S,sy.linNoC.S,1e-12);
            SPLIT_CHECK_SAME(r0.E,r.E);
            CHECK(r0.theta == r.theta);
            CHECK(r0.gamma == r.gamma);
            CHECK(r0.exactIterations >= 1);
            CHECK(!r0.exactFallback);
            CHECK(!r0.exactReclassified);
            CHECK_EQUAL(a==1,r0.guardFallback);

            CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());
            CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
        }
    }

    /// The provider S(E') = -S + C (E'-E) is solved by gamma = -gamma_lin < 0: the point is reclassified as taut
    /// and returns the unmodified S, C and E, with the wrinkle quantities undefined
    TEST(Exact_synthetic_reclassify_taut)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Taut);
        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitSynthetic::Result r = sy.run(sy.optEx,t==0,prov);
            CHECK_EQUAL(1,r.state);
            CHECK(r.exactReclassified);
            CHECK(!r.exactFallback);
            CHECK_EQUAL(1,sy.fx.split->numExactReclassified());
            CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());
            SPLIT_CHECK_SAME(r.S,sy.S);
            SPLIT_CHECK_SAME(r.E,sy.E);
            if (t==0)
                SPLIT_CHECK_SAME(r.C,sy.C);
            else
                CHECK_EQUAL(0,r.C.size());
            CHECK(std::isnan(r.theta) && std::isnan(r.gamma) && std::isnan(r.gammaHat) && std::isnan(r.residual));
        }
    }

    /// The provider with the wrinkle tension n4.S' shifted to -sw_lin converges to a stationary point with
    /// negative wrinkle tension and det Jac_0 < 0 whose unmodified strain is outside the slack band, i.e. a
    /// saddle: the linearised result is returned and counted
    TEST(Exact_synthetic_reclassify_slack)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Slack,2 * sy.swLin);
        for (index_t i = 0; i!=2; ++i)
            for (index_t t = 0; t!=2; ++t)
            {
                SplitSynthetic::Options opt = sy.optEx;
                opt.eta = 1e-3;
                opt.slackMultiplier = (i==0) ? 0.0 : 0.25;
                sy.fx.split->resetCounters();
                const SplitSynthetic::Result r = sy.run(opt,t==0,prov);
                SplitSynthetic::Options optL = opt;
                optL.exact = false;
                const SplitSynthetic::Result ref = sy.run(optL,t==0);
                CHECK_EQUAL(0,r.state);
                CHECK(r.exactFallback);
                CHECK(r.exactSaddle);
                CHECK(!r.exactWrongBranch);
                CHECK(!r.exactReclassified);
                CHECK_EQUAL(1,sy.fx.split->numExactSaddle());
                CHECK_EQUAL(1,sy.fx.split->numExactFallbacks());
                CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
                CHECK(std::isfinite(r.theta) && std::isfinite(r.gamma));
                SPLIT_CHECK_SAME(r.S,ref.S);
                SPLIT_CHECK_SAME(r.E,ref.E);
                if (t==0)
                    SPLIT_CHECK_SAME(r.C,ref.C);
                else
                    CHECK_EQUAL(0,r.C.size());
            }

        // Control: without the shift the wrinkle tension is positive and the point stays wrinkled
        sy.fx.split->resetCounters();
        const SplitSynthetic::Result rc = sy.run(sy.optEx,true,sy.provider(SplitSynthetic::Slack,0.0));
        CHECK_EQUAL(0,rc.state);
        CHECK(!rc.exactReclassified);
        CHECK(!rc.exactFallback);
        CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
    }

    /// The NaN provider ends the solve after its first call (non-convergence); the NegJ11 provider converges
    /// on its first call to a stationary point with n1^T C n1 < 0, which is a saddle; both return the
    /// linearised result bit for bit, with the fallback counted
    TEST(Exact_synthetic_nonconvergence)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Kind kinds[2] = { SplitSynthetic::NaN, SplitSynthetic::NegJ11 };
        for (index_t k = 0; k!=2; ++k)
        {
            const SplitSynthetic::Provider prov = sy.provider(kinds[k]);
            for (index_t t = 0; t!=2; ++t)
            {
                const SplitSynthetic::Result & ref = (t==0) ? sy.lin : sy.linNoC;
                sy.fx.split->resetCounters();
                const SplitSynthetic::Result r = sy.run(sy.optEx,t==0,prov);
                CHECK(r.exactFallback);
                CHECK(!r.exactReclassified);
                CHECK_EQUAL(1,r.exactIterations);
                CHECK_EQUAL(1,sy.fx.split->numExactFallbacks());
                CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
                CHECK_EQUAL(0,r.state);
                if (k==1)
                {
                    CHECK(r.exactSaddle);
                    CHECK(!r.exactWrongBranch);
                    CHECK_EQUAL(1,sy.fx.split->numExactSaddle());
                    CHECK_EQUAL(0,sy.fx.split->numExactNonconverged());
                }
                else
                {
                    CHECK(!r.exactSaddle);
                    CHECK(!r.exactWrongBranch);
                    CHECK_EQUAL(1,sy.fx.split->numExactNonconverged());
                    CHECK_EQUAL(0,sy.fx.split->numExactSaddle());
                }
                SPLIT_CHECK_SAME(r.S,ref.S);
                SPLIT_CHECK_SAME(r.E,ref.E);
                if (t==0)
                    SPLIT_CHECK_SAME(r.C,ref.C);
                else
                    CHECK_EQUAL(0,r.C.size());
                CHECK(std::isfinite(r.theta) && std::isfinite(r.gamma) && std::isfinite(r.gammaHat) && std::isfinite(r.residual));
                CHECK(r.theta == ref.theta);
                CHECK(r.gamma == ref.gamma);
                CHECK(r.gammaHat == ref.gammaHat);
                CHECK(r.residual == ref.residual);
            }
        }
    }

    /// The tangent C - beta w w^T with beta = M22 + sw_lin/gamma_lin leaves the solution (gamma_lin, theta_lin)
    /// unchanged but makes det Jac_0 = -M11 sw_lin - 2 gamma M12^2 < 0: a saddle, not a minimum of the relaxed
    /// energy. The result is the linearised one with the fallback counted. With beta = 0 the same solution is admissible.
    TEST(Exact_synthetic_saddle)
    {
        SplitSynthetic sy(false);
        const real_t beta = sy.m22() + sy.swLin / sy.gammaLin;
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Saddle,beta);
        for (index_t t = 0; t!=2; ++t)
        {
            const SplitSynthetic::Result & ref = (t==0) ? sy.lin : sy.linNoC;
            sy.fx.split->resetCounters();
            const SplitSynthetic::Result r = sy.run(sy.optEx,t==0,prov);
            CHECK(r.exactFallback);
            CHECK(!r.exactReclassified);
            CHECK(r.exactIterations >= 1);
            CHECK_EQUAL(1,sy.fx.split->numExactFallbacks());
            CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
            CHECK_EQUAL(0,r.state);
            SPLIT_CHECK_SAME(r.S,ref.S);
            SPLIT_CHECK_SAME(r.E,ref.E);
            if (t==0)
                SPLIT_CHECK_SAME(r.C,ref.C);
            else
                CHECK_EQUAL(0,r.C.size());
            CHECK(r.theta == ref.theta);
            CHECK(r.gamma == ref.gamma);
            CHECK(r.gammaHat == ref.gammaHat);
            CHECK(r.residual == ref.residual);
        }

        sy.fx.split->resetCounters();
        const SplitSynthetic::Result r0 = sy.run(sy.optEx,true,sy.provider(SplitSynthetic::Saddle,0.0));
        CHECK(!r0.exactFallback);
        CHECK(!r0.exactReclassified);
        CHECK_EQUAL(0,r0.state);
        CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());
    }

    /// A constant shift kappa ||S|| w of the provider stress moves only the residual component n2.S'. The scalar
    /// mode (closed-form theta accepted by the guard) converges in gamma with |R2| = kappa ||S|| above the tolerance,
    /// switches to the 2x2 mode and moves theta; with kappa = 0 theta stays at its linearised value.
    TEST(Exact_synthetic_switch_2x2)
    {
        SplitSynthetic sy(false);
        const real_t kappas[2] = { 1e-4, 0.0 };
        for (index_t i = 0; i!=2; ++i)
        {
            const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Switch,kappas[i]);
            sy.fx.split->resetCounters();
            const SplitSynthetic::Result r = sy.run(sy.optEx,true,prov);
            const SplitSynthetic::Result r0 = sy.run(sy.optEx,false,prov);
            CHECK(!r.guardFallback);
            CHECK(!r.exactFallback);
            CHECK(!r.exactReclassified);
            CHECK_EQUAL(0,r.state);
            const real_t dTheta = SplitThetaDiff(r.theta,sy.thetaLin);
            gsInfo << "  kappa " << kappas[i] << ": theta shift = " << dTheta << ", residual = " << r.residual / sy.sNorm << " ||S||\n";
            if (i==0)
                CHECK(dTheta > 1e-8);
            else
                CHECK(dTheta <= 1e-14);

            // Without the tangent the same solve runs
            CHECK(!r0.guardFallback);
            CHECK(!r0.exactFallback);
            CHECK(!r0.exactReclassified);
            CHECK_EQUAL(0,r0.state);
            CHECK_EQUAL(0,r0.C.size());
            if (i==0)
                CHECK(SplitThetaDiff(r0.theta,sy.thetaLin) > 1e-8);
            else
                CHECK(SplitThetaDiff(r0.theta,sy.thetaLin) <= 1e-14);
            SPLIT_CHECK_SAME(r0.S,r.S);
            SPLIT_CHECK_SAME(r0.E,r.E);
            CHECK(r0.theta == r.theta);
            CHECK(r0.gamma == r.gamma);

            const gsVector<real_t,3> n1 = SplitN1(r.theta), n2 = SplitN2(r.theta);
            const gsMatrix<real_t> Ep = sy.E + r.gamma * n1;
            gsMatrix<real_t> Sx, Cx;
            prov(Ep,Sx,Cx);
            CHECK(Sx.allFinite());
            CHECK(math::abs(n1.dot(Sx.col(0))) <= 1e-8 * sy.sNorm);
            CHECK(math::abs(n2.dot(Sx.col(0))) <= 1e-8 * sy.sNorm);
            CHECK(r.residual <= 1e-8 * sy.sNorm);
            SPLIT_CHECK_REL("S_eff vs provider at E'",r.S,Sx,1e-12);
            gsInfo << "  asymmetry of the exact tangent: " << SplitAsymmetry(r.C) << " (tol 1e-10)\n";
            CHECK(r.C.allFinite());
            CHECK(SplitAsymmetry(r.C) <= 1e-10);
            CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());
        }
    }

    /// The exact counters are incremented atomically. The taut provider takes the reclassification path, which
    /// increments without entering any critical section, so that the calls contend on the counter itself; the
    /// NaN provider ends in the fallback counter.
    TEST(OpenMP_exact_counters_atomic)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider taut = sy.provider(SplitSynthetic::Taut);
        const SplitSynthetic::Provider nan  = sy.provider(SplitSynthetic::NaN);
        gsInfo << "  omp_get_max_threads() = " << omp_get_max_threads() << "\n";
        const index_t N = 20000;

        sy.fx.split->resetCounters();
#       pragma omp parallel for
        for (index_t i = 0; i<N; ++i)
            sy.run(sy.optEx,false,taut);
        CHECK_EQUAL(N,sy.fx.split->numExactReclassified());
        CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());

        sy.fx.split->resetCounters();
#       pragma omp parallel for
        for (index_t i = 0; i<N; ++i)
            sy.run(sy.optEx,false,nan);
        CHECK_EQUAL(N,sy.fx.split->numExactFallbacks());
        CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// The outputs of all three states against independent oracles: eigenvalues of the metric pencil for the principal
    /// strains, stresses and Cauchy stresses, the wrapped determinant of the deformation gradient for the Cauchy
    /// outputs, and the elastic strain E' = E + gamma n1 (taut: E, slack: eta E). In exact mode the wrinkled stress also
    /// differs from the linearised one of a twin fixture.
    static void SplitOutputsCheck(const bool exact)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        const SplitStateSpec specs[3] = { SplitWrinkledSpec(), SplitTautSpec(), SplitSlackSpec() };
        const real_t etas[3] = { 0.0, 0.0, 1e-3 };
        const char * names[3] = { "wrinkled", "taut", "slack" };

        for (index_t s = 0; s!=3; ++s)
        {
            gsInfo << "state " << names[s] << "\n";
            SplitFixture<2,false> fx(SplitMatNHComp,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,pi/6);
            fx.split->options().setSwitch("Exact",exact);
            fx.split->options().setInt("ExactLocalSolve",0);
            fx.split->options().setReal("Eta",etas[s]);
            fx.split->options().setReal("SlackMultiplier",0.0);
            fx.split->resetCounters();
            CHECK_EQUAL(specs[s].state,fx.displayedState(pt,0));

            const gsMatrix<real_t> S  = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
            const gsMatrix<real_t> Sst = fx.split->eval3D_stress(0,pt,zmat,MaterialOutput::StressN);
            SPLIT_CHECK_SAME(Sst,S);

            // Elastic strain E' of the state
            gsMatrix<real_t> Ep = fx.EengDef();
            if (specs[s].state==-1)
                Ep *= etas[s];
            else if (specs[s].state==0)
            {
                const real_t theta = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
                const real_t gamma = fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);
                Ep += gamma * SplitN1(theta);
            }

            const gsMatrix<real_t> P = fx.split->eval3D_pstrain(0,pt,zmat);
            CHECK_EQUAL(3,P.rows());
            CHECK_EQUAL(1,P.cols());
            if (P.rows()==3 && P.cols()==1)
            {
                CHECK(P(2,0) == 0);
                const gsMatrix<real_t> lamE = SplitPencilEigenvalues(SplitSym2(Ep(0,0),0.5*Ep(2,0),Ep(1,0)),fx.Gref);
                SplitCheckPrincipal("pstrain",P,lamE,1e-12);
            }

            const gsMatrix<real_t> Scon = SplitSym2(S(0,0),S(2,0),S(1,0));
            const gsMatrix<real_t> pS = fx.split->eval3D_pstress(0,pt,zmat,MaterialOutput::Generic);
            CHECK_EQUAL(2,pS.rows());
            CHECK_EQUAL(1,pS.cols());
            if (pS.rows()==2 && pS.cols()==1)
                SplitCheckPrincipal("pstress",pS,SplitPencilEigenvalues(Scon,fx.Gref.inverse()),1e-12);

            const real_t detF = fx.wrapped->eval3D_detF(0,pt,zmat,MaterialOutput::CauchyStressN)(0,0);
            gsInfo << "  detF = " << detF << "\n";
            CHECK(math::abs(detF-1) > 1e-6);
            const gsMatrix<real_t> Sc = S / detF;
            const gsMatrix<real_t> Cs = fx.split->eval3D_CauchyStress(0,pt,zmat,MaterialOutput::CauchyStressN);
            const gsMatrix<real_t> Cv = fx.split->eval3D_CauchyVector(0,pt,zmat,MaterialOutput::CauchyVectorN);
            SPLIT_CHECK_REL("CauchyStress",Cs,Sc,1e-14);
            SPLIT_CHECK_REL("CauchyVector",Cv,Sc,1e-14);

            const real_t detFp = fx.wrapped->eval3D_detF(0,pt,zmat,MaterialOutput::PCauchyStressN)(0,0);
            const gsMatrix<real_t> pC = fx.split->eval3D_CauchyPStress(0,pt,zmat,MaterialOutput::PCauchyStressN);
            CHECK_EQUAL(2,pC.rows());
            CHECK_EQUAL(1,pC.cols());
            if (pC.rows()==2 && pC.cols()==1)
            {
                const gsMatrix<real_t> lamC = SplitPencilEigenvalues(Scon/detFp,fx.Gdef.inverse());
                SplitCheckPrincipal("CauchyPStress",pC,lamC,1e-12);
                if (specs[s].state==0)
                {
                    CHECK(pC.allFinite());
                    gsInfo << "  Cauchy principal stresses: " << pC(0,0) << ", " << pC(1,0) << "\n";
                    CHECK(math::min(math::abs(pC(0,0)),math::abs(pC(1,0))) <= 1e-8 * pC.cwiseAbs().maxCoeff());
                }
            }

            if (exact && specs[s].state==0)
            {
                SplitFixture<2,false> lin(SplitMatNHComp,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,pi/6);
                const gsMatrix<real_t> Sl = lin.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
                gsInfo << "  exact vs linearised S_eff: rel.diff = " << SplitRelErr(S,Sl) << " (must exceed 1e-4)\n";
                CHECK(Sl.allFinite());
                CHECK(SplitRelErr(S,Sl) > 1e-4);
            }
            CHECK_EQUAL(0,fx.split->numExactFallbacks());
            CHECK_EQUAL(0,fx.split->numExactReclassified());
        }
    }

    /// Output functions in linearised mode (see SplitOutputsCheck)
    TEST(Outputs_linearised)
    {
        SplitOutputsCheck(false);
    }

    /// Output functions in exact mode (see SplitOutputsCheck)
    TEST(Outputs_exact)
    {
        SplitOutputsCheck(true);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef detail_tftsplit::WrinkleOptions<real_t>        SplitXOpt;
    typedef detail_tftsplit::WrinkleResult<real_t>         SplitXRes;
    typedef SplitFixture<2,false>::SplitType::SC_provider  SplitXProv;

    /// Tolerance of the exact solve at an iterate: max(k ||S'||, floor), the floor being the round-off of the
    /// strain E' mapped to a stress, 64 eps ||C'||_F max(1,||E'||). \a k is ExactTol or ExactStagTol.
    static real_t SplitExactTolAt(const gsMatrix<real_t> & Sp, const gsMatrix<real_t> & Cp, const gsMatrix<real_t> & Eel, const real_t k)
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        const real_t floorR = 64 * eps * Cp.norm() * math::max(real_t(1),Eel.norm());
        return math::max(k * Sp.norm(),floorR);
    }

    /// Full counter vector of the exact solve: the expected value of each outcome counter, and the two legacy sums
    template<class SplitT>
    static void SplitCheckExactCounters(const SplitT & split, const index_t nonconv, const index_t saddle, const index_t wrong,
                                        const index_t taut, const index_t slack, const index_t skipped)
    {
        CHECK_EQUAL(nonconv,split.numExactNonconverged());
        CHECK_EQUAL(saddle,split.numExactSaddle());
        CHECK_EQUAL(wrong,split.numExactWrongBranch());
        CHECK_EQUAL(taut,split.numExactReclassifiedTaut());
        CHECK_EQUAL(slack,split.numExactReclassifiedSlack());
        CHECK_EQUAL(skipped,split.numExactSkipped());
        CHECK_EQUAL(split.numExactNonconverged()+split.numExactSaddle()+split.numExactWrongBranch(),split.numExactFallbacks());
        CHECK_EQUAL(split.numExactReclassifiedTaut()+split.numExactReclassifiedSlack(),split.numExactReclassified());
    }

    /// \a r equals the reference result \a ref bit for bit (the tangent only if requested)
    static void SplitCheckSameResult(const SplitXRes & r, const SplitXRes & ref, const bool tangent)
    {
        SPLIT_CHECK_SAME(r.S,ref.S);
        SPLIT_CHECK_SAME(r.E,ref.E);
        if (tangent)
            SPLIT_CHECK_SAME(r.C,ref.C);
        else
            CHECK_EQUAL(0,r.C.size());
        CHECK(r.theta == ref.theta);
        CHECK(r.gamma == ref.gamma);
        CHECK(r.gammaHat == ref.gammaHat);
        CHECK(r.residual == ref.residual);
    }

    /// Stress-like vectors dual to the directions at \a theta, with their defining identities checked
    /// (c = cos theta, s = sin theta, n4 = (s^2,c^2,-2cs)):
    ///   v  = (c^2, s^2, cs):       n1.v = 1,  n2.v = 0,  n4.v = 0
    ///   w  = (-2sc, 2sc, c^2-s^2): n1.w = 0,  n2.w = 1
    ///   Tm = (s^2, c^2, -sc):      n1.Tm = n2.Tm = 0, n4.Tm = 1
    static void SplitDualVectors(const real_t theta, gsMatrix<real_t> & v, gsMatrix<real_t> & w, gsMatrix<real_t> & Tm)
    {
        const real_t c = math::cos(theta), s = math::sin(theta);
        v.resize(3,1);
        v << c*c, s*s, c*s;
        w.resize(3,1);
        w << -2*s*c, 2*s*c, c*c-s*s;
        Tm.resize(3,1);
        Tm << s*s, c*c, -s*c;
        gsVector<real_t,3> n4;
        n4 << s*s, c*c, -2*c*s;
        const gsVector<real_t,3> n1 = SplitN1(theta), n2 = SplitN2(theta);
        CHECK(math::abs(n1.dot(v.col(0)) - 1) <= 1e-15);
        CHECK(math::abs(n2.dot(v.col(0))) <= 1e-15);
        CHECK(math::abs(n4.dot(v.col(0))) <= 1e-15);
        CHECK(math::abs(n1.dot(w.col(0))) <= 1e-15);
        CHECK(math::abs(n2.dot(w.col(0)) - 1) <= 1e-15);
        CHECK(math::abs(n1.dot(Tm.col(0))) <= 1e-15);
        CHECK(math::abs(n2.dot(Tm.col(0))) <= 1e-15);
        CHECK(math::abs(n4.dot(Tm.col(0)) - 1) <= 1e-15);
    }

    /// Provider of a scalar-shaped problem along the direction \a theta of the input (S0,C0,E0):
    ///   t = v.(E' - E0),  S(E') = S0 + C0 (E' - E0) + h(t) v,  C(E') = C0 + h'(t) v v^T.
    /// S is the gradient of an added potential H(t) = int h, so the tangent is consistent and symmetric. Along
    /// E' = E0 + gamma n1(theta) one has t = gamma and g(gamma) = n1.S = b (gamma - gamma_lin) + h(gamma), with
    /// b = n1.C0 n1 and gamma_lin = -n1.S0/b. h moves only n1.S: n2.v = n4.v = 0.
    static SplitXProv SplitScalarProvider(const gsMatrix<real_t> & S0, const gsMatrix<real_t> & C0, const gsMatrix<real_t> & E0,
                                          const real_t theta, const std::function<real_t(real_t)> & h,
                                          const std::function<real_t(real_t)> & hp)
    {
        gsMatrix<real_t> v, w, Tm;
        SplitDualVectors(theta,v,w,Tm);
        return [=](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            const gsMatrix<real_t> d = Ep - E0;
            const real_t t = v(0,0)*d(0,0) + v(1,0)*d(1,0) + v(2,0)*d(2,0);
            Se = S0 + C0 * d + h(t) * v;
            Ce = C0 + hp(t) * (v * v.transpose());
        };
    }

    /// Provider of the wrapped material of \a fx at through-thickness coordinate \a z. \a fx must outlive it.
    static SplitXProv SplitWrappedProvider(const SplitFixture<2,false> & fx, const real_t z)
    {
        const SplitFixture<2,false> * f = &fx;
        const gsVector<real_t> pt = SplitPoint();
        return [f,pt,z](const gsMatrix<real_t> & Eel, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            const gsMatrix<real_t> Cm = f->CmatFromE(Eel);
            Se = f->wrapped->eval3D_vector_C(Cm,0,pt,z,MaterialOutput::VectorN);
            Ce = f->wrapped->eval3D_matrix_C(Cm,0,pt,z,MaterialOutput::MatrixA);
        };
    }

    /// The pointwise core at the engineering strain \a Eeng of the fixture (Gcov = Gref) with the wrapped material as
    /// provider; \a S3out receives the wrapped stress S(Eeng)
    static SplitXRes SplitPointwise(const SplitFixture<2,false> & fx, const gsMatrix<real_t> & Eeng, const real_t z,
                                    const SplitXOpt & opt, const bool tangent, gsMatrix<real_t> * S3out = NULL)
    {
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> Cmat = fx.CmatFromE(Eeng);
        const gsMatrix<real_t> S3 = fx.wrapped->eval3D_vector_C(Cmat,0,pt,z,MaterialOutput::VectorN);
        const gsMatrix<real_t> C3 = fx.wrapped->eval3D_matrix_C(Cmat,0,pt,z,MaterialOutput::MatrixA);
        if (S3out)
            *S3out = S3;
        return fx.split->wrinkle(S3,C3,Eeng,fx.Gref,opt,tangent,SplitFixture<2,false>::SplitType::dC_provider(),SplitWrappedProvider(fx,z));
    }

    /// Exact-mode options of the pointwise core: mixed criterion, eta = 0, no slack multiplier
    static SplitXOpt SplitExactOptions()
    {
        SplitXOpt opt;
        opt.tensionField = 2;
        opt.eta = 0;
        opt.slackMultiplier = 0;
        opt.exact = true;
        return opt;
    }

    /// Runs \a f \a N times in an OpenMP parallel loop (no output inside)
    template<class F>
    static void SplitRepeatParallel(const index_t N, const F & f)
    {
#       pragma omp parallel for
        for (index_t i = 0; i<N; ++i)
            f();
    }

    /// Input of the pointwise core with both stresses negative and Emax = 0: E = (0,-0.01,0) rotated by \a rot, S = C E,
    /// C isotropic, Gcov = I. The linearised wrinkled state has n4.S' = 0 in exact arithmetic (the sign is round-off),
    /// so a slack shift of the provider decides the acceptance order.
    struct SplitSlackInput
    {
        SplitFixture<2,false> fx;
        gsMatrix<real_t> I2, S, C, E, Tm, v, w;
        SplitXOpt optEx;
        SplitXRes lin;
        real_t thetaLin, gammaLin, m22;

        explicit SplitSlackInput(const real_t rot)
        :
        fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0), I2(SplitIdentity2())
        {
            C = SplitIsoC(0.3);
            const real_t cr = math::cos(rot), sr = math::sin(rot);
            E.resize(3,1);
            E << -0.01*sr*sr, -0.01*cr*cr, 0.02*cr*sr;
            S = C * E;
            optEx = SplitExactOptions();
            SplitXOpt optLin = optEx;
            optLin.exact = false;
            lin = fx.split->wrinkle(S,C,E,I2,optLin,true);
            thetaLin = lin.theta;
            gammaLin = lin.gamma;
            SplitDualVectors(thetaLin,v,w,Tm);
            const gsVector<real_t,3> n2 = SplitN2(thetaLin);
            const gsVector<real_t,3> Cn2 = C * n2;
            m22 = n2.dot(Cn2);
        }

        /// Affine provider with the wrinkle tension shifted by -alpha: S0 + C0 (E'-E0) - alpha Tm
        SplitXProv provider(const real_t alpha) const
        {
            const gsMatrix<real_t> S0 = S, C0 = C, E0 = E, Tm0 = Tm;
            return [=](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
            {
                Se = S0 + C0 * (Ep - E0) - alpha * Tm0;
                Ce = C0;
            };
        }

        SplitXRes run(const SplitXOpt & opt, const bool tangent, const SplitXProv & prov) const
        {
            return fx.split->wrinkle(S,C,E,I2,opt,tangent,SplitFixture<2,false>::SplitType::dC_provider(),prov);
        }
    };

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// (a) The plain Newton step from gamma_lin overshoots below zero; the bracketed pass still finds the root.
    /// The provider is g(gamma) = A atan((gamma - gamma*)/l): root gamma* = gamma_lin/2, g(gamma_lin) > 0, g(0) < 0,
    /// and a Newton step from gamma_lin lands at about -2.6 gamma_lin.
    TEST(Exact_bracket_overshoot)
    {
        SplitSynthetic sy(false);
        CHECK(!sy.lin.guardFallback);
        const gsVector<real_t,3> n1l = SplitN1(sy.thetaLin);
        const gsVector<real_t,3> Cn1l = sy.C * n1l;
        const real_t b = n1l.dot(Cn1l);
        const real_t gl = sy.gammaLin, gs = 0.5 * gl, l = 0.1 * gl, A = b * gl;
        const std::function<real_t(real_t)> h = [=](const real_t g) { return A * math::atan((g-gs)/l) - b * (g-gl); };
        const std::function<real_t(real_t)> hp = [=](const real_t g) { const real_t x = (g-gs)/l; return A / (l * (1 + x*x)) - b; };
        const SplitXProv prov = SplitScalarProvider(sy.S,sy.C,sy.E,sy.thetaLin,h,hp);

        // Precondition: the plain Newton step from gamma_lin leaves (0,gamma_lin), and the bracket closes
        gsMatrix<real_t> Sx, Cx;
        prov(sy.E + gl * n1l,Sx,Cx);
        const real_t g1 = n1l.dot(Sx.col(0));
        const gsVector<real_t,3> Cxn1 = Cx * n1l;
        const real_t a1 = n1l.dot(Cxn1);
        const real_t gammaN = gl - g1 / a1;
        prov(sy.E,Sx,Cx);
        const real_t g0 = n1l.dot(Sx.col(0));
        gsInfo << "  gamma_lin = " << gl << ", root = " << gs << ", plain Newton step from gamma_lin -> " << gammaN
               << " (" << gammaN / gl << " gamma_lin), g(gamma_lin) = " << g1 << ", g(0) = " << g0 << "\n";
        CHECK(g1 > 0);
        CHECK(g0 < 0);
        CHECK(gammaN < 0);

        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.run(sy.optEx,t==0,prov);
            gsInfo << "  tangent " << (t==0) << ": exactIterations = " << r.exactIterations << ", gamma/gamma* - 1 = "
                   << r.gamma / gs - 1 << "\n";
            CHECK_EQUAL(0,r.state);
            CHECK(!r.exactFallback);
            CHECK(!r.exactReclassified);
            CHECK(!r.exactUsed2x2);
            CHECK(!r.exactStagnated);
            CHECK(r.gamma > 0);
            CHECK(math::abs(r.gamma - gs) <= 1e-10 * gs);

            gsMatrix<real_t> Sp, Cp;
            prov(r.E,Sp,Cp);
            const gsVector<real_t,3> n1 = SplitN1(r.theta);
            const real_t R1 = n1.dot(Sp.col(0));
            const real_t tol = SplitExactTolAt(Sp,Cp,r.E,1e-12);
            gsInfo << "  |n1.S'| = " << math::abs(R1) << " (tol " << tol << ")\n";
            CHECK(math::abs(R1) <= tol);
            SPLIT_CHECK_REL("S vs provider at E'",r.S,Sp,1e-12);
            SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,0);
            CHECK_EQUAL(0,sy.fx.split->numExactStagnated());
            CHECK_EQUAL(0,sy.fx.split->numExactUsed2x2());
        }
    }

    /// Surface law of the near-singular hard state of the TFT tangent diagnostics (3D wrapped material on an embedded
    /// flat unit square, so that G = I at every point)
    struct SplitHardState3D
    {
        typedef gsMaterialMatrixTFTSplit<3,real_t,false> SplitType;

        gsMultiPatch<real_t> mp, mp_def;
        gsConstantFunction<real_t> t, E, nu;
        gsMaterialMatrixBase<real_t>::uPtr wrapped;
        memory::unique_ptr<SplitType> split;

        static gsMultiPatch<real_t> EmbeddedSquare()
        {
            gsMultiPatch<real_t> m;
            m.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
            m.embed(3);
            return m;
        }

        SplitHardState3D()
        :
        mp(EmbeddedSquare()), mp_def(mp), t(0.05e-3,3), E(1e9,3), nu(0.5,3)
        {
            std::vector<gsFunctionSet<real_t>*> parameters(2);
            parameters[0] = &E;
            parameters[1] = &nu;
            gsOptionList options;
            options.addInt("Material","Material model: (0): SvK | (1): NH | (3): MR",1);
            options.addInt("Implementation","Implementation: (1): Analytical",1);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            wrapped = getMaterialMatrix<3,real_t>(mp,t,parameters,options);
            split = memory::make_unique(new SplitType(wrapped.get()));
            split->setUndeformed(&mp);
            split->setDeformed(&mp_def);
        }

        /// Prescribed metric (1 + 2 E1, 1 + 2 E2, E3) of an engineering Voigt strain for G = I
        static gsMatrix<real_t> CmatFromE(const gsMatrix<real_t> & Eeng)
        {
            gsMatrix<real_t> Cm(3,1);
            Cm << 1 + 2*Eeng(0,0), 1 + 2*Eeng(1,0), Eeng(2,0);
            return Cm;
        }
    };

    /// (b) The near-singular state E 1e9, nu 0.5, thickness
    /// 0.05e-3, incompressible Neo-Hookean, E = (8.72, -0.50, 0) with lambda_2^2 ~ 2e-5. The wrapped stress S(E) is of order
    /// 1e16 and the tangent is near-singular; the exact solve is accepted without any counter, with a residual far below
    /// the round-off floor of the iterate.
    TEST(Exact_hard_state)
    {
        SplitHardState3D hs;
        typedef SplitHardState3D::SplitType HType;
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> I2 = SplitIdentity2();
        gsMatrix<real_t> E3(3,1);
        E3 << 8.7176634982689478, -0.49998982751041243, 0;
        const gsMatrix<real_t> Cmat = SplitHardState3D::CmatFromE(E3);
        const gsMatrix<real_t> S3 = hs.wrapped->eval3D_vector_C(Cmat,0,pt,0,MaterialOutput::VectorN);
        const gsMatrix<real_t> C3 = hs.wrapped->eval3D_matrix_C(Cmat,0,pt,0,MaterialOutput::MatrixA);
        const HType::SC_provider prov = [&](const gsMatrix<real_t> & Eel, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            const gsMatrix<real_t> Cm = SplitHardState3D::CmatFromE(Eel);
            Se = hs.wrapped->eval3D_vector_C(Cm,0,pt,0,MaterialOutput::VectorN);
            Ce = hs.wrapped->eval3D_matrix_C(Cm,0,pt,0,MaterialOutput::MatrixA);
        };

        SplitXOpt opt = SplitExactOptions();
        opt.exact = false;
        const SplitXRes lin = hs.split->wrinkle(S3,C3,E3,I2,opt,true,HType::dC_provider(),prov);
        gsInfo << "  linearised: state " << lin.state << ", tftFailure " << lin.tftFailure << ", guardFallback " << lin.guardFallback
               << ", Smin " << lin.Smin << ", Smax " << lin.Smax << ", Emin " << lin.Emin << ", Emax " << lin.Emax << "\n";
        CHECK_EQUAL(0,lin.state);
        CHECK(!lin.tftFailure);

        opt.exact = true;
        hs.split->resetCounters();
        const SplitXRes res = hs.split->wrinkle(S3,C3,E3,I2,opt,true,HType::dC_provider(),prov);

        // Quantities at the strain E' of the result; for a fallback also at the linearised E'
        auto report = [&](const char * label, const gsMatrix<real_t> & Ee, const real_t theta)
        {
            gsMatrix<real_t> Sp, Cp;
            prov(Ee,Sp,Cp);
            const gsVector<real_t,3> n1 = SplitN1(theta);
            gsInfo << std::setprecision(12) << "  " << label << ": ||S'|| = " << Sp.norm() << ", ||C'||_F = " << Cp.norm()
                   << ", floorR = " << SplitExactTolAt(Sp,Cp,Ee,0) << ", |R1| = " << math::abs(n1.dot(Sp.col(0)))
                   << std::setprecision(6) << "\n";
        };
        gsInfo << std::setprecision(12) << "  ||S(E)|| = " << S3.norm() << ", state " << res.state << ", exactFallback " << res.exactFallback
               << ", saddle " << res.exactSaddle << ", wrongBranch " << res.exactWrongBranch << ", reclassified " << res.exactReclassified
               << ", stagnated " << res.exactStagnated << ", used2x2 " << res.exactUsed2x2 << ", skipped " << res.exactSkipped
               << ", exactIterations " << res.exactIterations << ", residual " << res.residual << std::setprecision(6) << "\n";
        gsInfo << "  counters: nonconv " << hs.split->numExactNonconverged() << ", saddle " << hs.split->numExactSaddle()
               << ", wrong " << hs.split->numExactWrongBranch() << ", taut " << hs.split->numExactReclassifiedTaut()
               << ", slack " << hs.split->numExactReclassifiedSlack() << ", skipped " << hs.split->numExactSkipped() << "\n";
        report("at E'",res.E,res.theta);

        CHECK_EQUAL(0,res.state);
        CHECK(!res.exactFallback);
        CHECK(!res.exactReclassified);
        CHECK(!res.exactSkipped);
        SplitCheckExactCounters(*hs.split,0,0,0,0,0,0);
        gsMatrix<real_t> Sp, Cp;
        prov(res.E,Sp,Cp);
        const real_t R1 = math::abs(SplitN1(res.theta).dot(Sp.col(0)));
        CHECK(R1 <= 1e-9 * Sp.norm() || R1 <= SplitExactTolAt(Sp,Cp,res.E,0));
    }

    /// (c) Wrinkled states of tiny amplitude |E| = m in {1e-3,1e-5,1e-7}, a rotated over-contracted strain
    /// m R diag(1,-0.75) R^T: the solve is accepted without any counter and its residual is at the stall tolerance
    TEST(Exact_near_natural)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t ca = math::cos(pi/6), sa = math::sin(pi/6);
        const SplitMaterial mats[2] = { SplitMatNHIncomp, SplitMatMRIncomp };
        const real_t ms[3] = { 1e-3, 1e-5, 1e-7 };
        for (index_t im = 0; im!=2; ++im)
            for (index_t k = 0; k!=3; ++k)
            {
                const real_t m = ms[k];
                gsInfo << "material " << im << ", m = " << m << "\n";
                SplitFixture<2,false> fx(mats[im],SplitIdentity2(),1.0,1.0,0.0);
                fx.split->resetCounters();
                gsMatrix<real_t> E3(3,1);
                E3 << m * (ca*ca - 0.75*sa*sa), m * (sa*sa - 0.75*ca*ca), 2 * m * 1.75 * ca * sa;

                SplitXOpt opt = SplitExactOptions();
                opt.exact = false;
                const SplitXRes lin = SplitPointwise(fx,E3,0,opt,true);
                CHECK_EQUAL(0,lin.state);
                CHECK(!lin.guardFallback);

                opt.exact = true;
                const SplitXRes res = SplitPointwise(fx,E3,0,opt,true);
                const SplitXProv prov = SplitWrappedProvider(fx,0);
                gsMatrix<real_t> Sp, Cp;
                prov(res.E,Sp,Cp);
                const real_t tolS = SplitExactTolAt(Sp,Cp,res.E,1e-9);
                const real_t tolC = SplitExactTolAt(Sp,Cp,res.E,1e-12);
                gsInfo << "  exactIterations = " << res.exactIterations << ", stagnated = " << res.exactStagnated
                       << ", residual = " << res.residual << " (tolS " << tolS << "), floorR / (1e-12 ||S'||) = "
                       << SplitExactTolAt(Sp,Cp,res.E,0) / (1e-12 * Sp.norm()) << ", tolC = " << tolC << "\n";
                CHECK_EQUAL(0,res.state);
                CHECK(!res.exactFallback);
                CHECK(!res.exactReclassified);
                CHECK(!res.exactSkipped);
                CHECK(res.residual <= tolS);
                SplitCheckExactCounters(*fx.split,0,0,0,0,0,0);
            }
    }

    /// (d) Compressible Neo-Hookean wrinkled states, at the mid-surface and off it: no fallback, no reclassification
    TEST(Exact_NH_Comp_no_fallback)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t deltas[2] = { 0.05, 0.2 };
        const real_t zs[2] = { 0.0, 0.4 * 0.01 };
        for (index_t id = 0; id!=2; ++id)
            for (index_t iz = 0; iz!=2; ++iz)
            {
                gsInfo << "delta " << deltas[id] << ", z " << zs[iz] << "\n";
                SplitFixture<2,false> fx(SplitMatNHComp,SplitSkewedJref(),1.2,(1-deltas[id])/math::sqrt(1.2),pi/6);
                fx.split->resetCounters();
                const gsMatrix<real_t> E3 = fx.EengDef();
                SplitXOpt opt = SplitExactOptions();
                opt.exact = false;
                const SplitXRes lin = SplitPointwise(fx,E3,zs[iz],opt,true);
                CHECK_EQUAL(0,lin.state);

                opt.exact = true;
                const SplitXRes res = SplitPointwise(fx,E3,zs[iz],opt,true);
                gsMatrix<real_t> Sp, Cp;
                SplitWrappedProvider(fx,zs[iz])(res.E,Sp,Cp);
                const real_t tolS = SplitExactTolAt(Sp,Cp,res.E,1e-9);
                gsInfo << "  exactIterations = " << res.exactIterations << ", stagnated = " << res.exactStagnated
                       << ", residual = " << res.residual << " (tolS " << tolS << ")\n";
                CHECK_EQUAL(0,res.state);
                CHECK(!res.exactFallback);
                CHECK(!res.exactReclassified);
                CHECK(!res.exactSkipped);
                CHECK(res.residual <= tolS);
                SplitCheckExactCounters(*fx.split,0,0,0,0,0,0);
            }
    }

    /// (e) Isotropic incompressible laws keep the wrinkle direction fixed: the scalar pass alone solves n2.S' = 0
    TEST(Exact_scalar_only_isotropic)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitMaterial mats[2] = { SplitMatNHIncomp, SplitMatMRIncomp };
        const real_t deltas[2] = { 0.05, 0.2 };
        for (index_t im = 0; im!=2; ++im)
            for (index_t ij = 0; ij!=2; ++ij)
                for (index_t id = 0; id!=2; ++id)
                {
                    gsInfo << "material " << im << ", Jref " << (ij==0 ? "identity" : "skewed") << ", delta " << deltas[id] << "\n";
                    const gsMatrix<real_t> Jref = (ij==0) ? SplitIdentity2() : SplitSkewedJref();
                    SplitFixture<2,false> fx(mats[im],Jref,1.2,(1-deltas[id])/math::sqrt(1.2),pi/6);
                    fx.split->resetCounters();
                    const gsMatrix<real_t> E3 = fx.EengDef();
                    SplitXOpt opt = SplitExactOptions();
                    opt.exact = false;
                    gsMatrix<real_t> S3;
                    const SplitXRes lin = SplitPointwise(fx,E3,0,opt,true,&S3);
                    gsInfo << "  Smax - Smin = " << lin.Smax - lin.Smin << " (>= " << 1e-2 * S3.norm() << "), guardFallback " << lin.guardFallback << "\n";
                    CHECK(lin.Smax - lin.Smin >= 1e-2 * S3.norm());
                    CHECK(!lin.guardFallback);

                    opt.exact = true;
                    const SplitXRes res = SplitPointwise(fx,E3,0,opt,true);
                    gsMatrix<real_t> Sp, Cp;
                    SplitWrappedProvider(fx,0)(res.E,Sp,Cp);
                    const real_t tol = SplitExactTolAt(Sp,Cp,res.E,res.exactStagnated ? 1e-9 : 1e-12);
                    const real_t R2 = SplitN2(res.theta).dot(res.S.col(0));
                    gsInfo << "  exactUsed2x2 " << res.exactUsed2x2 << ", |n2.S'| = " << math::abs(R2) << " (tol " << tol << ")\n";
                    CHECK(!res.exactUsed2x2);
                    CHECK_EQUAL(0,fx.split->numExactUsed2x2());
                    CHECK(!res.exactFallback);
                    CHECK_EQUAL(0,res.state);
                    CHECK(math::abs(R2) <= tol);
                    CHECK(!res.exactSkipped);
                    CHECK(res.exactIterations >= 1);
                    SplitCheckExactCounters(*fx.split,0,0,0,0,0,0);
                }
    }

    /// (f) SvK major-branch point (law S = Csvk E', lambda 2, mu 1): the closed-form theta and gamma_lin
    /// come from the classifier input, whose stress has its principal entries swapped. The scalar pass converges to the
    /// stationary point gamma = 0.04 with tension n4.S' = -0.30 and det Jac_0 = -0.88, a saddle outside the slack band,
    /// not a slack point.
    TEST(Exact_M8_major_branch_saddle)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const gsMatrix<real_t> I2 = SplitIdentity2();
        gsMatrix<real_t> Csvk(3,3);
        Csvk << 4, 2, 0,
                2, 4, 0,
                0, 0, 1;
        const gsMatrix<real_t> & Cin = Csvk;
        gsMatrix<real_t> Sin(3,1), Ein(3,1);
        Sin << -0.38, -0.16, 0;
        Ein << 0.01, -0.1, 0;
        const SplitXProv prov = [Csvk](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            Se = Csvk * Ep;
            Ce = Csvk;
        };
        SplitXOpt optEx = SplitExactOptions();
        SplitXOpt optLin = optEx;
        optLin.exact = false;

        const SplitXRes lin = fx.split->wrinkle(Sin,Cin,Ein,I2,optLin,true);
        gsInfo << "  linearised: state " << lin.state << ", guardFallback " << lin.guardFallback << ", tftFailure "
               << lin.tftFailure << ", theta " << lin.theta << ", gamma_lin " << lin.gamma << "\n";
        CHECK_EQUAL(0,lin.state);
        CHECK(!lin.guardFallback);
        CHECK(SplitThetaDiff(lin.theta,0) <= 1e-14);

        // The stationary point gamma = 0.04, from the law alone
        const gsVector<real_t,3> n1 = SplitN1(lin.theta), n2 = SplitN2(lin.theta);
        gsVector<real_t,3> n4;
        n4 << math::sin(lin.theta)*math::sin(lin.theta), math::cos(lin.theta)*math::cos(lin.theta), -2*math::cos(lin.theta)*math::sin(lin.theta);
        const real_t gam = 0.04;
        const gsMatrix<real_t> Ep = Ein + gam * n1;
        const gsMatrix<real_t> Sp = Csvk * Ep;
        const gsVector<real_t,3> Cn1 = Csvk * n1, Cn2 = Csvk * n2;
        const real_t sw = n4.dot(Sp.col(0));
        const real_t M11 = n1.dot(Cn1), M12 = n1.dot(Cn2), M21 = n2.dot(Cn1), M22 = n2.dot(Cn2);
        const real_t d = sw + 2 * gam * M22;
        const real_t det0 = M11 * d - 2 * gam * M12 * M21;
        gsInfo << "  sw = " << sw << ", d = " << d << ", det0 = " << det0 << "\n";
        CHECK(math::abs(sw + 0.30) <= 1e-12);
        CHECK(math::abs(det0 + 0.88) <= 1e-12);
        CHECK(math::abs(n1.dot(Sp.col(0))) <= 1e-12);
        CHECK(math::abs(n2.dot(Sp.col(0))) <= 1e-12);

        for (index_t t = 0; t!=2; ++t)
        {
            fx.split->resetCounters();
            const SplitXRes r = fx.split->wrinkle(Sin,Cin,Ein,I2,optEx,t==0,SplitFixture<2,false>::SplitType::dC_provider(),prov);
            const SplitXRes ref = fx.split->wrinkle(Sin,Cin,Ein,I2,optLin,t==0);
            gsInfo << "  tangent " << (t==0) << ": exactFallback " << r.exactFallback << ", saddle " << r.exactSaddle << ", wrongBranch "
                   << r.exactWrongBranch << ", reclassified " << r.exactReclassified << ", used2x2 " << r.exactUsed2x2
                   << ", exactIterations " << r.exactIterations << "\n";
            CHECK_EQUAL(0,r.state);
            CHECK(r.exactFallback);
            CHECK(r.exactSaddle);
            CHECK(!r.exactWrongBranch);
            CHECK(!r.exactReclassified);
            CHECK(!r.exactUsed2x2);
            SplitCheckExactCounters(*fx.split,0,1,0,0,0,0);
            SplitCheckSameResult(r,ref,t==0);
        }
    }

    /// (g) nonconverged: the NaN provider ends the solve after one call
    TEST(Exact_counter_nonconverged)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::NaN);
        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.run(sy.optEx,t==0,prov);
            CHECK_EQUAL(0,r.state);
            CHECK(r.exactFallback);
            CHECK(!r.exactSaddle);
            CHECK(!r.exactWrongBranch);
            CHECK(!r.exactReclassified);
            CHECK(!r.exactSkipped);
            CHECK_EQUAL(1,r.exactIterations);
            SplitCheckExactCounters(*sy.fx.split,1,0,0,0,0,0);
            SplitCheckSameResult(r,(t==0) ? sy.lin : sy.linNoC,t==0);
        }
    }

    /// (g) saddle: the stationary point of the affine problem with a tangent C - beta w w^T that is not a minimum
    TEST(Exact_counter_saddle)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Saddle,sy.m22() + sy.swLin/sy.gammaLin);
        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.run(sy.optEx,t==0,prov);
            CHECK_EQUAL(0,r.state);
            CHECK(r.exactFallback);
            CHECK(r.exactSaddle);
            CHECK(!r.exactWrongBranch);
            CHECK(!r.exactReclassified);
            SplitCheckExactCounters(*sy.fx.split,0,1,0,0,0,0);
            SplitCheckSameResult(r,(t==0) ? sy.lin : sy.linNoC,t==0);
        }
    }

    /// (g) wrong branch: the slack shift alpha = sw_lin + gamma_lin M22 keeps the solution at
    /// (gamma_lin,theta_lin) with sw = -gamma_lin M22 < 0 and det Jac_0 = M11 gamma_lin M22 > 0, a minimum with a
    /// non-positive wrinkle tension outside the slack band. Larger shifts make det Jac_0 negative: a saddle.
    TEST(Exact_counter_wrong_branch)
    {
        SplitSynthetic sy(false);
        const gsVector<real_t,3> n1 = SplitN1(sy.thetaLin), n2 = SplitN2(sy.thetaLin);
        gsVector<real_t,3> n4;
        n4 << math::sin(sy.thetaLin)*math::sin(sy.thetaLin), math::cos(sy.thetaLin)*math::cos(sy.thetaLin),
              -2*math::cos(sy.thetaLin)*math::sin(sy.thetaLin);
        const real_t alphas[5] = { sy.swLin + sy.gammaLin * sy.m22(), 1.01 * sy.swLin, 1.1 * sy.swLin, 1.5 * sy.swLin, 2 * sy.swLin };
        for (index_t k = 0; k!=5; ++k)
        {
            const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Slack,alphas[k]);
            gsMatrix<real_t> Sx, Cx;
            prov(sy.E + sy.gammaLin * n1,Sx,Cx);
            const gsVector<real_t,3> Cn1 = Cx * n1, Cn2 = Cx * n2;
            const real_t sw = n4.dot(Sx.col(0));
            const real_t M11 = n1.dot(Cn1), M12 = n1.dot(Cn2), M21 = n2.dot(Cn1), M22 = n2.dot(Cn2);
            const real_t det0 = M11 * (sw + 2 * sy.gammaLin * M22) - 2 * sy.gammaLin * M12 * M21;
            gsInfo << "  alpha = " << alphas[k] << ": sw = " << sw << ", det0 = " << det0 << "\n";
            CHECK(sw < 0);
            CHECK(M11 > 0);
            // sw = sw_lin - alpha and d = sw + 2 gamma_lin M22: det0 = M11 d > 0 iff alpha < sw_lin + 2 gamma_lin M22
            const bool minimum = det0 > 0;
            CHECK_EQUAL(k < 3,minimum);

            for (index_t t = 0; t!=2; ++t)
            {
                sy.fx.split->resetCounters();
                const SplitXRes r = sy.run(sy.optEx,t==0,prov);
                CHECK_EQUAL(0,r.state);
                CHECK(r.exactFallback);
                CHECK_EQUAL(minimum,r.exactWrongBranch);
                CHECK_EQUAL(!minimum,r.exactSaddle);
                CHECK(!r.exactReclassified);
                if (minimum)
                    SplitCheckExactCounters(*sy.fx.split,0,0,1,0,0,0);
                else
                    SplitCheckExactCounters(*sy.fx.split,0,1,0,0,0,0);
                SplitCheckSameResult(r,(t==0) ? sy.lin : sy.linNoC,t==0);
            }
        }
    }

    /// (g) stagnated: the jump g(gamma) = b (gamma - gs) + delta sign(gamma - gs) has no root and |g| >= delta, with delta
    /// between ExactTol and ExactStagTol relative to ||S||. The best iterate is accepted at the stall tolerance.
    TEST(Exact_counter_stagnated)
    {
        SplitSynthetic sy(false);
        const gsVector<real_t,3> n1l = SplitN1(sy.thetaLin);
        const gsVector<real_t,3> Cn1l = sy.C * n1l;
        const real_t b = n1l.dot(Cn1l);
        const real_t gl = sy.gammaLin, gs = 1.5 * gl, dlt = 1e-10 * sy.sNorm;
        const std::function<real_t(real_t)> h = [=](const real_t g) { return b * (gl - gs) + dlt * (g > gs ? 1 : -1); };
        const std::function<real_t(real_t)> hp = [](const real_t) { return real_t(0); };
        const SplitXProv prov = SplitScalarProvider(sy.S,sy.C,sy.E,sy.thetaLin,h,hp);

        gsMatrix<real_t> Sx, Cx;
        prov(sy.E + gl * n1l,Sx,Cx);
        const real_t g1 = n1l.dot(Sx.col(0));
        gsInfo << "  g(gamma_lin) = " << g1 << " (< 0: the upper end of the bracket is open), jump delta = " << dlt << "\n";
        CHECK(g1 < 0);

        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.run(sy.optEx,t==0,prov);
            gsMatrix<real_t> Sp, Cp;
            prov(r.E,Sp,Cp);
            const real_t R1 = math::abs(SplitN1(r.theta).dot(Sp.col(0)));
            const real_t tolC = SplitExactTolAt(Sp,Cp,r.E,1e-12), tolS = SplitExactTolAt(Sp,Cp,r.E,1e-9);
            gsInfo << "  exactIterations = " << r.exactIterations << ", |n1.S'| = " << R1 << " (tolC " << tolC << ", tolS " << tolS << ")\n";
            CHECK_EQUAL(0,r.state);
            CHECK(r.exactStagnated);
            CHECK(!r.exactFallback);
            CHECK(!r.exactReclassified);
            CHECK(R1 <= tolS);
            CHECK(R1 > tolC);
            SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,0);
            CHECK_EQUAL(1,sy.fx.split->numExactStagnated());
        }
    }

    /// (g) taut at gamma = 0: the provider -S + C (E'-E) has g(gamma_lin) > 0 and g(0) > 0
    TEST(Exact_counter_taut_at_zero)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Taut);
        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.run(sy.optEx,t==0,prov);
            CHECK_EQUAL(1,r.state);
            CHECK(r.exactReclassified);
            CHECK(!r.exactFallback);
            CHECK_EQUAL(2,r.exactIterations);
            SplitCheckExactCounters(*sy.fx.split,0,0,0,1,0,0);
            SPLIT_CHECK_SAME(r.S,sy.S);
            SPLIT_CHECK_SAME(r.E,sy.E);
            if (t==0)
                SPLIT_CHECK_SAME(r.C,sy.C);
            else
                CHECK_EQUAL(0,r.C.size());
            CHECK(std::isnan(r.theta) && std::isnan(r.gamma) && std::isnan(r.gammaHat) && std::isnan(r.residual));
        }
    }

    /// (g) gamma_lin <= 0 (strain-based criterion): one provider call at gamma = 0 decides taut or not converged
    TEST(Exact_counter_gamma_lin_nonpositive)
    {
        SplitSynthetic sy(false);
        gsMatrix<real_t> C = SplitIsoC(0.3), E(3,1);
        E << 0.02, -0.001, 0;
        const gsMatrix<real_t> S = C * E;
        SplitXOpt optEx = sy.optEx;
        optEx.tensionField = 0;
        SplitXOpt optLin = optEx;
        optLin.exact = false;

        const SplitXRes linT = sy.fx.split->wrinkle(S,C,E,sy.I2,optLin,true);
        const SplitXRes linN = sy.fx.split->wrinkle(S,C,E,sy.I2,optLin,false);
        gsInfo << "  linearised: state " << linT.state << ", gamma_lin " << linT.gamma << ", guardFallback " << linT.guardFallback << "\n";
        CHECK_EQUAL(0,linT.state);
        CHECK(linT.gamma < 0);
        CHECK(!linT.guardFallback);

        const SplitXProv affine = [C,E,S](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            Se = S + C * (Ep - E);
            Ce = C;
        };
        const SplitXProv negated = [C,E,S](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            Se = -S + C * (Ep - E);
            Ce = C;
        };
        for (index_t t = 0; t!=2; ++t)
        {
            // g(0) = n1.S > 0: taut at zero
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.fx.split->wrinkle(S,C,E,sy.I2,optEx,t==0,SplitFixture<2,false>::SplitType::dC_provider(),affine);
            CHECK_EQUAL(1,r.state);
            CHECK(r.exactReclassified);
            CHECK(!r.exactFallback);
            CHECK_EQUAL(1,r.exactIterations);
            SplitCheckExactCounters(*sy.fx.split,0,0,0,1,0,0);
            SPLIT_CHECK_SAME(r.S,S);
            SPLIT_CHECK_SAME(r.E,E);
            if (t==0)
                SPLIT_CHECK_SAME(r.C,C);
            else
                CHECK_EQUAL(0,r.C.size());
            CHECK(std::isnan(r.theta) && std::isnan(r.gamma) && std::isnan(r.gammaHat) && std::isnan(r.residual));

            // g(0) = -n1.S < 0: not converged
            sy.fx.split->resetCounters();
            const SplitXRes rn = sy.fx.split->wrinkle(S,C,E,sy.I2,optEx,t==0,SplitFixture<2,false>::SplitType::dC_provider(),negated);
            CHECK_EQUAL(0,rn.state);
            CHECK(rn.exactFallback);
            CHECK(!rn.exactReclassified);
            CHECK_EQUAL(1,rn.exactIterations);
            SplitCheckExactCounters(*sy.fx.split,1,0,0,0,0,0);
            SplitCheckSameResult(rn,(t==0) ? linT : linN,t==0);
        }
    }

    /// (g) slack reclassification: a negative wrinkle tension within the slack band, as a
    /// minimum (alpha = gamma_lin M22) and as a non-minimum (alpha = 3 gamma_lin M22), returns sigma (S,E,C)
    TEST(Exact_counter_slack_reclassification)
    {
        SplitSlackInput si(0.0);
        const gsVector<real_t,3> n1 = SplitN1(si.thetaLin), n2 = SplitN2(si.thetaLin);
        gsVector<real_t,3> n4;
        n4 << math::sin(si.thetaLin)*math::sin(si.thetaLin), math::cos(si.thetaLin)*math::cos(si.thetaLin),
              -2*math::cos(si.thetaLin)*math::sin(si.thetaLin);
        gsInfo << "  linearised: state " << si.lin.state << ", guardFallback " << si.lin.guardFallback << ", Emax " << si.lin.Emax
               << ", Emin " << si.lin.Emin << ", theta " << si.thetaLin << ", gamma_lin " << si.gammaLin << ", M22 " << si.m22 << "\n";
        CHECK_EQUAL(0,si.lin.state);
        CHECK(!si.lin.guardFallback);
        CHECK(si.lin.Emax <= 1e-8 * math::max(math::abs(si.lin.Emin),math::abs(si.lin.Emax)));

        const real_t alphas[2] = { si.gammaLin * si.m22, 3 * si.gammaLin * si.m22 };
        const real_t etas[2] = { 1e-3, 1e-3 };
        const real_t mults[2] = { 0.0, 0.25 };
        for (index_t k = 0; k!=2; ++k)
        {
            const SplitXProv prov = si.provider(alphas[k]);
            gsMatrix<real_t> Sx, Cx;
            prov(si.E + si.gammaLin * n1,Sx,Cx);
            const gsVector<real_t,3> Cn1 = Cx * n1, Cn2 = Cx * n2;
            const real_t sw = n4.dot(Sx.col(0));
            const real_t M11 = n1.dot(Cn1), M12 = n1.dot(Cn2), M21 = n2.dot(Cn1), M22 = n2.dot(Cn2);
            const real_t det0 = M11 * (sw + 2 * si.gammaLin * M22) - 2 * si.gammaLin * M12 * M21;
            gsInfo << "  route " << (k==0 ? "5d" : "5c") << ": alpha = " << alphas[k] << ", sw = " << sw << ", det0 = " << det0 << "\n";
            CHECK(sw < 0);
            CHECK_EQUAL(k==0,det0 > 0);

            for (index_t i = 0; i!=2; ++i)
                for (index_t t = 0; t!=2; ++t)
                {
                    SplitXOpt opt = si.optEx;
                    opt.eta = etas[i];
                    opt.slackMultiplier = mults[i];
                    const real_t sigma = (mults[i] > 0) ? mults[i] : etas[i];
                    si.fx.split->resetCounters();
                    const SplitXRes r = si.run(opt,t==0,prov);
                    CHECK_EQUAL(-1,r.state);
                    CHECK(r.exactReclassified);
                    CHECK(!r.exactFallback);
                    SplitCheckExactCounters(*si.fx.split,0,0,0,0,1,0);
                    const gsMatrix<real_t> Sexp = sigma * si.S, Eexp = sigma * si.E, Cexp = sigma * si.C;
                    SPLIT_CHECK_SAME(r.S,Sexp);
                    SPLIT_CHECK_SAME(r.E,Eexp);
                    if (t==0)
                        SPLIT_CHECK_SAME(r.C,Cexp);
                    else
                        CHECK_EQUAL(0,r.C.size());
                    CHECK(std::isnan(r.theta) && std::isnan(r.gamma) && std::isnan(r.gammaHat) && std::isnan(r.residual));
                }
        }
    }

    /// (g) skipped: an indefinite wrapped tangent along n1 is not attempted, and is not a fallback
    TEST(Exact_counter_skipped)
    {
        SplitSynthetic sy(false);
        const gsMatrix<real_t> negC = -sy.C;
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Affine);
        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.fx.split->wrinkle(sy.S,negC,sy.E,sy.I2,sy.optEx,t==0,SplitSynthetic::SplitType::dC_provider(),prov);
            CHECK(r.exactSkipped);
            CHECK(r.tftFailure);
            CHECK_EQUAL(0,r.exactIterations);
            CHECK(!r.exactFallback);
            CHECK(!r.exactReclassified);
            SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,1);
            CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());
        }
    }

    /// (g) 2x2 pass: a constant shift of n2.S' above the tolerance leaves the scalar pass unsolved
    TEST(Exact_counter_used2x2)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Switch,1e-4);
        for (index_t t = 0; t!=2; ++t)
        {
            sy.fx.split->resetCounters();
            const SplitXRes r = sy.run(sy.optEx,t==0,prov);
            gsInfo << "  exactIterations = " << r.exactIterations << ", stagnated = " << r.exactStagnated << "\n";
            CHECK_EQUAL(0,r.state);
            CHECK(r.exactUsed2x2);
            CHECK(!r.exactFallback);
            CHECK(!r.exactReclassified);
            SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,0);
            CHECK_EQUAL(1,sy.fx.split->numExactUsed2x2());
        }
    }

    /// Anisotropic input for which the isotropy guard rejects the closed form while the minor stress -r is inside or
    /// outside the taut band (1e-8 relative to the larger principal stress). With C = [[1,0.2,0],[0.2,1,c12],[0,c12,200]]
    /// and S = (1,-r,0) the minor direction is theta = pi/2, b = 1, gamma_lin = r and the closed-form residual is
    /// -c12 r. The provider -S + C (E'-E) has its 2x2 root at gamma <= 0, which is taut within the
    /// band and a wrong branch outside it.
    static void SplitBandTautCase(const real_t r, const real_t c12, const bool expectTaut)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const gsMatrix<real_t> I2 = SplitIdentity2();
        gsMatrix<real_t> C(3,3), S(3,1);
        C << 1, 0.2, 0,
             0.2, 1, c12,
             0, c12, 200;
        S << 1, -r, 0;
        const gsMatrix<real_t> E = C.inverse() * S;
        const SplitXProv prov = [C,E,S](const gsMatrix<real_t> & Ep, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
        {
            Se = -S + C * (Ep - E);
            Ce = C;
        };
        SplitXOpt optEx = SplitExactOptions();
        SplitXOpt optLin = optEx;
        optLin.exact = false;

        const SplitXRes lin = fx.split->wrinkle(S,C,E,I2,optLin,true);
        gsInfo << "  r " << r << ", c12 " << c12 << ": linearised state " << lin.state << ", guardFallback " << lin.guardFallback << ", tftFailure "
               << lin.tftFailure << ", Smin/max|S| " << lin.Smin / math::max(math::abs(lin.Smin),math::abs(lin.Smax)) << ", Emax " << lin.Emax << "\n";
        CHECK_EQUAL(0,lin.state);
        CHECK(lin.guardFallback);
        CHECK(!lin.tftFailure);
        const real_t rel = lin.Smin / math::max(math::abs(lin.Smin),math::abs(lin.Smax));
        CHECK_EQUAL(expectTaut,rel >= -1e-8);

        for (index_t t = 0; t!=2; ++t)
        {
            fx.split->resetCounters();
            const SplitXRes res = fx.split->wrinkle(S,C,E,I2,optEx,t==0,SplitFixture<2,false>::SplitType::dC_provider(),prov);
            const SplitXRes ref = fx.split->wrinkle(S,C,E,I2,optLin,t==0);
            gsInfo << "  tangent " << (t==0) << ": state " << res.state << ", exactUsed2x2 " << res.exactUsed2x2 << ", reclassified "
                   << res.exactReclassified << ", wrongBranch " << res.exactWrongBranch << ", saddle " << res.exactSaddle << ", fallback "
                   << res.exactFallback << ", exactIterations " << res.exactIterations << "\n";
            CHECK(res.exactUsed2x2);
            CHECK_EQUAL(1,fx.split->numExactUsed2x2());
            if (expectTaut)
            {
                CHECK_EQUAL(1,res.state);
                CHECK(res.exactReclassified);
                CHECK(!res.exactFallback);
                SplitCheckExactCounters(*fx.split,0,0,0,1,0,0);
            }
            else
            {
                CHECK_EQUAL(0,res.state);
                CHECK(res.exactFallback);
                CHECK(res.exactWrongBranch);
                CHECK(!res.exactReclassified);
                SplitCheckExactCounters(*fx.split,0,0,1,0,0,0);
                SplitCheckSameResult(res,ref,t==0);
            }
        }
    }

    /// (g) band-taut: gamma <= 0 with the minor stress inside the taut band is taut
    TEST(Exact_counter_band_taut)
    {
        SplitBandTautCase(0.5e-8,10,true);
    }

    /// (g) gamma <= 0 with the minor stress outside the taut band is a wrong branch
    TEST(Exact_counter_band_taut_wrong_branch)
    {
        SplitBandTautCase(1e-6,10,false);
    }

    /// (h) The counters are incremented atomically: one OpenMP loop per counter, each call ending in the same outcome
    TEST(OpenMP_exact_new_counters_atomic)
    {
        SplitSynthetic sy(false);
        SplitSlackInput si(0.0);
        const index_t N = 20000;
        gsInfo << "  omp_get_max_threads() = " << omp_get_max_threads() << "\n";

        const SplitSynthetic::Provider nan = sy.provider(SplitSynthetic::NaN);
        const SplitSynthetic::Provider saddle = sy.provider(SplitSynthetic::Saddle,sy.m22() + sy.swLin/sy.gammaLin);
        const SplitSynthetic::Provider wrong = sy.provider(SplitSynthetic::Slack,sy.swLin + sy.gammaLin * sy.m22());
        const SplitSynthetic::Provider taut = sy.provider(SplitSynthetic::Taut);
        const SplitSynthetic::Provider affine = sy.provider(SplitSynthetic::Affine);
        const SplitSynthetic::Provider sw2x2 = sy.provider(SplitSynthetic::Switch,1e-4);
        const gsMatrix<real_t> negC = -sy.C;

        const gsVector<real_t,3> n1l = SplitN1(sy.thetaLin);
        const gsVector<real_t,3> Cn1l = sy.C * n1l;
        const real_t b = n1l.dot(Cn1l);
        const real_t gl = sy.gammaLin, gs = 1.5 * gl, dlt = 1e-10 * sy.sNorm;
        const std::function<real_t(real_t)> h = [=](const real_t g) { return b * (gl - gs) + dlt * (g > gs ? 1 : -1); };
        const std::function<real_t(real_t)> hp = [](const real_t) { return real_t(0); };
        const SplitXProv jump = SplitScalarProvider(sy.S,sy.C,sy.E,sy.thetaLin,h,hp);

        SplitXOpt slackOpt = si.optEx;
        slackOpt.eta = 1e-3;
        slackOpt.slackMultiplier = 0.25;
        const SplitXProv slack = si.provider(si.gammaLin * si.m22);

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.run(sy.optEx,false,nan); });
        SplitCheckExactCounters(*sy.fx.split,N,0,0,0,0,0);

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.run(sy.optEx,false,saddle); });
        SplitCheckExactCounters(*sy.fx.split,0,N,0,0,0,0);

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.run(sy.optEx,false,wrong); });
        SplitCheckExactCounters(*sy.fx.split,0,0,N,0,0,0);

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.run(sy.optEx,false,taut); });
        SplitCheckExactCounters(*sy.fx.split,0,0,0,N,0,0);

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.fx.split->wrinkle(sy.S,negC,sy.E,sy.I2,sy.optEx,false,SplitSynthetic::SplitType::dC_provider(),affine); });
        SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,N);

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.run(sy.optEx,false,sw2x2); });
        SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,0);
        CHECK_EQUAL(N,sy.fx.split->numExactUsed2x2());

        sy.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { sy.run(sy.optEx,false,jump); });
        SplitCheckExactCounters(*sy.fx.split,0,0,0,0,0,0);
        CHECK_EQUAL(N,sy.fx.split->numExactStagnated());

        si.fx.split->resetCounters();
        SplitRepeatParallel(N,[&]() { si.run(slackOpt,false,slack); });
        SplitCheckExactCounters(*si.fx.split,0,0,0,0,N,0);
    }

    /// (i) The exact tangent against an explicit central finite difference of the exact stress with step 1e-6; the
    /// solve is polished to round-off, so the stress carries no tolerance-sized noise of order tol/h
    TEST(Exact_FD_C_h1e6)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t h = 1e-6;
        const SplitMaterial mats[2] = { SplitMatNHIncomp, SplitMatMRIncomp };
        const gsVector<real_t> pt = SplitPoint();
        for (index_t im = 0; im!=2; ++im)
        {
            gsInfo << "material " << im << "\n";
            SplitFixture<2,false> fx(mats[im],SplitSkewedJref(),1.2,0.8/math::sqrt(1.2),pi/6);
            fx.split->options().setSwitch("Exact",true);
            fx.split->options().setInt("ExactLocalSolve",0);
            fx.split->options().setReal("Eta",0);
            fx.split->resetCounters();
            const gsMatrix<real_t> Eeng = fx.EengDef();

            gsMatrix<real_t> J_FD(3,3);
            for (index_t j = 0; j!=3; ++j)
            {
                gsMatrix<real_t> Ep = Eeng, Em = Eeng;
                Ep(j,0) += h;
                Em(j,0) -= h;
                const gsMatrix<real_t> Cp = fx.CmatFromE(Ep), Cm = fx.CmatFromE(Em);
                CHECK_EQUAL(0,fx.stateAt(Cp,pt,0));
                CHECK_EQUAL(0,fx.stateAt(Cm,pt,0));
                const gsMatrix<real_t> Sp = fx.split->eval3D_vector_C(Cp,0,pt,0,MaterialOutput::Generic);
                const gsMatrix<real_t> Sm = fx.split->eval3D_vector_C(Cm,0,pt,0,MaterialOutput::Generic);
                J_FD.col(j) = (Sp - Sm) / (2*h);

                gsMatrix<real_t> S3p, S3m;
                const SplitXRes rp = SplitPointwise(fx,Ep,0,SplitExactOptions(),false,&S3p);
                const SplitXRes rm = SplitPointwise(fx,Em,0,SplitExactOptions(),false,&S3m);
                CHECK(rp.exactIterations >= 1);
                CHECK(rm.exactIterations >= 1);
                gsInfo << "  stencil e_" << j << ": residual / ||S||: " << rp.residual / S3p.norm() << " (+), " << rm.residual / S3m.norm() << " (-)\n";
            }
            const gsMatrix<real_t> Ceff = fx.split->eval3D_matrix_C(fx.CmatDef(),0,pt,0,MaterialOutput::Generic).reshape(3,3);
            SPLIT_CHECK_REL("exact C vs FD h 1e-6",Ceff,J_FD,1e-7);
            CHECK_EQUAL(0,fx.split->numExactFallbacks());
            CHECK_EQUAL(0,fx.split->numExactReclassified());
        }
    }

    /// (j) The exact tangent annihilates the wrinkle direction: C_t n1 = 0 and C_t^T n1 = 0 relative to ||C_t||
    TEST(Exact_Ct_n1_zero)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const SplitStateSpec w = SplitWrinkledSpec();
        const SplitMaterial mats[2] = { SplitMatNHIncomp, SplitMatMRIncomp };
        for (index_t im = 0; im!=2; ++im)
            for (index_t ij = 0; ij!=2; ++ij)
            {
                gsInfo << "material " << im << ", Jref " << (ij==0 ? "identity" : "skewed") << "\n";
                SplitFixture<2,false> fx(mats[im],(ij==0) ? SplitIdentity2() : SplitSkewedJref(),w.lam1,w.lam2,pi/6);
                fx.split->resetCounters();
                const SplitXRes res = SplitPointwise(fx,fx.EengDef(),0,SplitExactOptions(),true);
                CHECK_EQUAL(0,res.state);
                CHECK(!res.exactFallback);
                CHECK(!res.exactReclassified);
                CHECK(!res.exactSkipped);
                CHECK(res.exactIterations >= 1);
                SplitCheckExactCounters(*fx.split,0,0,0,0,0,0);
                CHECK_EQUAL(3,res.C.rows());
                CHECK_EQUAL(3,res.C.cols());
                gsMatrix<real_t> Sp, Cp;
                SplitWrappedProvider(fx,0)(res.E,Sp,Cp);
                const gsVector<real_t,3> n1 = SplitN1(res.theta);
                const gsVector<real_t,3> Cn1 = res.C * n1;
                const gsVector<real_t,3> CTn1 = res.C.transpose() * n1;
                const real_t Cnorm = res.C.norm();
                gsInfo << "  ||C'||_F / ||C_t||_F = " << Cp.norm() / Cnorm << ", ||C_t n1|| / ||C_t|| = " << Cn1.norm() / Cnorm
                       << ", ||C_t^T n1|| / ||C_t|| = " << CTn1.norm() / Cnorm << "\n";
                CHECK(Cn1.norm() <= 1e-12 * Cnorm);
                CHECK(CTn1.norm() <= 1e-12 * Cnorm);
            }
    }

    /// (j) Continuity of the exact stress as the wrinkle amplitude gamma -> 0 (wrinkled -> taut): the distance to the
    /// wrapped stress decays linearly in delta
    TEST(Exact_continuity_boundaries)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> zmat = SplitZ(0);
        {
            const real_t deltas[4] = { 1e-2, 1e-4, 1e-6, 1e-8 };
            real_t slope0 = 0;
            for (index_t i = 0; i!=4; ++i)
            {
                SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),1.2,(1-deltas[i])/math::sqrt(1.2),pi/6);
                fx.split->options().setSwitch("Exact",true);
                fx.split->options().setInt("ExactLocalSolve",0);
                fx.split->options().setReal("Eta",0);
                fx.split->options().setReal("SlackMultiplier",0);
                fx.split->resetCounters();
                CHECK_EQUAL(0,fx.displayedState(pt,0));
                const gsMatrix<real_t> Sex = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
                const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
                const real_t d = (Sex - Sw).norm() / Sw.norm();
                gsInfo << "  taut boundary, delta " << deltas[i] << ": d = " << d << ", d / delta = " << d / deltas[i] << "\n";
                if (i==0)
                    slope0 = d / deltas[i];
                else
                {
                    CHECK(d / deltas[i] <= 2 * slope0);
                    CHECK(d / deltas[i] >= 0.5 * slope0);
                }
                CHECK_EQUAL(0,fx.split->numExactFallbacks());
                CHECK_EQUAL(0,fx.split->numExactReclassified());
            }
        }
        {
            // Wrinkled -> slack: the strain along the wrinkles Emax = eps -> 0 and so does the whole stress
            const real_t epss[3] = { 1e-2, 1e-4, 1e-6 };
            real_t slope0 = 0;
            for (index_t i = 0; i!=3; ++i)
            {
                SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),math::sqrt(1 + 2*epss[i]),0.9,pi/6);
                fx.split->options().setSwitch("Exact",true);
                fx.split->options().setInt("ExactLocalSolve",0);
                fx.split->options().setReal("Eta",0);
                fx.split->options().setReal("SlackMultiplier",0);
                fx.split->resetCounters();
                CHECK_EQUAL(0,fx.displayedState(pt,0));
                const gsMatrix<real_t> Sex = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
                const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
                const real_t d = Sex.norm() / Sw.norm();
                gsInfo << "  slack boundary, eps " << epss[i] << ": d = " << d << ", d / eps = " << d / epss[i] << "\n";
                if (i==0)
                    slope0 = d / epss[i];
                else
                {
                    CHECK(d / epss[i] <= 2 * slope0);
                    CHECK(d / epss[i] >= 0.5 * slope0);
                }
                CHECK_EQUAL(0,fx.split->numExactFallbacks());
                CHECK_EQUAL(0,fx.split->numExactReclassified());
            }
        }
    }

    /// (k) The exact solve rejects options with ExactStagTol < ExactTol, ExactTol <= 0, ExactMaxIt < 1 or NaN, but only
    /// where it is used; the option list carries the documented defaults
    TEST(Exact_option_validation)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Affine);
        const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
        auto badOptions = [&](const index_t k)
        {
            SplitXOpt o = sy.optEx;
            switch (k)
            {
            case 0: o.exactTol = nan; break;
            case 1: o.exactStagTol = nan; break;
            case 2: o.exactTol = 1e-10; o.exactStagTol = 1e-12; break;
            case 3: o.exactTol = 0; break;
            case 4: o.exactTol = -1e-12; break;
            case 5: o.exactMaxIt = 0; break;
            default: o.exactMaxIt = -1; break;
            }
            return o;
        };
        auto runs = [&](const SplitXOpt & o, const gsMatrix<real_t> & S, const gsMatrix<real_t> & E, const bool tangent)
        {
            bool thrown = false;
            try { sy.fx.split->wrinkle(S,sy.C,E,sy.I2,o,tangent,SplitSynthetic::SplitType::dC_provider(),prov); }
            catch (...) { thrown = true; }
            return !thrown;
        };

        for (index_t k = 0; k!=7; ++k)
        {
            const SplitXOpt o = badOptions(k);
            CHECK_THROW(sy.run(o,true,prov), std::runtime_error);
            CHECK_THROW(sy.run(o,false,prov), std::runtime_error);

            // Not in the linearised mode
            SplitXOpt ol = o;
            ol.exact = false;
            CHECK(runs(ol,sy.S,sy.E,true));

            // Not at a taut point, which returns before the exact branch
            gsMatrix<real_t> E(3,1);
            E << 0.02, 0.01, 0.0;
            const gsMatrix<real_t> S = sy.C * E;
            CHECK(runs(o,S,E,true));
            CHECK_EQUAL(1,sy.fx.split->wrinkle(S,sy.C,E,sy.I2,o,true).state);
        }

        // The boundary values are valid
        SplitXOpt eq = sy.optEx;
        eq.exactTol = 1e-12;
        eq.exactStagTol = 1e-12;
        CHECK(runs(eq,sy.S,sy.E,true));
        SplitXOpt one = sy.optEx;
        one.exactMaxIt = 1;
        CHECK(runs(one,sy.S,sy.E,true));
        CHECK(runs(one,sy.S,sy.E,false));

        // Option list defaults and the options read by the pointwise core
        gsOptionList & ol = sy.fx.split->options();
        CHECK(ol.getReal("ExactTol") == real_t(1e-12));
        CHECK(ol.getReal("ExactStagTol") == real_t(1e-9));
        CHECK_EQUAL(50,ol.getInt("ExactMaxIt"));
        const SplitXOpt wd = sy.fx.split->wrinkleOptions();
        CHECK(wd.exactTol == real_t(1e-12));
        CHECK(wd.exactStagTol == real_t(1e-9));
        CHECK_EQUAL(50,wd.exactMaxIt);
        const SplitXOpt def;
        CHECK(def.exactTol == real_t(1e-12));
        CHECK(def.exactStagTol == real_t(1e-9));
        CHECK_EQUAL(50,def.exactMaxIt);

        ol.setReal("ExactTol",1e-10);
        ol.setReal("ExactStagTol",1e-8);
        ol.setInt("ExactMaxIt",7);
        const SplitXOpt ws = sy.fx.split->wrinkleOptions();
        CHECK(ws.exactTol == real_t(1e-10));
        CHECK(ws.exactStagTol == real_t(1e-8));
        CHECK_EQUAL(7,ws.exactMaxIt);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef SplitFixture<2,false>::SplitType::W_provider   SplitXWProv;
    typedef SplitFixture<2,false>::SplitType::dC_provider  SplitXdCProv;

    /**
     * Hyperelastic laws with closed-form W, S = dW/dE and C = dS/dE in engineering Voigt notation
     * E = (E_11, E_22, 2 E_12), S = (S^11, S^22, S^12). The evaluation allocates nothing for W, so that the
     * brute-force searches stay cheap. A non-finite intermediate result gives NaN entries.
     *
     * Ortho: W = 1/2 E^T D E, S = D E, C = D with D = Rot^T D0 Rot, D0 = [[10,3,0],[3,2,0],[0,0,1.5]]
     * rotated by phi0 (symmetric positive definite, quadratic).
     *
     * NHFibre: incompressible plane-stress neo-Hookean (mu = 1) on the reference metric G plus an exponential
     * fibre. With g = G + 2 E_t the deformed covariant metric, r = det G / det g:
     *   W_NH = mu/2 (tr(G^{-1} g) + r - 3),  S^{ab} = mu (G^{ab} - r g^{ab}),
     *   C^{abcd} = mu r (2 g^{ab} g^{cd} + g^{ac} g^{bd} + g^{ad} g^{bc}),
     * and, with the G-unit fibre a (a^T G a = 1), I4 = a^T g a, x = I4 - 1, A = (a_1^2, a_2^2, a_1 a_2),
     *   W_f = k1/(2 k2) (exp(k2 x^2) - 1),  S_f = 2 k1 x exp(k2 x^2) A,
     *   C_f = 4 k1 exp(k2 x^2) (1 + 2 k2 x^2) A A^T.
     * The law is strictly convex in E where g is positive definite, with W(0) = 0 and S(0) = 0.
     */
    struct SplitAniso
    {
        enum Kind { Ortho, NHFibre };

        Kind kind;
        real_t D[9];
        real_t G00, G01, G11, detG, Gi00, Gi01, Gi11;
        real_t mu, k1, k2, a0, a1;

        SplitAniso() : kind(Ortho), G00(1), G01(0), G11(1), detG(1), Gi00(1), Gi01(0), Gi11(1), mu(1), k1(0), k2(1), a0(1), a1(0)
        {
            for (index_t i = 0; i!=9; ++i)
                D[i] = 0;
        }

        /// Orthotropic quadratic law, material axes rotated by \a phi0
        static SplitAniso Orthotropic(const real_t phi0)
        {
            SplitAniso L;
            L.kind = Ortho;
            const real_t c = math::cos(phi0), s = math::sin(phi0);
            gsMatrix<real_t> Rt(3,3), D0(3,3);
            Rt << c*c,    s*s,   c*s,
                  s*s,    c*c,  -c*s,
                 -2*c*s,  2*c*s, c*c-s*s;
            D0 << 10, 3, 0,
                  3,  2, 0,
                  0,  0, 1.5;
            const gsMatrix<real_t> Dm = Rt.transpose() * D0 * Rt;
            for (index_t i = 0; i!=3; ++i)
                for (index_t j = 0; j!=3; ++j)
                    L.D[3*i+j] = Dm(i,j);
            return L;
        }

        /// Neo-Hooke on the reference metric G = Jref^T Jref plus a fibre along the Cartesian direction \a beta
        static SplitAniso NeoHookeFibre(const gsMatrix<real_t> & Jref, const real_t k1, const real_t k2, const real_t beta)
        {
            SplitAniso L;
            L.kind = NHFibre;
            L.k1 = k1;
            L.k2 = k2;
            const gsMatrix<real_t> G = Jref.transpose() * Jref;
            const gsMatrix<real_t> Gi = G.inverse();
            L.G00 = G(0,0); L.G01 = G(0,1); L.G11 = G(1,1);
            L.detG = G(0,0)*G(1,1) - G(0,1)*G(0,1);
            L.Gi00 = Gi(0,0); L.Gi01 = Gi(0,1); L.Gi11 = Gi(1,1);
            gsMatrix<real_t> dir(2,1);
            dir << math::cos(beta), math::sin(beta);
            const gsMatrix<real_t> a = Jref.inverse() * dir;
            L.a0 = a(0,0);
            L.a1 = a(1,0);
            return L;
        }

        real_t W(const real_t E0, const real_t E1, const real_t E2) const
        {
            const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
            if (kind==Ortho)
            {
                const real_t E[3] = { E0, E1, E2 };
                real_t w = 0;
                for (index_t i = 0; i!=3; ++i)
                    for (index_t j = 0; j!=3; ++j)
                        w += E[i] * D[3*i+j] * E[j];
                return real_t(0.5) * w;
            }
            const real_t g00 = G00 + 2*E0, g11 = G11 + 2*E1, g01 = G01 + E2;
            const real_t dg = g00*g11 - g01*g01;
            if (!(dg > 0 && g00 > 0))
                return nan;
            const real_t r = detG / dg;
            real_t w = real_t(0.5) * mu * (Gi00*g00 + 2*Gi01*g01 + Gi11*g11 + r - 3);
            if (k1 > 0)
            {
                const real_t x = a0*a0*g00 + 2*a0*a1*g01 + a1*a1*g11 - 1;
                w += k1 / (2*k2) * (std::exp(k2*x*x) - 1);
            }
            return math::isfinite(w) ? w : nan;
        }

        real_t W(const gsMatrix<real_t> & E) const { return W(E.data()[0],E.data()[1],E.data()[2]); }

        gsMatrix<real_t> S(const gsMatrix<real_t> & E) const
        {
            gsMatrix<real_t> s(3,1);
            if (kind==Ortho)
            {
                for (index_t i = 0; i!=3; ++i)
                    s(i,0) = D[3*i]*E.data()[0] + D[3*i+1]*E.data()[1] + D[3*i+2]*E.data()[2];
                return s;
            }
            const real_t g00 = G00 + 2*E.data()[0], g11 = G11 + 2*E.data()[1], g01 = G01 + E.data()[2];
            const real_t dg = g00*g11 - g01*g01;
            if (!(dg > 0 && g00 > 0))
            {
                s.setConstant(std::numeric_limits<real_t>::quiet_NaN());
                return s;
            }
            const real_t r = detG / dg;
            s(0,0) = mu * (Gi00 - r * g11 / dg);
            s(1,0) = mu * (Gi11 - r * g00 / dg);
            s(2,0) = mu * (Gi01 + r * g01 / dg);
            if (k1 > 0)
            {
                const real_t x = a0*a0*g00 + 2*a0*a1*g01 + a1*a1*g11 - 1;
                const real_t f = 2 * k1 * x * std::exp(k2*x*x);
                s(0,0) += f * a0*a0;
                s(1,0) += f * a1*a1;
                s(2,0) += f * a0*a1;
            }
            if (!s.allFinite())
                s.setConstant(std::numeric_limits<real_t>::quiet_NaN());
            return s;
        }

        gsMatrix<real_t> C(const gsMatrix<real_t> & E) const
        {
            gsMatrix<real_t> c(3,3);
            if (kind==Ortho)
            {
                for (index_t i = 0; i!=3; ++i)
                    for (index_t j = 0; j!=3; ++j)
                        c(i,j) = D[3*i+j];
                return c;
            }
            const real_t g00 = G00 + 2*E.data()[0], g11 = G11 + 2*E.data()[1], g01 = G01 + E.data()[2];
            const real_t dg = g00*g11 - g01*g01;
            if (!(dg > 0 && g00 > 0))
            {
                c.setConstant(std::numeric_limits<real_t>::quiet_NaN());
                return c;
            }
            const real_t r = detG / dg;
            const real_t gi[2][2] = { { g11/dg, -g01/dg }, { -g01/dg, g00/dg } };
            const index_t pa[3] = { 0, 1, 0 }, pb[3] = { 0, 1, 1 };
            for (index_t i = 0; i!=3; ++i)
                for (index_t j = 0; j!=3; ++j)
                {
                    const index_t a = pa[i], b = pb[i], cc = pa[j], d = pb[j];
                    c(i,j) = mu * r * (2*gi[a][b]*gi[cc][d] + gi[a][cc]*gi[b][d] + gi[a][d]*gi[b][cc]);
                }
            if (k1 > 0)
            {
                const real_t x = a0*a0*g00 + 2*a0*a1*g01 + a1*a1*g11 - 1;
                const real_t f = 4 * k1 * std::exp(k2*x*x) * (1 + 2*k2*x*x);
                const real_t A[3] = { a0*a0, a1*a1, a0*a1 };
                for (index_t i = 0; i!=3; ++i)
                    for (index_t j = 0; j!=3; ++j)
                        c(i,j) += f * A[i] * A[j];
            }
            if (!c.allFinite())
                c.setConstant(std::numeric_limits<real_t>::quiet_NaN());
            return c;
        }

        /// Provider of S and C at an elastic strain
        SplitXProv SCprov() const
        {
            const SplitAniso L = *this;
            return [L](const gsMatrix<real_t> & Eel, gsMatrix<real_t> & Se, gsMatrix<real_t> & Ce)
            {
                Se = L.S(Eel);
                Ce = L.C(Eel);
            };
        }

        /// Provider of W at an elastic strain
        SplitXWProv Wprov() const
        {
            const SplitAniso L = *this;
            return [L](const gsMatrix<real_t> & Eel) -> real_t { return L.W(Eel); };
        }
    };

    /// Orthotropic law (a) of the reference set (material axes rotated by pi/6) or the neo-Hooke-plus-fibre law (b)
    /// (k1 = 10, k2 = 20, fibre at 0.35 rad) on the reference Jacobian \a J
    static SplitAniso SplitAnisoLaw(const bool lawB, const gsMatrix<real_t> & J)
    {
        return lawB ? SplitAniso::NeoHookeFibre(J,10,20,0.35) : SplitAniso::Orthotropic(static_cast<real_t>(EIGEN_PI)/6);
    }

    /// Covariant engineering strain Jref^T eps Jref of the Cartesian strain Rot(alpha) diag(e1,e2) Rot(alpha)^T. Its
    /// principal values with respect to G = Jref^T Jref are e1 and e2.
    static gsMatrix<real_t> SplitAnisoE(const gsMatrix<real_t> & J, const real_t e1, const real_t e2, const real_t alpha)
    {
        gsMatrix<real_t> R(2,2), D(2,2);
        R << math::cos(alpha), -math::sin(alpha), math::sin(alpha), math::cos(alpha);
        D << e1, 0, 0, e2;
        const gsMatrix<real_t> eps = R * D * R.transpose();
        const gsMatrix<real_t> Et = J.transpose() * eps * J;
        gsMatrix<real_t> E(3,1);
        E << Et(0,0), Et(1,1), 2*Et(0,1);
        return E;
    }

    /// Ascending principal stresses of S with respect to G (pencil S x = lambda G^{-1} x)
    static gsMatrix<real_t> SplitAnisoStressPencil(const gsMatrix<real_t> & S, const gsMatrix<real_t> & G)
    {
        const gsMatrix<real_t> Gi = G.inverse();
        return SplitPencilEigenvalues(SplitSym2(S(0,0),S(2,0),S(1,0)),Gi);
    }

    /// Ascending principal strains of E with respect to G (pencil E x = lambda G x)
    static gsMatrix<real_t> SplitAnisoStrainPencil(const gsMatrix<real_t> & E, const gsMatrix<real_t> & G)
    {
        return SplitPencilEigenvalues(SplitSym2(E(0,0),E(2,0)/2,E(1,0)),G);
    }

    /// Criterion state of TensionField 2 from the closed-form stress: 1 if the minor principal stress is not
    /// negative, else -1 if the major principal strain is negative, else 0
    static index_t SplitAnisoClassify(const SplitAniso & law, const gsMatrix<real_t> & E, const gsMatrix<real_t> & G)
    {
        const gsMatrix<real_t> S = law.S(E);
        if (SplitAnisoStressPencil(S,G)(0,0) >= 0)
            return 1;
        return SplitAnisoStrainPencil(E,G)(1,0) < 0 ? -1 : 0;
    }

    /// Candidate strain of a state set. kind: 0 wrinkled candidate, 1 taut, 2 slack, 3 near-taut, 4 near-slack;
    /// cls is the closed-form classification
    struct SplitAnisoState
    {
        real_t e1, e2, alpha;
        gsMatrix<real_t> E;
        index_t kind, cls;
    };

    /// The taut/wrinkled boundary on the path E(s) = (1-s) E0 + s E1, E0 = E[(0.05,0.02,alpha)] and
    /// E1 = E[(0.05,e2,alpha)], by bisection to |ds| <= 1e-14 on the sign of the closed-form minor principal
    /// stress. e2 starts at -0.05 and is doubled until E1 is not taut.
    static real_t SplitAnisoTautBoundary(const SplitAniso & law, const gsMatrix<real_t> & J, const gsMatrix<real_t> & G,
                                         const real_t alpha, gsMatrix<real_t> & E0, gsMatrix<real_t> & E1)
    {
        E0 = SplitAnisoE(J,0.05,0.02,alpha);
        real_t e2 = -0.05;
        E1 = SplitAnisoE(J,0.05,e2,alpha);
        for (index_t k = 0; k!=6 && SplitAnisoClassify(law,E1,G)==1; ++k)
        {
            e2 *= 2;
            E1 = SplitAnisoE(J,0.05,e2,alpha);
        }
        CHECK_EQUAL(1,SplitAnisoClassify(law,E0,G));
        CHECK(SplitAnisoClassify(law,E1,G)!=1);
        real_t lo = 0, hi = 1;
        while (hi - lo > 1e-14)
        {
            const real_t mid = real_t(0.5) * (lo + hi);
            const gsMatrix<real_t> E = (1-mid) * E0 + mid * E1;
            if (SplitAnisoStressPencil(law.S(E),G)(0,0) >= 0)
                lo = mid;
            else
                hi = mid;
        }
        return real_t(0.5) * (lo + hi);
    }

    /// The fixed state set of a law and a metric: 16 wrinkled candidates, 4 taut, 4 slack, the near-taut state
    /// E(s_b + 1e-3) and the near-slack state (5e-4,-0.05,0.9)
    static std::vector<SplitAnisoState> SplitAnisoStates(const SplitAniso & law, const gsMatrix<real_t> & J, const gsMatrix<real_t> & G)
    {
        std::vector<SplitAnisoState> out;
        const real_t wr[4][2] = { { 0.05, -0.03 }, { 0.02, -0.05 }, { 0.08, -0.02 }, { 0.03, -0.06 } };
        const real_t wa[4] = { 0.2, 0.9, 1.4, 2.5 };
        const real_t ta[2][2] = { { 0.05, 0.02 }, { 0.04, 0.03 } };
        const real_t taa[2] = { 0.2, 1.4 };
        const real_t sl[2][2] = { { -0.02, -0.04 }, { -0.05, -0.01 } };
        const real_t sla[2] = { 0.9, 2.5 };
        auto add = [&](const real_t e1, const real_t e2, const real_t al, const index_t kind)
        {
            SplitAnisoState s;
            s.e1 = e1; s.e2 = e2; s.alpha = al; s.kind = kind;
            s.E = SplitAnisoE(J,e1,e2,al);
            s.cls = SplitAnisoClassify(law,s.E,G);
            out.push_back(s);
        };
        for (index_t i = 0; i!=4; ++i)
            for (index_t j = 0; j!=4; ++j)
                add(wr[i][0],wr[i][1],wa[j],0);
        for (index_t i = 0; i!=2; ++i)
            for (index_t j = 0; j!=2; ++j)
                add(ta[i][0],ta[i][1],taa[j],1);
        for (index_t i = 0; i!=2; ++i)
            for (index_t j = 0; j!=2; ++j)
                add(sl[i][0],sl[i][1],sla[j],2);

        gsMatrix<real_t> E0, E1;
        const real_t sb = SplitAnisoTautBoundary(law,J,G,0.9,E0,E1);
        SplitAnisoState nt;
        nt.e1 = nt.e2 = nt.alpha = std::numeric_limits<real_t>::quiet_NaN();
        nt.kind = 3;
        nt.E = (1 - (sb + 1e-3)) * E0 + (sb + 1e-3) * E1;
        nt.cls = SplitAnisoClassify(law,nt.E,G);
        out.push_back(nt);
        add(5e-4,-0.05,0.9,4);
        return out;
    }

    /// The state set with the counts asserted: at least 6 wrinkled candidates, near-boundary states wrinkled
    static std::vector<SplitAnisoState> SplitAnisoStatesChecked(const SplitAniso & law, const gsMatrix<real_t> & J, const gsMatrix<real_t> & G)
    {
        const std::vector<SplitAnisoState> st = SplitAnisoStates(law,J,G);
        index_t nw = 0;
        for (size_t i = 0; i!=st.size(); ++i)
        {
            if (st[i].kind==0 && st[i].cls==0)
                ++nw;
            if (st[i].kind>=3)
                CHECK_EQUAL(0,st[i].cls);
        }
        gsInfo << "  wrinkled candidates: " << nw << " of 16\n";
        CHECK(nw >= 6);
        return st;
    }

    static bool SplitAnisoWrinkledCand(const SplitAnisoState & s) { return s.kind==0 && s.cls==0; }

    /// Pencil Smin(S(E)) <= -1e-2 Sref and Emax(E) >= 1e-2 Eref: far from both interfaces
    static bool SplitAnisoAway(const SplitAniso & law, const gsMatrix<real_t> & E, const gsMatrix<real_t> & G)
    {
        const gsMatrix<real_t> lS = SplitAnisoStressPencil(law.S(E),G);
        const gsMatrix<real_t> lE = SplitAnisoStrainPencil(E,G);
        const real_t Sref = math::max(math::abs(lS(0,0)),math::abs(lS(1,0)));
        const real_t Eref = math::max(math::abs(lE(0,0)),math::abs(lE(1,0)));
        return lS(0,0) <= -1e-2 * Sref && lE(1,0) >= 1e-2 * Eref;
    }

    /// Exact-mode options with the local solver \a mode, regularisation \a eta and no slack multiplier
    static SplitXOpt SplitAnisoOpt(const index_t mode, const real_t eta = 0)
    {
        SplitXOpt opt = SplitExactOptions();
        opt.exactLocalSolve = mode;
        opt.slackMultiplier = 0;
        opt.eta = eta;
        return opt;
    }

    /// The pointwise core with S = S(E) of the law, the tangent argument \a Cin and the S/C provider of the law;
    /// the W provider only if \a withW
    static SplitXRes SplitAnisoCall(const SplitFixture<2,false> & fx, const SplitAniso & law, const gsMatrix<real_t> & E,
                                    const gsMatrix<real_t> & G, const SplitXOpt & opt, const bool tangent, const bool energy,
                                    const gsMatrix<real_t> & Cin, const bool withW)
    {
        const gsMatrix<real_t> S = law.S(E);
        const SplitXWProv Wp = withW ? law.Wprov() : SplitXWProv();
        return fx.split->wrinkle(S,Cin,E,G,opt,tangent,SplitXdCProv(),law.SCprov(),Wp,energy);
    }

    /// The pointwise core with the consistent tangent C(E) and the W provider
    static SplitXRes SplitAnisoRun(const SplitFixture<2,false> & fx, const SplitAniso & law, const gsMatrix<real_t> & E,
                                   const gsMatrix<real_t> & G, const index_t mode, const real_t eta = 0,
                                   const bool tangent = true, const bool energy = true)
    {
        return SplitAnisoCall(fx,law,E,G,SplitAnisoOpt(mode,eta),tangent,energy,law.C(E),true);
    }

    /// Round-off floor of the residual at an iterate: 64 eps ||C'||_F max(1,||E'||)
    static real_t SplitAnisoFloorR(const gsMatrix<real_t> & Cp, const gsMatrix<real_t> & Ep)
    {
        return 64 * std::numeric_limits<real_t>::epsilon() * Cp.norm() * math::max(real_t(1),Ep.norm());
    }

    /// Independent residual check at an accepted wrinkled result (Eta 0): the closed-form law at
    /// E' = E + gamma n1(theta) has n1.S' and n2.S' within the stall tolerance, and S_eff is S'
    static bool SplitAnisoResidualOK(const SplitAniso & law, const gsMatrix<real_t> & E, const SplitXRes & r)
    {
        if (!(r.state==0 && r.gamma > 0 && math::isfinite(r.theta)))
        {
            gsInfo << "  residual check: not a wrinkled result (state " << r.state << ")\n";
            return false;
        }
        const gsVector<real_t,3> n1 = SplitN1(r.theta), n2 = SplitN2(r.theta);
        const gsMatrix<real_t> Ep = E + r.gamma * n1;
        const gsMatrix<real_t> Sp = law.S(Ep), Cp = law.C(Ep);
        if (!(Sp.allFinite() && Cp.allFinite()))
            return false;
        const real_t R1 = n1.dot(Sp.col(0)), R2 = n2.dot(Sp.col(0));
        const real_t tolS = SplitExactTolAt(Sp,Cp,Ep,real_t(1e-9));
        const real_t floorR = SplitAnisoFloorR(Cp,Ep);
        const real_t dS = (r.S - Sp).cwiseAbs().maxCoeff();
        const bool ok = math::max(math::abs(R1),math::abs(R2)) <= tolS && dS <= 1e-12 * Sp.norm() + floorR;
        if (!ok)
            gsInfo << "  residual check failed: |R1| = " << math::abs(R1) << ", |R2| = " << math::abs(R2) << " (tolS " << tolS
                   << "), |S_eff - S'|_inf = " << dS << " (tol " << 1e-12 * Sp.norm() + floorR << ")\n";
        return ok;
    }

    /// KKT certificate of the global minimum of W(E + P) over the PSD cone at Eta 0: P* is 0 (taut), gamma n1(theta)
    /// (wrinkled) or -E (slack) according to the state, which must be \a expected. Checks that S' = S(E + P*) is
    /// positive semidefinite, P* is positive semidefinite and complementarity S'.P* = 0.
    static bool SplitAnisoKKT(const SplitAniso & law, const gsMatrix<real_t> & E, const gsMatrix<real_t> & G,
                              const SplitXRes & r, const index_t expected)
    {
        bool ok = (r.state == expected);
        gsMatrix<real_t> P(3,1);
        P.setZero();
        if (r.state==0)
        {
            if (r.gamma > 0)
            {
                const gsVector<real_t,3> n1 = SplitN1(r.theta);
                P = r.gamma * n1;
            }
            else
                ok = false;
        }
        else if (r.state==-1)
        {
            P = -E;
            ok = ok && SplitAnisoStrainPencil(E,G)(1,0) < 0;
        }
        const gsMatrix<real_t> Ep = E + P;
        const gsMatrix<real_t> Sp = law.S(Ep), Cp = law.C(Ep);
        if (!(Sp.allFinite() && Cp.allFinite()))
            return false;
        const real_t floorR = SplitAnisoFloorR(Cp,Ep);
        const gsMatrix<real_t> lam = SplitAnisoStressPencil(Sp,G);
        const real_t lmax = math::max(math::abs(lam(0,0)),math::abs(lam(1,0)));
        const bool psd = lam(0,0) >= -(1e-10 * lmax + floorR * G.norm());
        const real_t comp = math::abs(Sp.col(0).dot(P.col(0)));
        const bool compl0 = comp <= 1e-10 * Sp.norm() * P.norm() + floorR * P.norm();
        ok = ok && psd && compl0;
        if (!ok)
            gsInfo << "  KKT failed: state " << r.state << " (expected " << expected << "), lmin(S') = " << lam(0,0) << " (lmax " << lmax
                   << "), |S'.P*| = " << comp << " (||S'|| ||P*|| " << Sp.norm() * P.norm() << ")\n";
        return ok;
    }

    /// kappa = ||C'||_2 ||X||_2 with X = diag(1,2 gamma) Jac0^{-1}, and the tolerance tolC of the exact solve, from the
    /// closed-form law at the elastic strain of the wrinkled result \a r
    static void SplitAnisoKappa(const SplitAniso & law, const SplitXRes & r, real_t & kappa, real_t & tolC)
    {
        const real_t c = math::cos(r.theta), s = math::sin(r.theta);
        const gsVector<real_t,3> n1 = SplitN1(r.theta), n2 = SplitN2(r.theta);
        gsVector<real_t,3> n4;
        n4 << s*s, c*c, -2*c*s;
        const gsMatrix<real_t> Sp = law.S(r.E), Cp = law.C(r.E);
        const gsVector<real_t,3> Cn1 = Cp * n1, Cn2 = Cp * n2;
        const real_t sw = n4.dot(Sp.col(0));
        gsMatrix<real_t> Jac(2,2), Dg(2,2);
        Jac << n1.dot(Cn1), 2*r.gamma*n1.dot(Cn2), n2.dot(Cn1), sw + 2*r.gamma*n2.dot(Cn2);
        Dg << 1, 0, 0, 2*r.gamma;
        const gsMatrix<real_t> X = Dg * Jac.inverse();
        kappa = Cp.operatorNorm() * X.operatorNorm();
        tolC = SplitExactTolAt(Sp,Cp,r.E,real_t(1e-12));
    }

    static bool SplitAnisoSameScalar(const real_t a, const real_t b)
    {
        return a == b || (std::isnan(a) && std::isnan(b));
    }

    static bool SplitAnisoSameMatrix(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
    {
        if (A.rows()!=B.rows() || A.cols()!=B.cols())
            return false;
        for (index_t i = 0; i!=A.size(); ++i)
            if (!SplitAnisoSameScalar(A.data()[i],B.data()[i]))
                return false;
        return true;
    }

    /// Names of the fields in which two results differ bitwise (NaN equals NaN); empty if none. The energy is compared
    /// if \a withW and the outcome flags if \a withFlags.
    static std::string SplitAnisoDiff(const SplitXRes & a, const SplitXRes & b, const bool withW, const bool withFlags)
    {
        std::string d;
        if (a.state != b.state) d += " state";
        if (!SplitAnisoSameMatrix(a.S,b.S)) d += " S";
        if (!SplitAnisoSameMatrix(a.C,b.C)) d += " C";
        if (!SplitAnisoSameMatrix(a.E,b.E)) d += " E";
        if (!SplitAnisoSameScalar(a.theta,b.theta)) d += " theta";
        if (!SplitAnisoSameScalar(a.gamma,b.gamma)) d += " gamma";
        if (!SplitAnisoSameScalar(a.gammaHat,b.gammaHat)) d += " gammaHat";
        if (!SplitAnisoSameScalar(a.residual,b.residual)) d += " residual";
        if (withW && !SplitAnisoSameScalar(a.W,b.W)) d += " W";
        if (a.exactIterations != b.exactIterations) d += " exactIterations";
        if (withFlags)
        {
            if (a.guardFallback != b.guardFallback) d += " guardFallback";
            if (a.tftFailure != b.tftFailure) d += " tftFailure";
            if (a.exactFallback != b.exactFallback) d += " exactFallback";
            if (a.exactReclassified != b.exactReclassified) d += " exactReclassified";
            if (a.exactSaddle != b.exactSaddle) d += " exactSaddle";
            if (a.exactWrongBranch != b.exactWrongBranch) d += " exactWrongBranch";
            if (a.exactStagnated != b.exactStagnated) d += " exactStagnated";
            if (a.exactSkipped != b.exactSkipped) d += " exactSkipped";
            if (a.exactUsed2x2 != b.exactUsed2x2) d += " exactUsed2x2";
            if (a.exactMinimised != b.exactMinimised) d += " exactMinimised";
        }
        return d;
    }

    /// Minimum of W(E + P) over a grid and compass-search polish of the PSD matrices
    /// P = Rot(phi) diag(p1,p2) Rot(phi)^T (components, engineering Voigt (P_11,P_22,2 P_12) added to E). Grid:
    /// phi = k pi/36 (k = 0..35), p_i in {0} and pmax (j/12)^2 (j = 1..12) with
    /// pmax = 2 max(lambda_max(-E_t), \a gammaW, ||E_t||_2), E_t the component matrix of E. The polish starts
    /// from the 3 best grid points with p_i = s_i^2; steps pi/72 and sqrt(pmax)/24, halved on failure down to a
    /// phi step of 1e-10, at most 2000 evaluations per start. A NaN W counts as +inf.
    static real_t SplitAnisoBruteMin(const SplitAniso & law, const gsMatrix<real_t> & E, const real_t gammaW)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t inf = std::numeric_limits<real_t>::infinity();
        const real_t E0 = E(0,0), E1 = E(1,0), E2 = E(2,0);
        const real_t mid = real_t(0.5) * (E0 + E1);
        const real_t rad = math::sqrt(real_t(0.25) * (E0-E1) * (E0-E1) + real_t(0.25) * E2 * E2);
        const real_t lmin = mid - rad, lmax = mid + rad;
        const real_t pmax = 2 * math::max(math::max(-lmin,gammaW),math::max(math::abs(lmin),math::abs(lmax)));

        auto Wat = [&](const real_t phi, const real_t p1, const real_t p2) -> real_t
        {
            const real_t c = math::cos(phi), s = math::sin(phi);
            const real_t P00 = p1*c*c + p2*s*s, P11 = p1*s*s + p2*c*c, P01 = (p1-p2)*c*s;
            const real_t w = law.W(E0 + P00,E1 + P11,E2 + 2*P01);
            return math::isfinite(w) ? w : inf;
        };

        struct Pt { real_t w, phi, s1, s2; };
        Pt best[3];
        for (index_t i = 0; i!=3; ++i)
        {
            best[i].w = inf; best[i].phi = best[i].s1 = best[i].s2 = 0;
        }
        for (index_t k = 0; k!=36; ++k)
            for (index_t j1 = 0; j1!=13; ++j1)
                for (index_t j2 = 0; j2!=13; ++j2)
                {
                    const real_t q1 = real_t(j1)/12, q2 = real_t(j2)/12;
                    const real_t p1 = pmax * q1 * q1, p2 = pmax * q2 * q2;
                    Pt t;
                    t.phi = k * pi / 36;
                    t.s1 = math::sqrt(p1);
                    t.s2 = math::sqrt(p2);
                    t.w = Wat(t.phi,p1,p2);
                    for (index_t i = 0; i!=3; ++i)
                        if (t.w < best[i].w)
                        {
                            for (index_t m = 2; m > i; --m)
                                best[m] = best[m-1];
                            best[i] = t;
                            break;
                        }
                }

        real_t wmin = best[0].w;
        for (index_t m = 0; m!=3; ++m)
        {
            Pt x = best[m];
            if (!math::isfinite(x.w))
                continue;
            real_t dphi = pi / 72, ds = real_t(0.5) * math::sqrt(pmax) / 12;
            index_t nev = 0;
            while (dphi > 1e-10 && nev < 2000)
            {
                bool improved = false;
                for (index_t mv = 0; mv!=6 && nev < 2000; ++mv)
                {
                    Pt t = x;
                    const real_t sg = (mv%2==0) ? 1 : -1;
                    if (mv/2==0) t.phi += sg * dphi;
                    else if (mv/2==1) t.s1 += sg * ds;
                    else t.s2 += sg * ds;
                    t.w = Wat(t.phi,t.s1*t.s1,t.s2*t.s2);
                    ++nev;
                    if (t.w < x.w)
                    {
                        x = t;
                        improved = true;
                        break;
                    }
                }
                if (!improved)
                {
                    dphi *= real_t(0.5);
                    ds *= real_t(0.5);
                }
            }
            wmin = math::min(wmin,x.w);
        }
        return wmin;
    }

    /// Hessian B^T C' B + 2 Shat' of Phi(v) = W(E + v (x) v) for the stress \a Sp and tangent \a Cp at E + v (x) v
    static gsMatrix<real_t> SplitAnisoHes(const gsMatrix<real_t> & Cp, const gsMatrix<real_t> & Sp, const real_t v0, const real_t v1)
    {
        gsMatrix<real_t> B(3,2), Sh(2,2);
        B << 2*v0, 0, 0, 2*v1, 2*v1, 2*v0;
        Sh << Sp(0,0), Sp(2,0), Sp(2,0), Sp(1,0);
        const gsMatrix<real_t> H = B.transpose() * Cp * B + 2 * Sh;
        const gsMatrix<real_t> Hs = real_t(0.5) * (H + H.transpose());
        return Hs;
    }

    /// M11 = n1.C'n1 and Jac0 = [[M11, 2 gamma M12], [M21, sw + 2 gamma M22]] of the acceptance test at theta, gamma
    static gsMatrix<real_t> SplitAnisoJac0(const gsMatrix<real_t> & Cp, const gsMatrix<real_t> & Sp, const real_t theta,
                                           const real_t gamma, real_t & M11)
    {
        const real_t c = math::cos(theta), s = math::sin(theta);
        const gsVector<real_t,3> n1 = SplitN1(theta), n2 = SplitN2(theta);
        gsVector<real_t,3> n4;
        n4 << s*s, c*c, -2*c*s;
        const gsVector<real_t,3> Cn1 = Cp * n1, Cn2 = Cp * n2;
        M11 = n1.dot(Cn1);
        gsMatrix<real_t> Jac(2,2);
        Jac << M11, 2*gamma*n1.dot(Cn2), n2.dot(Cn1), n4.dot(Sp.col(0)) + 2*gamma*n2.dot(Cn2);
        return Jac;
    }

    static real_t SplitAnisoLambdaMin(const gsMatrix<real_t> & H)
    {
        const real_t mid = real_t(0.5) * (H(0,0) + H(1,1));
        const real_t rad = math::sqrt(real_t(0.25) * (H(0,0)-H(1,1)) * (H(0,0)-H(1,1)) + H(0,1)*H(0,1));
        return mid - rad;
    }

    /// Isotropic (pure neo-Hooke) and coaxial states, where the closed-form direction passes the isotropy guard and
    /// ExactLocalSolve 1 runs the scalar pass
    struct SplitAnisoCase
    {
        SplitAniso law;
        gsMatrix<real_t> G, E;
        real_t e1, e2, alpha;
        index_t id;
    };

    /// id 0: pure neo-Hooke on both metrics; id 1: law (b) with the fibre along e1 (beta = alpha), orthonormal G;
    /// id 2: law (a) with phi0 = 0 at alpha = 0, orthonormal G. Only states that classify as wrinkled.
    static std::vector<SplitAnisoCase> SplitAnisoIsoCases()
    {
        std::vector<SplitAnisoCase> out;
        const real_t wr[4][2] = { { 0.05, -0.03 }, { 0.02, -0.05 }, { 0.08, -0.02 }, { 0.03, -0.06 } };
        const real_t wa[4] = { 0.2, 0.9, 1.4, 2.5 };
        auto push = [&](const SplitAniso & law, const gsMatrix<real_t> & J, const real_t e1, const real_t e2, const real_t al, const index_t id)
        {
            SplitAnisoCase c;
            c.law = law;
            c.G = J.transpose() * J;
            c.E = SplitAnisoE(J,e1,e2,al);
            c.e1 = e1; c.e2 = e2; c.alpha = al; c.id = id;
            if (SplitAnisoClassify(law,c.E,c.G)==0)
                out.push_back(c);
        };
        for (index_t ig = 0; ig!=2; ++ig)
        {
            const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
            const SplitAniso nh = SplitAniso::NeoHookeFibre(J,0,20,0.35);
            for (index_t i = 0; i!=4; ++i)
                for (index_t j = 0; j!=4; ++j)
                    push(nh,J,wr[i][0],wr[i][1],wa[j],0);
        }
        const gsMatrix<real_t> I2 = SplitIdentity2();
        for (index_t i = 0; i!=4; ++i)
            for (index_t j = 0; j!=4; ++j)
                push(SplitAniso::NeoHookeFibre(I2,10,20,wa[j]),I2,wr[i][0],wr[i][1],wa[j],1);
        for (index_t i = 0; i!=4; ++i)
            push(SplitAniso::Orthotropic(0),I2,wr[i][0],wr[i][1],0,2);
        return out;
    }

    /// W(0) = 0, S = grad W and C = dS/dE by central differences (h = 1e-6) for both laws, both metrics, at a wrinkled,
    /// a taut and a slack candidate: the oracle laws are mutually consistent
    TEST(SynthLaws_FD_selfcheck)
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        const real_t h = 1e-6;
        const real_t st[3][3] = { { 0.02, -0.05, 0.9 }, { 0.05, 0.02, 0.2 }, { -0.02, -0.04, 0.9 } };
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const gsMatrix<real_t> Z = gsMatrix<real_t>::Zero(3,1);
                const real_t floor0 = 1e3 * eps * law.C(Z).cwiseAbs().maxCoeff();
                CHECK(math::abs(law.W(Z)) <= floor0);
                CHECK(law.S(Z).cwiseAbs().maxCoeff() <= floor0);
                for (index_t k = 0; k!=3; ++k)
                {
                    const gsMatrix<real_t> E = SplitAnisoE(J,st[k][0],st[k][1],st[k][2]);
                    const gsMatrix<real_t> S = law.S(E), C = law.C(E);
                    gsMatrix<real_t> Sfd(3,1), Cfd(3,3);
                    for (index_t j = 0; j!=3; ++j)
                    {
                        gsMatrix<real_t> Ep = E, Em = E;
                        Ep(j,0) += h;
                        Em(j,0) -= h;
                        Sfd(j,0) = (law.W(Ep) - law.W(Em)) / (2*h);
                        Cfd.col(j) = (law.S(Ep) - law.S(Em)) / (2*h);
                    }
                    const real_t floorC = 1e3 * eps * C.cwiseAbs().maxCoeff();
                    const real_t errS = (Sfd - S).cwiseAbs().maxCoeff(), errC = (Cfd - C).cwiseAbs().maxCoeff();
                    gsInfo << "  law " << (lb==1 ? "b" : "a") << ", G " << ig << ", state " << k << ": |S - FD W| = " << errS
                           << ", |C - FD S| = " << errC << "\n";
                    CHECK(S.allFinite() && C.allFinite());
                    CHECK(errS <= 1e-6 * S.cwiseAbs().maxCoeff() + floorC);
                    CHECK(errC <= 1e-6 * C.cwiseAbs().maxCoeff() + floorC);
                }
            }
    }

    /// Global-minimum certificate at Eta 0 for the exact mode: P* is 0, gamma n1 or -E; S(E + P*) and P* are positive
    /// semidefinite and complementary. A KKT point of the convex problem is its global minimum.
    static void SplitAnisoKKTRun(const bool lawB)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        index_t nWr = 0, nAll = 0;
        for (index_t ig = 0; ig!=2; ++ig)
        {
            const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
            const gsMatrix<real_t> G = J.transpose() * J;
            const SplitAniso law = SplitAnisoLaw(lawB,J);
            const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
            for (size_t i = 0; i!=st.size(); ++i)
            {
                const SplitXRes r = SplitAnisoRun(fx,law,st[i].E,G,1);
                ++nAll;
                CHECK_EQUAL(st[i].cls,r.state);
                CHECK(SplitAnisoKKT(law,st[i].E,G,r,st[i].cls));
                if (r.state==0)
                {
                    ++nWr;
                    CHECK(!r.exactFallback);
                    CHECK(SplitAnisoResidualOK(law,st[i].E,r));
                }
            }
            if (lawB)
            {
                // A failed direction (indefinite tangent argument): the minimisation starts away from the solution
                const real_t alphas[3] = { 0.2, 1.4, 2.5 };
                for (index_t a = 0; a!=3; ++a)
                {
                    const gsMatrix<real_t> E = SplitAnisoE(J,0.05,-0.03,alphas[a]);
                    const gsMatrix<real_t> negC = -law.C(E);
                    const SplitXRes r = SplitAnisoCall(fx,law,E,G,SplitAnisoOpt(1),true,true,negC,true);
                    ++nAll;
                    CHECK(SplitAnisoKKT(law,E,G,r,0));
                    CHECK(SplitAnisoResidualOK(law,E,r));
                }
            }
        }
        gsInfo << "  KKT certificates: " << nAll << " states, " << nWr << " wrinkled\n";
    }

    TEST(ExactMin_KKT_certificate_lawA)
    {
        SplitAnisoKKTRun(false);
    }

    TEST(ExactMin_KKT_certificate_lawB)
    {
        SplitAnisoKKTRun(true);
    }

    /// The accepted energy is not beaten by any PSD P found by a grid search with local polish
    static void SplitAnisoBruteRun(const bool lawB)
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        real_t minMargin = std::numeric_limits<real_t>::infinity();
        for (index_t ig = 0; ig!=2; ++ig)
        {
            const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
            const gsMatrix<real_t> G = J.transpose() * J;
            const SplitAniso law = SplitAnisoLaw(lawB,J);
            const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
            for (size_t i = 0; i!=st.size(); ++i)
            {
                const SplitXRes r = SplitAnisoRun(fx,law,st[i].E,G,1);
                CHECK_EQUAL(st[i].cls,r.state);
                CHECK(math::isfinite(r.W));
                const real_t Wmin = SplitAnisoBruteMin(law,st[i].E,r.state==0 ? r.gamma : real_t(0));
                const real_t tol = 1e-8 * math::abs(r.W) + 1e3 * eps * law.C(st[i].E).norm();
                minMargin = math::min(minMargin,Wmin - r.W);
                CHECK(Wmin >= r.W - tol);
            }
        }
        gsInfo << "  smallest (search minimum - W_eff) = " << minMargin << " (negative only within the tolerance)\n";
    }

    TEST(ExactMin_bruteForce_lawA)
    {
        SplitAnisoBruteRun(false);
    }

    TEST(ExactMin_bruteForce_lawB)
    {
        SplitAnisoBruteRun(true);
    }

    /// S_eff = grad W_eff: central differences of the energy of the minimisation, every evaluation a new direct call
    TEST(ExactMin_energy_gradient)
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        const real_t h = 1e-6;
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const real_t etas[2] = { 0, 1e-3 };
        index_t nChecked = 0;
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
                for (size_t i = 0; i!=st.size(); ++i)
                {
                    if (st[i].kind > 1)
                        continue;
                    for (index_t ie = 0; ie!=2; ++ie)
                    {
                        const SplitXRes r0 = SplitAnisoRun(fx,law,st[i].E,G,1,etas[ie]);
                        CHECK_EQUAL(st[i].cls,r0.state);
                        CHECK(!r0.exactFallback);
                        gsMatrix<real_t> g(3,1);
                        bool stable = true;
                        for (index_t j = 0; j!=3; ++j)
                        {
                            gsMatrix<real_t> Ep = st[i].E, Em = st[i].E;
                            Ep(j,0) += h;
                            Em(j,0) -= h;
                            const SplitXRes rp = SplitAnisoRun(fx,law,Ep,G,1,etas[ie],false,true);
                            const SplitXRes rm = SplitAnisoRun(fx,law,Em,G,1,etas[ie],false,true);
                            stable = stable && rp.state==r0.state && rm.state==r0.state && !rp.exactFallback && !rm.exactFallback;
                            g(j,0) = (rp.W - rm.W) / (2*h);
                        }
                        CHECK(stable);
                        const real_t err = (g - r0.S).cwiseAbs().maxCoeff();
                        const real_t tol = 1e-6 * r0.S.norm() + 1e3 * eps * law.C(st[i].E).norm();
                        if (!(err <= tol))
                            gsInfo << "  law " << (lb==1 ? "b" : "a") << ", G " << ig << ", state " << i << ", eta " << etas[ie]
                                   << ": |S_eff - FD W_eff| = " << err << " (tol " << tol << ")\n";
                        CHECK(err <= tol);
                        ++nChecked;
                    }
                }
            }

        // Synthetic affine and switched laws with the wrapped W of the SplitSynthetic input
        SplitSynthetic sy(true);
        for (index_t k = 0; k!=2; ++k)
            for (index_t ie = 0; ie!=2; ++ie)
            {
                const SplitSynthetic::Kind kind = k==0 ? SplitSynthetic::Affine : SplitSynthetic::Switch;
                const real_t par = k==0 ? real_t(0) : real_t(1e-2);
                const SplitSynthetic::Provider prov = sy.provider(kind,par);
                const gsMatrix<real_t> S0 = sy.S, C0 = sy.C, E0 = sy.E, w0 = sy.w;
                const real_t lin = par * sy.sNorm;
                const SplitXWProv Wp = [=](const gsMatrix<real_t> & Ee) -> real_t
                {
                    const gsMatrix<real_t> d = Ee - E0;
                    const gsMatrix<real_t> Cd = C0 * d;
                    return S0.col(0).dot(d.col(0)) + real_t(0.5) * d.col(0).dot(Cd.col(0)) + lin * w0.col(0).dot(d.col(0));
                };
                SplitXOpt opt = sy.optEx;
                opt.eta = etas[ie];
                opt.exactLocalSolve = 1;
                auto call = [&](const gsMatrix<real_t> & Ee) -> SplitXRes
                {
                    gsMatrix<real_t> S, C;
                    prov(Ee,S,C);
                    return sy.fx.split->wrinkle(S,C,Ee,sy.I2,opt,true,SplitXdCProv(),prov,Wp,true);
                };
                const SplitXRes r0 = call(E0);
                CHECK_EQUAL(0,r0.state);
                CHECK(!r0.exactFallback);
                gsMatrix<real_t> g(3,1);
                bool stable = true;
                for (index_t j = 0; j!=3; ++j)
                {
                    gsMatrix<real_t> Ep = E0, Em = E0;
                    Ep(j,0) += h;
                    Em(j,0) -= h;
                    const SplitXRes rp = call(Ep), rm = call(Em);
                    stable = stable && rp.state==r0.state && rm.state==r0.state && !rp.exactFallback && !rm.exactFallback;
                    g(j,0) = (rp.W - rm.W) / (2*h);
                }
                CHECK(stable);
                const real_t err = (g - r0.S).cwiseAbs().maxCoeff();
                const real_t tol = 1e-6 * r0.S.norm() + 1e3 * eps * C0.norm();
                gsInfo << "  synthetic " << (k==0 ? "affine" : "switch") << ", eta " << etas[ie] << ": |S_eff - FD W_eff| = " << err
                       << " (tol " << tol << ")\n";
                CHECK(err <= tol);
                ++nChecked;
            }
        gsInfo << "  gradient checks: " << nChecked << "\n";
    }

    /// J_ex = dS_eff/dE by central differences of S_eff, at wrinkled states away from both interfaces; symmetric, with
    /// the wrinkle direction in its kernel
    TEST(ExactMin_tangent_FD)
    {
        const real_t h = 1e-6;
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
                index_t nUsed = 0;
                for (size_t i = 0; i!=st.size(); ++i)
                {
                    if (!SplitAnisoWrinkledCand(st[i]) || !SplitAnisoAway(law,st[i].E,G))
                        continue;
                    ++nUsed;
                    const SplitXRes r = SplitAnisoRun(fx,law,st[i].E,G,1);
                    CHECK_EQUAL(0,r.state);
                    gsMatrix<real_t> Jfd(3,3);
                    bool stable = true;
                    for (index_t j = 0; j!=3; ++j)
                    {
                        gsMatrix<real_t> Ep = st[i].E, Em = st[i].E;
                        Ep(j,0) += h;
                        Em(j,0) -= h;
                        const SplitXRes rp = SplitAnisoRun(fx,law,Ep,G,1,0,false,false);
                        const SplitXRes rm = SplitAnisoRun(fx,law,Em,G,1,0,false,false);
                        stable = stable && rp.state==0 && rm.state==0 && !rp.exactFallback && !rm.exactFallback
                                 && !rp.exactStagnated && !rm.exactStagnated;
                        Jfd.col(j) = (rp.S - rm.S) / (2*h);
                    }
                    CHECK(stable);
                    const real_t Jmax = r.C.cwiseAbs().maxCoeff();
                    const real_t err = (r.C - Jfd).cwiseAbs().maxCoeff();
                    if (!(err <= 1e-6 * Jmax))
                        gsInfo << "  law " << (lb==1 ? "b" : "a") << ", G " << ig << ", state " << i << ": |J_ex - FD| = " << err
                               << " (tol " << 1e-6 * Jmax << ")\n";
                    CHECK(err <= 1e-6 * Jmax);
                    CHECK(SplitAsymmetry(r.C) <= 1e-12);
                    const gsVector<real_t,3> n1 = SplitN1(r.theta);
                    const gsMatrix<real_t> Jn1 = r.C * n1;
                    CHECK(Jn1.cwiseAbs().maxCoeff() <= 1e-10 * Jmax * n1.cwiseAbs().maxCoeff());
                }
                gsInfo << "  law " << (lb==1 ? "b" : "a") << ", G " << ig << ": " << nUsed << " states away from the interfaces\n";
                CHECK(nUsed >= 4);
            }
    }

    /// W_eff and S_eff are continuous across the taut and the slack interface (Eta 0, no slack multiplier)
    TEST(ExactMin_continuity_interfaces)
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        const real_t delta = 1e-6, dSlack = 1e-7;
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);

                gsMatrix<real_t> E0, E1;
                const real_t sb = SplitAnisoTautBoundary(law,J,G,0.9,E0,E1);
                const gsMatrix<real_t> dE = E1 - E0;
                const gsMatrix<real_t> Em = E0 + (sb - delta) * dE, Ep = E0 + (sb + delta) * dE, Eb = E0 + sb * dE;
                const SplitXRes rm = SplitAnisoRun(fx,law,Em,G,1);
                const SplitXRes rp = SplitAnisoRun(fx,law,Ep,G,1);
                CHECK_EQUAL(1,rm.state);
                CHECK_EQUAL(0,rp.state);
                CHECK(!rp.exactFallback);
                const real_t dS = (rp.S - rm.S).cwiseAbs().maxCoeff();
                const real_t dW = math::abs(rp.W - rm.W);
                const real_t boundS = 10 * law.C(Eb).norm() * dE.norm() * 2 * delta;
                const real_t boundW = 10 * law.S(Eb).norm() * dE.norm() * 2 * delta;
                gsInfo << "  law " << (lb==1 ? "b" : "a") << ", G " << ig << ", taut interface: |dS| = " << dS << " (bound " << boundS
                       << "), |dW| = " << dW << " (bound " << boundW << ")\n";
                CHECK(dS <= boundS);
                CHECK(dW <= boundW);

                const gsMatrix<real_t> Eneg = SplitAnisoE(J,-dSlack,-0.05,0.9), Epos = SplitAnisoE(J,dSlack,-0.05,0.9);
                const gsMatrix<real_t> Z = gsMatrix<real_t>::Zero(3,1);
                const real_t C0n = law.C(Z).norm();
                const SplitXRes sm = SplitAnisoRun(fx,law,Eneg,G,1);
                const SplitXRes sp = SplitAnisoRun(fx,law,Epos,G,1);
                CHECK_EQUAL(-1,sm.state);
                CHECK(sm.S.cwiseAbs().maxCoeff() == 0);
                CHECK(sm.W == 0);
                CHECK(sp.state==0 || sp.state==-1);
                CHECK(!sp.exactFallback);
                const real_t Sp = sp.S.cwiseAbs().maxCoeff();
                const real_t Wp = math::abs(sp.W);
                gsInfo << "  law " << (lb==1 ? "b" : "a") << ", G " << ig << ", slack interface: |S+| = " << Sp << " (bound "
                       << 10 * C0n * dSlack << "), |W+| = " << Wp << " (bound " << 10 * C0n * dSlack * dSlack + 1e3 * eps * C0n << ")\n";
                CHECK(Sp <= 10 * C0n * dSlack);
                CHECK(Wp <= 10 * C0n * dSlack * dSlack + 1e3 * eps * C0n);
            }
    }

    /// Anisotropic wrinkled states, where the isotropy guard rejects the closed-form direction: the minimisation runs
    /// and is accepted without fallback
    TEST(ExactMin_aniso_wrinkled_noFallback_minimised)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        index_t nIt3 = 0;
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
                for (size_t i = 0; i!=st.size(); ++i)
                {
                    if (!SplitAnisoWrinkledCand(st[i]) && st[i].kind!=3)
                        continue;
                    const SplitXRes r = SplitAnisoRun(fx,law,st[i].E,G,1);
                    CHECK_EQUAL(0,r.state);
                    CHECK(r.guardFallback);
                    CHECK(r.exactMinimised);
                    CHECK(!r.exactFallback);
                    CHECK(!r.exactSkipped);
                    CHECK(SplitAnisoResidualOK(law,st[i].E,r));
                    if (lb==1 && r.exactIterations >= 3)
                        ++nIt3;
                    gsInfo << "  law " << (lb==1 ? "b" : "a") << " G " << ig << " (" << st[i].e1 << "," << st[i].e2 << "," << st[i].alpha
                           << "): iterations " << r.exactIterations << ", gamma " << r.gamma << ", theta " << r.theta << "\n";
                }
            }
        gsInfo << "  law (b) states with at least 3 iterations: " << nIt3 << "\n";
        CHECK(nIt3 >= 4);
    }

    /// A tangent argument that is indefinite along every direction (-C) is a failed direction: mode 0 skips the exact
    /// solve and returns S, mode 1 minimises from the eigenvector start and agrees with the run on the consistent tangent
    TEST(ExactMin_hard_forcedFailedDirection_lawB)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const real_t alphas[3] = { 0.2, 1.4, 2.5 };
        for (index_t ig = 0; ig!=2; ++ig)
        {
            const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
            const gsMatrix<real_t> G = J.transpose() * J;
            const SplitAniso law = SplitAnisoLaw(true,J);
            for (index_t a = 0; a!=3; ++a)
            {
                const gsMatrix<real_t> E = SplitAnisoE(J,0.05,-0.03,alphas[a]);
                CHECK_EQUAL(0,SplitAnisoClassify(law,E,G));
                const gsMatrix<real_t> negC = -law.C(E);

                const SplitXRes r0 = SplitAnisoCall(fx,law,E,G,SplitAnisoOpt(0),true,true,negC,true);
                CHECK(r0.tftFailure);
                CHECK(r0.exactSkipped);
                CHECK(!r0.exactMinimised);
                CHECK(!r0.exactFallback);
                CHECK((r0.S - law.S(E)).cwiseAbs().maxCoeff() == 0);

                const SplitXRes r1 = SplitAnisoCall(fx,law,E,G,SplitAnisoOpt(1),true,true,negC,true);
                CHECK(r1.tftFailure);
                CHECK(!r1.exactSkipped);
                CHECK(r1.exactMinimised);
                CHECK(!r1.exactFallback);
                CHECK_EQUAL(0,r1.state);
                CHECK(r1.exactIterations >= 3);
                CHECK(SplitAnisoResidualOK(law,E,r1));

                const SplitXRes rc = SplitAnisoRun(fx,law,E,G,1);
                CHECK_EQUAL(0,rc.state);
                real_t kappa, kappa1, tol1, tol2;
                SplitAnisoKappa(law,rc,kappa,tol2);
                SplitAnisoKappa(law,r1,kappa1,tol1);
                const real_t dS = (r1.S - rc.S).cwiseAbs().maxCoeff();
                gsInfo << "  G " << ig << ", alpha " << alphas[a] << ": iterations " << r1.exactIterations << ", |S1 - S1consistent| = " << dS
                       << " (bound " << 4 * kappa * math::max(tol1,tol2) << ")\n";
                CHECK(dS <= 4 * kappa * math::max(tol1,tol2));
            }
        }
    }

    /// Stiff fibre with the tangent argument scaled by f (the direction search is scale-invariant, the linearised
    /// gamma is multiplied by 1/f): the states where mode 0 fails and mode 1 is accepted
    TEST(ExactMin_hard_mode0Fails_lawB)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const real_t fs[2] = { 0.01, 0.003 };
        index_t nHits = 0, nCand = 0;
        bool allAccepted = true;
        for (index_t ig = 0; ig!=2; ++ig)
        {
            const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
            const gsMatrix<real_t> G = J.transpose() * J;
            const SplitAniso law = SplitAniso::NeoHookeFibre(J,10,100,0.35);
            const std::vector<SplitAnisoState> st = SplitAnisoStates(law,J,G);
            for (size_t i = 0; i!=st.size(); ++i)
            {
                if (st[i].kind!=0)
                    continue;
                for (index_t k = 0; k!=2; ++k)
                {
                    const gsMatrix<real_t> Cin = fs[k] * law.C(st[i].E);
                    const SplitXRes r0 = SplitAnisoCall(fx,law,st[i].E,G,SplitAnisoOpt(0),true,true,Cin,true);
                    const SplitXRes r1 = SplitAnisoCall(fx,law,st[i].E,G,SplitAnisoOpt(1),true,true,Cin,true);
                    ++nCand;
                    const bool fail0 = r0.exactFallback && !r0.exactWrongBranch;
                    const bool acc1 = r1.state==0 && !r1.exactFallback && r1.exactMinimised && SplitAnisoResidualOK(law,st[i].E,r1);
                    if (st[i].cls==0)
                        allAccepted = allAccepted && acc1;
                    if (fail0 && acc1)
                        ++nHits;
                    gsInfo << "  G " << ig << " (" << st[i].e1 << "," << st[i].e2 << "," << st[i].alpha << ") f " << fs[k]
                           << ": class " << st[i].cls << ", mode 0: state " << r0.state << " fallback " << r0.exactFallback << " saddle " << r0.exactSaddle
                           << " wrong " << r0.exactWrongBranch << " it " << r0.exactIterations << "; mode 1: state " << r1.state
                           << " fallback " << r1.exactFallback << " min " << r1.exactMinimised << " it " << r1.exactIterations << " accepted " << acc1 << "\n";
                }
            }
        }
        gsInfo << "  hard states (mode 0 nonconverged or saddle, mode 1 accepted): " << nHits << " of " << nCand << " candidates\n";
        CHECK(allAccepted);
        CHECK(nHits >= 1);
    }

    /// Modes 1 and 2 are bitwise equal at the anisotropic states: the closed-form direction is guard-rejected, so mode 1
    /// skips the scalar/2x2 solve (cfOK is false) and runs the same start rule and minimisation as mode 2. At isotropic and
    /// coaxial states mode 1 runs the scalar pass and the two modes agree to the conditioning of the solve.
    TEST(ExactMin_mode1_vs_mode2)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        index_t nAniso = 0;
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
                for (size_t i = 0; i!=st.size(); ++i)
                {
                    if (!SplitAnisoWrinkledCand(st[i]))
                        continue;
                    const SplitXRes r1 = SplitAnisoRun(fx,law,st[i].E,G,1);
                    const SplitXRes r2 = SplitAnisoRun(fx,law,st[i].E,G,2);
                    const std::string d = SplitAnisoDiff(r1,r2,true,true);
                    if (!d.empty())
                        gsInfo << "  law " << (lb==1 ? "b" : "a") << " G " << ig << " state " << i << ": modes 1 and 2 differ in" << d << "\n";
                    CHECK(d.empty());
                    CHECK(r1.guardFallback);
                    CHECK(r1.exactMinimised && r2.exactMinimised);
                    ++nAniso;
                }
            }
        gsInfo << "  anisotropic states compared bitwise: " << nAniso << "\n";

        const std::vector<SplitAnisoCase> cases = SplitAnisoIsoCases();
        for (size_t i = 0; i!=cases.size(); ++i)
        {
            const SplitAnisoCase & c = cases[i];
            const SplitXRes r1 = SplitAnisoRun(fx,c.law,c.E,c.G,1);
            const SplitXRes r2 = SplitAnisoRun(fx,c.law,c.E,c.G,2);
            CHECK_EQUAL(0,r1.state);
            CHECK_EQUAL(0,r2.state);
            CHECK(!r1.exactFallback && !r2.exactFallback);
            CHECK(r2.exactMinimised);
            real_t kappa, tolC1, tolC2, k2;
            SplitAnisoKappa(c.law,r1,kappa,tolC1);
            SplitAnisoKappa(c.law,r2,k2,tolC2);
            const real_t dS = (r1.S - r2.S).cwiseAbs().maxCoeff();
            const real_t bound = 4 * kappa * math::max(tolC1,tolC2);
            gsInfo << "  case " << c.id << " (" << c.e1 << "," << c.e2 << "," << c.alpha << "): minimised in mode 1 " << r1.exactMinimised
                   << ", |S1 - S2| = " << dS << " (bound " << bound << ")\n";
            CHECK(dS <= bound);
        }
    }

    /// Without a W provider mode 1 is mode 0, bit for bit, including the skipped failed direction
    TEST(ExactMin_noW_mode1_equals_mode0)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        index_t nCmp = 0;
        auto compare = [&](const SplitAniso & law, const gsMatrix<real_t> & E, const gsMatrix<real_t> & G, const gsMatrix<real_t> & Cin) -> SplitXRes
        {
            const SplitXRes r0 = SplitAnisoCall(fx,law,E,G,SplitAnisoOpt(0),true,false,Cin,false);
            const SplitXRes r1 = SplitAnisoCall(fx,law,E,G,SplitAnisoOpt(1),true,false,Cin,false);
            const std::string d = SplitAnisoDiff(r1,r0,true,true);
            if (!d.empty())
                gsInfo << "  modes 1 and 0 differ in" << d << "\n";
            CHECK(d.empty());
            CHECK(!r1.exactMinimised);
            ++nCmp;
            return r1;
        };
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
                for (size_t i = 0; i!=st.size(); ++i)
                    if (SplitAnisoWrinkledCand(st[i]))
                        compare(law,st[i].E,G,law.C(st[i].E));
                if (lb==1)
                {
                    const gsMatrix<real_t> E = SplitAnisoE(J,0.05,-0.03,0.2);
                    const SplitXRes r = compare(law,E,G,-law.C(E));
                    CHECK(r.exactSkipped);
                }
            }
        gsInfo << "  compared states: " << nCmp << "\n";
    }

    /// Where mode 0 needs no 2x2 pass and has an admissible direction, mode 1 is mode 0, bit for bit
    TEST(ExactMin_mode1_equals_mode0_without2x2)
    {
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const std::vector<SplitAnisoCase> cases = SplitAnisoIsoCases();
        index_t nCmp = 0;
        for (size_t i = 0; i!=cases.size(); ++i)
        {
            const SplitAnisoCase & c = cases[i];
            const SplitXRes r0 = SplitAnisoRun(fx,c.law,c.E,c.G,0);
            if (r0.exactUsed2x2 || r0.tftFailure)
                continue;
            const SplitXRes r1 = SplitAnisoRun(fx,c.law,c.E,c.G,1);
            const std::string d = SplitAnisoDiff(r1,r0,true,true);
            if (!d.empty())
                gsInfo << "  case " << c.id << " state " << i << ": modes 1 and 0 differ in" << d << "\n";
            CHECK(d.empty());
            CHECK(!r1.exactMinimised);
            ++nCmp;
        }
        gsInfo << "  compared states: " << nCmp << " of " << cases.size() << "\n";
        CHECK(nCmp >= 4);
    }

    /// Hes of Phi(v) = W(E + v (x) v) is positive definite iff n1.C'n1 > 0 and det Jac0 > 0 (the minimum test), with
    /// det Hes = 8 gamma det Jac0 and (Q^T Hes Q)_11 = 4 gamma M11 at R = 0
    TEST(ExactMin_HesPD_iff_minimumTest)
    {
        // (A) accepted wrinkled states of the exact mode
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        index_t nA = 0;
        for (index_t lb = 0; lb!=2; ++lb)
            for (index_t ig = 0; ig!=2; ++ig)
            {
                const gsMatrix<real_t> J = ig==0 ? SplitIdentity2() : SplitSkewedJref();
                const gsMatrix<real_t> G = J.transpose() * J;
                const SplitAniso law = SplitAnisoLaw(lb==1,J);
                const std::vector<SplitAnisoState> st = SplitAnisoStatesChecked(law,J,G);
                for (size_t i = 0; i!=st.size(); ++i)
                {
                    if (!SplitAnisoWrinkledCand(st[i]) && st[i].kind < 3)
                        continue;
                    const SplitXRes r = SplitAnisoRun(fx,law,st[i].E,G,1);
                    CHECK_EQUAL(0,r.state);
                    if (r.state!=0)
                        continue;
                    const real_t c = math::cos(r.theta), s = math::sin(r.theta);
                    const real_t v0 = math::sqrt(r.gamma) * c, v1 = math::sqrt(r.gamma) * s;
                    gsMatrix<real_t> Pv(3,1);
                    Pv << v0*v0, v1*v1, 2*v0*v1;
                    const gsMatrix<real_t> Ep = st[i].E + Pv;
                    const gsMatrix<real_t> Sp = law.S(Ep), Cp = law.C(Ep);
                    const gsMatrix<real_t> Hes = SplitAnisoHes(Cp,Sp,v0,v1);
                    real_t M11;
                    const gsMatrix<real_t> Jac = SplitAnisoJac0(Cp,Sp,r.theta,r.gamma,M11);
                    const real_t detJ = Jac.determinant();
                    gsVector<real_t,2> nn;
                    nn << c, s;
                    const real_t q11 = nn.dot(Hes * nn);
                    const real_t detH = Hes.determinant();
                    const bool pd = SplitAnisoLambdaMin(Hes) > 0;
                    gsInfo << "  law " << (lb==1 ? "b" : "a") << " G " << ig << " state " << i << ": det Hes / (8 gamma det Jac0) = "
                           << detH / (8 * r.gamma * detJ) << ", Q11 / (4 gamma M11) = " << q11 / (4 * r.gamma * M11)
                           << ", Hes PD " << pd << ", minimum test " << (M11 > 0 && detJ > 0) << "\n";
                    CHECK(math::abs(detH - 8 * r.gamma * detJ) <= 1e-8 * math::abs(8 * r.gamma * detJ));
                    CHECK(math::abs(q11 - 4 * r.gamma * M11) <= 1e-8 * math::abs(4 * r.gamma * M11));
                    CHECK_EQUAL(M11 > 0 && detJ > 0,pd);
                    ++nA;
                }
            }
        gsInfo << "  accepted states checked: " << nA << "\n";
        CHECK(nA >= 20);

        // (B) algebraic sweep at R = 0: n1.S' = n2.S' = 0, S' = sw Tm, C' = C0 + t w w^T
        index_t nPD = 0, nNotPD = 0;
        const real_t thetas[2] = { 0.3, 1.1 };
        const real_t gammas[2] = { 1e-3, 0.1 };
        const real_t sws[2] = { 0.5, -0.5 };
        const real_t ts[4] = { 0, -0.3, -3, -30 };
        for (index_t it = 0; it!=2; ++it)
            for (index_t ig = 0; ig!=2; ++ig)
                for (index_t is = 0; is!=2; ++is)
                    for (index_t k = 0; k!=4; ++k)
                    {
                        gsMatrix<real_t> v, w, Tm;
                        SplitDualVectors(thetas[it],v,w,Tm);
                        const gsMatrix<real_t> Sp = sws[is] * Tm;
                        const gsMatrix<real_t> Cp = SplitIsoC() + ts[k] * (w * w.transpose());
                        const real_t gam = gammas[ig];
                        const real_t v0 = math::sqrt(gam) * math::cos(thetas[it]), v1 = math::sqrt(gam) * math::sin(thetas[it]);
                        const gsMatrix<real_t> Hes = SplitAnisoHes(Cp,Sp,v0,v1);
                        real_t M11;
                        const gsMatrix<real_t> Jac = SplitAnisoJac0(Cp,Sp,thetas[it],gam,M11);
                        const real_t detJ = Jac.determinant();
                        const bool pd = SplitAnisoLambdaMin(Hes) > 0;
                        CHECK(math::abs(Hes.determinant() - 8 * gam * detJ) <= 1e-8 * math::abs(8 * gam * detJ));
                        CHECK_EQUAL(M11 > 0 && detJ > 0,pd);
                        if (pd) ++nPD; else ++nNotPD;
                    }
        gsInfo << "  sweep: " << nPD << " positive definite, " << nNotPD << " not\n";
        CHECK(nPD > 0);
        CHECK(nNotPD > 0);
    }

    /// The truncated (linearised) wrinkled tangent is asymmetric at the taut interface by O(gamma) when the wrapped
    /// tangent varies with E, and symmetric where it does not
    TEST(Linearised_asymmetry_gammaScaling)
    {
        const real_t h = 1e-6;
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitIdentity2(),1.0,1.0,0.0);
        const gsMatrix<real_t> I2 = SplitIdentity2();
        SplitXOpt opt = SplitExactOptions();
        opt.exact = false;
        opt.eta = 0;
        opt.slackMultiplier = 0;

        // Orthotropic law with phi0 = 0 at alpha = 0: the closed-form direction passes the guard and C is constant
        {
            const SplitAniso law = SplitAniso::Orthotropic(0);
            gsMatrix<real_t> E0, E1;
            const real_t sb = SplitAnisoTautBoundary(law,I2,I2,0.0,E0,E1);
            const gsMatrix<real_t> Eb = E0 + sb * (E1 - E0);
            const gsMatrix<real_t> Delta = SplitAnisoE(I2,0,-1,0.0);
            const real_t ts[2] = { 0.005, 0.01 };
            for (index_t k = 0; k!=2; ++k)
            {
                const gsMatrix<real_t> E = Eb + ts[k] * Delta;
                const SplitXRes r = fx.split->wrinkle(law.S(E),law.C(E),E,I2,opt,true);
                gsInfo << "  orthotropic, t " << ts[k] << ": asymmetry " << SplitAsymmetry(r.C) << "\n";
                CHECK_EQUAL(0,r.state);
                CHECK(!r.guardFallback);
                CHECK(SplitAsymmetry(r.C) <= 1e-12);
            }
        }

        // Pure neo-Hooke: dC/dE by central differences of the closed-form tangent
        {
            const SplitAniso law = SplitAniso::NeoHookeFibre(I2,0,20,0.35);
            gsMatrix<real_t> E0, E1;
            const real_t sb = SplitAnisoTautBoundary(law,I2,I2,0.9,E0,E1);
            const gsMatrix<real_t> Eb = E0 + sb * (E1 - E0);
            const gsMatrix<real_t> Delta = SplitAnisoE(I2,0,-1,0.9);
            auto run = [&](const real_t t, real_t & gamma) -> real_t
            {
                const gsMatrix<real_t> E = Eb + t * Delta;
                const SplitAniso L = law;
                const SplitXdCProv dC = [L,E,h](std::vector<gsMatrix<real_t> > & dCs)
                {
                    dCs.resize(3);
                    for (index_t d = 0; d!=3; ++d)
                    {
                        gsMatrix<real_t> Ep = E, Em = E;
                        Ep(d,0) += h;
                        Em(d,0) -= h;
                        dCs[d] = (L.C(Ep) - L.C(Em)) / (2*h);
                    }
                };
                const SplitXRes r = fx.split->wrinkle(law.S(E),law.C(E),E,I2,opt,true,dC);
                CHECK_EQUAL(0,r.state);
                CHECK(!r.guardFallback);
                gamma = r.gamma;
                return SplitAsymmetry(r.C);
            };
            real_t g0;
            run(0.01,g0);
            real_t g1, g2;
            const real_t t1 = 0.01 * 5e-3 / g0, t2 = 0.01 * 1e-2 / g0;
            const real_t a1 = run(t1,g1), a2 = run(t2,g2);
            const real_t rel = (a2 / a1) / (g2 / g1);
            gsInfo << "  neo-Hooke: gamma " << g1 << ", " << g2 << "; asymmetry " << a1 << ", " << a2 << "; (a2/a1)/(g2/g1) = " << rel << "\n";
            CHECK(math::abs(rel - 1) <= 0.2);
            CHECK(a1 >= 1e-7);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Effective energy W_eff against a central finite difference of the effective stress S_eff,
    /// dW_eff/dE_i = S_eff,i (engineering Voigt strain), over the three criterion states and eta in {0,1e-3}.
    /// The slack state at eta = 0 (zero slack multiplier) is identically zero and is checked as such.
    template<bool linear>
    static void SplitEnergyFDCheck(const SplitMaterial mat)
    {
        const SplitStateSpec specs[3] = { SplitWrinkledSpec(), SplitTautSpec(), SplitSlackSpec() };
        const real_t etas[2] = { 0.0, 1e-3 };
        const gsVector<real_t> pt = SplitPoint();
        const real_t z = 0;
        const real_t h = 1e-6;
        const real_t tol = 1e-6;
        const char * names[3] = { "wrinkled", "taut", "slack" };

        for (index_t s = 0; s!=3; ++s)
            for (index_t e = 0; e!=2; ++e)
            {
                gsInfo << "state " << names[s] << ", eta " << etas[e] << "\n";
                SplitFixture<2,linear> fx(mat,SplitSkewedJref(),specs[s].lam1,specs[s].lam2,static_cast<real_t>(EIGEN_PI)/6);
                fx.split->options().setReal("Eta",etas[e]);
                fx.split->options().setReal("SlackMultiplier",0);
                if (!linear)
                    fx.split->options().setSwitch("Exact",true);
                fx.split->resetCounters();

                const gsMatrix<real_t> Cmat = fx.CmatDef();
                const gsMatrix<real_t> E0 = fx.EengDef();
                CHECK_EQUAL(specs[s].state, fx.stateAt(Cmat,pt,z));

                const gsMatrix<real_t> S = fx.split->eval3D_vector_C(Cmat,0,pt,z);
                const real_t W0 = fx.split->eval3D_energy_C(Cmat,0,pt,z)(0,0);
                CHECK(S.allFinite());
                CHECK(math::isfinite(W0));

                if (specs[s].state==-1 && etas[e]==0)
                {
                    CHECK(W0 == 0);
                    CHECK(S.cwiseAbs().maxCoeff() == 0);
                }
                else
                {
                    const real_t Smax = S.cwiseAbs().maxCoeff();
                    CHECK(Smax > 0);
                    for (index_t i = 0; i!=3; ++i)
                    {
                        gsMatrix<real_t> Ep = E0, Em = E0;
                        Ep(i,0) += h;
                        Em(i,0) -= h;
                        const gsMatrix<real_t> Cp = fx.CmatFromE(Ep), Cm = fx.CmatFromE(Em);
                        CHECK_EQUAL(specs[s].state, fx.stateAt(Cp,pt,z));
                        CHECK_EQUAL(specs[s].state, fx.stateAt(Cm,pt,z));
                        const real_t Wp = fx.split->eval3D_energy_C(Cp,0,pt,z)(0,0);
                        const real_t Wm = fx.split->eval3D_energy_C(Cm,0,pt,z)(0,0);
                        CHECK(math::isfinite(Wp));
                        CHECK(math::isfinite(Wm));
                        const real_t fd = (Wp - Wm) / (2*h);
                        gsInfo << "  component " << i << ": S = " << S(i,0) << ", FD = " << fd
                               << ", rel.err = " << math::abs(S(i,0) - fd) / Smax << " (tol " << tol << ")\n";
                        CHECK(math::abs(S(i,0) - fd) <= tol * Smax);
                    }
                }
                CHECK_EQUAL(0,fx.split->numTFTFailures());
                if (!linear && specs[s].state==0)
                {
                    CHECK_EQUAL(0,fx.split->numExactFallbacks());
                    CHECK_EQUAL(0,fx.split->numExactReclassified());
                    CHECK_EQUAL(0,fx.split->numExactSkipped());
                }
            }
    }

    TEST(Energy_FD_class_C)
    {
        gsInfo << "NH incompressible, exact\n";
        SplitEnergyFDCheck<false>(SplitMatNHIncomp);
        gsInfo << "SvK (linear)\n";
        SplitEnergyFDCheck<true>(SplitMatSvK);
    }

    TEST(Energy_NaN_linearised_class)
    {
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
        fx.split->options().setSwitch("Exact",false);
        const gsMatrix<real_t> Cmat = fx.CmatDef();
        CHECK_EQUAL(0, fx.stateAt(Cmat,pt,0));
        const real_t W = fx.split->eval3D_energy_C(Cmat,0,pt,0)(0,0);
        gsInfo << "  W_eff (linearised wrinkled) = " << W << "\n";
        CHECK(math::isnan(W));
    }

    TEST(ExactLocalSolve_mode1_equals_mode0_class)
    {
        const SplitMaterial mats[2] = { SplitMatNHIncomp, SplitMatNHComp };
        const SplitStateSpec specs[3] = { SplitWrinkledSpec(), SplitTautSpec(), SplitSlackSpec() };
        const gsVector<real_t> pt = SplitPoint();
        index_t nCmp = 0;

        for (index_t m = 0; m!=2; ++m)
            for (index_t s = 0; s!=3; ++s)
            {
                SplitFixture<2,false> fx(mats[m],SplitSkewedJref(),specs[s].lam1,specs[s].lam2,static_cast<real_t>(EIGEN_PI)/6);
                fx.split->options().setSwitch("Exact",true);
                fx.split->options().setReal("Eta",1e-3);
                const gsMatrix<real_t> Cc = fx.CmatDef();

                fx.split->options().setInt("ExactLocalSolve",0);
                fx.split->resetCounters();
                const gsMatrix<real_t> S0 = fx.split->eval3D_vector_C(Cc,0,pt,0);
                const gsMatrix<real_t> C0 = fx.split->eval3D_matrix_C(Cc,0,pt,0);
                const index_t u2 = fx.split->numExactUsed2x2();
                const index_t tf = fx.split->numTFTFailures();

                fx.split->options().setInt("ExactLocalSolve",1);
                fx.split->resetCounters();
                const gsMatrix<real_t> S1 = fx.split->eval3D_vector_C(Cc,0,pt,0);
                const gsMatrix<real_t> C1 = fx.split->eval3D_matrix_C(Cc,0,pt,0);
                CHECK_EQUAL(0,fx.split->numExactMinimised());

                gsInfo << "  material " << m << ", state " << specs[s].state << ": 2x2 pass " << u2
                       << ", TFT failures " << tf << "\n";
                if (m==0 && s==0)
                    CHECK_EQUAL(0,u2);
                if (u2==0 && tf==0)
                {
                    SPLIT_CHECK_SAME(S1,S0);
                    SPLIT_CHECK_SAME(C1,C0);
                    ++nCmp;
                }
            }
        gsInfo << "  bitwise-compared cases: " << nCmp << " of 6\n";
        CHECK(nCmp >= 5);
    }

    TEST(ExactLocalSolve_mode2_vs_mode1_class)
    {
        const SplitMaterial mats[2] = { SplitMatNHIncomp, SplitMatNHComp };
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> ptm = pt;
        const gsMatrix<real_t> zmat = SplitZ(0);

        for (index_t m = 0; m!=2; ++m)
            for (index_t path = 0; path!=2; ++path)
            {
                gsInfo << "material " << m << (path==0 ? ", prescribed metric\n" : ", geometry\n");
                SplitFixture<2,false> fx(mats[m],SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
                fx.split->options().setSwitch("Exact",true);
                fx.split->options().setReal("Eta",1e-3);
                const gsMatrix<real_t> Cc = fx.CmatDef();

                gsMatrix<real_t> S[3];
                for (index_t mode = 1; mode!=3; ++mode)
                {
                    fx.split->options().setInt("ExactLocalSolve",mode);
                    fx.split->resetCounters();
                    S[mode] = (path==0) ? fx.split->eval3D_vector_C(Cc,0,pt,0.0)
                                        : fx.split->eval3D_vector(0,ptm,zmat,MaterialOutput::VectorN);
                }
                CHECK(fx.split->numExactMinimised() > 0);
                CHECK_EQUAL(0,fx.split->numExactFallbacks());
                CHECK_EQUAL(0,fx.split->numExactSkipped());
                SPLIT_CHECK_REL("mode 2 vs mode 1",S[2],S[1],1e-6);
            }
    }

    TEST(ExactLocalSolve_invalid_class)
    {
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> ptm = pt;
        const gsMatrix<real_t> zmat = SplitZ(0);

        {
            SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
            fx.split->options().setSwitch("Exact",true);
            fx.split->options().setInt("ExactLocalSolve",3);
            CHECK_THROW(fx.split->eval3D_vector_C(fx.CmatDef(),0,pt,0.0), std::runtime_error);
        }
        {
            SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
            fx.split->options().setSwitch("Exact",false);
            fx.split->options().setInt("ExactLocalSolve",3);
            CHECK_THROW(fx.split->eval3D_vector_C(fx.CmatDef(),0,pt,0.0), std::runtime_error);
        }
        {
            SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
            fx.split->options().setSwitch("Exact",false);
            fx.split->options().setInt("ExactLocalSolve",3);
            CHECK_THROW(fx.split->eval3D_vector(0,ptm,zmat,MaterialOutput::VectorN), std::runtime_error);
        }
    }

    /// Eta outside [0,1] and a negative SlackMultiplier are not validated by setOptions; they throw at the first
    /// wrinkling evaluation (wrinkleOptions()), each with all other options at their defaults.
    TEST(WrinkleOptions_slackmultiplier_negative_throws)
    {
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> ptm = pt;
        const gsMatrix<real_t> zmat = SplitZ(0);

        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
        fx.split->options().setReal("SlackMultiplier",-1e-3);
        CHECK_THROW(fx.split->eval3D_vector(0,ptm,zmat,MaterialOutput::VectorN), std::runtime_error);
        CHECK_THROW(fx.split->wrinkleOptions(), std::runtime_error);
    }

    TEST(WrinkleOptions_eta_negative_throws)
    {
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> ptm = pt;
        const gsMatrix<real_t> zmat = SplitZ(0);

        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
        fx.split->options().setReal("Eta",-0.1);
        CHECK_THROW(fx.split->eval3D_vector(0,ptm,zmat,MaterialOutput::VectorN), std::runtime_error);
        CHECK_THROW(fx.split->wrinkleOptions(), std::runtime_error);
    }

    TEST(WrinkleOptions_eta_above_one_throws)
    {
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> ptm = pt;
        const gsMatrix<real_t> zmat = SplitZ(0);

        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
        fx.split->options().setReal("Eta",1.1);
        CHECK_THROW(fx.split->eval3D_vector(0,ptm,zmat,MaterialOutput::VectorN), std::runtime_error);
        CHECK_THROW(fx.split->wrinkleOptions(), std::runtime_error);
    }

    /// eval3D_pstretch and eval3D_detF forward to the wrapped law and never read the wrinkling options.
    TEST(WrinkleOptions_passthroughs_ignore_bad_options)
    {
        const SplitStateSpec spec = SplitWrinkledSpec();
        const gsVector<real_t> pt = SplitPoint();
        const gsMatrix<real_t> ptm = pt;
        const gsMatrix<real_t> zmat = SplitZ(0);

        SplitFixture<2,false> fx(SplitMatNHIncomp,SplitSkewedJref(),spec.lam1,spec.lam2,static_cast<real_t>(EIGEN_PI)/6);
        CHECK(fx.split->eval3D_pstretch(0,ptm,zmat).allFinite());
        CHECK(fx.split->eval3D_detF(0,ptm,zmat,MaterialOutput::Generic).allFinite());

        fx.split->options().setReal("SlackMultiplier",-1e-3);
        fx.split->options().setReal("Eta",1.1);
        CHECK(fx.split->eval3D_pstretch(0,ptm,zmat).allFinite());
        CHECK(fx.split->eval3D_detF(0,ptm,zmat,MaterialOutput::Generic).allFinite());
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#   undef SPLIT_CHECK_SAME
#   undef SPLIT_CHECK_REL
}
