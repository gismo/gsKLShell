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

    /// The provider with the wrinkle tension n4.S' shifted to -sw_lin is reclassified as slack and returns
    /// sigma (S, C, E), sigma = SlackMultiplier if positive, else Eta
    TEST(Exact_synthetic_reclassify_slack)
    {
        SplitSynthetic sy(false);
        const SplitSynthetic::Provider prov = sy.provider(SplitSynthetic::Slack,2 * sy.swLin);
        const real_t sigmas[2] = { 1e-3, 0.25 };
        for (index_t i = 0; i!=2; ++i)
            for (index_t t = 0; t!=2; ++t)
            {
                SplitSynthetic::Options opt = sy.optEx;
                opt.eta = 1e-3;
                opt.slackMultiplier = (i==0) ? 0.0 : 0.25;
                sy.fx.split->resetCounters();
                const SplitSynthetic::Result r = sy.run(opt,t==0,prov);
                CHECK_EQUAL(-1,r.state);
                CHECK(r.exactReclassified);
                CHECK(!r.exactFallback);
                CHECK_EQUAL(1,sy.fx.split->numExactReclassified());
                CHECK_EQUAL(0,sy.fx.split->numExactFallbacks());
                const gsMatrix<real_t> Se = sigmas[i] * sy.S;
                const gsMatrix<real_t> Ee = sigmas[i] * sy.E;
                SPLIT_CHECK_SAME(r.S,Se);
                SPLIT_CHECK_SAME(r.E,Ee);
                if (t==0)
                {
                    const gsMatrix<real_t> Ce = sigmas[i] * sy.C;
                    SPLIT_CHECK_SAME(r.C,Ce);
                }
                else
                    CHECK_EQUAL(0,r.C.size());
                CHECK(std::isnan(r.theta) && std::isnan(r.gamma));
            }

        // Control: without the shift the wrinkle tension is positive and the point stays wrinkled
        sy.fx.split->resetCounters();
        const SplitSynthetic::Result rc = sy.run(sy.optEx,true,sy.provider(SplitSynthetic::Slack,0.0));
        CHECK_EQUAL(0,rc.state);
        CHECK(!rc.exactReclassified);
        CHECK(!rc.exactFallback);
        CHECK_EQUAL(0,sy.fx.split->numExactReclassified());
    }

    /// A provider that returns NaN, or a tangent with n1^T C n1 < 0 (scalar mode), cannot be solved: the first
    /// provider call ends the solve and the result is the linearised one, bit for bit, with the fallback counted
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

#   undef SPLIT_CHECK_SAME
#   undef SPLIT_CHECK_REL
}
