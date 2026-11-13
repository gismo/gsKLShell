/** @file gsMaterialMatrixWrinkleOutputs_test.cpp

    @brief Unit tests of the post-processing outputs of gsMaterialMatrixTFT and gsMaterialMatrixTFTSplit:
    principal stress frame, wrinkle direction, wrinkling intensity, parametrisation-independent wrinkle
    strain, the zero/NaN conventions, the modified-stress outputs and gsShellStressFunction.

    The oracles do not reuse the helper code of the classes:
    - analytic states: affine deformations x = Q R diag(lam1,lam2) R^T X with a known principal strain
      E_min = (lam2^2 - 1)/2 and a known direction Q R e1 along the wrinkles, embedded in 3D by rigid rotations,
    - Jacobians of the patches (gsMapData) for the reference and deformed bases,
    - symmetric eigen solves of Cartesian stress tensors and generalized solves of metric pencils,
    - the wrapped material evaluated on its own at a single point,
    - the Voigt direction vectors n1 = (c^2,s^2,2cs), n2 = (-sc,sc,c^2-s^2) of the wrinkle covector.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

#include <cmath>
#include <iomanip>
#include <limits>

using namespace gismo;

SUITE(gsMaterialMatrixWrinkleOutputs_test)                 // The suite should have the same name as the file
{
    /// Largest entry-wise difference of \a A and \a B relative to the largest entry of the reference \a B
    static real_t WORelErr(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
    {
        return (A-B).cwiseAbs().maxCoeff() / B.cwiseAbs().maxCoeff();
    }

    /// Relative comparison of two matrices (reference \a B) with the achieved error printed
#   define WO_CHECK_REL(label,A,B,tol) \
    do { \
        gsInfo << "  " << (label) << ": rel.err = " << WORelErr((A),(B)) << " (tol " << (tol) << ")\n"; \
        CHECK((A).allFinite()); \
        CHECK((B).allFinite()); \
        CHECK(((A)-(B)).cwiseAbs().maxCoeff() <= (tol) * (B).cwiseAbs().maxCoeff()); \
    } while(0)

    /// Absolute difference of two angles modulo pi, folded into [0,pi/2]
    static real_t WOThetaDiff(const real_t a, const real_t b)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        real_t d = std::fmod(a - b, pi);
        if (d >  pi/2) d -= pi;
        if (d < -pi/2) d += pi;
        return math::abs(d);
    }

    /// 1x1 through-thickness coordinate matrix
    static gsMatrix<real_t> WOZ(const real_t z)
    {
        gsMatrix<real_t> zmat(1,1);
        zmat(0,0) = z;
        return zmat;
    }

    /// Distance of two vectors that are equal up to a sign
    static real_t WOSignedDiff(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
    {
        return math::min((a-b).norm(),(a+b).norm());
    }

    /// Rotation by \a angle about the coordinate axis \a axis (0: x, 1: y, 2: z)
    static gsMatrix<real_t> WORot(const index_t axis, const real_t angle)
    {
        const real_t c = math::cos(angle), s = math::sin(angle);
        gsMatrix<real_t> R(3,3);
        R.setIdentity();
        if (axis==0)      { R(1,1) = c; R(1,2) = -s; R(2,1) = s; R(2,2) = c; }
        else if (axis==1) { R(0,0) = c; R(0,2) =  s; R(2,0) = -s; R(2,2) = c; }
        else              { R(0,0) = c; R(0,1) = -s; R(1,0) = s; R(1,1) = c; }
        return R;
    }

    /// Planar rotation by \a angle
    static gsMatrix<real_t> WORot2(const real_t angle)
    {
        gsMatrix<real_t> R(2,2);
        R << math::cos(angle), -math::sin(angle), math::sin(angle), math::cos(angle);
        return R;
    }

    /// Rotation of the embedded reference patches in 3D
    static gsMatrix<real_t> WORref()
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        return WORot(2,20*pi/180) * WORot(0,35*pi/180);
    }

    /// Rotation of the embedded deformed patches in 3D: the deformed plane is tilted by 25 degrees w.r.t. the reference plane
    static gsMatrix<real_t> WORdef()
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        return WORot(0,-25*pi/180) * WORref();
    }

    /// Voigt direction n1 = (c^2,s^2,2cs) of the wrinkle covector (cos theta, sin theta)
    static gsVector<real_t,3> WON1(const real_t theta)
    {
        const real_t c = math::cos(theta), s = math::sin(theta);
        gsVector<real_t,3> n;
        n << c*c, s*s, 2*c*s;
        return n;
    }

    /// Voigt direction n2 = (-sc,sc,c^2-s^2) of the wrinkle covector (cos theta, sin theta)
    static gsVector<real_t,3> WON2(const real_t theta)
    {
        const real_t c = math::cos(theta), s = math::sin(theta);
        gsVector<real_t,3> n;
        n << -s*c, s*c, c*c-s*s;
        return n;
    }

    enum WOMaterial { WOMatSvK, WOMatNHComp, WOMatNHIncomp };

    static real_t WOPoisson(const WOMaterial m)
    {
        return (m==WOMatSvK || m==WOMatNHComp) ? 0.45 : 0.5;
    }

    /// Young's modulus from the shear modulus mu = 4.225e5
    static real_t WOYoung(const WOMaterial m)
    {
        return 2 * 4.225e5 * (1 + WOPoisson(m));
    }

    /// Principal stretches of an affine state and the criterion state they give (1 taut, 0 wrinkled, -1 slack)
    struct WOStateSpec
    {
        real_t lam1, lam2;
        index_t state;
    };

    static WOStateSpec WOWrinkledSpec() { WOStateSpec s = { 1.2, 0.8 / math::sqrt(1.2), 0 }; return s; }
    /// Taut although the minor Green-Lagrange strain is negative: the free contraction 1.2^(-1/2) = 0.913 is below lam2
    static WOStateSpec WOTautSpec()     { WOStateSpec s = { 1.2, 0.97, 1 }; return s; }
    static WOStateSpec WOSlackSpec()    { WOStateSpec s = { 0.95, 0.90, -1 }; return s; }

    /// Minor principal Green-Lagrange strain of an affine state: the rotations drop out of F^T F
    static real_t WOEmin(const WOStateSpec & s)
    {
        return 0.5 * (s.lam2*s.lam2 - 1);
    }

    /// Bilinear unit-square patch mapped by the 2x2 Jacobian \a J (x = J X)
    static gsMultiPatch<real_t> WOAffinePatch(const gsMatrix<real_t> & J)
    {
        gsMultiPatch<real_t> mp;
        mp.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        mp.patch(0).coefs() = mp.patch(0).coefs() * J.transpose();
        mp.addAutoBoundaries();
        return mp;
    }

    /// Reference Jacobian with g12 != 0 in the reference metric
    static gsMatrix<real_t> WOSkewedJref()
    {
        gsMatrix<real_t> Jref(2,2);
        Jref << 1.0, 0.3, 0.0, 0.8;
        return Jref;
    }

    static gsMatrix<real_t> WOIdentity2()
    {
        gsMatrix<real_t> I(2,2);
        I.setIdentity();
        return I;
    }

    /// Deformed Jacobian Q(beta) R(alpha) diag(lam1,lam2) R(alpha)^T Jref: principal stretches along the axes rotated by
    /// \a alpha, followed by the rigid rotation \a beta
    static gsMatrix<real_t> WOJdef(const gsMatrix<real_t> & Jref, const real_t lam1, const real_t lam2, const real_t alpha, const real_t beta)
    {
        gsMatrix<real_t> D(2,2);
        D << lam1, 0.0, 0.0, lam2;
        const gsMatrix<real_t> R = WORot2(alpha);
        return WORot2(beta) * R * D * R.transpose() * Jref;
    }

    /// Reference patch of an affine pair: planar in dim 2, embedded and rotated by WORref() in dim 3
    template<short_t dim>
    static gsMultiPatch<real_t> WOAffineRef(const gsMatrix<real_t> & Jref)
    {
        gsMultiPatch<real_t> mp = WOAffinePatch(Jref);
        if (dim==3)
        {
            mp.embed(3);
            mp.patch(0).coefs() = mp.patch(0).coefs() * WORref().transpose();
        }
        return mp;
    }

    /// Deformed patch of an affine pair: planar in dim 2, embedded and rotated by WORdef() in dim 3
    template<short_t dim>
    static gsMultiPatch<real_t> WOAffineDef(const gsMatrix<real_t> & Jref, const real_t lam1, const real_t lam2, const real_t alpha, const real_t beta)
    {
        gsMultiPatch<real_t> mp = WOAffinePatch(WOJdef(Jref,lam1,lam2,alpha,beta));
        if (dim==3)
        {
            mp.embed(3);
            mp.patch(0).coefs() = mp.patch(0).coefs() * WORdef().transpose();
        }
        return mp;
    }

    /// 3x2 basis (columns = tangent vectors) of the first patch at \a pt; padded with a zero row in dim 2
    static gsMatrix<real_t> WOJacobian(const gsMultiPatch<real_t> & mp, const gsVector<real_t> & pt)
    {
        gsMapData<real_t> map;
        map.flags = NEED_JACOBIAN | NEED_VALUE;
        map.points = pt;
        static_cast<const gsFunction<real_t>&>(mp.patch(0)).computeMap(map);
        const gsMatrix<real_t> J = map.jacobian(0);
        gsMatrix<real_t> B(3,2);
        B.setZero();
        B.topRows(J.rows()) = J;
        return B;
    }

    /// Cartesian 3x3 tensor sum_{ab} S^{ab} b_a b_b^T of the contravariant components \a S = (S11,S22,S12) in the basis \a basis
    static gsMatrix<real_t> WOCartesianTensor(const gsMatrix<real_t> & S, const gsMatrix<real_t> & basis)
    {
        gsMatrix<real_t> Smat(2,2);
        Smat << S(0,0), S(2,0), S(2,0), S(1,0);
        return basis * Smat * basis.transpose();
    }

    /// Unit normal g1 x g2 / |g1 x g2| of a 3x2 basis
    static gsVector<real_t,3> WONormal(const gsMatrix<real_t> & basis)
    {
        gsVector<real_t,3> g1 = basis.col(0), g2 = basis.col(1);
        return g1.cross(g2).normalized();
    }

    /// Single-point evaluation point of the affine tests
    static gsVector<real_t> WOPoint()
    {
        gsVector<real_t> pt(2);
        pt.setConstant(0.25);
        return pt;
    }

    /**
     * Owns a wrapped material, its parameter functions, the patches and both decorators acting on
     * the same wrapped object. The parameter functions are members because the nonlinear materials
     * keep references to them. Affine pairs also carry the metrics, the analytic wrinkle direction and
     * the reference and deformed unit normals.
     */
    template<short_t dim, bool linear>
    struct WOFixture
    {
        typedef gsMaterialMatrixTFTSplit<dim,real_t,linear> SplitType;
        typedef gsMaterialMatrixTFT<dim,real_t,linear>      TFTType;

        gsMultiPatch<real_t> mp, mp_def;
        gsConstantFunction<real_t> t, E, nu;
        gsMaterialMatrixBase<real_t>::uPtr wrapped;
        memory::unique_ptr<SplitType> split;
        memory::unique_ptr<TFTType>   tft;

        /// Affine pairs only: reference 2x2 metric, deformed 2x2 metric, reference (G11,G22,G12)
        gsMatrix<real_t> Gref, Gdef, gori;
        /// Affine pairs only: analytic Cartesian direction along the wrinkles on the deformed surface
        gsMatrix<real_t> dExact;
        /// Unit normals of the reference and deformed planes (affine pairs)
        gsMatrix<real_t> Nref, Ndef;

        WOFixture(const WOMaterial mat, const gsMultiPatch<real_t> & ref, const gsMultiPatch<real_t> & def, const real_t thickness)
        :
        mp(ref), mp_def(def),
        t(thickness,dim), E(WOYoung(mat),dim), nu(WOPoisson(mat),dim)
        {
            std::vector<gsFunctionSet<real_t>*> parameters(2);
            parameters[0] = &E;
            parameters[1] = &nu;

            gsOptionList options;
            options.addInt("Material","Material model: (0): SvK | (1): NH | (3): MR",mat==WOMatSvK ? 0 : 1);
            options.addInt("Implementation","Implementation: (1): Analytical",1);
            if (mat!=WOMatSvK)
                options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",mat==WOMatNHComp);
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

        /// Affine pair: Jref and the deformed Jacobian Q(beta) R(alpha) diag(lam1,lam2) R(alpha)^T Jref
        WOFixture(const WOMaterial mat, const gsMatrix<real_t> & Jref, const real_t lam1, const real_t lam2, const real_t alpha, const real_t beta, const real_t thickness = 0.01)
        :
        WOFixture(mat,WOAffineRef<dim>(Jref),WOAffineDef<dim>(Jref,lam1,lam2,alpha,beta),thickness)
        {
            const gsMatrix<real_t> Jdef = WOJdef(Jref,lam1,lam2,alpha,beta);
            Gref = Jref.transpose() * Jref;
            Gdef = Jdef.transpose() * Jdef;
            gori.resize(3,1);
            gori << Gref(0,0), Gref(1,1), Gref(0,1);

            gsMatrix<real_t> d2(2,1);
            d2 = WORot2(beta) * WORot2(alpha).col(0);
            dExact.setZero(3,1);
            dExact.topRows(2) = d2;
            Nref.setZero(3,1);
            Nref(2,0) = 1;
            Ndef = Nref;
            if (dim==3)
            {
                dExact = WORdef() * dExact;
                Nref = WORref() * Nref;
                Ndef = WORdef() * Ndef;
            }
        }

        /// Engineering Voigt strain of the deformed state (affine pairs)
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

        /// Decorator 0 = TFT, 1 = Split
        gsMaterialMatrixBase<real_t> * decorator(const index_t i) const
        {
            return i==0 ? static_cast<gsMaterialMatrixBase<real_t>*>(tft.get()) : static_cast<gsMaterialMatrixBase<real_t>*>(split.get());
        }
    };

    static const char * WOName(const index_t i)
    {
        return i==0 ? "TFT" : "Split";
    }

    static const char * WOMatName(const WOMaterial m)
    {
        return m==WOMatSvK ? "SvK" : (m==WOMatNHComp ? "NH comp" : "NH incomp");
    }

    /// Criterion state of the decorator at one point (1 taut, 0 wrinkled, -1 slack), read from eval3D_tensionfield
    static index_t WOStateOf(const gsMaterialMatrixBase<real_t> & m, const gsVector<real_t> & pt, const real_t z)
    {
        return static_cast<index_t>(m.eval3D_tensionfield(0,pt,WOZ(z),MaterialOutput::Generic)(0,0));
    }

    /// Orthonormality, orientation and the normal of the 3x3 frame stored column-major in the 9x1 vector \a D.
    /// Returns false (after failing a CHECK) if \a D has the wrong shape or non-finite entries.
    static bool WOUnpackFrame(const gsMatrix<real_t> & D, const gsMatrix<real_t> & basis, gsMatrix<real_t> & F, real_t & orthErr, real_t & detErr, real_t & nErr)
    {
        CHECK_EQUAL(9,D.rows());
        CHECK_EQUAL(1,D.cols());
        if (D.rows()!=9 || D.cols()!=1)
            return false;
        CHECK(D.allFinite());
        if (!D.allFinite())
            return false;
        F.resize(3,3);
        for (index_t l = 0; l!=3; ++l)
            for (index_t i = 0; i!=3; ++i)
                F(i,l) = D(i+3*l,0);
        gsMatrix<real_t> I3(3,3);
        I3.setIdentity();
        orthErr = (F.transpose()*F - I3).cwiseAbs().maxCoeff();
        detErr = math::abs(F.determinant() - 1);
        const gsVector<real_t,3> N = WONormal(basis);
        nErr = (F.col(2) - N).cwiseAbs().maxCoeff();
        return true;
    }

    /// Frame [d_min | d_max | N] of one evaluation (9x1 \a D) against the stress \a S (3x1 contravariant components in the
    /// reference basis \a basis) and its principal values \a P (2x1): orthonormal, right-handed, N the reference normal,
    /// T d = P d for the Cartesian tensor T of S (skipped for S = 0), and, at wrinkled points, a vanishing minor value
    /// whose direction stays tangent.
    static void WOCheckFrameColumn(const char * label, const gsMatrix<real_t> & D, const gsMatrix<real_t> & S, const gsMatrix<real_t> & P,
                                   const gsMatrix<real_t> & basis, const bool wrinkled)
    {
        gsMatrix<real_t> F;
        real_t orthErr = 0, detErr = 0, nErr = 0;
        if (!WOUnpackFrame(D,basis,F,orthErr,detErr,nErr))
            return;
        CHECK(P.allFinite());
        CHECK(S.allFinite());
        CHECK_EQUAL(2,P.rows());
        if (P.rows()!=2 || !P.allFinite() || !S.allFinite())
            return;
        gsInfo << "  " << label << ": orth " << orthErr << ", det " << detErr << ", N " << nErr;
        CHECK(orthErr <= 1e-12);
        CHECK(detErr <= 1e-12);
        CHECK(nErr <= 1e-12);
        CHECK(P(0,0) <= P(1,0));

        real_t pairErr = 0;
        if (S.norm() > 0)
        {
            const gsMatrix<real_t> T = WOCartesianTensor(S,basis);
            const real_t scale = P.cwiseAbs().maxCoeff();
            for (index_t l = 0; l!=2; ++l)
            {
                const gsMatrix<real_t> d = F.col(l);
                pairErr = math::max(pairErr,(T*d - P(l,0)*d).norm() / scale);
            }
            CHECK(pairErr <= 1e-10);
        }
        gsInfo << ", T d - P d " << pairErr;

        if (wrinkled)
        {
            const gsVector<real_t,3> N = WONormal(basis);
            const real_t minor = math::abs(P(0,0)) / math::abs(P(1,0));
            const real_t tang  = math::abs(F.col(0).dot(N));
            gsInfo << ", |P0|/|P1| " << minor << ", d_min.N " << tang;
            CHECK(minor <= 1e-8);
            CHECK(tang <= 1e-12);
        }
        gsInfo << "\n";
    }

    /// Frame, principal values and stress of decorator \a m at a single point
    static void WOFrameAtPoint(const char * label, const gsMaterialMatrixBase<real_t> & m, const gsMultiPatch<real_t> & mp,
                               const gsVector<real_t> & pt, const real_t z)
    {
        const gsMatrix<real_t> D = m.eval3D_pstressDir(0,pt,WOZ(z),MaterialOutput::Generic);
        const gsMatrix<real_t> P = m.eval3D_pstress(0,pt,WOZ(z),MaterialOutput::Generic);
        const gsMatrix<real_t> S = m.eval3D_stress(0,pt,WOZ(z),MaterialOutput::VectorN);
        WOCheckFrameColumn(label,D,S,P,WOJacobian(mp,pt),WOStateOf(m,pt,z)==0);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Frame cases of one material and state in 2D, for both decorators; \a eta and \a slack set the Split Eta and the TFT
    /// SlackMultiplier when \a setOptions is true
    template<bool linear>
    static void WOFrame2DCase(const WOMaterial mat, const WOStateSpec & st, const char * name, const bool setOptions, const real_t eta, const real_t slack)
    {
        WOFixture<2,linear> fx(mat,WOSkewedJref(),st.lam1,st.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        if (setOptions)
        {
            fx.split->options().setReal("Eta",eta);
            fx.tft->options().setReal("SlackMultiplier",slack);
        }
        const gsVector<real_t> pt = WOPoint();
        for (index_t c = 0; c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
            const index_t state = WOStateOf(m,pt,0);
            gsInfo << WOName(c) << ", " << WOMatName(mat) << ", " << name;
            if (setOptions)
                gsInfo << ", eta " << eta << ", slack multiplier " << slack;
            gsInfo << ": state " << state << "\n";
            CHECK_EQUAL(st.state,state);
            WOFrameAtPoint(name,m,fx.mp,pt,0);
        }
    }

    /// Principal stress frame of both decorators in 2D: wrinkled, taut with a negative minor strain, and slack
    TEST(Frame_2D)
    {
        const WOStateSpec specs[3] = { WOWrinkledSpec(), WOTautSpec(), WOSlackSpec() };
        const char * names[3] = { "wrinkled", "taut", "slack" };
        for (index_t s = 0; s!=3; ++s)
        {
            WOFrame2DCase<false>(WOMatNHComp,specs[s],names[s],false,0,0);
            WOFrame2DCase<false>(WOMatNHIncomp,specs[s],names[s],false,0,0);
            WOFrame2DCase<true>(WOMatSvK,specs[s],names[s],false,0,0);
        }
        // Non-zero slack stress: Split with a regularisation, TFT with a slack multiplier
        WOFrame2DCase<false>(WOMatNHComp,WOSlackSpec(),"slack",true,1e-3,0.25);
        WOFrame2DCase<false>(WOMatNHIncomp,WOSlackSpec(),"slack",true,1e-3,0.25);
        WOFrame2DCase<true>(WOMatSvK,WOSlackSpec(),"slack",true,1e-3,0.25);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Flat unit square embedded in 3D
    static gsMultiPatch<real_t> WOCurvedRef()
    {
        gsMultiPatch<real_t> ref;
        ref.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1) );
        ref.addAutoBoundaries();
        ref.embed(3);
        return ref;
    }

    /// The flat square stretched in x by 1.2, contracted in y by 0.75, with one corner lifted (the metric varies through the
    /// thickness), optionally followed by a rigid rotation about the x-axis
    static gsMultiPatch<real_t> WOCurvedDef(const bool rotate)
    {
        gsMultiPatch<real_t> def = WOCurvedRef();
        def.patch(0).coefs().col(0) *= 1.2;
        def.patch(0).coefs().col(1) *= 0.75;
        def.patch(0).coefs()(3,2) += 0.2;
        if (rotate)
            def.patch(0).coefs() = def.patch(0).coefs() * WORot(0,0.3).transpose();
        return def;
    }

    /// Four points and, per point, two fibre coordinates
    static void WOCurvedUZ(gsMatrix<real_t> & u, gsMatrix<real_t> & z)
    {
        u.resize(2,4);
        z.resize(2,4);
        u << 0.2, 0.7, 0.4, 0.85,
             0.3, 0.2, 0.8, 0.65;
        for (index_t k = 0; k!=4; ++k)
        {
            z(0,k) = -0.5 + 0.1 * k;
            z(1,k) =  0.5 - 0.1 * k;
        }
    }

    /// Column j*u.cols()+k of a multi-point frame equals the single-point call at (u.col(k), z(j,k)), and every column is a valid
    /// frame of its single-point stress (reference basis from the patch \a mp)
    static void WOMultipointFrame(const char * label, const gsMaterialMatrixBase<real_t> & m, const gsMultiPatch<real_t> & mp,
                                  const gsMatrix<real_t> & u, const gsMatrix<real_t> & z, const gsMatrix<real_t> * Nexp)
    {
        const gsMatrix<real_t> D = m.eval3D_pstressDir(0,u,z,MaterialOutput::Generic);
        CHECK_EQUAL(9,D.rows());
        CHECK_EQUAL(u.cols()*z.rows(),D.cols());
        if (D.rows()!=9 || D.cols()!=u.cols()*z.rows())
            return;
        real_t layoutErr = 0, nErr = 0;
        for (index_t j = 0; j!=z.rows(); ++j)
            for (index_t k = 0; k!=u.cols(); ++k)
            {
                const gsVector<real_t> uk = u.col(k);
                const index_t col = j*u.cols()+k;
                const gsMatrix<real_t> Ds = m.eval3D_pstressDir(0,uk,WOZ(z(j,k)),MaterialOutput::Generic);
                const gsMatrix<real_t> Ps = m.eval3D_pstress(0,uk,WOZ(z(j,k)),MaterialOutput::Generic);
                const gsMatrix<real_t> Ss = m.eval3D_stress(0,uk,WOZ(z(j,k)),MaterialOutput::VectorN);
                const gsMatrix<real_t> Dc = D.col(col);
                CHECK_EQUAL(9,Ds.rows());
                if (Ds.rows()!=9)
                    return;
                CHECK(Ds.allFinite());
                layoutErr = math::max(layoutErr,(Dc-Ds).cwiseAbs().maxCoeff());
                if (Nexp)
                    nErr = math::max(nErr,(Dc.bottomRows(3)-*Nexp).cwiseAbs().maxCoeff());
                WOCheckFrameColumn(label,Dc,Ss,Ps,WOJacobian(mp,uk),WOStateOf(m,uk,z(j,k))==0);
            }
        gsInfo << "  " << label << ": multi-point vs single-point " << layoutErr << ", N vs reference normal " << nErr << "\n";
        CHECK(layoutErr <= 1e-12);
        CHECK(nErr <= 1e-12);
    }

    /// The frame of both decorators on multi-point calls: N is the reference normal (the deformed plane is tilted
    /// w.r.t. the reference plane), each column equals its single-point call and is an eigenframe of its stress
    TEST(Frame_3D_multipoint)
    {
        const WOStateSpec w = WOWrinkledSpec();
        WOFixture<3,false> fx(WOMatNHComp,WOSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        const real_t dots = math::abs(fx.Nref.col(0).dot(fx.Ndef.col(0)));
        gsInfo << "  |N_ref . N_def| = " << dots << "\n";
        CHECK(dots < 1 - 1e-3);

        gsMatrix<real_t> u(2,3), z(2,3);
        u << 0.2, 0.5, 0.8,
             0.3, 0.6, 0.1;
        z << 0.0, 0.0, 0.0,
             0.25, 0.25, 0.25;
        for (index_t c = 0; c!=2; ++c)
            WOMultipointFrame((std::string("affine ") + WOName(c)).c_str(),*fx.decorator(c),fx.mp,u,z,&fx.Nref);

        // Curved pair, deformed plane rotated: the reference is flat, so its normal is the same everywhere
        WOFixture<3,false> cv(WOMatNHComp,WOCurvedRef(),WOCurvedDef(true),0.1);
        gsMatrix<real_t> uc, zc;
        WOCurvedUZ(uc,zc);
        gsMatrix<real_t> Nflat(3,1);
        Nflat << 0, 0, 1;
        for (index_t c = 0; c!=2; ++c)
            WOMultipointFrame((std::string("curved ") + WOName(c)).c_str(),*cv.decorator(c),cv.mp,uc,zc,&Nflat);
    }

    /// At a taut point S_eff is the wrapped stress, so the decorator frame equals the wrapped material's own frame
    /// (up to the sign of the eigenvectors; the normal exactly)
    template<short_t dim>
    static void WOTautWrappedCase(const WOMaterial mat)
    {
        const WOStateSpec st = WOTautSpec();
        WOFixture<dim,false> fx(mat,WOSkewedJref(),st.lam1,st.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        const gsVector<real_t> pt = WOPoint();
        const gsMatrix<real_t> Dw = fx.wrapped->eval3D_pstressDir(0,pt,WOZ(0),MaterialOutput::Generic);
        CHECK_EQUAL(9,Dw.rows());
        CHECK_EQUAL(1,Dw.cols());
        if (Dw.rows()!=9 || Dw.cols()!=1)
            return;
        for (index_t c = 0; c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
            CHECK_EQUAL(1,WOStateOf(m,pt,0));
            const gsMatrix<real_t> D = m.eval3D_pstressDir(0,pt,WOZ(0),MaterialOutput::Generic);
            CHECK_EQUAL(9,D.rows());
            if (D.rows()!=9)
                return;
            const real_t dmin = WOSignedDiff(D.block(0,0,3,1),Dw.block(0,0,3,1));
            const real_t dmax = WOSignedDiff(D.block(3,0,3,1),Dw.block(3,0,3,1));
            const real_t nn   = (D.block(6,0,3,1)-Dw.block(6,0,3,1)).cwiseAbs().maxCoeff();
            gsInfo << "  dim " << dim << ", " << WOMatName(mat) << ", " << WOName(c) << ": d_min " << dmin << ", d_max " << dmax << ", N " << nn << "\n";
            CHECK(dmin <= 1e-10);
            CHECK(dmax <= 1e-10);
            CHECK(nn <= 1e-12);
        }
    }

    TEST(Frame_taut_matches_wrapped)
    {
        WOTautWrappedCase<2>(WOMatNHComp);
        WOTautWrappedCase<3>(WOMatNHComp);
        WOTautWrappedCase<2>(WOMatSvK);
        WOTautWrappedCase<3>(WOMatSvK);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Wrinkle direction of an affine wrinkled state against the analytic direction Q R e1 (mapped to the deformed plane
    /// in 3D). The deformed basis is rotated w.r.t. the reference basis, so a reference-basis direction would differ.
    template<short_t dim, bool linear>
    static void WOWrinkleDirCase(const WOMaterial mat, const bool exact, const char * name)
    {
        const WOStateSpec w = WOWrinkledSpec();
        WOFixture<dim,linear> fx(mat,WOSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        fx.split->options().setInt("ExactLocalSolve",0);
        if (exact)
            fx.split->options().setSwitch("Exact",true);
        const gsVector<real_t> pt = WOPoint();
        for (index_t c = (exact ? 1 : 0); c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
            CHECK_EQUAL(0,WOStateOf(m,pt,0));
            const gsMatrix<real_t> d = m.eval3D_wrinkleDir(0,pt,WOZ(0),MaterialOutput::Generic);
            CHECK_EQUAL(3,d.rows());
            CHECK_EQUAL(1,d.cols());
            if (d.rows()!=3 || d.cols()!=1)
                return;
            CHECK(d.allFinite());
            const real_t lenErr  = math::abs(d.norm()-1);
            const real_t exactErr = WOSignedDiff(d,fx.dExact);
            const real_t tangErr = dim==3 ? math::abs(d.col(0).dot(fx.Ndef.col(0))) : math::abs(d(2,0));

            const real_t theta = m.eval3D_theta(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const gsMatrix<real_t> basis = WOJacobian(fx.mp,pt);
            const gsMatrix<real_t> vref = (-math::sin(theta)*basis.col(0) + math::cos(theta)*basis.col(1)).normalized();
            const real_t refDiff = WOSignedDiff(d,vref);
            gsInfo << "  dim " << dim << ", " << name << ", " << WOName(c) << ": |d|-1 " << lenErr << ", vs analytic " << exactErr
                   << ", tangency " << tangErr << ", vs reference-basis direction " << refDiff << "\n";
            CHECK(lenErr <= 1e-14);
            CHECK(exactErr <= 1e-8);
            CHECK(tangErr <= 1e-12);
            CHECK(refDiff > 1e-2);
        }
    }

    TEST(WrinkleDir_analytic_affine)
    {
        WOWrinkleDirCase<2,false>(WOMatNHComp,false,"NH comp");
        WOWrinkleDirCase<2,false>(WOMatNHIncomp,false,"NH incomp");
        WOWrinkleDirCase<2,false>(WOMatNHIncomp,true,"NH incomp exact");
        WOWrinkleDirCase<2,true>(WOMatSvK,false,"SvK");
        WOWrinkleDirCase<3,false>(WOMatNHComp,false,"NH comp");
        WOWrinkleDirCase<3,false>(WOMatNHIncomp,false,"NH incomp");
        WOWrinkleDirCase<3,false>(WOMatNHIncomp,true,"NH incomp exact");
        WOWrinkleDirCase<3,true>(WOMatSvK,false,"SvK");
    }

    /// The wrinkle direction on the curved, rotated pair against the eigenvectors of the wrapped material's Cauchy stress.
    /// For laws isotropic w.r.t. the reference metric, S and the deformed metric are coaxial w.r.t. G, so the Cauchy stress
    /// F S F^T / J has the eigenvectors F N_i; the wrinkle covector n is the minor one and d = t^a g_a satisfies d . F N_min = 0.
    /// The oracle uses neither the wrinkle angle nor the strain split of the classes.
    TEST(WrinkleDir_Cauchy_eigenvectors)
    {
        gsMatrix<real_t> u, z;
        WOCurvedUZ(u,z);
        const WOMaterial mats[2] = { WOMatNHComp, WOMatNHIncomp };
        for (index_t mi = 0; mi!=2; ++mi)
        {
            WOFixture<3,false> fx(mats[mi],WOCurvedRef(),WOCurvedDef(true),0.1);
            for (index_t c = 0; c!=2; ++c)
            {
                const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
                index_t used = 0, differs = 0;
                for (index_t k = 0; k!=u.cols(); ++k)
                {
                    const gsVector<real_t> pt = u.col(k);
                    if (WOStateOf(m,pt,0)!=0)
                        continue;
                    ++used;
                    const gsMatrix<real_t> basis = WOJacobian(fx.mp_def,pt);
                    const gsMatrix<real_t> Sig = fx.wrapped->eval3D_CauchyStress(0,pt,WOZ(0),MaterialOutput::CauchyStressN);
                    CHECK_EQUAL(3,Sig.rows());
                    if (Sig.rows()!=3)
                        return;
                    const gsMatrix<real_t> sigma = WOCartesianTensor(Sig,basis);
                    gsMatrix<real_t>::SelfAdjEigenSolver es(sigma);
                    const gsMatrix<real_t> ev = es.eigenvalues();
                    const gsMatrix<real_t> vmin = es.eigenvectors().col(0), vmax = es.eigenvectors().col(2);
                    CHECK(ev(0,0) < -1e-6 * ev(2,0));

                    const gsMatrix<real_t> d = m.eval3D_wrinkleDir(0,pt,WOZ(0),MaterialOutput::Generic);
                    CHECK_EQUAL(3,d.rows());
                    if (d.rows()!=3)
                        return;
                    const gsVector<real_t,3> Ndef = WONormal(basis);
                    const real_t lenErr = math::abs(d.norm()-1);
                    const real_t majErr = 1 - math::abs(d.col(0).dot(vmax.col(0)));
                    const real_t minDot = math::abs(d.col(0).dot(vmin.col(0)));
                    const real_t tang   = math::abs(d.col(0).dot(Ndef));
                    const gsMatrix<real_t> g1 = basis.col(0).normalized(), g2 = basis.col(1).normalized();
                    const real_t d1 = WOSignedDiff(d,g1), d2 = WOSignedDiff(d,g2);
                    gsInfo << "  " << WOMatName(mats[mi]) << ", " << WOName(c) << ", point " << k << ": sigma_min/sigma_max " << ev(0,0)/ev(2,0)
                           << ", 1-|d.v_max| " << majErr << ", |d.v_min| " << minDot << ", |d|-1 " << lenErr << ", tangency " << tang
                           << ", vs g1 " << d1 << ", vs g2 " << d2 << "\n";
                    CHECK(majErr <= 1e-10);
                    CHECK(minDot <= 1e-8);
                    CHECK(lenErr <= 1e-14);
                    CHECK(tang <= 1e-12);
                    if (d1 > 1e-2 && d2 > 1e-2)
                        ++differs;
                }
                gsInfo << "  " << WOMatName(mats[mi]) << ", " << WOName(c) << ": " << used << " wrinkled points, " << differs
                       << " with a direction different from both normalised g1 and g2\n";
                CHECK(used >= 1);
                CHECK(differs >= 1);
            }
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static const index_t WONumOutputs = 9;

    /// Row counts of the nine outputs of WOEvalOutput
    static index_t WOOutputRows(const index_t idx)
    {
        const index_t rows[WONumOutputs] = { 1, 1, 3, 1, 1, 2, 9, 3, 2 };
        return rows[idx];
    }

    static const char * WOOutputName(const index_t idx)
    {
        const char * names[WONumOutputs] = { "theta", "gamma", "wrinkleDir", "wrinkleIntensity", "gammaHat", "pstress",
                                             "pstressDir", "CauchyStress", "CauchyPStress" };
        return names[idx];
    }

    /// One of: theta, gamma, wrinkleDir, wrinkleIntensity, gammaHat, pstress, pstressDir, CauchyStress, CauchyPStress
    static gsMatrix<real_t> WOEvalOutput(const gsMaterialMatrixBase<real_t> & m, const index_t idx, const gsMatrix<real_t> & u, const gsMatrix<real_t> & z)
    {
        switch (idx)
        {
            case 0: return m.eval3D_theta(0,u,z,MaterialOutput::Generic);
            case 1: return m.eval3D_gamma(0,u,z,MaterialOutput::Generic);
            case 2: return m.eval3D_wrinkleDir(0,u,z,MaterialOutput::Generic);
            case 3: return m.eval3D_wrinkleIntensity(0,u,z,MaterialOutput::Generic);
            case 4: return m.eval3D_gammaHat(0,u,z,MaterialOutput::Generic);
            case 5: return m.eval3D_pstress(0,u,z,MaterialOutput::Generic);
            case 6: return m.eval3D_pstressDir(0,u,z,MaterialOutput::Generic);
            case 7: return m.eval3D_CauchyStress(0,u,z,MaterialOutput::CauchyStressN);
            default: return m.eval3D_CauchyPStress(0,u,z,MaterialOutput::Generic);
        }
    }

    /// Every output of TFT and Split, evaluated by one multi-point call on the curved, rotated pair, agrees
    /// column by column: the direction up to sign (TFT keeps the raw root of the wrinkle angle, Split folds it into [0,pi))
    TEST(Outputs_TFT_equals_Split)
    {
        WOFixture<3,false> fx(WOMatNHComp,WOCurvedRef(),WOCurvedDef(true),0.1);
        gsMatrix<real_t> u, z;
        WOCurvedUZ(u,z);
        const index_t n = u.cols()*z.rows();

        for (index_t idx = 0; idx!=WONumOutputs; ++idx)
        {
            const gsMatrix<real_t> A = WOEvalOutput(*fx.tft,idx,u,z);
            const gsMatrix<real_t> B = WOEvalOutput(*fx.split,idx,u,z);
            CHECK_EQUAL(WOOutputRows(idx),A.rows());
            CHECK_EQUAL(WOOutputRows(idx),B.rows());
            CHECK_EQUAL(n,A.cols());
            CHECK_EQUAL(n,B.cols());
            if (A.rows()!=WOOutputRows(idx) || B.rows()!=WOOutputRows(idx) || A.cols()!=n || B.cols()!=n)
                return;

            if (idx==0)
            {
                real_t worst = 0;
                for (index_t col = 0; col!=n; ++col)
                    worst = math::max(worst,WOThetaDiff(A(0,col),B(0,col)));
                gsInfo << "  theta: largest difference modulo pi " << worst << "\n";
                CHECK(A.allFinite() && B.allFinite());
                CHECK(worst <= 1e-8);
            }
            else if (idx==2)
            {
                real_t worst = 0;
                for (index_t col = 0; col!=n; ++col)
                {
                    const gsMatrix<real_t> a = A.col(col), b = B.col(col);
                    worst = math::max(worst,WOSignedDiff(a,b));
                    gsInfo << "  wrinkleDir column " << col << ": TFT.Split = " << a.col(0).dot(b.col(0)) << "\n";
                }
                gsInfo << "  wrinkleDir: largest signed difference " << worst << "\n";
                CHECK(A.allFinite() && B.allFinite());
                CHECK(worst <= 1e-8);
            }
            else if (idx==6)
            {
                real_t worstDirs = 0, worstN = 0;
                for (index_t col = 0; col!=n; ++col)
                {
                    worstDirs = math::max(worstDirs,math::max(WOSignedDiff(A.block(0,col,3,1),B.block(0,col,3,1)),
                                                             WOSignedDiff(A.block(3,col,3,1),B.block(3,col,3,1))));
                    worstN = math::max(worstN,(A.block(6,col,3,1)-B.block(6,col,3,1)).cwiseAbs().maxCoeff());
                }
                gsInfo << "  pstressDir: eigenvectors up to sign " << worstDirs << ", N " << worstN << "\n";
                CHECK(A.allFinite() && B.allFinite());
                CHECK(worstDirs <= 1e-8);
                CHECK(worstN <= 1e-12);
            }
            else if (idx!=1)
            {
                WO_CHECK_REL(WOOutputName(idx),A,B,1e-8);
            }
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Wrinkling intensity of an affine wrinkled state against the analytic minor Green-Lagrange strain
    /// max(0,-E_min) = (1 - lam2^2)/2; the reference parametrisation and the rigid rotations drop out.
    /// \a eta > 0 sets the Split regularisation.
    template<short_t dim, bool linear>
    static void WOIntensityCase(const WOMaterial mat, const bool exact, const real_t eta, const char * name, const gsMatrix<real_t> & Jref)
    {
        const WOStateSpec w = WOWrinkledSpec();
        WOFixture<dim,linear> fx(mat,Jref,w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        fx.split->options().setInt("ExactLocalSolve",0);
        if (exact)
            fx.split->options().setSwitch("Exact",true);
        if (eta > 0)
            fx.split->options().setReal("Eta",eta);
        const gsVector<real_t> pt = WOPoint();
        const real_t expected = -WOEmin(w);
        for (index_t c = ((exact || eta > 0) ? 1 : 0); c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
            CHECK_EQUAL(0,WOStateOf(m,pt,0));
            const gsMatrix<real_t> I = m.eval3D_wrinkleIntensity(0,pt,WOZ(0),MaterialOutput::Generic);
            CHECK_EQUAL(1,I.rows());
            CHECK_EQUAL(1,I.cols());
            if (I.rows()!=1 || I.cols()!=1)
                return;
            const real_t err = math::abs(I(0,0)-expected) / expected;
            gsInfo << "  dim " << dim << ", " << name << ", " << WOName(c) << ": intensity " << I(0,0) << " vs " << expected << ", rel.err " << err << "\n";
            CHECK(expected > 0.2);
            CHECK(err <= 1e-12);
        }
    }

    TEST(WrinkleIntensity_analytic)
    {
        const gsMatrix<real_t> Jrefs[2] = { WOSkewedJref(), WOIdentity2() };
        for (index_t j = 0; j!=2; ++j)
        {
            WOIntensityCase<2,false>(WOMatNHComp,false,0,"NH comp",Jrefs[j]);
            WOIntensityCase<2,false>(WOMatNHIncomp,false,0,"NH incomp",Jrefs[j]);
            WOIntensityCase<2,true>(WOMatSvK,false,0,"SvK",Jrefs[j]);
            WOIntensityCase<2,false>(WOMatNHIncomp,true,0,"NH incomp exact",Jrefs[j]);
            WOIntensityCase<2,false>(WOMatNHComp,false,1e-3,"NH comp eta",Jrefs[j]);
            WOIntensityCase<3,false>(WOMatNHComp,false,0,"NH comp",Jrefs[j]);
            WOIntensityCase<3,false>(WOMatNHIncomp,false,0,"NH incomp",Jrefs[j]);
            WOIntensityCase<3,true>(WOMatSvK,false,0,"SvK",Jrefs[j]);
            WOIntensityCase<3,false>(WOMatNHIncomp,true,0,"NH incomp exact",Jrefs[j]);
            WOIntensityCase<3,false>(WOMatNHComp,false,1e-3,"NH comp eta",Jrefs[j]);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// gamma-hat = gamma n^T G^{-1} n against two oracles. (1) The strain increment dE = E' - E returned by the decorator
    /// is gamma n n^T, so the non-zero eigenvalue of the pencil (dE, G) is gamma n^T G^{-1} n; this uses neither theta nor
    /// gamma. (2) The definition itself with theta and gamma of the class and the inverse metric computed here.
    template<short_t dim>
    static void WOGammaHatCase(const WOMaterial mat, const bool exact)
    {
        const WOStateSpec w = WOWrinkledSpec();
        WOFixture<dim,false> fx(mat,WOSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        fx.split->options().setInt("ExactLocalSolve",0);
        if (exact)
            fx.split->options().setSwitch("Exact",true);
        const gsVector<real_t> pt = WOPoint();
        for (index_t c = (exact ? 1 : 0); c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
            CHECK_EQUAL(0,WOStateOf(m,pt,0));
            const gsMatrix<real_t> Ep = m.eval3D_strain(0,pt,WOZ(0));
            CHECK_EQUAL(3,Ep.rows());
            if (Ep.rows()!=3)
                return;
            const gsMatrix<real_t> dE = Ep - fx.EengDef();
            gsMatrix<real_t> dEt(2,2);
            dEt << dE(0,0), 0.5*dE(2,0), 0.5*dE(2,0), dE(1,0);
            gsMatrix<real_t>::GenSelfAdjEigenSolver ges(dEt,fx.Gref);
            const gsMatrix<real_t> ev = ges.eigenvalues();
            const index_t big = math::abs(ev(0,0)) > math::abs(ev(1,0)) ? 0 : 1;

            const real_t gh = m.eval3D_gammaHat(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t theta = m.eval3D_theta(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t gamma = m.eval3D_gamma(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            gsMatrix<real_t> n(2,1);
            n << math::cos(theta), math::sin(theta);
            const real_t gh2 = gamma * (n.transpose() * fx.Gref.inverse() * n)(0,0);

            const real_t err1 = math::abs(ev(big,0)-gh) / math::abs(gh);
            const real_t rank1 = math::abs(ev(1-big,0)) / math::abs(gh);
            const real_t err2 = math::abs(gh2-gh) / math::abs(gh);
            gsInfo << "  dim " << dim << ", " << WOMatName(mat) << (exact ? " exact" : "") << ", " << WOName(c) << ": gammaHat " << gh
                   << ", gamma " << gamma << ", strain-pencil oracle " << err1 << ", other eigenvalue " << rank1 << ", definition oracle " << err2 << "\n";
            CHECK(std::isfinite(gh));
            CHECK(gh > 0);
            CHECK(err1 <= 1e-10);
            CHECK(rank1 <= 1e-10);
            CHECK(err2 <= 1e-12);
        }
    }

    TEST(GammaHat_definition)
    {
        WOGammaHatCase<2>(WOMatNHComp,false);
        WOGammaHatCase<2>(WOMatNHIncomp,false);
        WOGammaHatCase<2>(WOMatNHIncomp,true);
        WOGammaHatCase<3>(WOMatNHComp,false);
        WOGammaHatCase<3>(WOMatNHIncomp,false);
        WOGammaHatCase<3>(WOMatNHIncomp,true);
    }

    /// The same Cartesian state under two reference parametrisations (identity and skewed): gamma-hat, the intensity and the
    /// wrinkle direction are properties of the surface state and agree, gamma itself is a parametrisation-dependent component
    template<short_t dim>
    static void WOParametrisationCase()
    {
        const WOStateSpec w = WOWrinkledSpec();
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        WOFixture<dim,false> f1(WOMatNHComp,WOIdentity2(),w.lam1,w.lam2,pi/6,0.4);
        WOFixture<dim,false> f2(WOMatNHComp,WOSkewedJref(),w.lam1,w.lam2,pi/6,0.4);
        const gsVector<real_t> pt = WOPoint();
        for (index_t c = 0; c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m1 = *f1.decorator(c), & m2 = *f2.decorator(c);
            CHECK_EQUAL(0,WOStateOf(m1,pt,0));
            CHECK_EQUAL(0,WOStateOf(m2,pt,0));
            const real_t gh1 = m1.eval3D_gammaHat(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t gh2 = m2.eval3D_gammaHat(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t g1  = m1.eval3D_gamma(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t g2  = m2.eval3D_gamma(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t i1  = m1.eval3D_wrinkleIntensity(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const real_t i2  = m2.eval3D_wrinkleIntensity(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
            const gsMatrix<real_t> d1 = m1.eval3D_wrinkleDir(0,pt,WOZ(0),MaterialOutput::Generic);
            const gsMatrix<real_t> d2 = m2.eval3D_wrinkleDir(0,pt,WOZ(0),MaterialOutput::Generic);
            const real_t ghErr = math::abs(gh1-gh2) / math::abs(gh1);
            const real_t iErr  = math::abs(i1-i2) / math::abs(i1);
            const real_t dErr  = WOSignedDiff(d1,d2);
            const real_t gDiff = math::abs(g1-g2) / math::abs(g1);
            gsInfo << "  dim " << dim << ", " << WOName(c) << ": gammaHat " << gh1 << " vs " << gh2 << " (rel.err " << ghErr << "), intensity rel.err " << iErr
                   << ", direction " << dErr << ", gamma " << g1 << " vs " << g2 << " (relative difference " << gDiff << ")\n";
            CHECK(std::isfinite(gh1) && std::isfinite(gh2));
            CHECK(ghErr <= 1e-8);
            CHECK(iErr <= 1e-12);
            CHECK(dErr <= 1e-8);
            CHECK(gDiff > 1e-2);
        }
    }

    TEST(GammaHat_parametrisation_invariant)
    {
        WOParametrisationCase<2>();
        WOParametrisationCase<3>();
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Outputs at one non-wrinkled point: the direction and the intensity are exactly zero, gamma-hat, theta and gamma are NaN,
    /// and the principal frame stays orthonormal
    template<short_t dim>
    static void WOCheckNonWrinkled(const char * label, const gsMaterialMatrixBase<real_t> & m, const gsMultiPatch<real_t> & mp,
                                   const gsVector<real_t> & pt, const index_t state)
    {
        CHECK_EQUAL(state,WOStateOf(m,pt,0));
        const gsMatrix<real_t> d = m.eval3D_wrinkleDir(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> I = m.eval3D_wrinkleIntensity(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> gh = m.eval3D_gammaHat(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> th = m.eval3D_theta(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> ga = m.eval3D_gamma(0,pt,WOZ(0),MaterialOutput::Generic);
        CHECK_EQUAL(3,d.rows());
        CHECK_EQUAL(1,I.rows());
        CHECK_EQUAL(1,gh.rows());
        CHECK_EQUAL(1,th.rows());
        CHECK_EQUAL(1,ga.rows());
        if (d.rows()!=3 || I.rows()!=1 || gh.rows()!=1 || th.rows()!=1 || ga.rows()!=1)
            return;
        gsInfo << "  dim " << dim << ", " << label << ": state " << state << ", |wrinkleDir| " << d.cwiseAbs().maxCoeff() << ", intensity " << I(0,0)
               << ", gammaHat " << gh(0,0) << ", theta " << th(0,0) << ", gamma " << ga(0,0) << "\n";
        CHECK(d.cwiseAbs().maxCoeff() == 0);
        CHECK(I(0,0) == 0);
        CHECK(std::isnan(gh(0,0)));
        CHECK(std::isnan(th(0,0)));
        CHECK(std::isnan(ga(0,0)));

        gsMatrix<real_t> F;
        real_t orthErr = 0, detErr = 0, nErr = 0;
        if (WOUnpackFrame(m.eval3D_pstressDir(0,pt,WOZ(0),MaterialOutput::Generic),WOJacobian(mp,pt),F,orthErr,detErr,nErr))
        {
            CHECK(orthErr <= 1e-12);
            CHECK(detErr <= 1e-12);
            CHECK(nErr <= 1e-12);
        }
    }

    /// Affine wrinkled point with the Split display threshold S2crit lowered so that the point is displayed taut: the outputs
    /// follow the criterion state and stay those of a wrinkled point
    template<short_t dim>
    static void WOS2critCase()
    {
        const WOStateSpec w = WOWrinkledSpec();
        WOFixture<dim,false> fx(WOMatNHIncomp,WOSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6,0.4);
        const gsVector<real_t> pt = WOPoint();
        CHECK_EQUAL(0,WOStateOf(*fx.split,pt,0));
        fx.split->options().setReal("S2crit",-1e30);
        const index_t shown = WOStateOf(*fx.split,pt,0);
        const gsMatrix<real_t> d = fx.split->eval3D_wrinkleDir(0,pt,WOZ(0),MaterialOutput::Generic);
        const real_t I = fx.split->eval3D_wrinkleIntensity(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
        const real_t gh = fx.split->eval3D_gammaHat(0,pt,WOZ(0),MaterialOutput::Generic)(0,0);
        gsInfo << "  dim " << dim << ", S2crit -1e30: displayed state " << shown << ", |wrinkleDir| " << d.norm() << ", intensity " << I << ", gammaHat " << gh << "\n";
        CHECK_EQUAL(1,shown);
        CHECK_EQUAL(3,d.rows());
        if (d.rows()!=3)
            return;
        CHECK(math::abs(d.norm()-1) <= 1e-14);
        CHECK(I > 0);
        CHECK(std::isfinite(gh));
    }

    /// Non-affine 2D pair with wrinkled and taut points
    static gsMultiPatch<real_t> WOMixedDef()
    {
        const WOStateSpec w = WOWrinkledSpec();
        gsMultiPatch<real_t> def = WOAffinePatch(WOJdef(WOSkewedJref(),w.lam1,w.lam2,static_cast<real_t>(EIGEN_PI)/6,0.0));
        def.patch(0).coefs()(3,1) += 0.8;
        return def;
    }

    TEST(Outputs_nonwrinkled_conventions)
    {
        const real_t eta[4] = { 0, 1e-3, 0, 1e-3 }, slack[4] = { 0, 0, 0.25, 0.25 };
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const gsVector<real_t> pt = WOPoint();

        // The taut point has a negative minor Green-Lagrange strain, so a missing state gate on the intensity would show
        const WOStateSpec taut = WOTautSpec(), sl = WOSlackSpec();
        gsInfo << "  taut stretches (" << taut.lam1 << "," << taut.lam2 << "): E_min = " << WOEmin(taut) << "\n";
        CHECK(WOEmin(taut) < 0);

        {
            WOFixture<2,false> fx(WOMatNHIncomp,WOSkewedJref(),taut.lam1,taut.lam2,pi/6,0.4);
            for (index_t c = 0; c!=2; ++c)
                WOCheckNonWrinkled<2>((std::string("taut ") + WOName(c)).c_str(),*fx.decorator(c),fx.mp,pt,1);
        }
        {
            WOFixture<3,false> fx(WOMatNHIncomp,WOSkewedJref(),taut.lam1,taut.lam2,pi/6,0.4);
            for (index_t c = 0; c!=2; ++c)
                WOCheckNonWrinkled<3>((std::string("taut ") + WOName(c)).c_str(),*fx.decorator(c),fx.mp,pt,1);
        }
        for (index_t v = 0; v!=4; ++v)
        {
            WOFixture<2,false> f2(WOMatNHComp,WOSkewedJref(),sl.lam1,sl.lam2,pi/6,0.4);
            WOFixture<3,false> f3(WOMatNHComp,WOSkewedJref(),sl.lam1,sl.lam2,pi/6,0.4);
            f2.split->options().setReal("Eta",eta[v]);
            f2.split->options().setReal("SlackMultiplier",slack[v]);
            f2.tft->options().setReal("SlackMultiplier",slack[v]);
            f3.split->options().setReal("Eta",eta[v]);
            f3.split->options().setReal("SlackMultiplier",slack[v]);
            f3.tft->options().setReal("SlackMultiplier",slack[v]);
            for (index_t c = 0; c!=2; ++c)
            {
                const std::string label = std::string("slack ") + WOName(c) + ", eta " + std::to_string(eta[v]) + ", multiplier " + std::to_string(slack[v]);
                WOCheckNonWrinkled<2>(label.c_str(),*f2.decorator(c),f2.mp,pt,-1);
                WOCheckNonWrinkled<3>(label.c_str(),*f3.decorator(c),f3.mp,pt,-1);
            }
        }

        WOS2critCase<2>();
        WOS2critCase<3>();

        // Mixed pair, multi-point: the wrinkle angle is finite exactly where the direction is a unit vector, the intensity is
        // positive and gamma-hat is finite
        WOFixture<2,false> mx(WOMatNHIncomp,WOAffinePatch(WOSkewedJref()),WOMixedDef(),0.01);
        gsMatrix<real_t> u, z;
        WOCurvedUZ(u,z);
        const index_t n = u.cols()*z.rows();
        for (index_t c = 0; c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *mx.decorator(c);
            const gsMatrix<real_t> th = m.eval3D_theta(0,u,z,MaterialOutput::Generic);
            const gsMatrix<real_t> d  = m.eval3D_wrinkleDir(0,u,z,MaterialOutput::Generic);
            const gsMatrix<real_t> I  = m.eval3D_wrinkleIntensity(0,u,z,MaterialOutput::Generic);
            const gsMatrix<real_t> gh = m.eval3D_gammaHat(0,u,z,MaterialOutput::Generic);
            const gsMatrix<real_t> ga = m.eval3D_gamma(0,u,z,MaterialOutput::Generic);
            CHECK_EQUAL(n,th.cols());
            CHECK_EQUAL(n,d.cols());
            CHECK_EQUAL(n,I.cols());
            CHECK_EQUAL(n,gh.cols());
            CHECK_EQUAL(n,ga.cols());
            CHECK_EQUAL(3,d.rows());
            CHECK_EQUAL(1,th.rows());
            CHECK_EQUAL(1,I.rows());
            CHECK_EQUAL(1,gh.rows());
            CHECK_EQUAL(1,ga.rows());
            if (th.cols()!=n || d.cols()!=n || I.cols()!=n || gh.cols()!=n || ga.cols()!=n
                || d.rows()!=3 || th.rows()!=1 || I.rows()!=1 || gh.rows()!=1 || ga.rows()!=1)
                return;
            index_t wrinkled = 0, other = 0;
            for (index_t col = 0; col!=n; ++col)
            {
                const bool wr = std::isfinite(th(0,col));
                CHECK_EQUAL(wr,std::isfinite(ga(0,col)));
                CHECK_EQUAL(wr,std::isfinite(gh(0,col)));
                if (wr)
                {
                    ++wrinkled;
                    CHECK(math::abs(d.col(col).norm()-1) <= 1e-14);
                    CHECK(I(0,col) > 0);
                }
                else
                {
                    ++other;
                    CHECK(d.col(col).cwiseAbs().maxCoeff() == 0);
                    CHECK(I(0,col) == 0);
                    CHECK(std::isnan(gh(0,col)));
                    CHECK(std::isnan(ga(0,col)));
                }
            }
            gsInfo << "  mixed pair, " << WOName(c) << ": " << wrinkled << " wrinkled and " << other << " other columns\n";
            CHECK(wrinkled >= 1);
            CHECK(other >= 1);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// TFT outputs derived from the modified stress S' = S + gamma C n1 at an affine wrinkled point, and the coincidence of
    /// the stress outputs with those of the wrapped material at an affine taut point. Oracles: generalized eigenvalues of
    /// the metric pencils built here from S' and the metrics, the wrapped material's single-point stress, and the Jacobian
    /// J_w = S_w / sigma_w read off the wrapped material's own stress and Cauchy stress.
    static void WOModifiedStressCase(const WOMaterial mat)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const WOStateSpec w = WOWrinkledSpec();
        const gsVector<real_t> pt = WOPoint();
        WOFixture<2,false> fx(mat,WOSkewedJref(),w.lam1,w.lam2,pi/6,0.4);
        const gsMaterialMatrixBase<real_t> & tft = *fx.tft;
        const gsMaterialMatrixBase<real_t> & wr  = *fx.wrapped;
        CHECK_EQUAL(0,WOStateOf(tft,pt,0));

        const gsMatrix<real_t> S  = tft.eval3D_stress(0,pt,WOZ(0),MaterialOutput::VectorN);
        const gsMatrix<real_t> P  = tft.eval3D_pstress(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> Cs = tft.eval3D_CauchyStress(0,pt,WOZ(0),MaterialOutput::CauchyStressN);
        const gsMatrix<real_t> CP = tft.eval3D_CauchyPStress(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> Sw = wr.eval3D_stress(0,pt,WOZ(0),MaterialOutput::VectorN);
        const gsMatrix<real_t> Cw = wr.eval3D_CauchyStress(0,pt,WOZ(0),MaterialOutput::CauchyStressN);
        const gsMatrix<real_t> Pw = wr.eval3D_pstress(0,pt,WOZ(0),MaterialOutput::Generic);
        CHECK_EQUAL(3,S.rows());
        CHECK_EQUAL(2,P.rows());
        CHECK_EQUAL(3,Cs.rows());
        CHECK_EQUAL(2,CP.rows());
        CHECK_EQUAL(3,Sw.rows());
        CHECK_EQUAL(3,Cw.rows());
        CHECK_EQUAL(2,Pw.rows());
        if (S.rows()!=3 || P.rows()!=2 || Cs.rows()!=3 || CP.rows()!=2 || Sw.rows()!=3 || Cw.rows()!=3 || Pw.rows()!=2)
            return;

        gsMatrix<real_t> Smat(2,2);
        Smat << S(0,0), S(2,0), S(2,0), S(1,0);
        gsMatrix<real_t>::GenSelfAdjEigenSolver gp(Smat,fx.Gref.inverse());
        const gsMatrix<real_t> lam = gp.eigenvalues();
        const real_t pErr = (P-lam).cwiseAbs().maxCoeff() / lam.cwiseAbs().maxCoeff();
        const real_t pMinor = math::abs(P(0,0)) / math::abs(P(1,0));
        const real_t pContrast = Pw(0,0) / math::abs(Pw(1,0));
        gsInfo << "  " << WOMatName(mat) << ": pstress vs pencil " << pErr << ", |P0|/|P1| " << pMinor
               << ", unmodified wrapped P0/|P1| " << pContrast << "\n";
        CHECK(pErr <= 1e-12);
        CHECK(pMinor <= 1e-8);
        CHECK(pContrast < -1e-3);

        const real_t Jw = Sw(0,0) / Cw(0,0);
        const real_t J0sq = fx.Gdef.determinant() / fx.Gref.determinant();
        const gsMatrix<real_t> Sj = S / Jw;
        const real_t cErr = WORelErr(Cs,Sj);
        gsInfo << "  " << WOMatName(mat) << ": J_w " << Jw << ", J_w^2 " << Jw*Jw << " vs det g/det G " << J0sq << ", Cauchy vs S'/J_w " << cErr << "\n";
        CHECK(std::isfinite(Jw));
        CHECK(cErr <= 1e-12);
        if (mat==WOMatNHComp)
            CHECK(math::abs(Jw*Jw - J0sq) > 1e-6);
        else
            CHECK(math::abs(Jw - 1) <= 1e-12);

        gsMatrix<real_t>::GenSelfAdjEigenSolver gc(Smat/Jw,fx.Gdef.inverse());
        const gsMatrix<real_t> mu = gc.eigenvalues();
        const real_t cpErr = (CP-mu).cwiseAbs().maxCoeff() / mu.cwiseAbs().maxCoeff();
        const real_t cpMinor = math::abs(CP(0,0)) / CP.cwiseAbs().maxCoeff();
        gsInfo << "  " << WOMatName(mat) << ": CauchyPStress vs pencil " << cpErr << ", minor/max " << cpMinor << "\n";
        CHECK(cpErr <= 1e-12);
        CHECK(cpMinor <= 1e-8);

        // Taut control: the stress-based outputs are those of the wrapped material (principal Cauchy values: deformed-metric pencil)
        const WOStateSpec ta = WOTautSpec();
        WOFixture<2,false> tx(mat,WOSkewedJref(),ta.lam1,ta.lam2,pi/6,0.4);
        CHECK_EQUAL(1,WOStateOf(*tx.tft,pt,0));
        const gsMatrix<real_t> tP  = tx.tft->eval3D_pstress(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> tC  = tx.tft->eval3D_CauchyStress(0,pt,WOZ(0),MaterialOutput::CauchyStressN);
        const gsMatrix<real_t> tCP = tx.tft->eval3D_CauchyPStress(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> wP  = tx.wrapped->eval3D_pstress(0,pt,WOZ(0),MaterialOutput::Generic);
        const gsMatrix<real_t> wC  = tx.wrapped->eval3D_CauchyStress(0,pt,WOZ(0),MaterialOutput::CauchyStressN);
        const gsMatrix<real_t> wCP = tx.wrapped->eval3D_CauchyPStress(0,pt,WOZ(0),MaterialOutput::Generic);
        CHECK_EQUAL(wP.rows(),tP.rows());
        CHECK_EQUAL(wC.rows(),tC.rows());
        CHECK_EQUAL(wCP.rows(),tCP.rows());
        if (wP.rows()!=tP.rows() || wC.rows()!=tC.rows() || wCP.rows()!=tCP.rows())
            return;
        WO_CHECK_REL("taut pstress",tP,wP,1e-12);
        WO_CHECK_REL("taut CauchyStress",tC,wC,1e-12);
        // The wrapped material's own principal Cauchy values are not the oracle: they are the principal values of the Cauchy
        // components taken in the reference basis. The eigenvalues of (1/J) g S^con w.r.t. the deformed metric g are built here
        // from the wrapped single-point stress and Cauchy stress.
        const gsMatrix<real_t> tS = tx.wrapped->eval3D_stress(0,pt,WOZ(0),MaterialOutput::VectorN);
        const gsMatrix<real_t> tCw = tx.wrapped->eval3D_CauchyStress(0,pt,WOZ(0),MaterialOutput::CauchyStressN);
        const real_t tJ = tS(0,0) / tCw(0,0);
        gsMatrix<real_t> tSmat(2,2);
        tSmat << tS(0,0), tS(2,0), tS(2,0), tS(1,0);
        gsMatrix<real_t>::GenSelfAdjEigenSolver tg(tSmat/tJ,tx.Gdef.inverse());
        const gsMatrix<real_t> tmu = tg.eigenvalues();
        WO_CHECK_REL("taut CauchyPStress vs deformed-metric pencil",tCP,tmu,1e-12);
        gsInfo << "  taut, wrapped material: CauchyPStress " << wCP.transpose() << ", pstress " << wP.transpose()
               << "; deformed-metric pencil " << tmu.transpose() << "\n";
    }

    TEST(TFT_modified_stress_outputs)
    {
        WOModifiedStressCase(WOMatNHComp);
        WOModifiedStressCase(WOMatNHIncomp);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Largest entry-wise difference of \a A and \a B relative to the largest finite entry of \a B, infinite if the sizes or the
    /// NaN patterns differ
    static real_t WONaNAwareErr(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
    {
        if (A.rows()!=B.rows() || A.cols()!=B.cols())
            return std::numeric_limits<real_t>::infinity();
        real_t scale = 0, err = 0;
        for (index_t j = 0; j!=B.cols(); ++j)
            for (index_t i = 0; i!=B.rows(); ++i)
                if (std::isfinite(B(i,j)))
                    scale = math::max(scale,math::abs(B(i,j)));
        for (index_t j = 0; j!=B.cols(); ++j)
            for (index_t i = 0; i!=B.rows(); ++i)
            {
                if (std::isnan(A(i,j)) != std::isnan(B(i,j)))
                    return std::numeric_limits<real_t>::infinity();
                if (!std::isnan(A(i,j)))
                    err = math::max(err,math::abs(A(i,j)-B(i,j)));
            }
        return scale > 0 ? err / scale : err;
    }

    /// Multi-point calls of every output of both decorators on the curved pair, 2 fibres per point: each column
    /// j*u.cols()+k is the single-point call at (u.col(k), z(j,k))
    TEST(Multipoint_layout_outputs)
    {
        WOFixture<3,false> fx(WOMatNHComp,WOCurvedRef(),WOCurvedDef(false),0.1);
        gsMatrix<real_t> u, z;
        WOCurvedUZ(u,z);
        const index_t n = u.cols()*z.rows();

        for (index_t c = 0; c!=2; ++c)
        {
            const gsMaterialMatrixBase<real_t> & m = *fx.decorator(c);
            for (index_t idx = 0; idx!=WONumOutputs; ++idx)
            {
                const gsMatrix<real_t> M = WOEvalOutput(m,idx,u,z);
                CHECK_EQUAL(WOOutputRows(idx),M.rows());
                CHECK_EQUAL(n,M.cols());
                if (M.rows()!=WOOutputRows(idx) || M.cols()!=n)
                    return;
                real_t worst = 0;
                for (index_t j = 0; j!=z.rows(); ++j)
                    for (index_t k = 0; k!=u.cols(); ++k)
                    {
                        const gsVector<real_t> uk = u.col(k);
                        const gsMatrix<real_t> single = WOEvalOutput(m,idx,uk,WOZ(z(j,k)));
                        CHECK_EQUAL(WOOutputRows(idx),single.rows());
                        if (single.rows()!=WOOutputRows(idx))
                            return;
                        const gsMatrix<real_t> col = M.col(j*u.cols()+k);
                        if (idx==0)
                        {
                            CHECK_EQUAL(std::isnan(single(0,0)),std::isnan(col(0,0)));
                            if (!std::isnan(single(0,0)) && !std::isnan(col(0,0)))
                                worst = math::max(worst,WOThetaDiff(col(0,0),single(0,0)));
                        }
                        else
                            worst = math::max(worst,WONaNAwareErr(col,single));
                    }
                gsInfo << "  " << WOName(c) << ", " << WOOutputName(idx) << ": largest column difference to the single-point call " << worst << "\n";
                CHECK(worst <= 1e-12);
            }

            // The stress depends on the fibre coordinate: a column that used the coordinate of another column would differ
            real_t zEffect = 0;
            for (index_t j = 0; j!=z.rows(); ++j)
                for (index_t k = 1; k!=u.cols(); ++k)
                {
                    const gsVector<real_t> uk = u.col(k);
                    const gsMatrix<real_t> Pk = m.eval3D_pstress(0,uk,WOZ(z(j,k)),MaterialOutput::Generic);
                    const gsMatrix<real_t> P0 = m.eval3D_pstress(0,uk,WOZ(z(j,0)),MaterialOutput::Generic);
                    const gsMatrix<real_t> dk = m.eval3D_wrinkleDir(0,uk,WOZ(z(j,k)),MaterialOutput::Generic);
                    const gsMatrix<real_t> d0 = m.eval3D_wrinkleDir(0,uk,WOZ(z(j,0)),MaterialOutput::Generic);
                    zEffect = math::max(zEffect,math::max((Pk-P0).norm()/P0.norm(),(dk-d0).norm()));
                }
            gsInfo << "  " << WOName(c) << ": change of pstress or wrinkleDir with the fibre coordinate of another column " << zEffect << "\n";
            CHECK(zEffect > 1e-10);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// gsShellStressFunction of one stress type 100..104 on the decorator \a dec: target dimension and values at two points,
    /// equal to the decorator's own multi-point output at z = 0
    static gsMatrix<real_t> WOShellStress(const WOFixture<3,false> & fx, gsMaterialMatrixBase<real_t> * dec, const index_t type,
                                          const gsMatrix<real_t> & u, bool & ok)
    {
        const index_t expectedDim[5] = { 1, 1, 3, 1, 1 };
        gsMaterialMatrixContainer<real_t> cont(1);
        cont.set(0,dec);
        gsShellStressFunction<real_t> f(fx.mp,fx.mp_def,cont,0,static_cast<stress_type::type>(type));
        CHECK_EQUAL(expectedDim[type-100],f.targetDim());
        gsMatrix<real_t> r;
        f.eval_into(u,r);
        CHECK_EQUAL(expectedDim[type-100],r.rows());
        CHECK_EQUAL(u.cols(),r.cols());
        ok = (r.rows()==expectedDim[type-100] && r.cols()==u.cols());
        return r;
    }

    TEST(ShellStressFunction_wrinkle_types)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const WOStateSpec w = WOWrinkledSpec(), ta = WOTautSpec();
        WOFixture<3,false> wrk(WOMatNHComp,WOSkewedJref(),w.lam1,w.lam2,pi/6,0.4);
        WOFixture<3,false> tau(WOMatNHIncomp,WOSkewedJref(),ta.lam1,ta.lam2,pi/6,0.4);
        gsMatrix<real_t> u(2,2), z0(1,2);
        u << 0.5, 0.3,
             0.5, 0.7;
        z0.setZero();

        for (index_t c = 0; c!=2; ++c)
        {
            WOFixture<3,false> * fixtures[2] = { &wrk, &tau };
            for (index_t fi = 0; fi!=2; ++fi)
            {
                WOFixture<3,false> & fx = *fixtures[fi];
                gsMaterialMatrixBase<real_t> * dec = fx.decorator(c);
                for (index_t type = 100; type!=105; ++type)
                {
                    bool ok = false;
                    const gsMatrix<real_t> r = WOShellStress(fx,dec,type,u,ok);
                    if (!ok)
                        return;
                    const gsMatrix<real_t> ref = WOEvalOutput(*dec,type-100,u,z0);
                    const real_t err = WONaNAwareErr(r,ref);
                    gsInfo << "  " << WOName(c) << ", " << (fi==0 ? "wrinkled" : "taut") << ", stress type " << type << ": difference to the decorator output " << err << "\n";
                    CHECK(err <= 1e-12);

                    if (fi==0 && type==102)
                    {
                        for (index_t k = 0; k!=u.cols(); ++k)
                        {
                            const gsMatrix<real_t> rk = r.col(k);
                            const real_t dErr = WOSignedDiff(rk,wrk.dExact);
                            gsInfo << "    point " << k << ": direction vs analytic " << dErr << "\n";
                            CHECK(dErr <= 1e-8);
                        }
                    }
                    if (fi==0 && type==103)
                    {
                        const real_t expected = -WOEmin(w);
                        const real_t iErr = (r.row(0).array() - expected).abs().maxCoeff() / expected;
                        gsInfo << "    intensity vs analytic " << iErr << "\n";
                        CHECK(iErr <= 1e-12);
                    }
                    if (fi==1)
                    {
                        if (type==102)
                            CHECK(r.cwiseAbs().maxCoeff() == 0);
                        else if (type==103)
                            CHECK(r.cwiseAbs().maxCoeff() == 0);
                        else
                            CHECK(r.array().isNaN().all());
                    }
                }
            }

            gsMaterialMatrixBase<real_t> * dec = wrk.decorator(c);
            gsMatrix<real_t> pt = u.col(0);
            gsMaterialMatrixEvalSingle<real_t,MaterialOutput::WrinkleDir> eDir(0,dec,&wrk.mp_def,WOZ(0));
            gsMaterialMatrixEvalSingle<real_t,MaterialOutput::WrinkleIntensity> eInt(0,dec,&wrk.mp_def,WOZ(0));
            gsMaterialMatrixEvalSingle<real_t,MaterialOutput::GammaHat> eGh(0,dec,&wrk.mp_def,WOZ(0));
            CHECK_EQUAL(3,eDir.targetDim());
            CHECK_EQUAL(1,eInt.targetDim());
            CHECK_EQUAL(1,eGh.targetDim());
            const gsMatrix<real_t> dd = eDir.eval(pt);
            CHECK_EQUAL(3,dd.rows());
            if (dd.rows()==3)
                CHECK(math::abs(dd.col(0).norm()-1) <= 1e-14);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Exact mode of Split at an admissible eta > 0 point: S_eff = (1-eta) S_w(E') + eta S_w(E) with E' = E + gamma n1 the
    /// returned elastic strain, S_w the wrapped material on its own at the prescribed metric Cmat = 2 E + g_ori. The residual
    /// n1.S_w(E') = n2.S_w(E') = 0 pins E'; the blend is the independent value. An Exact=false twin must differ.
    static void WOExactEtaCase(const WOMaterial mat, const real_t tolS, const bool twinDiffers)
    {
        const real_t pi = static_cast<real_t>(EIGEN_PI);
        const real_t eta = 1e-3;
        const WOStateSpec w = WOWrinkledSpec();
        const gsVector<real_t> pt = WOPoint();
        const gsMatrix<real_t> zmat = WOZ(0);
        WOFixture<2,false> fx(mat,WOSkewedJref(),w.lam1,w.lam2,pi/6,0.0);
        WOFixture<2,false> lin(mat,WOSkewedJref(),w.lam1,w.lam2,pi/6,0.0);
        fx.split->options().setSwitch("Exact",true);
        fx.split->options().setInt("ExactLocalSolve",0);
        fx.split->options().setReal("Eta",eta);
        lin.split->options().setReal("Eta",eta);
        fx.split->resetCounters();

        const gsMatrix<real_t> Ep = fx.split->eval3D_strain(0,pt,zmat);
        const real_t theta = fx.split->eval3D_theta(0,pt,zmat,MaterialOutput::Generic)(0,0);
        const real_t gamma = fx.split->eval3D_gamma(0,pt,zmat,MaterialOutput::Generic)(0,0);
        CHECK_EQUAL(3,Ep.rows());
        if (Ep.rows()!=3)
            return;
        CHECK_EQUAL(0,WOStateOf(*fx.split,pt,0));
        CHECK(gamma > 0);

        const gsVector<real_t,3> n1 = WON1(theta), n2 = WON2(theta);
        const gsMatrix<real_t> gn1 = gamma * n1;
        const gsMatrix<real_t> dE = Ep - fx.EengDef();
        const real_t dErr = (dE - gn1).cwiseAbs().maxCoeff() / gn1.cwiseAbs().maxCoeff();
        CHECK(dE.allFinite());
        CHECK(dErr <= 1e-12);

        const gsMatrix<real_t> Sw  = fx.wrapped->eval3D_vector_C(fx.CmatFromE(fx.EengDef()),0,pt,0,MaterialOutput::VectorN);
        const gsMatrix<real_t> Swp = fx.wrapped->eval3D_vector_C(fx.CmatFromE(Ep),0,pt,0,MaterialOutput::VectorN);
        CHECK(Sw.allFinite());
        CHECK(Swp.allFinite());
        const real_t scale = Sw.norm();
        const real_t r1 = math::abs(n1.dot(Swp.col(0))), r2 = math::abs(n2.dot(Swp.col(0)));

        const gsMatrix<real_t> Sind = (1-eta) * Swp + eta * Sw;
        const gsMatrix<real_t> Seff = fx.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const gsMatrix<real_t> Sstr = fx.split->eval3D_stress(0,pt,zmat,MaterialOutput::StressN);
        const gsMatrix<real_t> Slin = lin.split->eval3D_vector(0,pt,zmat,MaterialOutput::VectorN);
        const real_t etaVisible = (Sind-Swp).norm() / Swp.norm();
        const real_t twin = WORelErr(Slin,Sind);
        gsInfo << "  " << WOMatName(mat) << ": gamma " << gamma << ", theta " << theta << ", E'-E vs gamma n1 " << dErr
               << ", residuals " << r1/scale << " " << r2/scale << ", eta contribution " << etaVisible << ", linearised twin vs blend " << twin << "\n";
        CHECK(r1 <= 1e-8 * scale);
        CHECK(r2 <= 1e-8 * scale);
        WO_CHECK_REL("S_eff vs blend",Seff,Sind,tolS);
        WO_CHECK_REL("StressN vs blend",Sstr,Sind,tolS);
        CHECK(etaVisible > 1e-5);
        if (twinDiffers)
            CHECK(twin > 1e-4);

        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
        CHECK_EQUAL(0,fx.split->numExactFallbacks());
        CHECK_EQUAL(0,fx.split->numExactReclassified());
    }

    TEST(Exact_eta_independent_value)
    {
        WOExactEtaCase(WOMatNHComp,1e-8,false);
        WOExactEtaCase(WOMatNHIncomp,1e-10,true);
    }

    /// Exact mode in 3D away from the mid-surface. The wrapped material's own fibre metrics G(z), g(z) (public BaseDim
    /// accessors, physical coordinate z*t) give the prescribed metric Cz = g(z) and the strain E = (g(z)-G(z))/2 independent of
    /// the decorator; the returned E' = E + gamma n1 must annihilate n1.S_w and n2.S_w at Cz' = Cz + 2 gamma (c^2,s^2,cs), and
    /// S_eff (eta = 0) is S_w(Cz').
    TEST(Exact_3D_fibre_z)
    {
        const real_t thickness = 0.1;
        WOFixture<3,false> fx(WOMatNHComp,WOCurvedRef(),WOCurvedDef(false),thickness);
        WOFixture<3,false> lin(WOMatNHComp,WOCurvedRef(),WOCurvedDef(false),thickness);
        fx.split->options().setSwitch("Exact",true);
        fx.split->options().setInt("ExactLocalSolve",0);
        fx.split->options().setReal("Eta",0);
        lin.split->options().setReal("Eta",0);
        fx.split->resetCounters();

        const gsMaterialMatrixBaseDim<3,real_t> * bd = dynamic_cast<const gsMaterialMatrixBaseDim<3,real_t>*>(fx.wrapped.get());
        CHECK(bd != nullptr);
        if (bd==nullptr)
            return;

        gsMatrix<real_t> pts(2,2);
        pts << 0.2, 0.7,
               0.3, 0.2;
        const real_t zs[3] = { 0, -0.4, 0.4 };
        real_t maxExactVsLin = 0, maxMetricChange = 0;
        for (index_t k = 0; k!=2; ++k)
        {
            const gsVector<real_t> pt = pts.col(k);
            gsMatrix<real_t> Cz0;
            for (index_t zi = 0; zi!=3; ++zi)
            {
                const real_t z = zs[zi];
                // The fibre metrics are read from the wrapped material's per-point cache right after it is filled
                bd->_computePoints(0,pts.col(k));
                const gsMatrix<real_t> Gz = bd->_getGcov_ori(0,z*thickness);
                const gsMatrix<real_t> gz = bd->_getGcov_def(0,z*thickness);
                gsMatrix<real_t> Cz(3,1);
                Cz << gz(0,0), gz(1,1), gz(0,1);
                gsMatrix<real_t> E(3,1);
                E << 0.5*(gz(0,0)-Gz(0,0)), 0.5*(gz(1,1)-Gz(1,1)), gz(0,1)-Gz(0,1);
                if (zi==0)
                    Cz0 = Cz;
                else
                    maxMetricChange = math::max(maxMetricChange,(Cz-Cz0).norm() / Cz0.norm());

                const gsMatrix<real_t> Sw = fx.wrapped->eval3D_vector_C(Cz,0,pt,z,MaterialOutput::VectorN);
                const gsMatrix<real_t> Sgeo = fx.wrapped->eval3D_vector(0,pt,WOZ(z),MaterialOutput::VectorN);
                const real_t selfCheck = WORelErr(Sw,Sgeo);
                gsInfo << "  point " << k << ", z " << z << ": prescribed-metric vs geometry stress of the wrapped material " << selfCheck << "\n";
                CHECK(selfCheck <= 1e-8);

                CHECK_EQUAL(0,WOStateOf(*fx.split,pt,z));
                const gsMatrix<real_t> Ep = fx.split->eval3D_strain(0,pt,WOZ(z));
                const real_t theta = fx.split->eval3D_theta(0,pt,WOZ(z),MaterialOutput::Generic)(0,0);
                const real_t gamma = fx.split->eval3D_gamma(0,pt,WOZ(z),MaterialOutput::Generic)(0,0);
                CHECK_EQUAL(3,Ep.rows());
                if (Ep.rows()!=3)
                    return;
                CHECK(gamma > 0);
                const gsVector<real_t,3> n1 = WON1(theta), n2 = WON2(theta);
                const gsMatrix<real_t> gn1 = gamma * n1;
                const real_t dErr = (Ep - E - gn1).cwiseAbs().maxCoeff() / gn1.cwiseAbs().maxCoeff();

                const real_t c = math::cos(theta), s = math::sin(theta);
                gsMatrix<real_t> Czp(3,1);
                Czp << Cz(0,0) + 2*gamma*c*c, Cz(1,0) + 2*gamma*s*s, Cz(2,0) + 2*gamma*c*s;
                const gsMatrix<real_t> Swp = fx.wrapped->eval3D_vector_C(Czp,0,pt,z,MaterialOutput::VectorN);
                const real_t r1 = math::abs(n1.dot(Swp.col(0))) / Sw.norm(), r2 = math::abs(n2.dot(Swp.col(0))) / Sw.norm();
                const gsMatrix<real_t> Seff = fx.split->eval3D_vector(0,pt,WOZ(z),MaterialOutput::VectorN);
                const real_t sErr = WORelErr(Seff,Swp);
                const gsMatrix<real_t> Slin = lin.split->eval3D_vector(0,pt,WOZ(z),MaterialOutput::VectorN);
                maxExactVsLin = math::max(maxExactVsLin,WORelErr(Slin,Seff));
                gsInfo << "    gamma " << gamma << ", E'-E vs gamma n1 " << dErr << ", residuals " << r1 << " " << r2 << ", S_eff vs S_w(Cz') " << sErr << "\n";
                CHECK(dErr <= 1e-10);
                CHECK(r1 <= 1e-8);
                CHECK(r2 <= 1e-8);
                CHECK(sErr <= 1e-8);
            }
        }
        gsInfo << "  fibre metric change with z " << maxMetricChange << ", exact vs linearised stress " << maxExactVsLin << "\n";
        CHECK(maxMetricChange > 1e-6);
        CHECK(maxExactVsLin > 1e-6);

        CHECK_EQUAL(0,fx.split->numThetaFallbacks());
        CHECK_EQUAL(0,fx.split->numTFTFailures());
        CHECK_EQUAL(0,fx.split->numExactFallbacks());
        CHECK_EQUAL(0,fx.split->numExactReclassified());
    }

}
