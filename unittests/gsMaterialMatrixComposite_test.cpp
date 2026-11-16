/** @file gsMaterialMatrixComposite_test.cpp

    @brief Oracle checks of the laminate matrices A, B, D and the stress resultants N, M of
    gsMaterialMatrixComposite.

    The laminate matrices are compared with a tensor-product oracle evaluated in this file: for a
    ply with fibre angle phi, the plane-stress tensor Q_{ijkl} given in the fibre frame is
    contracted with the fibre directions and the contravariant surface basis,
    C^{abcd} = sum Q_{ijkl} (f_i.a^a)(f_j.a^b)(f_k.a^c)(f_l.a^d), and integrated through the thickness
    (A = sum t C, B = sum t z C, D = sum (t z^2 + t^3/12) C). A single isotropic ply is
    also compared with gsMaterialMatrixLinear. The checks are run on an axis-aligned
    unit square and on a sheared, stretched and (in 3D) lifted patch on which the covariant basis is
    neither orthogonal nor of unit length. The stacking direction (ply 0 at the most negative z
    along the unit normal) is checked against the shell kinematics of the base class, with
    gsMaterialMatrixLinear stresses at physical heights as reference.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Hugo Verhelst
*/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

using namespace gismo;

SUITE(gsMaterialMatrixComposite_test)                 // The suite should have the same name as the file
{
    namespace {

    const real_t compPi = math::atan(1)*4;

    // Relative tolerance of the matrix and vector identities. Every compared quantity is a short
    // sum (at most a few hundred flops) of products of entries of size |C| t^k, so its round-off is
    // of order 1e-14 relative in double. 1e-10 leaves about four orders of margin over that, and is
    // still many orders below what a modelling error produces (a wrong fibre angle, metric, sign or
    // stacking order changes a result at the percent level or more). For float, 1e-10 lies below
    // round-off, so 1e5 eps (a fixed margin over machine precision) takes over; for double
    // 1e5 eps = 2.2e-11 < 1e-10 and the fixed value applies.
    real_t identityTol()
    {
        return std::max<real_t>(real_t(1e-10), real_t(1e5)*math::limits::epsilon());
    }

    // Unit square of degree 2 (3x3 control points, index 4 is the centre one), embedded in R^3 when
    // d == 3. With skewed == true the first coordinate becomes x <- 1.5 x + 0.4 y (control points are
    // transformed, the map is affine-invariant in them), and in 3D the centre control point is lifted
    // by 0.3, so that a1 = (1.5, 0, .) is longer than one, the angle between a1 and a2 is about 68 degrees
    // (cos about 0.37) and the surface is curved. The ply functions have domain dimension d because the
    // material evaluates them at the physical points. The material keeps non-owned pointers to the ply
    // functions, the reference patch and the deformed patch, so all of them are members, and the
    // containers are reserved before any address is taken.
    template <short_t d>
    struct Laminate
    {
        gsMultiPatch<> ref, def;
        std::vector<gsConstantFunction<> > thicks, Gfuns, phis;
        std::vector<gsFunctionSet<>*> Ts, Gs, Phis;
        real_t tTotal;
        memory::unique_ptr<gsMaterialMatrixComposite<d,real_t> > mm;

        // Plies of equal thickness tply with the fibre angles and the 3x3 ply matrices, both listed from
        // ply 0 (most negative z along the normal) upwards.
        Laminate(const std::vector<real_t> & angles, real_t tply, const std::vector<gsMatrix<> > & G3x3,
                 bool skewed)
        :
        tTotal(tply*angles.size())
        {
            ref.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
            ref.degreeElevate(1);
            if (d == 3)
                ref.embed(3);

            if (skewed)
            {
                gsMatrix<> & c = ref.patch(0).coefs();
                c.col(0) = 1.5*c.col(0) + 0.4*c.col(1);
                if (d == 3)
                    c(4,2) += 0.3;
            }
            def.addPatch( ref.patch(0).clone() );

            const size_t n = angles.size();
            thicks.reserve(n);
            Gfuns.reserve(n);
            phis.reserve(n);
            for (size_t i = 0; i != n; ++i)
            {
                gsMatrix<> Gmat = G3x3[i];
                Gmat.resize(Gmat.rows()*Gmat.cols(),1);
                thicks.push_back(gsConstantFunction<>(tply,d));
                Gfuns.push_back(gsConstantFunction<>(Gmat,d));
                phis.push_back(gsConstantFunction<>(angles[i],d));
            }
            for (size_t i = 0; i != n; ++i)
            {
                Ts.push_back(&thicks[i]);
                Gs.push_back(&Gfuns[i]);
                Phis.push_back(&phis[i]);
            }
            mm = memory::make_unique(new gsMaterialMatrixComposite<d,real_t>(ref,Ts,Gs,Phis));
            mm->setDeformed(&def);
        }

        // The same ply matrix for every ply.
        Laminate(const std::vector<real_t> & angles, real_t tply, const gsMatrix<> & G3x3, bool skewed)
        :
        Laminate(angles,tply,std::vector<gsMatrix<> >(angles.size(),G3x3),skewed)
        { }

    private:
        Laminate(const Laminate &);
        Laminate & operator=(const Laminate &);
    };

    // Orthotropic ply (Exx, Eyy, Gxy, nuxy, nuyx). Exx != Eyy, so a rotation of the fibre changes the
    // ply matrix. nuyx*Exx = 0.125*200 = nuxy*Eyy = 0.25*100 = 25 holds exactly in floating point, as
    // gsCompositeMatrix requires of a symmetric compliance. Gxy = 40 is a generic positive shear
    // modulus, not tied to Exx and Eyy by any isotropy relation, so the shear entry G22 enters the
    // rotation independently.
    gsMatrix<> orthotropicPly()
    {
        return gsCompositeMatrix<real_t>(200.0,100.0,40.0,0.25,0.125);
    }

    // Isotropic ply (E, nu): the plane-stress matrix of gsMaterialMatrixLinear.
    gsMatrix<> isotropicPly(real_t E, real_t nu)
    {
        return gsCompositeMatrix<real_t>(E,E,E/(2*(1+nu)),nu,nu);
    }

    gsMatrix<> pointsInterior()
    {
        gsMatrix<> u(2,3);
        u << 0.2, 0.5, 0.8,
             0.3, 0.6, 0.4;
        return u;
    }

    gsMatrix<> zeros(const gsMatrix<> & u)
    {
        return gsMatrix<>::Zero(1,u.cols());
    }

    // Covariant basis (columns a_1, a_2; d x 2) of patch 0 at the points u.
    std::vector<gsMatrix<> > covariantBases(const gsMultiPatch<> & g, const gsMatrix<> & u)
    {
        gsMapData<real_t> map;
        map.flags = NEED_JACOBIAN;
        map.points = u;
        static_cast<const gsFunction<real_t>&>(g.patch(0)).computeMap(map);
        std::vector<gsMatrix<> > bases(u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            bases[k] = map.jacobian(k);
        return bases;
    }

    // Covariant first and second fundamental forms of patch 0 at the points u (3D surfaces).
    void fundamentalForms(const gsMultiPatch<> & g, const gsMatrix<> & u,
                          std::vector<gsMatrix<> > & Acov, std::vector<gsMatrix<> > & Bcov)
    {
        gsMapData<real_t> map;
        map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
        map.points = u;
        static_cast<const gsFunction<real_t>&>(g.patch(0)).computeMap(map);
        gsMatrix<> normals = map.normals;
        normals.colwise().normalize();

        Acov.resize(u.cols());
        Bcov.resize(u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            gsMatrix<> acov = map.jacobian(k);
            Acov[k] = acov.transpose() * acov;
            gsMatrix<> deriv2 = map.deriv2(k);
            deriv2.resize(3,3);
            gsVector<> n = normals.col(k);
            Bcov[k].resize(2,2);
            Bcov[k](0,0) = deriv2.row(0).dot(n);
            Bcov[k](1,1) = deriv2.row(1).dot(n);
            Bcov[k](0,1) = Bcov[k](1,0) = deriv2.row(2).dot(n);
        }
    }

    // Engineering Voigt order [E00, E11, E01+E10], as folded by the material.
    gsVector<real_t,3> fold(const gsMatrix<> & E)
    {
        gsVector<real_t,3> v;
        v << E(0,0), E(1,1), E(0,1)+E(1,0);
        return v;
    }

    // Green-Lagrange membrane strain eps = fold(1/2 (A_def - A_ref)) and curvature change
    // kap = fold(B_ref - B_def) per point.
    void strains(const gsMultiPatch<> & ref, const gsMultiPatch<> & def, const gsMatrix<> & u,
                 std::vector<gsVector<real_t,3> > & eps, std::vector<gsVector<real_t,3> > & kap)
    {
        std::vector<gsMatrix<> > Ar, Br, Ad, Bd;
        fundamentalForms(ref,u,Ar,Br);
        fundamentalForms(def,u,Ad,Bd);
        eps.resize(u.cols());
        kap.resize(u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            eps[k] = fold(0.5*(Ad[k]-Ar[k]));
            kap[k] = fold(Br[k]-Bd[k]);
        }
    }

    // Deformed state with non-zero membrane strain and curvature change at every evaluation point:
    // a 5% stretch of the first coordinate, a 3% shear (second coordinate gains 3% of the first) and a
    // non-uniform out-of-plane lift of two control points. Both N and M are then non-trivial. The
    // isotropic-ply tests only need some deformed state, because their matrices do not depend on it.
    void stretchedBent(const gsMultiPatch<> & ref, gsMultiPatch<> & def)
    {
        gsMatrix<> & c = def.patch(0).coefs();
        c = ref.patch(0).coefs();
        c.col(0) *= 1.05;
        c.col(1) += 0.03*ref.patch(0).coefs().col(0);
        c(4,2) += 0.5;
        c(8,2) += 0.25;
    }

    // Smallest values over the points of the quantities that make a skewed fixture meaningful.
    struct GeometryGuards
    {
        real_t minCosAngle;      // |a1.a2|/(|a1||a2|): how far from orthogonal the basis is
        real_t minLengthDefect;  // | |a1| - 1 |: how far from unit length a1 is
        real_t minCurvature;     // Frobenius norm of the covariant second fundamental form (3D only)
    };

    GeometryGuards geometryGuards(const gsMultiPatch<> & ref, const gsMatrix<> & u, bool curved)
    {
        GeometryGuards g;
        g.minCosAngle = g.minLengthDefect = g.minCurvature = std::numeric_limits<real_t>::max();
        const std::vector<gsMatrix<> > bases = covariantBases(ref,u);
        for (size_t k = 0; k != bases.size(); ++k)
        {
            const gsVector<> a1 = bases[k].col(0), a2 = bases[k].col(1);
            g.minCosAngle     = std::min(g.minCosAngle,     math::abs(a1.dot(a2))/(a1.norm()*a2.norm()));
            g.minLengthDefect = std::min(g.minLengthDefect, math::abs(a1.norm()-1));
        }
        if (curved)
        {
            std::vector<gsMatrix<> > Acov, Bcov;
            fundamentalForms(ref,u,Acov,Bcov);
            for (size_t k = 0; k != Bcov.size(); ++k)
                g.minCurvature = std::min(g.minCurvature, Bcov[k].norm());
        }
        return g;
    }

    // Laminate matrices A, B, D (3x3 Voigt [11,22,12], contravariant components) at one point from
    // the covariant basis acov (d x 2). Every ply of the laminate is oriented at its own angle, the
    // tensor Q_{ijkl} (fibre frame) is contracted with the fibre directions f_i and the contravariant
    // basis a^a (a^a . a_b = delta^a_b), and the plies are stacked in +n from the most negative z.
    void tensorOracle(const gsMatrix<> & acov, const gsMatrix<> & G,
                      const std::vector<real_t> & thickness, const std::vector<real_t> & angles,
                      gsMatrix<> & A, gsMatrix<> & B, gsMatrix<> & D)
    {
        // Contravariant basis a^a: columns of a_cov (a_cov^T a_cov)^{-1}.
        const gsMatrix<> acon = acov * (acov.transpose()*acov).inverse();
        // Orthonormal tangent frame; a^2 is orthogonal to a_1 so e1 is a_1/|a_1| and e2 is a^2/|a^2|.
        const gsVector<> e1 = acov.col(0).normalized();
        const gsVector<> e2 = acon.col(1).normalized();

        // Q_{ijkl}: Q1111 = G00, Q2222 = G11, Q1122 = Q2211 = G01, Q1212 = Q1221 = Q2112 = Q2121 = G22.
        real_t Q[2][2][2][2];
        for (index_t i = 0; i != 2; ++i)
            for (index_t j = 0; j != 2; ++j)
                for (index_t k = 0; k != 2; ++k)
                    for (index_t l = 0; l != 2; ++l)
                        Q[i][j][k][l] = 0.0;
        Q[0][0][0][0] = G(0,0);
        Q[1][1][1][1] = G(1,1);
        Q[0][0][1][1] = Q[1][1][0][0] = G(0,1);
        Q[0][1][0][1] = Q[0][1][1][0] = Q[1][0][0][1] = Q[1][0][1][0] = G(2,2);

        // Voigt index -> tensor index pair.
        const index_t vi[3] = {0,1,0};
        const index_t vj[3] = {0,1,1};

        real_t tTot = 0.0;
        for (size_t p = 0; p != thickness.size(); ++p)
            tTot += thickness[p];

        A = gsMatrix<>::Zero(3,3);
        B = gsMatrix<>::Zero(3,3);
        D = gsMatrix<>::Zero(3,3);
        real_t tBelow = 0.0;
        for (size_t p = 0; p != thickness.size(); ++p)
        {
            const real_t t = thickness[p];
            const real_t c = math::cos(angles[p]), s = math::sin(angles[p]);
            // Fibre frame: f_1 along the fibre (phi from e1 towards e2), f_2 perpendicular to it.
            gsVector<> f[2];
            f[0] =  c*e1 + s*e2;
            f[1] = -s*e1 + c*e2;
            // P(i,a) = f_i . a^a
            real_t P[2][2];
            for (index_t i = 0; i != 2; ++i)
                for (index_t a = 0; a != 2; ++a)
                    P[i][a] = f[i].dot(acon.col(a));

            gsMatrix<> C(3,3);
            for (index_t I = 0; I != 3; ++I)
                for (index_t J = 0; J != 3; ++J)
                {
                    const index_t a = vi[I], b = vj[I], g = vi[J], h = vj[J];
                    real_t val = 0.0;
                    for (index_t i = 0; i != 2; ++i)
                        for (index_t j = 0; j != 2; ++j)
                            for (index_t k = 0; k != 2; ++k)
                                for (index_t l = 0; l != 2; ++l)
                                    val += Q[i][j][k][l]*P[i][a]*P[j][b]*P[k][g]*P[l][h];
                    C(I,J) = val;
                }

            // Height of the ply mid-surface above the laminate mid-surface.
            const real_t z = (tBelow + t/2) - tTot/2;
            A += t*C;
            B += t*z*C;
            D += (t*z*z + t*t*t/12)*C;
            tBelow += t;
        }
    }

    // Oracle matrices at all points: column k of the outputs holds the 3x3 matrix of point k as in
    // the material (9 x npts, column-major 3x3).
    void tensorOracleAtPoints(const gsMultiPatch<> & ref, const gsMatrix<> & u, const gsMatrix<> & G,
                              real_t tply, const std::vector<real_t> & angles,
                              gsMatrix<> & A, gsMatrix<> & B, gsMatrix<> & D)
    {
        const std::vector<gsMatrix<> > bases = covariantBases(ref,u);
        const std::vector<real_t> thickness(angles.size(),tply);
        A.resize(9,u.cols());
        B.resize(9,u.cols());
        D.resize(9,u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            gsMatrix<> Ak, Bk, Dk;
            tensorOracle(bases[k],G,thickness,angles,Ak,Bk,Dk);
            A.reshapeCol(k,3,3) = Ak;
            B.reshapeCol(k,3,3) = Bk;
            D.reshapeCol(k,3,3) = Dk;
        }
    }

    // Constant isotropic plane-stress matrix of gsMaterialMatrixLinear (9 x npts) on the reference patch.
    gsMatrix<> linearMatrix(const gsMultiPatch<> & ref, const gsMultiPatch<> & def, const gsMatrix<> & u,
                            real_t E, real_t nu, real_t t)
    {
        gsConstantFunction<> tfun(t,3), Efun(E,3), nufun(nu,3);
        gsMaterialMatrixLinear<3,real_t> lin(ref,tfun,Efun,nufun);
        lin.setDeformed(&def);
        return lin.eval3D_matrix(0,u,zeros(u));
    }

    // Matrices A, B, D of the material at the points u (each 9 x npts).
    template <short_t d>
    void compositeABD(const Laminate<d> & L, const gsMatrix<> & u,
                      gsMatrix<> & A, gsMatrix<> & B, gsMatrix<> & D)
    {
        A = L.mm->eval3D_matrix(0,u,zeros(u),MaterialOutput::MatrixA);
        B = L.mm->eval3D_matrix(0,u,zeros(u),MaterialOutput::MatrixB);
        D = L.mm->eval3D_matrix(0,u,zeros(u),MaterialOutput::MatrixD);
    }

    // Frobenius norm of the 3x3 matrix of point k.
    real_t normAt(const gsMatrix<> & M, index_t k)
    {
        return gsMatrix<>(M.reshapeCol(k,3,3)).norm();
    }

    // Frobenius norm of the difference of the 3x3 matrices of point k.
    real_t diffAt(const gsMatrix<> & M, const gsMatrix<> & Mref, index_t k)
    {
        return gsMatrix<>(M.reshapeCol(k,3,3) - Mref.reshapeCol(k,3,3)).norm();
    }

    // Frobenius norm of M - M^T at point k.
    real_t asymmetryAt(const gsMatrix<> & M, index_t k)
    {
        gsMatrix<> Mk = M.reshapeCol(k,3,3);
        return (Mk - gsMatrix<>(Mk.transpose())).norm();
    }

    // Three different angles, one of them negative (exercises the sign of the angle), and not mirror
    // symmetric, so that B != 0.
    std::vector<real_t> mixedAngles()
    {
        std::vector<real_t> a(3);
        a[0] = 0.0; a[1] = compPi/4; a[2] = -compPi/3;
        return a;
    }

    // The smallest unsymmetric laminate, so that B != 0.
    std::vector<real_t> twoPlyAngles()
    {
        std::vector<real_t> a(2);
        a[0] = 0.0; a[1] = compPi/4;
        return a;
    }

    // Mirror symmetric, so B = 0 up to round-off. This exercises the |A| t_tot scaling of the B check,
    // where a relative error |dB|/|B| would be undefined.
    std::vector<real_t> symmetricAngles()
    {
        std::vector<real_t> a(4);
        a[0] = 0.0; a[1] = compPi/4; a[2] = compPi/4; a[3] = 0.0;
        return a;
    }

    // Asserts the geometry of a skewed fixture before any value is compared: a vacuous skew (orthogonal
    // or unit-length basis, flat surface) would make the oracle comparison pass for the wrong reason.
    // Thresholds:
    //  - |cos(angle)| >= 0.1: the angle between a1 and a2 is at least about 6 degrees away from a
    //    right angle, so the off-diagonal metric term is at least 10% of sqrt(A11 A22) and code
    //    assuming an orthogonal basis would be off by far more than identityTol().
    //  - | |a1| - 1 | >= 0.1: a1 is at least 10% off unit length, so confusing the covariant basis with
    //    its normalised version would show at that order.
    //  - |b| >= 1e-3 (3D): the reference surface is really curved, so the normal and the second
    //    fundamental form vary over the patch and code assuming a fixed normal is exposed. A flat patch
    //    has |b| at round-off level (about 1e-16), so 1e-3 only confirms that the lift was applied.
    void checkSkewGuards(const GeometryGuards & g, bool curved)
    {
        CHECK(g.minCosAngle     >= 0.1);
        CHECK(g.minLengthDefect >= 0.1);
        if (curved)
            CHECK(g.minCurvature >= 1e-3);
    }

    } // namespace

    // One isotropic ply with angle 0 is rotation invariant, so on an orthogonal unit-length patch
    // A = t C_lin, D = t^3/12 C_lin and B = 0, with C_lin the plane-stress matrix of gsMaterialMatrixLinear
    // at z = 0 (an independent implementation of the same continuum tensor).
    TEST(Orthogonal_IsotropicPly_MatchesLinear)
    {
        // Set-up: one isotropic ply on the flat unit square (metric = identity). E only scales the
        // result; nu = 0.3 is generic (non-zero, so the lambda A^{ij}A^{kl} term of the plane-stress
        // tensor is exercised, and away from the incompressible limit 0.5); t = 0.01 is a thin shell,
        // and the identities checked are exact for any t. The deformation is arbitrary because neither
        // matrix depends on it.
        const real_t E = 100.0, nu = 0.3, t = 0.01, tol = identityTol();
        Laminate<3> L(std::vector<real_t>(1,0.0),t,isotropicPly(E,nu),false);
        const gsMatrix<> u = pointsInterior();
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);

        // Act: the laminate matrices of the material and, as oracle, the matrix of
        // gsMaterialMatrixLinear for the same E, nu and t.
        gsMatrix<> A, B, D;
        compositeABD(L,u,A,B,D);
        const gsMatrix<> Clin = linearMatrix(L.ref,L.def,u,E,nu,t);

        // Check: at every point A and D equal the Linear integrals and B vanishes, each to identityTol().
        for (index_t k = 0; k != u.cols(); ++k)
        {
            // A_ref = t C_lin, D_ref = t^3/12 C_lin (the integrals of a through-thickness constant).
            const gsMatrix<> Aref = t*Clin.reshapeCol(k,3,3);
            const gsMatrix<> Dref = t*t*t/12*Clin.reshapeCol(k,3,3);
            CHECK_CLOSE(0.0, gsMatrix<>(A.reshapeCol(k,3,3)-Aref).norm()/Aref.norm(), tol);
            CHECK_CLOSE(0.0, gsMatrix<>(D.reshapeCol(k,3,3)-Dref).norm()/Dref.norm(), tol);
            // B vanishes for a single ply centred on the mid-surface; scaled by |A| t_tot (B = 0 has no relative scale).
            CHECK_CLOSE(0.0, normAt(B,k)/(normAt(A,k)*t), tol);
        }
    }

    // Same identity as the control test on a sheared, stretched, curved patch where the covariant
    // basis is not orthonormal. The isotropic tensor is invariant under rotation of the ply, so the
    // gsMaterialMatrixLinear matrix is a valid oracle there too (it uses the true inverse metric).
    TEST(Skewed_IsotropicPly_MatchesLinear)
    {
        // Set-up: the isotropic ply of the control test (same E, nu, t and the same reasons) on the
        // skewed, lifted patch. The guards confirm that the skew and the curvature are real, so the
        // comparison with the Linear matrix is not the orthogonal case in disguise.
        const real_t E = 100.0, nu = 0.3, t = 0.01, tol = identityTol();
        Laminate<3> L(std::vector<real_t>(1,0.0),t,isotropicPly(E,nu),true);
        const gsMatrix<> u = pointsInterior();
        checkSkewGuards(geometryGuards(L.ref,u,true),true);
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);

        // Act: the laminate matrices of the material and the matrix of gsMaterialMatrixLinear, which
        // forms its tensor from the true inverse metric of the same patch.
        gsMatrix<> A, B, D;
        compositeABD(L,u,A,B,D);
        const gsMatrix<> Clin = linearMatrix(L.ref,L.def,u,E,nu,t);

        // Check: A = t C_lin, D = t^3/12 C_lin and B = 0 as in the control test, to identityTol()
        // relative to |A|, |D| and |A| t_tot respectively.
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const gsMatrix<> Aref = t*Clin.reshapeCol(k,3,3);
            const gsMatrix<> Dref = t*t*t/12*Clin.reshapeCol(k,3,3);
            CHECK_CLOSE(0.0, gsMatrix<>(A.reshapeCol(k,3,3)-Aref).norm()/Aref.norm(), tol);
            CHECK_CLOSE(0.0, gsMatrix<>(D.reshapeCol(k,3,3)-Dref).norm()/Dref.norm(), tol);
            CHECK_CLOSE(0.0, normAt(B,k)/(normAt(A,k)*t), tol);
        }
    }

    // The elasticity tensor has the major symmetry C^{abcd} = C^{cdab}, so the Voigt matrices of every
    // ply, and hence A and D (and B for a ply stack that is not mirror symmetric), are symmetric
    // on any patch, orthogonal or not.
    TEST(Skewed_ABD_Symmetric)
    {
        const real_t tol = identityTol();
        const gsMatrix<> u = pointsInterior();

        // Set-up: two orthotropic layups on the skewed, lifted patch, both with t_tot = 2*0.05 =
        // 4*0.025 = 0.1 so that the two runs differ only in the layup: the unsymmetric [0, pi/4]
        // (B != 0) and the mirror-symmetric [0, pi/4, pi/4, 0] (B = 0 up to round-off). The orthotropic
        // ply makes the matrices depend on the fibre angle; the guards confirm the skew is real.
        for (index_t layup = 0; layup != 2; ++layup)
        {
            const std::vector<real_t> angles = layup == 0 ? twoPlyAngles() : symmetricAngles();
            const real_t tply = layup == 0 ? 0.05 : 0.025;
            Laminate<3> L(angles,tply,orthotropicPly(),true);
            checkSkewGuards(geometryGuards(L.ref,u,true),true);

            // Act: the laminate matrices of the material; no oracle is needed because symmetry is
            // a property of each matrix alone.
            gsMatrix<> A, B, D;
            compositeABD(L,u,A,B,D);

            // Check: |M - M^T| vanishes relative to |A|, |D| and, for B, to the natural scale |A| t_tot
            // (a relative error |dB|/|B| would be undefined for the symmetric layup), to identityTol().
            for (index_t k = 0; k != u.cols(); ++k)
            {
                CHECK_CLOSE(0.0, asymmetryAt(A,k)/normAt(A,k), tol);
                CHECK_CLOSE(0.0, asymmetryAt(D,k)/normAt(D,k), tol);
                CHECK_CLOSE(0.0, asymmetryAt(B,k)/(normAt(A,k)*L.tTotal), tol);
            }
        }
    }

    // A and D of the laminate [0, pi/4, -pi/3] on the unit square: the metric is the identity, so
    // the tensor oracle reduces to the classical Voigt rotation and the fibre angle convention
    // (phi from e1 towards e2) is tested. D depends on z^2 only, so it does not see the stacking direction.
    TEST(Orthogonal_Laminate_AD_MatchesTensorOracle)
    {
        // Set-up: a thin three-ply orthotropic laminate (t_ply = 0.02, t_tot = 0.06; the oracle
        // identities are exact in t) with the three distinct angles of mixedAngles(), on the flat
        // unit square. No deformation is applied because A and D do not depend on it.
        const real_t tply = 0.02, tol = identityTol();
        const std::vector<real_t> angles = mixedAngles();
        Laminate<3> L(angles,tply,orthotropicPly(),false);
        const gsMatrix<> u = pointsInterior();

        // Act: the laminate matrices of the material and the tensor-product oracle for the same
        // plies, which rotates each ply matrix independently of the material code.
        gsMatrix<> A, B, D, Aor, Bor, Dor;
        compositeABD(L,u,A,B,D);
        tensorOracleAtPoints(L.ref,u,orthotropicPly(),tply,angles,Aor,Bor,Dor);

        // Check: A and D agree with the oracle to identityTol(), relative to the oracle norm.
        for (index_t k = 0; k != u.cols(); ++k)
        {
            CHECK_CLOSE(0.0, diffAt(A,Aor,k)/normAt(Aor,k), tol);
            CHECK_CLOSE(0.0, diffAt(D,Dor,k)/normAt(Dor,k), tol);
        }
    }

    // B of the unsymmetric laminate [0, pi/4, -pi/3] on the unit square. Ply 0 is at the most
    // negative z along the normal, so B = sum t_i z_i C_i with z_0 < 0. The oracle B must be
    // visibly nonzero, otherwise a sign error in B would pass.
    TEST(Orthogonal_Laminate_B_MatchesTensorOracle)
    {
        // Set-up: the laminate of the previous test (same angles, t_ply = 0.02, t_tot = 0.06, flat
        // unit square). The plies are not mirror symmetric, so B does not vanish.
        const real_t tply = 0.02, tol = identityTol();
        const std::vector<real_t> angles = mixedAngles();
        Laminate<3> L(angles,tply,orthotropicPly(),false);
        const gsMatrix<> u = pointsInterior();

        // Act: A and B of the material and the oracle B, whose plies are stacked from the most
        // negative z upwards by construction.
        gsMatrix<> A, B, D, Aor, Bor, Dor;
        compositeABD(L,u,A,B,D);
        tensorOracleAtPoints(L.ref,u,orthotropicPly(),tply,angles,Aor,Bor,Dor);

        // Check: B is non-vacuous and equal to the oracle, so a sign or stacking error cannot hide.
        for (index_t k = 0; k != u.cols(); ++k)
        {
            // Non-vacuity: the coupling is a sizeable fraction of its natural scale |A| t_tot.
            // |A| t_tot is the natural scale of B because |z| <= t_tot/2. For these angles with equal
            // plies at z in {-t_ply, 0, t_ply}, B = t_ply^2 (C_2 - C_0), non-zero because the two
            // outer plies differ. 1e-3 lies far above round-off and seven orders above the 1e-10
            // tolerance, while a sign or stacking error moves B by 2|B|, so it cannot hide below the
            // tolerance.
            CHECK(normAt(Bor,k) >= 1e-3*normAt(Aor,k)*L.tTotal);
            // B - B_or is measured in units of |A| t_tot, the natural scale of B, so that the same
            // check also holds for laminates where B is rounding noise.
            CHECK_CLOSE(0.0, diffAt(B,Bor,k)/(normAt(A,k)*L.tTotal), tol);
        }
    }

    // A, B and D of the unsymmetric laminate on the sheared, lifted patch against the tensor oracle,
    // which transforms the angle and the basis in one tensor contraction.
    TEST(Skewed_Laminate_ABD_MatchesTensorOracle)
    {
        // Set-up: the three-ply laminate of the orthogonal tests (t_ply = 0.02, t_tot = 0.06) on the
        // skewed, lifted patch, where the fibre angle and the non-orthogonal metric act together. The
        // guards confirm that the skew and the curvature are real.
        const real_t tply = 0.02, tol = identityTol();
        const std::vector<real_t> angles = mixedAngles();
        Laminate<3> L(angles,tply,orthotropicPly(),true);
        const gsMatrix<> u = pointsInterior();
        checkSkewGuards(geometryGuards(L.ref,u,true),true);

        // Act: the laminate matrices of the material and the tensor oracle on the same patch.
        gsMatrix<> A, B, D, Aor, Bor, Dor;
        compositeABD(L,u,A,B,D);
        tensorOracleAtPoints(L.ref,u,orthotropicPly(),tply,angles,Aor,Bor,Dor);

        // Check: A, B, D agree with the oracle to identityTol() (A and D relative to the oracle
        // norm, B relative to |A| t_tot), and B is non-vacuous (threshold as in the orthogonal B test).
        for (index_t k = 0; k != u.cols(); ++k)
        {
            CHECK_CLOSE(0.0, diffAt(A,Aor,k)/normAt(Aor,k), tol);
            CHECK_CLOSE(0.0, diffAt(D,Dor,k)/normAt(Dor,k), tol);
            CHECK(normAt(Bor,k) >= 1e-3*normAt(Aor,k)*L.tTotal);
            CHECK_CLOSE(0.0, diffAt(B,Bor,k)/(normAt(A,k)*L.tTotal), tol);
        }
    }

    // The stress resultants of a laminate with B != 0 are N = A eps + B kappa and M = B^T eps + D kappa
    // (classical lamination theory). Here eps and kappa are computed from the fundamental forms of the
    // reference and deformed patches, and A, B, D from the material itself, so the test isolates the
    // assembly of N and M from the correctness of the matrices.
    TEST(Skewed_NM_IncludeCoupling)
    {
        // Set-up: the unsymmetric two-ply laminate [0, pi/4] with t_ply = 0.05 (t_tot = 0.1), thick
        // enough that the coupling terms B kappa and B^T eps are not small against A eps and D kappa,
        // on the skewed, lifted patch. stretchedBent gives non-zero membrane strain and curvature
        // change at every point, so both N and M are non-trivial.
        const real_t tply = 0.05, tol = identityTol();
        Laminate<3> L(twoPlyAngles(),tply,orthotropicPly(),true);
        const gsMatrix<> u = pointsInterior();
        checkSkewGuards(geometryGuards(L.ref,u,true),true);
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);

        // Act: eps and kappa from the first and second fundamental forms of the two patches (not from
        // the material), the material's A, B, D, and its N and M.
        std::vector<gsVector<real_t,3> > eps, kap;
        strains(L.ref,L.def,u,eps,kap);

        gsMatrix<> A, B, D;
        compositeABD(L,u,A,B,D);
        const gsMatrix<> N = L.mm->eval3D_vector(0,u,zeros(u),MaterialOutput::VectorN);
        const gsMatrix<> M = L.mm->eval3D_vector(0,u,zeros(u),MaterialOutput::VectorM);

        // Check: N = A eps + B kappa and M = B^T eps + D kappa with all four terms present.
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const gsMatrix<> Ak = A.reshapeCol(k,3,3), Bk = B.reshapeCol(k,3,3), Dk = D.reshapeCol(k,3,3);
            const gsVector<> Bkap = Bk*kap[k];
            const gsVector<> Beps = Bk.transpose()*eps[k];
            const gsVector<> Aeps = Ak*eps[k];
            const gsVector<> Dkap = Dk*kap[k];

            // Non-vacuity: the coupling terms are not negligible against the direct terms. They must
            // be visibly present, not round-off: if either were dropped, N or M would change by at
            // least 1e-3 relative, far above the 1e-10 tolerance.
            CHECK(Bkap.norm() >= 1e-3*Aeps.norm());
            CHECK(Beps.norm() >= 1e-3*Dkap.norm());

            // Nscale and Mscale are scales of the same order as each summand: |A eps| <= |A||eps|, and
            // |B| is at most of order |A| t_tot because |z_i| <= t_tot/2 (an order-of-magnitude
            // scale, not a strict bound, since plies at different angles can partly cancel in the
            // sum of the C_i). Unlike |N| itself the scale cannot shrink through cancellation
            // between A eps and B kappa, so the relative error stays meaningful. Mscale is the same
            // construction for B^T eps + D kappa.
            const real_t Nscale = Ak.norm()*eps[k].norm() + Ak.norm()*L.tTotal*kap[k].norm();
            const real_t Mscale = Ak.norm()*L.tTotal*eps[k].norm() + Dk.norm()*kap[k].norm();
            CHECK_CLOSE(0.0, gsVector<>(N.col(k) - (Aeps + Bkap)).norm()/Nscale, tol);
            CHECK_CLOSE(0.0, gsVector<>(M.col(k) - (Beps + Dkap)).norm()/Mscale, tol);
        }
    }

    // A planar (2D) skewed patch with one ply at pi/4: A against the tensor oracle (which uses 2-vectors
    // here) and A = A^T.
    TEST(Skewed2D_A_MatchesTensorOracle)
    {
        // Set-up: a single orthotropic ply at pi/4 on the planar skewed patch. A is linear in t and
        // the checks are relative to |A|, so t only scales the result; 0.01 is a thin shell. With
        // Exx != Eyy every entry of the rotated ply matrix is non-zero at this angle, including the
        // normal-shear coupling entries (1,3) and (2,3), so the rotation is fully exercised. The
        // guards confirm the skew (no curvature in 2D).
        const real_t tply = 0.01, tol = identityTol();
        const std::vector<real_t> angles(1,compPi/4);
        Laminate<2> L(angles,tply,orthotropicPly(),true);
        const gsMatrix<> u = pointsInterior();
        checkSkewGuards(geometryGuards(L.ref,u,false),false);

        // Act: A of the material and A of the tensor oracle, which contracts with the 2-vectors of
        // the planar basis.
        const gsMatrix<> A = L.mm->eval3D_matrix(0,u,zeros(u),MaterialOutput::MatrixA);
        gsMatrix<> Aor, Bor, Dor;
        tensorOracleAtPoints(L.ref,u,orthotropicPly(),tply,angles,Aor,Bor,Dor);

        // Check: A agrees with the oracle and is symmetric, both to identityTol() relative to |A|.
        for (index_t k = 0; k != u.cols(); ++k)
        {
            CHECK_CLOSE(0.0, diffAt(A,Aor,k)/normAt(Aor,k), tol);
            CHECK_CLOSE(0.0, asymmetryAt(A,k)/normAt(A,k), tol);
        }
    }

    // The stacking direction is tied to the shell kinematics, independently of the laminate's own
    // z formula. Ply 0 lies at the most negative z along n = a1 x a2/|a1 x a2|. The base class
    // measures the strain at height z along that normal from the metric
    // Gcov(z) = A - 2 z B + z^2 n,^T n,, so E(z) = eps + z kappa + z^2 q exactly. The reference is the
    // through-thickness quadrature of gsMaterialMatrixLinear stresses evaluated at physical heights
    // z, minus the z^2 part, which lies outside the linear-in-z laminate model. The laminate must
    // match that reference, and the same plies stacked in the opposite order must miss it.
    TEST(Skewed_StackingDirection_MatchesKinematics)
    {
        // Set-up: two isotropic plies of different stiffness (the fibre angle is then irrelevant, and
        // B = (t_ply^2/2)(C_1 - C_0) != 0 only through the stiffness contrast), ply 0 = E[0] = 100 at
        // the bottom and ply 1 = E[1] = 300 on top; nu = 0.3 is generic and t_ply = 0.01 a thin shell.
        // Lmirror carries the plies in reverse order. Both sit on the skewed, lifted patch (guards)
        // with the deformation of stretchedBent, so eps and kappa are non-zero.
        const real_t tply = 0.01, nu = 0.3, tol = identityTol();
        std::vector<real_t> E(2);
        E[0] = 100.0; E[1] = 300.0;
        const std::vector<real_t> angles(2,0.0);
        std::vector<gsMatrix<> > G(2), Gmirror(2);
        for (index_t i = 0; i != 2; ++i)
        {
            G[i] = isotropicPly(E[i],nu);
            Gmirror[1-i] = G[i];
        }
        Laminate<3> L(angles,tply,G,true), Lmirror(angles,tply,Gmirror,true);
        const gsMatrix<> u = pointsInterior();
        checkSkewGuards(geometryGuards(L.ref,u,true),true);
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);
        stretchedBent(Lmirror.ref,Lmirror.def);
        Lmirror.mm->setDeformed(&Lmirror.def);

        // Act: the resultants N and M of the laminate and of its mirror image.
        const gsMatrix<> N  = L.mm->eval3D_vector(0,u,zeros(u),MaterialOutput::VectorN);
        const gsMatrix<> M  = L.mm->eval3D_vector(0,u,zeros(u),MaterialOutput::VectorM);
        const gsMatrix<> Nm = Lmirror.mm->eval3D_vector(0,u,zeros(u),MaterialOutput::VectorN);
        const gsMatrix<> Mm = Lmirror.mm->eval3D_vector(0,u,zeros(u),MaterialOutput::VectorM);

        // Reference: per ply, 2-point Gauss quadrature of the Linear stress at physical heights
        // (exact because S(z) = C(eps + z kappa + z^2 q) is quadratic in z, so the integrands of N
        // (degree 2) and M (degree 3) are integrated exactly), minus the z^2 stress C q, which
        // lies outside the laminate model and is extracted from the same Linear material by a central
        // second difference at 0 and +-h: C q = (S(h) + S(-h) - 2 S(0))/(2 h^2). Its integral
        // C q I2 (for N) or C q I3 (for M) is subtracted. The correction is even in z, so it
        // cannot fake the sign of the linear term. The second difference cancels terms of size |S|;
        // with h = t_tot/2 its absolute error is about eps |S|/h^2, which times I2 ~ t^3 is
        // ~ eps |S| t, round-off relative to the scales below. The heights are passed to Linear as
        // normalised heights zeta/t_tot because Linear is built with thickness t_tot. Only the
        // Stress output is used: the VectorN/VectorM outputs of Linear return the mid-surface
        // strains whatever the height. Linear uses the mid-surface constitutive tensor, the same
        // Kirchhoff-Love assumption as the laminate, so only the z-dependence of the kinematics is tested.
        const index_t npts = u.cols();
        const real_t tTot = L.tTotal, h = tTot/2, g = 1/math::sqrt(real_t(3));
        gsMatrix<> Nref = gsMatrix<>::Zero(3,npts), Mref = gsMatrix<>::Zero(3,npts);
        gsVector<> Nscale = gsVector<>::Zero(npts), Mscale = gsVector<>::Zero(npts);
        for (index_t i = 0; i != 2; ++i)
        {
            // Ply i occupies [-t_tot/2 + i t_ply, -t_tot/2 + (i+1) t_ply] along +n. The heights are the
            // two Gauss points of this interval (mid +- (t_ply/2)/sqrt(3)) followed by 0, +h and -h
            // for the second difference.
            const real_t zlo = -tTot/2 + i*tply, zhi = zlo + tply, zmid = (zlo + zhi)/2;
            gsVector<real_t,5> zeta;
            zeta << zmid - g*tply/2, zmid + g*tply/2, 0, h, -h;
            gsMatrix<> Z(5,npts);
            for (index_t j = 0; j != 5; ++j)
                Z.row(j).setConstant(zeta[j]/tTot);

            gsConstantFunction<> tfun(tTot,3), Efun(E[i],3), nufun(nu,3);
            gsMaterialMatrixLinear<3,real_t> lin(L.ref,tfun,Efun,nufun);
            lin.setDeformed(&L.def);
            const gsMatrix<> S = lin.eval3D_stress(0,u,Z,MaterialOutput::Stress);

            // I2 = int z^2 dz and I3 = int z^3 dz over the ply interval.
            const real_t I2 = (math::pow(zhi,3) - math::pow(zlo,3))/3;
            const real_t I3 = (math::pow(zhi,4) - math::pow(zlo,4))/4;
            for (index_t k = 0; k != npts; ++k)
            {
                const gsVector<> Cq = (S.col(3*npts+k) + S.col(4*npts+k) - 2*S.col(2*npts+k))/(2*h*h);
                for (index_t j = 0; j != 2; ++j)
                {
                    // The Gauss weights of a 2-point rule on an interval of length t_ply are t_ply/2 each.
                    // Nscale and Mscale accumulate the magnitudes of the summands (no cancellation).
                    const real_t w = tply/2;
                    Nref.col(k) += w*S.col(j*npts+k);
                    Mref.col(k) += w*zeta[j]*S.col(j*npts+k);
                    Nscale[k]   += w*S.col(j*npts+k).norm();
                    Mscale[k]   += w*math::abs(zeta[j])*S.col(j*npts+k).norm();
                }
                Nref.col(k) -= I2*Cq;
                Mref.col(k) -= I3*Cq;
            }
        }

        // Check: the laminate matches the reference to identityTol() relative to the sums of
        // magnitudes Nscale, Mscale (not to |N_ref|, which can shrink through cancellation). The
        // mirrored stack has the same A and D but B -> -B, so it misses the reference by 2 B kappa in
        // N and 2 B^T eps in M. Requiring that miss to be at least 1e-6 of the scale places it four
        // orders above the agreement tolerance in double (identityTol() = 1e-10), so the test decides
        // which ply order is right and not only that two computations are close. The bar is a fixed
        // number because the mirror miss is a property of the laminate (a fraction of the scale, at
        // most about 2), not of round-off; a bar tied to identityTol() would exceed that miss when
        // real_t is float, where identityTol() is much larger.
        for (index_t k = 0; k != npts; ++k)
        {
            CHECK_CLOSE(0.0, gsVector<>(N.col(k) - Nref.col(k)).norm()/Nscale[k], tol);
            CHECK_CLOSE(0.0, gsVector<>(M.col(k) - Mref.col(k)).norm()/Mscale[k], tol);
            CHECK(gsVector<>(Nm.col(k) - Nref.col(k)).norm()/Nscale[k] >= real_t(1e-6));
            CHECK(gsVector<>(Mm.col(k) - Mref.col(k)).norm()/Mscale[k] >= real_t(1e-6));
        }
    }
}
