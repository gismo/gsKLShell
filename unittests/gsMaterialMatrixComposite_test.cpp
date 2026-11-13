/** @file gsMaterialMatrixComposite_test.cpp

    @brief Checks the laminate strain energy of gsMaterialMatrixComposite<3>::eval3D_energy.

    The energy per unit reference area is W = 1/2 e.A.e + e.B.k + 1/2 k.D.k with the
    membrane strain e and curvature change k computed here from the geometry. It is
    compared with the laminate stiffness matrices (A, B, D), the stress resultants
    (N = A e, M = D k) and central differences along deformation paths
    (dW/dt = N.e' + M.k' for a symmetric layup, where B = 0).

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

const real_t compPi = math::atan(1)*4;

// Biquadratic, unscaled, axis-aligned unit square embedded in R^3. On it the
// covariant and contravariant metrics are the identity, so the laminate matrices
// are symmetric. The parameter functions are members because the material keeps
// (non-owned) pointers to them, to the reference patch and to the deformed patch.
struct CompositeLaminate
{
    gsMultiPatch<> ref, def;
    gsConstantFunction<> Gfun;
    std::vector<gsConstantFunction<> > thicks, phis;
    std::vector<gsFunctionSet<>*> Ts, Gs, Phis;
    real_t tTotal;
    memory::unique_ptr<gsMaterialMatrixComposite<3,real_t> > mm;

    // Plies of equal thickness tply with the given fibre angles, bottom to top.
    CompositeLaminate(const std::vector<real_t> & angles, real_t tply)
    :
    tTotal(tply*angles.size())
    {
        ref.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
        ref.degreeElevate(1);
        ref.embed(3);
        def.addPatch( ref.patch(0).clone() );

        // Orthotropic ply; nuyx*Exx == nuxy*Eyy holds exactly in floating point.
        gsMatrix<> Gmat = gsCompositeMatrix<real_t>(200.0,100.0,40.0,0.25,0.125);
        Gmat.resize(Gmat.rows()*Gmat.cols(),1);
        Gfun = gsConstantFunction<>(Gmat,3);

        const size_t n = angles.size();
        thicks.reserve(n);
        phis.reserve(n);
        for (size_t i = 0; i != n; ++i)
        {
            thicks.push_back(gsConstantFunction<>(tply,3));
            phis.push_back(gsConstantFunction<>(angles[i],3));
        }
        for (size_t i = 0; i != n; ++i)
        {
            Ts.push_back(&thicks[i]);
            Gs.push_back(&Gfun);
            Phis.push_back(&phis[i]);
        }
        mm = memory::make_unique(new gsMaterialMatrixComposite<3,real_t>(ref,Ts,Gs,Phis));
        mm->setDeformed(&def);
    }

private:
    CompositeLaminate(const CompositeLaminate &);
    CompositeLaminate & operator=(const CompositeLaminate &);
};

std::vector<real_t> symmetricAngles()
{
    std::vector<real_t> a(4);
    a[0] = 0.0; a[1] = compPi/4; a[2] = compPi/4; a[3] = 0.0;
    return a;
}

std::vector<real_t> unsymmetricAngles()
{
    std::vector<real_t> a(2);
    a[0] = 0.0; a[1] = compPi/4;
    return a;
}

// Covariant first and second fundamental forms of patch 0 at the points u,
// computed as in gsMaterialMatrixBaseDim.
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

gsVector<real_t,3> fold(const gsMatrix<> & E)
{
    gsVector<real_t,3> v;
    v << E(0,0), E(1,1), E(0,1)+E(1,0);
    return v;
}

// Green-Lagrange membrane strain e = fold(1/2 (A_def - A_ref)) and curvature
// change k = fold(B_ref - B_def) per point (engineering Voigt order).
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

gsMatrix<> pointsInterior()
{
    gsMatrix<> u(2,3);
    u << 0.2, 0.5, 0.8,
         0.3, 0.6, 0.4;
    return u;
}

// Stretch, shear and out-of-plane bending of the reference patch.
void stretchedBent(const gsMultiPatch<> & ref, gsMultiPatch<> & def)
{
    gsMatrix<> & c = def.patch(0).coefs();
    c = ref.patch(0).coefs();
    c.col(0) *= 1.05;
    c.col(1) += 0.03*ref.patch(0).coefs().col(0);
    c(4,2) += 0.5;
    c(8,2) += 0.25;
}

gsMatrix<> zeros(const gsMatrix<> & u)
{
    return gsMatrix<>::Zero(1,u.cols());
}

} // namespace

SUITE(gsMaterialMatrixComposite_test)                 // The suite should have the same name as the file
{
    TEST(Energy_SymmetricLayup_NM)
    {
        CompositeLaminate L(symmetricAngles(),0.025);
        const gsMatrix<> u = pointsInterior();
        const gsMatrix<> z = zeros(u);
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);

        std::vector<gsVector<real_t,3> > eps, kap;
        strains(L.ref,L.def,u,eps,kap);

        const gsMatrix<> W = L.mm->eval3D_energy(0,u,z);
        const gsMatrix<> N = L.mm->eval3D_vector(0,u,z,MaterialOutput::VectorN);
        const gsMatrix<> M = L.mm->eval3D_vector(0,u,z,MaterialOutput::VectorM);
        const gsMatrix<> A = L.mm->eval3D_matrix(0,u,z,MaterialOutput::MatrixA);
        const gsMatrix<> B = L.mm->eval3D_matrix(0,u,z,MaterialOutput::MatrixB);

        // The layup is symmetric: B vanishes up to roundoff.
        CHECK(B.cwiseAbs().maxCoeff() <= 1e-12 * A.cwiseAbs().maxCoeff() * L.tTotal);

        CHECK_EQUAL(1, W.rows());
        CHECK_EQUAL(u.cols(), W.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const real_t tn = eps[k].dot(N.col(k));
            const real_t tm = kap[k].dot(M.col(k));
            const real_t Wref = 0.5*(tn + tm);
            const real_t sum = math::abs(tn) + math::abs(tm);
            CHECK(math::abs(W(0,k) - Wref) <= 1e-10*sum + 1e-14);
            // Both contributions matter, so a wrong factor on either one is visible.
            CHECK(math::abs(tn) >= 1e-2*sum);
            CHECK(math::abs(tm) >= 1e-2*sum);
        }
    }

    TEST(Energy_UnsymmetricLayup_ABD)
    {
        CompositeLaminate L(unsymmetricAngles(),0.05);
        const gsMatrix<> u = pointsInterior();
        const gsMatrix<> z = zeros(u);
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);

        std::vector<gsVector<real_t,3> > eps, kap;
        strains(L.ref,L.def,u,eps,kap);

        const gsMatrix<> W  = L.mm->eval3D_energy(0,u,z);
        const gsMatrix<> Am = L.mm->eval3D_matrix(0,u,z,MaterialOutput::MatrixA);
        const gsMatrix<> Bm = L.mm->eval3D_matrix(0,u,z,MaterialOutput::MatrixB);
        const gsMatrix<> Dm = L.mm->eval3D_matrix(0,u,z,MaterialOutput::MatrixD);

        // The coupling must be present for the B term to be tested at all.
        CHECK(Bm.cwiseAbs().maxCoeff() >= 1e-3 * Am.cwiseAbs().maxCoeff() * L.tTotal);

        for (index_t k = 0; k != u.cols(); ++k)
        {
            const gsMatrix<> A = Am.reshapeCol(k,3,3);
            const gsMatrix<> B = Bm.reshapeCol(k,3,3);
            const gsMatrix<> D = Dm.reshapeCol(k,3,3);

            // The finite-difference tests rely on symmetric A and D.
            CHECK((A - A.transpose()).cwiseAbs().maxCoeff() <= 1e-12*A.cwiseAbs().maxCoeff());
            CHECK((D - D.transpose()).cwiseAbs().maxCoeff() <= 1e-12*D.cwiseAbs().maxCoeff());

            const real_t ta = 0.5*eps[k].dot(A*eps[k]);
            const real_t tb = eps[k].dot(B*kap[k]);
            const real_t td = 0.5*kap[k].dot(D*kap[k]);
            const real_t Wref = ta + tb + td;
            const real_t sum = math::abs(ta) + math::abs(tb) + math::abs(td);

            CHECK(math::abs(W(0,k) - Wref) <= 1e-10*sum + 1e-14);
            CHECK(math::abs(tb) >= 1e-2*math::abs(Wref));
            CHECK(math::abs(ta) >= 1e-2*sum);
            CHECK(math::abs(td) >= 1e-2*sum);
        }
    }

    TEST(Energy_FD_Membrane)
    {
        CompositeLaminate L(symmetricAngles(),0.025);
        const gsMatrix<> u = pointsInterior();
        const gsMatrix<> z = zeros(u);
        const real_t t0 = 1.0, h = 1e-6;

        // Flat in-plane path: z stays 0, hence k = 0 along it.
        gsMatrix<> base = L.ref.patch(0).coefs();
        base.col(0) *= 1.03;
        const gsMatrix<> dirX = L.ref.patch(0).coefs().col(0);
        auto setState = [&](real_t t)
        {
            L.def.patch(0).coefs() = base;
            L.def.patch(0).coefs().col(0) *= 1.0 + 0.03*t;
            L.def.patch(0).coefs().col(1) += 0.02*t*dirX;
            L.mm->setDeformed(&L.def);
        };

        std::vector<gsVector<real_t,3> > eps0, kap0, epsP, epsM, kapP, kapM;
        setState(t0);
        strains(L.ref,L.def,u,eps0,kap0);
        const gsMatrix<> N = L.mm->eval3D_vector(0,u,z,MaterialOutput::VectorN);
        const gsMatrix<> M = L.mm->eval3D_vector(0,u,z,MaterialOutput::VectorM);

        setState(t0+h);
        strains(L.ref,L.def,u,epsP,kapP);
        const gsMatrix<> Wp = L.mm->eval3D_energy(0,u,z);
        setState(t0-h);
        strains(L.ref,L.def,u,epsM,kapM);
        const gsMatrix<> Wm = L.mm->eval3D_energy(0,u,z);

        for (index_t k = 0; k != u.cols(); ++k)
        {
            const gsVector<real_t,3> de = (epsP[k]-epsM[k])/(2*h);
            const gsVector<real_t,3> dk = (kapP[k]-kapM[k])/(2*h);
            const real_t dW   = (Wp(0,k)-Wm(0,k))/(2*h);
            const real_t tn   = N.col(k).dot(de);
            const real_t tm   = M.col(k).dot(dk);
            const real_t pred = tn + tm;

            CHECK(math::abs(dW - pred) <= 1e-6*(math::abs(tn)+math::abs(tm)) + 1e-12);
            CHECK(math::abs(tn) > 1e-3);
            CHECK(kap0[k].cwiseAbs().maxCoeff() <= 1e-14);
        }
    }

    TEST(Energy_FD_Bending)
    {
        CompositeLaminate L(symmetricAngles(),0.025);
        const gsMatrix<> u = pointsInterior();
        const gsMatrix<> z = zeros(u);
        const real_t t0 = 1.0, h = 1e-6;

        // Out-of-plane path from a curved state, with a small in-plane stretch.
        gsMatrix<> base = L.ref.patch(0).coefs();
        base.col(0) *= 1.001;
        auto setState = [&](real_t t)
        {
            L.def.patch(0).coefs() = base;
            L.def.patch(0).coefs()(4,2) = t*0.05;
            L.mm->setDeformed(&L.def);
        };

        std::vector<gsVector<real_t,3> > eps0, kap0, epsP, epsM, kapP, kapM;
        setState(t0);
        strains(L.ref,L.def,u,eps0,kap0);
        const gsMatrix<> N = L.mm->eval3D_vector(0,u,z,MaterialOutput::VectorN);
        const gsMatrix<> M = L.mm->eval3D_vector(0,u,z,MaterialOutput::VectorM);

        setState(t0+h);
        strains(L.ref,L.def,u,epsP,kapP);
        const gsMatrix<> Wp = L.mm->eval3D_energy(0,u,z);
        setState(t0-h);
        strains(L.ref,L.def,u,epsM,kapM);
        const gsMatrix<> Wm = L.mm->eval3D_energy(0,u,z);

        for (index_t k = 0; k != u.cols(); ++k)
        {
            const gsVector<real_t,3> de = (epsP[k]-epsM[k])/(2*h);
            const gsVector<real_t,3> dk = (kapP[k]-kapM[k])/(2*h);
            const real_t dW   = (Wp(0,k)-Wm(0,k))/(2*h);
            const real_t tn   = N.col(k).dot(de);
            const real_t tm   = M.col(k).dot(dk);
            const real_t pred = tn + tm;

            CHECK(math::abs(dW - pred) <= 1e-6*(math::abs(tn)+math::abs(tm)) + 1e-12);
            // Bending dominates, so a wrong factor on the D term shows in dW/dt.
            CHECK(math::abs(tm) >= 0.5*(math::abs(tn)+math::abs(tm)));
        }
    }

    TEST(Energy_Undeformed_Zero)
    {
        const std::vector<std::vector<real_t> > layups = {symmetricAngles(),unsymmetricAngles()};
        for (size_t l = 0; l != layups.size(); ++l)
        {
            CompositeLaminate L(layups[l],0.05);
            const gsMatrix<> u = pointsInterior();
            const gsMatrix<> z = zeros(u);

            // Flat reference, then a curved reference with the same state.
            for (int curved = 0; curved != 2; ++curved)
            {
                if (curved)
                {
                    L.ref.patch(0).coefs()(4,2) += 0.1;
                    L.ref.patch(0).coefs()(8,2) += 0.05;
                }
                L.def.patch(0).coefs() = L.ref.patch(0).coefs();
                L.mm->setDeformed(&L.def);

                const gsMatrix<> W = L.mm->eval3D_energy(0,u,z);
                for (index_t k = 0; k != u.cols(); ++k)
                    CHECK(math::abs(W(0,k)) <= 1e-14);
            }
        }
    }

    TEST(Energy_OutputShape)
    {
        CompositeLaminate L(unsymmetricAngles(),0.05);
        const gsMatrix<> u = pointsInterior();
        const index_t npts = u.cols(), nz = 3;
        stretchedBent(L.ref,L.def);
        L.mm->setDeformed(&L.def);

        gsMatrix<> z(nz,npts);
        z.row(0).setConstant(-0.5);
        z.row(1).setConstant(0.0);
        z.row(2).setConstant(0.5);

        const gsMatrix<> W  = L.mm->eval3D_energy(0,u,z);
        const gsMatrix<> W0 = L.mm->eval3D_energy(0,u,zeros(u));

        CHECK_EQUAL(1, W.rows());
        CHECK_EQUAL(nz*npts, W.cols());
        CHECK_EQUAL(npts, W0.cols());
        for (index_t k = 0; k != npts; ++k)
        {
            CHECK(math::abs(W0(0,k)) > 0.0);
            CHECK_EQUAL(W0(0,k), W(0,k));
            for (index_t j = 0; j != nz; ++j)
                CHECK_EQUAL(W(0,k), W(0,j*npts+k));
        }
    }
}
