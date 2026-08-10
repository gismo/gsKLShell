/** @file gsMaterialMatrix3D.hpp

    @brief Implementation of the legacy-API adapter over 3D PFF material laws.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

#include <gsKLShell/src/gsMaterialMatrix3D.h>

#ifdef gsPhaseFieldFracture_ENABLED

namespace gismo
{

template <short_t dim, class T>
gsMaterialMatrix3D<dim,T>::gsMaterialMatrix3D(const gsFunctionSet<T> & mp,
                                              const gsFunctionSet<T> & thickness,
                                              const gsMaterialBase<T> & law)
:
Base(&mp,&thickness,nullptr),
m_law(&law),
m_condensation(&law)
{
    // No shell-side parameter functions: the law owns its own parameters, so
    // Base::m_pars stays empty and _computePoints fills a 0 x N m_parmat.
}

template <short_t dim, class T>
gsMaterialMatrix3D<dim,T>::gsMaterialMatrix3D(const gsFunctionSet<T> & mp,
                                              const gsFunctionSet<T> & thickness,
                                              const gsFunctionSet<T> & Density,
                                              const gsMaterialBase<T> & law)
:
Base(&mp,&thickness,&Density),
m_law(&law),
m_condensation(&law)
{
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrix3D<dim,T>::_localCartesianTriad(const gsMatrix<T> & gcov) const
{
    // Local ORTHONORMAL triad from the UNDEFORMED covariant basis (columns are
    // the covariant basis vectors g_1, g_2, g_3). e3 is the shell normal.
    // CORRECTNESS: the constitutive plane-stress structure (E13=E23=0, condense
    // along the 3-direction) only holds in a frame whose 3rd axis is the shell
    // normal; this is why we must NOT use the global identity as the Cartesian
    // basis (as eval3D_cov2cart does for mere post-processing).
    gsMatrix<T> triad(3,3);
    gsVector<T> e1 = gcov.col(0);
    e1.normalize();
    gsVector<T> e2 = gcov.col(1) - (gcov.col(1).dot(e1)) * e1; // Gram-Schmidt
    e2.normalize();
    gsVector<T> e3 = gcov.col(2);                              // shell normal
    e3.normalize();
    triad.col(0) = e1;
    triad.col(1) = e2;
    triad.col(2) = e3;
    return triad;
}

template <short_t dim, class T>
void gsMaterialMatrix3D<dim,T>::_computeBatch(const index_t patch,
                                              const gsMatrix<T> & u,
                                              const gsMatrix<T> & z) const
{
    // Metrics + thickness + (empty) parameters; task-11 same-input guard makes a
    // repeat call free.
    this->_computePoints(patch,u);

    // --- Cross-output cache probe (bump-on-set semantics only) ---------------
    Entry & e = m_cache.mine();
    if (   e.patch == patch && e.rev == Base::m_configRev
        && e.u.rows() == u.rows() && e.u.cols() == u.cols() && e.u == u
        && e.z.rows() == z.rows() && e.z.cols() == z.cols() && e.z == z)
    {
        gsMaterialMatrix3DIncrementHits();
        return;
    }
    gsMaterialMatrix3DIncrementSweeps();
    e.patch = patch; e.rev = Base::m_configRev; e.u = u; e.z = z;

    const index_t N  = u.cols();
    const index_t nz = z.rows();
    const index_t P  = N * nz;

    // --- Parameters at the midsurface, replicated per z-row ------------------
    // Parameters depend only on the in-plane point; the law evaluates its OWN
    // parameter functions at the midsurface points u.
    gsMaterialData<T> pdata;
    m_law->precomputeParameters(patch,u,pdata);
    const index_t npar = static_cast<index_t>(pdata.parameters.size());
    std::vector<gsMatrix<T>> params(npar);
    for (index_t v=0; v!=npar; ++v)
    {
        params[v].resize(1,P);
        for (index_t k=0; k!=N; ++k)
            for (index_t j=0; j!=nz; ++j)
                params[v](0, j*N+k) = pdata.parameters[v](0,k);
    }

    // --- First pass: covariant strain -> local Cartesian; cache back-transform
    gsMatrix<T> e2D(3,P);
    std::vector<gsMatrix<T>> Tback(P); // per point: cart -> contravariant curvilinear
    gsMatrix<T> ecov(3,1);
    for (index_t k=0; k!=N; ++k)
    {
        for (index_t j=0; j!=nz; ++j)
        {
            const index_t col = j*N+k;
            // PHYSICAL through-thickness height z*T (dimensionless z in [-1/2,1/2]).
            this->_getMetric(k, z(j,k) * m_data.mine().m_Tmat(0,k));

            // Covariant strain, shell Voigt-3 [E11, E22, E01+E10] (engineering),
            // pattern of gsMaterialMatrixBaseDim::eval3D_strain.
            const gsMatrix<T> Eblk = 0.5 * ( m_data.mine().m_Gcov_def.block(0,0,2,2)
                                           - m_data.mine().m_Gcov_ori.block(0,0,2,2) );
            ecov(0,0) = Eblk(0,0);
            ecov(1,0) = Eblk(1,1);
            ecov(2,0) = Eblk(0,1) + Eblk(1,0);

            const gsMatrix<T> triad = _localCartesianTriad(m_data.mine().m_gcov_ori);

            // CALIBRATION POINT (task 15 SvK oracle). The strain tensor is
            // E = E_ij g^i (x) g^j, so its covariant COMPONENTS live on the
            // CONTRAVARIANT in-plane basis g^1,g^2 (= m_gcon_ori cols 0,1).
            // Projecting onto the orthonormal Cartesian triad e1,e2 gives the
            // engineering-Voigt STRAIN transform R (e_cart = R e_cov), with
            // a_ia = g^i . e_a:
            //   E11cart = a11^2 E11 + a21^2 E22 + a11 a21 (2E12)
            //   E22cart = a12^2 E11 + a22^2 E22 + a12 a22 (2E12)
            //  2E12cart = 2a11a12 E11 + 2a21a22 E22 + (a11a22+a21a12)(2E12)
            // The contravariant factor a_ia ~ 1/|g| gives R ~ 1/|g|^2, so the
            // metric-scaled covariant strain (~|g|^2 * strain) maps to an O(strain)
            // Cartesian strain — feeding a well-scaled state to the plane-stress
            // condenser. By energy conjugacy the work-conjugate CONTRAVARIANT stress
            // and tangent transform with R^T: S^ij = R^T S_cart,
            // C^ijkl = R^T C_cart R (symmetric), matching the legacy contravariant
            // convention (gsMaterialMatrixLinear.hpp:554).
            //
            // Task-15 CALIBRATION FIX (both marked points, proven by the curved
            // z=0 SvK oracle SvK_curved_z0_pointwise): the original code applied
            // gsMaterialMatrixBaseDim::_transformation (which builds the engineering
            // STRESS transform T_sigma, factor-2 in the shear COLUMN) to the STRAIN
            // and used a metric-inflating gcov orientation. That (a) aborted the
            // condenser (det(2E+I)<0, ~|g|^2 inflation) and (b) after the basis flip
            // still left an O(1e-1) shear-convention error on a NON-orthogonal patch
            // (masked to ~1e-6 where g1 ~= g2-orthogonal). Strain and stress are
            // CONTRAGREDIENT (T_eps = R here, T_sigma = R^{-T} only for rotations —
            // NOT for non-orthonormal curvilinear bases), so R is built explicitly.
            const T a11 = m_data.mine().m_gcon_ori.col(0).dot(triad.col(0));
            const T a12 = m_data.mine().m_gcon_ori.col(0).dot(triad.col(1));
            const T a21 = m_data.mine().m_gcon_ori.col(1).dot(triad.col(0));
            const T a22 = m_data.mine().m_gcon_ori.col(1).dot(triad.col(1));
            gsMatrix<T> R(3,3);
            R(0,0) = a11*a11; R(0,1) = a21*a21; R(0,2) = a11*a21;
            R(1,0) = a12*a12; R(1,1) = a22*a22; R(1,2) = a12*a22;
            R(2,0) = 2*a11*a12; R(2,1) = 2*a21*a22; R(2,2) = a11*a22 + a21*a12;
            e2D.col(col) = R * ecov;
            Tback[col]   = R.transpose();
        }
    }

    // --- ONE batched plane-stress condensation over the whole grid -----------
    gsMatrix<T> S2Dc, C2Dc, E33;
    m_condensation.condense(e2D, params, S2Dc, C2Dc, E33);

    // --- Second pass: back-transform to contravariant curvilinear components -
    e.S.resize(3,P);
    e.C.resize(9,P);
    for (index_t col=0; col!=P; ++col)
    {
        const gsMatrix<T> & Tb = Tback[col];
        e.S.col(col) = Tb * S2Dc.col(col);
        const gsAsMatrix<T,Dynamic,Dynamic> Ccart = C2Dc.reshapeCol(col,3,3);
        gsAsMatrix<T,Dynamic,Dynamic> Ccurv = e.C.reshapeCol(col,3,3);
        Ccurv = Tb * Ccart * Tb.transpose();
    }
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrix3D<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    // NotIntegrated: the pointwise 9 x P tangent; the integrator's moment
    // weights build A/B/C/D from the SAME pointwise data.
    this->_computeBatch(patch,u,z);
    return m_cache.mine().C;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrix3D<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    // NotIntegrated: the pointwise 3 x P stress; the integrator's moment weights
    // build N/M from the SAME pointwise data.
    this->_computeBatch(patch,u,z);
    return m_cache.mine().S;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrix3D<dim,T>::eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput /*out*/) const
{
    // The stress IS the vector output's data (same condensed S^ij).
    this->_computeBatch(patch,u,z);
    return m_cache.mine().S;
}

template <short_t dim, class T>
std::ostream & gsMaterialMatrix3D<dim,T>::print(std::ostream & os) const
{
    os  <<"---------------------------------------------------------------------\n"
        <<"----------------gsMaterialMatrix3D (PFF adapter) Info----------------\n"
        <<"---------------------------------------------------------------------\n\n";
    os  <<"Wrapped 3D PFF law with "
        <<(m_law!=nullptr ? m_law->numParameters() : 0)<<" parameter(s).\n";
    return os;
}

} // namespace gismo

#endif // gsPhaseFieldFracture_ENABLED
