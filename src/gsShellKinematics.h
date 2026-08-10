/** @file gsShellKinematics.h

    @brief Stateless Kirchhoff-Love shell metric engine driven by an INJECTED
           gsMapData (no geometry ownership, no configuration state).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsMath.h>
#include <gsCore/gsDebug.h>

namespace gismo
{

/**
 * @brief      Shell metric (kinematics) engine: fundamental forms and
 *             through-thickness metric tensors of a Kirchhoff-Love shell,
 *             computed from a CALLER-SUPPLIED \ref gsMapData.
 *
 * This class is the kinematics core of the shell material provider. It
 * reproduces the fundamental-form / metric computations of
 * \ref gsMaterialMatrixBaseDim EXACTLY (the formulas are transplanted verbatim,
 * see the per-method @c \@note tags for the source lines), with two deliberate
 * differences:
 *
 * 1. **Injected geometry data.** It never calls @c computeMap itself: the
 *    caller computes ONE \ref gsMapData per configuration and hands it in.
 *    That is what allows several material laws (and several outputs of the same
 *    law) to share a single geometry evaluation instead of each re-running the
 *    map, which is the whole point of the provider.
 * 2. **Stateless by design.** There is NO configuration state here: no geometry
 *    pointers, no thickness, no material parameters, no options, no
 *    same-input guard and no configuration revision. The provider's correctness
 *    depends on there being *nothing to invalidate*: the only state is the
 *    batch of metric quantities produced by the last @c computeUndeformed /
 *    @c computeDeformed call, and those are overwritten wholesale. Consequently
 *    the class has no setters.
 *
 * ### Usage
 * \code
 * gsMapData<T> mapOri, mapDef;
 * mapOri.flags = mapDef.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
 * mapOri.points = mapDef.points = u;
 * // ... computeMap on the undeformed / deformed geometry ...
 * gsShellKinematics<3,T> kin;
 * kin.computeUndeformed(mapOri);      // batched: ALL points of the map
 * kin.computeDeformed  (mapDef);
 * typename gsShellKinematics<3,T>::PointMetrics m;
 * kin.getMetric(k, z*thickness, m);   // z*thickness == PHYSICAL height
 * \endcode
 *
 * ### Height convention
 * \c getMetric takes the **PHYSICAL** through-thickness coordinate, i.e. the
 * caller multiplies the dimensionless \f$ z\in[-\tfrac12,\tfrac12]\f$ by the
 * thickness itself (exactly as \ref gsMaterialMatrix3D does with
 * <tt>z(j,k)*m_Tmat(0,k)</tt>). Thickness handling is NOT this class' business.
 *
 * ### Required map flags
 * - BOTH dimensions: \c NEED_VALUE | \c NEED_JACOBIAN.
 * - \c dim==3 additionally: \c NEED_DERIV2 | \c NEED_NORMAL.
 *
 * \c NEED_VALUE is required in BOTH dimensions even though the metric formulas
 * never read the geometry values: <tt>map.values[0].cols()</tt> is the only
 * usable source for the number of points. \c map.points is NOT usable --
 * \c gsExprHelper::precompute swaps the quadrature points into it
 * (gsExprHelper.h:451) and straight back out after \c computeMap (:455), so it
 * is 0 x 0 for every element of an expression-assembler run. See \ref _nPoints.
 *
 * ### Dimensions
 * - \c dim==3 : surfaces in 3D (curved shells, \f$B\neq 0\f$, normal from the map).
 * - \c dim==2 : planar geometries (\f$B=0\f$, \f$n=e_3\f$); the \c dim==2 path
 *   never touches \c map.normals nor \c map.deriv2, so a planar caller may omit
 *   \c NEED_NORMAL / \c NEED_DERIV2 (but NOT \c NEED_VALUE, see above).
 *
 * @tparam     dim   Geometric dimension (2: planar, 3: surface)
 * @tparam     T     Real type
 *
 * @ingroup    KLShell
 */
template <short_t dim, class T>
class gsShellKinematics
{
public:

    /// Per-point, per-height metric quantities; the consumer contract of this class.
    struct PointMetrics
    {
        /// Covariant metric tensor \f$G_{ij}\f$ at height z, undeformed (3x3)
        gsMatrix<T> Gcov_ori;
        /// Covariant metric tensor \f$G_{ij}\f$ at height z, deformed (3x3)
        gsMatrix<T> Gcov_def;
        /// Covariant basis vectors \f$g_i\f$ at height z, undeformed (3x3, column 2 = normal)
        gsMatrix<T> gcov_ori;
        /// Contravariant basis vectors \f$g^i\f$ at height z, undeformed (3x3, column 2 = normal direction)
        gsMatrix<T> gcon_ori;
        /// \f$J_0^2 = \det(G_{cov,def})/\det(G_{cov,ori})\f$ at height z
        T J0_sq;

        PointMetrics() : J0_sq(0) { }
    };

public:

    gsShellKinematics()
    :
    m_nPointsOri(0),
    m_nPointsDef(0)
    { }

    /**
     * @brief      Computes the undeformed fundamental forms for ALL points of \a map (batched).
     *
     * @param[in]  map   Geometry map of the UNDEFORMED configuration; must have been
     *                   computed with NEED_VALUE | NEED_JACOBIAN (and, for \c dim==3,
     *                   additionally NEED_DERIV2 | NEED_NORMAL). The number of points
     *                   is taken from <tt>map.values[0]</tt>, see \ref _nPoints --
     *                   which is why NEED_VALUE is not optional.
     */
    void computeUndeformed(const gsMapData<T> & map)
    { _computeMetricUndeformed_impl<dim>(map); }

    /**
     * @brief      Computes the deformed fundamental forms for ALL points of \a map (batched).
     *
     * @param[in]  map   Geometry map of the DEFORMED configuration; same flag
     *                   requirements as \ref computeUndeformed.
     */
    void computeDeformed(const gsMapData<T> & map)
    { _computeMetricDeformed_impl<dim>(map); }

    /**
     * @brief      Evaluates the metric quantities at in-plane point \a k and
     *             PHYSICAL through-thickness height \a zPhysical.
     *
     * @param[in]  k          In-plane point index (column of the injected map)
     * @param[in]  zPhysical  Through-thickness coordinate in PHYSICAL units
     * @param[out] out        Filled with the quantities of \ref PointMetrics
     *
     * @note Transplanted from gsMaterialMatrixBaseDim<dim,T>::_getMetric(k,z),
     *       gsMaterialMatrixBaseDim.hpp:1066-1085 (the J0_sq part).
     */
    void getMetric(index_t k, T zPhysical, PointMetrics & out) const
    {
        GISMO_ENSURE(m_nPointsOri!=0,"Undeformed metric is not initialized; call computeUndeformed first.");
        GISMO_ENSURE(m_nPointsDef!=0,"Deformed metric is not initialized; call computeDeformed first.");
        GISMO_ENSURE(m_nPointsOri==m_nPointsDef,"Undeformed and deformed maps have a different number of points: "
                                                <<m_nPointsOri<<" != "<<m_nPointsDef);
        GISMO_ASSERT(k>=0 && k<m_nPointsOri,"Point index "<<k<<" out of range [0,"<<m_nPointsOri<<")");

        _getMetricDeformed_impl  <dim>(k,zPhysical,out);
        _getMetricUndeformed_impl<dim>(k,zPhysical,out);

        T ratio;
        T det_ori = out.Gcov_ori.determinant();
        T det_def = out.Gcov_def.determinant();

        if ((det_ori==0 && det_def==0) || (math::isnan(det_ori) && math::isnan(det_def)))
        {
            gsWarn<<"Jacobian determinant is undefined: J^2 = det(Gcov_def) / det(Gcov_ori) = "<<det_def<<"/"<<det_ori<<"! J^2 is set to 1";
            ratio = 1;
        }
        else
            ratio = det_def / det_ori;

        GISMO_ENSURE(ratio >= 0, "Jacobian determinant is negative! det(Gcov_def) = "<<det_def<<"; det(Gcov_ori) = "<<det_ori
                                 <<"\nGcov_def = "<<out.Gcov_def<<"\nGcov_ori = "<<out.Gcov_ori);
        out.J0_sq = ratio;
    }

    /// Number of in-plane points of the last injected map
    index_t nPoints() const { return m_nPointsOri; }

private:

    /**
     * @brief      Shape precondition on an injected map.
     *
     * \ref gsMaterialMatrixBaseDim enforces @c targetDim()==dim in its four
     * CONSTRUCTORS, because it owns the geometry. This class owns none, so that
     * guard has to move to the injection point: a map with the wrong shape (e.g. a
     * 3D->3D volume map handed to \c gsShellKinematics<3,T>) would make
     * <tt>acov = map.jacobian(k)</tt> 3x3 and the subsequent
     * <tt>reshapeCol(k,2,2) = Acov</tt> fail inside Eigen (or corrupt memory in
     * release). Checked once per batch, so the cost is irrelevant.
     */
    static void _checkShape(const gsMapData<T> & map)
    {
        GISMO_ENSURE(map.dim.first==2,
                     "The injected map must have a 2-dimensional parameter domain "
                     "(shell midsurface), but has "<<map.dim.first<<".");
        GISMO_ENSURE(map.dim.second==dim,
                     "The injected map's target dimension ("<<map.dim.second<<") does not "
                     "match the template dimension ("<<dim<<").");
    }

    /**
     * @brief      Number of in-plane points of an injected map.
     *
     * @warning    The count is taken from <tt>map.values[0]</tt> and NEVER from
     *             <tt>map.points</tt>. Under the expression framework
     *             \c gsExprHelper::precompute swaps the quadrature points INTO
     *             \c mapData.points (gsExprHelper.h:451) and straight back OUT
     *             after \c computeMap (:455), so \c map.points is 0 x 0 for every
     *             element of every assembler/evaluator run, while \c values,
     *             \c jacobians, \c deriv2, \c normals and \c patchId are valid.
     *             Sizing the batch on \c points therefore made
     *             \ref computeUndeformed / \ref computeDeformed compute NOTHING
     *             silently, and the next \ref getMetric threw.
     *             There is deliberately NO fallback to \c points: a silent
     *             fallback to the broken path is exactly the failure mode this
     *             guard exists to prevent.
     */
    static index_t _nPoints(const gsMapData<T> & map)
    {
        GISMO_ENSURE(!map.values.empty() && map.values[0].cols()>0,
                     "gsShellKinematics: the map has no values; NEED_VALUE must be requested. "
                     "Note map.points is unusable here -- gsExprHelper::precompute swaps the "
                     "quadrature points back out.");
        return map.values[0].cols();
    }

    /// Implementation of \ref computeDeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricDeformed_impl(const gsMapData<T> & map);

    /// Implementation of \ref computeDeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricDeformed_impl(const gsMapData<T> & map);

    /// Implementation of \ref computeUndeformed for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricUndeformed_impl(const gsMapData<T> & map);

    /// Implementation of \ref computeUndeformed for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricUndeformed_impl(const gsMapData<T> & map);

    /// Deformed part of \ref getMetric for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricDeformed_impl(index_t k, T z, PointMetrics & out) const;

    /// Deformed part of \ref getMetric for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricDeformed_impl(index_t k, T z, PointMetrics & out) const;

    /// Undeformed part of \ref getMetric for planar geometries
    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricUndeformed_impl(index_t k, T z, PointMetrics & out) const;

    /// Undeformed part of \ref getMetric for surface geometries
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricUndeformed_impl(index_t k, T z, PointMetrics & out) const;

private:

    // Batched fundamental forms / basis vectors, one COLUMN per in-plane point.
    // The layout mirrors the *_mat members of gsMaterialMatrixBaseDimData<dim,T>
    // (gsMaterialMatrixBaseDim.h:528-529) one-to-one, so that the transplanted
    // loop bodies below read exactly as they do in gsMaterialMatrixBaseDim.hpp:
    //   A*_mat, B*_mat : 4       x nPoints, reshapeCol(k,2,2)
    //   a*_mat, n*_mat : 2*dim   x nPoints, reshapeCol(k,dim,2)
    //   normal_*_mat   : 3       x nPoints (dim==3 only; the dim==2 normal is e3)
    gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Bcov_ori_mat;
    gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_ncov_ori_mat, m_normal_ori_mat;
    gsMatrix<T> m_Acov_def_mat, m_Acon_def_mat, m_Bcov_def_mat;
    gsMatrix<T> m_acov_def_mat, m_acon_def_mat, m_ncov_def_mat, m_normal_def_mat;

    index_t m_nPointsOri, m_nPointsDef;
};

// =============================================================================
// Batched fundamental forms (transplanted from gsMaterialMatrixBaseDim.hpp)
// =============================================================================

/// @note Transplanted VERBATIM from gsMaterialMatrixBaseDim<3,T>::_computeMetricDeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:470-528. Only the data source changes: the local
///       computeMap on m_defpatches is replaced by the injected \a map.
/// Complexity: O(nPoints) -- a fixed amount of 2x2/3x2 work per point.
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsShellKinematics<dim,T>::_computeMetricDeformed_impl(const gsMapData<T> & map)
{
    GISMO_ENSURE(map.flags & NEED_VALUE   ,"The injected map must be computed with NEED_VALUE");
    GISMO_ENSURE(map.flags & NEED_JACOBIAN,"The injected map must be computed with NEED_JACOBIAN");
    GISMO_ENSURE(map.flags & NEED_DERIV2  ,"The injected map must be computed with NEED_DERIV2");
    GISMO_ENSURE(map.flags & NEED_NORMAL  ,"The injected map must be computed with NEED_NORMAL");
    _checkShape(map);

    // The point count comes from map.values[0], NEVER from map.points: under the
    // expression framework gsExprHelper::precompute swaps the quadrature points
    // INTO mapData.points (gsExprHelper.h:451) and straight back OUT after
    // computeMap (:455), so map.points is 0 x 0 for every element. Sizing on it
    // makes this routine compute NOTHING, silently. See _nPoints().
    const index_t nPts = _nPoints(map);

    gsMatrix<T> deriv2(3,3), mixedB(2,2), acov(3,2), acon(3,2), ncov(3,2), Acov(2,2), Acon(2,2), Bcov(2,2);
    gsMatrix<T> normals;
    gsVector<T> normal;

    normals = map.normals;
    normals.colwise().normalize();
    m_normal_def_mat = normals;

    m_nPointsDef = nPts;

    m_Acov_def_mat.resize(4,nPts);    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,nPts);    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,nPts);    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*3,nPts);  m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*3,nPts);  m_acon_def_mat.setZero();
    m_ncov_def_mat.resize(2*3,nPts);  m_ncov_def_mat.setZero();

    for (index_t k=0; k!= nPts; k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2 = map.deriv2(k).reshaped(3,3);
        normal = normals.col(k);

        Bcov(0,0) = deriv2.row(0).dot(normal);
        Bcov(1,1) = deriv2.row(1).dot(normal);
        Bcov(0,1) = Bcov(1,0) = deriv2.row(2).dot(normal);

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = Acon(i,0)*Bcov(0,j) + Acon(i,1)*Bcov(1,j);

        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);

        // Assign members
        m_acov_def_mat.reshapeCol(k,3,2) = acov;
        m_acon_def_mat.reshapeCol(k,3,2) = acon;
        m_ncov_def_mat.reshapeCol(k,3,2) = ncov;
        m_Acov_def_mat.reshapeCol(k,2,2) = Acov;
        m_Acon_def_mat.reshapeCol(k,2,2) = Acon;
        m_Bcov_def_mat.reshapeCol(k,2,2) = Bcov;
    }
}

/// @note Transplanted VERBATIM from gsMaterialMatrixBaseDim<2,T>::_computeMetricDeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:533-574. Planar: B=0 and ncov=0; map.normals and
///       map.deriv2 are NOT touched.
/// Complexity: O(nPoints).
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsShellKinematics<dim,T>::_computeMetricDeformed_impl(const gsMapData<T> & map)
{
    GISMO_ENSURE(map.flags & NEED_VALUE   ,"The injected map must be computed with NEED_VALUE");
    GISMO_ENSURE(map.flags & NEED_JACOBIAN,"The injected map must be computed with NEED_JACOBIAN");
    _checkShape(map);

    // Sized from map.values[0], NEVER from map.points (0 x 0 after
    // gsExprHelper::precompute, which swaps the points in at gsExprHelper.h:451
    // and back out at :455). See _nPoints().
    const index_t nPts = _nPoints(map);

    gsMatrix<T> acov(2,2), acon(2,2), Acov(2,2), Acon(2,2);

    m_nPointsDef = nPts;

    m_Acov_def_mat.resize(4,nPts);    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,nPts);    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,nPts);    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*2,nPts);  m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*2,nPts);  m_acon_def_mat.setZero();
    m_ncov_def_mat.resize(2*2,nPts);  m_ncov_def_mat.setZero();

    // Kept verbatim from BaseDim (there the setZero above is unreliable through
    // gsThreaded); here it is merely a redundant, harmless explicit zero fill.
    gsMatrix<T> zero(2,2); zero.setZero();
    for (index_t k=0; k!= nPts; k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Assign members
        m_acov_def_mat.reshapeCol(k,2,2) = acov;
        m_acon_def_mat.reshapeCol(k,2,2) = acon;
        m_Acov_def_mat.reshapeCol(k,2,2) = Acov;
        m_Acon_def_mat.reshapeCol(k,2,2) = Acon;

        // Since setZero above does not work
        m_Bcov_def_mat.reshapeCol(k,2,2) = zero;
        m_ncov_def_mat.reshapeCol(k,2,2) = zero;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

/// @note Transplanted VERBATIM from gsMaterialMatrixBaseDim<3,T>::_computeMetricUndeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:588-646.
/// Complexity: O(nPoints).
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsShellKinematics<dim,T>::_computeMetricUndeformed_impl(const gsMapData<T> & map)
{
    GISMO_ENSURE(map.flags & NEED_VALUE   ,"The injected map must be computed with NEED_VALUE");
    GISMO_ENSURE(map.flags & NEED_JACOBIAN,"The injected map must be computed with NEED_JACOBIAN");
    GISMO_ENSURE(map.flags & NEED_DERIV2  ,"The injected map must be computed with NEED_DERIV2");
    GISMO_ENSURE(map.flags & NEED_NORMAL  ,"The injected map must be computed with NEED_NORMAL");
    _checkShape(map);

    // Sized from map.values[0], NEVER from map.points (0 x 0 after
    // gsExprHelper::precompute, which swaps the points in at gsExprHelper.h:451
    // and back out at :455). See _nPoints().
    const index_t nPts = _nPoints(map);

    gsMatrix<T> deriv2(3,3), mixedB(2,2), acov(3,2), acon(3,2), ncov(3,2), Acov(2,2), Acon(2,2), Bcov(2,2);
    gsMatrix<T> normals;
    gsVector<T> normal;

    normals = map.normals;
    normals.colwise().normalize();
    m_normal_ori_mat = normals;

    m_nPointsOri = nPts;

    m_Acov_ori_mat.resize(4,nPts);    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,nPts);    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,nPts);    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*3,nPts);  m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*3,nPts);  m_acon_ori_mat.setZero();
    m_ncov_ori_mat.resize(2*3,nPts);  m_ncov_ori_mat.setZero();

    for (index_t k=0; k!= nPts; k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2 = map.deriv2(k).reshaped(3,3);
        normal = normals.col(k);

        Bcov(0,0) = deriv2.row(0).dot(normal);
        Bcov(1,1) = deriv2.row(1).dot(normal);
        Bcov(0,1) = Bcov(1,0) = deriv2.row(2).dot(normal);

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = Acon(i,0)*Bcov(0,j) + Acon(i,1)*Bcov(1,j);

        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);

        // Assign members
        m_acov_ori_mat.reshapeCol(k,3,2) = acov;
        m_acon_ori_mat.reshapeCol(k,3,2) = acon;
        m_ncov_ori_mat.reshapeCol(k,3,2) = ncov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = Acov;
        m_Acon_ori_mat.reshapeCol(k,2,2) = Acon;
        m_Bcov_ori_mat.reshapeCol(k,2,2) = Bcov;
    }
}

/// @note Transplanted VERBATIM from gsMaterialMatrixBaseDim<2,T>::_computeMetricUndeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:651-691. Planar: B=0 and ncov=0; map.normals and
///       map.deriv2 are NOT touched.
/// Complexity: O(nPoints).
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsShellKinematics<dim,T>::_computeMetricUndeformed_impl(const gsMapData<T> & map)
{
    GISMO_ENSURE(map.flags & NEED_VALUE   ,"The injected map must be computed with NEED_VALUE");
    GISMO_ENSURE(map.flags & NEED_JACOBIAN,"The injected map must be computed with NEED_JACOBIAN");
    _checkShape(map);

    // Sized from map.values[0], NEVER from map.points (0 x 0 after
    // gsExprHelper::precompute, which swaps the points in at gsExprHelper.h:451
    // and back out at :455). See _nPoints().
    const index_t nPts = _nPoints(map);

    gsMatrix<T> acov(2,2), acon(2,2), Acov(2,2), Acon(2,2);

    m_nPointsOri = nPts;

    m_Acov_ori_mat.resize(4,nPts);    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,nPts);    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,nPts);    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*2,nPts);  m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*2,nPts);  m_acon_ori_mat.setZero();
    m_ncov_ori_mat.resize(2*2,nPts);  m_ncov_ori_mat.setZero();

    // See the dim==2 deformed counterpart: kept verbatim from BaseDim.
    gsMatrix<T> zero(2,2); zero.setZero();
    for (index_t k=0; k!= nPts; k++)
    {
        acov = map.jacobian(k);

        Acov = acov.transpose() * acov;
        Acon = Acov.inverse();

        // Construct basis
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = Acon(i,0)*acov.col(0) + Acon(i,1)*acov.col(1);

        // Assign members
        m_acov_ori_mat.reshapeCol(k,2,2) = acov;
        m_acon_ori_mat.reshapeCol(k,2,2) = acon;
        m_Acov_ori_mat.reshapeCol(k,2,2) = Acov;
        m_Acon_ori_mat.reshapeCol(k,2,2) = Acon;

        // Since setZero above does not work
        m_Bcov_ori_mat.reshapeCol(k,2,2) = zero;
        m_ncov_ori_mat.reshapeCol(k,2,2) = zero;
    }
}

// =============================================================================
// Per-point metrics at height z (transplanted from gsMaterialMatrixBaseDim.hpp)
// =============================================================================

/// @note Transplanted from gsMaterialMatrixBaseDim<3,T>::_getMetricDeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:1126-1181, restricted to Gcov_def:
///       G_def(z) = A_def - 2 z B_def + z^2 n_def^T n_def, G_def(2,2) = 1.
///       The deformed BASIS vectors (gcov_def/gcon_def, and hence the extra 3x3
///       inverse for Gcon_def) are NOT part of the consumer contract of
///       \ref PointMetrics and are deliberately not computed; add them to
///       PointMetrics if a consumer ever needs them.
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsShellKinematics<dim,T>::_getMetricDeformed_impl(index_t k, T z, PointMetrics & out) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_def, Bcov_def, ncov_def;

    // Get metric information
    Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    Bcov_def = m_Bcov_def_mat.reshapeCol(k,2,2);
    ncov_def = m_ncov_def_mat.reshapeCol(k,3,2);

    // Compute full metric
    out.Gcov_def.setZero(3,3);
    out.Gcov_def.block(0,0,2,2) = Acov_def - 2.0 * z * Bcov_def + z*z * ncov_def.transpose()*ncov_def;
    out.Gcov_def(2,2) = 1.0;
}

/// @note Transplanted from gsMaterialMatrixBaseDim<2,T>::_getMetricDeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:1186-1242, restricted to Gcov_def: planar, so
///       G_def = A_def (z-independent), G_def(2,2) = 1.
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsShellKinematics<dim,T>::_getMetricDeformed_impl(index_t k, T /*z*/, PointMetrics & out) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");

    gsMatrix<T> Acov_def;

    // Get metric information
    Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);

    // Compute full metric
    out.Gcov_def.setZero(3,3);
    out.Gcov_def.block(0,0,2,2) = Acov_def;
    out.Gcov_def(2,2) = 1.0;
}

//--------------------------------------------------------------------------------------------------------------------------------------

/// @note Transplanted from gsMaterialMatrixBaseDim<3,T>::_getMetricUndeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:1255-1312:
///       G_ori(z) = A_ori - 2 z B_ori + z^2 n_ori^T n_ori, G_ori(2,2) = 1;
///       g_cov = [a_1 + z n_1, a_2 + z n_2, normal];
///       g^c   = sum_i Gcon_ori(c,i) g_i.
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsShellKinematics<dim,T>::_getMetricUndeformed_impl(index_t k, T z, PointMetrics & out) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_normal_ori_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_ori, Bcov_ori, ncov_ori, Gcon_ori(3,3), acov_ori, normal(3,1);

    // Get metric information
    Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    Bcov_ori = m_Bcov_ori_mat.reshapeCol(k,2,2);
    ncov_ori = m_ncov_ori_mat.reshapeCol(k,3,2);

    // Compute full metric
    out.Gcov_ori.setZero(3,3);
    out.Gcov_ori.block(0,0,2,2) = Acov_ori - 2.0 * z * Bcov_ori + z*z * ncov_ori.transpose()*ncov_ori;
    out.Gcov_ori(2,2) = 1.0;
    Gcon_ori = out.Gcov_ori.inverse();

    // Get basis vectors
    acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    normal   = m_normal_ori_mat.reshapeCol(k,3,1);

    // Compute g_cov
    out.gcov_ori.setZero(3,3);
    out.gcov_ori.leftCols(2) = acov_ori + z * ncov_ori;
    out.gcov_ori.col(2) = normal;

    // Compute g_con
    out.gcon_ori.resize(3,3);
    for (index_t c = 0; c!=3; c++)
        out.gcon_ori.col(c) = Gcon_ori(c,0) * out.gcov_ori.col(0) + Gcon_ori(c,1) * out.gcov_ori.col(1) + Gcon_ori(c,2) * out.gcov_ori.col(2);
}

/// @note Transplanted from gsMaterialMatrixBaseDim<2,T>::_getMetricUndeformed_impl,
///       gsMaterialMatrixBaseDim.hpp:1317-1371. Planar: G_ori = A_ori (z-independent),
///       g_cov = [a_1, a_2, e_3] (block(0,0,2,2) = a, so row 2 of the in-plane columns
///       stays zero), g^c = sum_i Gcon_ori(c,i) g_i.
template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsShellKinematics<dim,T>::_getMetricUndeformed_impl(index_t k, T /*z*/, PointMetrics & out) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");

    gsMatrix<T> Acov_ori, Gcon_ori(3,3), acov_ori, normal(3,1);

    // Get metric information
    Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);

    // Compute full metric
    out.Gcov_ori.setZero(3,3);
    out.Gcov_ori.block(0,0,2,2) = Acov_ori;
    out.Gcov_ori(2,2) = 1.0;
    Gcon_ori = out.Gcov_ori.inverse();

    // Get basis vectors
    acov_ori = m_acov_ori_mat.reshapeCol(k,2,2);
    normal << 0,0,1;

    // Compute g_cov
    out.gcov_ori.setZero(3,3);
    out.gcov_ori.block(0,0,2,2) = acov_ori;
    out.gcov_ori.col(2) = normal;

    // Compute g_con
    out.gcon_ori.resize(3,3);
    for (index_t c = 0; c!=3; c++)
        out.gcon_ori.col(c) = Gcon_ori(c,0) * out.gcov_ori.col(0) + Gcon_ori(c,1) * out.gcov_ori.col(1) + Gcon_ori(c,2) * out.gcov_ori.col(2);
}

} // namespace gismo
