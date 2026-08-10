/** @file gsMaterialMatrix3D.h

    @brief Legacy-API adapter wrapping a 3D gsPhaseFieldFracture material law
           behind the classic gsMaterialMatrixBase interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

// gsMaterialMatrixBaseDim.h transitively pulls in gsCore/gsConfig.h, which
// (via gsConfigExt.h) defines gsPhaseFieldFracture_ENABLED. The whole content
// is only meaningful when the PFF module (which provides gsMaterialBase /
// gsMaterialData) is enabled; with PFF disabled this file is an empty TU.
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsPlaneStressCondensation.h>
#include <gsPhaseFieldFracture/materials/gsMaterialBase.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{

// ============================================================================
// Diagnostics counters for the cross-output cache (gsKLShell issue #28,
// constitutive side). One 3D constitutive SWEEP (metrics + triad transform +
// batched plane-stress condensation + back-transform) is executed per unique
// (patch, config-revision, u, z) request; every subsequent gsMaterialMatrix3D
// eval with the SAME request is a cache HIT. The legacy gsThinShellAssembler
// drives SIX gsMaterialMatrixIntegrate coefficients (MatrixA/B/C/D + VectorN/M)
// per element: matrix and vector paths share the same NumGauss z-grid, so all
// six collapse onto ONE sweep + FIVE hits per element.
//
// Declared GISMO_EXPORT here and DEFINED ONCE in gsMaterialMatrix3D_.cpp: with
// -fvisibility=hidden an inline function-local static would give the library
// and each consumer TU SEPARATE hidden instances (consumers would always read
// 0). A single exported symbol guarantees one counter across the .so boundary.
// NOT thread-exact (plain non-atomic globals); single-threaded measurement only.
// ============================================================================
/// Number of 3D constitutive SWEEPS (cache misses) since the last reset.
GISMO_EXPORT size_t gsMaterialMatrix3DSweeps();
/// Resets the sweep counter to zero.
GISMO_EXPORT void   gsMaterialMatrix3DResetSweeps();
/// Internal: increments the sweep counter.
GISMO_EXPORT void   gsMaterialMatrix3DIncrementSweeps();

/// Number of cache HITS (reused sweeps) since the last reset.
GISMO_EXPORT size_t gsMaterialMatrix3DHits();
/// Resets the hit counter to zero.
GISMO_EXPORT void   gsMaterialMatrix3DResetHits();
/// Internal: increments the hit counter.
GISMO_EXPORT void   gsMaterialMatrix3DIncrementHits();

/**
 * @brief   Legacy-API adapter wrapping ONE 3D PFF material law behind the
 *          classic \ref gsMaterialMatrixBase interface.
 *
 * This class lets the UNTOUCHED \ref gsMaterialMatrixIntegrate machinery and
 * \ref gsThinShellAssembler drive a gsPhaseFieldFracture (PFF) 3D constitutive
 * law today. It derives \ref gsMaterialMatrixBaseDim to inherit the shell metric
 * engine (curvilinear kinematics, thickness, same-input guard) and delegates the
 * pointwise constitutive response to \ref gsPlaneStressCondensation, which
 * enforces plane stress (\f$S_{33}=0\f$) on the 3D law.
 *
 * ### Integration mode
 * NotIntegrated: \ref eval3D_matrix / \ref eval3D_vector return the pointwise
 * through-thickness value at each dimensionless \f$z\in[-\tfrac12,\tfrac12]\f$
 * node; the integrator's moment weights build A/B/C/D and N/M from the SAME
 * pointwise data. The material scales \f$z\f$ by the physical thickness itself
 * (\f$z\cdot T\f$), exactly as the legacy nonlinear material does.
 * \ref eval3D_stress is overridden too and returns that same pointwise
 * condensed \f$S^{ij}\f$ (it is the vector output's data, not a moment).
 *
 * ### Cross-output cache (the 6->1 mechanism, issue #28 constitutive side)
 * A per-thread single-entry cache keyed on (patch, config-revision, u, z) makes
 * the legacy 6-integrator pattern cost ONE constitutive sweep + 5 cache hits per
 * element. Correctness rests on bump-on-set semantics only: \ref setDeformed /
 * \ref setUndeformed bump \c m_configRev unconditionally (even re-setting the
 * same pointer with Newton-mutated coefficients), which auto-invalidates the
 * cache; pointers are never compared (task-11 lesson). This complements the
 * task-11 metric same-input guard on the geometric side.
 *
 * ### Out of scope (Phase-4/5 work)
 *  - eval3D_pstress / pstretch / CauchyVector / dmatrix / eval3D_matrix_C: NOT
 *    implemented here; the Base GISMO_NO_IMPLEMENTATION defaults remain.
 *  - No \ref getMaterialMatrix registration: callers construct this adapter
 *    directly (deliberate Phase-3 deviation; registration revisited at parity).
 *  - Does NOT use the standalone \ref gsShellKinematics: deriving BaseDim IS
 *    the delegation for this adapter. The standalone class now exists (it is
 *    the metric engine of \ref gsShellMaterialProvider), but it is a different
 *    contract: it is stateless and driven by a CALLER-INJECTED gsMapData -- no
 *    geometry pointers, no thickness, no configuration state, no setters --
 *    whereas the legacy \ref gsMaterialMatrixBase interface this adapter must
 *    satisfy is exactly setUndeformed / setDeformed on owned geometry plus
 *    patch-based eval, which is what BaseDim provides. Its consumer,
 *    \ref gsShellMaterialProvider, is a gsFunctionSet and not a BaseDim, so it
 *    cannot serve the legacy interface in this adapter's place.
 *  - TFT / Composite interop untested.
 *
 * @tparam dim The shell parametric embedding (2 = planar, 3 = surface)
 * @tparam T   Real type
 *
 * @ingroup KLShell
 */
template <short_t dim, class T>
class gsMaterialMatrix3D : public gsMaterialMatrixBaseDim<dim,T>
{
public:

    GISMO_CLONE_FUNCTION(gsMaterialMatrix3D);

    using Base = gsMaterialMatrixBaseDim<dim,T>;

    typedef typename Base::function_ptr function_ptr;

    /**
     * @brief      Constructor without density.
     *
     * @param[in]  mp         Original (undeformed) geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  law        The 3D PFF material law (non-owning; caller owns it)
     */
    gsMaterialMatrix3D(const gsFunctionSet<T> & mp,
                       const gsFunctionSet<T> & thickness,
                       const gsMaterialBase<T> & law);

    /**
     * @brief      Full constructor with density.
     *
     * @param[in]  mp         Original (undeformed) geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  Density    Density function
     * @param[in]  law        The 3D PFF material law (non-owning; caller owns it)
     */
    gsMaterialMatrix3D(const gsFunctionSet<T> & mp,
                       const gsFunctionSet<T> & thickness,
                       const gsFunctionSet<T> & Density,
                       const gsMaterialBase<T> & law);

    /// Destructor
    virtual ~gsMaterialMatrix3D() {}

public:

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isMatIntegrated() const override { return MatIntegration::NotIntegrated; }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isVecIntegrated() const override { return MatIntegration::NotIntegrated; }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_matrix;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    using Base::eval3D_vector;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;
    using Base::eval3D_stress;

    /// See \ref gsMaterialMatrixBase for details
    std::ostream & print(std::ostream & os) const override;

public:
    /// Shared pointer for gsMaterialMatrix3D
    typedef memory::shared_ptr< gsMaterialMatrix3D > Ptr;
    /// Unique pointer for gsMaterialMatrix3D
    typedef memory::unique_ptr< gsMaterialMatrix3D > uPtr;

protected:

    /**
     * @brief   Cross-output cache entry: the pointwise condensed response
     *          (\a C = 9 x P tangent, \a S = 3 x P stress) of one
     *          (patch, config-revision, u, z) request over P = z.rows()*u.cols()
     *          grid points, colIdx = j*u.cols()+k.
     */
    struct Entry
    {
        index_t     patch = -1;
        uint64_t    rev   = 0;
        gsMatrix<T> u, z, C, S;
    };

    /// Shared engine: fills (or reuses) the per-thread cache for (patch,u,z).
    void _computeBatch(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z) const;

    /// Builds a local orthonormal triad from the UNDEFORMED covariant basis.
    gsMatrix<T> _localCartesianTriad(const gsMatrix<T> & gcov) const;

    using Base::m_data;

    const gsMaterialBase<T>       * m_law;         ///< Non-owning 3D PFF law.
    gsPlaneStressCondensation<T>    m_condensation;///< Plane-stress condenser.
    mutable util::gsThreaded<Entry> m_cache;       ///< Per-thread single-entry cache.

}; // class gsMaterialMatrix3D

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrix3D.hpp)
#endif

#endif // gsPhaseFieldFracture_ENABLED
