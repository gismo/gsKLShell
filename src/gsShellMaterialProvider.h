/** @file gsShellMaterialProvider.h

    @brief Batched, stateless provider of ALL SIX through-thickness moment
           matrices (A/B/C/D/N/M) of a 3D PFF material law, in ONE per-element
           sweep.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

// The whole content is only meaningful when the gsPhaseFieldFracture module
// (which provides gsMaterialBase / gsMaterialData) is enabled. The macro is set
// in the generated gsCore/gsConfigExt.h, pulled in by any gismo header above
// the guard (pattern of gsPlaneStressCondensation.h:17-24).
#include <gsCore/gsLinearAlgebra.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsCore/gsFunction.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsMemory.h>
#include <gsUtils/gsThreaded.h>
#include <gsAssembler/gsGaussRule.h>

#include <gsKLShell/src/gsShellKinematics.h>
#include <gsKLShell/src/gsPlaneStressCondensation.h>

#include <gsPhaseFieldFracture/materials/gsMaterialBase.h>
#include <gsPhaseFieldFracture/materials/gsMaterialContainer.h>

namespace gismo
{

template <short_t dim, class T>
class gsShellMaterialProviderSingle;

/**
 * @brief   Batched, CONFIGURATION-STATELESS provider of the six through-thickness
 *          moment matrices of a Kirchhoff-Love shell whose constitutive response
 *          is a 3D \ref gsMaterialBase (gsPhaseFieldFracture) law under plane
 *          stress.
 *
 * ONCE per element (per thread) this class
 *  1. evaluates the thickness,
 *  2. builds the shell metrics through \ref gsShellKinematics from the TWO
 *     injected geometry maps (undeformed + deformed),
 *  3. maps the covariant strain of every (in-plane point, z-node) pair into a
 *     local Cartesian frame,
 *  4. runs ONE batched \ref gsPlaneStressCondensation over the whole
 *     \f$P = N\cdot n_z\f$ grid,
 *  5. back-transforms and accumulates the z-Gauss moments
 *
 * into a single 42-row per-point column. It replaces the legacy pattern of SIX
 * independent \ref gsMaterialMatrixIntegrate coefficients (MatrixA/B/C/D +
 * VectorN/M), each of which redoes the full geometric precomputation
 * (gsKLShell issue #28).
 *
 * ### Why there is no configuration state
 * The class holds NO deformed geometry, no \c setDeformed, no configuration
 * revision and therefore nothing that could go stale. The deformed
 * configuration arrives EXCLUSIVELY through the map injected at parse time
 * (\ref bindKinematics). That absence of state is the point of this class: the
 * legacy adapter's \c assembleVector regression came from the invalidate/rebuild
 * machinery around \c setDeformed, which simply does not exist here.
 *
 * ### Result layout (42 rows per in-plane point)
 * Each column \a k of \ref cachedMoments is
 * \f$[\,A(9);\,B(9);\,C(9);\,D(9);\,N(3);\,M(3)\,]\f$ with the 3x3 blocks stored
 * COLUMN-MAJOR, so a view can map any block zero-copy:
 * \code
 * gsAsConstMatrix<T>(prov.cachedMoments().data() + k*TargetDim + OffsetA, 3,3);
 * \endcode
 * Offsets: \ref OffsetA = 0, \ref OffsetB = 9, \ref OffsetC = 18,
 * \ref OffsetD = 27, \ref OffsetN = 36, \ref OffsetM = 39.
 *
 * ### B == C
 * The legacy integrator uses moment 0 for A, moment 1 for BOTH B and C and
 * moment 2 for D, and all four read the SAME pointwise tangent
 * (gsMaterialMatrixIntegrate.h:246-260 -- MatrixC is annotated "must be 1" --
 * and gsMaterialMatrixIntegrate.hpp:408-409, where all four MatrixX outputs
 * share one \c _eval3D overload). The condensed tangent is symmetric, so B and C
 * are IDENTICAL pointwise. The moment-1 tangent is therefore accumulated ONCE
 * into rows [9,18) and COPIED into rows [18,27); it is never integrated twice.
 *
 * ### Evaluation protocol (ONE mode: map-bound, deferred)
 * Exactly the mechanism of the solids \ref gsMaterialDataProvider:
 *  - \ref bindKinematics injects the calling thread's TWO geometry-map data
 *    addresses. It MUST re-run on EVERY parse: \c gsExprHelper::cleanUp destroys
 *    and re-allocates the map data, so cached addresses go dangling.
 *  - A piece's \c eval_into (\ref gsShellMaterialProviderSingle) is the
 *    element-boundary signal only: it stores the EXPLICIT quadrature points and
 *    sets \c cacheValid=false. It must NOT read the maps -- the expression
 *    framework may still evaluate them AFTER this provider in the same
 *    per-element loop.
 *  - The FIRST view evaluation on that element calls \ref ensureFilled, which
 *    runs the real work with provably fresh map data; further views reuse it.
 * There is NO unbound fallback: \ref ensureFilled asserts that both slots are
 * bound.
 *
 * ### Two hard constraints on the consuming views
 * 1. Register this provider as a PLAIN variable (\c gsExprHelper::getVar(prov)),
 *    NEVER as a composition (\c getVar(prov,G)). A plain variable lands in
 *    \c m_fdata and is computed at the PARAMETRIC quadrature points
 *    (gsExprHelper.h:458-463), which is what \ref beginElement must store,
 *    because those points are forwarded to \c precomputeParameters. A
 *    composition lands in \c m_cdata and is computed at \c values[0]
 *    (:465-470), i.e. at PHYSICAL coordinates -- invisible for constant material
 *    parameters, wrong for varying ones.
 * 2. Request \c NEED_VALUE ONLY on this function set. \c gsFunction's default
 *    \c deriv_into is a FINITE DIFFERENCE that re-enters \c eval_into at
 *    perturbed points; with \c NEED_DERIV set, every element would call
 *    \ref beginElement several times with perturbed points and thrash the cache.
 *
 * ### Requirements on the two injected maps
 * Both must be computed on the SAME quadrature points, in the same order, with
 * at least
 * \c NEED_VALUE|NEED_JACOBIAN|NEED_DERIV2|NEED_NORMAL for \c dim==3 and
 * \c NEED_VALUE|NEED_JACOBIAN for \c dim==2.
 * \c NEED_VALUE is required on BOTH maps: on the undeformed one it also supplies
 * the physical points of the thickness evaluation, and on BOTH it is what
 * \ref gsShellKinematics sizes its batches from.
 *
 * The reason is the central trap of this file: \c mapData.points is EMPTY (0 x 0)
 * after \c gsExprHelper::precompute, which swaps the quadrature points into it
 * (gsExprHelper.h:451) and straight back out after \c computeMap (:455). Only
 * \c values, \c jacobians, \c deriv2, \c normals and \c patchId survive. Hence
 * the thickness is evaluated at \c mapOri->values[0], the material parameters at
 * the EXPLICITLY stored quadrature points (\ref beginElement), and the metric
 * engine counts points through <tt>map.values[0].cols()</tt>
 * (\c gsShellKinematics::_nPoints). Nothing in this class may read
 * \c mapData.points.
 *
 * ### ori == def aliasing is legal
 * A linear (undeformed) assembly registers the same multipatch twice and
 * \c gsExprHelper::getMap deduplicates, so BOTH bound pointers may be the very
 * same address. The two \ref gsShellKinematics compute calls then fill the
 * undeformed and deformed member sets from identical data, the strain is
 * identically zero and the plane-stress Newton converges at zero strain. Never
 * assume the two addresses differ.
 *
 * ### Thickness vs. parameters: a deliberate asymmetry
 * The thickness is evaluated at PHYSICAL points (\c mapOri->values[0]), exactly
 * as the legacy shell engine does (gsMaterialMatrixBaseDim.hpp:447), while the
 * PFF law's own parameters are evaluated at PARAMETRIC points, exactly as the
 * legacy adapter does (\c precomputeParameters, gsMaterialMatrix3D.hpp:100).
 * The two conventions disagree for spatially varying data and agree for the
 * constant parameters used so far; this is inherited behavior, reproduced on
 * purpose so that a parity test against the legacy path is meaningful.
 *
 * ### Output-request protocol
 * A fill computes ONLY the moments the consuming assembly routine will read (see
 * \ref setRequested, default \ref ShellReq_All). Non-requested rows are ZERO,
 * not stale: the 42 x N buffer is fully zeroed before accumulation. A
 * mask/view mismatch is caught by the views, which assert \c requested(bit)
 * before reading their block. When no matrix moment is requested, the condensed
 * tangent itself is skipped (\c wantC2D=false, see
 * \ref gsPlaneStressCondensation::condense); the Newton loop still evaluates
 * both stress and tangent per sweep -- that is structural, not an oversight.
 *
 * ### Complexity
 * Per element: O(N) geometry work, O(N n_z) frame transforms, and ONE batched
 * condensation over P = N n_z points (whose cost is (#Newton sweeps) x (one
 * batched stress + one batched tangent evaluation)). The legacy six-coefficient
 * pattern pays the geometry work six times.
 *
 * @tparam dim The shell embedding dimension (2 = planar, 3 = surface)
 * @tparam T   Real type
 *
 * @ingroup KLShell
 */
template <short_t dim, class T>
class gsShellMaterialProvider : public gsFunctionSet<T>
{
    typedef typename gsMaterialBase<T>::function_ptr function_ptr;

public:

    GISMO_CLONE_FUNCTION(gsShellMaterialProvider)

    /// Output-request bits (bitwise OR into the request mask, see \ref setRequested)
    enum ShellRequest : unsigned
    {
        ShellReq_A = 1,  ///< Membrane stiffness A (moment 0 of the tangent)
        ShellReq_B = 2,  ///< Coupling stiffness B (moment 1 of the tangent)
        ShellReq_C = 4,  ///< Coupling stiffness C (moment 1 of the tangent; == B)
        ShellReq_D = 8,  ///< Bending stiffness D  (moment 2 of the tangent)
        ShellReq_N = 16, ///< Normal force N       (moment 0 of the stress)
        ShellReq_M = 32, ///< Bending moment M     (moment 1 of the stress)
        /// All four tangent moments; gates the condensed tangent itself
        ShellReq_MatrixMoments = ShellReq_A | ShellReq_B | ShellReq_C | ShellReq_D,
        /// Both stress moments
        ShellReq_StressMoments = ShellReq_N | ShellReq_M,
        /// Everything (default)
        ShellReq_All = ShellReq_MatrixMoments | ShellReq_StressMoments
    };

    /// Row layout of one result column. Unnamed enum (not static const members):
    /// no out-of-line definition is needed when a header-only view odr-uses them.
    enum
    {
        OffsetA   =  0, ///< Row offset of the 3x3 A block (column-major)
        OffsetB   =  9, ///< Row offset of the 3x3 B block (column-major)
        OffsetC   = 18, ///< Row offset of the 3x3 C block (column-major)
        OffsetD   = 27, ///< Row offset of the 3x3 D block (column-major)
        OffsetN   = 36, ///< Row offset of the 3x1 N block
        OffsetM   = 39, ///< Row offset of the 3x1 M block
        TargetDim = 42  ///< Total number of rows: 4*9 + 2*3
    };

    /**
     * @brief   Constructor from a container of (per-patch) 3D material laws.
     *
     * @param[in] materials  Per-patch 3D PFF laws (non-owning; caller owns them)
     * @param[in] undeformed The UNDEFORMED midsurface. The deformed configuration
     *                       is NOT a constructor argument by design: it arrives
     *                       only through \ref bindKinematics.
     * @param[in] thickness  Shell thickness, evaluated at PHYSICAL points, so its
     *                       domain dimension must be \a dim.
     * @param[in] nz         Number of through-thickness Gauss nodes (legacy option
     *                       "NumGauss", default 4).
     */
    gsShellMaterialProvider( const gsMaterialContainer<T> & materials,
                             const gsFunctionSet<T>       & undeformed,
                             const gsFunctionSet<T>       & thickness,
                             index_t                        nz = 4)
    :
    m_materials(materials),
    m_undeformed(memory::make_shared_not_owned(&undeformed)),
    m_thickness(memory::make_shared_not_owned(&thickness)),
    m_nz(nz)
    {
        _init();
    }

    /**
     * @brief   Constructor from a single 3D material law applied to all patches.
     *
     * See the container constructor for the parameter contract.
     */
    gsShellMaterialProvider( const gsMaterialBase<T>      * material,
                             const gsFunctionSet<T>       & undeformed,
                             const gsFunctionSet<T>       & thickness,
                             index_t                        nz = 4)
    :
    m_materials(undeformed.nPieces()),
    m_undeformed(memory::make_shared_not_owned(&undeformed)),
    m_thickness(memory::make_shared_not_owned(&thickness)),
    m_nz(nz)
    {
        GISMO_ENSURE(material!=nullptr,"gsShellMaterialProvider: null material law.");
        for (index_t p = 0; p!=undeformed.nPieces(); ++p)
            m_materials.set(p,const_cast<gsMaterialBase<T>*>(material));
        _init();
    }

    /// Copy constructor. The pieces are OWNED raw pointers, so they must be
    /// rebuilt (a default memberwise copy would share them and double-free);
    /// the per-thread injection slots are deliberately NOT copied -- a copy has
    /// to be re-bound by its own parse.
    gsShellMaterialProvider(const gsShellMaterialProvider & other)
    :
    gsFunctionSet<T>(other),
    m_materials(other.m_materials),
    m_undeformed(other.m_undeformed),
    m_thickness(other.m_thickness),
    m_nz(other.m_nz),
    m_z(other.m_z),
    m_w(other.m_w),
    m_requested(other.m_requested)
    {
        m_mapOri.mine() = nullptr;
        m_mapDef.mine() = nullptr;
        this->_makePieces();
    }

    ~gsShellMaterialProvider()
    {
        freeAll(m_pieces);
    }

    /// Implementation of domainDim, see \ref gsFunctionSet (shell midsurface: 2)
    short_t domainDim() const override { return 2; }

    /// Implementation of targetDim, see \ref gsFunctionSet (the 42 stacked moment rows)
    short_t targetDim() const override { return (short_t)TargetDim; }

    /// Implementation of piece, see \ref gsFunctionSet
    const gsFunction<T> & piece(const index_t p) const override
    {
        GISMO_ASSERT(p>=0 && p<(index_t)m_pieces.size(),"Patch index "<<p<<" out of range.");
        return *m_pieces[p];
    }

    /// Implementation of nPieces, see \ref gsFunctionSet
    index_t nPieces() const override { return (index_t)m_pieces.size(); }

    /// Not evaluated directly; use the per-patch pieces (see \ref piece)
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    { GISMO_UNUSED(u); GISMO_UNUSED(result); GISMO_NO_IMPLEMENTATION; }

    // ------------------------------------------------------------------------
    // Parse-time injection (map-bound mode)
    // ------------------------------------------------------------------------

    /**
     * @brief   Binds the CALLING thread's undeformed and deformed geometry-map data.
     *
     * Called from a view expression's \c parse (once per parse, per thread). The
     * addresses MUST be resolved from the caller's own per-thread geometry-map
     * copies (a \c gsGeometryMap nests by value, so each OpenMP thread's
     * expression copy owns a distinct map whose \c data() resolves to that
     * thread's \c gsMapData slot).
     *
     * @warning Must re-run on EVERY parse: \c gsExprHelper::cleanUp destroys and
     *          re-allocates the map data, invalidating any address cached here.
     *
     * @note \a ori and \a def MAY be the same address (see the class doc on
     *       aliasing); that is the undeformed/linear case, not an error.
     */
    void bindKinematics(const gsMapData<T>* ori, const gsMapData<T>* def) const
    { m_mapOri.mine() = ori; m_mapDef.mine() = def; }

    /// Current thread's injected undeformed map data (nullptr if unbound)
    const gsMapData<T>* mapOri() const { return m_mapOri.mine(); }
    /// Current thread's injected deformed map data (nullptr if unbound)
    const gsMapData<T>* mapDef() const { return m_mapDef.mine(); }

    // ------------------------------------------------------------------------
    // Output requests
    // ------------------------------------------------------------------------

    /**
     * @brief   Sets which of A/B/C/D/N/M the consuming assembly routine will read
     *          (bitwise OR of \ref ShellRequest values).
     *
     * Everything NOT requested is not integrated at all; if no matrix moment is
     * requested the condensed tangent is not even formed. Call this ONCE on the
     * main thread BEFORE assembly starts -- the mask is a PLAIN member (not
     * per-thread state) that the fills only read.
     *
     * Parse-time OR-ing is deliberately NOT used: no view call can distinguish a
     * new parse cycle from a sibling view of the same cycle, so an accumulated
     * mask could never be reset safely.
     */
    void setRequested(unsigned m) const
    {
        // GISMO_ENSURE, not GISMO_ASSERT: the mask is CALLER-supplied on a public
        // API, so these are external contracts and must hold in every build. The
        // cost is nil -- setRequested() runs ONCE per parse, on the main thread.
        GISMO_ENSURE(m!=0,"An empty output-request mask would make every fill a no-op.");
        GISMO_ENSURE((m & ~(unsigned)ShellReq_All)==0,"Unknown bits in the output-request mask.");
        m_requested = m;
    }

    /// True if @a bit (a \ref ShellRequest value) is part of the request mask
    bool requested(unsigned bit) const { return (m_requested & bit)!=0; }

    /// The current output-request mask
    unsigned requestMask() const { return m_requested; }

    // ------------------------------------------------------------------------
    // Element protocol
    // ------------------------------------------------------------------------

    /**
     * @brief   Element-boundary signal: records this element's EXPLICIT quadrature
     *          points, patch and material, and INVALIDATES the cache.
     *
     * Called from a piece's \c eval_into. Does NOT read the injected maps: they
     * may still be recomputed after this provider in the same per-element loop.
     */
    void beginElement(const gsMatrix<T> & points, index_t patch,
                      const gsMaterialBase<T> * material) const
    {
        Cache & c = m_cache.mine();
        c.points     = points;
        c.patch      = patch;
        c.material   = material;
        c.cacheValid = false;
    }

    /**
     * @brief   Deferred fill: if the cache is invalid, run the whole per-element
     *          sweep and store the 42 x N moment block; a no-op otherwise.
     *
     * Called from the FIRST view evaluation of the element, i.e. at a point where
     * all registered map data is provably fresh.
     */
    void ensureFilled() const
    {
        Cache & c = m_cache.mine();
        // Checked BEFORE the early return: a fill reached without a preceding
        // beginElement() has an EMPTY cache, so returning here would silently
        // serve nothing to the views. cacheValid defaults to FALSE precisely so
        // that this cannot pass unnoticed (task-17 lesson, solids provider).
        //
        // GISMO_ENSURE, not GISMO_ASSERT: this is a CALLER-SEQUENCING contract,
        // and under -DNDEBUG an ASSERT would not merely stay silent -- it would
        // vanish, letting ensureFilled() walk into _fill with a null material
        // (a SIGSEGV, observed in the Release unittest run of task 41).
        //
        // COST: one always-false, perfectly predictable pointer compare, run
        // ONCE PER QUADRATURE POINT PER VIEW -- not once per element. Measured:
        // ensureFilled() 20,160 calls vs _fill() 512, i.e. ~39x per fill
        // (callers: var2deriv2dot_expr::eval, flatdot_expr::eval,
        // flatdot2_expr::eval, _eval::quadrature). The cost was measured on
        // that real per-point path and is still negligible: +103,808 Ir
        // (+0.00083%), net whole-program -0.0005%, per-point view kernels
        // bit-identical. Contrast :423 below, which IS once per element.
        GISMO_ENSURE(c.material!=nullptr,
                     "gsShellMaterialProvider: fill without a preceding beginElement(): the "
                     "element-boundary signal (a piece's eval_into) never ran, so this thread's "
                     "cache is empty.");
        if (c.cacheValid) return;

        const gsMapData<T> * mapOri = m_mapOri.mine();
        const gsMapData<T> * mapDef = m_mapDef.mine();
        // GISMO_ENSURE for the same reason as the check above: a missing bind is
        // a caller-sequencing error, and without the guard the next line
        // dereferences a null gsMapData. Runs only on ACTUAL fills (once per
        // element, after the cacheValid early return).
        GISMO_ENSURE(mapOri!=nullptr && mapDef!=nullptr,
                     "gsShellMaterialProvider: the geometry maps are not bound on this thread. "
                     "bindKinematics() must run in EVERY parse (gsExprHelper::cleanUp "
                     "invalidates the addresses).");

        _fill(c,*mapOri,*mapDef);

        c.cacheValid = true;
        ++m_fillCount;
    }

    /// The current element's moment block, 42 x N (valid after \ref ensureFilled).
    /// Column \a k is [A(9);B(9);C(9);D(9);N(3);M(3)], 3x3 blocks column-major.
    const gsMatrix<T> & cachedMoments() const { return m_cache.mine().M; }

    /// Number of in-plane quadrature points of the current element
    index_t cachedPoints() const { return m_cache.mine().M.cols(); }

    // ------------------------------------------------------------------------
    // Observables. Deliberately MEMBER variables with accessors, never function
    // -local statics in exported free functions: under -fvisibility=hidden the
    // library and each consumer TU would get SEPARATE hidden instances (a trap
    // that already cost a repair round in task 10).
    // ------------------------------------------------------------------------

    /// Number of per-element sweeps performed (approximate under OpenMP; plain
    /// counter with benign races -- a positive value proves the fill fired).
    size_t fillCount() const { return m_fillCount; }
    /// Number of sweeps that integrated the tangent moments A/B/C/D
    size_t matrixMomentFills() const { return m_matrixMomentFills; }
    /// Number of sweeps that integrated the stress moments N/M
    size_t stressMomentFills() const { return m_stressMomentFills; }
    /// Resets all three counters (test helper)
    void resetCounters() const
    { m_fillCount = m_matrixMomentFills = m_stressMomentFills = 0; }

    /// Number of through-thickness Gauss nodes
    index_t numGauss() const { return m_nz; }

protected:

    /// Shared constructor body: z-grid + pieces.
    void _init()
    {
        GISMO_ENSURE(m_nz>0,"gsShellMaterialProvider: NumGauss must be positive, got "<<m_nz<<".");
        GISMO_ENSURE(m_undeformed->domainDim()==2,
                     "gsShellMaterialProvider: the midsurface must have a 2-dimensional parameter "
                     "domain, but has "<<m_undeformed->domainDim()<<".");
        GISMO_ENSURE(m_undeformed->targetDim()==dim,
                     "gsShellMaterialProvider: the midsurface target dimension ("
                     <<m_undeformed->targetDim()<<") does not match the template dimension ("<<dim<<").");
        GISMO_ENSURE(m_thickness->domainDim()==dim,
                     "gsShellMaterialProvider: the thickness is evaluated at PHYSICAL points, so its "
                     "domain dimension must be "<<dim<<", but it is "<<m_thickness->domainDim()<<".");

        // z-grid: DIMENSIONLESS nodes in [-1/2,1/2]. The legacy integrator rebuilds
        // this rule for every in-plane point (gsMaterialMatrixIntegrate.hpp:301-311)
        // although the interval is the constant [-1/2,1/2] and the nodes are
        // therefore identical for all points; here it is built ONCE.
        gsGaussRule<T> gauss(m_nz);
        gsMatrix<T> quNodes(1,m_nz);
        gsVector<T> quWeights(m_nz);
        gauss.mapTo((T)(-0.5),(T)(0.5),quNodes,quWeights);
        m_z = quNodes.transpose();
        m_w = quWeights;

        // gsThreaded value-initializes to nullptr under OpenMP; the single-slot
        // (no-OpenMP) case is default-initialized, so set it explicitly here.
        m_mapOri.mine() = nullptr;
        m_mapDef.mine() = nullptr;

        this->_makePieces();
    }

    void _makePieces()
    {
        m_pieces.resize(m_undeformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces[p] = new gsShellMaterialProviderSingle<dim,T>((index_t)p,m_materials.piece((index_t)p),this);
    }

    /**
     * @brief   Local ORTHONORMAL triad from the UNDEFORMED covariant basis.
     *
     * @note Copied VERBATIM from gsMaterialMatrix3D<dim,T>::_localCartesianTriad
     *       (gsMaterialMatrix3D.hpp:50-69).
     *
     * CORRECTNESS: the constitutive plane-stress structure (E13=E23=0, condense
     * along the 3-direction) only holds in a frame whose 3rd axis is the shell
     * NORMAL; the global identity is NOT an admissible Cartesian basis here.
     */
    static void _localCartesianTriad(const gsMatrix<T> & gcov, gsMatrix<T> & triad)
    {
        triad.resize(3,3);
        gsVector<T> e1 = gcov.col(0);
        e1.normalize();
        gsVector<T> e2 = gcov.col(1) - (gcov.col(1).dot(e1)) * e1; // Gram-Schmidt
        e2.normalize();
        gsVector<T> e3 = gcov.col(2);                              // shell normal
        e3.normalize();
        triad.col(0) = e1;
        triad.col(1) = e2;
        triad.col(2) = e3;
    }

protected:

    gsMaterialContainer<T> m_materials;
    function_ptr           m_undeformed;   ///< UNDEFORMED midsurface (no deformed twin: by design)
    function_ptr           m_thickness;    ///< Thickness, evaluated at PHYSICAL points
    index_t                m_nz;           ///< Through-thickness Gauss nodes (legacy "NumGauss")
    gsVector<T>            m_z;            ///< DIMENSIONLESS z nodes in [-1/2,1/2] (size m_nz)
    gsVector<T>            m_w;            ///< z-quadrature weights (size m_nz)

    mutable std::vector<gsShellMaterialProviderSingle<dim,T> *> m_pieces;

    // Per-thread injection slots (unset == nullptr). NEVER store a gsGeometryMap
    // or an expression symbol here: only the resolved per-thread data addresses.
    mutable util::gsThreaded<const gsMapData<T>*> m_mapOri, m_mapDef;

    /// Per-thread element cache. One element belongs to one patch and every
    /// element change goes through beginElement(), so one cache per thread is
    /// enough. All buffers are reused across elements (resized in place).
    struct Cache
    {
        /// false = must (re)fill on the next view eval. Defaults to FALSE: a
        /// freshly constructed cache holds nothing, so claiming validity would
        /// let an ensureFilled() preceding the first beginElement() serve zeros.
        bool cacheValid = false;
        index_t patch = -1;
        const gsMaterialBase<T> * material = nullptr;

        gsMatrix<T> points;   ///< EXPLICIT parametric quadrature points, 2 x N
        gsMatrix<T> M;        ///< The result: 42 x N moment block

        // Scratch, reused across elements (no per-element allocation churn).
        gsMatrix<T>  Tmat;    ///< thickness, 1 x N
        gsMatrix<T>  e2D;     ///< local-Cartesian strain, 3 x P
        gsMatrix<T>  Tback;   ///< per-column R^T, 9 x P (reshapeCol(col,3,3))
        gsMatrix<T>  S2Dc;    ///< condensed Cartesian stress,  3 x P
        gsMatrix<T>  C2Dc;    ///< condensed Cartesian tangent, 9 x P
        gsMatrix<T>  E33;     ///< converged out-of-plane strain, 1 x P
        std::vector<gsMatrix<T>> params; ///< per-parameter 1 x P rows
        gsMaterialData<T>        pdata;  ///< midsurface parameter buffer
        gsShellKinematics<dim,T> kin;    ///< the metric engine (stateless w.r.t. config)
    };
    mutable util::gsThreaded<Cache> m_cache;

    /// The per-element sweep. See the class documentation for the contract.
    /// Declared AFTER \ref Cache: a member declaration is not a complete-class
    /// context, so the type must already be visible here.
    void _fill(Cache & c, const gsMapData<T> & mapOri, const gsMapData<T> & mapDef) const;

    /// Output-request mask: plain member, set once on the main thread before
    /// assembly, read-only during the fills. Default: everything.
    mutable unsigned m_requested = (unsigned)ShellReq_All;

    // Observables (approximate under OpenMP; benign races).
    mutable size_t m_fillCount = 0, m_matrixMomentFills = 0, m_stressMomentFills = 0;

}; // class gsShellMaterialProvider

/**
 * @brief   Per-patch piece of \ref gsShellMaterialProvider (a \ref gsFunction).
 *
 * Its \c eval_into is ONLY the element-boundary signal: it stores the element's
 * quadrature points in the parent's per-thread cache and invalidates it. The
 * returned buffer is deliberately left uninitialised -- the views read the
 * parent cache (\ref gsShellMaterialProvider::cachedMoments), never this result.
 *
 * @tparam dim The shell embedding dimension (2 = planar, 3 = surface)
 * @tparam T   Real type
 *
 * @ingroup KLShell
 */
template <short_t dim, class T>
class gsShellMaterialProviderSingle : public gsFunction<T>
{
public:

    GISMO_CLONE_FUNCTION(gsShellMaterialProviderSingle)

    gsShellMaterialProviderSingle(index_t patch,
                                  const gsMaterialBase<T> * material,
                                  const gsShellMaterialProvider<dim,T> * parent)
    :
    m_pIndex(patch),
    m_material(material),
    m_parent(parent)
    {
        GISMO_ENSURE(m_parent  !=nullptr,"gsShellMaterialProviderSingle: null parent provider.");
        GISMO_ENSURE(m_material!=nullptr,"gsShellMaterialProviderSingle: no material on patch "<<patch<<".");
    }

    /// Shell midsurface parameter domain
    short_t domainDim() const override { return 2; }

    /// The 42 stacked moment rows [A(9);B(9);C(9);D(9);N(3);M(3)]
    short_t targetDim() const override
    { return (short_t)gsShellMaterialProvider<dim,T>::TargetDim; }

    /// Element-boundary signal only -- see the class documentation.
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        m_parent->beginElement(u,m_pIndex,m_material);
        // Resize only, no setZero: in this mode the views read the parent cache,
        // never this result buffer.
        result.resize((index_t)gsShellMaterialProvider<dim,T>::TargetDim,u.cols());
    }

protected:
    index_t m_pIndex;
    const gsMaterialBase<T> * m_material;
    const gsShellMaterialProvider<dim,T> * m_parent;
};

// =============================================================================
// The per-element sweep
// =============================================================================

/**
 * @note The body reproduces gsMaterialMatrix3D<dim,T>::_computeBatch
 *       (gsMaterialMatrix3D.hpp:72-188) step for step -- that sequence is the
 *       task-15 frame-calibrated reference -- with two differences: the metrics
 *       come from the INJECTED maps through \ref gsShellKinematics instead of
 *       gsMaterialMatrixBaseDim, and the pointwise response is integrated into
 *       through-thickness moments instead of being returned pointwise.
 *
 * Complexity: O(N) geometry + O(N*nz) frame transforms + ONE batched
 * condensation over P = N*nz points.
 */
template <short_t dim, class T>
void gsShellMaterialProvider<dim,T>::_fill(Cache & c,
                                           const gsMapData<T> & mapOri,
                                           const gsMapData<T> & mapDef) const
{
    const gsMatrix<T> & u = c.points;     // EXPLICIT parametric quadrature points
    const index_t N  = u.cols();
    const index_t nz = m_nz;
    const index_t P  = N * nz;

    GISMO_ENSURE(mapOri.flags & NEED_VALUE,
                 "gsShellMaterialProvider: the undeformed map must be computed with NEED_VALUE "
                 "(the thickness is evaluated at PHYSICAL points, and gsShellKinematics takes "
                 "its point count from values[0] -- mapData.points is empty here).");
    // The deformed map never carries the thickness, but it too is sized from
    // values[0] by gsShellKinematics; check it HERE so the diagnostic names the
    // deformed map instead of surfacing as a generic kinematics error.
    GISMO_ENSURE(mapDef.flags & NEED_VALUE,
                 "gsShellMaterialProvider: the deformed map must be computed with NEED_VALUE "
                 "(gsShellKinematics takes its point count from values[0] -- mapData.points is "
                 "empty after gsExprHelper::precompute).");
    GISMO_ASSERT(mapOri.values[0].cols()==N,
                 "The injected undeformed map has "<<mapOri.values[0].cols()<<" points but the "
                 "element cached "<<N<<".");
    GISMO_ASSERT(mapDef.values[0].cols()==N,
                 "The injected deformed map has "<<mapDef.values[0].cols()<<" points but the "
                 "element cached "<<N<<".");
    // The metrics come from the bound maps while the thickness and the material
    // parameters are taken on c.patch: a patch mismatch would silently MIX patch
    // data into plausible-looking numbers. Mirrors the solids provider's guard
    // (gsMaterialDataProvider.h:479). The ori==def aliasing case satisfies both.
    GISMO_ASSERT(mapOri.patchId==c.patch,
                 "Injected undeformed map patchId "<<mapOri.patchId<<" does not match the element's "
                 "patch "<<c.patch<<".");
    GISMO_ASSERT(mapDef.patchId==c.patch,
                 "Injected deformed map patchId "<<mapDef.patchId<<" does not match the element's "
                 "patch "<<c.patch<<".");

    // --- (1) Thickness at PHYSICAL points ------------------------------------
    // Legacy convention (gsMaterialMatrixBaseDim.hpp:447). mapData.points is
    // EMPTY (0 x 0) after gsExprHelper::precompute (gsExprHelper.h:451/455);
    // .values are valid, and are the only admissible source here.
    m_thickness->piece(c.patch).eval_into(mapOri.values[0],c.Tmat);
    GISMO_ASSERT(c.Tmat.rows()==1 && c.Tmat.cols()==N,
                 "The thickness must evaluate to 1 x "<<N<<", got "<<c.Tmat.rows()<<" x "<<c.Tmat.cols()<<".");

    // --- (2) Shell metrics from the two injected maps ------------------------
    // ori and def may be the SAME address (see the class doc): the two calls fill
    // disjoint member sets, so aliasing simply yields a zero strain state.
    c.kin.computeUndeformed(mapOri);
    c.kin.computeDeformed  (mapDef);

    // --- (6) Parameters at PARAMETRIC points, replicated per z-row -----------
    // The parameters depend only on the in-plane point; the law evaluates its OWN
    // parameter functions (gsMaterialMatrix3D.hpp:96-109). NOTE the deliberate
    // asymmetry with the thickness above -- documented in the class header.
    c.material->precomputeParameters(c.patch,u,c.pdata);
    const index_t npar = (index_t)c.pdata.parameters.size();
    c.params.resize(npar);
    for (index_t v=0; v!=npar; ++v)
    {
        c.params[v].resize(1,P);
        for (index_t k=0; k!=N; ++k)
            for (index_t j=0; j!=nz; ++j)
                c.params[v](0, j*N+k) = c.pdata.parameters[v](0,k);
    }

    // --- (3)-(5) covariant strain -> local Cartesian; cache the back-transform
    // Column index col = j*N+k, matching both the legacy adapter and the legacy
    // integrator (gsMaterialMatrixIntegrate.hpp:320: vals(i, j*u.cols()+k)).
    c.e2D.resize(3,P);
    c.Tback.resize(9,P);
    gsMatrix<T> ecov(3,1), R(3,3), triad(3,3), Eblk;
    typename gsShellKinematics<dim,T>::PointMetrics pm;
    for (index_t k=0; k!=N; ++k)
    {
        const T t = c.Tmat(0,k);
        for (index_t j=0; j!=nz; ++j)
        {
            const index_t col = j*N+k;
            // getMetric takes the PHYSICAL height: dimensionless z in [-1/2,1/2]
            // times the thickness (gsMaterialMatrix3D.hpp:121).
            c.kin.getMetric(k, m_z(j)*t, pm);

            // (3) Covariant strain, shell Voigt-3 [E11, E22, E01+E10] (ENGINEERING
            // shear), gsMaterialMatrix3D.hpp:125-129.
            Eblk = 0.5 * ( pm.Gcov_def.block(0,0,2,2) - pm.Gcov_ori.block(0,0,2,2) );
            ecov(0,0) = Eblk(0,0);
            ecov(1,0) = Eblk(1,1);
            ecov(2,0) = Eblk(0,1) + Eblk(1,0);

            // (4) Local orthonormal triad, e3 = shell normal.
            _localCartesianTriad(pm.gcov_ori,triad);

            // (5) The engineering-STRAIN transform R, copied VERBATIM from
            // gsMaterialMatrix3D.hpp:160-167.
            //
            // The strain tensor is E = E_ij g^i (x) g^j, so its covariant
            // COMPONENTS live on the CONTRAVARIANT in-plane basis g^1,g^2
            // (= gcon_ori cols 0,1). Projecting onto the orthonormal triad gives
            // e_cart = R e_cov with a_ia = g^i . e_a.
            //
            // TRAP (this is the task-15 bug): do NOT use
            // gsMaterialMatrixBaseDim::_transformation -- that builds the
            // engineering STRESS transform T_sigma (factor 2 in the shear COLUMN).
            // Strain and stress are CONTRAGREDIENT: T_eps = R^{-T}_sigma only for
            // ROTATIONS, not for non-orthonormal curvilinear bases. Using the
            // stress transform on the strain is an O(1) error on a non-orthogonal
            // patch, masked to ~1e-6 on near-orthogonal ones.
            //
            // By energy conjugacy the work-conjugate CONTRAVARIANT stress and
            // tangent come back with R^T: S^ij = R^T S_cart, C^ijkl = R^T C_cart R.
            const T a11 = pm.gcon_ori.col(0).dot(triad.col(0));
            const T a12 = pm.gcon_ori.col(0).dot(triad.col(1));
            const T a21 = pm.gcon_ori.col(1).dot(triad.col(0));
            const T a22 = pm.gcon_ori.col(1).dot(triad.col(1));
            R(0,0) = a11*a11; R(0,1) = a21*a21; R(0,2) = a11*a21;
            R(1,0) = a12*a12; R(1,1) = a22*a22; R(1,2) = a12*a22;
            R(2,0) = 2*a11*a12; R(2,1) = 2*a21*a22; R(2,2) = a11*a22 + a21*a12;

            c.e2D.col(col) = R * ecov;
            c.Tback.reshapeCol(col,3,3) = R.transpose();
        }
    }

    // --- Output-request gating ----------------------------------------------
    const bool needA   = this->requested(ShellReq_A);
    const bool needD   = this->requested(ShellReq_D);
    const bool needN   = this->requested(ShellReq_N);
    const bool needM   = this->requested(ShellReq_M);
    // B and C are the SAME moment-1 integral of the SAME tangent (see the class
    // doc): one accumulation serves both.
    const bool needBC  = this->requested(ShellReq_B) || this->requested(ShellReq_C);
    const bool needMat = (m_requested & (unsigned)ShellReq_MatrixMoments)!=0;
    const bool needStr = needN || needM;

    // --- (7) ONE batched plane-stress condensation over ALL P columns --------
    // The condenser is a pointer plus two scalars, so constructing it here (rather
    // than caching one per patch) is free and keeps this class free of per-material
    // state. wantC2D=false skips the FINAL tangent evaluation + Schur condensation;
    // the Newton loop still needs both stress and tangent per sweep.
    gsPlaneStressCondensation<T> condensation(c.material);
    condensation.condense(c.e2D,c.params,c.S2Dc,c.C2Dc,c.E33,nullptr,nullptr,needMat);

    // --- (8)-(9) back-transform (gated) + through-thickness moments ----------
    // Legacy weight, gsMaterialMatrixIntegrate.hpp:320:
    //   res += w(j,k) * math::pow(z(j,k)*Tmat(0,k), moment) * vals(i, j*N+k) * Tmat(0,k)
    // The multiplication ORDER is kept exactly as written there (((w*pow)*val)*t),
    // so the result is bit-identical to the legacy integrator; do not fold the
    // trailing thickness into the weight prefix.
    // Moments: A=0, B=1, C=1, D=2 (gsMaterialMatrixIntegrate.h:246-260),
    //          N=0, M=1 (:234-244).
    // Rows that are not requested stay ZERO (the buffer is fully zeroed here), so
    // no stale data from a previous element can masquerade as a result.
    c.M.setZero((index_t)TargetDim,N);
    gsMatrix<T> Ccurv(3,3), Scurv(3,1);
    for (index_t k=0; k!=N; ++k)
    {
        const T t = c.Tmat(0,k);
        for (index_t j=0; j!=nz; ++j)
        {
            const index_t col = j*N+k;
            const T zp   = m_z(j)*t;                              // PHYSICAL height
            const T wp0  = m_w(j) * math::pow(zp,(index_t)0);
            const T wp1  = m_w(j) * math::pow(zp,(index_t)1);
            const T wp2  = m_w(j) * math::pow(zp,(index_t)2);
            const gsAsMatrix<T,Dynamic,Dynamic> Tb = c.Tback.reshapeCol(col,3,3);

            if (needMat)
            {
                // C^ijkl = R^T C_cart R, i.e. Tb * C_cart * Tb^T with Tb = R^T
                // (gsMaterialMatrix3D.hpp:184-186).
                const gsAsMatrix<T,Dynamic,Dynamic> Ccart = c.C2Dc.reshapeCol(col,3,3);
                Ccurv = Tb * Ccart * Tb.transpose();
                const T * Cd = Ccurv.data();          // 3x3, column-major
                for (index_t i=0; i!=9; ++i)
                {
                    const T v = Cd[i];
                    if (needA ) c.M(OffsetA+i,k) += wp0 * v * t;
                    if (needBC) c.M(OffsetB+i,k) += wp1 * v * t;
                    if (needD ) c.M(OffsetD+i,k) += wp2 * v * t;
                }
            }
            if (needStr)
            {
                // S^ij = R^T S_cart (gsMaterialMatrix3D.hpp:183).
                Scurv = Tb * c.S2Dc.col(col);
                for (index_t i=0; i!=3; ++i)
                {
                    const T v = Scurv(i,0);
                    if (needN) c.M(OffsetN+i,k) += wp0 * v * t;
                    if (needM) c.M(OffsetM+i,k) += wp1 * v * t;
                }
            }
        }
    }
    // B == C pointwise AND as an integral (same integrand, same moment 1): the
    // moment-1 tangent was accumulated ONCE into rows [9,18) and is copied here
    // (disjoint row ranges, no aliasing).
    // The copy is gated on C alone, and the B rows are cleared again when B was
    // not requested, so the class-doc invariant "non-requested rows are ZERO"
    // holds for EVERY mask -- a C-only mask must not leave the B block populated
    // just because the shared accumulator lives in the B rows.
    if (needBC)
    {
        if ( this->requested(ShellReq_C))
            c.M.middleRows((index_t)OffsetC,9) = c.M.middleRows((index_t)OffsetB,9);
        if (!this->requested(ShellReq_B))
            c.M.middleRows((index_t)OffsetB,9).setZero();
    }

    if (needMat) ++m_matrixMomentFills;
    if (needStr) ++m_stressMomentFills;
}

} // namespace gismo

#endif // gsPhaseFieldFracture_ENABLED
