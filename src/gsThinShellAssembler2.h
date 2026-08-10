/** @file gsThinShellAssembler2.h

    @brief Kirchhoff-Love shell assembler built ENTIRELY on the batched
           \ref gsShellMaterialProvider and its zero-copy moment views.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

// The whole content is only meaningful when the gsPhaseFieldFracture module
// (which provides gsMaterialBase / gsMaterialContainer) is enabled. The macro is
// set in the generated gsCore/gsConfigExt.h, pulled in by any gismo header above
// the guard (pattern of gsMaterialMatrix3D_.cpp / gsShellMaterialProvider.h).
#include <gsCore/gsLinearAlgebra.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsShellMaterialProvider.h>
#include <gsKLShell/src/gsShellMaterialExpr.h>

// ThinShellAssemblerStatus (and shell_coupling) live in the legacy header. That
// header is NOT modified by this file in any way; it is included only for the
// status enum, so that a caller can drive the legacy and the new assembler with
// the same status handling. NOTE the include order: gsShellMaterialExpr.h enters
// the gsMaterialMatrixUtils <-> gsMaterialMatrixBase include cycle through Base,
// which is the only order that compiles (see gsShellMaterialExpr.h:28-34).
#include <gsKLShell/src/gsThinShellAssembler.h>

#include <gsPde/gsPointLoads.h>
#include <gsAssembler/gsExprAssembler.h>

namespace gismo
{

/**
 * @brief   Assembles the linear system of a 2D/3D Kirchhoff-Love shell whose
 *          constitutive response is a 3D \ref gsMaterialBase law under plane
 *          stress, using ONE \ref gsShellMaterialProvider per assembly routine
 *          instead of the legacy SIX \c gsMaterialMatrixIntegrate coefficients.
 *
 * This is the parallel assembler of gsKLShell issue #28: the algebra, the option
 * handling, the Dirichlet/Neumann/point-load machinery and the solution
 * construction are reproduced from \ref gsThinShellAssembler VERBATIM, so that a
 * parity suite can compare the two term by term. The legacy assembler is never
 * touched.
 *
 * ### The ONE structural difference
 * Legacy declares six \c gsMaterialMatrixIntegrate<T,MaterialOutput::X>
 * coefficients (MatrixA/B/C/D + VectorN/M), registers each with \c getCoeff and
 * wraps each in \c reshape(mm,3,3). Each of those SIX independently redoes the
 * full per-element geometric precomputation. Here, ONE
 * \ref gsShellMaterialProvider is constructed as a LOCAL in the assembly routine
 * and registered ONCE as a plain variable; the individual moments are read
 * through \ref expr::shellMaterialView expressions that map a block of the
 * provider's single per-element cache ZERO-COPY.
 *
 * The views are natively 3x3 / 3x1 and must therefore NEVER be wrapped in
 * \c reshape (\c reshape_expr hard-asserts the total size). The only surviving
 * \c reshape in the algebra is the one on \c mult2t, which is a genuine 9-row
 * plain variable.
 *
 * ### Conventions of the input functions (three DIFFERENT domains, on purpose)
 *  - the 3D material law's own parameter functions live on the 2D PARAMETRIC
 *    domain (they are evaluated at the parametric quadrature points, see
 *    \ref gsShellMaterialProvider);
 *  - the THICKNESS is evaluated at PHYSICAL points, so \c domainDim()==d
 *    (enforced by the provider);
 *  - the DENSITY is likewise evaluated at PHYSICAL points (it is registered as a
 *    composition with the undeformed map in \ref assembleMass), so
 *    \c domainDim()==d as well (enforced in the constructor).
 *
 * ### Per-routine output-request masks
 * Each assembly routine builds its OWN provider and, as the FIRST statement
 * after the construction, tells it which of the six moments it will read
 * (\c gsShellMaterialProvider::setRequested). Everything else is not integrated
 * at all, and when no tangent moment is requested the plane-stress condensation
 * does not even form the condensed tangent:
 *
 * | routine                       | mask                                   |
 * |-------------------------------|----------------------------------------|
 * | \ref assemble bending         | A\|B\|C\|D                              |
 * | \ref assemble membrane        | A                                       |
 * | \ref assembleMatrix bending   | ALL SIX (N and M enter \c m_Em_der2 / \c m_Ef_der2) |
 * | \ref assembleMatrix membrane  | A\|N                                    |
 * | \ref assembleVector bending   | N\|M                                    |
 * | \ref assembleVector membrane  | N                                       |
 * | \ref assembleMass             | -- (no provider at all)                 |
 *
 * The gating is observable from outside through \ref materialFills (one sweep
 * per element and routine) and, for the two DEFORMED-configuration routines,
 * through \ref matrixMomentFills / \ref stressMomentFills: a
 * \ref assembleVector leaves \ref matrixMomentFills untouched.
 *
 * ### What is NOT ported (later phases of issue #28)
 * Weak boundary conditions and weak interface coupling, the elastic foundation,
 * the standalone pressure entry points (\c assemblePressureMatrix /
 * \c assemblePressureVector, which assemble the pressure ALONE into an
 * otherwise empty system), the DWR assembler, the stress post-processing
 * (\c constructStress / principal stretches), the interface error measures and
 * the (deformed, previous, update) MIP overload. Consequently the \c try blocks
 * of \ref assembleMatrix / \ref assembleVector reproduce the \c try bodies of
 * legacy's \c assembleMatrix_impl and \c assembleVector_impl
 * (\c gsThinShellAssembler.hpp:1966-1988 and \c :2231-2262, in the \c d==3
 * bending instantiations) MINUS
 * \c _assembleFoundation, \c _assembleWeakBCs and \c _assembleWeakIfc.
 *
 * ### DELIBERATE DIVERGENCE from the legacy assembler
 *  1. \ref _initialize does NOT reset the pressure / foundation indicators.
 *     Legacy resets them (\c gsThinShellAssembler.hpp:267-269), and since
 *     \c _getOptions re-runs \c _initialize whenever the "Continuity" option
 *     changes, a \c setPressure issued BEFORE an option change is silently
 *     dropped there. Here the two indicators are initialised in the constructor
 *     and nowhere else.
 *  2. \c m_status is initialised in the constructor (legacy leaves it
 *     indeterminate until the first assembly -- \c gsThinShellAssembler.h:704 has
 *     no default member initialiser, so \c status() before any assembly is still
 *     an uninitialised read there). Task 44 made legacy's \ref assembleMass set
 *     the status explicitly on success, so THAT particular divergence is gone.
 *  3. \ref assembler returns a REFERENCE (legacy returns a copy, and
 *     \c gsExprAssembler's copy constructor does not copy the expression
 *     helper, which makes the copy nearly useless).
 *
 * ### Known inherited hazards (reproduced on purpose, for parity)
 *  1. STALE \c m_rhs once the point loads are gone. \c m_rhs is only written
 *     when point loads are present (\c gsThinShellAssembler2.hpp:704 and
 *     \c gsThinShellAssembler2.hpp:795) and is never cleared. Repeated
 *     \ref assemble calls WITH point
 *     loads are therefore CORRECT: \c m_rhs is ASSIGNED from the expression
 *     assembler's rhs, not accumulated, so every loaded assembly rebuilds it
 *     (measured on the Scordelis-Lo roof: \c ||r||=57540.3 for a point load of
 *     -1e5, and 355467 after changing that same load to -7e5, i.e. it tracks the
 *     new load exactly). The hazard is the MIRROR IMAGE: once any loaded
 *     assembly has run, \ref rhs keeps returning that \c m_rhs for every later
 *     call that does not repopulate it -- after the point loads have been
 *     removed, or after a routine that writes only the expression assembler's
 *     own rhs. \ref assembleVector carries the SAME load-guarded \c m_rhs write
 *     as legacy (the \c if(m_pLoads.numLoads()!=0) blocks of
 *     \c assembleVector_impl, \c gsThinShellAssembler.hpp:2255-2259 for the
 *     bending branch and \c :2323-2327 for the non-bending one), so it
 *     writes only the expression assembler's rhs precisely when there are NO
 *     point loads; \ref assembleMatrix writes neither object (its
 *     \c initSystem() zeroes the expression assembler's rhs and nothing refills
 *     it). A Newton loop on a LOADED problem is therefore consistent, while a
 *     Newton loop that follows a loaded \ref assemble after the loads were
 *     dropped reads the stale linear rhs. Measured: with the loads removed, \ref rhs
 *     still returns 355467 instead of the true load-free 17413.6. \c m_rhs is
 *     protected and has no clear/reset entry point, so the only escape today is
 *     a FRESH assembler instance. This is legacy behaviour
 *     (legacy writes \c m_rhs at exactly the same four load-guarded places,
 *     \c gsThinShellAssembler.hpp:1819-1823, \c :1895-1899, \c :2255-2259 and
 *     \c :2323-2327, and its \c rhs() -- \c gsThinShellAssembler.h:402 -- falls
 *     back to the expression assembler's rhs only while \c m_rhs is still
 *     EMPTY, i.e. never again once a loaded assembly has run); it is
 *     kept so that
 *     the parity suite compares like with like.
 *  2. ~~\ref assembleMass leaves the trial space HOMOGENIZED.~~ FIXED (task 44),
 *     in BOTH assemblers. \ref assembleMass still sets the space up with
 *     \c dirichlet::homogeneous -- the mass operator carries no Dirichlet lifting
 *     -- but it now RESTORES the \c l2Projection setup (through
 *     \c _assembleDirichlet plus the \c m_ddofs / \c m_mapper refresh, i.e. the
 *     body of \ref updateBCs) before returning, on the error path as well. The
 *     sequence \c assemble() -> \c assembleMass() -> \c assemble() therefore
 *     returns the SAME rhs as \c assemble() alone. Before the fix it silently
 *     returned the HOMOGENIZED rhs (measured 17274.8 instead of 1.01007e+08, four
 *     orders of magnitude off, the stiffness being unaffected because the
 *     Dirichlet lifting enters the rhs only). The old remedy -- calling
 *     \ref updateBCs after the mass assembly -- is now redundant but still
 *     harmless, and the tests keep exercising it.
 *  3. ~~\ref assembleMass with \a lumped = true does NOT lump.~~ FIXED (task 44),
 *     in BOTH assemblers. The old code assembled
 *     \c m_space.rowSum() and then read \c m_assembler.matrix(). That is not a
 *     lumping operation and never was: \c rowsum_expr declares
 *     \c Space \c = \c E::Space (\c gsExpressions/rowsum_expr.h:32), so the
 *     expression stays VECTOR-valued and \c gsExprAssembler dispatches it into the
 *     RHS at compile time (\c push<E::isMatrix()>), leaving the system matrix
 *     untouched. \c matrix() then returned an unmanaged CACHE that no reset path
 *     clears -- empty, correctly-sized-but-zero, or the STALE matrix of an earlier
 *     call, depending on the caller's history. The reported "bit-identical to the
 *     consistent matrix, 14608 non-zeros" was the stale case, not a degenerate
 *     lumping. \ref assembleMass now assembles the CONSISTENT matrix
 *     unconditionally and row-sums it when \a lumped is true (the
 *     \c gsUtils/gsProjection.hpp:59-69 pattern), which conserves the grand sum by
 *     construction on any fixture.
 *  4. \ref assembleVector CANNOT BE THE FIRST ASSEMBLY on a fresh instance. It
 *     runs \c initVector(1) only -- which is exactly what preserves a tangent
 *     matrix assembled by a preceding \ref assembleMatrix, so it must stay --
 *     but \c initVector does NOT size the system matrix
 *     (\c gsExprAssembler.h:449-453), while \c gsExprAssembler::assemble asserts
 *     \c m_fmatrix.cols()==numDofs() (\c :1213, first line of \c assemble()). On a fresh instance that assert
 *     is caught by the routine's own \c catch and returns
 *     \c ThinShellAssemblerStatus::AssemblyError. Measured on the Scordelis-Lo
 *     roof: the LEGACY assembler and this one both return \c AssemblyError for a
 *     first-call \c assembleVector, i.e. the behaviour is IDENTICAL and
 *     inherited, not introduced here. Any real Newton loop sizes the system
 *     first (through \ref assemble or \ref assembleMatrix), which is why the
 *     defect has never surfaced.
 *
 * @tparam d        The dimension (2 = planar, 3 = surface)
 * @tparam T        Real type
 * @tparam bending  True: assemble the bending terms; False: membrane only
 *
 * @ingroup KLShell
 */
template <short_t d, class T, bool bending>
class gsThinShellAssembler2
{
public:

    /// The batched six-moment provider that replaces the legacy six integrators
    typedef gsShellMaterialProvider<d,T> Provider;

public:

    /**
     * @brief      Constructor from a single 3D material law (all patches).
     *
     * @param[in]  patches        The undeformed midsurface geometry
     * @param[in]  basis          The discretization basis
     * @param[in]  bconditions    The boundary conditions (needs \c setGeoMap)
     * @param[in]  surface_force  The surface force, \c targetDim()==d
     * @param[in]  thickness      The shell thickness. Evaluated at PHYSICAL
     *                            points, so \c domainDim() must be \a d
     * @param[in]  law            The 3D material law. Its OWN parameter
     *                            functions live on the 2D PARAMETRIC domain
     */
    gsThinShellAssembler2(const gsMultiPatch<T>          & patches,
                          const gsMultiBasis<T>          & basis,
                          const gsBoundaryConditions<T>  & bconditions,
                          const gsFunctionSet<T>         & surface_force,
                          const gsFunctionSet<T>         & thickness,
                          const gsMaterialBase<T>        & law);

    /**
     * @brief      Constructor from a single 3D material law, with a density.
     *
     * @param[in]  density  The mass density. Evaluated at PHYSICAL points
     *                      (a composition with the undeformed map), so
     *                      \c domainDim() must be \a d. Required by
     *                      \ref assembleMass.
     *
     * See the other constructor for the remaining parameters.
     */
    gsThinShellAssembler2(const gsMultiPatch<T>          & patches,
                          const gsMultiBasis<T>          & basis,
                          const gsBoundaryConditions<T>  & bconditions,
                          const gsFunctionSet<T>         & surface_force,
                          const gsFunctionSet<T>         & thickness,
                          const gsMaterialBase<T>        & law,
                          const gsFunctionSet<T>         & density);

    /**
     * @brief      Constructor from a container of per-patch 3D material laws.
     *
     * The container must hold one law per patch. See the single-law constructor
     * for the domain conventions of @a thickness and of the laws' parameters.
     */
    gsThinShellAssembler2(const gsMultiPatch<T>          & patches,
                          const gsMultiBasis<T>          & basis,
                          const gsBoundaryConditions<T>  & bconditions,
                          const gsFunctionSet<T>         & surface_force,
                          const gsFunctionSet<T>         & thickness,
                          const gsMaterialContainer<T>   & laws);

    /**
     * @brief      Constructor from a container of per-patch laws, with a density.
     *
     * See the other constructors for the parameter contract.
     */
    gsThinShellAssembler2(const gsMultiPatch<T>          & patches,
                          const gsMultiBasis<T>          & basis,
                          const gsBoundaryConditions<T>  & bconditions,
                          const gsFunctionSet<T>         & surface_force,
                          const gsFunctionSet<T>         & thickness,
                          const gsMaterialContainer<T>   & laws,
                          const gsFunctionSet<T>         & density);

    //--------------------- OPTIONS ------------------------------------------//

    /// Returns the option list ("ExprAssembler" group + "Continuity" + "NumGauss")
    gsOptionList & options() { return m_options; }

    /// Overwrites the options; re-initializes the space if "Continuity" changed
    void setOptions(gsOptionList & options);

    /// Returns the internal expression assembler (by REFERENCE, see the class doc)
    gsExprAssembler<T> & assembler() { return m_assembler; }

    /// Returns the internal expression assembler (by REFERENCE, see the class doc)
    const gsExprAssembler<T> & assembler() const { return m_assembler; }

    //--------------------- PROBLEM FORMULATION ------------------------------//

    /// Registers a \ref gsPointLoads object for point loads acting on the shell
    void setPointLoads(const gsPointLoads<T> & pLoads) { m_pLoads = pLoads; }

    /// Gets the registered point loads
    const gsPointLoads<T> & getPointLoads() const { return m_pLoads; }

    /**
     * @brief   Registers a FOLLOWER pressure (acting on the DEFORMED normal).
     *
     * @param[in] pressure  Scalar pressure function, \c targetDim()==1,
     *                      evaluated at PHYSICAL points (registered as a
     *                      composition with the undeformed map, exactly as
     *                      legacy does).
     *
     * All three assembly routines use it, exactly as legacy does:
     * \ref assembleMatrix adds the follower stiffness, \ref assembleVector the
     * follower load, and the LINEAR \ref assemble adds the load evaluated on
     * the UNDEFORMED normal (\c usn(m_def) in legacy's
     * \c _assemblePressure_impl<3,!_matrix>(pressFun),
     * \c gsThinShellAssembler.hpp:386-403; the linear routine has NO stiffness
     * contribution, its matrix body -- \c :378-384 -- is empty by design). Has no effect for \a d != 3: a pressure
     * works out-of-plane.
     */
    void setPressure(const gsFunction<T> & pressure)
    { m_pressFun = &pressure; m_pressInd = true; }

    /// Overwrites the boundary conditions
    void updateBCs(const gsBoundaryConditions<T> & bconditions)
    {
        m_bcs = bconditions;
        space m_space = m_assembler.trialSpace(0);
        this->_assembleDirichlet();

        m_ddofs = m_space.fixedPart();
        m_mapper = m_space.mapper();
    }

    /// Sets the Dirichlet BCs to zero
    void homogenizeDirichlet();

    /// Returns the number of degrees of freedom in the assembler
    index_t numDofs() const { return m_assembler.numDofs(); }

    /// Returns the assembler status
    ThinShellAssemblerStatus status() const { return m_status; }

    //--------------------- SYSTEM ASSEMBLY ----------------------------------//

    /// Assembles the LINEAR system matrix and the corresponding right-hand side
    ThinShellAssemblerStatus assemble();

    /// Assembles the mass matrix (thickness x density); if @a lumped is true, a
    /// lumped mass matrix is constructed
    ThinShellAssemblerStatus assembleMass(const bool lumped = false);

    /**
     * @brief   Assembles the TANGENT stiffness matrix on the deformed
     *          configuration @a deformed (the Newton Jacobian).
     *
     * Request mask: all six moments for the bending branch (the stress moments
     * N and M enter the geometric terms \c m_Em_der2 / \c m_Ef_der2), A|N for
     * the membrane branch. The Dirichlet dofs are HOMOGENIZED, the external
     * force, the Neumann data and the point loads are NOT assembled (they carry
     * no dependence on the unknown) -- that is legacy's split between
     * \c assembleMatrix and \ref assembleVector, reproduced verbatim.
     *
     * @warning Writes the system matrix only. Its \c initSystem() ZEROES the
     *          expression assembler's rhs and nothing refills it, so \ref rhs is
     *          meaningless after this call (see hazard 1 in the class doc).
     */
    ThinShellAssemblerStatus assembleMatrix(const gsFunctionSet<T> & deformed);

    /// Same, with the deformed configuration built from @a solVector
    ThinShellAssemblerStatus assembleMatrix(const gsMatrix<T> & solVector);

    /**
     * @brief   Assembles the RESIDUAL vector (external force + Neumann + point
     *          loads - internal force) on the deformed configuration.
     *
     * Request mask: N|M for the bending branch, N for the membrane branch -- no
     * tangent moment is integrated at all, and the plane-stress condensation
     * skips the condensed tangent entirely.
     *
     * @param[in] deformed    The deformed configuration
     * @param[in] homogenize  true (Newton default): homogenize the Dirichlet
     *                        dofs; false: recompute the Dirichlet data
     *
     * @warning Init sequence \c initVector(1) ONLY: the system MATRIX is not
     *          touched (so a tangent from a preceding \ref assembleMatrix
     *          survives) and not SIZED either -- see hazard 4 of the class doc:
     *          this routine cannot be the first assembly on a fresh instance.
     * @warning Writes \c m_rhs only when point loads are present (legacy
     *          verbatim); see hazard 1 of the class doc.
     */
    ThinShellAssemblerStatus assembleVector(const gsFunctionSet<T> & deformed,
                                            const bool homogenize = true);

    /// Same, with the deformed configuration built from @a solVector
    ThinShellAssemblerStatus assembleVector(const gsMatrix<T> & solVector,
                                            const bool homogenize = true);

    //--------------------- SYSTEM ACCESS ------------------------------------//

    /// Returns a reference to the assembled system matrix
    const gsSparseMatrix<T> & matrix() const { return m_assembler.matrix(); }

    /// Returns a reference to the assembled mass matrix
    gsSparseMatrix<T> & massMatrix() { return m_mass; }

    /// Returns a reference to the assembled mass matrix
    const gsSparseMatrix<T> & massMatrix() const { return m_mass; }

    /// Returns a reference to the assembled right-hand side
    const gsMatrix<T> & rhs() const
    { return m_rhs.size()==0 ? m_assembler.rhs() : m_rhs; }

    /**
     * @brief   Accumulated number of per-element material sweeps.
     *
     * The benchmark observable of issue #28: the sum of
     * \c gsShellMaterialProvider::fillCount over every assembly routine that ran
     * so far. A single \ref assemble on a mesh of \a N elements adds exactly
     * \a N (single-threaded; the provider's counter is a plain, non-atomic
     * member, so under OpenMP it is only approximate). \ref assembleMass uses no
     * material and adds nothing.
     */
    size_t materialFills() const { return m_materialFills; }

    /**
     * @brief   Accumulated number of per-element sweeps that integrated the
     *          TANGENT moments A/B/C/D.
     *
     * The direct observable of the per-routine request masks: it grows by one
     * per element in \ref assembleMatrix and stays UNCHANGED across a
     * \ref assembleVector, which is the end-to-end proof that the residual path
     * does no matrix-moment work. Same OpenMP caveat as \ref materialFills.
     *
     * @note Only the DEFORMED-configuration routines (\ref assembleMatrix and
     *       \ref assembleVector) contribute. The LINEAR \ref assemble was left
     *       byte-identical for the parity suite and feeds \ref materialFills
     *       only, so these two counters are NOT a total over all routines.
     */
    size_t matrixMomentFills() const { return m_matrixMomentFills; }

    /// Accumulated number of per-element sweeps that integrated the STRESS
    /// moments N/M (see \ref matrixMomentFills)
    size_t stressMomentFills() const { return m_stressMomentFills; }

    /// Resets all three accumulated material-sweep counters
    void resetMaterialFills()
    { m_materialFills = m_matrixMomentFills = m_stressMomentFills = 0; }

    //--------------------- SOLUTION CONSTRUCTION ----------------------------//

    /// Constructs a multipatch from the free/eliminated coefficients of @a solVector
    gsMultiPatch<T> constructMultiPatch(const gsMatrix<T> & solVector) const;

    /// Overwrites @a mp with the deformed geometry of @a solVector
    void updateMultiPatch(const gsMatrix<T> & solVector, gsMultiPatch<T> & mp) const;

    /// Constructs the DEFORMED geometry from the solution vector
    gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const;

    /// Constructs the DEFORMED geometry from the solution vector into @a deformed
    void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    /// Constructs the DISPLACEMENT field from the solution vector
    gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const;

    /// Constructs the DISPLACEMENT field from the solution vector into @a deformed
    void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    /// The reverse direction: the free coefficients of a displacement multipatch
    gsVector<T> constructSolutionVector(const gsMultiPatch<T> & displacements) const;

protected:

    /// Initializes the space, the Dirichlet dofs and the mapper
    void _initialize();

    /// Registers the default options
    void _defaultOptions();

    /// Reads the options into the members; re-initializes if "Continuity" changed
    void _getOptions();

    /// Computes the Dirichlet dofs by L2 projection
    void _assembleDirichlet();

    /// Assembles the Neumann boundary contribution (rhs only)
    void _assembleNeumann();

    /// Adds the point loads to \c m_rhs
    void _applyLoads();

    /// Assembles the FOLLOWER pressure contribution on the deformed
    /// configuration: the stiffness for @a _matrix, the load otherwise
    template<bool _matrix>
    void _assemblePressure(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    /// Assembles the FOLLOWER pressure contribution for the LINEAR routine,
    /// i.e. on the UNDEFORMED configuration (the deformed geometry IS the
    /// undeformed one there). There is no stiffness contribution in that case,
    /// so the @a _matrix instantiation is an explicit no-op
    template<bool _matrix>
    void _assemblePressure(const gsFunction<T> & pressFun);

    /// Shared constructor body (options, checks, initialization)
    void _construct(const gsFunctionSet<T> & surface_force);

    /// The deformed geometry of @a solVector, on top of the undeformed patches
    gsMultiPatch<T> _constructSolution(const gsMatrix<T> & solVector) const;

private:

    /// Specialisation of assemble() for surfaces (3D) WITH bending
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type assemble_impl();

    /// Specialisation of assemble() for the membrane case (planar, or 3D without bending)
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type assemble_impl();

    /// Implementation of the Neumann contribution for surfaces (3D)
    template<short_t _d>
    typename std::enable_if<(_d==3), void>::type _assembleNeumann_impl();

    /// Implementation of the Neumann contribution for planar geometries (2D)
    template<short_t _d>
    typename std::enable_if<!(_d==3), void>::type _assembleNeumann_impl();

    /// Specialisation of assembleMatrix() for surfaces (3D) WITH bending
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
    assembleMatrix_impl(const gsFunctionSet<T> & deformed);

    /// Specialisation of assembleMatrix() for the membrane case
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
    assembleMatrix_impl(const gsFunctionSet<T> & deformed);

    /// Specialisation of assembleVector() for surfaces (3D) WITH bending
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
    assembleVector_impl(const gsFunctionSet<T> & deformed, const bool homogenize);

    /// Specialisation of assembleVector() for the membrane case
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
    assembleVector_impl(const gsFunctionSet<T> & deformed, const bool homogenize);

    /// Follower-pressure STIFFNESS, surfaces (3D)
    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    /// Follower-pressure LOAD, surfaces (3D)
    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    /// Follower pressure, planar geometries (2D): a NO-OP for both values of
    /// @a _matrix, since a pressure works out-of-plane (as in legacy, ONE body
    /// serves both instantiations)
    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3), void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    /// Follower-pressure STIFFNESS on the UNDEFORMED configuration, surfaces
    /// (3D): EMPTY by design, the linear case has no matrix contribution
    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun);

    /// Follower-pressure LOAD on the UNDEFORMED configuration, surfaces (3D)
    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun);

    /// Follower pressure on the UNDEFORMED configuration, planar geometries
    /// (2D): a NO-OP for both values of @a _matrix (as in legacy, ONE body
    /// serves both instantiations)
    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3), void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun);

protected:

    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space       space;
    typedef typename gsExprAssembler<T>::solution    solution;
    typedef typename gsExprAssembler<T>::variable    variable;

    gsDofMapper m_mapper;

    gsExprAssembler<T> m_assembler;

    gsMultiPatch<T> m_patches;
    mutable gsMultiBasis<T> m_basis;
    const gsFunctionSet<T> * m_spaceBasis;
    gsBoundaryConditions<T> m_bcs;

    mutable gsMatrix<T> m_ddofs;

    gsSparseMatrix<T> m_mass;

    const gsFunctionSet<T> * m_forceFun;
    bool m_parametricForce;

    /// Thickness, evaluated at PHYSICAL points (domainDim()==d)
    const gsFunctionSet<T> * m_thickFun;
    /// Mass density, evaluated at PHYSICAL points (domainDim()==d); may be null
    const gsFunctionSet<T> * m_densityFun;

    /// One 3D material law per patch (non-owning)
    gsMaterialContainer<T> m_materials;

    gsPointLoads<T> m_pLoads;

    mutable gsMatrix<T> m_solvector;

    gsMatrix<T> m_rhs;

    mutable gsOptionList m_options;

    /// The FOLLOWER pressure function (non-owning, null unless \ref setPressure
    /// was called), evaluated at PHYSICAL points
    const gsFunction<T> * m_pressFun;

    /// Pressure / foundation indicators. Initialised in the constructor and
    /// changed only by \ref setPressure: unlike legacy, \ref _initialize does not
    /// reset them (see the class doc). \c m_foundInd stays false -- the elastic
    /// foundation is not ported -- and is kept so that the legacy line numbering
    /// of the assembly routines still reads across.
    mutable bool m_foundInd;
    mutable bool m_pressInd;

    mutable index_t m_continuity;

    /// Number of through-thickness Gauss nodes, forwarded to the provider
    mutable index_t m_numGauss;

    /// Accumulated provider fill count over all assemblies (see \ref materialFills)
    mutable size_t m_materialFills;

    /// Accumulated per-moment-group sweep counts (see \ref matrixMomentFills)
    mutable size_t m_matrixMomentFills;
    mutable size_t m_stressMomentFills;

    mutable ThinShellAssemblerStatus m_status;
};

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssembler2.hpp)
#endif

#endif // gsPhaseFieldFracture_ENABLED
