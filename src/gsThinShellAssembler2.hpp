/** @file gsThinShellAssembler2.hpp

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

#include <gsCore/gsLinearAlgebra.h>

#ifdef gsPhaseFieldFracture_ENABLED

#include <gsKLShell/src/gsThinShellAssembler2.h>
#include <gsKLShell/src/gsThinShellUtils.h>

#include <gsPde/gsBoundaryConditions.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsMSplines/gsMappedBasis.h>

namespace gismo
{

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
gsThinShellAssembler2<d, T, bending>::gsThinShellAssembler2(
                                        const gsMultiPatch<T>         & patches,
                                        const gsMultiBasis<T>         & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunctionSet<T>        & surface_force,
                                        const gsFunctionSet<T>        & thickness,
                                        const gsMaterialContainer<T>  & laws)
:
m_patches(patches),
m_basis(basis),
m_spaceBasis(&m_basis),
m_bcs(bconditions),
m_forceFun(&surface_force),
m_thickFun(&thickness),
m_densityFun(nullptr),
m_materials(laws)
{
    this->_construct(surface_force);
}

template<short_t d, class T, bool bending>
gsThinShellAssembler2<d, T, bending>::gsThinShellAssembler2(
                                        const gsMultiPatch<T>         & patches,
                                        const gsMultiBasis<T>         & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunctionSet<T>        & surface_force,
                                        const gsFunctionSet<T>        & thickness,
                                        const gsMaterialContainer<T>  & laws,
                                        const gsFunctionSet<T>        & density)
:
m_patches(patches),
m_basis(basis),
m_spaceBasis(&m_basis),
m_bcs(bconditions),
m_forceFun(&surface_force),
m_thickFun(&thickness),
m_densityFun(&density),
m_materials(laws)
{
    this->_construct(surface_force);
}

template<short_t d, class T, bool bending>
gsThinShellAssembler2<d, T, bending>::gsThinShellAssembler2(
                                        const gsMultiPatch<T>         & patches,
                                        const gsMultiBasis<T>         & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunctionSet<T>        & surface_force,
                                        const gsFunctionSet<T>        & thickness,
                                        const gsMaterialBase<T>       & law)
:
m_patches(patches),
m_basis(basis),
m_spaceBasis(&m_basis),
m_bcs(bconditions),
m_forceFun(&surface_force),
m_thickFun(&thickness),
m_densityFun(nullptr),
m_materials((index_t)patches.nPatches())
{
    // Non-owning: the caller owns the law (same convention as the provider's
    // single-law constructor, gsShellMaterialProvider.h:249-263).
    for (size_t p = 0; p!=m_patches.nPatches(); ++p)
        m_materials.set((index_t)p,const_cast<gsMaterialBase<T>*>(&law));

    this->_construct(surface_force);
}

template<short_t d, class T, bool bending>
gsThinShellAssembler2<d, T, bending>::gsThinShellAssembler2(
                                        const gsMultiPatch<T>         & patches,
                                        const gsMultiBasis<T>         & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunctionSet<T>        & surface_force,
                                        const gsFunctionSet<T>        & thickness,
                                        const gsMaterialBase<T>       & law,
                                        const gsFunctionSet<T>        & density)
:
m_patches(patches),
m_basis(basis),
m_spaceBasis(&m_basis),
m_bcs(bconditions),
m_forceFun(&surface_force),
m_thickFun(&thickness),
m_densityFun(&density),
m_materials((index_t)patches.nPatches())
{
    for (size_t p = 0; p!=m_patches.nPatches(); ++p)
        m_materials.set((index_t)p,const_cast<gsMaterialBase<T>*>(&law));

    this->_construct(surface_force);
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_construct(const gsFunctionSet<T> & surface_force)
{
    // surface forces defined in the parametric domain (ONLY WORKS FOR 3D)
    m_parametricForce = (surface_force.domainDim()==2 && (d==2||d==3));

    // Pressure and foundation are off by default. DELIBERATE DIVERGENCE from
    // gsThinShellAssembler.hpp:266-269: this is the ONLY place where the two
    // indicators are INITIALISED (setPressure is the only other writer), so that
    // an option change (which re-runs _initialize through _getOptions) cannot
    // silently drop them.
    m_pressFun = nullptr;
    m_foundInd = false;
    m_pressInd = false;

    m_materialFills     = 0;
    m_matrixMomentFills = 0;
    m_stressMomentFills = 0;
    // Legacy leaves m_status indeterminate until the first assembly.
    m_status = ThinShellAssemblerStatus::Success;

    // The thickness is evaluated at PHYSICAL points, so its domain dimension is
    // the EMBEDDING dimension d, not 2. The provider enforces the same thing at
    // its own construction (gsShellMaterialProvider.h:482-484); checking it here
    // reports the mistake at the assembler's construction instead of inside the
    // first assembly.
    GISMO_ENSURE(m_thickFun!=nullptr,"gsThinShellAssembler2: no thickness function was given.");
    GISMO_ENSURE(m_thickFun->domainDim()==d,
                 "gsThinShellAssembler2: the thickness is evaluated at PHYSICAL points, so its "
                 "domain dimension must be "<<d<<", but it is "<<m_thickFun->domainDim()<<".");
    GISMO_ENSURE(m_thickFun->targetDim()==1,
                 "gsThinShellAssembler2: the thickness must be scalar, but has target dimension "
                 <<m_thickFun->targetDim()<<".");

    // The density is registered as a COMPOSITION with the undeformed map in
    // assembleMass, i.e. it is evaluated at PHYSICAL points as well.
    if (m_densityFun!=nullptr)
    {
        GISMO_ENSURE(m_densityFun->domainDim()==d,
                     "gsThinShellAssembler2: the density is evaluated at PHYSICAL points, so its "
                     "domain dimension must be "<<d<<", but it is "<<m_densityFun->domainDim()<<".");
        GISMO_ENSURE(m_densityFun->targetDim()==1,
                     "gsThinShellAssembler2: the density must be scalar, but has target dimension "
                     <<m_densityFun->targetDim()<<".");
    }

    GISMO_ENSURE(m_materials.size()==(index_t)m_patches.nPatches(),
                 "gsThinShellAssembler2: the material container holds "<<m_materials.size()
                 <<" laws but the geometry has "<<m_patches.nPatches()<<" patches.");
    for (index_t p = 0; p!=m_materials.size(); ++p)
        GISMO_ENSURE(m_materials.piece(p)!=nullptr,
                     "gsThinShellAssembler2: no material law on patch "<<p<<".");

    this->_defaultOptions();
    this->_getOptions();
    this->_initialize();
}

// ---------------------------------------------------------------------------
// Options and initialization (gsThinShellAssembler.hpp:195-279)
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_defaultOptions()
{
    m_options.addInt("Continuity","Set the continuity for the space",-1);
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);

    // Assembler options
    gsOptionList assemblerOptions = m_assembler.defaultOptions().wrapIntoGroup("ExprAssembler");
    m_options.update(assemblerOptions,gsOptionList::addIfUnknown);

    m_continuity = -1;
    m_numGauss   = 4;
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_getOptions()
{
    // If the continuity changed, we need to re-initialize the space.
    index_t continuity = m_continuity;
    m_continuity = m_options.getInt("Continuity");
    if (continuity != m_options.getInt("Continuity"))
        this->_initialize();

    // Forwarded to the provider at every assembly (the provider is a local, so
    // there is no stale copy of this number anywhere).
    m_numGauss = m_options.getInt("NumGauss");
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::setOptions(gsOptionList & options)
{
    // Check if the continuity option changed
    // Get old continuity
    index_t continuity = m_options.getInt("Continuity");

    m_options.update(options,gsOptionList::ignoreIfUnknown);

    // If the continuity changed, we need to re-initialize the space.
    if (continuity != m_options.getInt("Continuity"))
        this->_initialize();
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_initialize()
{
    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    GISMO_ASSERT(m_bcs.hasGeoMap(),"No geometry map was assigned to the boundary conditions. Use bc.setGeoMap to assign one!");

    // Set the discretization space
    space m_space = m_assembler.getSpace(*m_spaceBasis, d, 0); // last argument is the space ID

    this->_assembleDirichlet();

    m_ddofs = m_space.fixedPart();
    m_mapper = m_space.mapper();

    // DELIBERATE DIVERGENCE: gsThinShellAssembler.hpp:266-269 resets the
    // foundation and pressure indicators HERE. It must not happen here: this
    // function re-runs on every "Continuity" change, so a load registered before
    // that change would be silently dropped. The two indicators are set in the
    // constructor and nowhere else.

    GISMO_ASSERT(m_forceFun->targetDim()==d,"Force must have " << d<<" dimensions but has "<<m_forceFun->targetDim());
}

// ---------------------------------------------------------------------------
// Boundary conditions and loads (gsThinShellAssembler.hpp:341-368, 1304-1385)
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_assembleDirichlet()
{
    this->_getOptions();
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    m_space.setup(m_bcs, dirichlet::l2Projection, m_continuity);
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::homogenizeDirichlet()
{
    this->_getOptions();
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    m_space.setup(m_bcs, dirichlet::homogeneous, m_continuity);
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_assembleNeumann()
{
    _assembleNeumann_impl<d>();
}

template <short_t d, class T, bool bending>
template <short_t _d>
typename std::enable_if<(_d==3), void>::type
gsThinShellAssembler2<d, T, bending>::_assembleNeumann_impl()
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);
    m_assembler.assembleBdr(m_bcs.get("Neumann"),m_space * g_N * meas(m_ori));
}

template <short_t d, class T, bool bending>
template <short_t _d>
typename std::enable_if<!(_d==3), void>::type
gsThinShellAssembler2<d, T, bending>::_assembleNeumann_impl()
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);
    m_assembler.assembleBdr(m_bcs.get("Neumann"), m_space * g_N * meas(m_ori));
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::_applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts,globalActs;

    space       m_space = m_assembler.trialSpace(0);
    m_mapper = m_space.mapper();

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        GISMO_ASSERT(m_pLoads[i].value.size()==d,"Point load has wrong dimension "<<m_pLoads[i].value.size()<<" instead of "<<d<<"\n");
        GISMO_ASSERT((size_t)m_pLoads[i].patch<m_patches.nPatches(),"Point load is defined on a patch with index "<<m_pLoads[i].patch<<" while the geometry has "<<m_patches.nPatches()<<" patches\n");
        // Compute actives and values of basis functions on point load location.
        if ( m_pLoads[i].parametric )   // in parametric space
        {
            if (const gsMappedBasis<2,T> * mbasis = dynamic_cast<const gsMappedBasis<2,T> * >(m_spaceBasis))
            {
                mbasis->active_into(m_pLoads[i].patch,m_pLoads[i].point, acts );
                mbasis->eval_into  (m_pLoads[i].patch,m_pLoads[i].point, bVals );
            }
            else if (const gsMultiBasis<T> * mbasis = dynamic_cast<const gsMultiBasis<T> * >(m_spaceBasis))
            {
                mbasis->basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts);
                mbasis->basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
            }
            else
                GISMO_ERROR("Basis type not understood");
        }
        else                            // in physical space
        {
            gsMatrix<T> forcePoint;
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);

            if (const gsMappedBasis<2,T> * mbasis = dynamic_cast<const gsMappedBasis<2,T> * >(m_spaceBasis))
            {
                mbasis->active_into(m_pLoads[i].patch,forcePoint, acts );
                mbasis->eval_into  (m_pLoads[i].patch,forcePoint, bVals );
            }
            else if (const gsMultiBasis<T> * mbasis = dynamic_cast<const gsMultiBasis<T> * >(m_spaceBasis))
            {
                mbasis->basis(m_pLoads[i].patch).active_into( forcePoint, acts);
                mbasis->basis(m_pLoads[i].patch).eval_into  ( forcePoint, bVals);
            }
            else
                GISMO_ERROR("Basis type not understood");
        }

        // Add the point load values in the right entries in the global RHS
        for (size_t j = 0; j< d; ++j)
        {
            if (m_pLoads[i].value[j] != 0.0)
            {
                m_mapper.localToGlobal(acts, m_pLoads[i].patch, globalActs,j);
                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    if (m_mapper.is_free_index(globalActs(k,0)))
                        m_rhs(globalActs(k,0), 0) += bVals(k,0) * m_pLoads[i].value[j];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Follower pressure on the DEFORMED configuration
// (gsThinShellAssembler.hpp:413-467). No material is involved, so there is no
// provider and no view here: the four bodies are byte-for-byte the legacy
// algebra. The UNDEFORMED-geometry overload used by the LINEAR routine follows
// below (:370-411); the standalone assemblePressure* entry points, which
// assemble the pressure alone into an otherwise empty system, are NOT ported.
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler2<d, T, bending>::_assemblePressure(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    this->_getOptions();
    _assemblePressure_impl<d,_matrix>(pressFun,deformed);
}

// assembles eq 3.26 from http://resolver.tudelft.nl/uuid:56c0cc91-643d-4817-9702-93fedce5fd78
template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler2<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_pressure = m_assembler.getCoeff(pressFun, m_ori);
    GISMO_ASSERT(pressFun.targetDim()==1,"Pressure function has dimension "<<pressFun.targetDim()<<", but expected 1");

    m_assembler.assemble(
                            -m_pressure.val() * m_space * var1(m_space,m_def).tr()* meas(m_ori)
                        );
}

// assembles eq 3.25 from http://resolver.tudelft.nl/uuid:56c0cc91-643d-4817-9702-93fedce5fd78
template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler2<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_pressure = m_assembler.getCoeff(pressFun, m_ori);
    GISMO_ASSERT(pressFun.targetDim()==1,"Pressure function has dimension "<<pressFun.targetDim()<<", but expected 1");

    // Assemble vector
    m_assembler.assemble(
                  m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
                  );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3), void>::type
gsThinShellAssembler2<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & /*pressFun*/, const gsFunctionSet<T> & /*deformed*/)
{
    // Since pressure works out-of-plane, this function has no effect. ONE body
    // for both values of _matrix, exactly as gsThinShellAssembler.hpp:461-467.
}

// ---------------------------------------------------------------------------
// Follower pressure on the UNDEFORMED configuration, i.e. the LINEAR routine
// (gsThinShellAssembler.hpp:370-411). Same three-body structure as above: one
// EMPTY matrix body, one load body, one shared 2D no-op.
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler2<d, T, bending>::_assemblePressure(const gsFunction<T> & pressFun)
{
    this->_getOptions();
    _assemblePressure_impl<d,_matrix>(pressFun);
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler2<d, T, bending>::_assemblePressure_impl(const gsFunction<T> &)
{
    // No matrix contribution for the linear case
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler2<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & pressFun)
{
    // As everywhere in the linear routine, the "deformed" geometry IS the
    // undeformed one, so m_def and m_ori resolve to the SAME gsMapData.
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(defpatches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_pressure = m_assembler.getCoeff(pressFun, m_ori);
    GISMO_ASSERT(pressFun.targetDim()==1,"Pressure function has dimension "<<pressFun.targetDim()<<", but expected 1");

    m_assembler.assemble(
        m_pressure.val() * m_space * usn(m_def) * meas(m_ori)
        );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3), void>::type
gsThinShellAssembler2<d, T, bending>::_assemblePressure_impl(const gsFunction<T> &)
{
    // Since pressure works out-of-plane, this function has no effect. ONE body
    // for both values of _matrix, exactly as gsThinShellAssembler.hpp:405-411.
}

// ---------------------------------------------------------------------------
// Mass matrix (legacy assembleMass, gsThinShellAssembler.hpp:1605-1698;
// its point-mass helper _applyMass is at :1479-1602)
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler2<d, T, bending>::assembleMass(const bool lumped)
{
    GISMO_ENSURE(m_densityFun!=nullptr,
                 "gsThinShellAssembler2::assembleMass needs a density function; use one of the "
                 "constructor overloads that takes one.");

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);

    // Initialize stystem
    m_assembler.initSystem();

    // NO provider here: the legacy Density output of gsMaterialMatrixIntegrate is
    // exactly thickness(x) * rho(x) evaluated at PHYSICAL points
    // (gsMaterialMatrixBaseDim.hpp:46-60), with no constitutive law involved.
    // Both are therefore registered as COMPOSITIONS with the undeformed map --
    // which is what makes them physical-point evaluations. (The "plain variable,
    // never a composition" rule of gsShellMaterialProvider applies to the
    // PROVIDER only; using plain variables here would move the two functions to
    // the parametric points and break parity with the legacy Density path.)
    auto m_thick   = m_assembler.getCoeff(*m_thickFun  , m_ori);
    auto m_density = m_assembler.getCoeff(*m_densityFun, m_ori);

    space       m_space = m_assembler.trialSpace(0);
    // The mass operator carries no Dirichlet lifting, so the space is set up
    // homogeneously here. This OVERWRITES the space's fixed part, which is
    // shared state: it is restored after the try/catch below.
    m_space.setup(m_bcs, dirichlet::homogeneous, m_continuity);

    try
    {
        // The CONSISTENT mass matrix is assembled UNCONDITIONALLY, also when a
        // lumped matrix is asked for -- see gsThinShellAssembler.hpp for the full
        // reasoning. In short: rowSum() is transparent to the Space trait
        // (gsExpressions/rowsum_expr.h:32), so m_space.rowSum() is a VECTOR-valued
        // expression that gsExprAssembler dispatches into the rhs at compile time;
        // m_assembler.matrix() then returns an unmanaged (empty, zero, or stale)
        // cache instead of a lumped matrix.
        m_assembler.assemble(m_thick.val()*m_density.val()*m_space*m_space.tr()*meas(m_ori));
        m_mass = m_assembler.matrix();

        if (lumped)
        {
            // Row-sum lumping of the ASSEMBLED matrix: the gsProjection pattern,
            // src/gsUtils/gsProjection.hpp:59-69. Since sum_i sum_j M_ij IS the
            // grand sum of M, the total mass is conserved to machine precision on
            // ANY fixture, with or without eliminated Dirichlet dofs.
            //
            // COST, DELIBERATE -- DO NOT "OPTIMISE" THIS AWAY: a full consistent
            // assembly plus one O(nnz) sparse mat-vec, where the vector-valued
            // rowSum() expression plus m_assembler.rhs() would be far cheaper. The
            // cheap route is NOT equivalent: its push<false> path has no column
            // loop, so it also deposits the contributions of the ELIMINATED
            // columns and the resulting diagonal does not sum to the consistent
            // grand sum on a constrained problem.
            gsMatrix<T> ones = gsMatrix<T>::Ones(m_mass.cols(),1);
            gsMatrix<T> rowSums = m_mass * ones;
            gsSparseEntries<T> entries;
            entries.reserve(rowSums.rows());
            for (index_t i = 0; i != rowSums.rows(); ++i)
                entries.add(i,i,rowSums(i,0));
            gsSparseMatrix<T> lumpedMass(m_mass.rows(),m_mass.cols());
            lumpedMass.setFrom(entries);
            m_mass = give(lumpedMass);
        }

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }

    // Undo the homogenization performed above. The trial space is shared state:
    // leaving it homogenized makes every later assemble() on this instance
    // silently return the rhs WITHOUT the Dirichlet lifting. This is exactly what
    // updateBCs(m_bcs) does; it runs on the error path too so that a failed mass
    // assembly does not leave the instance in a homogenized state.
    this->_assembleDirichlet();
    m_ddofs  = m_space.fixedPart();
    m_mapper = m_space.mapper();

    return m_status;
}

// ---------------------------------------------------------------------------
// Linear system (legacy assemble / assemble_impl,
// gsThinShellAssembler.hpp:1731-1832 bending and :1834-1909 non-bending)
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler2<d, T, bending>::assemble()
{
    return assemble_impl<d, bending>();
}

/**
    @brief Assembles the Kirchhoff-Love shell equations including the bending terms.

    The material coefficients are the ONLY difference w.r.t.
    gsThinShellAssembler::assemble_impl: FOUR gsMaterialMatrixIntegrate objects
    (MatrixA/B/C/D), each wrapped in reshape(mm,3,3), are replaced by ONE
    gsShellMaterialProvider and four zero-copy views onto its single per-element
    cache. The views are natively 3x3 and are therefore NOT reshaped.
*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler2<d, T, bending>::assemble_impl()
{
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // Linear assembly: deformed and undeformed geometries are the same. Both
    // maps therefore resolve to the SAME gsMapData (gsExprHelper::getMap
    // deduplicates); the provider supports that aliasing by design (zero strain).
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(defpatches);

    // Initialize stystem
    m_assembler.initSystem();
    m_assembler.initVector(1);

    // ------------------------------------------------------------------------
    // THE structural difference w.r.t. the legacy assembler.
    // The provider is a LOCAL declared BEFORE every expression object below:
    // locals are destroyed in reverse order, so it outlives the views and hence
    // the m_assembler.assemble(...) calls that evaluate them.
    // ------------------------------------------------------------------------
    Provider provider(m_materials,m_patches,*m_thickFun,m_numGauss);
    // FIRST statement after construction: the fill computes only what is asked
    // for. Linear stiffness needs the four tangent moments and no stress moment
    // (this mask is Provider::ShellReq_MatrixMoments, spelled out).
    provider.setRequested( (unsigned)Provider::ShellReq_A | (unsigned)Provider::ShellReq_B
                         | (unsigned)Provider::ShellReq_C | (unsigned)Provider::ShellReq_D );

    // Registered as a PLAIN VARIABLE, never as a composition: a composition is
    // evaluated at PHYSICAL points, which would then be forwarded to
    // precomputeParameters -- which expects PARAMETRIC ones
    // (gsShellMaterialProvider.h:109-117). No NEED_DERIV is ever requested on
    // this symbol either (its deriv_into would be a finite difference).
    auto pv = m_assembler.getCoeff(provider);

    // The four tangent views. They are NATIVELY 3x3 and must NOT be wrapped in
    // reshape (reshape_expr hard-asserts the total size). Both geometry maps
    // handed to them are also parsed by the integrand below (meas(m_ori),
    // jac(m_def)), which is what keeps them bound.
    auto mmA = expr::shellMaterialView<MaterialOutput::MatrixA>(pv,m_ori,m_def,&provider);
    auto mmB = expr::shellMaterialView<MaterialOutput::MatrixB>(pv,m_ori,m_def,&provider);
    auto mmC = expr::shellMaterialView<MaterialOutput::MatrixC>(pv,m_ori,m_def,&provider);
    auto mmD = expr::shellMaterialView<MaterialOutput::MatrixD>(pv,m_ori,m_def,&provider);

    // mult2t is a genuine 9-row plain variable, so this reshape STAYS.
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);

    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) );
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3);

    auto m_N_der    = m_Em_der * mmA + m_Ef_der * mmB;
    auto m_M_der    = m_Em_der * mmC + m_Ef_der * mmD;

    try
    {
        // FIRST in the try block, for BOTH values of _matrix and gated on
        // m_pressInd -- gsThinShellAssembler.hpp:1793-1797 verbatim. The rhs
        // ACCUMULATES, so this ordering is also what keeps the summation order
        // (and hence the last bits) identical to legacy's.
        if (m_pressInd)
        {
            this->_assemblePressure<true>(*m_pressFun);
            this->_assemblePressure<false>(*m_pressFun);
        }

        m_assembler.assemble(
            (
                m_N_der * m_Em_der.tr()
                +
                m_M_der * m_Ef_der.tr()
            ) * meas(m_ori)
            );

        if (m_parametricForce)  // Assemble the force defined in the parameter domain
            m_assembler.assemble(m_space * m_parforce  * meas(m_ori));
        else                    // Assemble the force defined in the physical domain
            m_assembler.assemble(m_space * m_physforce * meas(m_ori));

        this->_assembleNeumann();

        // Assemble the loads
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        // The benchmark observable: exactly one sweep per element (the force and
        // Neumann assemblies re-parse WITHOUT the provider, so they add nothing).
        m_materialFills += provider.fillCount();
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

/**
    @brief Assembles the membrane-only Kirchhoff-Love shell equations (planar
    geometries, or surfaces without bending terms).

    Only the membrane stiffness A is needed, so the provider is asked for
    ShellReq_A alone: the B/C/D moments are not integrated at all.
*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
gsThinShellAssembler2<d, T, bending>::assemble_impl()
{
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // Linear assembly: deformed and undeformed geometries are the same
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(defpatches);

    // Initialize stystem. initVector(1) is a no-op right after initSystem()
    // (gsExprAssembler.h:329-334 vs :449-453); it is spelled out so that both
    // branches share one documented init sequence.
    m_assembler.initSystem();
    m_assembler.initVector(1);

    // ONE provider instead of the single legacy MatrixA integrator; see the
    // bending branch for the lifetime and registration contract.
    Provider provider(m_materials,m_patches,*m_thickFun,m_numGauss);
    provider.setRequested( (unsigned)Provider::ShellReq_A );

    auto pv = m_assembler.getCoeff(provider);

    // Natively 3x3: NO reshape.
    auto mmA = expr::shellMaterialView<MaterialOutput::MatrixA>(pv,m_ori,m_def,&provider);

    space       m_space = m_assembler.trialSpace(0);
    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    auto jacG       = jac(m_def);
    auto m_Em_der   = flat( jacG.tr() * jac(m_space) ) ; //[checked]
    auto m_N_der    = m_Em_der * mmA;

    try
    {
        // See the bending branch: gsThinShellAssembler.hpp:1793-1797 verbatim.
        // For d==2 both instantiations are the shared no-op (a pressure works
        // out-of-plane), so this block only bites for <3,T,false>.
        if (m_pressInd)
        {
            this->_assemblePressure<true>(*m_pressFun);
            this->_assemblePressure<false>(*m_pressFun);
        }

        m_assembler.assemble(
            (
                m_N_der * m_Em_der.tr()
            ) * meas(m_ori)
            );

        if (m_parametricForce)  // Assemble the force defined in the parameter domain
            m_assembler.assemble(m_space * m_parforce  * meas(m_ori));
        else                    // Assemble the force defined in the physical domain
            m_assembler.assemble(m_space * m_physforce * meas(m_ori));

        this->_assembleNeumann();

        // Assemble the loads
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        m_materialFills += provider.fillCount();
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

// ---------------------------------------------------------------------------
// Tangent stiffness on a DEFORMED configuration
// (legacy assembleMatrix / assembleMatrix_impl,
//  gsThinShellAssembler.hpp:1911-1995 bending and :1997-2053 non-bending)
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler2<d, T, bending>::assembleMatrix(const gsFunctionSet<T> & deformed)
{
    return assembleMatrix_impl<d, bending>(deformed);
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler2<d, T, bending>::assembleMatrix(const gsMatrix<T> & solVector)
{
    gsMultiPatch<T> def;
    constructSolution(solVector, def);
    return assembleMatrix(def);
}

/**
    @brief Assembles the tangent stiffness INCLUDING the bending terms.

    Legacy declares SIX gsMaterialMatrixIntegrate coefficients here
    (gsThinShellAssembler.hpp:1933-1938), four of them wrapped in
    reshape(mm,3,3); here ONE gsShellMaterialProvider feeds six zero-copy views
    of its single per-element cache. The views are natively 3x3 (A/B/C/D) and
    3x1 (N/M) and are therefore NOT reshaped.

    Request mask: ALL SIX. The stress moments are not decoration -- N enters the
    geometric membrane term m_Em_der2 and M the geometric bending term
    m_Ef_der2, so a matrix-moments-only mask would trip the view's request
    assert.
*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler2<d, T, bending>::assembleMatrix_impl(const gsFunctionSet<T> & deformed)
{
    // DIVERGENCE (necessary): legacy reaches _getOptions only through
    // homogenizeDirichlet() at gsThinShellAssembler.hpp:1951, i.e. AFTER the
    // material objects are built. The provider needs "NumGauss" at CONSTRUCTION,
    // so the options are read first -- exactly as assemble_impl above does.
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize matrix
    m_assembler.initSystem();
    m_assembler.initMatrix();

    // ------------------------------------------------------------------------
    // THE structural difference w.r.t. the legacy assembler. The provider is a
    // LOCAL declared BEFORE every expression object below: locals are destroyed
    // in reverse order, so it outlives the views and hence the
    // m_assembler.assemble(...) calls that evaluate them.
    // ------------------------------------------------------------------------
    Provider provider(m_materials,m_patches,*m_thickFun,m_numGauss);
    // FIRST statement after construction. ALL SIX: the two stress moments are
    // read by m_Em_der2 / m_Ef_der2 below (this mask is Provider::ShellReq_All).
    provider.setRequested( (unsigned)Provider::ShellReq_All );

    // Registered as a PLAIN VARIABLE, never as a composition: a composition is
    // evaluated at PHYSICAL points, which would then be forwarded to
    // precomputeParameters -- which expects PARAMETRIC ones
    // (gsShellMaterialProvider.h:109-117). No NEED_DERIV is ever requested on
    // this symbol either (its deriv_into would be a finite difference).
    auto pv = m_assembler.getCoeff(provider);

    // The six views. NATIVELY 3x3 (A/B/C/D) and 3x1 (N/M): NO reshape
    // (reshape_expr hard-asserts the total size). Both geometry maps handed to
    // them are also parsed by the integrand below (meas(m_ori), jac(m_def)),
    // which is what keeps them bound.
    auto mmA = expr::shellMaterialView<MaterialOutput::MatrixA>(pv,m_ori,m_def,&provider);
    auto mmB = expr::shellMaterialView<MaterialOutput::MatrixB>(pv,m_ori,m_def,&provider);
    auto mmC = expr::shellMaterialView<MaterialOutput::MatrixC>(pv,m_ori,m_def,&provider);
    auto mmD = expr::shellMaterialView<MaterialOutput::MatrixD>(pv,m_ori,m_def,&provider);
    auto S0  = expr::shellMaterialView<MaterialOutput::VectorN>(pv,m_ori,m_def,&provider);
    auto S1  = expr::shellMaterialView<MaterialOutput::VectorM>(pv,m_ori,m_def,&provider);

    // mult2t is a genuine 9-row plain variable, so this reshape STAYS.
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);

    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Ef_der2  = -(flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
                            + var2deriv2(m_space,m_space,m_def, m_M ));

    auto m_N_der    = m_Em_der * mmA + m_Ef_der * mmB;
    auto m_M_der    = m_Em_der * mmC + m_Ef_der * mmD;

    try
    {
        // Legacy also calls _assembleFoundation<true> here (:1762) and
        // _assembleWeakBCs/_assembleWeakIfc after the assembly (:1778-1779);
        // neither is ported (see the class doc).
        if (m_pressInd) this->_assemblePressure<true>(*m_pressFun,deformed);

        // Assemble matrix
        m_assembler.assemble(
                (
                    m_N_der * m_Em_der.tr()
                    +
                    m_Em_der2
                    +
                    m_M_der * m_Ef_der.tr()
                    +
                    m_Ef_der2
                ) * meas(m_ori)
            );

        // Exactly one sweep per element: the pressure assembly re-parses WITHOUT
        // the provider, so it adds nothing.
        m_materialFills     += provider.fillCount();
        m_matrixMomentFills += provider.matrixMomentFills();
        m_stressMomentFills += provider.stressMomentFills();
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

/**
    @brief Assembles the MEMBRANE tangent stiffness (planar geometries, or
    surfaces without bending terms).

    Request mask: A|N -- the membrane tangent plus the normal force that the
    geometric term m_Em_der2 contracts with. B, C, D and M are not integrated.
*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
gsThinShellAssembler2<d, T, bending>::assembleMatrix_impl(const gsFunctionSet<T> & deformed)
{
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize matrix
    m_assembler.initSystem();
    m_assembler.initMatrix();

    // ONE provider instead of legacy's MatrixA + VectorN integrators; see the
    // bending branch for the lifetime and registration contract.
    Provider provider(m_materials,m_patches,*m_thickFun,m_numGauss);
    provider.setRequested( (unsigned)Provider::ShellReq_A | (unsigned)Provider::ShellReq_N );

    auto pv = m_assembler.getCoeff(provider);

    // Natively 3x3 / 3x1: NO reshape.
    auto mmA = expr::shellMaterialView<MaterialOutput::MatrixA>(pv,m_ori,m_def,&provider);
    auto S0  = expr::shellMaterialView<MaterialOutput::VectorN>(pv,m_ori,m_def,&provider);

    space       m_space = m_assembler.trialSpace(0);

    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_N_der    = m_Em_der * mmA;

    // Assemble matrix
    try
    {
        if (m_pressInd) this->_assemblePressure<true>(*m_pressFun,deformed);

        m_assembler.assemble(
                (
                    m_N_der * m_Em_der.tr()
                    +
                    m_Em_der2
                ) * meas(m_ori)
            );

        m_materialFills     += provider.fillCount();
        m_matrixMomentFills += provider.matrixMomentFills();
        m_stressMomentFills += provider.stressMomentFills();
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...) // add specific cases?
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

// ---------------------------------------------------------------------------
// Residual vector on a DEFORMED configuration
// (legacy assembleVector / assembleVector_impl,
//  gsThinShellAssembler.hpp:2188-2269 bending and :2271-2338 non-bending;
//  the solution-vector overload is at :2654-2660)
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler2<d, T, bending>::assembleVector(const gsFunctionSet<T> & deformed, const bool homogenize)
{
    return assembleVector_impl<d, bending>(deformed,homogenize);
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler2<d, T, bending>::assembleVector(const gsMatrix<T> & solVector, const bool homogenize)
{
    gsMultiPatch<T> def;
    constructSolution(solVector, def);
    return assembleVector(def,homogenize);
}

/**
    @brief Assembles the residual INCLUDING the bending terms.

    Request mask: N|M. NO tangent moment is integrated, and since no matrix
    moment is requested the batched plane-stress condensation does not even form
    the condensed tangent (gsShellMaterialProvider.h:774-783). That saving is the
    measurable point of the per-routine masks; it is observable through
    \ref matrixMomentFills, which this routine leaves untouched.

    Note the init sequence: initVector(1) ONLY -- no initSystem(), no
    initMatrix() -- so the system MATRIX is untouched by this routine.
*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler2<d, T, bending>::assembleVector_impl(const gsFunctionSet<T> & deformed, const bool homogenize)
{
    // See assembleMatrix_impl: the provider needs "NumGauss" at construction.
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize vector
    m_assembler.initVector(1);

    // ONE provider instead of legacy's VectorN + VectorM integrators.
    Provider provider(m_materials,m_patches,*m_thickFun,m_numGauss);
    // FIRST statement after construction: the two stress moments and NOTHING
    // else. The four tangent moments are skipped entirely.
    provider.setRequested( (unsigned)Provider::ShellReq_N | (unsigned)Provider::ShellReq_M );

    auto pv = m_assembler.getCoeff(provider);

    // Natively 3x1: NO reshape.
    auto S0  = expr::shellMaterialView<MaterialOutput::VectorN>(pv,m_ori,m_def,&provider);
    auto S1  = expr::shellMaterialView<MaterialOutput::VectorM>(pv,m_ori,m_def,&provider);

    // mult2t is a genuine 9-row plain variable, so this reshape STAYS.
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m_m2 = m_assembler.getCoeff(mult2t);

    space m_space       = m_assembler.trialSpace(0);
    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    if (homogenize) this->homogenizeDirichlet();
    else            this->_assembleDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    try
    {
        // Legacy also calls _assembleFoundation<false> here (:2027) and
        // _assembleWeakBCs/_assembleWeakIfc at :2044-2045; neither is ported.
        if (m_pressInd) this->_assemblePressure<false>(*m_pressFun,deformed);

        // Assemble vector
        ////// External force
        if (m_parametricForce)  // Assemble the force defined in the parameter domain
            m_assembler.assemble(m_space * m_parforce  * meas(m_ori));
        else                    // Assemble the force defined in the physical domain
            m_assembler.assemble(m_space * m_physforce * meas(m_ori));

        ////// Internal force
        m_assembler.assemble(
            (
                 - ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
            )
            );

        this->_assembleNeumann();

        // Assemble the loads. NOTE (hazard 1 of the class doc): m_rhs is written
        // ONLY when point loads exist -- reproduced from
        // gsThinShellAssembler.hpp:2255-2259 -- so on an UNLOADED problem this
        // routine writes the expression assembler's rhs alone and rhs() returns
        // a stale m_rhs if a loaded assembly ran earlier on this instance.
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        // One sweep per element: the force, pressure and Neumann assemblies
        // re-parse WITHOUT the provider, so they add nothing.
        m_materialFills     += provider.fillCount();
        m_matrixMomentFills += provider.matrixMomentFills();
        m_stressMomentFills += provider.stressMomentFills();
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

/**
    @brief Assembles the MEMBRANE residual (planar geometries, or surfaces
    without bending terms).

    Request mask: N alone -- the single moment this form reads.
*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
gsThinShellAssembler2<d, T, bending>::assembleVector_impl(const gsFunctionSet<T> & deformed, const bool homogenize)
{
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize vector
    m_assembler.initVector(1);

    Provider provider(m_materials,m_patches,*m_thickFun,m_numGauss);
    provider.setRequested( (unsigned)Provider::ShellReq_N );

    auto pv = m_assembler.getCoeff(provider);

    // Natively 3x1: NO reshape.
    auto S0  = expr::shellMaterialView<MaterialOutput::VectorN>(pv,m_ori,m_def,&provider);

    space m_space       = m_assembler.trialSpace(0);
    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    if (homogenize) this->homogenizeDirichlet();
    else            this->_assembleDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    try
    {
        if (m_pressInd) this->_assemblePressure<false>(*m_pressFun,deformed);

        // Assemble vector
        ////// External force
        if (m_parametricForce)  // Assemble the force defined in the parameter domain
            m_assembler.assemble(m_space * m_parforce  * meas(m_ori));
        else                    // Assemble the force defined in the physical domain
            m_assembler.assemble(m_space * m_physforce * meas(m_ori));

        ////// Internal force
        m_assembler.assemble(
            (
                - ( ( m_N * m_Em_der.tr() ) * meas(m_ori) ).tr()
            )
            );

        this->_assembleNeumann();

        // Assemble the loads (see the bending branch on m_rhs;
        // gsThinShellAssembler.hpp:2323-2327)
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        m_materialFills     += provider.fillCount();
        m_matrixMomentFills += provider.matrixMomentFills();
        m_stressMomentFills += provider.stressMomentFills();
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

// ---------------------------------------------------------------------------
// Solution construction (legacy _constructSolution / constructSolution at
// gsThinShellAssembler.hpp:2685-2701, constructMultiPatch at :2831-2886)
// No material is involved in any of these.
// ---------------------------------------------------------------------------

template<short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler2<d, T, bending>::constructMultiPatch(const gsMatrix<T> & solVector) const
{
    m_solvector = solVector;
    space m_space = m_assembler.trialSpace(0);
    m_space.setup(m_bcs, dirichlet::l2Projection, m_continuity);
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart() = m_ddofs;

    if (const gsMappedBasis<2,T> * mbasis = dynamic_cast<const gsMappedBasis<2,T> * >(m_spaceBasis))
    {
        gsMatrix<T> tmp;
        const index_t dim = m_space.dim();
        GISMO_ASSERT(static_cast<size_t>(dim*mbasis->size())==m_mapper.mapSize(),"Something is wrong in the sizes, basis size = "<<mbasis->size()<<" mapper size = "<<m_mapper.mapSize());

        gsMatrix<T> cc(mbasis->size(),d);
        cc.setZero();

        for ( index_t p =0; p!=m_patches.nPieces(); ++p) // Deform the geometry
        {
            for (index_t c = 0; c!=dim; c++) // for all components
            {
                // loop over all basis functions (even the eliminated ones)
                for (size_t i = 0; i < m_mapper.patchSize(p,c); ++i)
                {
                    const index_t ii = m_mapper.index(i, p, c);
                    if ( m_mapper.is_free_index(ii) ) // DoF value is in the solVector
                        cc(i,c) = m_solvector.at(ii);
                    else // eliminated DoF: fill with Dirichlet data
                        cc(i,c) =  m_ddofs.at( m_mapper.global_to_bindex(ii) );
                }
            }
        }
        mbasis->global_coef_to_local_coef(cc,tmp);
        return mbasis->exportToPatches(tmp);
    }
    else
    {
        gsMultiPatch<T> result;
        // Solution vector and solution variable
        solution m_solution = m_assembler.getSolution(m_space, m_solvector);

        gsMatrix<T> cc;
        for ( index_t p =0; p!=m_patches.nPieces(); ++p) // Deform the geometry
        {
            m_solution.extract(cc, p);
            result.addPatch(m_basis.basis(p).makeGeometry( give(cc) ));
        }
        return result;
    }
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler2<d, T, bending>::_constructSolution(const gsMatrix<T> & solVector) const
{
    gsMultiPatch<T> mp = m_patches;
    gsMultiPatch<T> displacement = constructDisplacement(solVector);
    for ( size_t k =0; k!=displacement.nPatches(); ++k) // Deform the geometry
        mp.patch(k).coefs() += displacement.patch(k).coefs();

    return mp;
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler2<d, T, bending>::constructSolution(const gsMatrix<T> & solVector) const
{
    return _constructSolution(solVector);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = _constructSolution(solVector);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::updateMultiPatch(const gsMatrix<T> & solVector, gsMultiPatch<T> & mp) const
{
    mp = _constructSolution(solVector);
}

template<short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler2<d, T, bending>::constructDisplacement(const gsMatrix<T> & solVector) const
{
    return constructMultiPatch(solVector);
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler2<d, T, bending>::constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructDisplacement(solVector);
}

template<short_t d, class T, bool bending>
gsVector<T> gsThinShellAssembler2<d, T, bending>::constructSolutionVector(const gsMultiPatch<T> & displacements) const
{
    gsVector<T> result(m_mapper.freeSize());

    for (size_t p=0; p!=displacements.nPatches(); p++)
    {
        for (size_t dim = 0; dim!=d; dim++)
        {
            for (size_t k=0; k!=m_mapper.patchSize(p,dim); k++)
            {
                if (m_mapper.is_free(k,p,dim))
                    result.at(m_mapper.index(k,p,dim)) = displacements.patch(p).coefs()(k,dim);
            }
        }
    }
    return result;
}

} // namespace gismo

#endif // gsPhaseFieldFracture_ENABLED
