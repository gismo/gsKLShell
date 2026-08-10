/** @file gsThinShellAssembler.hpp

    @brief Provides linear and nonlinear assemblers for thin shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>

#include <gsPde/gsBoundaryConditions.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsPiecewiseFunction.h>

#include <unordered_set>

namespace gismo
{

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>::gsThinShellAssembler(const gsMultiPatch<T> & patches,
                                                          const gsMultiBasis<T> & basis,
                                                          const gsBoundaryConditions<T> & bconditions,
                                                          const gsFunctionSet<T> & surface_force,
                                                          const gsMaterialMatrixContainer<T> & materialMatrices
                                                          )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_spaceBasis(&m_basis),
                                        m_bcs(bconditions),
                                        m_forceFun(&surface_force),
                                        m_materialMatrices(materialMatrices)
{
    // surface forces defined in the parametric domain (ONLY WORKS FOR 3D)
    m_parametricForce = (surface_force.domainDim()==2 && (d==2||d==3));

    this->_defaultOptions();
    this->_getOptions();
    this->_initialize();
}

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>::gsThinShellAssembler(const gsMultiPatch<T> & patches,
                                                          const gsMultiBasis<T> & basis,
                                                          const gsBoundaryConditions<T> & bconditions,
                                                          const gsFunctionSet<T> & surface_force,
                                                          typename gsMaterialMatrixBase<T>::uPtr & materialMatrix
                                                          )
:
gsThinShellAssembler<d, T, bending>(patches,basis,bconditions,surface_force,materialMatrix.get())
{

}

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>::gsThinShellAssembler(const gsMultiPatch<T> & patches,
                                                          const gsMultiBasis<T> & basis,
                                                          const gsBoundaryConditions<T> & bconditions,
                                                          const gsFunctionSet<T> & surface_force,
                                                          gsMaterialMatrixBase<T> * materialMatrix
                                                          )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_spaceBasis(&basis),
                                        m_bcs(bconditions),
                                        m_forceFun(&surface_force)
{
    m_materialMatrices = gsMaterialMatrixContainer<T>(m_patches.nPatches());
    GISMO_ASSERT(materialMatrix!=nullptr,"Material matrix is incomplete!");
    GISMO_ASSERT(materialMatrix->initialized(),"Material matrix is incomplete!");
    for (size_t p=0; p!=m_patches.nPatches(); p++)
        m_materialMatrices.set(p,materialMatrix);

    // surface forces defined in the parametric domain (ONLY WORKS FOR 3D)
    m_parametricForce = (surface_force.domainDim()==2 && (d==2||d==3));

    this->_defaultOptions();
    this->_getOptions();
    this->_initialize();
}

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>& gsThinShellAssembler<d, T, bending>::operator=( const gsThinShellAssembler& other )
{
    if (this!=&other)
    {
        m_mapper=other.m_mapper;
        m_assembler=other.m_assembler;
        m_evaluator=other.m_evaluator;
        m_patches=other.m_patches;
        m_itpatches=other.m_itpatches;
        m_basis=other.m_basis;
        m_spaceBasis=other.m_spaceBasis;
        m_bcs=other.m_bcs;
        m_ddofs=other.m_ddofs;
        m_mass=other.m_mass;
        m_parametricForce=other.m_parametricForce;
        m_forceFun=other.m_forceFun;
        m_foundFun=other.m_foundFun;
        m_pressFun=other.m_pressFun;
        m_materialMatrices=other.m_materialMatrices;
        m_pLoads=other.m_pLoads;
        m_pMass=other.m_pMass;
        m_solvector=other.m_solvector;
        m_rhs=other.m_rhs;
        m_options=other.m_options;
        m_foundInd=other.m_foundInd;
        m_pressInd=other.m_pressInd;
        m_continuity=other.m_continuity;
        m_alpha_d_bc=other.m_alpha_d_bc;
        m_alpha_r_bc=other.m_alpha_r_bc;
        m_alpha_d_ifc=other.m_alpha_d_ifc;
        m_alpha_r_ifc=other.m_alpha_r_ifc;
        m_IfcDefault=other.m_IfcDefault;
        m_inPlane=other.m_inPlane;
        m_outPlane=other.m_outPlane;
        m_uncoupled=other.m_uncoupled;
        m_strongC0=other.m_strongC0;
        m_weakC0=other.m_weakC0;
        m_strongC1=other.m_strongC1;
        m_weakC1=other.m_weakC1;
        m_unassigned=other.m_unassigned;

        // To do: make copy constructor for the gsExprAssembler
        m_assembler.setIntegrationElements(m_basis);
        GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
        m_assembler.setOptions(m_options.getGroup("ExprAssembler"));
    }
    return *this;
}

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>& gsThinShellAssembler<d, T, bending>::operator=( gsThinShellAssembler&& other )
{
    m_mapper=give(other.m_mapper);
    m_assembler=give(other.m_assembler);
    m_evaluator=give(other.m_evaluator);
    m_patches=give(other.m_patches);
    m_itpatches=give(other.m_itpatches);
    m_basis=give(other.m_basis);
    m_spaceBasis=give(other.m_spaceBasis);
    m_bcs=give(other.m_bcs);
    m_ddofs=give(other.m_ddofs);
    m_mass=give(other.m_mass);
    m_parametricForce=give(other.m_parametricForce);
    m_forceFun=give(other.m_forceFun);
    m_foundFun=give(other.m_foundFun);
    m_pressFun=give(other.m_pressFun);
    m_materialMatrices=give(other.m_materialMatrices);
    m_pLoads=give(other.m_pLoads);
    m_pMass=give(other.m_pMass);
    m_solvector=give(other.m_solvector);
    m_rhs=give(other.m_rhs);
    m_options=give(other.m_options);
    m_foundInd=give(other.m_foundInd);
    m_pressInd=give(other.m_pressInd);
    m_continuity=give(other.m_continuity);
    m_alpha_d_bc=give(other.m_alpha_d_bc);
    m_alpha_r_bc=give(other.m_alpha_r_bc);
    m_alpha_d_ifc=give(other.m_alpha_d_ifc);
    m_alpha_r_ifc=give(other.m_alpha_r_ifc);
    m_IfcDefault=give(other.m_IfcDefault);
    m_inPlane=give(other.m_inPlane);
    m_outPlane=give(other.m_outPlane);
    m_uncoupled=give(other.m_uncoupled);
    m_strongC0=give(other.m_strongC0);
    m_weakC0=give(other.m_weakC0);
    m_strongC1=give(other.m_strongC1);
    m_weakC1=give(other.m_weakC1);
    m_unassigned=give(other.m_unassigned);
    // To do: make copy constructor for the gsExprAssembler
    m_assembler.setIntegrationElements(m_basis);
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));
    return *this;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_defaultOptions()
{
    m_options.addReal("WeakDirichlet","Penalty parameter weak dirichlet conditions",1e3);
    m_options.addReal("WeakClamped","Penalty parameter weak clamped conditions",1e3);
    m_options.addInt("Continuity","Set the continuity for the space",-1);

    m_options.addReal("IfcPenalty","Penalty parameter weak coupling conditions on the interface",1e3);
    m_options.addInt("IfcDefault","Default weak(!) interface coupling; C^k, k={-1,0,1}",1);
    m_options.addString("Solver","Sparse linear solver", "CGDiagonal");

    // Assembler options
    gsOptionList assemblerOptions = m_assembler.defaultOptions().wrapIntoGroup("ExprAssembler");
    m_options.update(assemblerOptions,gsOptionList::addIfUnknown);

    m_continuity = -1;
}


template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_getOptions()
{
    // If the continuity changed, we need to re-initialize the space.
    index_t continuity = m_continuity;
    m_continuity = m_options.getInt("Continuity");
    if (continuity != m_options.getInt("Continuity"))
        this->_initialize();

    m_alpha_d_bc = m_options.getReal("WeakDirichlet");
    m_alpha_r_bc = m_options.getReal("WeakClamped");
    m_alpha_d_ifc = m_alpha_r_ifc = m_options.getReal("IfcPenalty");
    m_IfcDefault = m_options.getInt("IfcDefault");
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::setOptions(gsOptionList & options)
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
void gsThinShellAssembler<d, T, bending>::_initialize()
{
    //gsInfo<<"Active options:\n"<< m_assembler.options() <<"\n";

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    GISMO_ASSERT(m_bcs.hasGeoMap(),"No geometry map was assigned to the boundary conditions. Use bc.setGeoMap to assign one!");

    // Initialize the geometry maps
    // geometryMap m_ori   = m_assembler.getMap(m_patches);
    // geometryMap m_def   = m_assembler.getMap(*m_defpatches);

    // Set the discretization space
    space m_space = m_assembler.getSpace(*m_spaceBasis, d, 0); // last argument is the space ID

    this->_assembleDirichlet();

    m_ddofs = m_space.fixedPart();
    m_mapper = m_space.mapper();

    // foundation is off by default
    m_foundInd = false;
    // pressure is off by default
    m_pressInd = false;

    GISMO_ASSERT(m_forceFun->targetDim()==d,"Force must have " << d<<" dimensions but has "<<m_forceFun->targetDim());

    // test interfaces on in-plane and out-of-plane connection and put them in respective containers
    // _ifcTest();
    // match interfaces where needed
    // todo
    // Put the interfaces in the right container depending on the

}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::addStrongC0(const gsBoxTopology::ifContainer & interfaces)
{
    m_strongC0 = interfaces;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::addStrongC1(const gsBoxTopology::ifContainer & interfaces)
{
    m_strongC1 = interfaces;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::addWeakC0(const gsBoxTopology::ifContainer & interfaces)
{
    m_weakC0 = interfaces;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::addWeakC1(const gsBoxTopology::ifContainer & interfaces)
{
    m_weakC1 = interfaces;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::addUncoupled(const gsBoxTopology::ifContainer & interfaces)
{
    m_uncoupled = interfaces;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::initInterfaces()
{
    this->_getOptions();
    // Find unassigned interfaces and add them to the right containers
    for (gsBoxTopology::const_iiterator it = m_patches.topology().iBegin(); it!=m_patches.topology().iEnd(); it++)
    {
        if (
                std::find(m_strongC0.begin(), m_strongC0.end(), *it) == m_strongC0.end() // m_strongC0 does not contain *it
            &&  std::find(m_strongC1.begin(), m_strongC1.end(), *it) == m_strongC1.end() // m_strongC1 does not contain *it
            &&  std::find(m_weakC0.begin(), m_weakC0.end(), *it) == m_weakC0.end() // m_weakC0 does not contain *it
            &&  std::find(m_weakC1.begin(), m_weakC1.end(), *it) == m_weakC1.end() // m_weakC1 does not contain *it
            &&  std::find(m_uncoupled.begin(), m_uncoupled.end(), *it) == m_uncoupled.end() // m_uncoupled does not contain *it
                )
        {
            if (m_IfcDefault==-1)
                continue;
            else if (m_IfcDefault==0)
                m_weakC0.push_back(*it);
            else if (m_IfcDefault==1)
                m_weakC1.push_back(*it);
            else
                GISMO_ERROR("Option unknown");
        }
    }

    // Set strong C0 using the setup function.
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_assembleNeumann()
{
    _assembleNeumann_impl<d>();
}

template <short_t d, class T, bool bending>
template <short_t _d>
typename std::enable_if<(_d==3), void>::type
gsThinShellAssembler<d, T, bending>::_assembleNeumann_impl()
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);
    m_assembler.assembleBdr(m_bcs.get("Neumann"),m_space * g_N * meas(m_ori));
}

template <short_t d, class T, bool bending>
template <short_t _d>
typename std::enable_if<!(_d==3), void>::type
gsThinShellAssembler<d, T, bending>::_assembleNeumann_impl()
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);
    m_assembler.assembleBdr(m_bcs.get("Neumann"), m_space * g_N * meas(m_ori));
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assemblePressure(const gsFunction<T> & pressFun)
{
    this->_getOptions();
    _assemblePressure_impl<d,_matrix>(pressFun);
}

template <short_t d, class T, bool bending>
template <short_t _d, bool matrix>
typename std::enable_if<(_d==3) && matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> &)
{
    // No matrix contribution for the linear case
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & pressFun)
{
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
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> &)
{
    // Since pressure works out-of-plane, this function has no effect
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assemblePressure(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    this->_getOptions();
    _assemblePressure_impl<d,_matrix>(pressFun,deformed);
}

// assembles eq 3.26 from http://resolver.tudelft.nl/uuid:56c0cc91-643d-4817-9702-93fedce5fd78
template <short_t d, class T, bool bending>
template <short_t _d, bool matrix>
typename std::enable_if<(_d==3) && matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_pressure = m_assembler.getCoeff(pressFun, m_ori);
    GISMO_ASSERT(pressFun.targetDim()==1,"Pressure function has dimension "<<pressFun.targetDim()<<", but expected 1");

    m_assembler.assemble(
                            -m_pressure.val() * m_space * var1(m_space,m_def).tr()* meas(m_ori)
                            //-m_pressure.val() * jac(m_space) * sn(m_def).normalized() * meas(m_ori)
                        );
}

// assembles eq 3.25 from http://resolver.tudelft.nl/uuid:56c0cc91-643d-4817-9702-93fedce5fd78
template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
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
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & /*presFun*/, const gsFunctionSet<T> & /*deformed*/)
{
    // Since pressure works out-of-plane, this function has no effect
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assembleFoundation(const gsFunction<T> & foundFun)
{
    this->_getOptions();
    _assembleFoundation_impl<d,_matrix>(foundFun);
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & foundFun)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_foundation = m_assembler.getCoeff(foundFun, m_ori);
    GISMO_ASSERT(foundFun.targetDim()==3,"Foundation function has dimension "<<foundFun.targetDim()<<", but expected 3");

    m_assembler.assemble(
        m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori)
        );
}

// assembles eq 3.27 from http://resolver.tudelft.nl/uuid:56c0cc91-643d-4817-9702-93fedce5fd78
template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & /* foundFun */)
{
    // No rhs contribution for the linear case
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3), void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & )
{
    // Since foundation works out-of-plane, this function has no effect
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assembleFoundation(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed)
{
    this->_getOptions();
    _assembleFoundation_impl<d,_matrix>(foundFun,deformed);
}

// assembles eq 3.28 from http://resolver.tudelft.nl/uuid:56c0cc91-643d-4817-9702-93fedce5fd78
template <short_t d, class T, bool bending>
template <short_t _d, bool matrix>
typename std::enable_if<(_d==3) && matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & foundFun, const gsFunctionSet<T> & /*deformed*/)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_foundation = m_assembler.getCoeff(foundFun, m_ori);
    GISMO_ASSERT(foundFun.targetDim()==3,"Foundation function has dimension "<<foundFun.targetDim()<<", but expected 3");

    m_assembler.assemble(
            m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori)
        );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    auto m_foundation = m_assembler.getCoeff(foundFun, m_ori);
    GISMO_ASSERT(foundFun.targetDim()==3,"Foundation function has dimension "<<foundFun.targetDim()<<", but expected 3");

    // Assemble vector
    m_assembler.assemble(
                  m_space * m_foundation.asDiag() * (m_def - m_ori) * meas(m_ori) // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
                );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3), void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & , const gsFunctionSet<T> & )
{
    // Since foundation works out-of-plane, this function has no effect
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assembleWeakBCs()
{
    this->_getOptions();
    _assembleWeakBCs_impl<d,_matrix>();
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl()
{
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmD = m_assembler.getCoeff(m_mmD);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);
    auto alpha_r = m_alpha_r_bc * reshape(mmDcart,9,1).max().val() / el.area(m_ori);

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        -(alpha_d * m_space * m_space.tr()) * meas(m_ori)
    );

    // for weak clamped
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Clamped")
        ,
        (
            alpha_r * ( ( var1(m_space,m_ori) * unv(m_ori) ) * ( var1(m_space,m_ori) * unv(m_ori) ).tr() )
        ) * meas(m_ori)
    );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl()
{
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);

    // Weak BCs

    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        -(alpha_d * m_space * g_N         ) * meas(m_ori)
    );

    // for weak clamped
    // do nothing
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl()
{
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        -(alpha_d * m_space * m_space.tr()) * meas(m_ori)
    );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl()
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        (m_alpha_d_bc * (m_space * (m_ori - m_ori) - m_space * (g_N) )) * meas(m_ori)
    );
}

template <short_t d, class T, bool bending>
template <bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assembleWeakBCs(const gsFunctionSet<T> & deformed)
{
    this->_getOptions();
    _assembleWeakBCs_impl<d,_matrix>(deformed);
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmD = m_assembler.getCoeff(m_mmD);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);
    auto alpha_r = m_alpha_r_bc * reshape(mmDcart,9,1).max().val() / el.area(m_ori);


    auto du  = m_def - m_ori;
    auto dnN = ( usn(m_def).tr()*unv(m_ori) - usn(m_ori).tr()*unv(m_ori) ).val();

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        -alpha_d * m_space * m_space.tr() * meas(m_ori)
    );

    // for weak clamped
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Clamped")
        ,
        (
            alpha_r * dnN * ( var2deriv2(m_space,m_space,m_def,unv(m_ori).tr()) )
            +
            alpha_r * ( ( var1(m_space,m_def) * unv(m_ori) ) * ( var1(m_space,m_def) * unv(m_ori) ).tr() )
        ) * meas(m_ori)
    );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmD = m_assembler.getCoeff(m_mmD);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);
    auto alpha_r = m_alpha_r_bc * reshape(mmDcart,9,1).max().val() / el.area(m_ori);

    auto du  = m_def - m_ori;
    auto dnN = ( usn(m_def).tr()*nv(m_ori) - usn(m_ori).tr()*nv(m_ori) ).val();

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        alpha_d * (m_space * du - m_space * (g_N) ) * meas(m_ori)
    );

    // for weak clamped
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Clamped")
        ,
        (
            - alpha_r * dnN * ( var1(m_space,m_def) * usn(m_ori) )
        ) * meas(m_ori)
    );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        -alpha_d * m_space * m_space.tr() * meas(m_ori)
    );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el = m_assembler.getElement();
    auto alpha_d = m_alpha_d_bc * reshape(mmAcart,9,1).max().val() / el.area(m_ori);

    // Weak BCs
    m_assembler.assembleBdr
    (
        m_bcs.get("Weak Dirichlet")
        ,
        alpha_d * (m_space * (m_def - m_ori) - m_space * (g_N) ) * meas(m_ori)
    );
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assembleWeakIfc()
{
    this->_getOptions();
    if (m_weakC0.size()==0 && m_weakC1.size()==0)
        return;
    _assembleWeakIfc_impl<d,_matrix>();
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl()
{
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // auto g_N = m_assembler.getBdrFunction(m_ori);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmD = m_assembler.getCoeff(m_mmD);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);

    element el   = m_assembler.getElement();
    auto h       = (el.area(m_ori.left()) + el.area(m_ori.right())) / 2;
    auto alpha_d = m_alpha_d_ifc * reshape(mmAcart,9,1).max().val() / h;
    auto alpha_r = m_alpha_r_ifc * reshape(mmDcart,9,1).max().val() / h;

    // C^0 coupling
    m_assembler.assembleIfc(m_weakC0,
                     alpha_d * m_space.left() * m_space.left().tr() * meas(m_ori)
                    ,
                    -alpha_d * m_space.right()* m_space.left() .tr() * meas(m_ori)
                    ,
                    -alpha_d * m_space.left() * m_space.right().tr() * meas(m_ori)
                    ,
                     alpha_d * m_space.right()* m_space.right().tr() * meas(m_ori)
                     );

    // C^1 coupling
    // Penalty of out-of-plane coupling
    // dW^pr / du_r --> second line
    m_assembler.assembleIfc(m_weakC1,
                     alpha_r * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ) * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ).tr() * meas(m_ori)    // left left
                    ,
                     alpha_r * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ) * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ).tr() * meas(m_ori)   // left right
                    ,
                     alpha_r * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ) * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ).tr() * meas(m_ori)   // right left
                    ,
                     alpha_r * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ) * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ).tr() * meas(m_ori)  // right right
                    ,
                    // Symmetry
                     alpha_r * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ) * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ).tr() * meas(m_ori)    // right right
                    ,
                     alpha_r * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ) * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ).tr() * meas(m_ori)   // right left
                    ,
                     alpha_r * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ) * ( var1(m_space.right(),m_ori.right()) * usn(m_ori.left()) ).tr() * meas(m_ori)   // left right
                    ,
                     alpha_r * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ) * ( var1(m_space.left(),m_ori.left()) * usn(m_ori.right()) ).tr() * meas(m_ori)  // left left
                     );

    // Penalty of in-plane coupling
    // dW^pr / du_r --> fourth line
    m_assembler.assembleIfc(m_weakC1,
                     alpha_r * ( ovar1(m_space.left() ,m_ori.left() ) * usn(m_ori.right()) ) * ( ovar1(m_space.left() ,m_ori.left() ) * usn(m_ori.right()) ).tr() * meas(m_ori) // left left
                    + // Symmetry
                     alpha_r * (  var1(m_space.left() ,m_ori.left() ) * unv(m_ori.right()) ) * (  var1(m_space.left() ,m_ori.left() ) * unv(m_ori.right()) ).tr() * meas(m_ori) // left left
                    ,
                     alpha_r * ( ovar1(m_space.left() ,m_ori.left() ) * usn(m_ori.right()) ) * (  var1(m_space.right(),m_ori.right()) * unv(m_ori.left() ) ).tr() * meas(m_ori) // left right
                    + // Symmetry
                     alpha_r * (  var1(m_space.left() ,m_ori.left() ) * unv(m_ori.right()) ) * ( ovar1(m_space.right(),m_ori.right()) * usn(m_ori.left() ) ).tr() * meas(m_ori) // left right
                    ,
                     alpha_r * (  var1(m_space.right(),m_ori.right()) * unv(m_ori.left() ) ) * ( ovar1(m_space.left() ,m_ori.left() ) * usn(m_ori.right()) ).tr() * meas(m_ori) // right left
                    + // Symmetry
                     alpha_r * ( ovar1(m_space.right(),m_ori.right()) * usn(m_ori.left() ) ) * (  var1(m_space.left() ,m_ori.left() ) * unv(m_ori.right()) ).tr() * meas(m_ori) // right left
                    ,
                     alpha_r * (  var1(m_space.right(),m_ori.right()) * unv(m_ori.left() ) ) * (  var1(m_space.right(),m_ori.right()) * unv(m_ori.left() ) ).tr() * meas(m_ori) // right right
                    + // Symmetry
                     alpha_r * ( ovar1(m_space.right(),m_ori.right()) * usn(m_ori.left() ) ) * ( ovar1(m_space.right(),m_ori.right()) * usn(m_ori.left() ) ).tr() * meas(m_ori) // right right
                     );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl()
{
/*
    empty
 */
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl()
{
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el   = m_assembler.getElement();
    auto h       = (el.area(m_ori.left()) + el.area(m_ori.right())) / 2;
    auto alpha_d = m_alpha_d_ifc * reshape(mmAcart,9,1).max().val() / h;

    // C^0 coupling
    m_assembler.assembleIfc(m_weakC0,
                     alpha_d * m_space.left() * m_space.left().tr() * meas(m_ori) * meas(m_ori)
                    ,
                    -alpha_d * m_space.right()* m_space.left() .tr() * meas(m_ori) * meas(m_ori)
                    ,
                    -alpha_d * m_space.left() * m_space.right().tr() * meas(m_ori) * meas(m_ori)
                    ,
                     alpha_d * m_space.right()* m_space.right().tr() * meas(m_ori) * meas(m_ori)
                     );

    // C^1 coupling DOES NOT CONTRIBUTE IN 2D PROBLEMS
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl()
{
/*
    empty
 */
}

template<short_t d, class T, bool bending>
template<bool _matrix>
void gsThinShellAssembler<d, T, bending>::_assembleWeakIfc(const gsFunctionSet<T> & deformed)
{
    this->_getOptions();
    _assembleWeakIfc_impl<d,_matrix>(deformed);
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmD = m_assembler.getCoeff(m_mmD);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);

    element el   = m_assembler.getElement();
    auto h       = (el.area(m_ori.left()) + el.area(m_ori.right())) / 2;
    auto alpha_d = m_alpha_d_ifc * reshape(mmAcart,9,1).max().val() / h;
    auto alpha_r = m_alpha_r_ifc * reshape(mmDcart,9,1).max().val() / h;

    auto du = ((m_def.left()-m_ori.left()) - (m_def.right()-m_ori.right()));

    auto dN_lr = (usn(m_def.left()).tr()*usn(m_def.right())
                    - usn(m_ori.left()).tr()*usn(m_ori.right())).val();

    auto dN_rl = (usn(m_def.right()).tr()*usn(m_def.left())
                    - usn(m_ori.right()).tr()*usn(m_ori.left())).val();

    auto dnN_lr= (unv(m_def.left()).tr()*usn(m_def.right())
                    - unv(m_ori.left()).tr()*usn(m_ori.right())).val();

    auto dnN_rl= (unv(m_def.right()).tr()*usn(m_def.left())
                    - unv(m_ori.right()).tr()*usn(m_ori.left())).val();

    // C^0 coupling
    m_assembler.assembleIfc(m_weakC0,
                     alpha_d * m_space.left() * m_space.left().tr() * meas(m_ori)
                    ,
                    -alpha_d * m_space.right()* m_space.left() .tr() * meas(m_ori)
                    ,
                    -alpha_d * m_space.left() * m_space.right().tr() * meas(m_ori)
                    ,
                     alpha_d * m_space.right()* m_space.right().tr() * meas(m_ori)
                     );

    // C^1 coupling
    // Penalty of out-of-plane coupling
    // dW^pr / du_r --> first line
    m_assembler.assembleIfc(m_weakC1,
                     alpha_r * dN_lr * var2(m_space.left() ,m_space.left() ,m_def.left() ,usn(m_def.right()).tr() ) * meas(m_ori)      // left left
                     +//Symmetry
                     alpha_r * dN_rl * var2( m_space.left(),m_space.left(),m_def.left(),usn(m_def.right() ).tr() ) * meas(m_ori)     // left left
                    ,
                     alpha_r * dN_lr * ( var1(m_space.left() ,m_def.left() ) * var1(m_space.right(),m_def.right()).tr() ) * meas(m_ori)// left right
                     +//Symmetry
                     alpha_r * dN_rl * ( var1(m_space.left(),m_def.left()) * var1(m_space.right() ,m_def.right() ).tr() ) * meas(m_ori)// left right
                    ,
                     alpha_r * dN_lr * ( var1(m_space.right(),m_def.right()) * var1(m_space.left() ,m_def.left() ).tr() ) * meas(m_ori)// right left
                     +//Symmetry
                     alpha_r * dN_rl * ( var1(m_space.right() ,m_def.right() ) * var1(m_space.left(),m_def.left()).tr() ) * meas(m_ori)// right left
                    ,
                     alpha_r * dN_lr * var2( m_space.right(),m_space.right(),m_def.right(),usn(m_def.left() ).tr() ) * meas(m_ori)     // right right
                     +//Symmetry
                     alpha_r * dN_rl * var2(m_space.right() ,m_space.right() ,m_def.right() ,usn(m_def.left()).tr() ) * meas(m_ori)      // right right
                     );

    // Penalty of out-of-plane coupling
    // dW^pr / du_r --> second line
    m_assembler.assembleIfc(m_weakC1,
                     alpha_r * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ) * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ).tr() * meas(m_ori)    // left left
                    ,
                     alpha_r * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ) * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ).tr() * meas(m_ori)   // left right
                    ,
                     alpha_r * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ) * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ).tr() * meas(m_ori)   // right left
                    ,
                     alpha_r * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ) * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ).tr() * meas(m_ori)  // right right
                    ,
                    // Symmetry
                     alpha_r * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ) * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ).tr() * meas(m_ori)    // right right
                    ,
                     alpha_r * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ) * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ).tr() * meas(m_ori)   // right left
                    ,
                     alpha_r * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ) * ( var1(m_space.right(),m_def.right()) * usn(m_def.left()) ).tr() * meas(m_ori)   // left right
                    ,
                     alpha_r * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ) * ( var1(m_space.left(),m_def.left()) * usn(m_def.right()) ).tr() * meas(m_ori)  // left left
                     );

    // Penalty of in-plane coupling
    // dW^pr / du_r --> third line
    m_assembler.assembleIfc(m_weakC1,
                     alpha_r * dnN_lr * ovar2(m_space.left(),m_space.left(),m_def.left(),usn(m_def.right()).tr()) * meas(m_ori) // left left
                    + // Symmetry
                     alpha_r * dnN_rl * ovar2(m_space.left(),m_space.left(),m_def.left(),usn(m_def.right()).tr()) * meas(m_ori) // left left
                    ,
                     alpha_r * dnN_lr * ( ovar1(m_space.left() ,m_def.left() ) * var1(m_space.right(),m_def.right()).tr() ) * meas(m_ori) // left right
                    + // Symmetry
                     alpha_r * dnN_rl * ( ovar1(m_space.left() ,m_def.left() ) * var1(m_space.right(),m_def.right()).tr() ) * meas(m_ori) // right left
                    ,
                     alpha_r * dnN_lr * ( ovar1(m_space.right(),m_def.right()) * var1(m_space.left() ,m_def.left() ).tr() ) * meas(m_ori) // right left
                    + // Symmetry
                     alpha_r * dnN_rl * ( ovar1(m_space.right(),m_def.right()) * var1(m_space.left() ,m_def.left() ).tr() ) * meas(m_ori) // right left
                    ,
                     alpha_r * dnN_lr * ovar2(m_space.right(),m_space.right(),m_def.right(),usn(m_def.left()).tr()) * meas(m_ori) // right right
                    + // Symmetry
                     alpha_r * dnN_rl * ovar2(m_space.right(),m_space.right(),m_def.right(),usn(m_def.left()).tr()) * meas(m_ori) // right right
                     );

    // Penalty of in-plane coupling
    // dW^pr / du_r --> fourth line
    m_assembler.assembleIfc(m_weakC1,
                     alpha_r * ( ovar1(m_space.left() ,m_def.left() ) * usn(m_def.right()) ) * ( ovar1(m_space.left() ,m_def.left() ) * usn(m_def.right()) ).tr() * meas(m_ori) // left left
                    + // Symmetry
                     alpha_r * (  var1(m_space.left() ,m_def.left() ) * unv(m_def.right()) ) * (  var1(m_space.left() ,m_def.left() ) * unv(m_def.right()) ).tr() * meas(m_ori) // left left
                    ,
                     alpha_r * ( ovar1(m_space.left() ,m_def.left() ) * usn(m_def.right()) ) * (  var1(m_space.right(),m_def.right()) * unv(m_def.left() ) ).tr() * meas(m_ori) // left right
                    + // Symmetry
                     alpha_r * (  var1(m_space.left() ,m_def.left() ) * unv(m_def.right()) ) * ( ovar1(m_space.right(),m_def.right()) * usn(m_def.left() ) ).tr() * meas(m_ori) // left right
                    ,
                     alpha_r * (  var1(m_space.right(),m_def.right()) * unv(m_def.left() ) ) * ( ovar1(m_space.left() ,m_def.left() ) * usn(m_def.right()) ).tr() * meas(m_ori) // right left
                    + // Symmetry
                     alpha_r * ( ovar1(m_space.right(),m_def.right()) * usn(m_def.left() ) ) * (  var1(m_space.left() ,m_def.left() ) * unv(m_def.right()) ).tr() * meas(m_ori) // right left
                    ,
                     alpha_r * (  var1(m_space.right(),m_def.right()) * unv(m_def.left() ) ) * (  var1(m_space.right(),m_def.right()) * unv(m_def.left() ) ).tr() * meas(m_ori) // right right
                    + // Symmetry
                     alpha_r * ( ovar1(m_space.right(),m_def.right()) * usn(m_def.left() ) ) * ( ovar1(m_space.right(),m_def.right()) * usn(m_def.left() ) ).tr() * meas(m_ori) // right right
                     );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmD = m_assembler.getCoeff(m_mmD);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);

    element el   = m_assembler.getElement();
    auto h       = (el.area(m_ori.left()) + el.area(m_ori.right())) / 2;
    auto alpha_d = m_alpha_d_ifc * reshape(mmAcart,9,1).max().val() / h;
    auto alpha_r = m_alpha_r_ifc * reshape(mmDcart,9,1).max().val() / h;

    auto du = ((m_def.left()-m_ori.left()) - (m_def.right()-m_ori.right()));

    auto dN_lr = (usn(m_def.left()).tr()*usn(m_def.right())
                    - usn(m_ori.left()).tr()*usn(m_ori.right())).val();

    auto dN_rl = (usn(m_def.right()).tr()*usn(m_def.left())
                    - usn(m_ori.right()).tr()*usn(m_ori.left())).val();

    auto dnN_lr= (unv(m_def.left()).tr()*usn(m_def.right())
                    - unv(m_ori.left()).tr()*usn(m_ori.right())).val();

    auto dnN_rl= (unv(m_def.right()).tr()*usn(m_def.left())
                    - unv(m_ori.right()).tr()*usn(m_ori.left())).val();

    // C^0 coupling
    m_assembler.assembleIfc(m_weakC0,
                    -alpha_d * m_space.left() * du * meas(m_ori)
                    ,
                     alpha_d * m_space.right()* du * meas(m_ori)
                     );

   // C^1 coupling
    m_assembler.assembleIfc(m_weakC1,
                    -alpha_r * dN_lr * var1(m_space.left(),m_def.left())   * usn(m_def.right()) * meas(m_ori)
                    ,
                    -alpha_r * dN_lr * var1(m_space.right(),m_def.right()) * usn(m_def.left() ) * meas(m_ori)
                    ,
                    // Symmetry
                    -alpha_r * dN_rl * var1(m_space.right(),m_def.right())   * usn(m_def.left()) * meas(m_ori)
                    ,
                    -alpha_r * dN_rl * var1(m_space.left(),m_def.left()) * usn(m_def.right() ) * meas(m_ori)
                     );

    // Penalty of in-plane coupling
    // dW^pr / du_r --> second line
    m_assembler.assembleIfc(m_weakC1,
                    -alpha_r * dnN_lr* ovar1(m_space.left(),m_def.left())  * usn(m_def.right()) * meas(m_ori)
                    ,
                    -alpha_r * dnN_lr* var1(m_space.right(),m_def.right()) * unv(m_def.left() ) * meas(m_ori)
                    ,
                    // Symmetry
                    -alpha_r * dnN_rl* ovar1(m_space.right(),m_def.right())  * usn(m_def.left()) * meas(m_ori)
                    ,
                    -alpha_r * dnN_rl* var1(m_space.left(),m_def.left()) * unv(m_def.right() ) * meas(m_ori)
                     );
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && _matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el   = m_assembler.getElement();
    auto h       = (el.area(m_ori.left()) + el.area(m_ori.right())) / 2;
    auto alpha_d = m_alpha_d_ifc * reshape(mmAcart,9,1).max().val() / h;

    auto du = ((m_def.left()-m_ori.left()) - (m_def.right()-m_ori.right()));

    // C^0 coupling
    m_assembler.assembleIfc(m_weakC0,
                     alpha_d * m_space.left() * m_space.left().tr() * meas(m_ori)
                    ,
                    -alpha_d * m_space.right()* m_space.left() .tr() * meas(m_ori)
                    ,
                    -alpha_d * m_space.left() * m_space.right().tr() * meas(m_ori)
                    ,
                     alpha_d * m_space.right()* m_space.right().tr() * meas(m_ori)
                     );

    // C^1 coupling DOES NOT CONTRIBUTE IN 2D PROBLEMS
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<!(_d==3) && !_matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakIfc_impl(const gsFunctionSet<T> & deformed)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);

    auto cart2cov = cartcov(m_ori);
    auto con2cartI = cartcon(m_ori).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);

    element el   = m_assembler.getElement();
    auto h       = (el.area(m_ori.left()) + el.area(m_ori.right())) / 2;
    auto alpha_d = m_alpha_d_ifc * reshape(mmAcart,9,1).max().val() / h;

    auto du = ((m_def.left()-m_ori.left()) - (m_def.right()-m_ori.right()));

    // C^0 coupling
     m_assembler.assembleIfc(m_weakC0,
                      alpha_d * m_space.left() * du * meas(m_ori)
                     ,
                     -alpha_d * m_space.right()* du * meas(m_ori)
                      );

    // C^1 coupling DOES NOT CONTRIBUTE IN 2D PROBLEMS
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_assembleDirichlet()
{
    this->_getOptions();
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // if statement
    m_space.setup(m_bcs, dirichlet::l2Projection, m_continuity);
    // m_assembler.initSystem();
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::homogenizeDirichlet()
{
    this->_getOptions();
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    m_space.setup(m_bcs, dirichlet::homogeneous, m_continuity);
    // space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart().setZero();
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_applyLoads()
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

/*
    Adds the point masses of m_pMass to the ALREADY ASSEMBLED mass matrix m_mass.

    Each point mass deposits the CONSISTENT (rank-one) block
    value * N_k(p) * N_l(p) into every one of the d component blocks, restricted to
    the FREE dofs -- hence the k x l double loop below, and hence the fact that this
    routine must run BEFORE any row-sum lumping (see assembleMass).

    CONTRACT -- five measured behaviours; the first two are handled here, the other
    three are LOUD already and are documented rather than converted:

    (a) The point must lie in the parameter domain of the patch's basis. An
        out-of-domain point used to NaN the ENTIRE mass matrix while assembleMass()
        still reported Success: silent corruption under a success status, which is
        what every dynamic-analysis driver branches on. It is now rejected by the
        GISMO_ENSURE below. ENSURE, not ASSERT: the point is a CALLER-SUPPLIED value
        of the public setPointMass() API, and this routine runs once per point mass
        (not per quadrature point), so the guard is live under -DNDEBUG as well. The
        std::runtime_error it throws is caught by assembleMass()'s own catch(...)
        and surfaces to the caller as AssemblyError -- with the full diagnostic
        (the point and the domain) already on std::cerr, because GISMO_ENSURE prints
        condition, message, file and line BEFORE it throws
        (src/gsCore/gsDebug.h:120-124); only what() is the bare "GISMO_ENSURE".
        CONTRACT WIDENING this implies: a point mass whose VALUE is zero at an
        out-of-domain point used to be a silent no-op -- the "value != 0" test below
        gated the write, so the NaN basis values were never consumed -- and now
        raises AssemblyError as well. Deliberate: the input is wrong either way, and
        distinguishing the two would make the guard value-dependent.

        The bound is read from gsBasis::support(). TWO residues follow from that,
        and they are NOT the same thing:
          - the TOLERANCE band (see the guard below): points just outside the domain
            are accepted and then CLAMPED onto it, so their mass is deposited at the
            boundary instead of being lost. Closed.
          - the BOUNDING BOX: support() is documented as *a bounding box for* the
            domain, so a point in box \ domain is admitted. The clamp does NOT
            relocate such a point (it is already inside the box), and it would
            evaluate to sum(N) = 0, i.e. a dropped mass under Success (a NaN one on
            a rational basis, see the clamp comment below). This residue
            is EMPTY for every basis this assembler can reach: for a tensor /
            B-spline / NURBS patch the parameter domain IS its box, and
            gsHTensorBasis::support() (gsHTensorBasis.hpp:122-126) returns the
            LEVEL-0 box, which hierarchical refinement does not shrink. It is
            documented rather than guarded because closing it in general needs a
            per-basis domain-membership query that gsBasis does not offer.
        A basis type that does not implement support() raises
        GISMO_NO_IMPLEMENTATION, again an AssemblyError rather than a silent NaN.

    (b) The actives must be numbered in the SPACE basis, not in the integration
        basis. m_mapper is read from m_space (:1485), and _initialize() builds
        m_space from *m_spaceBasis (:259) -- m_basis is only the integration basis
        (:248). The two diverge as soon as the public setSpaceBasis() is used, and
        this routine used to resolve the basis as m_basis.front().basis(patch): the
        actives were numbered in one basis and resolved in the other. MEASURED on the
        Scordelis-Lo roof with m_basis = coarse and setSpaceBasis(refined), a 7.5
        point mass at the interior point (0.5,0.5): the mass landed on an entirely
        different dof set (rows 44 45 54 55 363 364 372 651 652 660 instead of
        152 153 170 171 459 460 475 476 747 748 763 764) and the total deposited was
        15.9375 instead of d*value = 22.5 -- silently, under Success, because the
        grand sum is blind to WHERE the mass lands and some of the mis-resolved
        indices hit ELIMINATED dofs. The reverse split (m_basis finer than the space
        basis) reads m_dofs[c] past patchSize(0,c) and returns an unrelated index
        that is_free_index() accepts. Repaired by dispatching over *m_spaceBasis
        with the SAME type ladder _applyLoads carries (:1337-1369).
        SIDE EFFECT, measured: a point mass on patch > 0 used to THROW through
        gsBasis::piece()'s GISMO_ENSURE(0==k) (src/gsCore/gsBasis.h:105-109), because
        m_basis.front() is patch 0's basis and .basis(patch) then asks it for a piece
        it does not have. gsMultiBasis::basis(patch) indexes patches properly, so that
        throw is gone and a multipatch point mass is now applied to its own patch --
        verified on a two-patch DISCONNECTED fixture (see the report of task 63; the
        deposited block is identical to the patch-0 one and sits on patch 1's dofs).
        A glued multipatch was not probed.
        That same throw was, however, ALSO the only thing stopping an OUT-OF-RANGE
        patch index in a Release build, so the patch test below is a GISMO_ENSURE
        and not the GISMO_ASSERT _applyLoads carries -- the one deliberate deviation
        from that routine's literal form, and the reason is written at the line.

    (c) On a problem with ZERO free dofs the GISMO_ASSERT(m_mass.rows()!=0) below
        fires even when there is no point mass at all, because assembleMass() calls
        this routine unconditionally; assembleMass() then returns AssemblyError.
        Under -DNDEBUG that assert is absent, the empty loop is skipped and the same
        call returns Success -- Debug and Release diverge on that (pathological)
        fixture. Documented, not guarded: both branches are visible, neither is
        silently wrong.

    (d) A NEGATIVE mass is PERMITTED and is not validated. It produces a negative
        contribution and hence a NEGATIVE lumped diagonal entry (measured: a
        -1e6 point mass gives a minimum diagonal of -2.5e5), which makes the lumped
        operator indefinite and unusable in a generalized eigenproblem. This is
        deliberate: negative point masses are legitimate input in model updating.
*/
template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_applyMass()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts,globalActs;

    space       m_space = m_assembler.trialSpace(0);
    m_mapper = m_space.mapper();

    GISMO_ASSERT(m_mass.rows()!=0,"Mass matrix must be assembled first");

    for (size_t i = 0; i< m_pMass.numLoads(); ++i )
    {
        GISMO_ASSERT(m_pMass[i].value.size()==1,"Mass should be one-dimensional");
        // ENSURE, not the ASSERT _applyLoads uses -- see contract (b): until this
        // routine dispatched over m_spaceBasis, EVERY patch != 0 was stopped by
        // gsBasis::piece()'s own GISMO_ENSURE, which is live under -DNDEBUG. Making
        // the valid patches work removes that net, and an ASSERT here would leave
        // patch >= nPatches() to fall through into gsMultiBasis::basis(patch) ->
        // *m_bases[patch] in a Release build, i.e. an out-of-range deref instead of
        // a throw. Same argument as (a): a caller-supplied value of a public setter,
        // tested once per point mass.
        GISMO_ENSURE((size_t)m_pMass[i].patch<m_patches.nPatches(),"Point mass is defined on a patch with index "<<m_pMass[i].patch<<" while the geometry has "<<m_patches.nPatches()<<" patches\n");

        // Resolve the PARAMETRIC point at which this mass acts. Both branches feed
        // the same active_into/eval_into below, so the domain guard covers the
        // physical-space input too -- an inversion that does not converge lands
        // outside the domain or returns NaN, and NaN fails both comparisons.
        gsMatrix<T> parPoint;
        if ( m_pMass[i].parametric )   // in parametric space
            parPoint = m_pMass[i].point;
        else                            // in physical space
            m_patches.patch(m_pMass[i].patch).invertPoints(m_pMass[i].point,parPoint);

        // Contract (b): the actives are resolved through m_mapper, which comes from
        // the space, which is built from *m_spaceBasis -- so they must be NUMBERED
        // in *m_spaceBasis, not in the integration basis m_basis. Same type ladder,
        // same order and same fallback as _applyLoads (:1337-1369); the two casts
        // are hoisted out of it only because the domain guard has to sit BETWEEN
        // support() and active_into(), and duplicating the guard in each branch
        // would be worse than hoisting the dispatch.
        const gsMappedBasis<2,T> * mappedBasis = dynamic_cast<const gsMappedBasis<2,T> * >(m_spaceBasis);
        const gsMultiBasis<T>    * multiBasis  = dynamic_cast<const gsMultiBasis<T>    * >(m_spaceBasis);
        GISMO_ENSURE(mappedBasis!=nullptr || multiBasis!=nullptr,"Basis type not understood");

        // Contract (a): out-of-domain points used to NaN the WHOLE mass matrix
        // while assembleMass() still reported Success. The bounds are READ FROM THE
        // SAME BASIS the actives come from -- the parameter domain is not
        // guaranteed to be [0,1]^d, and reading it from a different basis would
        // re-create contract (b) in a new place.
        const gsMatrix<T> supp = (mappedBasis!=nullptr)
                               ? mappedBasis->getMappedSingleBasis(m_pMass[i].patch).support()
                               : multiBasis->basis(m_pMass[i].patch).support();
        // This dimension test MUST precede the first .col(0) access below.
        GISMO_ENSURE(parPoint.rows()==supp.rows(),
                     "Point mass "<<i<<" has "<<parPoint.rows()<<" coordinate(s), but the basis of patch "
                     <<m_pMass[i].patch<<" has parameter dimension "<<supp.rows());
        // Tolerance in units of eps, on the extent of the domain, so that a float
        // or multiprecision T scales with the arithmetic instead of degenerating
        // into an exact test (for an exact-arithmetic T with epsilon() == 0 it IS
        // an exact test, which is correct). The 1e3 is a round allowance for the
        // round-off accumulated by whatever produced the point -- a coordinate
        // computed from a knot vector, or a physical point run through
        // invertPoints; it is NOT there to cover a Newton iterate landing outside
        // the domain, which cannot happen: gsGeometry::invertPoints ->
        // gsFunction::newtonRaphson(...,withSupport=true) clamps EVERY iterate with
        // arg.cwiseMax(supp.col(0)).cwiseMin(supp.col(1)) against this same
        // support(), and writes +inf on failure.
        // Note the real_t scaling of the BAND, which is why the clamp below is not
        // optional: in double it is ~2.2e-13*extent, but for real_t = float it is
        // ~1.2e-4*extent, i.e. wider than one element on a mesh finer than that.
        // A point admitted by tol is therefore RELOCATED onto the domain, never
        // silently dropped, whatever the width of the band.
        const T tol = 1e3 * std::numeric_limits<T>::epsilon()
                          * (supp.col(1)-supp.col(0)).cwiseAbs().maxCoeff();
        GISMO_ENSURE( ((parPoint.col(0).array()-supp.col(0).array()) >= -tol).all() &&
                      ((parPoint.col(0).array()-supp.col(1).array()) <=  tol).all(),
                      "Point mass "<<i<<" lies outside the parameter domain of patch "
                      <<m_pMass[i].patch<<": the point is ("<<parPoint.col(0).transpose()
                      <<") while the domain is ("<<supp.col(0).transpose()<<") x ("
                      <<supp.col(1).transpose()<<")");

        // ... and CLAMP what the tolerance admitted. Without this, a point ACCEPTED
        // inside the band was NOT harmless: at u = 1+1e-14 on a [0,1] knot vector
        // active_into returns a SHIFTED active set (measured on a B-spline basis:
        // 12..26 instead of 15..29) on which every N vanishes. On a POLYNOMIAL basis
        // the k x l loop then deposits exactly nothing -- the mass is silently lost
        // under Success. On a RATIONAL one (this class' own Scordelis-Lo fixture is
        // a NURBS) it is worse: gsRationalBasis divides by sum_i w_i N_i, which is
        // ALSO zero, so eval_into returns NaN and the WHOLE mass matrix goes NaN --
        // measured, still under Success. So the tolerance band re-opened, 2.2e-13
        // wide, the very failure shape (a) exists to remove; the clamp closes it.
        // Clamping is exact and a no-op for any point genuinely inside.
        parPoint.col(0) = parPoint.col(0).cwiseMax(supp.col(0)).cwiseMin(supp.col(1));

        // Compute actives and values of basis functions on point load location.
        if (mappedBasis!=nullptr)
        {
            mappedBasis->active_into(m_pMass[i].patch,parPoint, acts );
            mappedBasis->eval_into  (m_pMass[i].patch,parPoint, bVals);
        }
        else
        {
            multiBasis->basis(m_pMass[i].patch).active_into( parPoint, acts );
            multiBasis->basis(m_pMass[i].patch).eval_into  ( parPoint, bVals);
        }

        // Add the point load values in the right entries in the global RHS
        for (size_t j = 0; j< d; ++j)
        {
            if (m_pMass[i].value[0] != 0.0)
            {
                m_mapper.localToGlobal(acts, m_pMass[i].patch, globalActs,j);
                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    for (index_t l=0; l < globalActs.rows(); ++l)
                    {
                        if (m_mapper.is_free_index(globalActs(k,0)) && m_mapper.is_free_index(globalActs(l,0)))
                            m_mass(globalActs(k,0), globalActs(l,0)) += bVals(k,0) * bVals(l,0) * m_pMass[i].value[0];
                    }
                }
            }
        }
    }
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleMass(const bool lumped)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);

    // Initialize stystem
    m_assembler.initSystem();

    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMatrices,&m_patches);
    auto mm0 = m_assembler.getCoeff(m_mm);

    space       m_space = m_assembler.trialSpace(0);
    // The mass operator carries no Dirichlet lifting, so the space is set up
    // homogeneously here. This OVERWRITES the space's fixed part, which is
    // shared state: it is restored after the try/catch below.
    m_space.setup(m_bcs, dirichlet::homogeneous, m_continuity);

    try
    {
        // The CONSISTENT mass matrix is assembled UNCONDITIONALLY, also when a
        // lumped matrix is asked for. Assembling mm0.val()*m_space.rowSum()*meas()
        // instead does NOT lump: rowSum() is transparent to the Space trait
        // (gsExpressions/rowsum_expr.h:32 declares Space = E::Space), so that
        // expression is VECTOR-valued and gsExprAssembler dispatches it -- at
        // compile time, through push<E::isMatrix()> -- into the rhs. The system
        // matrix is then never written and m_assembler.matrix() returns an
        // unmanaged cache (empty, zero, or the STALE matrix of an earlier call).
        m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
        m_mass = m_assembler.matrix();

        // Point masses are CONSISTENT (rank-one) contributions: _applyMass writes
        // value*N_k*N_l over the FULL k x l double loop (:1591-1598 above), i.e.
        // genuine OFF-DIAGONAL entries. It must therefore run BEFORE the row-sum
        // lumping below -- applying it afterwards would re-introduce off-diagonals,
        // so the result would not be diagonal at all. That ordering is also what
        // makes the point masses reach the lumped diagonal (the dead block's own
        // "To do: add point masses in lumped case").
        this->_applyMass();

        if (lumped)
        {
            // Row-sum lumping of the ASSEMBLED matrix: the gsProjection pattern,
            // src/gsUtils/gsProjection.hpp:59-69. Since sum_i sum_j M_ij IS the
            // grand sum of M, the total mass is conserved to machine precision on
            // ANY fixture, with or without eliminated Dirichlet dofs.
            //
            // COST, DELIBERATE -- DO NOT "OPTIMISE" THIS AWAY: this pays a full
            // consistent assembly plus one O(nnz) sparse mat-vec, where assembling
            // the vector-valued rowSum() expression and reading m_assembler.rhs()
            // would be far cheaper. The cheap route is NOT equivalent: its
            // push<false> path has no column loop, so it also deposits the
            // contributions of the ELIMINATED columns, and the resulting diagonal
            // does not sum to the consistent grand sum on a constrained problem.
            // The two coincide only on a fixture without Dirichlet elimination.
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
    // updateBCs(m_bcs) does, and it runs on the error path too so that a failed
    // mass assembly does not leave the instance in a homogenized state. Sitting
    // after the try/catch is the ONLY load-bearing part of its position:
    // its order relative to _applyMass() is NOT load-bearing. (An earlier comment
    // here claimed that _applyMass "re-reads m_mapper from the space", which is
    // true -- :1485 -- but is not a constraint: setup(bc,homogeneous,cont) and
    // setup(bc,l2Projection,cont) build an IDENTICAL mapper and differ only in the
    // gsDirichletValues call that writes fixedPart, so _applyMass sees the same
    // free/eliminated classification either way.)
    this->_assembleDirichlet();
    m_ddofs  = m_space.fixedPart();
    m_mapper = m_space.mapper();

    return m_status;
}

// legacy
template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleFoundation()
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);

    // Initialize stystem
    m_assembler.initSystem();
    auto    m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
    GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

    space       m_space = m_assembler.trialSpace(0);

    try
    {
        m_assembler.assemble(m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori));
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemble()
{
    return assemble_impl<d, bending>();
}

/**
    @brief Assembles the Kirchhoff-Love shell equations including the bending terms.
    Optionally, pressure is included via \a p * n * u
    Optionally, foundation stiffness is included via k_x v_x v_x + k_y v_y v_y + k_z v_z v_z
    Since the variational energy of the foundation force k_i u_i is equal to k_i u_i v_i where i denotes any direction, u_i are displacemets and v_i are spaces.

*/
template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assemble_impl()
{
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // Linear assembly: deformed and undeformed geometries are the same
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(defpatches);

    // Initialize stystem
    m_assembler.initSystem();
    m_assembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> m_mmB(m_materialMatrices,&m_patches,&defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> m_mmC(m_materialMatrices,&m_patches,&defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmB = m_assembler.getCoeff(m_mmB);
    auto mmC = m_assembler.getCoeff(m_mmC);
    auto mmD = m_assembler.getCoeff(m_mmD);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);

    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) );
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3);

    auto m_N_der    = m_Em_der * reshape(mmA,3,3) + m_Ef_der * reshape(mmB,3,3);
    auto m_M_der    = m_Em_der * reshape(mmC,3,3) + m_Ef_der * reshape(mmD,3,3);

    try
    {
        if (m_foundInd)
        {
            this->_assembleFoundation<true>(*m_foundFun);
            this->_assembleFoundation<false>(*m_foundFun);
        }
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

        this->_assembleWeakBCs<true>();
        this->_assembleWeakBCs<false>();
        this->_assembleWeakIfc<true>();
        this->_assembleWeakIfc<false>();
        this->_assembleNeumann();

        // Assemble the loads
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assemble_impl()
{
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // Linear assembly: deformed and undeformed geometries are the same
    gsMultiPatch<T> & defpatches = m_patches;
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(defpatches);

    // Initialize stystem
    m_assembler.initSystem();

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&defpatches);
    auto mmA = m_assembler.getCoeff(m_mmA);

    space       m_space = m_assembler.trialSpace(0);
    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    auto jacG       = jac(m_def);
    auto m_Em_der   = flat( jacG.tr() * jac(m_space) ) ; //[checked]
    auto m_N_der    = m_Em_der * reshape(mmA,3,3);

    try
    {
        if (m_foundInd)
        {
            this->_assembleFoundation<true>(*m_foundFun);
            this->_assembleFoundation<false>(*m_foundFun);
        }
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

        this->_assembleWeakBCs<true>();
        this->_assembleWeakBCs<false>();
        this->_assembleWeakIfc<true>();
        this->_assembleWeakIfc<false>();
        this->_assembleNeumann();

        // Assemble the loads
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleMatrix(const gsFunctionSet<T> & deformed)
{
    return assembleMatrix_impl<d, bending>(deformed);
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assembleMatrix_impl(const gsFunctionSet<T> & deformed)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize matrix
    m_assembler.initSystem();
    m_assembler.initMatrix();

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> m_mmB(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> m_mmC(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmB = m_assembler.getCoeff(m_mmB);
    auto mmC = m_assembler.getCoeff(m_mmC);
    auto mmD = m_assembler.getCoeff(m_mmD);
    auto S0  = m_assembler.getCoeff(m_S0);
    auto S1  = m_assembler.getCoeff(m_S1);

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

    auto m_N_der    = m_Em_der * reshape(mmA,3,3) + m_Ef_der * reshape(mmB,3,3);
    auto m_M_der    = m_Em_der * reshape(mmC,3,3) + m_Ef_der * reshape(mmD,3,3);

    try
    {
        if (m_foundInd) this->_assembleFoundation<true>(*m_foundFun,deformed);
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

        this->_assembleWeakBCs<true>(deformed);
        this->_assembleWeakIfc<true>(deformed);

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assembleMatrix_impl(const gsFunctionSet<T> & deformed)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize matrix
    m_assembler.initSystem();
    m_assembler.initMatrix();

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto S0  = m_assembler.getCoeff(m_S0);

    space       m_space = m_assembler.trialSpace(0);


    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_N_der    = m_Em_der * reshape(mmA,3,3);

    // Assemble matrix
    try
    {
        if (m_foundInd) this->_assembleFoundation<true>(*m_foundFun,deformed);
        if (m_pressInd) this->_assemblePressure<true>(*m_pressFun,deformed);

        m_assembler.assemble(
                (
                    m_N_der * m_Em_der.tr()
                    +
                    m_Em_der2
                ) * meas(m_ori)
            );
        this->_assembleWeakBCs<true>(deformed);
        this->_assembleWeakIfc<true>(deformed);

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...) // add specific cases?
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleMatrix(const gsMatrix<T> & solVector)
{
    gsMultiPatch<T> def;
    constructSolution(solVector, def);
    return assembleMatrix(def);
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleMatrix(const gsFunctionSet<T> & deformed, const gsFunctionSet<T> & previous, gsMatrix<T> & update)
{
    return assembleMatrix_impl<d, bending>(deformed, previous, update);
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assembleMatrix_impl(const gsFunctionSet<T> & deformed, const gsFunctionSet<T> & previous, gsMatrix<T> & update)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);
    geometryMap m_prev  = m_assembler.getMap(previous);
    // Initialize matrix
    m_assembler.initSystem();
    m_assembler.initMatrix();

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> m_mmB(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> m_mmC(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,&m_patches,&deformed);
    auto mmA = m_assembler.getCoeff(m_mmA);
    auto mmB = m_assembler.getCoeff(m_mmB);
    auto mmC = m_assembler.getCoeff(m_mmC);
    auto mmD = m_assembler.getCoeff(m_mmD);
    // auto S0  = m_assembler.getCoeff(m_S0);
    // auto S1  = m_assembler.getCoeff(m_S1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmAd(m_materialMatrices,&m_patches,&previous);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> m_mmBd(m_materialMatrices,&m_patches,&previous);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> m_mmCd(m_materialMatrices,&m_patches,&previous);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmDd(m_materialMatrices,&m_patches,&previous);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0d(m_materialMatrices,&m_patches,&previous);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1d(m_materialMatrices,&m_patches,&previous);
    auto mmAp = m_assembler.getCoeff(m_mmAd);
    auto mmBp = m_assembler.getCoeff(m_mmBd);
    auto mmCp = m_assembler.getCoeff(m_mmCd);
    auto mmDp = m_assembler.getCoeff(m_mmDd);
    auto S0  = m_assembler.getCoeff(m_S0d);
    auto S1  = m_assembler.getCoeff(m_S1d);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    solution    m_du = m_assembler.getSolution(m_space,update);

    this->homogenizeDirichlet();

    auto m_E_mc = flat( jac(m_prev).tr() * grad(m_du) ) ; //[checked]
    auto m_E_fc = -( deriv2(m_du,sn(m_prev).normalized().tr() ) + deriv2(m_prev,var1(m_du,m_prev) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_N_c  = m_E_mc * reshape(mmAp,3,3) + m_E_fc * reshape(mmBp,3,3);
    auto m_M_c  = m_E_mc * reshape(mmCp,3,3) + m_E_fc * reshape(mmDp,3,3);

    auto m_N        = S0.tr() + m_N_c;
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_M        = S1.tr() + m_M_c; // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Ef_der2  = -(flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M ).symmetrize()
                            + var2deriv2(m_space,m_space,m_def, m_M ));

    auto m_N_der    = m_Em_der * reshape(mmA,3,3) + m_Ef_der * reshape(mmB,3,3);
    auto m_M_der    = m_Em_der * reshape(mmC,3,3) + m_Ef_der * reshape(mmD,3,3);

    try
    {
        if (m_foundInd) this->_assembleFoundation<true>(*m_foundFun,deformed);
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
        this->_assembleWeakBCs<true>(deformed);
        this->_assembleWeakIfc<true>(deformed);

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

// template <short_t d, typename T, bool bending>
// template <short_t _d, bool _bending>
// typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
// gsThinShellAssembler<d, T, bending>::assembleMatrix_impl(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> & previous, gsMatrix<T> & update)
// {
//     GISMO_NO_IMPLEMENTATION;
// }

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleMatrix(const gsMatrix<T> & solVector, const gsMatrix<T> & prevVector)
{
    // gsMultiPatch<T> deformed;
    // constructSolution(solVector, deformed);
    // assembleMatrix(deformed);

    gsMultiPatch<T> def, it;
    constructSolution(solVector, def);
    constructSolution(prevVector, it);
    gsMatrix<T> update = solVector - prevVector;
    return assembleMatrix(def,it,update);
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleVector(const gsFunctionSet<T> & deformed, const bool homogenize)
{
  return assembleVector_impl<d, bending>(deformed,homogenize);
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assembleVector_impl(const gsFunctionSet<T> & deformed, const bool homogenize)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize vector
    m_assembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,&m_patches,&deformed);
    auto S0  = m_assembler.getCoeff(m_S0);
    auto S1  = m_assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
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
        if (m_foundInd) this->_assembleFoundation<false>(*m_foundFun,deformed);
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

        this->_assembleWeakBCs<false>(deformed);
        this->_assembleWeakIfc<false>(deformed);
        this->_assembleNeumann();

        // Assemble the loads
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
gsThinShellAssembler<d, T, bending>::assembleVector_impl(const gsFunctionSet<T> & deformed, const bool homogenize)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    // Initialize vector
    m_assembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
    auto S0  = m_assembler.getCoeff(m_S0);

    space m_space       = m_assembler.trialSpace(0);
    auto m_physforce = m_assembler.getCoeff(*m_forceFun,m_ori); // force defined in physical domain
    auto m_parforce  = m_assembler.getCoeff(*m_forceFun); // force defined in parametric domain

    if (homogenize) this->homogenizeDirichlet();
    else            this->_assembleDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    try
    {
        if (m_foundInd) this->_assembleFoundation<false>(*m_foundFun,deformed);
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

        this->_assembleWeakBCs<false>(deformed);
        this->_assembleWeakIfc<false>(deformed);
        this->_assembleNeumann();

        // Assemble the loads
        if ( m_pLoads.numLoads() != 0 )
        {
            m_rhs = m_assembler.rhs();
            _applyLoads();
        }

        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureMatrix(const gsFunction<T> & pressFun)
{
    try
    {
        this->_assemblePressure<true>(pressFun);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureMatrix(const T pressure)
{
    gsConstantFunction<T> pressFun(pressure,d);
    try
    {
        this->_assemblePressure<true>(pressFun);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureMatrix(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    try
    {
        this->_assemblePressure<true>(pressFun,deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureMatrix(const T pressure, const gsFunctionSet<T> & deformed)
{
    gsConstantFunction<T> pressFun(pressure,d);
    try
    {
        this->_assemblePressure<true>(pressFun, deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureVector(const gsFunction<T>   & pressFun )
{
    try
    {
        this->_assemblePressure<false>(pressFun);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureVector(const T pressure)
{
    gsConstantFunction<T> pressFun(pressure,d);
    try
    {
        this->_assemblePressure<false>(pressFun);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureVector(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed)
{
    try
    {
        this->_assemblePressure<false>(pressFun,deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureVector(const T pressure, const gsFunctionSet<T> & deformed)
{
    gsConstantFunction<T> pressFun(pressure,d);
    try
    {
        this->_assemblePressure<false>(pressFun, deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemblePressureVector(const gsFunctionSet<T> & deformed)
{
    try
    {
        this->assemblePressureVector(*m_pressFun,deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleFoundationVector(const gsFunctionSet<T> & deformed)
{
    try
    {
        this->assembleFoundationVector(*m_foundFun,deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleFoundationVector(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed)
{
    try
    {
        this->_assembleFoundation<false>(foundFun,deformed);
        m_status = ThinShellAssemblerStatus::Success;
    }
    catch (...)
    {
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
    return m_status;
}

template<short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::boundaryForce(const gsFunctionSet<T> & deformed,  const std::vector<patchSide> & patchSides) const
{
    return boundaryForce_impl<d, bending>(deformed,patchSides);
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<(_d==3) && _bending, gsMatrix<T> >::type
gsThinShellAssembler<d, T, bending>::boundaryForce_impl(const gsFunctionSet<T> & deformed, const std::vector<patchSide> & patchSides) const
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(*m_spaceBasis, d, 0); // last argument is the space ID

    gsBoundaryConditions<T> bc;
    u.setup(bc, dirichlet::l2Projection, m_continuity);

    gsVector<T> F(d);
    F.setZero();
    if (const gsMultiBasis<T> * mbasis = dynamic_cast<const gsMultiBasis<T>*>(&u.source()))
    {
        // Collect indices of the functions on the selected boundaries
        std::vector<std::unordered_set<index_t>> indices(d);
        gsMatrix<index_t> boundary;
        for (std::vector<patchSide>::const_iterator bdr = patchSides.begin(); bdr != patchSides.end(); bdr++)
        {
            boundary = mbasis->basis(bdr->patch).boundary(bdr->side());
            for (index_t k=0; k!=boundary.rows(); k++)
                for (index_t c=0; c!=d; c++)
                    indices[c].insert(u.mapper().index(boundary.at(k),bdr->patch,c));
        }

        assembler.initSystem();

        geometryMap m_ori   = assembler.getMap(m_patches);
        geometryMap m_def   = assembler.getMap(deformed);

        // Initialize vector
        // m_assembler.initVector(1);

        gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
        gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,&m_patches,&deformed);
        auto S0  = assembler.getCoeff(m_S0);
        auto S1  = assembler.getCoeff(m_S1);

        gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
        auto m_m2 = assembler.getCoeff(mult2t);

        // this->homogenizeDirichlet();

        auto m_N        = S0.tr();
        auto m_Em_der   = flat( jac(m_def).tr() * jac(u) ) ;

        auto m_M        = S1.tr(); // output is a column
        auto m_Ef_der   = -( deriv2(u,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(u,m_def) ) ) * reshape(m_m2,3,3); //[checked]

        // Assemble vector (slow?)
        try
        {
            assembler.assemble(
                          - ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                        );
        }
        catch (...)
        {
            GISMO_ERROR("Assembly of the force vector failed.");
        }

        // Grab and sum control point forces on boundary indices
        for (index_t c = 0; c != d; c++)
            for (std::unordered_set<index_t>::const_iterator it = indices[c].begin(); it!=indices[c].end(); it++)
                F[c] += assembler.rhs().at(*it);
    }
    else
        GISMO_ERROR("The basis is not a gsMultiBasis!");

    return F;
}

template <short_t d, typename T, bool bending>
template <short_t _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), gsMatrix<T> >::type
gsThinShellAssembler<d, T, bending>::boundaryForce_impl(const gsFunctionSet<T> & deformed, const std::vector<patchSide> & patchSides) const
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(*m_spaceBasis, d, 0); // last argument is the space ID

    gsBoundaryConditions<T> bc;
    u.setup(bc, dirichlet::l2Projection, m_continuity);

    gsVector<T> F(d);
    F.setZero();
    if (const gsMultiBasis<T> * mbasis = dynamic_cast<const gsMultiBasis<T>*>(&u.source()))
    {
        // Collect indices of the functions on the selected boundaries
        std::vector<std::unordered_set<index_t>> indices(d);
        gsMatrix<index_t> boundary;
        for (std::vector<patchSide>::const_iterator bdr = patchSides.begin(); bdr != patchSides.end(); bdr++)
        {
            boundary = mbasis->basis(bdr->patch).boundary(bdr->side());
            for (index_t k=0; k!=boundary.rows(); k++)
                for (index_t c=0; c!=d; c++)
                    indices[c].insert(u.mapper().index(boundary.at(k),bdr->patch,c));
        }

        assembler.initSystem();

        geometryMap m_ori   = assembler.getMap(m_patches);
        geometryMap m_def   = assembler.getMap(deformed);

        // Initialize vector
        // m_assembler.initVector(1);

        gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&m_patches,&deformed);
        auto S0  = assembler.getCoeff(m_S0);

        // this->homogenizeDirichlet();

        auto m_N        = S0.tr();
        auto m_Em_der   = flat( jac(m_def).tr() * jac(u) ) ;

        // Assemble vector (slow?)
        try
        {
            assembler.assemble(
                          - ( ( m_N * m_Em_der.tr() ) * meas(m_ori) ).tr()
                        );
        }
        catch (...)
        {
            GISMO_ERROR("Assembly of the force vector failed.");
        }

        // Grab and sum control point forces on boundary indices
        for (index_t c = 0; c != d; c++)
            for (std::unordered_set<index_t>::const_iterator it = indices[c].begin(); it!=indices[c].end(); it++)
                F[c] += assembler.rhs().at(*it);
    }
    else
        GISMO_ERROR("The basis is not a gsMultiBasis!");

    return F;
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assembleVector(const gsMatrix<T> & solVector, const bool homogenize)
{
    gsMultiPatch<T> def;
    constructSolution(solVector, def);
    return assembleVector(def,homogenize);
}

template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemble(const gsFunctionSet<T> & deformed,
                                                   const bool Matrix, const bool homogenize)
{
    ThinShellAssemblerStatus status;
    if (Matrix)
    {
        status = assembleMatrix(deformed);
        if (status!=ThinShellAssemblerStatus::Success)
            return status;
    }

    return assembleVector(deformed,homogenize);
}
template<short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssembler<d, T, bending>::assemble(const gsMatrix<T> & solVector,
                                                   const bool Matrix, const bool homogenize)
{
    gsMultiPatch<T> def;
    constructSolution(solVector, def);
    return assemble(def,Matrix,homogenize);
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::_constructSolution(const gsMatrix<T> & solVector, const gsMultiPatch<T> & /*undeformed*/) const
{
    gsMultiPatch<T> mp = m_patches;
    gsMultiPatch<T> displacement = constructDisplacement(solVector);
    for ( size_t k =0; k!=displacement.nPatches(); ++k) // Deform the geometry
        mp.patch(k).coefs() += displacement.patch(k).coefs();;  // defG points to mp_def, therefore updated

    return mp;
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::constructSolution(const gsMatrix<T> & solVector) const
{
    return _constructSolution(solVector,m_patches);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = _constructSolution(solVector,m_patches);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::updateMultiPatch(const gsMatrix<T> & solVector, gsMultiPatch<T> & mp) const
{
    mp = _constructSolution(solVector,mp);
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::getArea(const gsFunctionSet<T> & geometry)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap G = m_assembler.getMap(geometry);

    gsExprEvaluator<T> evaluator(m_assembler);
    T result = evaluator.integral(meas(G));
    return result;
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::getDisplacementNorm(const gsFunctionSet<T> & deformed)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    auto u   = m_def - m_ori;

    gsExprEvaluator<T> evaluator(m_assembler);
    T result = evaluator.integral( u.tr() * u * meas(m_def));
    T area = evaluator.integral(meas(m_ori));

    return std::pow(result/area,0.5);
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::getElasticEnergy(const gsFunctionSet<T> & deformed)
{
    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(deformed);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,&deformed);
    auto S0  = m_assembler.getCoeff(m_S0);
    auto S1  = m_assembler.getCoeff(m_S1);
    auto u   = m_def - m_ori;

    auto m_N        = S0.tr();
    auto m_M        = S1.tr(); // output is a column

    gsExprEvaluator<T> evaluator(m_assembler);
    T result = evaluator.integral(0.5 * ( u.tr() * ( m_N + m_M ).tr() ) * meas(m_def));
    return result;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::plotSolution(std::string string, const gsMatrix<T> & solVector)
{
    m_solvector = solVector;
    space m_space = m_assembler.trialSpace(0);
    solution m_solution = m_assembler.getSolution(m_space, m_solvector);
    geometryMap G = m_assembler.getMap(m_patches);
    gsExprEvaluator<T> ev(m_assembler);
    ev.options().setSwitch("plot.elements", false);
    ev.options().setInt   ("plot.npts"    , 500);
    ev.writeParaview( m_solution, G, string);
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::interfaceErrorC0(const gsFunctionSet<T> & deformed, const ifContainer & iFaces)
{
    geometryMap G = m_assembler.getMap(deformed);
    gsExprEvaluator<T> ev(m_assembler);
    ev.integralInterface( ( G.left() - G.right() ).sqNorm() , iFaces);
    ev.calcSqrt();
    return ev.value();
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::interfaceErrorG1(const gsFunctionSet<T> & deformed, const ifContainer & iFaces)
{
    geometryMap G = m_assembler.getMap(deformed);
    gsExprEvaluator<T> ev(m_assembler);
    ev.integralInterface( (sn(G.left()).normalized()-sn(G.right()).normalized()).sqNorm() , iFaces);
    ev.calcSqrt();
    return ev.value();
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::interfaceErrorNormal(const gsFunctionSet<T> & deformed, const ifContainer & iFaces)
{
    geometryMap G = m_assembler.getMap(deformed);
    gsExprEvaluator<T> ev(m_assembler);
    ev.maxInterface( (sn(G.left())-sn(G.right())).norm() , iFaces);
    return ev.value();
}

template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::interfaceErrorGaussCurvature(const gsFunctionSet<T> & deformed, const ifContainer & iFaces)
{
    geometryMap G = m_assembler.getMap(deformed);
    gsExprEvaluator<T> ev(m_assembler);
    ev.maxInterface( abs( (fform(G.left() ).inv()*fform2nd(G.left() )).det() -
                          (fform(G.right()).inv()*fform2nd(G.right())).det() ) , iFaces);
    return ev.value();
}
template<short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::interfaceErrorMeanCurvature(const gsFunctionSet<T> & deformed, const ifContainer & iFaces)
{
    geometryMap G = m_assembler.getMap(deformed);
    gsExprEvaluator<T> ev(m_assembler);
    ev.maxInterface( abs( (fform(G.left() ).inv()*fform2nd(G.left() )).trace().val() -
                          (fform(G.right()).inv()*fform2nd(G.right())).trace().val() ) , iFaces);
    return ev.value();
}
template<short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::constructMultiPatch(const gsMatrix<T> & solVector) const
{
    m_solvector = solVector;
    space m_space = m_assembler.trialSpace(0);
    m_space.setup(m_bcs, dirichlet::l2Projection, m_continuity);
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart() = m_ddofs; //CHECK FIXEDPART

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
                    {
                        cc(i,c) =  m_ddofs.at( m_mapper.global_to_bindex(ii) );
                    }
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
            result.addPatch(m_basis.basis(p).makeGeometry( give(cc) ));  // defG points to mp_def, therefore updated
        }
        return result;
    }
}

template<short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::constructDisplacement(const gsMatrix<T> & solVector) const
{
    return constructMultiPatch(solVector);
}

template<short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::fullSolutionVector(const gsMatrix<T> & vector) const
{
    gsMatrix<T> solVector = vector;
    space m_space = m_assembler.trialSpace(0);
    m_space.setup(m_bcs, dirichlet::l2Projection, m_continuity);
    solution m_solution = m_assembler.getSolution(m_space, solVector);
    gsMatrix<T> result;
    m_solution.extractFull(result);
    return result.col(0);
}

template<short_t d, class T, bool bending>
gsVector<T> gsThinShellAssembler<d, T, bending>::constructSolutionVector(const gsMultiPatch<T> & displacements) const
{
    gsVector<T> result(m_mapper.freeSize());

    for (size_t p=0; p!=displacements.nPatches(); p++)
    {
        for (size_t dim = 0; dim!=d; dim++)
        {
            for (size_t k=0; k!=m_mapper.patchSize(p,dim); k++)
            {
                if (m_mapper.is_free(k,p,dim))
                {
                    result.at(m_mapper.index(k,p,dim)) = displacements.patch(p).coefs()(k,dim);
                }
            }
        }

    }
    return result;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructDisplacement(solVector);
}

// template<short_t d, class T, bool bending>
// void gsThinShellAssembler<d, T, bending>::constructStresses(const gsMultiPatch<T> & deformed,
//                                                     gsPiecewiseFunction<T> & result,
//                                                     stress_type::type type) const
// {
//     deformed = constructDisplacement(solVector);
// }

template<short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::computePrincipalStretches(const gsMatrix<T> & u, const gsFunctionSet<T> & deformed, const T z)
{
    // gsDebug<<"Warning: Principle Stretch computation of gsThinShellAssembler is depreciated...\n";
    gsMatrix<T> Z(1,1);
    Z.setZero();
    gsMatrix<T> result(3,u.cols());
    result.setZero();
    gsMatrix<T> zmat(1,1);
    zmat<<z;
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // geometryMap m_ori   = m_assembler.getMap(m_patches);
    // geometryMap m_def   = m_assembler.getMap(*m_defpatches);
    // m_assembler.initSystem();

    gsMaterialMatrixEval<T,MaterialOutput::Stretch> m_mm(m_materialMatrices,&deformed,zmat);
    auto mm0 = m_assembler.getCoeff(m_mm);

    gsExprEvaluator<T> evaluator(m_assembler);

    for (index_t k = 0; k != u.cols(); ++k)
        result.col(k) = evaluator.eval(mm0,u.col(k));
    return result;
}

template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::computePrincipalStresses(const gsMatrix<T> & u, const gsFunctionSet<T> & deformed, const T z)
{
    // gsDebug<<"Warning: Principle Stretch computation of gsThinShellAssembler is depreciated...\n";
    gsMatrix<T> Z(1,1);
    Z.setZero();
    gsMatrix<T> result(2,u.cols());
    result.setZero();
    gsMatrix<T> zmat(1,1);
    zmat<<z;
    this->_getOptions();

    m_assembler.cleanUp();
    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // geometryMap m_ori   = m_assembler.getMap(m_patches);
    // geometryMap m_def   = m_assembler.getMap(*m_defpatches);
    // m_assembler.initSystem();

    gsMaterialMatrixEval<T,MaterialOutput::PStress> m_mm(m_materialMatrices,&deformed,zmat);
    auto mm0 = m_assembler.getCoeff(m_mm);

    gsExprEvaluator<T> evaluator(m_assembler);

    for (index_t k = 0; k != u.cols(); ++k)
    {
        result.col(k) = evaluator.eval(mm0,u.col(k));
    }
    return result;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructStress(const gsFunctionSet<T> & deformed,
                                                    gsPiecewiseFunction<T> & result,
                                                    stress_type::type type)
{
    constructStress(m_patches,deformed,result,type);
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructStress(
                                                    const gsFunctionSet<T> & original,
                                                    const gsFunctionSet<T> & deformed,
                                                    gsPiecewiseFunction<T> & result,
                                                    stress_type::type type)
{
    result.clear();

    for (size_t p = 0; p < m_patches.nPatches(); ++p )
        result.addPiecePointer(new gsShellStressFunction<T>(original,deformed,m_materialMatrices,p,type));

}

// template<short_t d, class T, bool bending>
// gsField<T> gsThinShellAssembler<d, T, bending>::constructStress(const gsMultiPatch<T> & deformed,
//                                                     stress_type::type type)
// {
//     gsPiecewiseFunction<T> result;
//     result.clear();

//     for (size_t p = 0; p < m_patches.nPatches(); ++p )
//         result.addPiecePointer(new gsShellStressFunction<d, T, bending>(m_patches,deformed,m_materialMatrices,p,type,m_assembler));

//     gsField<T> stressField(m_patches,result, true);
//     return stressField;

// }

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::projectL2_into(const gsFunction<T> & /*fun*/, gsMatrix<T>& /*result*/)
{
    // /// todo: make a projection with BCs?
    // /// todo: test
    // // this->_getOptions();

    // m_assembler.cleanUp();
    // GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    // m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    // geometryMap m_ori   = m_assembler.getMap(m_patches);

    // // Initialize stystem
    // m_assembler.initSystem();

    // space       m_space = m_assembler.trialSpace(0);
    // auto    function = m_assembler.getCoeff(fun, m_ori);
    // // auto    function = m_assembler.getCoeff(fun);

    // // assemble system
    // m_assembler.assemble(m_space*m_space.tr()*meas(m_ori),m_space * function*meas(m_ori));

    // gsSparseSolver<>::uPtr solver = gsSparseSolver<T>::get( m_options.askString("Solver","CGDiagonal") );
    // solver->compute(m_assembler.matrix());
    // result = solver->solve(m_assembler.rhs());
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::projectL2_into(const gsFunction<T> & fun, gsMultiPatch<T>& mp)
{
    /// todo: make a projection with BCs?
    /// todo: test
    gsMatrix<T> tmp = projectL2(fun);
    mp = m_patches;

    // Solution vector and solution variable
    space m_space = m_assembler.trialSpace(0);
    m_space.setup(m_bcs, dirichlet::l2Projection, m_continuity);
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart() = m_ddofs;

    solution m_solution = m_assembler.getSolution(m_space, tmp);

    gsMatrix<T> cc;
    for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
    {
        // // extract deformed geometry
        m_solution.extract(cc, k);
        mp.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }
}


template<short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::projectL2(const gsFunction<T> & fun)
{
    /// todo: make a projection with BCs?
    /// todo: test
    gsMatrix<T> result;
    this->projectL2_into(fun,result);
    return result;
}

template <short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::deformationNorm(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> & original)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    geometryMap ori = m_assembler.getMap(original);
    geometryMap def = m_assembler.getMap(deformed);

    gsExprEvaluator<T> evaluator(m_assembler);
    T result = evaluator.integral(def.sqNorm() * meas(ori));
    return result;
}

template<short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_ifcTest(const T tol)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    gsExprEvaluator<T> ev(m_assembler);

    m_inPlane.clear();
    m_outPlane.clear();

    for (gsBoxTopology::const_iiterator it = m_patches.topology().iBegin(); it!=m_patches.topology().iEnd(); it++)
    {
        // G1 condition
        ev.integralInterface( (sn(m_ori.left()).normalized()-sn(m_ori.right()).normalized()).sqNorm() );
        ev.calcSqrt();

        // // Continuous normal condition
        // ev.maxInterface( (sn(m_ori.left())-sn(m_ori.right())).norm() );
        // ev.calcSqrt();

        if (ev.value() < tol)
            m_inPlane.push_back(*it);
        else
            m_outPlane.push_back(*it);
    }
}

template<short_t d, class T, bool bending>
bool gsThinShellAssembler<d, T, bending>::_isInPlane(const boundaryInterface & /*ifc*/, const T tol)
{
    geometryMap m_ori   = m_assembler.getMap(m_patches);
    gsExprEvaluator<T> ev(m_assembler);

    // G1 condition
    ev.integralInterface( (sn(m_ori.left()).normalized()-sn(m_ori.right()).normalized()).sqNorm() );
    ev.calcSqrt();

    // // Continuous normal condition
    // ev.maxInterface( (sn(m_ori.left())-sn(m_ori.right())).norm() );
    // ev.calcSqrt();

    return (ev.value() < tol);
}


}// namespace gismo
