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
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>

#include <gsPde/gsBoundaryConditions.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsPiecewiseFunction.h>

#include <gsIO/gsWriteParaview.h>

#include <gsMSplines/gsWeightMapper.h>

#include <unordered_set>

namespace gismo
{

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>::gsThinShellAssembler(const gsMultiPatch<T> & patches,
                                                          const gsMultiBasis<T> & basis,
                                                          const gsBoundaryConditions<T> & bconditions,
                                                          const gsFunction<T> & surface_force,
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
    this->_defaultOptions();
    this->_getOptions();
    this->_initialize();
}

template<short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>::gsThinShellAssembler(const gsMultiPatch<T> & patches,
                                                          const gsMultiBasis<T> & basis,
                                                          const gsBoundaryConditions<T> & bconditions,
                                                          const gsFunction<T> & surface_force,
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
gsThinShellAssembler<d, T, bending>::_assemblePressure_impl(const gsFunction<T> & , const gsFunctionSet<T> & )
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
    // No matrix contribution for the linear case
}

template <short_t d, class T, bool bending>
template <short_t _d, bool _matrix>
typename std::enable_if<(_d==3) && !_matrix, void>::type
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

template <short_t d, class T, bool bending>
template <short_t _d, bool matrix>
typename std::enable_if<(_d==3) && matrix, void>::type
gsThinShellAssembler<d, T, bending>::_assembleFoundation_impl(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed)
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

        // Compute actives and values of basis functions on point load location.
        if ( m_pMass[i].parametric )   // in parametric space
        {
            m_basis.front().basis(m_pMass[i].patch).active_into( m_pMass[i].point, acts );
            m_basis.front().basis(m_pMass[i].patch).eval_into  ( m_pMass[i].point, bVals);
        }
        else                            // in physical space
        {
            gsMatrix<T> forcePoint;
            m_patches.patch(m_pMass[i].patch).invertPoints(m_pMass[i].point,forcePoint);
            m_basis.front().basis(m_pMass[i].patch).active_into( forcePoint, acts );
            m_basis.front().basis(m_pMass[i].patch).eval_into  ( forcePoint, bVals);
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
    m_space.setup(m_bcs, dirichlet::homogeneous, m_continuity);
    // this->homogenizeDirichlet();

    gsExprEvaluator<T> ev(m_assembler);
    gsVector<> pt(2);
    pt.setConstant(.25);

    try
    {
        // assemble system
        if (!lumped)
            m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
        else
            m_assembler.assemble(mm0.val()*(m_space.rowSum())*meas(m_ori));

/*        // assemble system
        if (!lumped)
        {
            m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
            m_mass = m_assembler.matrix();
            this->_applyMass();
        }
        else
        {
            // To do: add point masses in lumped case
            m_assembler.assemble(mm0.val()*(m_space.rowSum())*meas(m_ori));
            m_rhs = m_assembler.rhs();
        }
        m_mass = m_assembler.matrix();

        this->_applyMass();
        m_status = ThinShellAssemblerStatus::Success;*/
    }
    catch (...)
    {
        m_assembler.cleanUp();
        m_status = ThinShellAssemblerStatus::AssemblyError;
    }
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

    auto m_force = m_assembler.getCoeff(*m_forceFun, m_ori);


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
            ,
            m_space * m_force * meas(m_ori)
            );

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
    auto m_force = m_assembler.getCoeff(*m_forceFun, m_ori);

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
            ,
            m_space * m_force * meas(m_ori)
            );

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

    gsVector<T> pt(2);
    pt.setConstant(0.25);
    gsExprEvaluator<T> ev(m_assembler);

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
    auto m_force = m_assembler.getCoeff(*m_forceFun, m_ori);


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
        m_assembler.assemble(m_space * m_force * meas(m_ori) -
                                ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
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
    auto m_force = m_assembler.getCoeff(*m_forceFun, m_ori);

    if (homogenize) this->homogenizeDirichlet();
    else            this->_assembleDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    try
    {
        if (m_foundInd) this->_assembleFoundation<false>(*m_foundFun,deformed);
        if (m_pressInd) this->_assemblePressure<false>(*m_pressFun,deformed);

        // Assemble vector
        m_assembler.assemble(m_space * m_force * meas(m_ori) -
                    ( ( m_N * m_Em_der.tr() ) * meas(m_ori) ).tr()
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
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::_constructSolution(const gsMatrix<T> & solVector, const gsMultiPatch<T> & undeformed) const
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
void gsThinShellAssembler<d, T, bending>::projectL2_into(const gsFunction<T> & fun, gsMatrix<T>& result)
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
bool gsThinShellAssembler<d, T, bending>::_isInPlane(const boundaryInterface & ifc, const T tol)
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
