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

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>

namespace gismo
{

template <short_t d, class T, bool bending>
gsThinShellAssembler<d, T, bending>::gsThinShellAssembler(const gsMultiPatch<T> & patches,
                                                          const gsMultiBasis<T> & basis,
                                                          const gsBoundaryConditions<T> & bconditions,
                                                          const gsFunction<T> & surface_force,
                                                          gsMaterialMatrixBase<T> * materialmatrix
                                                          )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_bcs(bconditions),
                                        m_forceFun(&surface_force),
                                        m_materialMat(materialmatrix)
{
    this->_defaultOptions();
    this->_initialize();
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_defaultOptions()
{
    m_options.addReal("WeakDirichlet","Penalty parameter weak dirichlet conditions",1e3);
    m_options.addReal("WeakClamped","Penalty parameter weak clamped conditions",1e3);

    // Assembler options
    m_options.addInt("DirichletStrategy","Method for enforcement of Dirichlet BCs [11..14]",11);
    m_options.addInt("DirichletValues","Method for computation of Dirichlet DoF values [100..103]",101);
    m_options.addInt("InterfaceStrategy","Method of treatment of patch interfaces [0..3]",1);
    m_options.addReal("bdA","Estimated nonzeros per column of the matrix: bdA*deg + bdB",2);
    m_options.addInt("bdB","Estimated nonzeros per column of the matrix: bdA*deg + bdB",1);
    m_options.addReal("bdO","Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]",0.333);
    m_options.addReal("quA","Number of quadrature points: quA*deg + quB",1);
    m_options.addInt("quB","Number of quadrature points: quA*deg + quB",1);
    m_options.addInt("quRule","Quadrature rule [1:GaussLegendre,2:GaussLobatto]",1);
}


template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_getOptions() const
{
    m_alpha_d = m_options.getReal("WeakDirichlet");
    m_alpha_r = m_options.getReal("WeakClamped");
}


template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_initialize()
{
    //gsInfo<<"Active options:\n"<< m_assembler.options() <<"\n";
    m_defpatches = m_patches;

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);
    m_assembler.setOptions(m_options);

    // Initialize the geometry maps
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Set the discretization space
    space m_space = m_assembler.getSpace(m_basis, d, 0); // last argument is the space ID

    this->_assembleDirichlet();

    m_ddofs = m_space.fixedPart();
    m_mapper = m_space.mapper();

    // foundation is off by default
    m_foundInd = false;
    // pressure is off by default
    m_pressInd = false;

    GISMO_ASSERT(m_forceFun->targetDim()==d,"Force must have " << d<<" dimensions but has "<<m_forceFun->targetDim());
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_assembleNeumann()
{
    _assembleNeumann_impl<d,bending>();
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssembler<d, T, bending>::_assembleNeumann_impl()
{
    m_assembler.getMap(m_patches);           // this map is used for integrals
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();
    m_assembler.assembleRhsBc(m_space * g_N * tv(m_ori).norm(), m_bcs.neumannSides() );
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssembler<d, T, bending>::_assembleNeumann_impl()
{
    m_assembler.getMap(m_patches);           // this map is used for integrals
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();
    m_assembler.assembleRhsBc(m_space * g_N * tv(m_ori).norm(), m_bcs.neumannSides() );
}


template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_assembleWeakBCs()
{
    this->_getOptions();
    _assembleWeakBCs_impl<d,bending>();
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl()
{
    m_assembler.getMap(m_patches);           // this map is used for integrals

    geometryMap m_ori   = m_assembler.exprData()->getMap();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();

    // Weak BCs
    m_assembler.assembleLhsRhsBc
    (
        -(m_alpha_d * m_space * m_space.tr()) * meas(m_ori)
        ,
        (m_alpha_d * (m_space * (m_ori - m_ori) - m_space * (g_N) )) * meas(m_ori)
        ,
        m_bcs.container("Weak Dirichlet")
    );

    // for weak clamped
    m_assembler.assembleLhsRhsBc
    (
        (
            m_alpha_r * ( sn(m_ori).tr()*nv(m_ori) - sn(m_ori).tr()*nv(m_ori) ).val() * ( var2(m_space,m_space,m_ori,nv(m_ori).tr()) )
            +
            m_alpha_r * ( ( var1(m_space,m_ori) * nv(m_ori) ) * ( var1(m_space,m_ori) * nv(m_ori) ).tr() )
        ) * meas(m_ori)
        ,
        // THIS LINE SHOULD BE ZERO!
        ( m_alpha_r * ( sn(m_ori).tr()*sn(m_ori) - sn(m_ori).tr()*sn(m_ori) ).val() * ( var1(m_space,m_ori) * sn(m_ori) ) ) * meas(m_ori)
        ,
        m_bcs.container("Weak Clamped")
    );
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl()
{
    m_assembler.getMap(m_patches);           // this map is used for integrals

    geometryMap m_ori   = m_assembler.exprData()->getMap();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();

    // Weak BCs
    m_assembler.assembleLhsRhsBc
    (
        -(m_alpha_d * m_space * m_space.tr()) * meas(m_ori)
        ,
        (m_alpha_d * (m_space * (m_ori - m_ori) - m_space * (g_N) )) * meas(m_ori)
        ,
        m_bcs.container("Weak Dirichlet")
    );
}


template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_assembleWeakBCs(const gsMultiPatch<T> & deformed)
{
    this->_getOptions();
    _assembleWeakBCs_impl<d,bending>(deformed);
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl(const gsMultiPatch<T> & deformed)
{
    m_defpatches = deformed;
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();

    // Weak BCs
    m_assembler.assembleLhsRhsBc
    (
        -m_alpha_d * m_space * m_space.tr() * meas(m_ori)
        ,
        m_alpha_d * (m_space * (m_def - m_ori) - m_space * (g_N) ) * meas(m_ori)
        ,
        m_bcs.container("Weak Dirichlet")
    );

    // for weak clamped
    m_assembler.assembleLhsRhsBc
    (
        (
            m_alpha_r * ( sn(m_def).tr()*nv(m_ori) - sn(m_ori).tr()*nv(m_ori) ).val() * ( var2(m_space,m_space,m_def,nv(m_ori).tr()) )
            +
            m_alpha_r * ( ( var1(m_space,m_def) * nv(m_ori) ) * ( var1(m_space,m_def) * nv(m_ori) ).tr() )
        ) * meas(m_ori)
        ,
        ( m_alpha_r * ( sn(m_def).tr()*sn(m_ori) - sn(m_ori).tr()*sn(m_ori) ).val() * ( var1(m_space,m_def) * sn(m_ori) ) ) * meas(m_ori)
        ,
        m_bcs.container("Weak Clamped")
    );
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssembler<d, T, bending>::_assembleWeakBCs_impl(const gsMultiPatch<T> & deformed)
{
    m_defpatches = deformed;
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();

    // Weak BCs
    m_assembler.assembleLhsRhsBc
    (
        -m_alpha_d * m_space * m_space.tr() * meas(m_ori)
        ,
        m_alpha_d * (m_space * (m_def - m_ori) - m_space * (g_N) ) * meas(m_ori)
        ,
        m_bcs.container("Weak Dirichlet")
    );
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_assembleDirichlet()
{
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // if statement
    m_space.setup(m_bcs, dirichlet::interpolation, 0);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::homogenizeDirichlet()
{
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    m_space.setup(m_bcs, dirichlet::homogeneous, 0);
    // space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart().setZero();
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::_applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts,globalActs;

    space       m_space = m_assembler.trialSpace(0);
    m_mapper = m_space.mapper();

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if (m_pLoads[i].value.size()!=d)
            gsWarn<<"Point load has wrong dimension "<<m_pLoads[i].value.size()<<" instead of "<<d<<"\n";
        // Compute actives and values of basis functions on point load location.
        if ( m_pLoads[i].parametric )   // in parametric space
        {
            m_basis.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_basis.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else                            // in physical space
        {
            gsMatrix<T> forcePoint;
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);
            m_basis.front().basis(m_pLoads[i].patch).active_into( forcePoint, acts );
            m_basis.front().basis(m_pLoads[i].patch).eval_into  ( forcePoint, bVals);
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

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assembleMass(bool lumped)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    m_assembler.getMap(m_patches);           // this map is used for integrals

    // Initialize stystem
    m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,m_patches);
    variable mm0 = m_assembler.getCoeff(m_mm);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    // assemble system
    if (!lumped)
        m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
    else
        m_assembler.assemble(mm0.val()*(m_space.rowSum())*meas(m_ori));
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assembleFoundation()
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    m_assembler.getMap(m_patches);           // this map is used for integrals
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    // Initialize stystem
    m_assembler.initSystem(false);
    variable    m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
    GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

    space       m_space = m_assembler.trialSpace(0);

    m_assembler.assemble(m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori));
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assemble()
{
    assemble_impl<d, bending>();
}

/**
    @brief Assembles the Kirchhoff-Love shell equations including the bending terms.
    Optionally, pressure is included via \a p * n * u
    Optionally, foundation stiffness is included via k_x v_x v_x + k_y v_y v_y + k_z v_z v_z
    Since the variational energy of the foundation force k_i u_i is equal to k_i u_i v_i where i denotes any direction, u_i are displacemets and v_i are spaces.

*/
template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssembler<d, T, bending>::assemble_impl()
{
    // this->_getOptions();

    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize stystem
    m_assembler.initSystem(false);
    m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> m_mmB(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> m_mmC(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMat,m_defpatches);
    variable mmA = m_assembler.getCoeff(m_mmA);
    variable mmB = m_assembler.getCoeff(m_mmB);
    variable mmC = m_assembler.getCoeff(m_mmC);
    variable mmD = m_assembler.getCoeff(m_mmD);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);

    /*
        NOTE: after update of the exprParser, we can use
        geometryMap m_ori   = m_assembler.getMap(m_patches);
        geometryMap m_def   = m_assembler.getMap(m_defpatches);
     */
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);


    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) );
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3);

    auto m_N_der    = m_Em_der * reshape(mmA,3,3) + m_Ef_der * reshape(mmB,3,3);
    auto m_M_der    = m_Em_der * reshape(mmC,3,3) + m_Ef_der * reshape(mmD,3,3);

    if (m_foundInd)
    {
        variable m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        m_assembler.assemble(
            m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori)
            );

    }
    if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        m_assembler.assemble(
            m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
            );
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

    this->_assembleWeakBCs();
    this->_assembleNeumann();

    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        _applyLoads();
    }
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssembler<d, T, bending>::assemble_impl()
{
    // this->_getOptions();

    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize stystem
    m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMat,m_defpatches);
    variable mmA = m_assembler.getCoeff(m_mmA);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);

    auto jacG       = jac(m_def);
    auto m_Em_der   = flat( jacG.tr() * jac(m_space) ) ; //[checked]
    auto m_N_der    = m_Em_der * reshape(mmA,3,3);

    if (m_foundInd)
    {
        variable m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        m_assembler.assemble(
            m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori)
            );
    }
    if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        // Assemble vector
        m_assembler.assemble(
            m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
            );
    }

    m_assembler.assemble(
        (
            m_N_der * m_Em_der.tr()
        ) * meas(m_ori)
        ,
        m_space * m_force * meas(m_ori)
        );

    this->_assembleWeakBCs();
    this->_assembleNeumann();

    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        _applyLoads();
    }
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assembleMatrix(const gsMultiPatch<T> & deformed)
{
    assembleMatrix_impl<d, bending>(deformed);
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssembler<d, T, bending>::assembleMatrix_impl(const gsMultiPatch<T> & deformed)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize matrix
    m_assembler.initMatrix(false);
    // m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> m_mmB(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> m_mmC(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> m_mmD(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMat,m_defpatches);
    variable mmA = m_assembler.getCoeff(m_mmA);
    variable mmB = m_assembler.getCoeff(m_mmB);
    variable mmC = m_assembler.getCoeff(m_mmC);
    variable mmD = m_assembler.getCoeff(m_mmD);
    variable S0  = m_assembler.getCoeff(m_S0);
    variable S1  = m_assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();


    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]


    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Ef_der2  = -(flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
                            + var2(m_space,m_space,m_def, m_M ));

    auto m_N_der    = m_Em_der * reshape(mmA,3,3) + m_Ef_der * reshape(mmB,3,3);
    auto m_M_der    = m_Em_der * reshape(mmC,3,3) + m_Ef_der * reshape(mmD,3,3);

    if (m_foundInd)
    {
        variable m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        m_assembler.assemble(
                m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori)
            );
    }
    if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        m_assembler.assemble(
                                -m_pressure.val() * m_space * var1(m_space,m_def).tr()* meas(m_ori)
                            );
    }
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
    this->_assembleWeakBCs(deformed);
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssembler<d, T, bending>::assembleMatrix_impl(const gsMultiPatch<T> & deformed)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize matrix
    m_assembler.initMatrix(false);
    // m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> m_mmA(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    variable mmA = m_assembler.getCoeff(m_mmA);
    variable S0  = m_assembler.getCoeff(m_S0);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();


    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_N_der    = m_Em_der * reshape(mmA,3,3);

    // Assemble matrix
    m_assembler.assemble(
            (
                m_N_der * m_Em_der.tr()
                +
                m_Em_der2
            ) * meas(m_ori)
        );

    if (m_foundInd)
    {
        variable m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        m_assembler.assemble(
                m_space * m_foundation.asDiag() * m_space.tr() * meas(m_ori)
            );
    }
    if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        m_assembler.assemble(
                -m_pressure.val() * m_space * var1(m_space,m_def).tr()* meas(m_ori)
            );
    }
    this->_assembleWeakBCs(deformed);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assembleMatrix(const gsMatrix<T> & solVector)
{
    // gsMultiPatch<T> deformed;
    // constructSolution(solVector, deformed);
    // assembleMatrix(deformed);

    constructSolution(solVector, m_defpatches);
    assembleMatrix(m_defpatches);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assembleVector(const gsMultiPatch<T> & deformed)
{
  assembleVector_impl<d, bending>(deformed);
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssembler<d, T, bending>::assembleVector_impl(const gsMultiPatch<T> & deformed)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize vector
    m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMat,m_defpatches);
    variable S0  = m_assembler.getCoeff(m_S0);
    variable S1  = m_assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space m_space       = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);


    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    if (m_foundInd)
    {
        variable m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        // Assemble vector
        m_assembler.assemble(
                      m_space * m_foundation.asDiag() * (m_def - m_ori) * meas(m_ori) // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
                    );
    }
    if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        // Assemble vector
        m_assembler.assemble(
                      m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
                      );
    }

        // Assemble vector
    m_assembler.assemble(m_space * m_force * meas(m_ori) -
                            ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                        );

    this->_assembleWeakBCs(deformed);
    this->_assembleNeumann();

    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        _applyLoads();
    }
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssembler<d, T, bending>::assembleVector_impl(const gsMultiPatch<T> & deformed)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize vector
    m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    variable S0  = m_assembler.getCoeff(m_S0);

    space m_space       = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);

    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    if (m_foundInd)
    {
        variable m_foundation = m_assembler.getCoeff(*m_foundFun, m_ori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        // Assemble vector
        m_assembler.assemble(
                      m_space * m_foundation.asDiag() * (m_def - m_ori) * meas(m_ori) // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
                    );
    }
    if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        // Assemble vector
        m_assembler.assemble(
                      m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
                    );
    }

    // Assemble vector
    m_assembler.assemble(m_space * m_force * meas(m_ori) -
                ( ( m_N * m_Em_der.tr() ) * meas(m_ori) ).tr()
                );

    this->_assembleWeakBCs(deformed);
    this->_assembleNeumann();

    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        _applyLoads();
    }
}

template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::boundaryForce(const gsMultiPatch<T> & deformed, patchSide& ps)
{
    return boundaryForce_impl<d, bending>(deformed,ps);
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, gsMatrix<T> >::type
gsThinShellAssembler<d, T, bending>::boundaryForce_impl(const gsMultiPatch<T> & deformed, patchSide& ps)
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(m_basis, d, 0); // last argument is the space ID

    assembler.initSystem();

    m_defpatches = deformed;

    assembler.getMap(m_patches);           // this map is used for integrals
    assembler.getMap(m_defpatches);

    // Initialize vector
    // m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMat,m_defpatches);
    variable S0  = assembler.getCoeff(m_S0);
    variable S1  = assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = assembler.getCoeff(mult2t);

    geometryMap m_ori   = assembler.exprData()->getMap();
    geometryMap m_def   = assembler.exprData()->getMap2();


    // this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(u) ) ;

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(u,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(u,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    // Assemble vector
    assembler.assemble(
                  - ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                );

    gsMatrix<T> Fint = assembler.rhs();
    gsVector<T> result(d);
    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>*>(&u.source());
    gsMatrix<index_t> boundary;

    gsBoundaryConditions<real_t>::bcContainer container;
    m_bcs.getConditionFromSide(ps,container);

    for ( index_t com = 0; com!=d; com++)
    {
        const gsBasis<T> & basis = mbasis[ps.patch];
        boundary = basis.boundary(ps.side());

        T offset = u.mapper().offset(ps.patch);
        T size = u.mapper().size(com);
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            index_t ii = offset + size*com + boundary.at(l);
            result.at(com) += Fint.at(ii);
        }
    }

    return result;
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), gsMatrix<T> >::type
gsThinShellAssembler<d, T, bending>::boundaryForce_impl(const gsMultiPatch<T> & deformed, patchSide& ps)
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(m_basis, d, 0); // last argument is the space ID

    assembler.initSystem();

    m_defpatches = deformed;

    assembler.getMap(m_patches);           // this map is used for integrals
    assembler.getMap(m_defpatches);

    // Initialize vector
    // m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    variable S0  = assembler.getCoeff(m_S0);

    geometryMap m_ori   = assembler.exprData()->getMap();
    geometryMap m_def   = assembler.exprData()->getMap2();


    // this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(u) ) ;

    // Assemble vector
    assembler.assemble(
                  - ( ( m_N * m_Em_der.tr() ) * meas(m_ori) ).tr()
                );

    gsMatrix<T> Fint = assembler.rhs();
    gsVector<T> result(d);
    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>*>(&u.source());
    gsMatrix<index_t> boundary;


    for ( index_t com = 0; com!=d; com++)
    {
        const gsBasis<T> & basis = mbasis[ps.patch];
        boundary = basis.boundary(ps.side());

        T offset = u.mapper().offset(ps.patch);
        T size = u.mapper().size(com);
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            index_t ii = offset + size*com + boundary.at(l);
            result.at(com) += Fint.at(ii);
        }
    }

    return result;
}

template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::boundaryForceVector(const gsMultiPatch<T> & deformed, patchSide& ps, index_t com)
{
    return boundaryForceVector_impl<d, bending>(deformed,ps,com);
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, gsMatrix<T> >::type
gsThinShellAssembler<d, T, bending>::boundaryForceVector_impl(const gsMultiPatch<T> & deformed, patchSide& ps, index_t com)
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(m_basis, d, 0); // last argument is the space ID

    assembler.initSystem();

    m_defpatches = deformed;

    assembler.getMap(m_patches);           // this map is used for integrals
    assembler.getMap(m_defpatches);

    // Initialize vector
    // m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMat,m_defpatches);
    variable S0  = assembler.getCoeff(m_S0);
    variable S1  = assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = assembler.getCoeff(mult2t);

    geometryMap m_ori   = assembler.exprData()->getMap();
    geometryMap m_def   = assembler.exprData()->getMap2();


    // this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(u) ) ;

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(u,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(u,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    // Assemble vector
    assembler.assemble(
                  - ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                );

    gsMatrix<T> Fint = assembler.rhs();
    gsMatrix<T> result;
    const gsMultiBasis<T> & mbasis =
        *dynamic_cast<const gsMultiBasis<T>*>(&u.source());
    gsMatrix<index_t> boundary;

    typedef gsBoundaryConditions<T> bcList;

    gsBoundaryConditions<real_t>::bcContainer container;
    m_bcs.getConditionFromSide(ps,container);

    // If the boundary is not a dirichlet boundary, then there are forces applied????
    for ( typename bcList::const_iterator it =  container.begin();
          it != container.end() ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;
        if( (it->unkComponent()!=com) && (it->unkComponent()!=-1) ) continue;

        const int k = it->patch();
        const gsBasis<T> & basis = mbasis[k];

        // Get dofs on this boundary
        boundary = basis.boundary(it->side());
        result.resize(boundary.size(),1);

        T offset = u.mapper().offset(k);
        T size = u.mapper().size(com);
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            index_t ii = offset + size*com + boundary.at(l);
            result.at(l) = Fint.at(ii);
        }
    }

    // gsInfo<<"Neumann forces\n";
    // for ( typename bcList::const_iterator it =  m_bcs.begin("Neumann");
    //       it != m_bcs.end("Neumann") ; ++it )
    // {
    //     if( it->unknown()!=u.id() ) continue;
    //     //
    //     for (index_t r = 0; r!=u.dim(); ++r)
    //     {
    //         const int k = it->patch();
    //         const gsBasis<T> & basis = mbasis[k];

    //         // Get dofs on this boundary
    //         boundary = basis.boundary(it->side());

    //         T offset = u.mapper().offset(0);
    //         T size = u.mapper().size(r);

    //         for (index_t l=0; l!= boundary.size(); ++l)
    //         {
    //             index_t ii = offset + size*r + boundary.at(l);
    //             gsInfo<<result.at(ii)<<"\n";
    //         }
    //     }
    // }
    return result;
}

template<int d, typename T, bool bending>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), gsMatrix<T> >::type
gsThinShellAssembler<d, T, bending>::boundaryForceVector_impl(const gsMultiPatch<T> & deformed, patchSide& ps, index_t com)
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(m_basis, d, 0); // last argument is the space ID

    assembler.initSystem();

    m_defpatches = deformed;

    assembler.getMap(m_patches);           // this map is used for integrals
    assembler.getMap(m_defpatches);

    // Initialize vector
    // m_assembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMat,m_defpatches);
    variable S0  = assembler.getCoeff(m_S0);

    geometryMap m_ori   = assembler.exprData()->getMap();
    geometryMap m_def   = assembler.exprData()->getMap2();

    // this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(u) ) ;

    // Assemble vector
    assembler.assemble(
                  - ( ( m_N * m_Em_der.tr() ) * meas(m_ori) ).tr()
                );

    gsMatrix<T> Fint = assembler.rhs();
    gsMatrix<T> result;
    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>*>(&u.source());
    gsMatrix<index_t> boundary;

    typedef gsBoundaryConditions<T> bcList;

    gsBoundaryConditions<real_t>::bcContainer container;
    m_bcs.getConditionFromSide(ps,container);

    for ( typename bcList::const_iterator it =  container.begin();
          it != container.end() ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;

        if( (it->unkComponent()!=com) && (it->unkComponent()!=-1) ) continue;

        const int k = it->patch();
        const gsBasis<T> & basis = mbasis[k];

        // if (it->type()==condition_type::dirichlet)
        //     gsDebug<<"Dirichlet\n";
        // else if (it->type()==condition_type::neumann)
        //     gsDebug<<"Neumann\n";
        // else
        //     GISMO_ERROR("Type unknown");

        // Get dofs on this boundary
        boundary = basis.boundary(it->side());
        result.resize(boundary.size(),1);

        T offset = u.mapper().offset(k);
        T size = u.mapper().size(com);
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            index_t ii = offset + size*com + boundary.at(l);
            result.at(l) = Fint.at(ii);
        }
    }

    // gsInfo<<"Neumann forces\n";
    // for ( typename bcList::const_iterator it =  m_bcs.begin("Neumann");
    //       it != m_bcs.end("Neumann") ; ++it )
    // {
    //     if( it->unknown()!=u.id() ) continue;
    //     //
    //     for (index_t r = 0; r!=u.dim(); ++r)
    //     {
    //         const int k = it->patch();
    //         const gsBasis<T> & basis = mbasis[k];

    //         // Get dofs on this boundary
    //         boundary = basis.boundary(it->side());

    //         T offset = u.mapper().offset(0);
    //         T size = u.mapper().size(r);

    //         for (index_t l=0; l!= boundary.size(); ++l)
    //         {
    //             index_t ii = offset + size*r + boundary.at(l);
    //             gsInfo<<result.at(ii)<<"\n";
    //         }
    //     }
    // }
    return result;
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assembleVector(const gsMatrix<T> & solVector)
{
    constructSolution(solVector, m_defpatches);
    assembleVector(m_defpatches);
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assemble(const gsMultiPatch<T> & deformed,
                                                   bool Matrix)
{
    if (Matrix)
        assembleMatrix(deformed);

    assembleVector(deformed);
}
template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::assemble(const gsMatrix<T> & solVector,
                                                   bool Matrix)
{
    constructSolution(solVector, m_defpatches);
    assemble(m_defpatches,Matrix);
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::constructSolution(const gsMatrix<T> & solVector) const
{
    gsMultiPatch<T> mp = m_patches;
    gsMultiPatch<T> displacement = constructDisplacement(solVector);
    for ( size_t k =0; k!=displacement.nPatches(); ++k) // Deform the geometry
        mp.patch(k).coefs() += displacement.patch(k).coefs();;  // defG points to mp_def, therefore updated

    return mp;
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructSolution(solVector);
}

template <short_t d, class T, bool bending>
T gsThinShellAssembler<d, T, bending>::getArea(const gsMultiPatch<T> & geometry)
{
    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(geometry);           // this map is used for integrals

    // Initialize vector
    geometryMap defG   = m_assembler.exprData()->getMap2();

    gsExprEvaluator<> evaluator(m_assembler);
    T result = evaluator.integral(meas(defG));
    return result;
}


template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::constructMultiPatch(const gsMatrix<T> & solVector) const
{
    m_solvector = solVector;
    gsMultiPatch<T> result;

    // Solution vector and solution variable
    space m_space = m_assembler.trialSpace(0);
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart() = m_ddofs;

    solution m_solution = m_assembler.getSolution(m_space, m_solvector);

    gsMatrix<T> cc;
    for ( size_t k =0; k!=m_defpatches.nPatches(); ++k) // Deform the geometry
    {
        // extract deformed geometry
        m_solution.extract(cc, k);
        result.addPatch(m_basis.basis(k).makeGeometry( give(cc) ));  // defG points to mp_def, therefore updated
    }

    return result;
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssembler<d, T, bending>::constructDisplacement(const gsMatrix<T> & solVector) const
{
    return constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending>
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

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructDisplacement(solVector);
}

// template <short_t d, class T, bool bending>
// void gsThinShellAssembler<d, T, bending>::constructStresses(const gsMultiPatch<T> & deformed,
//                                                     gsPiecewiseFunction<T> & result,
//                                                     stress_type::type type) const
// {
//     deformed = constructDisplacement(solVector);
// }

template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::computePrincipalStretches(const gsMatrix<T> & u, const gsMultiPatch<T> & deformed, const T z)
{
    // gsDebug<<"Warning: Principle Stretch computation of gsThinShellAssembler is depreciated...\n";
    gsMatrix<T> result(3,u.cols());
    result.setZero();
    // this->_getOptions();

    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);
    // m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::Stretch> m_mm(m_materialMat,m_defpatches);
    variable mm0 = m_assembler.getCoeff(m_mm);

    gsExprEvaluator<> evaluator(m_assembler);

    for (index_t k = 0; k != u.cols(); ++k)
    {
        result.col(k) = evaluator.eval(mm0,u.col(k));
    }
    return result;
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::constructStress(const gsMultiPatch<T> & deformed,
                                                    gsPiecewiseFunction<T> & result,
                                                    stress_type::type type)
{
    result.clear();

    for (size_t p = 0; p < m_patches.nPatches(); ++p )
        result.addPiecePointer(new gsShellStressFunction<T>(m_patches,deformed,m_materialMat,p,type,m_assembler));

}

// template <short_t d, class T, bool bending>
// gsField<T> gsThinShellAssembler<d, T, bending>::constructStress(const gsMultiPatch<T> & deformed,
//                                                     stress_type::type type)
// {
//     gsPiecewiseFunction<T> result;
//     result.clear();

//     for (size_t p = 0; p < m_patches.nPatches(); ++p )
//         result.addPiecePointer(new gsShellStressFunction<d, T, bending>(m_patches,deformed,m_materialMat,p,type,m_assembler));

//     gsField<T> stressField(m_patches,result, true);
//     return stressField;

// }

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::projectL2_into(const gsFunction<T> & fun, gsMatrix<T>& result)
{
    // this->_getOptions();

    m_assembler.cleanUp();
    m_assembler.setOptions(m_options);

    m_assembler.getMap(m_patches);           // this map is used for integrals

    // Initialize stystem
    m_assembler.initSystem(false);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    variable    function = m_assembler.getCoeff(fun, m_ori);
    // variable    function = m_assembler.getCoeff(fun);

    // assemble system
    m_assembler.assemble(m_space*m_space.tr()*meas(m_ori),m_space * function*meas(m_ori));
    m_solver.compute(m_assembler.matrix());
    result = m_solver.solve(m_assembler.rhs());
}

template <short_t d, class T, bool bending>
void gsThinShellAssembler<d, T, bending>::projectL2_into(const gsFunction<T> & fun, gsMultiPatch<T>& mp)
{
    gsMatrix<T> tmp = projectL2(fun);
    mp = m_patches;

    // Solution vector and solution variable
    space m_space = m_assembler.trialSpace(0);
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart() = m_ddofs;

    solution m_solution = m_assembler.getSolution(m_space, tmp);

    GISMO_ASSERT(m_defpatches.nPatches()==mp.nPatches(),"The number of patches of the result multipatch is not equal to that of the geometry!");

    gsMatrix<T> cc;
    for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
    {
        // // extract deformed geometry
        m_solution.extract(cc, k);
        mp.patch(k).coefs() = cc;  // defG points to mp_def, therefore updated
    }
}


template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssembler<d, T, bending>::projectL2(const gsFunction<T> & fun)
{
    gsMatrix<T> result;
    this->projectL2_into(fun,result);
    return result;
}

}// namespace gismo
