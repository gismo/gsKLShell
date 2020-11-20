/** @file gsMembraneAssembler.hpp

    @brief Provides linear and nonlinear elasticity systems for thin shells.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsThinPlateAssembler.h>
#include <gsKLShell/gsMaterialMatrix.h>

namespace gismo
{

// template<class T>
// gsMembraneAssembler<T>::gsMembraneAssembler(  const gsMultiPatch<T> & patches,
//                                         const gsMultiBasis<T> & basis,
//                                         const gsBoundaryConditions<T> & bconditions,
//                                         const gsFunction<T> & surface_force,
//                                         const gsFunction<T> & thickness,
//                                         T YoungsModulus,
//                                         T PoissonsRatio
//                                         )
//                                         :
//                                         m_patches(patches),
//                                         m_basis(basis),
//                                         m_bcs(bconditions),
//                                         m_forceFun(&surface_force),
//                                         m_thickFun(&thickness)
// {
//     m_YoungsModulus = memory::make_shared(new gsConstantFunction(YoungsModulus,3));
//     m_PoissonsRatio = memory::make_shared(new gsConstantFunction(PoissonsRatio,3));
//     this->initialize();
// }

template<class T>
gsMembraneAssembler<T>::gsMembraneAssembler(  const gsMultiPatch<T> & patches,
                                        const gsMultiBasis<T> & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunction<T> & surface_force,
                                        const gsMaterialMatrix<T> & materialmatrix
                                        )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_bcs(bconditions),
                                        m_forceFun(&surface_force),
                                        m_materialMat(materialmatrix)
{
    this->initialize();
}

template <class T>
void gsMembraneAssembler<T>::defaultOptions()
{
    m_options.addInt("NonlinearLoads","Nonlinear Loads: 0: off, 1: on",nl_loads::off);
}

template <class T>
void gsMembraneAssembler<T>::getOptions() const
{
    m_nl_loads = m_options.getInt("NonlinearLoads");
}

/*
    TO INITIALIZE
        basis
        boundary conditions
*/

template <class T>
void gsMembraneAssembler<T>::initialize()
{
    this->defaultOptions();

    //gsInfo<<"Active options:\n"<< m_assembler.options() <<"\n";
    m_defpatches = m_patches;

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);

    // Initialize the geometry maps
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Set the discretization space
    space m_space = m_assembler.getSpace(m_basis, 3, 0); // last argument is the space ID
    // m_space.setup(m_bcs, dirichlet::interpolation, 0);

    // Define fields as variables:
    // ... surface force
    //m_force = &
    // m_assembler.getCoeff(*m_forceFun,m_ori);
    // // ... thickness
    // m_assembler.getCoeff(*m_thickFun,m_ori);

    // call a defineComponents() depending on model

    this->assembleDirichlet();

    m_ddofs = m_space.fixedPart();

    m_mapper = m_space.mapper();
    m_dim = m_space.dim();

    // foundation is off by default
    m_foundInd = false;
    // pressure is off by default
    m_pressInd = false;
}


template <class T>
void gsMembraneAssembler<T>::assembleNeumann()
{
    m_assembler.getMap(m_patches);           // this map is used for integrals
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    variable g_N = m_assembler.getBdrFunction();
    m_assembler.assembleRhsBc(m_space * g_N * otangent(m_ori).norm(), m_bcs.neumannSides() );
}

template <class T>
void gsMembraneAssembler<T>::assembleDirichlet()
{
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // if statement
    m_space.setup(m_bcs, dirichlet::interpolation, 0);
}

template <class T>
void gsMembraneAssembler<T>::homogenizeDirichlet()
{
    // space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // m_space.setup(m_bcs, dirichlet::homogeneous, 0);
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart().setZero();
}

template<class T>
void gsMembraneAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts,globalActs;

    space       m_space = m_assembler.trialSpace(0);
    m_mapper = m_space.mapper();

    // /*
    //     NOTE!!
    //     This does not work yet. We need the dofMappers per degree of freedom!!
    // */
    // // -----------------------------------
    // m_dofMappers.resize(3);
    // m_dofMappers.at(0) = m_space.mapper();
    // m_dofMappers.at(1) = m_space.mapper();
    // m_dofMappers.at(2) = m_space.mapper();
    // // -----------------------------------


    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
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
        for (size_t j = 0; j< 3; ++j)
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

template<class T>
void gsMembraneAssembler<T>::assembleMass()
{
    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals

    // Initialize stystem
    m_assembler.initSystem(false);
    gsMaterialMatrix m_mm = m_materialMat;
    m_mm.makeDensity();
    variable mm0 = m_assembler.getCoeff(m_mm);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    // gsVector<> pt(2);
    // pt.setConstant(0.5);
    // gsExprEvaluator<T> evaluator(m_assembler);
    // gsDebug<<evaluator.eval(mm0,pt)<<"\n";

    // assemble system
    m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
}

template<class T>
void gsMembraneAssembler<T>::assembleFoundation()
{
    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    // Initialize stystem
    m_assembler.initSystem(false);
    variable    m_found = m_assembler.getCoeff(*m_foundFun, m_ori);
    space       m_space = m_assembler.trialSpace(0);

    m_assembler.assemble(m_found.val()*m_space*m_space.tr()*meas(m_ori));
}

template<class T>
void gsMembraneAssembler<T>::assemble()
{
    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize stystem
    m_assembler.initSystem(false);

    gsMaterialMatrix m_mmA = m_materialMat;
    m_mmA.makeMatrix(0);
    gsMaterialMatrix m_mmB = m_materialMat;
    m_mmB.makeMatrix(1);
    gsMaterialMatrix m_mmC = m_materialMat;
    m_mmC.makeMatrix(2);
    gsMaterialMatrix m_mmD = m_materialMat;
    m_mmD.makeMatrix(3);
    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector(0);
    gsMaterialMatrix m_S1 = m_materialMat;
    m_S1.makeVector(1);

    variable mmA = m_assembler.getCoeff(m_mmA);
    variable mmB = m_assembler.getCoeff(m_mmB);
    variable mmC = m_assembler.getCoeff(m_mmC);
    variable mmD = m_assembler.getCoeff(m_mmD);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    // auto m_Sm_der   = m_Em_der * reshape(m_materialMat,3,3);
    // auto m_N_der    = m_thick.val() * m_Sm_der;
    // auto m_N_der    = m_Em_der * reshape(mm0,3,3);

    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    // auto m_Sf_der   = m_Ef_der * reshape(m_materialMat,3,3);
    // auto m_M_der    = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf_der;
    // auto m_M_der    = m_Ef_der * reshape(mm2,3,3);

    // auto m_M_der    = pow(m_thick.val(),3) / 12.0 * m_Sf_der;

    auto m_N_der    = m_Em_der * reshape(mmA,3,3) + m_Ef_der * reshape(mmB,3,3);
    auto m_M_der    = m_Em_der * reshape(mmC,3,3) + m_Ef_der * reshape(mmD,3,3);

    if (m_foundInd) // no foundation
    {
        GISMO_ERROR("Not implemented");
    }
    else if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);

        // gsVector<> pt(2);
        // pt.setConstant(0.25);
        // gsExprEvaluator<> evaluator(m_assembler);

        // gsDebug<<evaluator.eval(m_pressure.val(),pt)<<"\n";
        // gsDebug<<evaluator.eval(m_space,pt)<<"\n";
        // gsDebug<<evaluator.eval(sn(m_def).normalized(),pt)<<"\n";

        // Assemble vector
        m_assembler.assemble(
            (
                m_N_der * m_Em_der.tr()
                +
                m_M_der * m_Ef_der.tr()
            ) * meas(m_ori)
            ,
            m_space * m_force * meas(m_ori) + m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
            );

        this->assembleNeumann();
    }

    else // no foundation, no pressure
    {
        // assemble system
        m_assembler.assemble(
            (
                m_N_der * m_Em_der.tr()
                +
                m_M_der * m_Ef_der.tr()
            ) * meas(m_ori)
            ,
            m_space * m_force * meas(m_ori)
            );

        this->assembleNeumann();
    }
    // else
    // {
    //     variable    m_found = m_assembler.getCoeff(*m_foundFun, m_ori);
    //     // assemble system
    //     m_assembler.assemble(
    //         (
    //             m_N_der * m_Em_der.tr()
    //             +
    //             m_M_der * m_Ef_der.tr()
    //             +
    //             m_found.val()*m_space*m_space.tr()
    //         ) * meas(m_ori)
    //         ,
    //         m_space * m_force * meas(m_ori)
    //         );
    // }


    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        applyLoads();
    }
    // Neumann

    // gsVector<> pt(2);
    // pt.setConstant(0.25);
    // gsExprEvaluator<> evaluator(m_assembler);
    // gsDebug<<evaluator.eval(reshape(mmA,3,3),pt)<<"\n";
    // gsDebug<<evaluator.eval(reshape(mmB,3,3),pt)<<"\n";
    // gsDebug<<evaluator.eval(reshape(mmC,3,3),pt)<<"\n";
    // gsDebug<<evaluator.eval(reshape(mmD,3,3),pt)<<"\n";

    // gsDebug<<evaluator.eval((m_N_der * m_Em_der.tr()
    //             +
    //             m_M_der * m_Ef_der.tr()
    //         ) * meas(m_ori),pt)<<"\n";

}

// TO DO
// template <class T>
// bool gsMembraneAssembler<T>::assemble(const gsMatrix<T> & solutionVector,
//                                         const std::vector<gsMatrix<T> > & fixedDoFs,
//                                         bool assembleMatrix)
// {
//     gsMultiPatch<T> displacement;
//     constructSolution(solutionVector,fixedDoFs,displacement);
//     if (checkSolution(displacement) != -1)
//         return false;
//     assemble(displacement,assembleMatrix);
//     return true;
// }

template <class T>
void gsMembraneAssembler<T>::assembleMatrix(const gsMultiPatch<T> & deformed)
{
    m_assembler.cleanUp();
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize matrix
    m_assembler.initMatrix(false);
    // m_assembler.initSystem(false);

    gsMaterialMatrix m_mmA = m_materialMat;
    m_mmA.makeMatrix(0);
    gsMaterialMatrix m_mmB = m_materialMat;
    m_mmB.makeMatrix(1);
    gsMaterialMatrix m_mmC = m_materialMat;
    m_mmC.makeMatrix(2);
    gsMaterialMatrix m_mmD = m_materialMat;
    m_mmD.makeMatrix(3);
    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector(0);
    gsMaterialMatrix m_S1 = m_materialMat;
    m_S1.makeVector(1);

    variable mmA = m_assembler.getCoeff(m_mmA);
    variable mmB = m_assembler.getCoeff(m_mmB);
    variable mmC = m_assembler.getCoeff(m_mmC);
    variable mmD = m_assembler.getCoeff(m_mmD);
    variable S0 = m_assembler.getCoeff(m_S0);
    variable S1 = m_assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

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

    if (m_foundInd) // no foundation
    {
        GISMO_ERROR("not implemented");
    }
    else if (m_pressInd)
    {
        // gsVector<> pt(2);
        // pt.setConstant(0.25);
        // gsExprEvaluator<> evaluator(m_assembler);
        // gsDebug<<evaluator.eval(m_def,pt)<<"\n";

        // gsDebug<<evaluator.eval(m_pressure.val(),pt)<<"\n";
        // gsDebug<<evaluator.eval(m_space,pt)<<"\n";
        // gsDebug<<evaluator.eval(sn(m_def).normalized(),pt)<<"\n";


        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);
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
                -
                m_pressure.val() * m_space * var1(m_space,m_def).tr()* meas(m_ori)
            );
    }
    else // no foundation, no pressure
    {
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
    }
    // else
    // {
    //     variable    m_found = m_assembler.getCoeff(*m_foundFun, m_ori);
    //     // Assemble matrix
    //     m_assembler.assemble(
    //             (
    //                  m_N_der * m_Em_der.tr()
    //                  +
    //                  m_Em_der2
    //                  +
    //                  m_M_der * m_Ef_der.tr()
    //                  -
    //                  m_Ef_der2
    //                  +
    //                 m_found.val()*m_space*m_space.tr()
    //             ) * meas(m_ori)
    //         );
    // }
    // gsVector<> pt(2);
    // pt.setConstant(0.5);
    // gsExprEvaluator<> evaluator(m_assembler);
    // gsDebug<<"\n"<<evaluator.eval(reshape(mmA,3,3),pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(reshape(mmB,3,3),pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(reshape(mmC,3,3),pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(reshape(mmD,3,3),pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(S0,pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(S2,pt)<<"\n";

}
template<class T>
void gsMembraneAssembler<T>::assembleMatrix(const gsMatrix<T> & solVector)
{
    // gsMultiPatch<T> deformed;
    // constructSolution(solVector, deformed);
    // assembleMatrix(deformed);

    constructSolution(solVector, m_defpatches);
    assembleMatrix(m_defpatches);
}

template <class T>
void gsMembraneAssembler<T>::assembleVector(const gsMultiPatch<T> & deformed)
{
    m_assembler.cleanUp();
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize vector
    m_assembler.initVector(1,false);

    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector(0);
    gsMaterialMatrix m_S1 = m_materialMat;
    m_S1.makeVector(1);

    variable S0 = m_assembler.getCoeff(m_S0);
    variable S1 = m_assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space m_space       = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    this->homogenizeDirichlet();

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    if (m_foundInd) // no foundation
    {
        //     variable    m_found = m_assembler.getCoeff(*m_foundFun, m_ori);
        //     // Assemble vector
        //     m_assembler.assemble(
        //                 m_space * m_force * meas(m_ori)
        //                 -
        //                 ( ( m_N * m_Em_der.tr() - m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
        //                 +

        //                 );
        GISMO_ERROR("Not implemented");
    }
    else if (m_pressInd)
    {
        variable m_pressure = m_assembler.getCoeff(*m_pressFun, m_ori);

        // Assemble vector
        m_assembler.assemble(
                        m_space * m_force * meas(m_ori)
                      + m_pressure.val() * m_space * sn(m_def).normalized() * meas(m_ori)
                      - ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                    );

        this->assembleNeumann();
    }
    else
    {
        // Assemble vector
        m_assembler.assemble(m_space * m_force * meas(m_ori) -
                    ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                    );

        this->assembleNeumann();
    }
    // else
    // {

    // }


    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        applyLoads();
    }
    // gsVector<> pt(2);
    // pt.setConstant(0.5);
    // gsExprEvaluator<> evaluator(m_assembler);
    // gsDebug<<"\n"<<evaluator.eval(S0,pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(S1,pt)<<"\n";
}

template <class T>
gsMatrix<T> gsMembraneAssembler<T>::boundaryForceVector(const gsMultiPatch<T> & deformed, patchSide& ps, int com)
{
    gsExprAssembler<T> assembler;
    assembler.setIntegrationElements(m_basis);
    space u = assembler.getSpace(m_basis, 3, 0); // last argument is the space ID

    assembler.initSystem();

    m_defpatches = deformed;

    assembler.getMap(m_patches);           // this map is used for integrals
    assembler.getMap(m_defpatches);

    // Initialize vector
    // m_assembler.initVector(1,false);

    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector(0);
    gsMaterialMatrix m_S1 = m_materialMat;
    m_S1.makeVector(1);

    variable S0 = assembler.getCoeff(m_S0);
    variable S1 = assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = assembler.getCoeff(mult2t);

    geometryMap m_ori   = assembler.exprData()->getMap();
    geometryMap m_def   = assembler.exprData()->getMap2();
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

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
template<class T>
void gsMembraneAssembler<T>::assembleVector(const gsMatrix<T> & solVector)
{
    // gsMultiPatch<T> deformed;
    // constructSolution(solVector, deformed);
    // assembleVector(deformed);

    constructSolution(solVector, m_defpatches);
    assembleVector(m_defpatches);
}

template<class T>
void gsMembraneAssembler<T>::assemble(const gsMultiPatch<T> & deformed,
                                        bool Matrix)
{
    if (Matrix)
        assembleMatrix(deformed);

    assembleVector(deformed);
}
template<class T>
void gsMembraneAssembler<T>::assemble(const gsMatrix<T> & solVector,
                                        bool Matrix)
{
    // gsMultiPatch<T> deformed;
    // constructSolution(solVector, deformed);
    // assemble(deformed,Matrix);

    constructSolution(solVector, m_defpatches);
    assemble(m_defpatches,Matrix);
}

template <class T>
gsMultiPatch<T> gsMembraneAssembler<T>::constructSolution(const gsMatrix<T> & solVector) const
{
    m_solvector = solVector;
    gsMultiPatch<T> mp = m_patches;

    // Solution vector and solution variable
    space m_space = m_assembler.trialSpace(0);
    const_cast<expr::gsFeSpace<T> & >(m_space).fixedPart() = m_ddofs;

    solution m_solution = m_assembler.getSolution(m_space, m_solvector);

    GISMO_ASSERT(m_defpatches.nPatches()==mp.nPatches(),"The number of patches of the result multipatch is not equal to that of the geometry!");

    gsMatrix<T> cc;
    for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
    {
        // // extract deformed geometry
        m_solution.extract(cc, k);
        mp.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    return mp;
}

template <class T>
void gsMembraneAssembler<T>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructSolution(solVector);
}

template <class T>
T gsMembraneAssembler<T>::getArea(const gsMultiPatch<T> & geometry)
{
    // if (deformed)
    //     m_assembler.getMap(m_defpatches);           // this map is used for integrals
    // else
    //     m_assembler.getMap(m_patches);           // this map is used for integrals

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(geometry);           // this map is used for integrals

    // Initialize vector
    geometryMap defG   = m_assembler.exprData()->getMap2();

    gsExprEvaluator<> evaluator(m_assembler);
    T result = evaluator.integral(meas(defG));
    return result;
}


template <class T>
gsMultiPatch<T> gsMembraneAssembler<T>::constructDisplacement(const gsMatrix<T> & solVector) const
{
    gsMultiPatch<T> displacement = constructSolution(solVector);
    for ( size_t k =0; k!=displacement.nPatches(); ++k) // Deform the geometry
    {
        displacement.patch(k).coefs() -= m_patches.patch(k).coefs();;  // defG points to mp_def, therefore updated
    }

    return displacement;
}

template <class T>
void gsMembraneAssembler<T>::constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructDisplacement(solVector);
}

// template <class T>
// void gsMembraneAssembler<T>::constructStresses(const gsMultiPatch<T> & deformed,
//                                                     gsPiecewiseFunction<T> & result,
//                                                     stress_type::type type) const
// {
//     deformed = constructDisplacement(solVector);
// }

template <class T>
gsMatrix<T> gsMembraneAssembler<T>::computePrincipalStretches(const gsMatrix<T> & u, const gsMultiPatch<T> & deformed, const T z)
{
    // gsDebug<<"Warning: Principle Stretch computation of gsMembraneAssembler is depreciated...\n";
    gsMatrix<T> result(3,u.cols());
    result.setZero();
    // this->getOptions();

    m_assembler.cleanUp();
    m_defpatches = deformed;

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();

    // m_assembler.initSystem(false);

    gsMaterialMatrix m_mm = m_materialMat;
    m_mm.makeStretch();
    variable mm0 = m_assembler.getCoeff(m_mm);

    gsExprEvaluator<> evaluator(m_assembler);

    for (index_t k = 0; k != u.cols(); ++k)
    {
        result.col(k) = evaluator.eval(mm0,u.col(k));
    }
    return result;
}

template <class T>
void gsMembraneAssembler<T>::constructStress(const gsMultiPatch<T> & deformed,
                                                    gsPiecewiseFunction<T> & result,
                                                    stress_type::type type)
{
    result.clear();

    for (size_t p = 0; p < m_patches.nPatches(); ++p )
        result.addPiecePointer(new gsShellStressFunction<T>(m_patches,deformed,m_materialMat,p,type,m_assembler));

}

// template <class T>
// gsField<T> gsMembraneAssembler<T>::constructStress(const gsMultiPatch<T> & deformed,
//                                                     stress_type::type type)
// {
//     gsPiecewiseFunction<T> result;
//     result.clear();

//     for (size_t p = 0; p < m_patches.nPatches(); ++p )
//         result.addPiecePointer(new gsShellStressFunction<T>(m_patches,deformed,m_materialMat,p,type,m_assembler));

//     gsField<T> stressField(m_patches,result, true);
//     return stressField;

// }

}// namespace gismo