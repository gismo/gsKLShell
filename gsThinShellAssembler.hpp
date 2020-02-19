/** @file gsThinShellAssembler.hpp

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

#include <gsThinShell2/gsThinShellAssembler.h>
#include <gsThinShell2/gsMaterialMatrix.h>

namespace gismo
{

// template<class T>
// gsThinShellAssembler<T>::gsThinShellAssembler(  const gsMultiPatch<T> & patches,
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
gsThinShellAssembler<T>::gsThinShellAssembler(  const gsMultiPatch<T> & patches,
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
void gsThinShellAssembler<T>::defaultOptions()
{
    m_options.addInt("NonlinearLoads","Nonlinear Loads: 0: off, 1: on",nl_loads::off);
}

template <class T>
void gsThinShellAssembler<T>::getOptions() const
{
    m_nl_loads = m_options.getInt("NonlinearLoads");
}

/*
    TO INITIALIZE
        basis
        boundary conditions
*/

template <class T>
void gsThinShellAssembler<T>::initialize()
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
    m_space.setInterfaceCont(0); // todo: 1 (smooth basis)

    // Define fields as variables:
    // ... surface force
    //m_force = &
    // m_assembler.getCoeff(*m_forceFun,m_ori);
    // // ... thickness
    // m_assembler.getCoeff(*m_thickFun,m_ori);

    // call a defineComponents() depending on model

    assembleDirichlet();

    m_mapper = m_space.mapper();
    m_dim = m_space.dim();

    // foundation is off by default
    m_foundInd = false;

}


template <class T>
void gsThinShellAssembler<T>::assembleNeumann()
{

}

template <class T>
void gsThinShellAssembler<T>::assembleDirichlet()
{
    space m_space = m_assembler.trialSpace(0); // last argument is the space ID
    // if statement
    m_space.addBc( m_bcs.get("Dirichlet") ); // (!) must be called only once
}

template <class T>
void gsThinShellAssembler<T>::assembleClamped()
{

}

template<class T>
void gsThinShellAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<unsigned> acts,globalActs;

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

// template<class T>
// void gsShellAssembler<T>::assembleClamped()
// {
//     typedef gsBoundaryConditions<T>::const_iterator cIterator;

//     for( cIterator iit = m_bcs.begin("clamped"); iit!=m_bcs.end("clamped"); ++iit)
//     {
//         const boundary_condition<T> * it = &iit->get();

//         gsDofMapper & mapper  = m_dofMappers[it->unknown()];
//         const patchSide & cur = it->side();
//         // Get boundary dofs
//         gsMatrix<unsigned> bDofs = m_basis[0][cur.patch].boundary(cur);

//         // Cast to tensor b-spline basis
//         const gsTensorBSplineBasis<2,T> * tp =
//             dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&m_basis[0][cur.patch]);

//         if ( tp != NULL) // clamp adjacent dofs
//         {
//             const int str = tp->stride( cur.direction() );
//             if ( cur.parameter() )
//             {
//                 for ( index_t k=0; k<bDofs.size(); ++k)
//                     mapper.matchDof( cur.patch, (bDofs)(k,0),
//                                      cur.patch, (bDofs)(k,0) - str );
//             }
//             else
//             {
//                 for ( index_t k=0; k<bDofs.size(); ++k)
//                     mapper.matchDof( cur.patch, (bDofs)(k,0),
//                                      cur.patch, (bDofs)(k,0) + str );
//             }
//         }
//         else
//             gsWarn<<"Unable to apply clamped condition.\n";

//     }

//     // SET MAPPERS IN SPACE!
// }


// template <class T>
// void gsThinShellAssembler<T>::defineComponents()
// {

//     gsMaterialMatrix mm = m_materialMat;
//     // gsMaterialMatrix materialMat(m_patches, *m_thickFun, *m_YoungsModulus, *m_PoissonsRatio);
//     variable m_materialMat = m_assembler.getCoeff(materialMat);

//     gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
//     variable m_m2 = m_assembler.getCoeff(mult2t);

//     /*
//         We provide the following functions:                                 checked with previous assembler
//             m_Em            membrane strain tensor.                             V
//             m_Em_der        first variation of E_m.                             V
//             m_Em_der2       second variation of E_m MULTIPLIED BY S_m.          V
//             m_Ef            flexural strain tensor.                             V
//             m_Ef_der        second variation of E_f.                            V
//             m_Ef_der2       second variation of E_f MULTIPLIED BY S_f.          V

//         Where:
//             m_ori           geometry map of the the undeformed geometry,
//             m_def           geometry map of the deformed geometry,
//             m_materialMat   the material matrix,
//             m_m2            an auxillary matrix to multiply the last row of a tensor with 2
//             m_space         solution space
//     **/

//     space m_space = m_assembler.trialSpace(0);
//     geometryMap m_ori = m_assembler.exprData()->getMap();
//     geometryMap m_def = m_assembler.exprData()->getMap2();
//     variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);
//     variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

//     auto m_Em       = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
//     auto m_Sm      = m_Em * reshape(m_materialMat,3,3);
//     auto m_N        = m_thick.val() * m_Sm;
//     auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
//     auto m_Sm_der  = m_Em_der * reshape(m_materialMat,3,3);
//     auto m_N_der    = m_thick.val() * m_Sm_der;
//     auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

//     auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
//     auto m_Sf      = m_Ef * reshape(m_materialMat,3,3);
//     auto m_M        = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf;
//     auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
//     auto m_Sf_der  = m_Ef_der * reshape(m_materialMat,3,3);
//     auto m_M_der    = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf_der;
//     auto m_Ef_der2  = flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
//                         + var2(m_space,m_space,m_def,m_Ef * reshape(m_materialMat,3,3) );

//     //auto m_ff       = m_force;

// }

template<class T>
void gsThinShellAssembler<T>::assembleMass()
{
    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals

    // Initialize stystem
    m_assembler.initSystem();
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
void gsThinShellAssembler<T>::assembleFoundation()
{
    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    geometryMap m_ori   = m_assembler.exprData()->getMap();

    // Initialize stystem
    m_assembler.initSystem();
    variable    m_found = m_assembler.getCoeff(*m_foundFun, m_ori);
    space       m_space = m_assembler.trialSpace(0);

    m_assembler.assemble(m_found.val()*m_space*m_space.tr()*meas(m_ori));
}

template<class T>
void gsThinShellAssembler<T>::assemble()
{
    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize stystem
    m_assembler.initSystem();

    gsMaterialMatrix m_mm0 = m_materialMat;
    m_mm0.makeMatrix(0);
    gsMaterialMatrix m_mm1 = m_materialMat;
    m_mm1.makeMatrix(1);
    gsMaterialMatrix m_mm2 = m_materialMat;
    m_mm2.makeMatrix(2);

    // gsMaterialMatrix materialMat(m_patches, *m_thickFun, *m_YoungsModulus, *m_PoissonsRatio);
    variable mm0 = m_assembler.getCoeff(m_mm0);
    variable mm1 = m_assembler.getCoeff(m_mm1);
    variable mm2 = m_assembler.getCoeff(m_mm2);

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

    auto m_N_der    = m_Em_der * reshape(mm0,3,3) + m_Ef_der * reshape(mm1,3,3);
    auto m_M_der    = m_Em_der * reshape(mm1,3,3) + m_Ef_der * reshape(mm2,3,3);

    if (!m_foundInd) // no foundation
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
    // pt.setConstant(0.5);
    // gsExprEvaluator<> evaluator(m_assembler);
    // gsDebug<<evaluator.eval(reshape(mm1,3,3),pt)<<"\n";

}

// TO DO
// template <class T>
// bool gsThinShellAssembler<T>::assemble(const gsMatrix<T> & solutionVector,
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
void gsThinShellAssembler<T>::assembleMatrix(const gsMultiPatch<T> & deformed)
{
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize matrix
    // m_assembler.initMatrix();
    m_assembler.initSystem();

    gsMaterialMatrix m_mm0 = m_materialMat;
    m_mm0.makeMatrix(0);
    gsMaterialMatrix m_mm1 = m_materialMat;
    m_mm1.makeMatrix(1);
    gsMaterialMatrix m_mm2 = m_materialMat;
    m_mm2.makeMatrix(2);
    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector(0);
    gsMaterialMatrix m_S1 = m_materialMat;
    m_S1.makeVector(1);

    variable mm0 = m_assembler.getCoeff(m_mm0);
    variable mm1 = m_assembler.getCoeff(m_mm1);
    variable mm2 = m_assembler.getCoeff(m_mm2);
    variable S0 = m_assembler.getCoeff(m_S0);
    variable S1 = m_assembler.getCoeff(m_S1);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]


    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Ef_der2  = -(flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
                            + var2(m_space,m_space,m_def, m_M ));

    auto m_N_der    = m_Em_der * reshape(mm0,3,3) + m_Ef_der * reshape(mm1,3,3);
    auto m_M_der    = m_Em_der * reshape(mm1,3,3) + m_Ef_der * reshape(mm2,3,3);

    if (!m_foundInd) // no foundation
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
    // gsDebug<<"\n"<<evaluator.eval(reshape(mm0,3,3),pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(S0,pt)<<"\n";
    // gsDebug<<"\n"<<evaluator.eval(S2,pt)<<"\n";

}
template<class T>
void gsThinShellAssembler<T>::assembleMatrix(const gsMatrix<T> & solVector)
{
    constructSolution(solVector, m_defpatches);
    assembleMatrix(m_defpatches);
}

template <class T>
void gsThinShellAssembler<T>::assembleVector(const gsMultiPatch<T> & deformed)
{
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize vector
    m_assembler.initVector();

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

    auto m_N        = S0.tr();
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    auto m_M        = S1.tr(); // output is a column
    auto m_Ef_der   = -( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    if (!m_foundInd) // no foundation
    {
        // Assemble vector
        m_assembler.assemble(m_space * m_force * meas(m_ori) -
                    ( ( m_N * m_Em_der.tr() + m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                    );
    }
    // else
    // {
    //     variable    m_found = m_assembler.getCoeff(*m_foundFun, m_ori);
    //     // Assemble vector
    //     m_assembler.assemble(
    //                 m_space * m_force * meas(m_ori)
    //                 -
    //                 ( ( m_N * m_Em_der.tr() - m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
    //                 +

    //                 );
    // }


    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        applyLoads();
    }
}
template<class T>
void gsThinShellAssembler<T>::assembleVector(const gsMatrix<T> & solVector)
{
    constructSolution(solVector, m_defpatches);
    assembleVector(m_defpatches);
}

template<class T>
void gsThinShellAssembler<T>::assemble(const gsMultiPatch<T> & deformed,
                                        bool Matrix)
{
    m_defpatches = deformed;

    if (Matrix)
    {
        m_assembler.cleanUp();
        assembleMatrix(deformed);
    }

    m_assembler.cleanUp();
    assembleVector(deformed);
}
template<class T>
void gsThinShellAssembler<T>::assemble(const gsMatrix<T> & solVector,
                                        bool Matrix)
{
    constructSolution(solVector, m_defpatches);
    gsDebugVar(m_defpatches);
    assemble(m_defpatches,Matrix);
}

template <class T>
gsMultiPatch<T> gsThinShellAssembler<T>::constructSolution(const gsMatrix<T> & solVector) const
{
    m_solvector = solVector;
    gsMultiPatch<T> mp = m_patches;

    // Solution vector and solution variable
    space m_space = m_assembler.trialSpace(0);
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
void gsThinShellAssembler<T>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructSolution(solVector);
}

template <class T>
gsMultiPatch<T> gsThinShellAssembler<T>::constructDisplacement(const gsMatrix<T> & solVector) const
{
    gsMultiPatch<T> displacement = constructSolution(solVector);
    for ( size_t k =0; k!=displacement.nPatches(); ++k) // Deform the geometry
    {
        displacement.patch(k).coefs() -= m_patches.patch(k).coefs();;  // defG points to mp_def, therefore updated
    }

    return displacement;
}

template <class T>
void gsThinShellAssembler<T>::constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructDisplacement(solVector);
}

template <class T>
gsMatrix<T> gsThinShellAssembler<T>::computePrincipalStretches(const gsMatrix<T> & u, const gsMultiPatch<T> & deformed, const T z)
{
    gsMatrix<T> result(3,u.cols());

    this->getOptions();

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(deformed);

    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();

    m_assembler.initSystem();

    auto expr       = jac(m_def);
    auto expr2      = sn(m_def);

    gsExprEvaluator<> evaluator(m_assembler);

    gsMatrix<> metricA(3,3);
    gsMatrix<> evs;
    metricA.setZero();
    for (index_t k=0; k!=u.cols(); k++)
    {
        metricA.block(0,0,2,2) = evaluator.eval(expr,u.col(k)).block(0,0,2,2);
        metricA.col(2) = evaluator.eval(expr2,u.col(k));
        metricA = metricA.transpose()*metricA;
        Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;
        eigSolver.compute(metricA);
        evs = eigSolver.eigenvalues();
        for (index_t r=0; r!=3; r++)
            result(r,k) = math::sqrt(evs(r,0));
    }

    // result.col(0) = evaluator.eval(expr,pt);
    return result;
}

template <class T>
void gsThinShellAssembler<T>::constructStress(const gsMultiPatch<T> & deformed,
                                                    gsPiecewiseFunction<T> & result,
                                                    stress_type::type type)
{
    result.clear();

    for (index_t p = 0; p < m_patches.nPatches(); ++p )
        result.addPiecePointer(new gsShellStressFunction<T>(m_patches,deformed,m_materialMat,p,type,m_assembler));
}

/*
    To do; make warnings
*/

}// namespace gismo