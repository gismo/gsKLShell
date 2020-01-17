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
gsOptionList gsThinShellAssembler<T>::defaultOptions()
{
    gsOptionList & opt = m_assembler.options();
    // to do
    // opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::saint_venant_kirchhoff);
    return opt;
}

/*
    TO INITIALIZE
        basis
        boundary conditions
*/

template <class T>
void gsThinShellAssembler<T>::initialize()
{
    //gsInfo<<"Active options:\n"<< m_assembler.options() <<"\n";
    m_defpatches = m_patches;

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);

    // // initialize expression evaluator
    // gsExprEvaluator<T> evaluator(m_assembler);
    // m_evaluator = evaluator;

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

// template <class T>
// void gsThinShellAssembler<T>::createBendingStrips(gsMultiPatch<T> & mp, gsMultiPatch<T> & strips, gsStripIndices & indices)
// {
//     index_t bsize1, bsize2, par1, par2, dir1, dir2;
//     gsMatrix<unsigned> tempIndices;
//     gsMatrix<T> coefs;
//     std::vector<index_t> patches;
//     gsKnotVector<> kv0, kv1;
//     gsTensorBSplineBasis<2,real_t> basis;
//     gsTensorBSpline<2,real_t> shape;

//     indices.clear();

//     GISMO_ASSERT(( mp.nPatches()==1 || mp.nInterfaces()!=0 ),"Multiple patches, but no interfaces?? nPatches="<<mp.nPatches()<<", nInterfaces="<<mp.nInterfaces());
//     if (mp.nInterfaces()!=0)
//     {
//         for (typename gsMultiPatch<>::const_iiterator it = mp.iBegin(); it<mp.iEnd(); it++)
//         {
//             bsize1 = mp.basis(it->first().patch).boundary(it->first()).size();
//             bsize2 = mp.basis(it->second().patch).boundary(it->second()).size();

//             GISMO_ASSERT(bsize2==bsize1, "Number of control points along boundary are not equal.");

//             // Collect indices of patch
//             tempIndices = gsMatrix<unsigned>(bsize1,3);


//             tempIndices.col(0) = mp.basis(it->first().patch).boundaryOffset(it->first(),1);
//             tempIndices.col(1) = mp.basis(it->second().patch).boundaryOffset(it->second(),0);
//             tempIndices.col(2) = mp.basis(it->second().patch).boundaryOffset(it->second(),1);

//             par1 = it->first().parameter();
//             par2 = it->second().parameter();
//             dir1 = it->first().direction();
//             dir2 = it->second().direction();

//             // matchInterface and then shift with stride index. NOTE: this does not work since the stride is not implemented for gsBasis and some TensorBases...
//             // gsMatrix<unsigned> bndThis, bndOther;
//             // mp.basis(it->first().patch).matchWith(*it, mp.basis(it->second().patch), bndThis,bndOther);
//             // gsDebugVar(bndThis);
//             // gsDebugVar(bndOther);

//             if (((par1 != par2) && ( dir1 != dir2 ))) // || ((par1 == par2) && ( dir1 == dir2 )) //then the directions should be different
//                 tempIndices.col(0).reverseInPlace();

//             // Collect control point coordinates
//             coefs = gsMatrix<T>(3*bsize1,mp.targetDim());

//             // Helper vector for patch numbers
//             patches = std::vector<index_t>{it->first().patch, it->second().patch, it->second().patch};
//             for (index_t j=0; j!=3; ++j)
//             for (index_t i=0; i!=bsize1; ++i)
//                       coefs.row(i*3 + j) = mp.patch(patches[j]).coef(tempIndices(i,j));

//             // Make a basis
//             kv0.initUniform(0,1,0,3,1);
//             kv1.initUniform(0,1,bsize1-2,2,1);  // TODO: WHAT IF BSIZE1<2?

//             basis = gsTensorBSplineBasis<2,real_t>(kv0,kv1);
//             shape = gsTensorBSpline<2,real_t>(basis,coefs);

//             strips.addPatch(shape);
//             indices.push_back(tempIndices);
//         }
//     }
// }

// template <class T>
// void gsThinShellAssembler<T>::assembleBendingStrips()
// {
//     // make strips on m_patches and on m_defpatches
//     this->createBendingStrips(m_patches,m_strips,m_stripIndices);

//     gsMaterialMatrix m_mm0 = m_materialMat;
//     gsMaterialMatrix m_mm2 = m_materialMat;
//     m_mm2.setMoment(2);
//     variable mm0 = m_stripAssembler.getCoeff(m_mm0);
//     variable mm2 = m_stripAssembler.getCoeff(m_mm2);

//     gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
//     variable m_m2 = m_stripAssembler.getCoeff(mult2t);

//     for (index_t k=0; k!=m_strips.nPatches(); ++k) // for each strip
//     {
//          geometryMap S   = m_stripAssembler.getMap(m_strips.patch(k));

//         // set assembler for the strips (over all strips together)
//         m_stripAssembler.setIntegrationElements(m_strips.basis(k));
//         space stripSpace = m_stripAssembler.getSpace(m_strips.basis(k), 3, 0); // last argument is the space ID

//         auto m_Ef_der   = ( deriv2(stripSpace,sn(S).normalized().tr() ) + deriv2(S,var1(stripSpace,S) ) ) * reshape(m_m2,3,3); //[checked]
//         auto m_M_der    = m_Ef_der * reshape(mm2,3,3);

//         m_stripAssembler.initSystem();
//         m_stripAssembler.assemble( ( m_M_der * m_Ef_der.tr() ) * meas(S) );

//         // for (index_t i=0; i!=m_stripIndices[k].rows(); ++i)
//         // {
//         //     index_t gidx = stripSpace.index(  )
//         //     (this->matrix())(stripSpace.index)
//         // }
//         // make assembler
//         // assemble object
//         // add the entries to the systrem matrix BUT THAT IS A REFERENCE TO AN OBJECT IN gsExprAssembler!!
//     }


//     // after writing in sparsematrix;
//     // m_matrix.makeCompressed();
// }

// template <class T>
// void gsThinShellAssembler<T>::assembleBendingStrips()
// {
//     // make strips on m_patches and on m_defpatches
//     this->createBendingStrips(m_patches,m_strips,m_stripIndices);
//     this->createBendingStrips(m_defpatches,m_defstrips,m_defstripIndices);


//     geometryMap S;
//     geometryMap defS;

//     gsMaterialMatrix m_mm0 = m_materialMat;
//     gsMaterialMatrix m_mm2 = m_materialMat;
//     m_mm2.setMoment(2);
//     variable mm0 = m_stripAssembler.getCoeff(m_mm0);
//     variable mm2 = m_stripAssembler.getCoeff(m_mm2);

//     gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
//     variable m_m2 = m_stripAssembler.getCoeff(mult2t);

//     auto m_Ef       = ( deriv2(m_ori,sn(S).normalized().tr()) - deriv2(m_def,sn(defS).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
//     auto m_M        = m_Ef * reshape(mm,3,3);

//     auto m_Ef_der   = ( deriv2(m_space,sn(defS).normalized().tr() ) + deriv2(m_def,var1(m_space,defS) ) ) * reshape(m_m2,3,3); //[checked]
//     auto m_M_der    = m_Ef_der * reshape(mm2,3,3);

//     for (index_t k=0; k!=m_strips.nPatches(); ++k)
//     {
//         S   = A.getMap(m_strips.patch(k));
//         defS= A.getMap(m_defstrips.patch(k));

//         // set assembler for the strips (over all strips together)
//         m_stripAssembler.setIntegrationElements(m_strips.basis(k));
//         space stripSpace = m_stripAssembler.getSpace(m_strips.basis(k), 3, 0); // last argument is the space ID


//         m_stripAssembler.initSystem();

//         m_stripAssembler.assemble(
//             ( m_M_der * m_Ef_der.tr() ) * meas(m_ori)
//             ,
//             - ( - m_M * m_Ef_der.tr() * meas(m_ori) ).tr()
//             );


//         // make assembler
//         // assemble object
//         // add the entries to the systrem matrix BUT THAT IS A REFERENCE TO AN OBJECT IN gsExprAssembler!!
//     }


//     // after writing in sparsematrix;
//     // m_matrix.makeCompressed();
// }

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
void gsThinShellAssembler<T>::assemble()
{
    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize stystem
    m_assembler.initSystem();

    gsMaterialMatrix m_mm0 = m_materialMat;
    gsMaterialMatrix m_mm2 = m_materialMat;
    m_mm2.setMoment(2);
    // gsMaterialMatrix materialMat(m_patches, *m_thickFun, *m_YoungsModulus, *m_PoissonsRatio);
    variable mm0 = m_assembler.getCoeff(m_mm0);
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
    auto m_N_der    = m_Em_der * reshape(mm0,3,3);

    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    // auto m_Sf_der   = m_Ef_der * reshape(m_materialMat,3,3);
    // auto m_M_der    = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf_der;
    auto m_M_der    = m_Ef_der * reshape(mm2,3,3);

    // auto m_M_der    = pow(m_thick.val(),3) / 12.0 * m_Sf_der;

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

    // Assemble the loads
    if ( m_pLoads.numLoads() != 0 )
    {
        m_rhs = m_assembler.rhs();
        applyLoads();
    }
    // Neumann
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
    gsMaterialMatrix m_mm2 = m_materialMat;
    m_mm2.setMoment(2);
    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector();
    gsMaterialMatrix m_S2 = m_materialMat;
    m_S2.makeVector();
    m_S2.setMoment(2);

    variable mm0 = m_assembler.getCoeff(m_mm0);
    variable mm2 = m_assembler.getCoeff(m_mm2);
    variable S0 = m_assembler.getCoeff(m_S0);
    variable S2 = m_assembler.getCoeff(m_S2);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    auto m_Em       = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
    // auto m_Sm       = m_Em * reshape(mm0,3,3);
    // auto m_N2        = m_thick.val() * m_Sm;
    auto m_N        = S0;
    // auto m_N        = m_Em * reshape(mm0,3,3);

    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    // auto m_Sm_der   = m_Em_der * reshape(m_materialMat,3,3);
    // auto m_N_der    = m_thick.val() * m_Sm_der;
    auto m_N_der    = m_Em_der * reshape(mm0,3,3);

    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N.tr() ); //[checked]


    auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
    // auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) ) ;//* reshape(m_m2,3,3) ; //[checked]
    // auto m_Sf       = m_Ef * reshape(m_materialMat,3,3);
    // auto m_M        = pow(m_thick.val(),3) / 12.0 * m_Sf;
    auto m_M        = S2; // output is a column
    // auto m_M        = m_Ef * reshape(mm2,3,3);

    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    // auto m_Sf_der   = m_Ef_der * reshape(m_materialMat,3,3);
    // auto m_M_der    = pow(m_thick.val(),3) / 12.0 * m_Sf_der;
    auto m_M_der    = m_Ef_der * reshape(mm2,3,3);

    auto m_Ef_der2  = flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M.tr()  ).symmetrize()
                        + var2(m_space,m_space,m_def, m_M.tr() );

    // Assemble matrix
    m_assembler.assemble(
                (
                    m_N_der * m_Em_der.tr()
                    +
                    m_Em_der2
                    +
                    m_M_der * m_Ef_der.tr()
                    -
                    m_Ef_der2
                    ) * meas(m_ori)
                );
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

    gsMaterialMatrix m_mm0 = m_materialMat;
    gsMaterialMatrix m_mm2 = m_materialMat;
    m_mm2.setMoment(2);
    // gsMaterialMatrix materialMat(m_patches, *m_thickFun, *m_YoungsModulus, *m_PoissonsRatio);
    variable mm0 = m_assembler.getCoeff(m_mm0);
    variable mm2 = m_assembler.getCoeff(m_mm2);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space m_space       = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun, m_ori);
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    auto m_Em       = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
    // auto m_Sm       = m_Em * reshape(m_materialMat,3,3);
    // auto m_N        = m_thick.val() * m_Sm;
    auto m_N    = m_Em * reshape(mm0,3,3);

    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ;

    auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
    // auto m_Sf       = m_Ef * reshape(m_materialMat,3,3);
    // auto m_M        = pow(m_thick.val(),3) / 12.0 * m_Sf;
    auto m_M        = m_Ef * reshape(mm2,3,3);

    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    // Assemble vector
    m_assembler.assemble(m_space * m_force * meas(m_ori) -
                ( ( m_N * m_Em_der.tr() - m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                );

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

    // m_assembler.cleanUp();

    // assemble();
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

/*
    To do; make warnings
*/
template <class T>
void gsThinShellAssembler<T>::constructStress(const gsMultiPatch<T> & deformed,
                                                    gsPiecewiseFunction<T> & result,
                                                    stress_type::type type)
{
    result.clear();

    for (index_t p = 0; p < m_patches.nPatches(); ++p )
        result.addPiecePointer(new gsShellStressFunction<T>(m_patches,deformed,m_materialMat,p,type,m_assembler));
}

}// namespace gismo