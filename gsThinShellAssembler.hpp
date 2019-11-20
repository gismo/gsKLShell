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

template<class T>
gsThinShellAssembler<T>::gsThinShellAssembler(  const gsMultiPatch<T> & patches,
                                        const gsMultiBasis<T> & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunction<T> & surface_force,
                                        const gsFunction<T> & thickness,
                                        T YoungsModulus,
                                        T PoissonsRatio
                                        )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_bcs(bconditions),
                                        m_forceFun(&surface_force),
                                        m_thickFun(&thickness)
{
    m_YoungsModulus = memory::make_shared(new gsConstantFunction(YoungsModulus,3));
    m_PoissonsRatio = memory::make_shared(new gsConstantFunction(PoissonsRatio,3));
    this->initialize();
}

template<class T>
gsThinShellAssembler<T>::gsThinShellAssembler(  const gsMultiPatch<T> & patches,
                                        const gsMultiBasis<T> & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunction<T> & surface_force,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & YoungsModulus,
                                        const gsFunction<T> & PoissonsRatio
                                        )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_bcs(bconditions),
                                        m_forceFun(&surface_force),
                                        m_thickFun(&thickness),
                                        m_YoungsModulus(memory::make_shared_not_owned(&YoungsModulus)),
                                        m_PoissonsRatio(memory::make_shared_not_owned(&PoissonsRatio))
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

    // Set the geometry map
    geometryMap m_ori = m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Set the discretization space
    space m_space = m_assembler.getSpace(m_basis, 3, 0); // last argument is the space ID
    m_space.setInterfaceCont(0); // todo: 1 (smooth basis)
    m_space.addBc( m_bcs.get("Dirichlet") ); // (!) must be called only once

    // Solution vector and solution variable
    m_solution = m_assembler.getSolution(m_space, m_solvector);

    // Define fields as variables:
    // ... surface force
    //m_force = &
    m_assembler.getCoeff(*m_forceFun,m_ori);
    // ... thickness
    m_assembler.getCoeff(*m_thickFun,m_ori);

    // call a defineComponents() depending on model

}

template <class T>
void gsThinShellAssembler<T>::defineComponents()
{
    gsMaterialMatrix materialMat(m_patches, *m_YoungsModulus, *m_PoissonsRatio);
    variable m_materialMat = m_assembler.getCoeff(materialMat);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    /*
        We provide the following functions:                                 checked with previous assembler
            m_Em            membrane strain tensor.                             V
            m_Em_der        first variation of E_m.                             V
            m_Em_der2       second variation of E_m MULTIPLIED BY S_m.          V
            m_Ef            flexural strain tensor.                             V
            m_Ef_der        second variation of E_f.                            V
            m_Ef_der2       second variation of E_f MULTIPLIED BY S_f.          V

        Where:
            m_ori           geometry map of the the undeformed geometry,
            m_def           geometry map of the deformed geometry,
            m_materialMat   the material matrix,
            m_m2            an auxillary matrix to multiply the last row of a tensor with 2
            m_space         solution space
    **/

    space m_space = m_assembler.trialSpace(0);
    geometryMap m_ori = m_assembler.exprData()->getMap();
    geometryMap m_def = m_assembler.exprData()->getMap2();
    variable m_force = m_assembler.getCoeff(*m_forceFun);
    variable m_thick = m_assembler.getCoeff(*m_thickFun);

    auto m_Em       = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
    auto m_Sm      = m_Em * reshape(m_materialMat,3,3);
    auto m_N        = m_thick.val() * m_Sm;
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Sm_der  = m_Em_der * reshape(m_materialMat,3,3);
    auto m_N_der    = m_thick.val() * m_Sm_der;
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
    auto m_Sf      = m_Ef * reshape(m_materialMat,3,3);
    auto m_M        = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf;
    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Sf_der  = m_Ef_der * reshape(m_materialMat,3,3);
    auto m_M_der    = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf_der;
    auto m_Ef_der2  = flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
                        + var2(m_space,m_space,m_def,m_Ef * reshape(m_materialMat,3,3) );

    auto m_ff       = m_force;
    // Mat & Rhs
}

template<class T>
void gsThinShellAssembler<T>::assemble()
{
    // Initialize stystem
    m_assembler.initSystem();

    gsMaterialMatrix materialMat(m_patches, *m_YoungsModulus, *m_PoissonsRatio);
    variable m_materialMat = m_assembler.getCoeff(materialMat);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force    = m_assembler.getCoeff(*m_forceFun);
    variable m_thick    = m_assembler.getCoeff(*m_thickFun);

    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Sm_der   = m_Em_der * reshape(m_materialMat,3,3);
    auto m_N_der    = m_thick.val() * m_Sm_der;
    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Sf_der   = m_Ef_der * reshape(m_materialMat,3,3);
    auto m_M_der    = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf_der;

    // assemble system
    m_assembler.assemble(
        ( m_N_der * m_Em_der.tr() + m_M_der * m_Ef_der.tr() ) * meas(m_ori)
        ,m_space * m_force * meas(m_ori)
        );
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

    m_defpatches = deformed;

    // Initialize matrix
    m_assembler.initMatrix();

    gsMaterialMatrix materialMat(m_patches, *m_YoungsModulus, *m_PoissonsRatio);
    variable m_materialMat = m_assembler.getCoeff(materialMat);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_thick    = m_assembler.getCoeff(*m_thickFun);

    auto m_Em       = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
    auto m_Sm       = m_Em * reshape(m_materialMat,3,3);
    auto m_N        = m_thick.val() * m_Sm;
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Sm_der   = m_Em_der * reshape(m_materialMat,3,3);
    auto m_N_der    = m_thick.val() * m_Sm_der;
    auto m_Em_der2  = flatdot( jac(m_space),jac(m_space).tr(), m_N ); //[checked]

    auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
    auto m_Sf       = m_Ef * reshape(m_materialMat,3,3);
    auto m_M        = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf;
    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Sf_der   = m_Ef_der * reshape(m_materialMat,3,3);
    auto m_M_der    = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf_der;
    auto m_Ef_der2  = flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
                        + var2(m_space,m_space,m_def, m_M );

    // Assemble matrix
    m_assembler.assemble(
                ( m_N_der * m_Em_der.tr() + m_Em_der2 + m_M_der * m_Ef_der.tr() - m_Ef_der2 ) * meas(m_ori)
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
    m_defpatches = deformed;

    // Initialize vector
    m_assembler.initVector();

    gsMaterialMatrix materialMat(m_patches, *m_YoungsModulus, *m_PoissonsRatio);
    variable m_materialMat = m_assembler.getCoeff(materialMat);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space m_space       = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    variable m_force    = m_assembler.getCoeff(*m_forceFun);
    variable m_thick    = m_assembler.getCoeff(*m_thickFun);

    auto m_Em       = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
    auto m_Sm       = m_Em * reshape(m_materialMat,3,3);
    auto m_N        = m_thick.val() * m_Sm;
    auto m_Em_der   = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]

    auto m_Ef       = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
    auto m_Sf       = m_Ef * reshape(m_materialMat,3,3);
    auto m_M        = m_thick.val() * m_thick.val() * m_thick.val() / 12.0 * m_Sf;
    auto m_Ef_der   = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(m_space,m_def) ) ) * reshape(m_m2,3,3); //[checked]

    // Assemble vector
    m_assembler.assemble(m_space * m_force * meas(m_ori) -
                ( ( m_N * m_Em_der.tr() - m_M * m_Ef_der.tr() ) * meas(m_ori) ).tr()
                );
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
        assembleMatrix(deformed);
    }
    assembleVector(deformed);
}
template<class T>
void gsThinShellAssembler<T>::assemble(const gsMatrix<T> & solVector,
                                        bool Matrix)
{
    constructSolution(solVector, m_defpatches);
    assemble(m_defpatches,Matrix);
}

template <class T>
gsMultiPatch<T> gsThinShellAssembler<T>::constructSolution(const gsMatrix<T> & solVector) const
{
    gsMultiPatch<T> mp = m_patches;

    m_solution.setSolutionVector(solVector);

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
        displacement.patch(k).coefs() -= m_patches.patch(0).coefs();;  // defG points to mp_def, therefore updated
    }

    return displacement;
}

template <class T>
void gsThinShellAssembler<T>::constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const
{
    deformed = constructDisplacement(solVector);
}

}// namespace gismo
