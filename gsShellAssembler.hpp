/** @file gsShellAssembler.hpp

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

#include <gsThinShell2/gsShellAssembler.h>
#include <gsThinShell2/gsMaterialMatrix.h>
#include <gsThinShell2/gsShellUtils.h>

namespace gismo
{

template<class T>
gsShellAssembler<T>::gsShellAssembler(  const gsMultiPatch<T> & patches,
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
                                        m_surfforce(surface_force),
                                        m_thickness(thickness),
{
    m_Youngs = gsConstantFunction(YoungsModulus,3);
    m_Youngs = gsConstantFunction(PoissonsRatio,3);
}

template<class T>
gsShellAssembler<T>::gsShellAssembler(  const gsMultiPatch<T> & patches,
                                        const gsMultiBasis<T> & basis,
                                        const gsBoundaryConditions<T> & bconditions,
                                        const gsFunction<T> & surface_force,
                                        const gsFunction<T> & thickness,
                                        const gsFunction<T> & YoungsModulus,
                                        const gsFunction<T> & PoissonsRatio,
                                        )
                                        :
                                        m_patches(patches),
                                        m_basis(basis),
                                        m_bcs(bconditions),
                                        m_surfforce(surface_force),
                                        m_thickness(thickness),
                                        m_Youngs(YoungsModulus),
                                        m_Poisson(PoissonsRatio),
{

}

template <class T>
gsOptionList gsElasticityAssembler<T>::defaultOptions()
{
    gsOptionList & opt = m_assembler.options();
    // opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    // opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    // opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::saint_venant_kirchhoff);
    return opt;
}

/*
    TO INITIALIZE
        basis
        boundary conditions
*/

template <class T>
void gsShellAssembler<T>::initialize()
{
    //gsInfo<<"Active options:\n"<< m_assembler.options() <<"\n";


    // Elements used for numerical integration
    m_assembler.setIntegrationElements(dbasis);
    m_evaluator = gsExprEvaluator<> ev(m_assembler);

    // Set the geometry map
    m_ori = m_assembler.getMap(mp);           // this map is used for integrals
    m_def = m_assembler.getMap(mp_def);

    // Set the discretization space
    m_space = m_assembler.getSpace(dbasis, 3);
    m_space.setInterfaceCont(0); // todo: 1 (smooth basis)
    m_space.addBc( bc.get("Dirichlet") ); // (!) must be called only once

    // Solution vector and solution variable
    m_solution = m_assembler.getSolution(m_space, m_solvector);

    // Define fields as variables:
    // ... surface force
    m_force = m_assembler.getCoeff(m_forceFun, m_ori);
    // ... thickness
    m_thick = m_assembler.getCoeff(m_thickFun,m_ori);
}

template <class T>
void gsShellAssembler<T>::defineComponents()
{
    gsMaterialMatrix materialMat(m_patches, m_young, m_poiss);
    m_materialMat = m_assembler.getCoeff(materialMat);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    m_m2 = m_assembler.getCoeff(mult2t, m_ori);

    // Define components
    // m_def is geometry map of deformed shell
    // m_ori is geometry map of undeformed shell
    // m_space is the solution space
    auto m_Em = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) ; //[checked]
    auto m_Em_der = flat( jac(m_def).tr() * jac(m_space) ) ; //[checked]
    auto m_Em_der2 = flatdot( jac(m_space),jac(m_space).tr(), m_Em * reshape(m_materialMat,3,3) ); //[checked]

    auto m_Ef = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) ; //[checked]
    auto m_Ef_der = ( deriv2(m_space,sn(m_def).normalized().tr() ) + deriv2(m_def,var1(u,m_def) ) ) * reshape(m_m2,3,3); //[checked]
    auto m_Ef_der2 = flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_Ef * reshape(m_materialMat,3,3)  ).symmetrize()
                    + var2(m_space,m_space,m_def,m_Ef * reshape(m_materialMat,3,3) );
}

template<class T>
void gsShellAssembler<T>::assemble()
{
    // Initialize stystem
    m_assembler.initSystem();

    // assemble system
    m_assembler.assemble(
        (
            (m_thick.val()) * (m_Em_der * reshape(mm,3,3) * m_Em_der.tr())
            +
            (m_thick.val() * m_thick.val() * m_thick.val())/3.0 * (m_Ef_der * reshape(m_materialMat,3,3) * m_Ef_der.tr())
        ) * meas(m_ori)
        ,m_space * m_force * meas(m_ori)
        );
}

// TO DO
// template <class T>
// bool gsShellAssembler<T>::assemble(const gsMatrix<T> & solutionVector,
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
bool gsShellAssembler<T>::assembleMatrix()
{
    // Initialize matrix
    m_assembler.initMatrix();

    // Assemble matrix
    m_assembler.assemble(
                (
                    (m_thick.val()) * (m_Em_der * reshape(mm,3,3) * m_Em_der.tr() + m_Em_der2)
                    +
                    (m_thick.val() * m_thick.val() * m_thick.val())/3.0 * (m_Ef_der * reshape(m_materialMat,3,3) * m_Ef_der.tr() - m_Ef_der2)
                ) * meas(m_ori)
                );
}

template <class T>
bool gsShellAssembler<T>::assembleVector()
{
    // Initialize vector
    m_assembler.initVector();

    // Assemble vector
    m_assembler.assemble(m_space * m_force * meas(m_ori) -
                (
                 (
                    (m_thick.val()) *(m_Em * reshape(m_materialMat,3,3) * m_Em_der.tr()) -
                    (m_thick.val() * m_thick.val() * m_thick.val())/3.0 * (m_Ef * reshape(m_materialMat,3,3) * m_Ef_der.tr())
                 ) * meas(m_ori)
                ).tr()
                );
}

// TO DO
template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & displacement,
                                        bool assembleMatrix)
{
    constructSolution(displacement);

    if (assembleMatrix)
    {
        assembleMatrix(displacement);
    }
    assembleVector(displacement);
}

template <class T>
void gsElasticityAssembler<T>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) const
{

}


}// namespace gismo
