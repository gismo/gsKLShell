/** @file gsThinShellAssemblerDWR.hpp

    @brief Provides DWR assembly routines for the Kirchhoff-Love shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/gsMaterialMatrixEval.h>

namespace gismo
{

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::gsThinShellAssemblerDWR(
                                                                const gsMultiPatch<T> & patches,
                                                                const gsMultiBasis<T> & basisL,
                                                                const gsMultiBasis<T> & basisH,
                                                                const gsBoundaryConditions<T> & bconditions,
                                                                const gsFunction<T> & surface_force,
                                                                gsMaterialMatrixBase<T> * materialmatrix
                                                            )
                                                            :
                                                            Base(patches,basisL,bconditions,surface_force,materialmatrix),
                                                            m_basisL(basisL),
                                                            m_basisH(basisH)
{
    m_assemblerL = new gsThinShellAssembler<d,T,bending>(patches,basisL,bconditions,surface_force,materialmatrix);
    m_assemblerH = new gsThinShellAssembler<d,T,bending>(patches,basisH,bconditions,surface_force,materialmatrix);

    m_dL = gsVector<T>::Zero(m_assemblerL->numDofs());
    m_dH = gsVector<T>::Zero(m_assemblerH->numDofs());
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
void gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_setBasis(const gsMultiBasis<T> & basis)
{
    Base::m_basis = basis;
    Base::_initialize();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assembleMass(gsThinShellAssemblerBase<T> * assembler, bool lumped)
{
    assembler->assembleMass(lumped);
    return assembler->matrix();
}

// template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
// gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assembleMassLM(gsThinShellAssemblerBase<T> * assembler, bool lumped)
// {
//     gsExprAssembler<T> assembler(1,1);
//     // Elements used for numerical integration

//     assembler.setIntegrationElements(m_basisH);
//     assembler.setOptions(m_assemblerL->options());

//     // Initialize the geometry maps
//     assembler.getMap(m_patches);           // this map is used for integrals

//     // Set the discretization space
//     space spaceL = assembler.getSpace(m_basisL, d, 0); // last argument is the space ID
//     space spaceH = assembler.getTestSpace(spaceL , m_basisH);

//     this->_assembleDirichlet();

//     // Initialize stystem
//     assembler.initSystem(false);

//     gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,m_patches);
//     variable mm0 = assembler.getCoeff(m_mm);

//     space       m_space = assembler.trialSpace(0);
//     geometryMap m_ori   = assembler.exprData()->getMap();

//     // assemble system
//     assembler.assemble(mm0.val()*spaceL*spaceH.tr()*meas(m_ori));
//     return assembler.matrix();

    // m_assembler.cleanUp();
    // m_assembler.setOptions(m_options);

    // m_assembler.getMap(m_patches);           // this map is used for integrals

    // // Initialize stystem
    // m_assembler.initSystem(false);

    // gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,m_patches);
    // variable mm0 = m_assembler.getCoeff(m_mm);

    // space       m_space = m_assembler.trialSpace(0);
    // geometryMap m_ori   = m_assembler.exprData()->getMap();

    // // assemble system
    // if (!lumped)
    //     m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
    // else
    //     m_assembler.assemble(mm0.val()*(m_space.rowSum())*meas(m_ori));
// }

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assembleMatrix(gsThinShellAssemblerBase<T> * assembler)
{
    assembler->assemble();
    return assembler->matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assembleMatrix(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed)
{
    assembler->assemble();
    assembler->assembleMatrix(deformed);
    return assembler->matrix(); //Base::matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assemblePrimal(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed)
{
    assembler->assembleVector(deformed);
    return assembler->rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assemblePrimal(gsThinShellAssemblerBase<T> * assembler)
{
    assembler->assemble();
    return assembler->rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space   space   = exprAssembler.trialSpace(0);
    variable usol   = exprAssembler.getCoeff(primal);
    geometryMap Gori= exprAssembler.exprData()->getMap();

    auto expr = 2 * space * usol * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());
    exprAssembler.getMap(m_patches);           // this map is used for integrals
    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = space * gismo::expr::uv(_comp,3) * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
    // variable Lambda = exprAssembler.getCoeff(lambdaf);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m2,3,3);
    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m2,3,3);
    auto lambda     = Cm     * reshape(Tmat,3,3).tr();
    auto lambda_der = Cm_der * reshape(Tmat,3,3).tr();

    gsExprEvaluator<T> ev(exprAssembler);

    auto expr = 2 * lambda_der * lambda.tr() * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");
    GISMO_ASSERT(_comp!=2 ,"Not available for comp = 2");

    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m2,3,3);
    auto lambda_der = Cm_der * reshape(Tmat,3,3).tr();

    gsExprEvaluator<T> ev(exprAssembler);
    gsVector<T> pt(2);
    pt.setConstant(0.25);
    gsDebugVar(ev.eval(flat( jac(Gdef).tr() * jac(space) ),pt));
    gsDebugVar(ev.eval(m2,pt));
    gsDebugVar(ev.eval(Cm_der,pt));
    gsDebugVar(ev.eval(reshape(Tmat,3,3).tr(),pt));

    auto expr = lambda_der * gismo::expr::uv(_comp,3) * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto expr = 2 * Em_der * Em.tr() * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto expr = Em_der *  gismo::expr::uv(_comp,3) * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);


    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto Emp    = Em * reshape(Tmat,3,3).tr();
    auto Emp_der= Em_der * reshape(Tmat,3,3).tr();


    auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);

    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);


    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto Emp    = Em * reshape(Tmat,3,3).tr();
    auto Emp_der= Em_der * reshape(Tmat,3,3).tr();


    // auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
    auto expr = Emp_der * gismo::expr::uv(_comp,3) * meas(Gori);

    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable S0 = exprAssembler.getCoeff(S0f);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der     = Em_der * reshape(mmA,3,3);
    auto expr = 2 * Sm_der * S0 * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    // gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    // variable S0 = exprAssembler.getCoeff(S0f);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der     = Em_der * reshape(mmA,3,3);
    auto expr = Sm_der * gismo::expr::uv(_comp,3) * meas(Gori);
    exprAssembler.assemble(expr);
    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable P0  = exprAssembler.getCoeff(P0f);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // gsFunctionExpr<> Tmatf("1","0","0","0","2","0","0","0","3",2);
    // variable Tmat = exprAssembler.getCoeff(Tmatf);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der     = Em_der * reshape(mmA,3,3);
    // auto expr = (Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(_comp,3) * meas(Gori);
    auto expr = 2 * (Sm_der * reshape(Tmat,3,3).tr()) * P0 * meas(Gori);
    exprAssembler.assemble(expr);

    // gsExprEvaluator<T> ev(exprAssembler);
    // gsVector<T> pt(2);
    // pt.setConstant(0.25);
    // gsDebugVar(ev.eval((Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(0,3) * meas(Gori),pt));
    // gsDebugVar(ev.eval((Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(1,3) * meas(Gori),pt));
    // gsDebugVar(ev.eval((Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(2,3) * meas(Gori),pt));

    return exprAssembler.rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // gsFunctionExpr<> Tmatf("1","0","0","0","2","0","0","0","3",2);
    // variable Tmat = exprAssembler.getCoeff(Tmatf);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der     = Em_der * reshape(mmA,3,3);
    auto expr = (Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(_comp,3) * meas(Gori);
    exprAssembler.assemble(expr);

    // gsExprEvaluator<T> ev(exprAssembler);
    // gsVector<T> pt(2);
    // pt.setConstant(0.25);
    // gsDebugVar(ev.eval((Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(0,3) * meas(Gori),pt));
    // gsDebugVar(ev.eval((Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(1,3) * meas(Gori),pt));
    // gsDebugVar(ev.eval((Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(2,3) * meas(Gori),pt));

    return exprAssembler.rhs();
}

// template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
// template<enum GoalFunction _GF, short_t _comp>
// typename std::enable_if<_GF==GoalFunction::FlexuralStress, void>::type
// gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
// {
//     this->_setBasis(basis);

//     gsMatrix<T> Z(1,1);
//     Z.setZero();

//     exprAssembler.cleanUp();
//     exprAssembler.setOptions(assembler->options());

//     m_defpatches = m_patches;
//     for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
//         m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

//     exprAssembler.getMap(m_patches);           // this map is used for integrals
//     exprAssembler.getMap(m_defpatches);

//     // Initialize vector
//     exprAssembler.initSystem(false);
//     exprAssembler.initVector(1,false);

//     gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(m_materialMat,m_defpatches,Z);
//     gsMaterialMatrixEval<T,MaterialOutput::MatrixB> mmBf(m_materialMat,m_defpatches,Z);
//     variable mmA = exprAssembler.getCoeff(mmAf);
//     variable mmB = exprAssembler.getCoeff(mmBf);

//     gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
//     variable m2 = exprAssembler.getCoeff(mult2t);

//     space  space       = exprAssembler.trialSpace(0);
//     geometryMap Gori   = exprAssembler.exprData()->getMap();
//     geometryMap Gdef   = exprAssembler.exprData()->getMap2();

//     auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);

//     auto Sf_der     = Ef_der * reshape(mmB,3,3);
//     auto expr = Sm_der * gismo::expr::uv(_comp,3) * meas(Gori);
//     exprAssembler.assemble(expr);

//     // _assembleDual_expr(expr,deformed,basis);
// }

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf(assembler->material(),m_defpatches);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable mmB = exprAssembler.getCoeff(mmBf);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);
    // auto N_der    = Em_der * reshape(mmA,3,3)
    //                 +
    //                 Ef_der * reshape(mmB,3,3);

    gsWarn<<"The MembraneForce goal function seems to be problematic (Nder)\n";
    auto N_der1     = Em_der * reshape(mmA,3,3);
    auto N_der2     = Ef_der * reshape(mmB,3,3);
    auto N_der      = N_der1 + N_der2;

    // auto expr = N_der * gismo::expr::uv(_comp,3) * meas(Gori);
    auto expr = N_der * gismo::expr::uv(_comp,3) * meas(Gori);
    auto expr1 = N_der1 * gismo::expr::uv(_comp,3) * meas(Gori);
    auto expr2 = N_der2 * gismo::expr::uv(_comp,3) * meas(Gori);
    exprAssembler.assemble(expr);


    gsDebugVar(exprAssembler.rhs().rows());


    gsVector<T> pt(2);
    pt.setConstant(0.25);
    gsExprEvaluator<T> ev(exprAssembler);
    gsDebugVar(ev.eval(Gdef,pt));
    gsDebugVar(ev.eval(Em_der,pt));
    gsDebugVar(ev.eval(Ef_der,pt));
    gsDebugVar(ev.eval(Em_der * reshape(mmA,3,3),pt));
    gsDebugVar(ev.eval(Ef_der * reshape(mmB,3,3),pt));
    gsDebugVar(ev.eval(N_der,pt));
    gsDebugVar(ev.eval(expr,pt));
    gsDebugVar(ev.eval(expr1,pt));
    gsDebugVar(ev.eval(expr2,pt));

    gsDebugVar(exprAssembler.rhs());

    return exprAssembler.rhs();
}

// template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
// gsVector<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::_assembleDual_point(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMatrix<T> & tmp)
// {
//     /* SINGLE PATCH IMPLEMENTATION */
//     index_t pIndex = 0;
//     gsVector<T> result;
//     gsMatrix<index_t> actives, globalActives;

//     space  space       = exprAssembler.trialSpace(0);

//     result.resize(exprAssembler.numDofs());
//     result.setZero();

//     for (index_t k = 0; k!=u.cols(); k++)
//     {
//         tmp = ev.eval(expr,u.col(k));
//         basis.front().basis(pIndex).active_into( u.col(k), actives ); // note: takes the first patch of pids as patch id. Works for non-overlapping patches.
//         for (index_t j = 0; j< space.dim(); ++j)
//         {
//             space.mapper().localToGlobal(actives, pIndex, globalActives,j);
//             for (index_t i=0; i < globalActives.rows(); ++i)
//                 if (v.mapper().is_free_index(globalActs(i,0)))
//                     result.at(globalActives(i,0)) += tmp(i,0);
//         }

//     }

//     return result;
// }

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space   space   = exprAssembler.trialSpace(0);
    variable usol   = exprAssembler.getCoeff(primal);
    geometryMap Gori= exprAssembler.exprData()->getMap();

    auto expr = 2 * space * usol * meas(Gori);
    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());
    exprAssembler.getMap(m_patches);           // this map is used for integrals
    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = space * gismo::expr::uv(_comp,3) * meas(Gori);
    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
    variable lambda = exprAssembler.getCoeff(lambdaf);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Cm_der     = 2*flat( jac(Gdef).tr() * jac(space) ) * reshape(m2,3,3);
    auto lambda_der = Cm_der * reshape(Tmat,3,3).tr();

    auto expr = 2* lambda_der * lambda * meas(Gori);//<<<--------NOTE: Square missing?
    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");
    GISMO_ASSERT(_comp!=2 ,"Not available for comp = 2");

    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Cm_der     = 2*flat( jac(Gdef).tr() * jac(space) ) * reshape(m2,3,3);
    auto lambda_der = Cm_der * reshape(Tmat,3,3).tr();

    auto expr = lambda_der * gismo::expr::uv(_comp,3) * meas(Gori);
    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto expr = 2 * Em_der * Em.tr() * meas(Gori);
    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto expr = Em_der * gismo::expr::uv(_comp,3) * meas(Gori);
    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto Emp    = Em * reshape(Tmat,3,3).tr();
    auto Emp_der= Em_der * reshape(Tmat,3,3).tr();

    auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto Emp    = Em * reshape(Tmat,3,3).tr();
    auto Emp_der= Em_der * reshape(Tmat,3,3).tr();

    // auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
    auto expr = Emp_der * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable S0 = exprAssembler.getCoeff(S0f);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der    = Em_der * reshape(mmA,3,3);

    auto expr = 2 * Sm_der * S0 * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    // gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    // variable S0 = exprAssembler.getCoeff(S0f);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der    = Em_der * reshape(mmA,3,3);

    // auto expr = 2 * Sm_der * S0 * meas(Gori);
    auto expr = Sm_der * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable P0  = exprAssembler.getCoeff(P0f);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // gsFunctionExpr<> Tmatf("1","0","0","0","2","0","0","0","3",2);
    // variable Tmat = exprAssembler.getCoeff(Tmatf);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der    = Em_der * reshape(mmA,3,3);

    auto expr = (Sm_der * reshape(Tmat,3,3).tr()) * P0 * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(assembler->material(),m_defpatches,Z);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // gsFunctionExpr<> Tmatf("1","0","0","0","2","0","0","0","3",2);
    // variable Tmat = exprAssembler.getCoeff(Tmatf);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Sm_der    = Em_der * reshape(mmA,3,3);

    auto expr = (Sm_der * reshape(Tmat,3,3).tr()) * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp==9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf(assembler->material(),m_defpatches);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable mmB = exprAssembler.getCoeff(mmBf);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->material(),m_defpatches,Z);
    variable S0 = exprAssembler.getCoeff(S0f);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);
    auto N_der    = Em_der * reshape(mmA,3,3) + Ef_der * reshape(mmB,3,3);

    auto expr = 2 * N_der * S0 * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp!=9, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::assembleDual_impl(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf(assembler->material(),m_defpatches);
    variable mmA = exprAssembler.getCoeff(mmAf);
    variable mmB = exprAssembler.getCoeff(mmBf);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    space  space       = exprAssembler.trialSpace(0);
    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto Em_der   = flat( jac(Gdef).tr() * jac(space) ) ;
    auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);
    auto N_der    = Em_der * reshape(mmA,3,3) + Ef_der * reshape(mmB,3,3);

    auto expr = N_der * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    result.resize(exprAssembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        exprAssembler.integrationElements().basis(pIndex).active_into( points.col(k), actives );
        for (index_t j = 0; j< space.dim(); ++j)
        {
            space.mapper().localToGlobal(actives, pIndex, globalActives,j);
            for (index_t i=0; i < globalActives.rows(); ++i)
                if (space.mapper().is_free_index(globalActives(i,0)))
                    result.at(globalActives(i,0)) += tmp(i+j*actives.rows(),0);
        }
    }
    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
void gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    result = m_assemblerL->constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
void gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    result = m_assemblerH->constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
void gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed)
{
    m_assemblerL->constructSolution(solVector,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
void gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed)
{
    m_assemblerH->constructSolution(solVector,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructSolutionL(const gsMatrix<T> & solVector)
{
    return m_assemblerL->constructSolution(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructSolutionH(const gsMatrix<T> & solVector)
{
    return m_assemblerH->constructSolution(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructDisplacementL(const gsMatrix<T> & solVector)
{
    return m_assemblerL->constructDisplacement(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF, comp>::constructDisplacementH(const gsMatrix<T> & solVector)
{
    return m_assemblerH->constructDisplacement(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
T gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori= exprAssembler.exprData()->getMap();
    variable F      = exprAssembler.getCoeff(*m_forceFun, Gori);
    variable zsolL  = exprAssembler.getCoeff(dualL);
    variable zsolH  = exprAssembler.getCoeff(dualH);

    auto expr = (zsolH-zsolL).tr() * F * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );

    // _assembleDual_expr(expr,deformed,basis);
}


template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
T gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    return computeError_impl<d,bending>(dualL,dualH,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());
    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches);
    variable S0  = exprAssembler.getCoeff(S0f);
    variable S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    geometryMap Gori  = exprAssembler.exprData()->getMap();
    geometryMap Gdef  = exprAssembler.exprData()->getMap2();
    variable F      = exprAssembler.getCoeff(*m_forceFun, Gori);
    variable zsolL  = exprAssembler.getCoeff(dualL);
    variable zsolH  = exprAssembler.getCoeff(dualH);
        // variable m_thick = exprAssembler.getCoeff(*m_thickFun, m_ori);

    this->homogenizeDirichlet();

    auto N        = S0.tr();
    // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
    auto Em_der   = flat( jac(Gdef).tr() * (fjac(zsolH) - fjac(zsolL)) ) ;

    auto M        = S1.tr(); // output is a column
    // auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);
    auto Ef_der   = -( deriv2(zsolH,sn(Gdef).normalized().tr() )
                    -  deriv2(zsolL,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(zsolH,Gdef) )
                                                                 - deriv2(Gdef,var1(zsolL,Gdef) ) ) * reshape(m2,3,3);

    auto Fext = (zsolH-zsolL).tr() * F * meas(Gori);
    auto Fint = ( N * Em_der.tr() + M * Ef_der.tr() ) * meas(Gori);

    auto expr = ( Fext - Fint ) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T integral = ev.integral(expr);

    gsVector<T> pt(2);
    pt.setConstant(0.25);

    gsDebugVar(ev.eval(Gori.tr(),pt));
    gsDebugVar(ev.eval(Gdef.tr(),pt));
    gsDebugVar(ev.eval(zsolL.tr(),pt));
    gsDebugVar(ev.eval(zsolH.tr(),pt));

    gsDebugVar(ev.eval(N.tr(),pt));
    gsDebugVar(ev.eval(M.tr(),pt));

    gsDebugVar(ev.eval(Em_der.tr(),pt));
    gsDebugVar(ev.eval(Ef_der.tr(),pt));

    real_t fint_m = ev.integral(( N * Em_der.tr() ) * meas(Gori) );
    real_t fint_f = ev.integral(( M * Ef_der.tr() ) * meas(Gori) );
    real_t fint = fint_m - fint_f;
    real_t fext = ev.integral( Fext  );

    gsDebugVar(fint_m);
    gsDebugVar(fint_f);

    gsDebugVar(fint);
    gsDebugVar(fext);
    real_t Res = fext-fint;
    gsDebugVar(Res);

    if (m_foundInd)
    {
        variable foundation = exprAssembler.getCoeff(*m_foundFun, Gori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        integral += ev.integral( ( zsolH - zsolL ) * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    }
    if (m_pressInd)
    {
        variable pressure = exprAssembler.getCoeff(*m_pressFun, Gori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        integral += ev.integral( pressure.val() * ( zsolH - zsolL ) * sn(Gdef).normalized() * meas(Gori) );
    }

    return integral;

}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());
    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches);
    variable S0  = exprAssembler.getCoeff(S0f);

    geometryMap Gori  = exprAssembler.exprData()->getMap();
    geometryMap Gdef  = exprAssembler.exprData()->getMap2();
    variable F      = exprAssembler.getCoeff(*m_forceFun, Gori);

    variable zsolL  = exprAssembler.getCoeff(dualL);
    variable zsolH  = exprAssembler.getCoeff(dualH);
        // variable m_thick = exprAssembler.getCoeff(*m_thickFun, m_ori);

    this->homogenizeDirichlet();

    auto N        = S0.tr();
    // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
    auto Em_der   = flat( jac(Gdef).tr() * (fjac(zsolH) - fjac(zsolL)) ) ;

    auto Fext = (zsolH-zsolL).tr() * F * meas(Gori);
    auto Fint = ( N * Em_der.tr() ) * meas(Gori);

    auto expr = ( Fext - Fint ) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T integral = ev.integral(expr);

    if (m_foundInd)
    {
        variable foundation = exprAssembler.getCoeff(*m_foundFun, Gori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        integral += ev.integral( ( zsolH - zsolL ) * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    }
    if (m_pressInd)
    {
        variable pressure = exprAssembler.getCoeff(*m_pressFun, Gori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        integral += ev.integral( pressure.val() * ( zsolH - zsolL ) * sn(Gdef).normalized() * meas(Gori) );
    }

    return integral;

}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF>
typename std::enable_if<(_GF==GoalFunction::Buckling), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeErrorEig_impl( const T evPrimalL,
                                                                        const T evDualL,
                                                                        const T evDualH,
                                                                        const gsMultiPatch<T> & dualL,
                                                                        const gsMultiPatch<T> & dualH,
                                                                        const gsMultiPatch<T> & primal,
                                                                        const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;
    // for ( size_t k =0; k!=deformed.nPatches(); ++k) // Deform the geometry
    //     m_defpatches.patch(k).coefs() += deformed.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_L(m_assemblerL->material(),m_patches);
    variable mmA_L    = exprAssembler.getCoeff(mmAf_L);
    variable mmB_L    = exprAssembler.getCoeff(mmBf_L);
    variable mmC_L    = exprAssembler.getCoeff(mmCf_L);
    variable mmD_L    = exprAssembler.getCoeff(mmDf_L);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_NL(m_assemblerL->material(),m_defpatches);
    variable mmA_NL    = exprAssembler.getCoeff(mmAf_NL);
    variable mmB_NL    = exprAssembler.getCoeff(mmBf_NL);
    variable mmC_NL    = exprAssembler.getCoeff(mmCf_NL);
    variable mmD_NL    = exprAssembler.getCoeff(mmDf_NL);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerL->material(),m_defpatches);
    variable S0  = exprAssembler.getCoeff(S0f);
    variable S1  = exprAssembler.getCoeff(S1f);
    variable rho  = exprAssembler.getCoeff(rhof);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    variable usol   = exprAssembler.getCoeff(primal);
    variable zsolL  = exprAssembler.getCoeff(dualL);
    variable zsolH  = exprAssembler.getCoeff(dualH);

    geometryMap Gori  = exprAssembler.exprData()->getMap();
    geometryMap Gdef  = exprAssembler.exprData()->getMap2();


    // // Matrix K_L
    // auto N      = S0.tr();
    // auto Em_der = flat( jac(Gori).tr() * (fjac(zsolH) - fjac(zsolL)) ) ;

    // auto M      = S1.tr(); // output is a column
    // auto Ef_der = -(deriv2(zsolH, sn(Gori).normalized().tr()) - deriv2(zsolL, sn(Gori).normalized().tr()) + deriv2(Gori, var1(zsolH, Gori)) - deriv2(Gori, var1(zsolL, Gori))) * reshape(m2, 3, 3);

    // auto N_der  = Em_der * reshape(mmA,3,3) + Ef_der * reshape(mmB,3,3);
    // auto M_der  = Em_der * reshape(mmC,3,3) + Ef_der * reshape(mmD,3,3);

    // auto Fint   = ( N_der * Em_der.tr() + M_der * Ef_der.tr() ) * meas(Gori);

    // // Matrix K_NL
    // auto Em_der2  = flatdot( jac(Gdef),jac(zsolH), N ) - flatdot( jac(Gdef),jac(zsolL), N );
    // auto Ef_der2  = -(flatdot2( deriv2(m_space), var1(m_space,m_def).tr(), m_M  ).symmetrize()
    //                         + var2dot(m_space,m_space,m_def, m_M ))

    /*
        Weak form:
            Wint(u,v) = int  eps(u,v)'*n(u) + kappa(u,v)'*m(u) dint

        Jacobian
            Wint(u,v,w) = int  eps(u,v)'*n(u,w)' + eps''(u,v,w)*n(u) + kappa(u,v)'*m'(u,w) + kappa''(u,v,w)*m(u) dint

     */


    auto N_L        = S0.tr();
    auto Em_der_L   = flat(  jac(Gori).tr() * ( fjac(usol )               ) ) ;
    auto Em_derd_L  = flat(  jac(Gori).tr() * ( fjac(zsolH) - fjac(zsolL) ) ) ;

    auto M_L        = S1.tr(); // output is a column
    auto Ef_der_L   = -( deriv2(usol , sn(Gori).normalized().tr() ) - deriv2(Gori, var1(usol , Gori) ) ) * reshape(m2, 3, 3);
    auto Ef_derd_L  = -(
                         deriv2(zsolH, sn(Gori).normalized().tr() ) - deriv2(Gori, var1(zsolH, Gori) )
                        -deriv2(zsolL, sn(Gori).normalized().tr() ) - deriv2(Gori, var1(zsolL, Gori) )
                       ) * reshape(m2, 3, 3);

    auto N_der_L    = Em_der_L * reshape(mmA_L,3,3) + Ef_der_L * reshape(mmB_L,3,3);
    auto M_der_L    = Em_der_L * reshape(mmC_L,3,3) + Ef_der_L * reshape(mmD_L,3,3);

    auto N_derd_L   = Em_derd_L * reshape(mmA_L,3,3) + Ef_derd_L * reshape(mmB_L,3,3);
    auto M_derd_L   = Em_derd_L * reshape(mmC_L,3,3) + Ef_derd_L * reshape(mmD_L,3,3);

    auto Fint_L     = ( N_derd_L * Em_derd_L.tr() + M_derd_L * Ef_derd_L.tr() ) * meas(Gori);

    auto N_NL       = S0.tr();
    auto Em_der_NL  = flat(  jac(Gdef).tr() * ( fjac(usol )               ) ) ;
    auto Em_derd_NL = flat(  jac(Gdef).tr() * ( fjac(zsolH) - fjac(zsolL) ) ) ;

    auto Em_der2    = flat( fjac(usol).tr() * ( fjac(usol)                ) ) * N_NL.tr();
    auto Emd_der2   = flat( fjac(usol).tr() * ( fjac(zsolH) - fjac(zsolL) ) ) * N_NL.tr();

    auto M_NL       = S1.tr(); // output is a column
    auto Ef_der_NL  = -( deriv2(usol , sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(usol , Gdef) ) ) * reshape(m2, 3, 3);
    auto Ef_derd_NL = -(
                         deriv2(zsolH, sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(zsolH, Gdef) )
                        -deriv2(zsolL, sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(zsolL, Gdef) )
                       ) * reshape(m2, 3, 3);

    auto Ef_der2    = -(
                          (  ( deriv2(usol )                 )  * ( var1(usol ,Gdef)                    ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (  ( deriv2(usol )                 )  * ( var1(usol ,Gdef)                    ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (    deriv2(Gdef)                     * ( var2(usol ,usol ,Gdef)              )       ).tr()      * reshape(m2, 3, 3)
                       ) * M_NL.tr();

    auto Efd_der2   = -(
                          (  ( deriv2(usol )                 )  * ( var1(zsolH,Gdef) - var1(zsolL,Gdef)             ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (  ( deriv2(zsolH) - deriv2(zsolL) )  * ( var1(usol ,Gdef)                                ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (    deriv2(Gdef)                     * ( var2(usol ,zsolH,Gdef) - var2(usol ,zsolL,Gdef) )       ).tr()      * reshape(m2, 3, 3)
                       ) * M_NL.tr();


    auto N_der_NL   = Em_der_NL  * reshape(mmA_NL,3,3) + Ef_der_NL  * reshape(mmB_NL,3,3);
    auto M_der_NL   = Em_der_NL  * reshape(mmC_NL,3,3) + Ef_der_NL  * reshape(mmD_NL,3,3);

    auto N_derd_NL  = Em_derd_NL * reshape(mmA_NL,3,3) + Ef_derd_NL * reshape(mmB_NL,3,3);
    auto M_derd_NL  = Em_derd_NL * reshape(mmC_NL,3,3) + Ef_derd_NL * reshape(mmD_NL,3,3);

    auto Fint_NL    = ( N_derd_NL * Em_der_NL.tr() + M_derd_NL * Ef_der_NL.tr() + Em_der2 + Ef_der2 ) * meas(Gori);

    auto K_L        = ( N_der_L   * Em_der_L.tr()  + M_der_L   * Ef_der_L.tr()                          ) * meas(Gori);
    auto Kd_L       = ( N_derd_L  * Em_der_L.tr()  + M_derd_L  * Ef_der_L.tr()                          ) * meas(Gori);
    auto K_NL       = ( N_der_NL  * Em_der_NL.tr() + M_der_NL  * Ef_der_NL.tr() + Em_der2  + Ef_der2    ) * meas(Gori);
    auto Kd_NL      = ( N_derd_NL * Em_der_NL.tr() + M_derd_NL * Ef_der_NL.tr() + Emd_der2 + Efd_der2   ) * meas(Gori);

    auto A          = Kd_L;
    auto Bdiff      = Kd_NL - Kd_L;
    auto Bprimal    = K_NL  - K_L ;

    // gsVector<T> pt(2);
    // pt.setConstant(0.25);

    gsExprEvaluator<T> ev(exprAssembler);

    // gsDebugVar(ev.eval(Gori.tr(), pt));
    // gsDebugVar(ev.eval(Gdef.tr(), pt));

    // gsDebug<<"E_mder_L:\n";
    // gsDebugVar("\n"<<ev.eval(jac(Gori).tr(), pt));
    // gsDebugVar("\n"<<ev.eval(fjac(usol), pt));
    // gsDebugVar("\n"<<ev.eval(jac(Gori).tr() * fjac(usol), pt));
    // gsDebugVar("\n"<<ev.eval(flat(jac(Gori).tr() * fjac(usol)), pt));

    // gsDebug<<"E_mder_NL:\n";
    // gsDebugVar("\n"<<ev.eval(jac(Gdef).tr(), pt));
    // gsDebugVar("\n"<<ev.eval(fjac(usol), pt));
    // gsDebugVar("\n"<<ev.eval(jac(Gdef).tr() * fjac(usol), pt));
    // gsDebugVar("\n"<<ev.eval(flat(jac(Gdef).tr() * fjac(usol)), pt));

    // gsDebug<<"E_fder_L:\n";
    // gsDebugVar("\n"<<ev.eval(deriv2(usol, sn(Gori).normalized().tr()), pt));
    // gsDebugVar("\n"<<ev.eval(deriv2(Gori, var1(usol, Gori)), pt));
    // gsDebugVar("\n"<<ev.eval(-(deriv2(usol, sn(Gori).normalized().tr()) - deriv2(Gori, var1(usol, Gori))) * reshape(m2, 3, 3), pt));

    // gsDebug<<"E_fder_NL:\n";
    // gsDebugVar("\n"<<ev.eval(deriv2(usol, sn(Gdef).normalized().tr()), pt));
    // gsDebugVar("\n"<<ev.eval(deriv2(Gdef, var1(usol, Gdef)), pt));
    // gsDebugVar("\n"<<ev.eval(-(deriv2(usol, sn(Gdef).normalized().tr()) - deriv2(Gdef, var1(usol, Gdef))) * reshape(m2, 3, 3), pt));

    // gsDebug<<"E_mder2:\n";
    // gsDebugVar("\n"<<ev.eval(flat( fjac(usol).tr() * fjac(usol) ), pt));
    // gsDebugVar("\n"<<ev.eval(N_NL.tr(), pt));
    // gsDebugVar("\n"<<ev.eval(flat( fjac(usol).tr() * fjac(usol) ) * N_NL.tr(), pt));

    // gsDebug<<"E_fder2:\n";
    // gsDebugVar("\n"<<ev.eval((deriv2(usol) * var1(usol,Gdef).tr()).tr() * reshape(m2, 3, 3), pt));
    // gsDebugVar("\n"<<ev.eval(N_NL, pt));
    // gsDebugVar("\n"<<ev.eval((deriv2(usol) * var1(usol,Gdef).tr()).tr() * reshape(m2, 3, 3) * N_NL.tr(), pt));

    // gsDebugVar("\n"<<ev.eval((deriv2(usol) * var2(usol,usol,Gdef) ).tr() * reshape(m2, 3, 3), pt));
    // gsDebugVar("\n"<<ev.eval(M_NL.tr(), pt));
    // gsDebugVar("\n"<<ev.eval((deriv2(usol) * var2(usol,usol,Gdef) ).tr() * reshape(m2, 3, 3) * M_NL.tr(), pt));

    // gsDebugVar("Em_der_L = "<< ev.eval(Em_der_L, pt));
    // gsDebugVar("Ef_der_L = "<< ev.eval(Ef_der_L, pt));

    // gsDebugVar("Em_der_NL = "<< ev.eval(Em_der_NL, pt));
    // gsDebugVar("Ef_der_NL = "<< ev.eval(Ef_der_NL, pt));

    // gsDebugVar("N_der_L = "<< ev.eval(N_der_L, pt));
    // gsDebugVar("M_der_L = "<< ev.eval(M_der_L, pt));

    // gsDebugVar("N_der_NL = "<< ev.eval(N_der_NL, pt));
    // gsDebugVar("M_der_NL = "<< ev.eval(M_der_NL, pt));

    // gsDebugVar("Em_der2 = "<< ev.eval(Em_der2, pt));
    // gsDebugVar("Ef_der2 = "<< ev.eval(Ef_der2, pt));

    // gsDebugVar(ev.eval((deriv2(usol) * var2(usol,usol,Gdef) ).tr() * reshape(m2, 3, 3) * M_NL.tr(), pt));



    // gsDebugVar(ev.eval(Fint_NL, pt));

    // gsInfo << "Fint_m = " << ev.integral((N_der_L * Em_der_L.tr()) * meas(Gori)));
    // gsInfo << "Fint_f = " << ev.integral((M_der_L * Ef_der_L.tr()) * meas(Gori)));

    // gsInfo << "Fint_m = " << ev.integral((N_der_NL * Em_der_NL.tr() + Em_der2) * meas(Gori)));
    // gsInfo << "Fint_f = " << ev.integral((M_der_NL * Ef_der_NL.tr() + Ef_der2) * meas(Gori)));



    // gsInfo << "Fint_L  = " << ev.integral(Fint_L));
    // gsInfo << "Fint_NL = " << ev.integral(Fint_NL));



    // auto A      = Fint;
    // auto Bdiff  = rho.val() * usol.tr() * (zsolH - zsolL);
    // auto Bprimal= rho.val() * usol.tr() * usol;

    T integral = 0.0;
    integral += evPrimalL * ev.integral( Bdiff * meas(Gori));
    integral += -ev.integral( A  * meas(Gori) );
    integral += (evDualH - evDualL) * (ev.integral( Bprimal * meas(Gori)) - 1.0);

    return integral;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF>
typename std::enable_if<(_GF==GoalFunction::Modal), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeErrorEig_impl( const T evPrimalL,
                                                                        const T evDualL,
                                                                        const T evDualH,
                                                                        const gsMultiPatch<T> & dualL,
                                                                        const gsMultiPatch<T> & dualH,
                                                                        const gsMultiPatch<T> & primal)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // defG points to mp_def, therefore updated

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf(m_assemblerL->material(),m_patches);
    variable mmA    = exprAssembler.getCoeff(mmAf);
    variable mmB    = exprAssembler.getCoeff(mmBf);
    variable mmC    = exprAssembler.getCoeff(mmCf);
    variable mmD    = exprAssembler.getCoeff(mmDf);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerL->material(),m_patches);
    variable S0  = exprAssembler.getCoeff(S0f);
    variable S1  = exprAssembler.getCoeff(S1f);
    variable rho  = exprAssembler.getCoeff(rhof);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    variable usol   = exprAssembler.getCoeff(primal);
    variable zsolL  = exprAssembler.getCoeff(dualL);
    variable zsolH  = exprAssembler.getCoeff(dualH);

    geometryMap Gori  = exprAssembler.exprData()->getMap();
    geometryMap Gdef  = exprAssembler.exprData()->getMap2();

    auto N      = S0.tr();
    auto Em_derd= flat( jac(Gori).tr() * (fjac(zsolH) - fjac(zsolL)) ) ;
    auto Em_der = flat( jac(Gori).tr() * (fjac(usol) ) ) ;

    auto M      = S1.tr(); // output is a column
    auto Ef_derd= -(deriv2(zsolH, sn(Gori).normalized().tr()) - deriv2(zsolL, sn(Gori).normalized().tr()) + deriv2(Gori, var1(zsolH, Gori)) - deriv2(Gori, var1(zsolL, Gori))) * reshape(m2, 3, 3);
    auto Ef_der = -(deriv2(usol, sn(Gori).normalized().tr()) + deriv2(Gori, var1(usol, Gori)) ) * reshape(m2, 3, 3);

    auto N_der  = Em_derd * reshape(mmA,3,3) + Ef_derd * reshape(mmB,3,3);
    auto M_der  = Em_derd * reshape(mmC,3,3) + Ef_derd * reshape(mmD,3,3);

    auto Fint   = ( N_der * Em_der.tr() + M_der * Ef_der.tr() ) * meas(Gori);

    gsVector<T> pt(2);
    pt.setConstant(0.25);

    gsExprEvaluator<T> ev(exprAssembler);

    auto A      = Fint;
    auto Bdiff  = rho.val() * usol.tr() * (zsolH - zsolL);
    auto Bprimal= rho.val() * usol.tr() * usol;

    T integral = 0.0;
    integral += evPrimalL * ev.integral( Bdiff * meas(Gori));
    integral += -ev.integral( A  * meas(Gori) );
    integral += (evDualH - evDualL) * (ev.integral( Bprimal * meas(Gori)) - 1.0);

    // gsDebugVar(ev.eval(Gori.tr(), pt));
    // gsDebugVar(ev.eval(Gdef.tr(), pt));
    // gsDebugVar(ev.eval(zsolL.tr(), pt));
    // gsDebugVar(ev.eval(zsolH.tr(), pt));
    // gsDebugVar(ev.eval(usol.tr(), pt));
    // gsDebugVar(ev.eval(N_der, pt));
    // gsDebugVar(ev.eval(M_der, pt));
    // gsDebugVar(ev.eval(Fint, pt));
    // gsDebugVar(ev.eval(Fint * meas(Gori), pt));
    // gsInfo << "Fint_m = " << ev.integral((N_der * Em_der.tr()) * meas(Gori)));
    // gsInfo << "Fint_f = " << ev.integral((M_der * Ef_der.tr()) * meas(Gori)));

    return integral;
}

//============================================================
template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto expr = (Gdef - Gori).tr() * (Gdef - Gori) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto expr = (Gdef-Gori).tr() * gismo::expr::uv(_comp,3)*meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    /*
        NOTE:
        - Cm * Transformation does only give the in-plane stretches, since C is constructed with the in-plane jacobians (Gdef) only. In fact, Cm should have a third entry with the Jacobian Determinant J0, but this one is hard to compute the variations for.
     */

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto cov_ori = jac(Gori).tr()*jac(Gori);
    auto cov_def = jac(Gdef).tr()*jac(Gdef);

    auto Cm     = flat(jac(Gdef).tr()*jac(Gdef)) * reshape(m2,3,3);
    auto stretch= reshape(Tmat,3,3) * Cm.tr();

    auto expr = stretch.tr() * stretch * meas(Gori);

    // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
    // variable lambda = exprAssembler.getCoeff(lambdaf);
    // auto expr1 = gismo::expr::pow(lambda.tr() * lambda,2) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");
    GISMO_ASSERT(_comp!=2 ,"Not available for comp = 2");

    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Cm     = flat(jac(Gdef).tr()*jac(Gdef)) * reshape(m2,3,3);
    auto stretch= reshape(Tmat,3,3) * Cm.tr();

    auto expr = stretch.tr() * gismo::expr::uv(_comp,3) * meas(Gori);

    // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
    // variable lambda = exprAssembler.getCoeff(lambdaf);
    // auto expr = gismo::expr::pow(lambda.tr() * gismo::expr::uv(_comp,3),2) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    auto expr = Em * Em.tr() * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    auto expr = Em * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Emp    = Em * reshape(Tmat,3,3).tr();

    // auto expr = Emp * Emp.tr() * meas(Gori);
    auto expr = Emp * Emp.tr() * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Emp    = Em * reshape(Tmat,3,3).tr();

    // auto expr = Emp * Emp.tr() * meas(Gori);
    auto expr = Emp * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches,Z);
    variable S0  = exprAssembler.getCoeff(S0f);

    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = S0.tr() * S0 * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches,Z);
    variable S0  = exprAssembler.getCoeff(S0f);

    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = S0.tr() * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_materialMat,m_defpatches,Z);
    variable P0     = exprAssembler.getCoeff(P0f);

    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = P0.tr() * P0 * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    T result = ev.integral( expr  );

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_materialMat,m_defpatches,Z);
    variable P0     = exprAssembler.getCoeff(P0f);

    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = P0.tr() * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    T result = ev.integral( expr  );

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    variable N  = exprAssembler.getCoeff(Nf);

    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = N.tr() * N * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    variable N  = exprAssembler.getCoeff(Nf);

    geometryMap Gori   = exprAssembler.exprData()->getMap();

    auto expr = N.tr() * gismo::expr::uv(_comp,3) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);

    return ev.integral( expr  );
}

//============================================================
template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto expr = (Gdef - Gori).tr() * (Gdef - Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Displacement && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    auto expr = (Gdef-Gori).tr() * gismo::expr::uv(_comp,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    // exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    // gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    // variable m2 = exprAssembler.getCoeff(mult2t);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    // geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    // gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    // variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // auto Cm     = flat(jac(Gdef).tr()*jac(Gdef)) * reshape(m2,3,3);
    // auto stretch= reshape(Tmat,3,3) * Cm.tr();

    // auto expr = stretch.tr() * gismo::expr::uv(_comp,3) * meas(Gori);

    gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
    variable lambda = exprAssembler.getCoeff(lambdaf);
    auto expr = gismo::expr::pow(lambda.tr() * lambda,2) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::Stretch && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    // ! Only works with Z=0 for now
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    // exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    // To compensate for the flat, since flat does [E11,E22,2*E12]
    // gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","0.5",2);
    // variable m2 = exprAssembler.getCoeff(mult2t);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    // geometryMap Gdef   = exprAssembler.exprData()->getMap2();

    // gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_assemblerL->material(),m_defpatches,Z);
    // variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // auto Cm     = flat(jac(Gdef).tr()*jac(Gdef)) * reshape(m2,3,3);
    // auto stretch= reshape(Tmat,3,3) * Cm.tr();

    // auto expr = stretch.tr() * gismo::expr::uv(_comp,3) * meas(Gori);

    gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
    variable Lambda = exprAssembler.getCoeff(lambdaf);
    auto expr = gismo::expr::pow(Lambda.tr() * gismo::expr::uv(_comp,3),2) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    auto expr = Em.tr() * Em;

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStrain && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    auto expr = Em.tr() * gismo::expr::uv(_comp,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_materialMat,m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Emp    = Em * reshape(Tmat,3,3).tr();

    auto expr = Emp.tr() * Emp;

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStrain && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_materialMat,m_defpatches,Z);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    geometryMap Gori   = exprAssembler.exprData()->getMap();
    geometryMap Gdef   = exprAssembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Emp    = Em * reshape(Tmat,3,3).tr();

    // auto expr = Emp.tr() * Emp;
    auto expr = Emp.tr() * gismo::expr::uv(_comp,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_materialMat,m_defpatches,Z);
    variable S0  = exprAssembler.getCoeff(S0f);

    auto expr = S0.tr() * S0;

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneStress && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_materialMat,m_defpatches,Z);
    variable S0  = exprAssembler.getCoeff(S0f);

    auto expr = S0.tr() * gismo::expr::uv(_comp,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_materialMat,m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_materialMat,m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_materialMat,m_defpatches,Z);
    variable S0  = exprAssembler.getCoeff(S0f);
    variable P0  = exprAssembler.getCoeff(P0f);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // gsFunctionExpr<> Tmatf("1","0","0","0","2","0","0","0","3",2);
    // variable Tmat = exprAssembler.getCoeff(Tmatf);

    // auto expr = (S0.tr() * reshape(Tmat,3,3).tr()) * (S0.tr() * reshape(Tmat,3,3).tr()).tr();
    auto expr = P0.tr() * P0;

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembranePStress && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_materialMat,m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_materialMat,m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::Transformation>  Tmatf(m_materialMat,m_defpatches,Z);
    variable S0  = exprAssembler.getCoeff(S0f);
    variable P0  = exprAssembler.getCoeff(P0f);
    variable Tmat   = exprAssembler.getCoeff(Tmatf);

    // gsFunctionExpr<> Tmatf("1","0","0","0","2","0","0","0","3",2);
    // variable Tmat = exprAssembler.getCoeff(Tmatf);

    // auto expr = (S0.tr() * reshape(Tmat,3,3).tr()) * gismo::expr::uv(_comp,3);
    auto expr = P0.tr() * gismo::expr::uv(_comp,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp==9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    variable N  = exprAssembler.getCoeff(Nf);

    auto expr = N.tr() * N;

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template<enum GoalFunction _GF, short_t _comp>
typename std::enable_if<_GF==GoalFunction::MembraneForce && _comp!=9, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    GISMO_ASSERT(_comp < d,"Component is out of bounds");

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    variable N  = exprAssembler.getCoeff(Nf);

    auto expr = N.tr() * gismo::expr::uv(_comp,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval( expr, points.col(k));
        GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
        result += tmp(0,0);
    }

    return result;
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template <enum GoalFunction _GF>
typename std::enable_if<_GF == GoalFunction::Modal, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::matrixNorm_impl(const gsMultiPatch<T> &primal, const gsMultiPatch<T> &other) const
{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerH->options());

    exprAssembler.getMap(m_patches); // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1, false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerH->material(),m_defpatches);
    variable rho = exprAssembler.getCoeff(rhof);

    variable usol = exprAssembler.getCoeff(primal);
    variable zsol = exprAssembler.getCoeff(other);

    geometryMap Gori = exprAssembler.exprData()->getMap();

    gsExprEvaluator<T> ev(exprAssembler);
    return ev.integral(rho.val() * usol.tr() * zsol * meas(Gori));
}

template <short_t d, class T, bool bending, enum GoalFunction GF, short_t comp>
template <enum GoalFunction _GF>
typename std::enable_if<_GF == GoalFunction::Buckling, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF, comp>::matrixNorm_impl(const gsMultiPatch<T> &primal, const gsMultiPatch<T> &other, const gsMultiPatch<T> &deformed) const
{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerH->options());

    gsMultiPatch<T> defpatches = deformed;

    exprAssembler.getMap(m_patches);           // this map is used for integrals
    exprAssembler.getMap(defpatches);

    // Initialize vector
    exprAssembler.initSystem(false);
    exprAssembler.initVector(1,false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_L(m_assemblerL->material(),m_patches);
    variable mmA_L    = exprAssembler.getCoeff(mmAf_L);
    variable mmB_L    = exprAssembler.getCoeff(mmBf_L);
    variable mmC_L    = exprAssembler.getCoeff(mmCf_L);
    variable mmD_L    = exprAssembler.getCoeff(mmDf_L);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_NL(m_assemblerL->material(),m_defpatches);
    variable mmA_NL    = exprAssembler.getCoeff(mmAf_NL);
    variable mmB_NL    = exprAssembler.getCoeff(mmBf_NL);
    variable mmC_NL    = exprAssembler.getCoeff(mmCf_NL);
    variable mmD_NL    = exprAssembler.getCoeff(mmDf_NL);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->material(),defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->material(),defpatches);
    variable S0  = exprAssembler.getCoeff(S0f);
    variable S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = exprAssembler.getCoeff(mult2t);

    geometryMap Gori  = exprAssembler.exprData()->getMap();
    geometryMap Gdef  = exprAssembler.exprData()->getMap2();
    variable usol = exprAssembler.getCoeff(primal);
    variable zsol = exprAssembler.getCoeff(other);

    auto Em_der_L   = flat(  jac(Gori).tr() * ( fjac(usol )               ) ) ;
    auto Ef_der_L   = -( deriv2(usol , sn(Gori).normalized().tr() ) - deriv2(Gori, var1(usol , Gori) ) ) * reshape(m2, 3, 3);

    auto Em_derd_L  = flat(  jac(Gori).tr() * ( fjac(zsol )               ) ) ;
    auto Ef_derd_L  = -( deriv2(zsol , sn(Gori).normalized().tr() ) - deriv2(Gori, var1(zsol , Gori) ) ) * reshape(m2, 3, 3);


    auto N_der_L     = Em_derd_L * reshape(mmA_L,3,3) + Ef_derd_L * reshape(mmB_L,3,3);
    auto M_der_L     = Em_derd_L * reshape(mmC_L,3,3) + Ef_derd_L * reshape(mmD_L,3,3);

    auto N_NL       = S0.tr();
    auto Em_der_NL  = flat(  jac(Gdef).tr() * ( fjac(usol )               ) ) ;
    auto Ef_der_NL  = -( deriv2(usol , sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(usol , Gdef) ) ) * reshape(m2, 3, 3);

    auto Em_derd_NL = flat(  jac(Gdef).tr() * ( fjac(zsol )               ) ) ;
    auto Ef_derd_NL = -( deriv2(zsol , sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(zsol , Gdef) ) ) * reshape(m2, 3, 3);

    auto Em_der2    = flat( fjac(usol).tr() * ( fjac(zsol )               ) ) * N_NL.tr();

    auto M_NL       = S1.tr(); // output is a column
    auto Ef_der2    = -(
                          (  ( deriv2(zsol )                 )  * ( var1(usol ,Gdef)                    ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (  ( deriv2(usol )                 )  * ( var1(zsol ,Gdef)                    ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (    deriv2(Gdef)                     * ( var2(usol ,zsol ,Gdef)              )       ).tr()      * reshape(m2, 3, 3)
                       ) * M_NL.tr();

    auto N_der_NL   = Em_derd_NL  * reshape(mmA_NL,3,3) + Ef_derd_NL  * reshape(mmB_NL,3,3);
    auto M_der_NL   = Em_derd_NL  * reshape(mmC_NL,3,3) + Ef_derd_NL  * reshape(mmD_NL,3,3);

    auto K_L        = ( N_der_L   * Em_der_L.tr()  + M_der_L   * Ef_der_L.tr()                          ) * meas(Gori);
    auto K_NL       = ( N_der_NL  * Em_der_NL.tr() + M_der_NL  * Ef_der_NL.tr() + Em_der2  + Ef_der2    ) * meas(Gori);

    auto Bprimal    = K_NL - K_L ;

    auto expr = ( Bprimal ) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T integral = ev.integral(expr);

    if (m_foundInd)
    {
        variable foundation = exprAssembler.getCoeff(*m_foundFun, Gori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        integral += ev.integral( zsol * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    }
    if (m_pressInd)
    {
        variable pressure = exprAssembler.getCoeff(*m_pressFun, Gori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        integral += ev.integral( pressure.val() * zsol * sn(Gdef).normalized() * meas(Gori) );
    }

    return integral;
}

} // namespace gismo
