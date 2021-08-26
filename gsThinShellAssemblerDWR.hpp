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

#include <gsCore/gsFunctionExpr.h>

namespace gismo
{

template <short_t d, class T, bool bending>
gsThinShellAssemblerDWR<d, T, bending>::gsThinShellAssemblerDWR(
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
                                                            m_basisH(basisH),
                                                            m_bcs(bconditions)
{
    m_assemblerL = new gsThinShellAssembler<d,T,bending>(patches,basisL,bconditions,surface_force,materialmatrix);
    m_assemblerH = new gsThinShellAssembler<d,T,bending>(patches,basisH,bconditions,surface_force,materialmatrix);

    m_dL = gsVector<T>::Zero(m_assemblerL->numDofs());
    m_dH = gsVector<T>::Zero(m_assemblerH->numDofs());
 }

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::_setBasis(const gsMultiBasis<T> & basis)
{
    Base::m_basis = basis;
    Base::_initialize();
}

template <short_t d, class T, bool bending>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleMass(gsThinShellAssemblerBase<T> * assembler, bool lumped)
{
    assembler->assembleMass(lumped);
    return assembler->matrix();
}

// template <short_t d, class T, bool bending>
// gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleMassLM(gsThinShellAssemblerBase<T> * assembler, bool lumped)
// {
//     gsExprAssembler<T> assembler(1,1);
//     // Elements used for numerical integration

//     assembler.setIntegrationElements(m_basisH);
//     assembler.setOptions(m_assemblerL->options());

//     // Initialize the geometry maps
//     geometryMap m_ori   = assembler.getMap(m_patches);           // this map is used for integrals

//     // Set the discretization space
//     space spaceL = assembler.getSpace(m_basisL, d, 0); // last argument is the space ID
//     space spaceH = assembler.getTestSpace(spaceL , m_basisH);

//     this->_assembleDirichlet();

//     // Initialize stystem
//     assembler.initSystem();

//     gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,m_patches);
//     auto mm0 = assembler.getCoeff(m_mm);

//     space       m_space = assembler.trialSpace(0);

//     // assemble system
//     assembler.assemble(mm0.val()*spaceL*spaceH.tr()*meas(m_ori));
//     return assembler.matrix();

    // m_assembler.cleanUp();
    // m_assembler.setOptions(m_options);

    // m_assembler.getMap(m_patches);           // this map is used for integrals

    // // Initialize stystem
    // m_assembler.initSystem();

    // gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,m_patches);
    // auto mm0 = m_assembler.getCoeff(m_mm);

    // space       m_space = m_assembler.trialSpace(0);
    // geometryMap m_ori   = m_assembler.exprData()->getMap();

    // // assemble system
    // if (!lumped)
    //     m_assembler.assemble(mm0.val()*m_space*m_space.tr()*meas(m_ori));
    // else
    //     m_assembler.assemble(mm0.val()*(m_space.rowSum())*meas(m_ori));
// }

template <short_t d, class T, bool bending>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleMatrix(gsThinShellAssemblerBase<T> * assembler)
{
    assembler->assemble();
    return assembler->matrix();
}

template <short_t d, class T, bool bending>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleMatrix(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed)
{
    assembler->assemble();
    assembler->assembleMatrix(deformed);
    return assembler->matrix(); //Base::matrix();
}

template <short_t d, class T, bool bending>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::_assemblePrimal(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed)
{
    assembler->assembleVector(deformed);
    return assembler->rhs();
}

template <short_t d, class T, bool bending>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::_assemblePrimal(gsThinShellAssemblerBase<T> * assembler)
{
    assembler->assemble();
    return assembler->rhs();
}

template <short_t d, class T, bool bending>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    GISMO_ASSERT(m_component < d || m_component==9,"Component is out of bounds");

    // For the stretches
    // ! Only works with Z=0 for now, since principal stresses and directions are implemented for Z=0, since principal stresses and directions are implemented for Z=0
    gsMatrix<T> Z(1,1);
    Z.setZero();
    ////

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Space
    space   space   = exprAssembler.trialSpace(0);
    // Homogenize Dirichlet
    space.setup(m_bcs, dirichlet::homogeneous, 0);

    // Solution
    auto usol   = exprAssembler.getCoeff(primal);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::CovTransform>  Tcovf(m_assemblerL->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::ConTransform>  Tconf(m_assemblerL->material(),m_defpatches,Z);
    auto Tcov = exprAssembler.getCoeff(Tcovf);
    auto Tcon = exprAssembler.getCoeff(Tconf);

    // Material tensors at Z
    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::MatrixD> mmDf(assembler->material(),m_defpatches,Z);
    auto mmA = exprAssembler.getCoeff(mmAf);
    auto mmD = exprAssembler.getCoeff(mmDf);

    // Material tensors integrated
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAfi(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBfi(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCfi(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(assembler->material(),m_defpatches);
    auto mmAi = exprAssembler.getCoeff(mmAfi);
    auto mmBi = exprAssembler.getCoeff(mmBfi);
    auto mmCi = exprAssembler.getCoeff(mmCfi);
    auto mmDi = exprAssembler.getCoeff(mmDfi);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_assemblerL->material(),m_defpatches,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Force and moment tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerL->material(),m_defpatches);
    auto N = exprAssembler.getCoeff(Nf);
    auto M = exprAssembler.getCoeff(Mf);

    // Helper matrix for flexural components
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient and its first variation
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m12,3,3);

    // Principal stretch and its first variation
    auto lambda     = Cm     * reshape(Tcon,3,3).tr();
    auto lambda_der = Cm_der * reshape(Tcon,3,3).tr();

    // Membrane strain and its first variation
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    // Principal membrane strain and its first variation
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tcon,3,3).tr();
    auto Emp_der= (Em_der * reshape(m12,3,3)) * reshape(Tcon,3,3).tr();

    // Flexural strain and its first variation
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;
    auto Ef_der = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);

    // Membrane stress (its first variation)
    auto Sm_der     = Em_der * reshape(mmA,3,3);
    auto Sf_der     = Ef_der * reshape(mmD,3,3);

    // Membrane force (its first variation)
    auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
    auto M_der    = Em_der * reshape(mmCi,3,3) + Ef_der * reshape(mmDi,3,3);

    gsExprEvaluator<T> ev(exprAssembler);
    gsVector<T> pt(2);
    pt.setConstant(0.25);

    if (m_goalFunction == GoalFunction::Displacement)
    {
        if (m_component==9)
        {
            auto expr = 2 * space * usol * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = space * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        if (m_component==9)
        {
            auto expr = 2 * lambda_der * lambda.tr() * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = lambda_der * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Em_der * Em.tr() * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = Em_der *  gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::MembranePStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = Emp_der * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sm_der * S0 * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = Sm_der * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::MembranePStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * (Sm_der * reshape(Tcov,3,3).tr()) * P0 * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = (Sm_der * reshape(Tcov,3,3).tr()) * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            auto expr = 2 * N_der * N * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            gsWarn<<"The MembraneForce goal function seems to be problematic (Nder)\n";
            auto N_der1     = Em_der * reshape(mmAi,3,3);
            auto N_der2     = Ef_der * reshape(mmBi,3,3);
            auto N_der0a     = N_der1;// + N_der2; // -->> commenting the second term out will give correct results.
            auto N_der0b     = N_der2;// + N_der2; // -->> commenting the second term out will give correct results.
            auto N_der0c     = N_der1 + N_der2; // -->> commenting the second term out will give correct results.

            auto expr0a = N_der0a * gismo::expr::uv(m_component,3) * meas(Gori);
            auto expr0b = N_der0b * gismo::expr::uv(m_component,3) * meas(Gori);
            auto expr0c = N_der0c * gismo::expr::uv(m_component,3) * meas(Gori);

            exprAssembler.initVector(1);
            exprAssembler.assemble(expr0a);
            gsDebugVar(exprAssembler.rhs());

            exprAssembler.initVector(1);
            exprAssembler.assemble(expr0b);
            gsDebugVar(exprAssembler.rhs());

            exprAssembler.initVector(1);
            exprAssembler.assemble(expr0c);
            gsDebugVar(exprAssembler.rhs());


            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Ef_der * Ef.tr() * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = Ef_der *  gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sf_der * S1 * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = Sf_der * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            auto expr = 2 * M_der * M * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
        else
        {
            auto expr = M_der * gismo::expr::uv(m_component,3) * meas(Gori);
            exprAssembler.assemble(expr);
            return exprAssembler.rhs();
        }
    }
    else
        GISMO_ERROR("Goal function not known");
}


template <short_t d, class T, bool bending>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;

    GISMO_ASSERT(m_component < d || m_component==9,"Component is out of bounds");

    // For the stretches
    // ! Only works with Z=0 for now, since principal stresses and directions are implemented for Z=0
    gsMatrix<T> Z(1,1);
    Z.setZero();
    ////

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(assembler->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Space
    space   space   = exprAssembler.trialSpace(0);
    // Homogenize Dirichlet
    space.setup(m_bcs, dirichlet::homogeneous, 0);

    // Solution
    auto usol   = exprAssembler.getCoeff(primal);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::CovTransform>  Tcovf(m_assemblerL->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::ConTransform>  Tconf(m_assemblerL->material(),m_defpatches,Z);
    auto Tcov = exprAssembler.getCoeff(Tcovf);
    auto Tcon = exprAssembler.getCoeff(Tconf);

    // Material tensors at Z
    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::MatrixD> mmDf(assembler->material(),m_defpatches,Z);
    auto mmA = exprAssembler.getCoeff(mmAf);
    auto mmD = exprAssembler.getCoeff(mmDf);

    // Material tensors integrated
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAfi(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBfi(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCfi(assembler->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(assembler->material(),m_defpatches);
    auto mmAi = exprAssembler.getCoeff(mmAfi);
    auto mmBi = exprAssembler.getCoeff(mmBfi);
    auto mmCi = exprAssembler.getCoeff(mmCfi);
    auto mmDi = exprAssembler.getCoeff(mmDfi);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S1f(assembler->material(),m_defpatches,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(assembler->material(),m_defpatches,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Force and moment tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerL->material(),m_defpatches);
    auto N = exprAssembler.getCoeff(Nf);
    auto M = exprAssembler.getCoeff(Mf);

    // Helper matrix for flexural components
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient and its first variation
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m12,3,3);

    // Principal stretch and its first variation
    auto lambda     = Cm     * reshape(Tcon,3,3).tr();
    auto lambda_der = Cm_der * reshape(Tcon,3,3).tr();

    // Membrane strain and its first variation
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    // Principal membrane strain and its first variation
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tcon,3,3).tr();
    auto Emp_der= (Em_der * reshape(m12,3,3)) * reshape(Tcon,3,3).tr();

    // Membrane stress (its first variation)
    auto Sm_der     = Em_der * reshape(mmA,3,3);

    // Principal membrane stress (its first variation)
    auto Smp_der     = Sm_der * reshape(Tcov,3,3).tr();

    // Flexural strain and its first variation
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;
    auto Ef_der = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);

    // Flexural stress (its first variation)
    auto Sf_der     = Ef_der * reshape(mmD,3,3);

    // Membrane force (its first variation)
    // auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
    auto M_der    = Em_der * reshape(mmCi,3,3) + Ef_der * reshape(mmDi,3,3);


    gsExprEvaluator<T> ev(exprAssembler);
    gsVector<T> pt(2);
    pt.setConstant(0.25);

    result.resize(exprAssembler.numDofs());
    result.setZero();
    if (m_goalFunction == GoalFunction::Displacement)
    {
        if (m_component==9)
        {
            auto expr = 2 * space * usol * meas(Gori);
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
        else
        {
            auto expr = space * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        if (m_component==9)
        {
            auto expr = 2* lambda_der * lambda * meas(Gori);//<<<--------NOTE: Square missing?
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
        else
        {
            auto expr = lambda_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Em_der * Em.tr() * meas(Gori);
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
        else
        {
            auto expr = Em_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::MembranePStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
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
        else
        {
            auto expr = Emp_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sm_der * S0 * meas(Gori);
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
        else
        {
            auto expr = Sm_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::MembranePStress)
    {
        if (m_component==9)
        {
            auto expr = Smp_der * P0 * meas(Gori);
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
        else
        {
            auto expr = Sm_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            gsWarn<<"The MembraneForce goal function seems to be problematic (Nder)\n";
            // remove N_der from this scipe
            auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
            auto expr = 2 * N_der * N * meas(Gori);
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
        else
        {
            // remove N_der from this scipe
            gsWarn<<"The MembraneForce goal function seems to be problematic (Nder)\n";
            auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
            auto expr = N_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Ef_der * Ef.tr() * meas(Gori);
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
        else
        {
            auto expr = Ef_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sf_der * S1 * meas(Gori);
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
        else
        {
            auto expr = Sf_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            gsWarn<<"The FlexuralMoment goal function seems to be problematic (Mder)\n";
            // remove N_der from this scipe
            auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
            auto expr = 2 * M_der * M * meas(Gori);
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
        else
        {
            // remove N_der from this scipe
            gsWarn<<"The FlexuralMoment goal function seems to be problematic (Mder)\n";
            auto expr = M_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
    }
    else
        GISMO_ERROR("Goal function not known");


};


// template <short_t d, class T, bool bending>
// gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleDual_point(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMatrix<T> & tmp)
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

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    result = m_assemblerL->constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    result = m_assemblerH->constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed)
{
    m_assemblerL->constructSolution(solVector,deformed);
}

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed)
{
    m_assemblerH->constructSolution(solVector,deformed);
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending>::constructSolutionL(const gsMatrix<T> & solVector)
{
    return m_assemblerL->constructSolution(solVector);
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending>::constructSolutionH(const gsMatrix<T> & solVector)
{
    return m_assemblerH->constructSolution(solVector);
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending>::constructDisplacementL(const gsMatrix<T> & solVector)
{
    return m_assemblerL->constructDisplacement(solVector);
}

template <short_t d, class T, bool bending>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending>::constructDisplacementH(const gsMatrix<T> & solVector)
{
    return m_assemblerH->constructDisplacement(solVector);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH)
{
    std::vector<T> result = computeError_impl<0>(dualL,dualH);
    GISMO_ASSERT(result.size()==1,"Must be scalar!");
    return result[0];
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH)
{
    return computeError_impl<1>(dualL,dualH);
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH)
{
    return computeError_impl<2>(dualL,dualH);
}

template <short_t d, class T, bool bending>
template <int _elWise>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    auto F      = exprAssembler.getCoeff(*m_forceFun, Gori);
    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);

    auto g_N = exprAssembler.getBdrFunction();

    auto expr = (zsolH-zsolL).tr() * F * meas(Gori);
    auto bexpr= (zsolH-zsolL).tr() * g_N * tv(Gori).norm();

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    std::vector<T> result;
    if (_elWise == 0)
    {
        // Get neumann sides, put them in a container of patchSide s.
        std::vector<patchSide> bContainer;
        typename gsBoundaryConditions<T>::bcContainer bcContainer = m_bcs.neumannSides();
        for (typename gsBoundaryConditions<T>::bcContainer::iterator it = bcContainer.begin(); it!=bcContainer.end(); ++it)
            bContainer.push_back(patchSide(it->patch(),it->side()));

        bintegral = ev.integralBdr( bexpr, bContainer);
        integral = ev.integral(expr);

        gsDebugVar(bintegral);
        gsDebugVar(integral);

        result.push_back(integral+bintegral);
    }
    else if (_elWise == 1) // element-wise
    {
        integral = ev.integralElWise(expr);
        gsDebugVar(integral);
        result = ev.elementwise();
    }
    else if (_elWise == 2) // function-wise
    {
        integral = 0;

        gsVector<T> pt(2);
        pt.setConstant(0.25);

        space   space   = exprAssembler.trialSpace(0);

        gsDebugVar(ev.eval(zsolL.tr(),pt));
        gsDebugVar(ev.eval(space,pt));
        gsDebugVar(ev.eval( space * zsolL,pt));

        // result =
    }
    else
        GISMO_ERROR("Unknown");

    return result;
}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    std::vector<T> result = computeError_impl<d,bending,0>(dualL,dualH,deformed);
    GISMO_ASSERT(result.size()==1,"Must be scalar!");
    return result[0];
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    return computeError_impl<d,bending,1>(dualL,dualH,deformed);
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    return computeError_impl<d,bending,2>(dualL,dualH,deformed);
}

template <short_t d, class T, bool bending>
template<int _d, bool _bending, int _elWise>
typename std::enable_if<_d==3 && _bending, std::vector<T>>::type
gsThinShellAssemblerDWR<d, T, bending>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());
    m_defpatches = deformed;

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    auto F      = exprAssembler.getCoeff(*m_forceFun, Gori);
    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);
        // auto m_thick = exprAssembler.getCoeff(*m_thickFun, m_ori);

    Base::homogenizeDirichlet();

    auto N        = S0.tr();
    // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
    auto Em_der   = flat( jac(Gdef).tr() * (jac(zsolH) - jac(zsolL)) ) ;

    auto M        = S1.tr(); // output is a column
    // auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);
    auto Ef_der   = -( deriv2(zsolH,sn(Gdef).normalized().tr() )
                    -  deriv2(zsolL,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(zsolH,Gdef) )
                                                                 - deriv2(Gdef,var1(zsolL,Gdef) ) ) * reshape(m2,3,3);

    auto Fext = (zsolH-zsolL).tr() * F;
    auto Fint = ( N * Em_der.tr() + M * Ef_der.tr() );

    auto g_N = exprAssembler.getBdrFunction();

    auto expr = ( Fext - Fint ) * meas(Gori);
    auto bexpr= (zsolH-zsolL).tr() * g_N * tv(Gori).norm();

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    std::vector<T> result;
    if (_elWise == 0)
    {
        // Get neumann sides, put them in a container of patchSide s.
        std::vector<patchSide> bContainer;
        typename gsBoundaryConditions<T>::bcContainer bcContainer = m_bcs.neumannSides();
        for (typename gsBoundaryConditions<T>::bcContainer::iterator it = bcContainer.begin(); it!=bcContainer.end(); ++it)
            bContainer.push_back(patchSide(it->patch(),it->side()));

        bintegral = ev.integralBdr( bexpr, bContainer);
        integral = ev.integral(expr);
    }
    else if (_elWise == 1) // element-wise
    {
        integral = ev.integralElWise(expr);
    }
    else if (_elWise == 2) // function-wise
    {
        gsVector<T> pt(2);
        pt.setConstant(0.25);

        space   space   = exprAssembler.trialSpace(0);
        auto Em_der_u = flat( jac(Gdef).tr() * ( (jac(zsolH) - jac(zsolL)) ) * space + jac(space) * (zsolH - zsolL) );
        // auto Ef_der_u = -( deriv2(zsolH,sn(Gdef).normalized().tr() )
        //                 -  deriv2(zsolL,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(zsolH,Gdef) )
        //                                                          - deriv2(Gdef,var1(zsolL,Gdef) ) ) * reshape(m2,3,3);



        gsDebugVar(ev.eval(Em_der_u,pt));
        gsDebugVar(ev.eval(zsolL * space,pt));

        exprAssembler.assemble( ( Fext - Fint ) * meas(Gori));
        std::vector<T> vec( exprAssembler.rhs().data(),
                            exprAssembler.rhs().data() + exprAssembler.rhs().rows() * exprAssembler.rhs().cols());
        integral = 0;
    }
    else
        GISMO_ERROR("Unknown");

    // if (m_foundInd)
    // {
    //     auto foundation = exprAssembler.getCoeff(*m_foundFun, Gori);
    //     GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

    //     integral += ev.integral( ( zsolH - zsolL ) * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    // }
    // if (m_pressInd)
    // {
    //     auto pressure = exprAssembler.getCoeff(*m_pressFun, Gori);
    //     GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

    //     integral += ev.integral( pressure.val() * ( zsolH - zsolL ) * sn(Gdef).normalized() * meas(Gori) );
    // }

    if (_elWise == 0)
        result.push_back(integral+bintegral);
    else if (_elWise == 1)
        result = ev.elementwise();
    else if (_elWise == 2)
    {

        // result =
    }
    else
        GISMO_ERROR("Unknown");

    return result;

}

template <short_t d, class T, bool bending>
template<int _d, bool _bending, int _elWise>
typename std::enable_if<!(_d==3 && _bending), std::vector<T>>::type
gsThinShellAssemblerDWR<d, T, bending>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());
    m_defpatches = deformed;

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches);
    auto S0  = exprAssembler.getCoeff(S0f);

    auto F      = exprAssembler.getCoeff(*m_forceFun, Gori);

    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);
        // auto m_thick = exprAssembler.getCoeff(*m_thickFun, m_ori);

    auto N        = S0.tr();
    // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
    auto Em_der   = flat( jac(Gdef).tr() * (jac(zsolH) - jac(zsolL)) ) ;

    auto Fext = (zsolH-zsolL).tr() * F;
    auto Fint = ( N * Em_der.tr() );

    auto g_N = exprAssembler.getBdrFunction();

    auto expr = ( Fext - Fint ) * meas(Gori);
    auto bexpr= (zsolH-zsolL).tr() * g_N * tv(Gori).norm();

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    std::vector<T> result;
    if (_elWise == 0)
    {
        // Get neumann sides, put them in a container of patchSide s.
        std::vector<patchSide> bContainer;
        typename gsBoundaryConditions<T>::bcContainer bcContainer = m_bcs.neumannSides();
        for (typename gsBoundaryConditions<T>::bcContainer::iterator it = bcContainer.begin(); it!=bcContainer.end(); ++it)
            bContainer.push_back(patchSide(it->patch(),it->side()));

        bintegral = ev.integralBdr( bexpr, bContainer);
        integral = ev.integral(expr);
    }
    else if (_elWise == 1) // element-wise
    {
        integral = ev.integralElWise(expr);
    }
    else if (_elWise == 2) // function-wise
    {
        gsVector<T> pt(2);
        pt.setConstant(0.25);

        space   space   = exprAssembler.trialSpace(0);
        auto Em_der_u = flat( jac(Gdef).tr() * ( (jac(zsolH) - jac(zsolL)) ) * space + jac(space) * (zsolH - zsolL) );
        // auto Ef_der_u = -( deriv2(zsolH,sn(Gdef).normalized().tr() )
        //                 -  deriv2(zsolL,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(zsolH,Gdef) )
        //                                                          - deriv2(Gdef,var1(zsolL,Gdef) ) ) * reshape(m2,3,3);



        gsDebugVar(ev.eval(zsolL * space,pt));
        gsDebugVar(ev.eval(Em_der_u,pt));

        // exprAssembler.assemble( ( Fext - Fint ) * meas(Gori));
        std::vector<T> vec( exprAssembler.rhs().data(),
                            exprAssembler.rhs().data() + exprAssembler.rhs().rows() * exprAssembler.rhs().cols());
        integral = 0;
    }
    else
        GISMO_ERROR("Unknown");
    // if (m_foundInd)
    // {
    //     auto foundation = exprAssembler.getCoeff(*m_foundFun, Gori);
    //     GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

    //     integral += ev.integral( ( zsolH - zsolL ) * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    // }
    // if (m_pressInd)
    // {
    //     auto pressure = exprAssembler.getCoeff(*m_pressFun, Gori);
    //     GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

    //     integral += ev.integral( pressure.val() * ( zsolH - zsolL ) * sn(Gdef).normalized() * meas(Gori) );
    // }

    if (_elWise == 0)
        result.push_back(integral+bintegral);
    else if (_elWise == 1)
        result = ev.elementwise();
    else if (_elWise == 2)
    {

        // result =
    }
    else
        GISMO_ERROR("Unknown");

    return result;
}


/**
 * @brief      Computes the eigenvalue error estimate for Buckling Analysis
 *
 * @param[in]  evPrimalL  The primal eigenvalue on the low order basis
 * @param[in]  evDualL    The dual eigenvalue on the low order basis
 * @param[in]  evDualH    The dual eigenvalue on the high order basis
 * @param[in]  dualL      The dual modeshape on the low order basis
 * @param[in]  dualH      The dual modeshape on the high order basis
 * @param[in]  primal     The physical mode shape
 * @param[in]  deformed   The deformed geometry around which the buckling solver is linearized
 *
 * @tparam     d          Dimension
 * @tparam     T          Real type
 * @tparam     bending    Bending flag
 *
 * @return     The error estimate
 */
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig(         const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                const gsMultiPatch<T> & deformed)
{
    std::vector<T> result = computeErrorEig_impl<0>(evPrimalL, evDualL, evDualH, dualL, dualH, primal, deformed);
    GISMO_ASSERT(result.size()==1,"Must be scalar!");
    return result[0];
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigElements( const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                const gsMultiPatch<T> & deformed)
{
    return computeErrorEig_impl<1>(evPrimalL, evDualL, evDualH, dualL, dualH, primal, deformed);
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigDofs( const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                const gsMultiPatch<T> & deformed)
{
    return computeErrorEig_impl<2>(evPrimalL, evDualL, evDualH, dualL, dualH, primal, deformed);
}


template <short_t d, class T, bool bending>
template <int _elWise>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig_impl(    const T evPrimalL,
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
    //     m_defpatches.patch(k).coefs() += deformed.patch(k).coefs();  // Gdef points to mp_def, therefore updated

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_L(m_assemblerL->material(),m_patches);
    auto mmA_L    = exprAssembler.getCoeff(mmAf_L);
    auto mmB_L    = exprAssembler.getCoeff(mmBf_L);
    auto mmC_L    = exprAssembler.getCoeff(mmCf_L);
    auto mmD_L    = exprAssembler.getCoeff(mmDf_L);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_NL(m_assemblerL->material(),m_defpatches);
    auto mmA_NL    = exprAssembler.getCoeff(mmAf_NL);
    auto mmB_NL    = exprAssembler.getCoeff(mmBf_NL);
    auto mmC_NL    = exprAssembler.getCoeff(mmCf_NL);
    auto mmD_NL    = exprAssembler.getCoeff(mmDf_NL);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerL->material(),m_defpatches);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);
    // auto rho  = exprAssembler.getCoeff(rhof);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    auto usol   = exprAssembler.getCoeff(primal);
    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);

    // // Matrix K_L
    // auto N      = S0.tr();
    // auto Em_der = flat( jac(Gori).tr() * (jac(zsolH) - jac(zsolL)) ) ;

    // auto M      = S1.tr(); // output is a column
    // auto Ef_der = -(deriv2(zsolH, sn(Gori).normalized().tr()) - deriv2(zsolL, sn(Gori).normalized().tr()) + deriv2(Gori, var1(zsolH, Gori)) - deriv2(Gori, var1(zsolL, Gori))) * reshape(m2, 3, 3);

    // auto N_der  = Em_der * reshape(mmA,3,3) + Ef_der * reshape(mmB,3,3);
    // auto M_der  = Em_der * reshape(mmC,3,3) + Ef_der * reshape(mmD,3,3);

    // auto Fint   = ( N_der * Em_der.tr() + M_der * Ef_der.tr() );

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
    auto Em_der_L   = flat(  jac(Gori).tr() * ( jac(usol )               ) ) ;
    auto Em_derd_L  = flat(  jac(Gori).tr() * ( jac(zsolH) - jac(zsolL) ) ) ;

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

    auto Fint_L     = ( N_derd_L * Em_derd_L.tr() + M_derd_L * Ef_derd_L.tr() );

    auto N_NL       = S0.tr();
    auto Em_der_NL  = flat(  jac(Gdef).tr() * ( jac(usol )               ) ) ;
    auto Em_derd_NL = flat(  jac(Gdef).tr() * ( jac(zsolH) - jac(zsolL) ) ) ;

    auto Em_der2    = flat( jac(usol).tr() * ( jac(usol)                ) ) * N_NL.tr();
    auto Emd_der2   = flat( jac(usol).tr() * ( jac(zsolH) - jac(zsolL) ) ) * N_NL.tr();

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

    auto Fint_NL    = ( N_derd_NL * Em_der_NL.tr() + M_derd_NL * Ef_der_NL.tr() + Em_der2 + Ef_der2 );

    auto K_L        = ( N_der_L   * Em_der_L.tr()  + M_der_L   * Ef_der_L.tr()                          );
    auto Kd_L       = ( N_derd_L  * Em_der_L.tr()  + M_derd_L  * Ef_der_L.tr()                          );
    auto K_NL       = ( N_der_NL  * Em_der_NL.tr() + M_der_NL  * Ef_der_NL.tr() + Em_der2  + Ef_der2    );
    auto Kd_NL      = ( N_derd_NL * Em_der_NL.tr() + M_derd_NL * Ef_der_NL.tr() + Emd_der2 + Efd_der2   );

    auto A          = Kd_L;
    auto Bdiff      = Kd_NL - Kd_L;
    auto Bprimal    = K_NL  - K_L ;

    gsExprEvaluator<T> ev(exprAssembler);

    std::vector<T> result;

    T integral = 0.0;
    integral += ev.integral( A  * meas(Gori) );
    integral += -evPrimalL * ev.integral( Bdiff * meas(Gori));
    integral += (evDualH - evDualL) * (ev.integral( Bprimal * meas(Gori)) - 1.0);

    if (_elWise == 0)
        result.push_back(integral);
    else if (_elWise == 1)
        result.push_back(integral);
    else if (_elWise == 2)
        result.push_back(integral);
    else
        GISMO_ERROR("Unknown");

    return result;
}


/**
 * @brief      Computes the eigenvalue error estimate for Modal Analysis
 *
 * @param[in]  evPrimalL  The primal eigenvalue on the low order basis
 * @param[in]  evDualL    The dual eigenvalue on the low order basis
 * @param[in]  evDualH    The dual eigenvalue on the high order basis
 * @param[in]  dualL      The dual modeshape on the low order basis
 * @param[in]  dualH      The dual modeshape on the high order basis
 * @param[in]  primal     The physical mode shape
 *
 * @tparam     d          Dimension
 * @tparam     T          Real type
 * @tparam     bending    Bending flag
 *
 * @return     The error estimate
 */
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig(                 const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal)
{
    std::vector<T> result = computeErrorEig_impl<0>(evPrimalL, evDualL, evDualH, dualL, dualH, primal);
    GISMO_ASSERT(result.size()==1,"Must be scalar!");
    return result[0];
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigElements(    const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal)
{
    return computeErrorEig_impl<1>(evPrimalL, evDualL, evDualH, dualL, dualH, primal);
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigDofs(    const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal)
{
    return computeErrorEig_impl<2>(evPrimalL, evDualL, evDualH, dualL, dualH, primal);
}

template <short_t d, class T, bool bending>
template <int _elWise>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig_impl(    const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal)
{
    // Everything is evaluated in the lower basis L
    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        m_defpatches.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf(m_assemblerL->material(),m_patches);
    auto mmA    = exprAssembler.getCoeff(mmAf);
    auto mmB    = exprAssembler.getCoeff(mmBf);
    auto mmC    = exprAssembler.getCoeff(mmCf);
    auto mmD    = exprAssembler.getCoeff(mmDf);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerL->material(),m_patches);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);
    auto rho  = exprAssembler.getCoeff(rhof);

    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    auto usol   = exprAssembler.getCoeff(primal);
    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);

    auto N      = S0.tr();
    auto Em_derd= flat( jac(Gori).tr() * (jac(zsolH) - jac(zsolL)) ) ;
    auto Em_der = flat( jac(Gori).tr() * (jac(usol) ) ) ;

    auto M      = S1.tr(); // output is a column
    auto Ef_derd= -( deriv2(zsolH, sn(Gori).normalized().tr() )
                  -  deriv2(zsolL, sn(Gori).normalized().tr() ) + deriv2(Gori, var1(zsolH, Gori) )
                                                                - deriv2(Gori, var1(zsolL, Gori) ) ) * reshape(m2, 3, 3);
    auto Ef_der = -( deriv2(usol, sn(Gori).normalized().tr() ) + deriv2(Gori, var1(usol, Gori)) ) * reshape(m2, 3, 3);

    auto N_der  = Em_derd * reshape(mmA,3,3) + Ef_derd * reshape(mmB,3,3);
    auto M_der  = Em_derd * reshape(mmC,3,3) + Ef_derd * reshape(mmD,3,3);

    auto Fint   = ( N_der * Em_der.tr() + M_der * Ef_der.tr() );

//    space   u   = exprAssembler.trialSpace(0);
    gsVector<T> pt(2);
    pt.setConstant(0.25);

    gsExprEvaluator<T> ev(exprAssembler);

    std::vector<T> result;

    auto A      = Fint;
    auto Bdiff  = rho.val() * usol.tr() * (zsolH - zsolL);
    auto Bprimal= rho.val() * usol.tr() * usol;
    gsConstantFunction<T> onefun(1,2);
    auto one    =exprAssembler.getCoeff(onefun);

    T integral = 0;

    if (_elWise == 0)
    {
        integral = ev.integral( A  * meas(Gori) - evPrimalL * Bdiff * meas(Gori) + (evDualH - evDualL) * ( Bprimal * meas(Gori)  - one) );
        result.push_back(integral);
    }
    else if (_elWise == 1)
    {
        integral = ev.integralElWise( A  * meas(Gori) - evPrimalL * Bdiff * meas(Gori) + (evDualH - evDualL) * ( Bprimal * meas(Gori)  - one) );
        result = ev.elementwise();
    }
    else if (_elWise == 2)
        result.push_back(integral);
    else
        GISMO_ERROR("Unknown");


    space   u2= exprAssembler.getSpace(m_basisL,d);
    gsBoundaryConditions<> bc;
    u2.setup(bc, dirichlet::interpolation, 0);

    gsDebugVar(ev.integral( Gori.tr() * Gori));
    gsDebugVar(ev.eval( Gori.tr() * Gori, pt));
    gsDebugVar(ev.eval( u2.rowSum(), pt));

    gsDebugVar((u2).rows());
    gsDebugVar((u2).cols());
    gsDebugVar((u2.rowSum()).rows());
    gsDebugVar((u2.rowSum()).cols());
    gsDebugVar(( Gori.tr() * Gori).rows());
    gsDebugVar(( Gori.tr() * Gori).cols());


    gsDebugVar(ev.eval( Gori.tr() * Gori* u2.rowSum(), pt));

    exprAssembler.assemble(Gori.tr() * Gori * u2.rowSum());
    gsDebugVar(exprAssembler.rhs().sum());





    return result;
}

//============================================================
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeGoal(const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::ConTransform>  Tconf(m_assemblerL->material(),m_defpatches,Z);
    auto Tcon = exprAssembler.getCoeff(Tconf);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Force tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerL->material(),m_defpatches);
    auto N  = exprAssembler.getCoeff(Nf);
    auto M  = exprAssembler.getCoeff(Mf);

    ///////// TEMPORARY
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0if(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1if(m_assemblerL->material(),m_defpatches);
    auto S0i = exprAssembler.getCoeff(S0if);
    auto S1i = exprAssembler.getCoeff(S1if);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmAfi(m_assemblerL->material(),m_defpatches);
    auto mmAi = exprAssembler.getCoeff(mmAfi);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(m_assemblerL->material(),m_defpatches);
    auto mmDi = exprAssembler.getCoeff(mmDfi);
    ///////// TEMPORARY

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_assemblerL->material(),m_defpatches,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Helper matrix for flexural components
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);

    // Principal stretch
    auto lambda     = reshape(Tcon,3,3) * Cm.tr();

    // Membrane strain
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    // Principal membrane strain
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tcon,3,3).tr();

    // Flexural strain
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;


    // Membrane force (its first variation)
    // auto N_der    = Em_der * reshape(mmA,3,3) + Ef_der * reshape(mmB,3,3);

    gsExprEvaluator<T> ev(exprAssembler);
    if (m_goalFunction == GoalFunction::Displacement)
    {
        if (m_component==9)
        {
           auto expr = (Gdef - Gori).tr() * (Gdef - Gori) * meas(Gori);
           return ev.integral( expr  );
        }
        else
        {
            auto expr = (Gdef - Gori).tr() * gismo::expr::uv(m_component,3)*meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        /*
            NOTE:
            - Cm * Transformation does only give the in-plane stretches, since C is constructed with the in-plane jacobians (Gdef) only. In fact, Cm should have a third entry with the Jacobian Determinant J0, but this one is hard to compute the variations for.
         */
        if (m_component==9)
        {
            // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
            // auto lambda = exprAssembler.getCoeff(lambdaf);
            // auto expr1 = gismo::expr::pow(lambda.tr() * lambda,2) * meas(Gori);
            auto expr = lambda.tr() * lambda * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = lambda.tr() * gismo::expr::uv(m_component,3)*meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            // auto expr = Em * Em.tr() * meas(Gori);
            auto expr = Em.sqNorm() * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = Em * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembranePStrain)
    {
        if (m_component==9)
        {
            auto expr = Emp * Emp.tr() * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = Emp * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = S0.tr() * S0 * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = S0.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembranePStress)
    {
        if (m_component==9)
        {
            auto expr = P0.tr() * P0 * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = P0.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            auto expr = N.tr() * N * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = N.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            auto expr = Ef * Ef.tr() * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            gsVector<T> pt(2);
            pt.setConstant(0.25);
            gsDebugVar(ev.eval(Em,pt));
            gsDebugVar(ev.eval(S0i,pt));
            gsDebugVar(ev.eval(reshape(mmAi,3,3)*Em.tr(),pt));

            gsDebugVar(ev.eval(Ef,pt));
            gsDebugVar(ev.eval(S1i,pt));
            gsDebugVar(ev.eval(reshape(mmDi,3,3)*Ef.tr(),pt));

            auto expr = Ef * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = S1.tr() * S1 * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = S1.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            auto expr = M.tr() * M * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            auto expr = M.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integral( expr  );
        }
    }
    else
        GISMO_ERROR("Goal function not known");
}
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerL->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerL->options());

    m_defpatches = deformed;

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(m_defpatches);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::ConTransform>  Tconf(m_assemblerL->material(),m_defpatches,Z);
    auto Tcon = exprAssembler.getCoeff(Tconf);
    // Material tensors

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerL->material(),m_defpatches,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(m_assemblerL->material(),m_defpatches,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Force tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerL->material(),m_defpatches);
    auto N  = exprAssembler.getCoeff(Nf);
    auto M  = exprAssembler.getCoeff(Mf);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStressN> P0f(m_assemblerL->material(),m_defpatches,Z);
    auto P0 = exprAssembler.getCoeff(P0f);

    // // Helper matrix for flexural components
    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm      = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);

    // Principal stretch
    auto lambda  = reshape(Tcon,3,3).tr() * Cm.tr();

    // Membrane strain
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    // Principal membrane strain
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tcon,3,3).tr();

    // Flexural strain
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;

    // Membrane force (its first variation)
    // auto N_der    = Em_der * reshape(mmA,3,3) + Ef_der * reshape(mmB,3,3);

    gsExprEvaluator<T> ev(exprAssembler);
    T result = 0;
    gsMatrix<T> tmp;
    if (m_goalFunction == GoalFunction::Displacement)
    {
        if (m_component==9)
        {
           auto expr = (Gdef - Gori).tr() * (Gdef - Gori) * meas(Gori);
           for (index_t k = 0; k!=points.cols(); k++)
           {
               tmp = ev.eval( expr, points.col(k));
               GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
               result += tmp(0,0);
           }
           return result;
        }
        else
        {
            auto expr = (Gdef - Gori).tr() * gismo::expr::uv(m_component,3)*meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        /*
            NOTE:
            - Cm * Transformation does only give the in-plane stretches, since C is constructed with the in-plane jacobians (Gdef) only. In fact, Cm should have a third entry with the Jacobian Determinant J0, but this one is hard to compute the variations for.
         */
        if (m_component==9)
        {
            // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerL->material(),m_defpatches,Z);
            // auto lambda = exprAssembler.getCoeff(lambdaf);
            // auto expr1 = gismo::expr::pow(lambda.tr() * lambda,2) * meas(Gori);
            auto expr = lambda.tr() * lambda * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = (Gdef - Gori).tr() * gismo::expr::uv(m_component,3)*meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            auto expr = Em * Em.tr() * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = Em * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::MembranePStrain)
    {
        if (m_component==9)
        {
            auto expr = Emp * Emp.tr() * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = Emp * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = S0.tr() * S0 * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = S0.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::MembranePStress)
    {
        if (m_component==9)
        {
            auto expr = P0.tr() * P0 * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = P0.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            auto expr = N.tr() * N * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = N.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            auto expr = Ef * Ef.tr() * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = Ef * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = S1.tr() * S1 * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = S1.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            auto expr = M.tr() * M * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
        else
        {
            auto expr = M.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            for (index_t k = 0; k!=points.cols(); k++)
            {
                tmp = ev.eval( expr, points.col(k));
                GISMO_ASSERT(tmp.cols()==tmp.rows() && tmp.cols()==1,"Goal function must be scalar!");
                result += tmp(0,0);
            }
            return result;
        }
    }
    else
        GISMO_ERROR("Goal function not known");

}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::matrixNorm(const gsMultiPatch<T> &primal, const gsMultiPatch<T> &other) const
{
    GISMO_ASSERT(m_goalFunction==GoalFunction::Modal,"Only available for modal goal function");

    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerH->options());

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches); // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerH->material(),m_defpatches);
    auto rho = exprAssembler.getCoeff(rhof);

    auto usol = exprAssembler.getCoeff(primal);
    auto zsol = exprAssembler.getCoeff(other);

    gsExprEvaluator<T> ev(exprAssembler);
    return ev.integral(rho.val() * usol.tr() * zsol * meas(Gori));
}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::matrixNorm(const gsMultiPatch<T> &primal, const gsMultiPatch<T> &other, const gsMultiPatch<T> &deformed) const
{
    GISMO_ASSERT(m_goalFunction==GoalFunction::Buckling,"Only available for buckling goal function");
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    exprAssembler.setOptions(m_assemblerH->options());

    gsMultiPatch<T> defpatches = deformed;

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_L(m_assemblerL->material(),m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_L(m_assemblerL->material(),m_patches);
    auto mmA_L    = exprAssembler.getCoeff(mmAf_L);
    auto mmB_L    = exprAssembler.getCoeff(mmBf_L);
    auto mmC_L    = exprAssembler.getCoeff(mmCf_L);
    auto mmD_L    = exprAssembler.getCoeff(mmDf_L);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_NL(m_assemblerL->material(),m_defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_NL(m_assemblerL->material(),m_defpatches);
    auto mmA_NL    = exprAssembler.getCoeff(mmAf_NL);
    auto mmB_NL    = exprAssembler.getCoeff(mmBf_NL);
    auto mmC_NL    = exprAssembler.getCoeff(mmCf_NL);
    auto mmD_NL    = exprAssembler.getCoeff(mmDf_NL);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->material(),defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->material(),defpatches);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    auto usol = exprAssembler.getCoeff(primal);
    auto zsol = exprAssembler.getCoeff(other);

    auto Em_der_L   = flat(  jac(Gori).tr() * ( jac(usol )               ) ) ;
    auto Ef_der_L   = -( deriv2(usol , sn(Gori).normalized().tr() ) - deriv2(Gori, var1(usol , Gori) ) ) * reshape(m2, 3, 3);

    auto Em_derd_L  = flat(  jac(Gori).tr() * ( jac(zsol )               ) ) ;
    auto Ef_derd_L  = -( deriv2(zsol , sn(Gori).normalized().tr() ) - deriv2(Gori, var1(zsol , Gori) ) ) * reshape(m2, 3, 3);


    auto N_der_L     = Em_derd_L * reshape(mmA_L,3,3) + Ef_derd_L * reshape(mmB_L,3,3);
    auto M_der_L     = Em_derd_L * reshape(mmC_L,3,3) + Ef_derd_L * reshape(mmD_L,3,3);

    auto N_NL       = S0.tr();
    auto Em_der_NL  = flat(  jac(Gdef).tr() * ( jac(usol )               ) ) ;
    auto Ef_der_NL  = -( deriv2(usol , sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(usol , Gdef) ) ) * reshape(m2, 3, 3);

    auto Em_derd_NL = flat(  jac(Gdef).tr() * ( jac(zsol )               ) ) ;
    auto Ef_derd_NL = -( deriv2(zsol , sn(Gdef).normalized().tr() ) - deriv2(Gdef, var1(zsol , Gdef) ) ) * reshape(m2, 3, 3);

    auto Em_der2    = flat( jac(usol).tr() * ( jac(zsol )               ) ) * N_NL.tr();

    auto M_NL       = S1.tr(); // output is a column
    auto Ef_der2    = -(
                          (  ( deriv2(zsol )                 )  * ( var1(usol ,Gdef)                    ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (  ( deriv2(usol )                 )  * ( var1(zsol ,Gdef)                    ).tr()  ).tr()      * reshape(m2, 3, 3)
                        + (    deriv2(Gdef)                     * ( var2(usol ,zsol ,Gdef)              )       ).tr()      * reshape(m2, 3, 3)
                       ) * M_NL.tr();

    auto N_der_NL   = Em_derd_NL  * reshape(mmA_NL,3,3) + Ef_derd_NL  * reshape(mmB_NL,3,3);
    auto M_der_NL   = Em_derd_NL  * reshape(mmC_NL,3,3) + Ef_derd_NL  * reshape(mmD_NL,3,3);

    auto K_L        = ( N_der_L   * Em_der_L.tr()  + M_der_L   * Ef_der_L.tr()                          );
    auto K_NL       = ( N_der_NL  * Em_der_NL.tr() + M_der_NL  * Ef_der_NL.tr() + Em_der2  + Ef_der2    );

    auto Bprimal    = K_NL - K_L ;

    auto expr = ( Bprimal ) * meas(Gori);

    gsExprEvaluator<T> ev(exprAssembler);
    T integral = ev.integral(expr);

    if (m_foundInd)
    {
        auto foundation = exprAssembler.getCoeff(*m_foundFun, Gori);
        GISMO_ASSERT(m_foundFun->targetDim()==3,"Foundation function has dimension "<<m_foundFun->targetDim()<<", but expected 3");

        integral += ev.integral( zsol * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    }
    if (m_pressInd)
    {
        auto pressure = exprAssembler.getCoeff(*m_pressFun, Gori);
        GISMO_ASSERT(m_pressFun->targetDim()==1,"Pressure function has dimension "<<m_pressFun->targetDim()<<", but expected 1");

        integral += ev.integral( pressure.val() * zsol * sn(Gdef).normalized() * meas(Gori) );
    }

    return integral;
}

} // namespace gismo
