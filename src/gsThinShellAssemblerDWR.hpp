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

#include <gsKLShell/src/gsThinShellAssemblerDWR.h>
#include <gsKLShell/src/gsMaterialMatrixNonlinear.h>
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsAssembler/gsExprEvaluator.h>

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
                                                            m_bcs(bconditions)
{
    m_assemblerL = new gsThinShellAssembler<d,T,bending>(patches,basisL,bconditions,surface_force,materialmatrix);
    m_assemblerH = new gsThinShellAssembler<d,T,bending>(patches,basisH,bconditions,surface_force,materialmatrix);

    m_dL = gsVector<T>::Zero(m_assemblerL->numDofs());
    m_dH = gsVector<T>::Zero(m_assemblerH->numDofs());
 }

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assembleMass(gsThinShellAssemblerBase<T> * assembler, gsSparseMatrix<T> & result, bool lumped)
{
    ThinShellAssemblerStatus status =  assembler->assembleMass(lumped);
    result = assembler->matrix();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    return status;
}

// template <short_t d, class T, bool bending>
// gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::_assembleMassLM(gsThinShellAssemblerBase<T> * assembler, bool lumped)
// {
//     gsExprAssembler<T> assembler(1,1);
//     // Elements used for numerical integration

//     assembler.setIntegrationElements(m_assemblerH->getBasis());
    // GISMO_ENSURE(m_assemblerL->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    // exprAssembler.setOptions(m_assemblerL->options().getGroup("ExprAssembler"));


//     // Initialize the geometry maps
//     geometryMap m_ori   = assembler.getMap(m_patches);           // this map is used for integrals

//     // Set the discretization space
//     space spaceL = assembler.getSpace(m_assemblerL->getBasis(), d, 0); // last argument is the space ID
//     space spaceH = assembler.getTestSpace(spaceL , m_assemblerH->getBasis());

//     this->_assembleDirichlet();

//     // Initialize stystem
//     assembler.initSystem();

//     gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,&m_patches);
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

    // gsMaterialMatrixIntegrate<T,MaterialOutput::Density> m_mm(m_materialMat,&m_patches);
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
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::assembleL()
{
    std::pair<gsSparseMatrix<T>,gsVector<T>> result;
    ThinShellAssemblerStatus status = _assemble(m_assemblerL,result);
    m_matrixL = result.first;
    m_pL = result.second;
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::assembleH()
{
    std::pair<gsSparseMatrix<T>,gsVector<T>> result;
    ThinShellAssemblerStatus status = _assemble(m_assemblerH,result);
    m_matrixH = result.first;
    m_pH = result.second;
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assemble(gsThinShellAssemblerBase<T> * assembler, std::pair<gsSparseMatrix<T>,gsVector<T>> & result)
{
    ThinShellAssemblerStatus status = assembler->assemble();
    result = std::make_pair(assembler->matrix(),assembler->rhs());
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assembleMatrix(gsThinShellAssemblerBase<T> * assembler, gsSparseMatrix<T> & result)
{
    ThinShellAssemblerStatus status = assembler->assemble();
    result = assembler->matrix();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assembleMatrix(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed, gsSparseMatrix<T> & result)
{
    ThinShellAssemblerStatus status = assembler->assembleMatrix(deformed);
    result = assembler->matrix();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assemblePrimal(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed, gsVector<T> & result)
{
    ThinShellAssemblerStatus status = assembler->assembleVector(deformed);
    result = assembler->rhs();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assemblePrimal(gsThinShellAssemblerBase<T> * assembler, gsVector<T> & result)
{
    ThinShellAssemblerStatus status = assembler->assemble();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    result = assembler->rhs();
    return status;
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed, gsVector<T> & result)
{
    GISMO_ASSERT(m_component < d || m_component==9,"Component is out of bounds");

    // For the stretches
    // ! Only works with Z=0 for now, since principal stresses and directions are implemented for Z=0, since principal stresses and directions are implemented for Z=0
    gsMatrix<T> Z(1,1);
    Z.setZero();
    ////

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(assembler->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(assembler->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

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
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>  Tstretchf(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressTransform>  Tpstressf(assembler->materials(),&deformed,Z);
    auto Tstretch = exprAssembler.getCoeff(Tstretchf);
    auto Tpstress = exprAssembler.getCoeff(Tpstressf);

    // Material tensors at Z
    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::MatrixD> mmDf(assembler->materials(),&deformed,Z);
    auto mmA = exprAssembler.getCoeff(mmAf);
    auto mmD = exprAssembler.getCoeff(mmDf);

    // Material tensors integrated
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(assembler->materials(),&deformed);
    auto mmAi = exprAssembler.getCoeff(mmAfi);
    auto mmBi = exprAssembler.getCoeff(mmBfi);
    auto mmCi = exprAssembler.getCoeff(mmCfi);
    auto mmDi = exprAssembler.getCoeff(mmDfi);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(assembler->materials(),&deformed,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStress> P0f(assembler->materials(),&deformed,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Force and moment tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(assembler->materials(),&deformed);
    auto N = exprAssembler.getCoeff(Nf);
    auto M = exprAssembler.getCoeff(Mf);

    // Helper matrix for flexural components
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient and its first variation
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<T> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m12,3,3);

    // Principal stretch and its first variation
    auto lambda     = Cm     * reshape(Tstretch,3,3).tr();
    auto lambda_der = Cm_der * reshape(Tstretch,3,3).tr();

    // Membrane strain and its first variation
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    // Principal membrane strain and its first variation
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();
    auto Emp_der= (Em_der * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();

    // Flexural strain and its first variation
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;
    auto Ef_der = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);

    // Membrane stress (its first variation)
    auto Sm_der     = Em_der * reshape(mmA,3,3);

    // Principal membrane stress (its first variation)
    auto Smp_der     = Sm_der * reshape(Tpstress,3,3).tr();

    // Flexural stress (its first variation)
    auto Sf_der     = Ef_der * reshape(mmD,3,3);

    // Membrane force (its first variation)
    auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
    auto M_der    = Em_der * reshape(mmCi,3,3) + Ef_der * reshape(mmDi,3,3);

    if (m_goalFunction == GoalFunction::Displacement)
    {
        if (m_component==9)
        {
            auto expr = 2 * space * usol * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = space * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        if (m_component==9)
        {
            auto expr = 2 * lambda_der * lambda.tr() * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = lambda_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::PStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Emp_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::PStress)
    {
        // to convert P0 from a vector of length 2 to 3, this matrix is [[1,0],[0,1],[0,0]]
        gsMatrix<T> convMat(3,2);
        convMat.setZero();
        convMat(0,0) = convMat(1,1) = 1;
        gsConstantFunction<T> convFun(convMat.reshape(6,1),2);
        auto conv = exprAssembler.getCoeff(convFun);
        if (m_component==9)
        {
            auto expr = 2 * Smp_der * reshape(conv,3,2) * P0 * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            GISMO_ASSERT(m_component < 2,"Can only select principle stress component 0 or 1, but "<<m_component<<" selected.");
            auto expr = Smp_der * reshape(conv,3,2) * gismo::expr::uv(m_component,2) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Em_der * Em.tr() * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Em_der *  gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sm_der * S0 * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Sm_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            auto expr = 2 * N_der * N * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = N_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Ef_der * Ef.tr() * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Ef_der *  gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sf_der * S1 * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Sf_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            auto expr = 2 * M_der * M * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = M_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assemble(expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else
        GISMO_ERROR("Goal function not known");
}

template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assembleDual(const bContainer & bnds, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed, gsVector<T> & result)
{
    GISMO_ASSERT(m_component < d || m_component==9,"Component is out of bounds");

    // For the stretches
    // ! Only works with Z=0 for now, since principal stresses and directions are implemented for Z=0, since principal stresses and directions are implemented for Z=0
    gsMatrix<T> Z(1,1);
    Z.setZero();
    ////

    gsExprAssembler<T> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(assembler->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(assembler->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);
    if (bnds.size()==0)
    {
        result = exprAssembler.rhs();
        return ThinShellAssemblerStatus::Success;
    }

    // Space
    space   space   = exprAssembler.trialSpace(0);
    // Homogenize Dirichlet
    space.setup(m_bcs, dirichlet::homogeneous, 0);

    // Solution
    auto usol   = exprAssembler.getCoeff(primal);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>  Tstretchf(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressTransform>  Tpstressf(assembler->materials(),&deformed,Z);
    auto Tstretch = exprAssembler.getCoeff(Tstretchf);
    auto Tpstress = exprAssembler.getCoeff(Tpstressf);

    // Material tensors at Z
    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::MatrixD> mmDf(assembler->materials(),&deformed,Z);
    auto mmA = exprAssembler.getCoeff(mmAf);
    auto mmD = exprAssembler.getCoeff(mmDf);

    // Material tensors integrated
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(assembler->materials(),&deformed);
    auto mmAi = exprAssembler.getCoeff(mmAfi);
    auto mmBi = exprAssembler.getCoeff(mmBfi);
    auto mmCi = exprAssembler.getCoeff(mmCfi);
    auto mmDi = exprAssembler.getCoeff(mmDfi);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(assembler->materials(),&deformed,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStress> P0f(assembler->materials(),&deformed,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Force and moment tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(assembler->materials(),&deformed);
    auto N = exprAssembler.getCoeff(Nf);
    auto M = exprAssembler.getCoeff(Mf);

    // Helper matrix for flexural components
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient and its first variation
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<T> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m12,3,3);

    // Principal stretch and its first variation
    auto lambda     = Cm     * reshape(Tstretch,3,3).tr();
    auto lambda_der = Cm_der * reshape(Tstretch,3,3).tr();

    // Membrane strain and its first variation
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    // Principal membrane strain and its first variation
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();
    auto Emp_der= (Em_der * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();

    // Flexural strain and its first variation
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;
    auto Ef_der = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);

    // Membrane stress (its first variation)
    auto Sm_der     = Em_der * reshape(mmA,3,3);

    // Principal membrane stress (its first variation)
    auto Smp_der     = Sm_der * reshape(Tpstress,3,3).tr();

    // Flexural stress (its first variation)
    auto Sf_der     = Ef_der * reshape(mmD,3,3);

    // Membrane force (its first variation)
    auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
    auto M_der    = Em_der * reshape(mmCi,3,3) + Ef_der * reshape(mmDi,3,3);

    if (m_goalFunction == GoalFunction::Displacement)
    {
        if (m_component==9)
        {
            auto expr = 2 * space * usol * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = space * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        if (m_component==9)
        {
            auto expr = 2 * lambda_der * lambda.tr() * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = lambda_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::PStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Emp_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::PStress)
    {
        // to convert P0 from a vector of length 2 to 3, this matrix is [[1,0],[0,1],[0,0]]
        gsMatrix<T> convMat(3,2);
        convMat.setZero();
        convMat(0,0) = convMat(1,1) = 1;
        gsConstantFunction<T> convFun(convMat.reshape(6,1),2);
        auto conv = exprAssembler.getCoeff(convFun);
        if (m_component==9)
        {
            auto expr = 2 * Smp_der * reshape(conv,3,2) * P0 * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            GISMO_ASSERT(m_component < 2,"Can only select principle stress component 0 or 1, but "<<m_component<<" selected.");
            auto expr = Smp_der * reshape(conv,3,2) * gismo::expr::uv(m_component,2) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            auto expr = 2 * Em_der * Em.tr() * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Em_der *  gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sm_der * S0 * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Sm_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            auto expr = 2 * N_der * N * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = N.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            gsFunctionExpr<T> mult110t("1","0","0","0","1","0","0","0","0",2);
            auto m110 = exprAssembler.getCoeff(mult2t);
            auto expr = 2*Ef_der * reshape(m110,3,3) *  Ef.tr() * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Ef_der *  gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Sf_der * S1 * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = Sf_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            auto expr = 2 * M_der * M * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
        else
        {
            auto expr = M_der * gismo::expr::uv(m_component,3) * meas(Gori);
            try
            {
                exprAssembler.assembleBdr(bnds,expr);
                result = exprAssembler.rhs();
                return ThinShellAssemblerStatus::Success;
            }
            catch (...)
            {
                exprAssembler.cleanUp();
                return ThinShellAssemblerStatus::AssemblyError;
            }
        }
    }
    else
        GISMO_ERROR("Goal function not known");
}


template <short_t d, class T, bool bending>
ThinShellAssemblerStatus gsThinShellAssemblerDWR<d, T, bending>::_assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed, gsVector<T> & result)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
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
    GISMO_ENSURE(assembler->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(assembler->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);
    if (points.cols()==0)
    {
        result = exprAssembler.rhs();
        return ThinShellAssemblerStatus::Success;
    }

    // Space
    space   space   = exprAssembler.trialSpace(0);
    // Homogenize Dirichlet
    space.setup(m_bcs, dirichlet::homogeneous, 0);

    // Solution
    auto usol   = exprAssembler.getCoeff(primal);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>  Tstretchf(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::PStressTransform>  Tpstressf(assembler->materials(),&deformed,Z);
    auto Tstretch = exprAssembler.getCoeff(Tstretchf);
    auto Tpstress = exprAssembler.getCoeff(Tpstressf);

    // Material tensors at Z
    gsMaterialMatrixEval<T,MaterialOutput::MatrixA> mmAf(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::MatrixD> mmDf(assembler->materials(),&deformed,Z);
    auto mmA = exprAssembler.getCoeff(mmAf);
    auto mmD = exprAssembler.getCoeff(mmDf);

    // Material tensors integrated
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCfi(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(assembler->materials(),&deformed);
    auto mmAi = exprAssembler.getCoeff(mmAfi);
    auto mmBi = exprAssembler.getCoeff(mmBfi);
    auto mmCi = exprAssembler.getCoeff(mmCfi);
    auto mmDi = exprAssembler.getCoeff(mmDfi);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(assembler->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(assembler->materials(),&deformed,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStress> P0f(assembler->materials(),&deformed,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Force and moment tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(assembler->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(assembler->materials(),&deformed);
    auto N = exprAssembler.getCoeff(Nf);
    auto M = exprAssembler.getCoeff(Mf);

    // Helper matrix for flexural components
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient and its first variation
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<T> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
    auto Cm_der     = 2* flat( jac(Gdef).tr() * jac(space) ) * reshape(m12,3,3);

    // Principal stretch and its first variation
    auto lambda     = Cm     * reshape(Tstretch,3,3).tr();
    auto lambda_der = Cm_der * reshape(Tstretch,3,3).tr();

    // Membrane strain and its first variation
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    // Principal membrane strain and its first variation
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();
    auto Emp_der= (Em_der * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();

    // Flexural strain and its first variation
    auto Ef     = ( deriv2(Gori,sn(Gori).normalized().tr()) - deriv2(Gdef,sn(Gdef).normalized().tr()) ) * reshape(m2,3,3) ;
    auto Ef_der = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);

    // Membrane stress (its first variation)
    auto Sm_der     = Em_der * reshape(mmA,3,3);

    // Principal membrane stress (its first variation)
    auto Smp_der     = Sm_der * reshape(Tpstress,3,3).tr();

    // Flexural stress (its first variation)
    auto Sf_der     = Ef_der * reshape(mmD,3,3);

    // Membrane force (its first variation)
    auto N_der    = Em_der * reshape(mmAi,3,3) + Ef_der * reshape(mmBi,3,3);
    auto M_der    = Em_der * reshape(mmCi,3,3) + Ef_der * reshape(mmDi,3,3);


    gsExprEvaluator<T> ev(exprAssembler);
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
        }
    }
    else if (m_goalFunction == GoalFunction::Stretch)
    {
        if (m_component==9)
        {
            auto expr = 2* lambda_der * lambda.tr() * meas(Gori);//<<<--------NOTE: Square missing?
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
        }
    }
    else if (m_goalFunction == GoalFunction::PStrain)
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
        }
    }
    else if (m_goalFunction == GoalFunction::PStress)
    {
        if (m_component==9)
        {
            auto expr = 2 * Smp_der * P0 * meas(Gori);
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
        }
        else
        {
            auto expr = Smp_der * gismo::expr::uv(m_component,3) * meas(Gori);
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
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
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
        }
        else
        {
            // remove N_der from this scipe
            auto expr = N.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
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
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
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
        }
        else
        {
            // remove N_der from this scipe
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
        }
    }
    else
        GISMO_ERROR("Goal function not known");

    return ThinShellAssemblerStatus::Success;
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
void gsThinShellAssemblerDWR<d, T, bending>::updateMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    m_assemblerL->updateMultiPatch(solVector,result);
}

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::updateMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    m_assemblerH->updateMultiPatch(solVector,result);
}

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

template <short_t d, class T, bool bending>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::constructSolutionVectorL(const gsMultiPatch<T> & deformed)
{
    return m_assemblerL->constructSolutionVector(deformed);
}

template <short_t d, class T, bool bending>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending>::constructSolutionVectorH(const gsMultiPatch<T> & deformed)
{
    return m_assemblerH->constructSolutionVector(deformed);
}

template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::projectL2_L(const gsFunction<T> &fun)
{
    return m_assemblerL->projectL2(fun);
}

template <short_t d, class T, bool bending>
gsMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::projectL2_H(const gsFunction<T> &fun)
{
    return m_assemblerH->projectL2(fun);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeError_impl<0>(dualL,dualH,withLoads,filename,np,parametric,mesh);
    return m_error;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeError_impl<1>(dualL,dualH,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeError_impl<2>(dualL,dualH,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
template <index_t _elWise>
void gsThinShellAssemblerDWR<d, T, bending>::computeError_impl( const gsMultiPatch<T> & dualL,
                                                                const gsMultiPatch<T> & dualH,
                                                                bool withLoads,
                                                                std::string filename,
                                                                unsigned np,
                                                                bool parametric,
                                                                bool mesh)
{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);

    auto F      = exprAssembler.getCoeff(*m_forceFun, Gori);
    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);

    auto g_N = exprAssembler.getBdrFunction(Gori);

    auto expr = (zsolH-zsolL).tr() * F;
    auto bexpr= (zsolH-zsolL).tr() * g_N;

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    if (_elWise == 0)
    {
        integral = ev.integral(expr * meas(Gori)); // this one before otherwise it gives an error (?)
        bintegral = ev.integralBdrBc( m_bcs.get("Neumann"), bexpr * tv(Gori).norm());
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,integral);

        m_error = integral+bintegral;
    }
    else if (_elWise == 1) // element-wise
    {
        m_error = ev.integralElWise(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
        m_errors = ev.elementwise();
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToElWiseError(dualL,dualH,m_errors);
    }
    else if (_elWise == 2) // function-wise
    {
        integral = 0;
    }
    else
        GISMO_ERROR("Unknown");

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeError_impl<d,bending,0>(dualL,dualH,deformed,withLoads,filename,np,parametric,mesh);
    return m_error;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeError_impl<d,bending,1>(dualL,dualH,deformed,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeError_impl<d,bending,2>(dualL,dualH,deformed,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
template<short_t _d, bool _bending, index_t _elWise>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssemblerDWR<d, T, bending>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                            // bool withLoads,
                                                            std::string filename, unsigned np, bool parametric, bool mesh)

{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&deformed);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
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

    auto g_N = exprAssembler.getBdrFunction(Gori);

    auto expr = ( Fext - Fint );
    auto bexpr= (zsolH-zsolL).tr() * g_N;

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral = 0;
    if (_elWise == 0)
    {
        if (withLoads)
            bintegral = ev.integralBdrBc( m_bcs.get("Neumann"), bexpr * tv(Gori).norm()); // goes wrong with sizes?
        integral = ev.integral(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,integral);

        m_error = integral+bintegral;
    }
    else if (_elWise == 1) // element-wise
    {
        m_error = ev.integralElWise(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
        m_errors = ev.elementwise();
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToElWiseError(dualL,dualH,m_errors);
    }
    else if (_elWise == 2) // function-wise
    {
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

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

template <short_t d, class T, bool bending>
template<short_t _d, bool _bending, index_t _elWise>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssemblerDWR<d, T, bending>::computeError_impl(  const gsMultiPatch<T> & dualL,
                                                            const gsMultiPatch<T> & dualH,
                                                            const gsMultiPatch<T> & deformed,
                                                            bool withLoads,
                                                            std::string filename,
                                                            unsigned np,
                                                            bool parametric,
                                                            bool mesh)
{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed);
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

    auto g_N = exprAssembler.getBdrFunction(Gori);

    auto expr = ( Fext - Fint );
    auto bexpr= (zsolH-zsolL).tr() * g_N;

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    if (_elWise == 0)
    {
        bintegral = ev.integralBdrBc( m_bcs.get("Neumann"), bexpr * tv(Gori).norm());
        integral = ev.integral(expr * meas(Gori));
        m_error = bintegral+integral;
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
    }
    else if (_elWise == 1) // element-wise
    {
        m_error = ev.integralElWise(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
        m_errors = ev.elementwise();
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToElWiseError(dualL,dualH,m_errors);
    }
    else if (_elWise == 2) // function-wise
    {
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

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeSquaredError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeSquaredError_impl<0>(dualL,dualH,withLoads,filename,np,parametric,mesh);
    return m_error;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeSquaredErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeSquaredError_impl<1>(dualL,dualH,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeSquaredErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeSquaredError_impl<2>(dualL,dualH,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
template <index_t _elWise>
void gsThinShellAssemblerDWR<d, T, bending>::computeSquaredError_impl( const gsMultiPatch<T> & dualL,
                                                                const gsMultiPatch<T> & dualH,
                                                                bool withLoads,
                                                                std::string filename,
                                                                unsigned np,
                                                                bool parametric,
                                                                bool mesh)
{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);

    auto F      = exprAssembler.getCoeff(*m_forceFun, Gori);
    auto zsolL  = exprAssembler.getCoeff(dualL);
    auto zsolH  = exprAssembler.getCoeff(dualH);

    auto g_N = exprAssembler.getBdrFunction(Gori);

    auto expr = ((zsolH-zsolL).tr() * F).sqNorm();
    auto bexpr= (zsolH-zsolL).tr() * g_N;

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    if (_elWise == 0)
    {
        integral = ev.integral(expr * meas(Gori)); // this one before otherwise it gives an error (?)
        bintegral = ev.integralBdrBc( m_bcs.get("Neumann"), bexpr * tv(Gori).norm());
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,integral);

        m_error = integral+bintegral;
    }
    else if (_elWise == 1) // element-wise
    {
        m_error = ev.integralElWise(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
        m_errors = ev.elementwise();
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToElWiseError(dualL,dualH,m_errors);
    }
    else if (_elWise == 2) // function-wise
    {
        integral = 0;
    }
    else
        GISMO_ERROR("Unknown");

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeSquaredError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeSquaredError_impl<d,bending,0>(dualL,dualH,deformed,withLoads,filename,np,parametric,mesh);
    return m_error;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeSquaredErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeSquaredError_impl<d,bending,1>(dualL,dualH,deformed,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeSquaredErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                        std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeSquaredError_impl<d,bending,2>(dualL,dualH,deformed,withLoads,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
template<short_t _d, bool _bending, index_t _elWise>
typename std::enable_if<_d==3 && _bending, void>::type
gsThinShellAssemblerDWR<d, T, bending>::computeSquaredError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                                                            // bool withLoads,
                                                            std::string filename, unsigned np, bool parametric, bool mesh)

{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&deformed);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
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

    auto g_N = exprAssembler.getBdrFunction(Gori);

    auto expr = ( Fext - Fint ).sqNorm();
    auto bexpr= (zsolH-zsolL).tr() * g_N;

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    if (_elWise == 0)
    {
        bintegral = ev.integralBdrBc( m_bcs.get("Neumann"), bexpr * tv(Gori).norm());
        integral = ev.integral(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,integral);

        m_error = integral+bintegral;
    }
    else if (_elWise == 1) // element-wise
    {
        m_error = ev.integralElWise(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
        m_errors = ev.elementwise();
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToElWiseError(dualL,dualH,m_errors);
    }
    else if (_elWise == 2) // function-wise
    {
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

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

template <short_t d, class T, bool bending>
template<short_t _d, bool _bending, index_t _elWise>
typename std::enable_if<!(_d==3 && _bending), void>::type
gsThinShellAssemblerDWR<d, T, bending>::computeSquaredError_impl(  const gsMultiPatch<T> & dualL,
                                                            const gsMultiPatch<T> & dualH,
                                                            const gsMultiPatch<T> & deformed,
                                                            bool withLoads,
                                                            std::string filename,
                                                            unsigned np,
                                                            bool parametric,
                                                            bool mesh)
{
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed);
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

    auto g_N = exprAssembler.getBdrFunction(Gori);

    auto expr = ( Fext - Fint ).sqNorm();
    auto bexpr= (zsolH-zsolL).tr() * g_N;

    gsExprEvaluator<T> ev(exprAssembler);

    T integral, bintegral;
    if (_elWise == 0)
    {
        bintegral = ev.integralBdrBc( m_bcs.get("Neumann"), bexpr * tv(Gori).norm());
        integral = ev.integral(expr * meas(Gori));
        m_error = bintegral+integral;
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
    }
    else if (_elWise == 1) // element-wise
    {
        m_error = ev.integralElWise(expr * meas(Gori));
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToError(dualL,dualH,m_error);
        m_errors = ev.elementwise();
        if (m_pLoads.numLoads()!=0 && withLoads)
            _applyLoadsToElWiseError(dualL,dualH,m_errors);
    }
    else if (_elWise == 2) // function-wise
    {
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

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
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
T gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig(  const T evPrimalL,
                                                            const T evDualL,
                                                            const T evDualH,
                                                            const gsMultiPatch<T> & dualL,
                                                            const gsMultiPatch<T> & dualH,
                                                            const gsMultiPatch<T> & primal,
                                                            const gsMultiPatch<T> & deformed,
                                                            std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeErrorEig_impl<0>(evPrimalL, evDualL, evDualH, dualL, dualH, primal, deformed,filename,np,parametric,mesh);
    return m_error;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigElements( const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                const gsMultiPatch<T> & deformed,
                                                                                std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeErrorEig_impl<1>(evPrimalL, evDualL, evDualH, dualL, dualH, primal, deformed,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigDofs( const T evPrimalL,
                                                                            const T evDualL,
                                                                            const T evDualH,
                                                                            const gsMultiPatch<T> & dualL,
                                                                            const gsMultiPatch<T> & dualH,
                                                                            const gsMultiPatch<T> & primal,
                                                                            const gsMultiPatch<T> & deformed,
                                                                            std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeErrorEig_impl<2>(evPrimalL, evDualL, evDualH, dualL, dualH, primal, deformed,filename,np,parametric,mesh);
    return m_errors;
}


template <short_t d, class T, bool bending>
template <index_t _elWise>
void gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig_impl(    const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                const gsMultiPatch<T> & deformed,
                                                                                std::string filename,
                                                                                unsigned np,
                                                                                bool parametric,
                                                                                bool mesh)
{
    // Everything is evaluated in the lower basis L
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_L(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_L(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_L(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_L(m_assemblerH->materials(),&m_patches);
    auto mmA_L    = exprAssembler.getCoeff(mmAf_L);
    auto mmB_L    = exprAssembler.getCoeff(mmBf_L);
    auto mmC_L    = exprAssembler.getCoeff(mmCf_L);
    auto mmD_L    = exprAssembler.getCoeff(mmDf_L);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_NL(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_NL(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_NL(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_NL(m_assemblerH->materials(),&deformed);
    auto mmA_NL    = exprAssembler.getCoeff(mmAf_NL);
    auto mmB_NL    = exprAssembler.getCoeff(mmBf_NL);
    auto mmC_NL    = exprAssembler.getCoeff(mmCf_NL);
    auto mmD_NL    = exprAssembler.getCoeff(mmDf_NL);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&deformed);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerH->materials(),deformed);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);
    // auto rho  = exprAssembler.getCoeff(rhof);

    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
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

    auto K_L        = ( N_der_L   * Em_der_L.tr()  + M_der_L   * Ef_der_L.tr()                          );
    auto Kd_L       = ( N_derd_L  * Em_der_L.tr()  + M_derd_L  * Ef_der_L.tr()                          );
    auto K_NL       = ( N_der_NL  * Em_der_NL.tr() + M_der_NL  * Ef_der_NL.tr() + Em_der2  + Ef_der2    );
    auto Kd_NL      = ( N_derd_NL * Em_der_NL.tr() + M_derd_NL * Ef_der_NL.tr() + Emd_der2 + Efd_der2   );

    auto A          = Kd_L;
    auto Bdiff      = Kd_NL - Kd_L;
    auto Bprimal    = K_NL  - K_L ;

    auto expr   =  A - evPrimalL * Bdiff + (evDualH - evDualL) * Bprimal;

    gsExprEvaluator<T> ev(exprAssembler);
    if (_elWise == 0)
    {
        m_error = ev.integral(expr * meas(Gori)) - (evDualH - evDualL);
    }
    else if (_elWise == 1)
    {
        m_error = ev.integralElWise(expr * meas(Gori)) - (evDualH - evDualL);
        m_errors = ev.elementwise();
    }
    else if (_elWise == 2)
    {
        m_error = 0;
    }
    else
        GISMO_ERROR("Unknown");

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

// template <short_t d, class T, bool bending>
// gsMatrix<T> gsThinShellAssemblerDWR<d, T, bending>::evalErrorEig(   gsMatrix<T> & u,
//                                                                     const T evPrimalL,
//                                                                     const T evDualL,
//                                                                     const T evDualH,
//                                                                     const gsMultiPatch<T> & dualL,
//                                                                     const gsMultiPatch<T> & dualH,
//                                                                     const gsMultiPatch<T> & primal,
//                                                                     const gsMultiPatch<T> & deformed)
// {
//     gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
//     exprAssembler.cleanUp();
    // GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    // exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

//     gsExprEvaluator<T> ev(exprAssembler);

//     return ev.eval(expr,u);
// }

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
T gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig(  const T evPrimalL,
                                                            const T evDualL,
                                                            const T evDualH,
                                                            const gsMultiPatch<T> & dualL,
                                                            const gsMultiPatch<T> & dualH,
                                                            const gsMultiPatch<T> & primal,
                                                            std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeErrorEig_impl<0>(evPrimalL, evDualL, evDualH, dualL, dualH, primal,filename,np,parametric,mesh);
    return m_error;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigElements( const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeErrorEig_impl<1>(evPrimalL, evDualL, evDualH, dualL, dualH, primal,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
std::vector<T> gsThinShellAssemblerDWR<d, T, bending>::computeErrorEigDofs(    const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                std::string filename, unsigned np, bool parametric, bool mesh)
{
    computeErrorEig_impl<2>(evPrimalL, evDualL, evDualH, dualL, dualH, primal,filename,np,parametric,mesh);
    return m_errors;
}

template <short_t d, class T, bool bending>
template <index_t _elWise>
void gsThinShellAssemblerDWR<d, T, bending>::computeErrorEig_impl(    const T evPrimalL,
                                                                                const T evDualL,
                                                                                const T evDualH,
                                                                                const gsMultiPatch<T> & dualL,
                                                                                const gsMultiPatch<T> & dualH,
                                                                                const gsMultiPatch<T> & primal,
                                                                                std::string filename,
                                                                                unsigned np,
                                                                                bool parametric,
                                                                                bool mesh)
{
    // Everything is evaluated in the lower basis L
    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    gsMultiPatch<T> deformed = m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf(m_assemblerH->materials(),&m_patches);
    auto mmA    = exprAssembler.getCoeff(mmAf);
    auto mmB    = exprAssembler.getCoeff(mmBf);
    auto mmC    = exprAssembler.getCoeff(mmCf);
    auto mmD    = exprAssembler.getCoeff(mmDf);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerH->materials(),&m_patches);
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

    auto Em_derdw= flat( jac(Gori).tr() * (jac(zsolH) - jac(zsolL)) ) ;


    auto M      = S1.tr(); // output is a column
    auto Ef_derd= -( deriv2(zsolH, sn(Gori).normalized().tr() )
                  -  deriv2(zsolL, sn(Gori).normalized().tr() ) + deriv2(Gori, var1(zsolH, Gori) )
                                                                - deriv2(Gori, var1(zsolL, Gori) ) ) * reshape(m2, 3, 3);
    auto Ef_der = -( deriv2(usol, sn(Gori).normalized().tr() ) + deriv2(Gori, var1(usol, Gori)) ) * reshape(m2, 3, 3);

    auto N_der  = Em_derd * reshape(mmA,3,3) + Ef_derd * reshape(mmB,3,3);
    auto M_der  = Em_derd * reshape(mmC,3,3) + Ef_derd * reshape(mmD,3,3);

    auto Fint   = ( N_der * Em_der.tr() + M_der * Ef_der.tr() );

    gsExprEvaluator<T> ev(exprAssembler);

    auto A      = Fint;
    auto Bdiff  = rho.val() * usol.tr() * (zsolH - zsolL);
    auto Bprimal= rho.val() * usol.tr() * usol;
    auto expr   = A - evPrimalL * Bdiff + (evDualH - evDualL) * Bprimal;

    if (_elWise == 0)
        m_error = ev.integral(expr * meas(Gori)) - (evDualH - evDualL);
    else if (_elWise == 1)
    {
        m_error = ev.integralElWise(expr * meas(Gori) ) - (evDualH - evDualL);
        m_errors = ev.elementwise();
    }
    else if (_elWise == 2)
        m_error = 0;
    else
        GISMO_ERROR("Unknown");

    if (!filename.empty())
    {
        ev.options().setSwitch("plot.elements",mesh);
        ev.options().setInt("plot.npts",np);
        // if (parametric)
        //     ev.writeParaview(expr,filename);
        // else
            ev.writeParaview(expr,Gori,filename);
    }
}

//============================================================
template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeGoal(const gsMultiPatch<T> & deformed)
{
    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>  Tstretchf(m_assemblerH->materials(),&deformed,Z);
    auto Tstretch = exprAssembler.getCoeff(Tstretchf);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&deformed,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Force tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerH->materials(),&deformed);
    auto N  = exprAssembler.getCoeff(Nf);
    auto M  = exprAssembler.getCoeff(Mf);

    ///////// TEMPORARY
    // gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0if(m_assemblerH->materials(),&deformed);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1if(m_assemblerH->materials(),&deformed);
    // auto S0i = exprAssembler.getCoeff(S0if);
    // auto S1i = exprAssembler.getCoeff(S1if);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmAfi(m_assemblerH->materials(),&deformed);
    // auto mmAi = exprAssembler.getCoeff(mmAfi);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(m_assemblerH->materials(),&deformed);
    // auto mmDi = exprAssembler.getCoeff(mmDfi);
    ///////// TEMPORARY

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStress> P0f(m_assemblerH->materials(),&deformed,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Helper matrix for flexural components
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<T> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);

    // Principal stretch
    auto lambda     = reshape(Tstretch,3,3) * Cm.tr();

    // Membrane strain
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    // Principal membrane strain
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();

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
            // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerH->materials(),&deformed,Z);
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
    else if (m_goalFunction == GoalFunction::PStrain)
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
    else if (m_goalFunction == GoalFunction::PStress)
    {
        if (m_component==9)
        {
            auto expr = P0.tr() * P0 * meas(Gori);
            return ev.integral( expr  );
        }
        else
        {
            GISMO_ASSERT(m_component < 2,"Can only select principle stress component 0 or 1, but "<<m_component<<" selected.");
            auto expr = P0.tr() * gismo::expr::uv(m_component,2) * meas(Gori);
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
T gsThinShellAssemblerDWR<d, T, bending>::computeGoal(const bContainer & bnds, const gsMultiPatch<T> & deformed)
{
    if (bnds.size()==0) return 0;

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>  Tstretchf(m_assemblerH->materials(),&deformed,Z);
    auto Tstretch = exprAssembler.getCoeff(Tstretchf);

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&deformed,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Force tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerH->materials(),&deformed);
    auto N  = exprAssembler.getCoeff(Nf);
    auto M  = exprAssembler.getCoeff(Mf);

    ///////// TEMPORARY
    // gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0if(m_assemblerH->materials(),&deformed);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1if(m_assemblerH->materials(),&deformed);
    // auto S0i = exprAssembler.getCoeff(S0if);
    // auto S1i = exprAssembler.getCoeff(S1if);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmAfi(m_assemblerH->materials(),&deformed);
    // auto mmAi = exprAssembler.getCoeff(mmAfi);
    // gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDfi(m_assemblerH->materials(),&deformed);
    // auto mmDi = exprAssembler.getCoeff(mmDfi);
    ///////// TEMPORARY

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStress> P0f(m_assemblerH->materials(),&deformed,Z);
    auto P0  = exprAssembler.getCoeff(P0f);

    // Helper matrix for flexural components
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<T> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);

    // Principal stretch
    auto lambda     = reshape(Tstretch,3,3) * Cm.tr();

    // Membrane strain
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    // Principal membrane strain
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();

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
           return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = (Gdef - Gori).tr() * gismo::expr::uv(m_component,3)*meas(Gori);
            return ev.integralBdr( expr,bnds  );
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
            // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerH->materials(),&deformed,Z);
            // auto lambda = exprAssembler.getCoeff(lambdaf);
            // auto expr1 = gismo::expr::pow(lambda.tr() * lambda,2) * meas(Gori);
            auto expr = lambda.tr() * lambda * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = lambda.tr() * gismo::expr::uv(m_component,3)*meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::PStrain)
    {
        if (m_component==9)
        {
            auto expr = Emp * Emp.tr() * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = Emp * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::PStress)
    {
        if (m_component==9)
        {
            auto expr = P0.tr() * P0 * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            GISMO_ASSERT(m_component < 2,"Can only select principle stress component 0 or 1, but "<<m_component<<" selected.");
            auto expr = P0.tr() * gismo::expr::uv(m_component,2) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStrain)
    {
        if (m_component==9)
        {
            // auto expr = Em * Em.tr() * meas(Gori);
            auto expr = Em.sqNorm() * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = Em * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneStress)
    {
        if (m_component==9)
        {
            auto expr = S0.tr() * S0 * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = S0.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::MembraneForce)
    {
        if (m_component==9)
        {
            auto expr = N.tr() * N * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = N.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStrain)
    {
        if (m_component==9)
        {
            auto expr = Ef * Ef.tr() * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = Ef * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralStress)
    {
        if (m_component==9)
        {
            auto expr = S1.tr() * S1 * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = S1.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else if (m_goalFunction == GoalFunction::FlexuralMoment)
    {
        if (m_component==9)
        {
            auto expr = M.tr() * M * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
        else
        {
            auto expr = M.tr() * gismo::expr::uv(m_component,3) * meas(Gori);
            return ev.integralBdr( expr,bnds  );
        }
    }
    else
        GISMO_ERROR("Goal function not known");
}

template <short_t d, class T, bool bending>
T gsThinShellAssemblerDWR<d, T, bending>::computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    if (points.cols()==0) return 0;

    gsMatrix<T> Z(1,1);
    Z.setZero();

    gsExprAssembler<T> exprAssembler = m_assemblerH->assembler();
    exprAssembler.cleanUp();
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));


    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(deformed);

    // Transformation for stretches
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>  Tstretchf(m_assemblerH->materials(),&deformed,Z);
    auto Tstretch = exprAssembler.getCoeff(Tstretchf);
    // Material tensors

    // Stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&deformed,Z);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&deformed,Z);
    auto S0 = exprAssembler.getCoeff(S0f);
    auto S1 = exprAssembler.getCoeff(S1f);

    // Force tensors
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> Nf(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> Mf(m_assemblerH->materials(),&deformed);
    auto N  = exprAssembler.getCoeff(Nf);
    auto M  = exprAssembler.getCoeff(Mf);

    // Principal stress tensors
    gsMaterialMatrixEval<T,MaterialOutput::PStress> P0f(m_assemblerH->materials(),&deformed,Z);
    auto P0 = exprAssembler.getCoeff(P0f);

    // // Helper matrix for flexural components
    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
    auto m2 = exprAssembler.getCoeff(mult2t);

    // Deformation gradient
        // To compensate for the flat, since flat does [E11,E22,2*E12]
    gsFunctionExpr<T> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    auto m12 = exprAssembler.getCoeff(mult12t);

    auto Cm      = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);

    // Principal stretch
    auto lambda  = reshape(Tstretch,3,3).tr() * Cm.tr();

    // Membrane strain
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;

    // Principal membrane strain
    auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();

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
            // gsMaterialMatrixEval<T,MaterialOutput::Stretch>  lambdaf(m_assemblerH->materials(),&deformed,Z);
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
    else if (m_goalFunction == GoalFunction::PStrain)
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
    else if (m_goalFunction == GoalFunction::PStress)
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
            GISMO_ASSERT(m_component < 2,"Can only select principle stress component 0 or 1, but "<<m_component<<" selected.");
            auto expr = P0.tr() * gismo::expr::uv(m_component,2) * meas(Gori);
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
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches); // this map is used for integrals

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::Density> rhof(m_assemblerH->materials(),&m_patches); // note: is this right?
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
    GISMO_ENSURE(m_assemblerH->options().hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    exprAssembler.setOptions(m_assemblerH->options().getGroup("ExprAssembler"));

    gsMultiPatch<T> defpatches = deformed;

    // Geometries
    geometryMap Gori = exprAssembler.getMap(m_patches);           // this map is used for integrals
    geometryMap Gdef = exprAssembler.getMap(defpatches);

    // Initialize vector
    exprAssembler.initSystem();
    exprAssembler.initVector(1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_L(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_L(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_L(m_assemblerH->materials(),&m_patches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_L(m_assemblerH->materials(),&m_patches);
    auto mmA_L    = exprAssembler.getCoeff(mmAf_L);
    auto mmB_L    = exprAssembler.getCoeff(mmBf_L);
    auto mmC_L    = exprAssembler.getCoeff(mmCf_L);
    auto mmD_L    = exprAssembler.getCoeff(mmDf_L);

    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixA> mmAf_NL(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixB> mmBf_NL(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixC> mmCf_NL(m_assemblerH->materials(),&deformed);
    gsMaterialMatrixIntegrate<T,MaterialOutput::MatrixD> mmDf_NL(m_assemblerH->materials(),&deformed);
    auto mmA_NL    = exprAssembler.getCoeff(mmAf_NL);
    auto mmB_NL    = exprAssembler.getCoeff(mmBf_NL);
    auto mmC_NL    = exprAssembler.getCoeff(mmCf_NL);
    auto mmD_NL    = exprAssembler.getCoeff(mmDf_NL);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_assemblerH->materials(),&defpatches);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_assemblerH->materials(),&defpatches);
    auto S0  = exprAssembler.getCoeff(S0f);
    auto S1  = exprAssembler.getCoeff(S1f);

    gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",2);
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

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::_applyLoadsToElWiseError(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, std::vector<T> & errors) const
{

    for (size_t i = 0; i<m_pLoads.numLoads(); ++i )
    {
        if (m_pLoads[i].value.size()!=d)
            gsWarn<<"Point load has wrong dimension "<<m_pLoads[i].value.size()<<" instead of "<<d<<"\n";
        // Compute actives and values of basis functions on point load location.
        gsMatrix<T> forcePoint;
        if ( m_pLoads[i].parametric )   // in parametric space
            forcePoint = m_pLoads[i].point;
        else                            // in physical space
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);

        std::vector<index_t> elements;

        #pragma omp parallel
        {
            std::vector<index_t> elements_private;
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
            index_t patch_cnt = 0;
#endif

            index_t c = 0;
            for (unsigned patchInd=0; patchInd < m_assemblerH->getBasis().nBases(); ++patchInd)
            {
                // Initialize domain element iterator
                typename gsBasis<T>::domainIter domIt =
                    m_assemblerH->getBasis().piece(patchInd).makeDomainIterator();

#ifdef _OPENMP
                c = patch_cnt + tid;
                patch_cnt += domIt->numElements();// a bit costy
                for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
                for (; domIt->good(); domIt->next() )
#endif
                {
                    bool test = !(((((forcePoint - domIt->lowerCorner()).array() >= 0) &&
                                    (domIt->upperCorner() - forcePoint).array() >= 0)   ).array() == 0).any();
                    if (test)
                        elements_private.push_back(c);
                    #ifdef _OPENMP
                                    c += nt;
                    #else
                                    c++;
                    #endif
                }
#pragma omp critical
                elements.insert(elements.end(), elements_private.begin(), elements_private.end());
            }
        }

        // Compute load value on the point
        gsMatrix<T> Lres, Hres;
        dualL.patch(m_pLoads[i].patch).eval_into(forcePoint,Lres);
        dualH.patch(m_pLoads[i].patch).eval_into(forcePoint,Hres);
        gsMatrix<T> result = (Hres-Lres).transpose() * m_pLoads[i].value;
        // Distribute load over all the elements
        GISMO_ASSERT(result.rows()==result.cols() && result.rows()==1,"Result must be scalar");
        for (size_t k=0; k!=elements.size(); k++)
            errors.at(elements.at(k)) += result(0,0)/elements.size();
    }
}

template <short_t d, class T, bool bending>
void gsThinShellAssemblerDWR<d, T, bending>::_applyLoadsToError(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, T & error) const
{

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if (m_pLoads[i].value.size()!=d)
            gsWarn<<"Point load has wrong dimension "<<m_pLoads[i].value.size()<<" instead of "<<d<<"\n";
        // Compute actives and values of basis functions on point load location.
        gsMatrix<T> forcePoint;
        if ( m_pLoads[i].parametric )   // in parametric space
            forcePoint = m_pLoads[i].point;
        else                            // in physical space
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);

        // Compute load value on the point
        gsMatrix<T> Lres, Hres;
        dualL.patch(m_pLoads[i].patch).eval_into(forcePoint,Lres);
        dualH.patch(m_pLoads[i].patch).eval_into(forcePoint,Hres);
        gsMatrix<T> result = (Hres-Lres).transpose() * m_pLoads[i].value;
        // Add load tot total error
        error += result(0,0);
    }
}


} // namespace gismo
