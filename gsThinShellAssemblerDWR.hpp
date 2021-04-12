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
#include <gsKLShell/gsMaterialMatrixEval.h>

namespace gismo
{

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsThinShellAssemblerDWR<d, T, bending, GF>::gsThinShellAssemblerDWR(
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
    _setBasis(m_basisL);
    m_dL = gsVector<T>::Zero(Base::m_assembler.numDofs());
    _setBasis(m_basisH);
    m_dH = gsVector<T>::Zero(Base::m_assembler.numDofs());
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::_setBasis(const gsMultiBasis<T> & basis)
{
    Base::m_basis = basis;
    Base::_initialize();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF>::_assembleMass(const gsMultiBasis<T> & basis, bool lumped)
{
    this->_setBasis(basis);
    Base::assembleMass(lumped);
    return Base::matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF>::_assembleMatrix(const gsMultiBasis<T> & basis)
{
    this->_setBasis(basis);
    Base::assemble();
    return Base::matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsSparseMatrix<T> gsThinShellAssemblerDWR<d, T, bending, GF>::_assembleMatrix(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(basis);
    Base::assemble(deformed);
    return Base::matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending, GF>::_assemblePrimal(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(basis);
    Base::assemble(deformed);
    return Base::rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsVector<T> gsThinShellAssemblerDWR<d, T, bending, GF>::_assemblePrimal(const gsMultiBasis<T> & basis)
{
    this->_setBasis(basis);
    Base::assemble();
    return Base::rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::_assembleDual_expr(auto expr, const gsMultiPatch<T> & deformed, const gsMultiBasis<T> & basis)
{
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);
    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initVector(1,false);

    Base::homogenizeDirichlet();

    // Assemble vector
    Base::m_assembler.assemble(expr);

}


template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::DisplacementNorm, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);
    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space  space       = Base::m_assembler.trialSpace(0);
    geometryMap Gori   = Base::m_assembler.exprData()->getMap();

    auto expr = space * gismo::expr::uv(2,3) * meas(Gori);
    Base::m_assembler.assemble(expr);

    // _assembleDual_expr(expr,deformed,basis);

}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::Displacement, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space   space   = Base::m_assembler.trialSpace(0);
    variable usol   = Base::m_assembler.getCoeff(primal);
    geometryMap Gori= Base::m_assembler.exprData()->getMap();

    auto expr = 2 * space * usol * meas(Gori);
    Base::m_assembler.assemble(expr);

}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStrain, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = Base::m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        Base::m_defpatches.patch(k).coefs() += primal.patch(k).coefs();;  // defG points to mp_def, therefore updated

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space  space       = Base::m_assembler.trialSpace(0);
    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ; //[checked]
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto expr = 2 * Em_der * Em.tr() * meas(Gori);
    Base::m_assembler.assemble(expr);

    // _assembleDual_expr(expr,deformed,basis);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStress, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = Base::m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        Base::m_defpatches.patch(k).coefs() += primal.patch(k).coefs();;  // defG points to mp_def, therefore updated

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(Base::m_materialMat,Base::m_defpatches);
    variable S0  = Base::m_assembler.getCoeff(S0f);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();

    auto N        = S0.tr();

    auto expr = N * gismo::expr::uv(0,3) * meas(Gori);
    Base::m_assembler.assemble(expr);

    // _assembleDual_expr(expr,deformed,basis);
}

// template <short_t d, class T, bool bending, enum GoalFunction GF>
// gsVector<T> gsThinShellAssemblerDWR<d, T, bending, GF>::_assembleDual_point(const gsMatrix<T> & points, const gsMultiBasis<T> & basis, const gsMatrix<T> & tmp)
// {
//     /* SINGLE PATCH IMPLEMENTATION */
//     index_t pIndex = 0;
//     gsVector<T> result;
//     gsMatrix<index_t> actives, globalActives;

//     space  space       = Base::m_assembler.trialSpace(0);

//     result.resize(Base::m_assembler.numDofs());
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::DisplacementNorm, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMatrix<T> & points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);
    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space  space       = Base::m_assembler.trialSpace(0);
    geometryMap Gori   = Base::m_assembler.exprData()->getMap();

    auto expr = space * gismo::expr::uv(2,3) * meas(Gori);
    gsExprEvaluator<T> ev(Base::m_assembler);

    result.resize(Base::m_assembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        basis.basis(pIndex).active_into( points.col(k), actives );
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::Displacement, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMatrix<T> & points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space   space   = Base::m_assembler.trialSpace(0);
    variable usol   = Base::m_assembler.getCoeff(primal);
    geometryMap Gori= Base::m_assembler.exprData()->getMap();

    auto expr = 2 * space * usol * meas(Gori);
    gsExprEvaluator<T> ev(Base::m_assembler);

    result.resize(Base::m_assembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        basis.basis(pIndex).active_into( points.col(k), actives );
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStrain, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMatrix<T> & points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = Base::m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        Base::m_defpatches.patch(k).coefs() += primal.patch(k).coefs();;  // defG points to mp_def, therefore updated

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space  space       = Base::m_assembler.trialSpace(0);
    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();

    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ; //[checked]
    auto Em_der = flat( jac(Gdef).tr() * jac(space) ) ;

    auto expr = 2 * Em_der * Em.tr() * meas(Gori);
    gsExprEvaluator<T> ev(Base::m_assembler);

    result.resize(Base::m_assembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        basis.basis(pIndex).active_into( points.col(k), actives );
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStress, gsVector<T>>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDual_impl(const gsMatrix<T> & points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
{
    /* SINGLE PATCH IMPLEMENTATION */
    index_t pIndex = 0;
    gsVector<T> result;
    gsMatrix<T> tmp;
    gsMatrix<index_t> actives, globalActives;
    this->_setBasis(basis);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = Base::m_patches;
    for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
        Base::m_defpatches.patch(k).coefs() += primal.patch(k).coefs();;  // defG points to mp_def, therefore updated

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    space  space       = Base::m_assembler.trialSpace(0);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(Base::m_materialMat,Base::m_defpatches);
    variable S0  = Base::m_assembler.getCoeff(S0f);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();

    auto N        = S0.tr();

    auto expr = N * gismo::expr::uv(0,3) * meas(Gori);
    gsExprEvaluator<T> ev(Base::m_assembler);

    result.resize(Base::m_assembler.numDofs());
    result.setZero();

    for (index_t k = 0; k!=points.cols(); k++)
    {
        tmp = ev.eval(expr,points.col(k));
        basis.basis(pIndex).active_into( points.col(k), actives );
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    this->_setBasis(m_basisL);
    result = Base::constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result)
{
    this->_setBasis(m_basisH);
    result = Base::constructMultiPatch(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);
    Base::constructSolution(solVector,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisH);
    Base::constructSolution(solVector,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF>::constructSolutionL(const gsMatrix<T> & solVector)
{
    this->_setBasis(m_basisL);
    return Base::constructSolution(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF>::constructSolutionH(const gsMatrix<T> & solVector)
{
    this->_setBasis(m_basisH);
    return Base::constructSolution(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF>::constructDisplacementL(const gsMatrix<T> & solVector)
{
    this->_setBasis(m_basisL);
    return Base::constructDisplacement(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsMultiPatch<T> gsThinShellAssemblerDWR<d, T, bending, GF>::constructDisplacementH(const gsMatrix<T> & solVector)
{
    this->_setBasis(m_basisH);
    return Base::constructDisplacement(solVector);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
T gsThinShellAssemblerDWR<d, T, bending, GF>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori= Base::m_assembler.exprData()->getMap();
    variable F      = Base::m_assembler.getCoeff(*Base::m_forceFun, Gori);
    variable zsolL  = Base::m_assembler.getCoeff(dualL);
    variable zsolH  = Base::m_assembler.getCoeff(dualH);

    auto expr = (zsolH-zsolL).tr() * F * meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);

    return ev.integral( expr  );

    // _assembleDual_expr(expr,deformed,basis);
}


template <short_t d, class T, bool bending, enum GoalFunction GF>
T gsThinShellAssemblerDWR<d, T, bending, GF>::computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    return computeError_impl<d,bending>(dualL,dualH,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<int _d, bool _bending>
typename std::enable_if<_d==3 && _bending, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);
    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);
    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(Base::m_materialMat,Base::m_defpatches);
    gsMaterialMatrixEval<T,MaterialOutput::VectorM> S1f(Base::m_materialMat,Base::m_defpatches);
    variable S0  = Base::m_assembler.getCoeff(S0f);
    variable S1  = Base::m_assembler.getCoeff(S1f);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m2 = Base::m_assembler.getCoeff(mult2t);

    geometryMap Gori  = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef  = Base::m_assembler.exprData()->getMap2();
    variable F      = Base::m_assembler.getCoeff(*Base::m_forceFun, Gori);
    variable zsolL  = Base::m_assembler.getCoeff(dualL);
    variable zsolH  = Base::m_assembler.getCoeff(dualH);
        // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    this->homogenizeDirichlet();

    auto N        = S0.tr();
    // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
    auto Em_der   = flat( jac(Gdef).tr() * (fjac(zsolH) - fjac(zsolL)) ) ; //[checked]

    auto M        = S1.tr(); // output is a column
    // auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3); //[checked]
    auto Ef_der   = -( deriv2(zsolH,sn(Gdef).normalized().tr() )
                    -  deriv2(zsolL,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(zsolH,Gdef) )
                                                                 - deriv2(Gdef,var1(zsolL,Gdef) ) ) * reshape(m2,3,3); //[checked]

    auto Fext = (zsolH-zsolL).tr() * F * meas(Gori);
    auto Fint = ( N * Em_der.tr() + M * Ef_der.tr() ) * meas(Gori);

    auto expr = ( Fext - Fint ) * meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);
    T integral = ev.integral(expr);


    if (Base::m_foundInd)
    {
        variable foundation = Base::m_assembler.getCoeff(*Base::m_foundFun, Gori);
        GISMO_ASSERT(Base::m_foundFun->targetDim()==3,"Foundation function has dimension "<<Base::m_foundFun->targetDim()<<", but expected 3");

        integral += ev.integral( ( zsolH - zsolL ) * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    }
    if (Base::m_pressInd)
    {
        variable pressure = Base::m_assembler.getCoeff(*Base::m_pressFun, Gori);
        GISMO_ASSERT(Base::m_pressFun->targetDim()==1,"Pressure function has dimension "<<Base::m_pressFun->targetDim()<<", but expected 1");

        integral += ev.integral( pressure.val() * ( zsolH - zsolL ) * sn(Gdef).normalized() * meas(Gori) );
    }

    return integral;

}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<int _d, bool _bending>
typename std::enable_if<!(_d==3 && _bending), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);
    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);
    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(Base::m_materialMat,Base::m_defpatches);
    variable S0  = Base::m_assembler.getCoeff(S0f);

    geometryMap Gori  = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef  = Base::m_assembler.exprData()->getMap2();
    variable F      = Base::m_assembler.getCoeff(*Base::m_forceFun, Gori);

    variable zsolL  = Base::m_assembler.getCoeff(dualL);
    variable zsolH  = Base::m_assembler.getCoeff(dualH);
        // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    this->homogenizeDirichlet();

    auto N        = S0.tr();
    // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
    auto Em_der   = flat( jac(Gdef).tr() * (fjac(zsolH) - fjac(zsolL)) ) ; //[checked]

    auto Fext = (zsolH-zsolL).tr() * F * meas(Gori);
    auto Fint = ( N * Em_der.tr() ) * meas(Gori);

    auto expr = ( Fext - Fint ) * meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);
    T integral = ev.integral(expr);

    if (Base::m_foundInd)
    {
        variable foundation = Base::m_assembler.getCoeff(*Base::m_foundFun, Gori);
        GISMO_ASSERT(Base::m_foundFun->targetDim()==3,"Foundation function has dimension "<<Base::m_foundFun->targetDim()<<", but expected 3");

        integral += ev.integral( ( zsolH - zsolL ) * foundation.asDiag() * (Gdef - Gori) * meas(Gori) ); // [v_x,v_y,v_z] diag([k_x,k_y,k_z]) [u_x; u_y; u_z]
    }
    if (Base::m_pressInd)
    {
        variable pressure = Base::m_assembler.getCoeff(*Base::m_pressFun, Gori);
        GISMO_ASSERT(Base::m_pressFun->targetDim()==1,"Pressure function has dimension "<<Base::m_pressFun->targetDim()<<", but expected 1");

        integral += ev.integral( pressure.val() * ( zsolH - zsolL ) * sn(Gdef).normalized() * meas(Gori) );
    }

    return integral;

}

template <short_t d, class T, bool bending, enum GoalFunction GF>
T gsThinShellAssemblerDWR<d, T, bending, GF>::computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed)
{
    return computeErrorEig_impl<GF>(evPrimalL,evDualL,evDualH,dualL,dualH,deformed);
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<!(_GF==GoalFunction::Buckling || _GF==GoalFunction::Modal), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed)
{
    GISMO_ERROR("Eigenvalue error estimation not available for this goal function");
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<(_GF==GoalFunction::Buckling), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed)
{
    // this->_setBasis(m_basisL);
    // Base::m_assembler.cleanUp();
    // Base::m_assembler.setOptions(Base::m_options);
    // Base::m_defpatches = deformed;

    // Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    // Base::m_assembler.getMap(Base::m_defpatches);

    // // Initialize vector
    // Base::m_assembler.initSystem(false);
    // Base::m_assembler.initVector(1,false);

    // geometryMap Gori  = Base::m_assembler.exprData()->getMap();

    // gsExprEvaluator<T> ev(Base::m_assembler);
    // T First = evPrimalL * ev.integral( mass * meas(mapL) );
    // // T Second = -ev.integral( Fint  * meas(mapL) );
    // T Third = (evDualH-evDualL)*(ev.integral( mass2 * meas(mapL)) -1.0 );


    // T integral = ev.integral(expr);


    GISMO_NO_IMPLEMENTATION;

}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<(_GF==GoalFunction::Modal), T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);
    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);
    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::Density> rhof(Base::m_materialMat,Base::m_defpatches);
    variable rho  = Base::m_assembler.getCoeff(rhof);

    variable zsolL  = Base::m_assembler.getCoeff(dualL);
    variable zsolH  = Base::m_assembler.getCoeff(dualH);

    geometryMap Gori  = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef  = Base::m_assembler.exprData()->getMap2();

    gsExprEvaluator<T> ev(Base::m_assembler);
    T integral = 0.0;
    integral += evPrimalL * ev.integral( rho.val() * Gdef * ( zsolH - zsolL ) * meas(Gori) );
    // T Second = -ev.integral( Fint  * meas(Gori) );
    integral += (evDualH-evDualL)*(ev.integral( rho.val() * Gdef * Gdef.tr() * meas(Gori)) -1.0 );

    return integral;
}

//============================================================
template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::DisplacementNorm, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();

    auto expr = (Gdef-Gori).tr() * gismo::expr::uv(2,3)*meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::Displacement, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();

    auto expr = (Gdef - Gori).tr() * (Gdef - Gori) * meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStrain, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ; //[checked]

    auto expr = Em.tr() * Em * meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);

    return ev.integral( expr  );
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStress, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(Base::m_materialMat,Base::m_defpatches);
    variable S0  = Base::m_assembler.getCoeff(S0f);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();

    auto expr = S0.tr() * S0 * meas(Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);

    return ev.integral( expr  );
}

//============================================================
template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::DisplacementNorm, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();

    auto expr = (Gdef-Gori).tr() * gismo::expr::uv(2,3);

    gsExprEvaluator<T> ev(Base::m_assembler);
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::Displacement, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();

    auto expr = (Gdef - Gori).tr() * (Gdef - Gori);

    gsExprEvaluator<T> ev(Base::m_assembler);
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStrain, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    geometryMap Gori   = Base::m_assembler.exprData()->getMap();
    geometryMap Gdef   = Base::m_assembler.exprData()->getMap2();
    auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ; //[checked]

    auto expr = Em.tr() * Em;

    gsExprEvaluator<T> ev(Base::m_assembler);
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

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStress, T>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
{
    this->_setBasis(m_basisL);

    Base::m_assembler.cleanUp();
    Base::m_assembler.setOptions(Base::m_options);

    Base::m_defpatches = deformed;

    Base::m_assembler.getMap(Base::m_patches);           // this map is used for integrals
    Base::m_assembler.getMap(Base::m_defpatches);

    // Initialize vector
    Base::m_assembler.initSystem(false);
    Base::m_assembler.initVector(1,false);

    gsMaterialMatrixEval<T,MaterialOutput::VectorN> S0f(Base::m_materialMat,Base::m_defpatches);
    variable S0  = Base::m_assembler.getCoeff(S0f);

    auto expr = S0.tr() * S0;

    gsExprEvaluator<T> ev(Base::m_assembler);
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
}// namespace gismo
