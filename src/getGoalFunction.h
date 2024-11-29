/** @file gsThinShellUtils.h

    @brief Utilities for gsThinShellAssembler. Mainly expressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

//! [Include namespace]
#include <gsKLShell/src/gsGoalFunction.h>

#  define return_t  auto

namespace gismo
{

enum class GoalFunction : short_t
{
    DisplacementNorm        = 10,
    DisplacementComponent   = 11,
    StretchNorm             = 20,
    StretchComponent        = 21,
    PStrainNorm             = 30,
    PStrainComponent        = 31,
    PStressNorm             = 40,
    PStressComponent        = 41,
    MembraneStrainNorm      = 50,
    MembraneStrainComponent = 51,
    MembraneStressNorm      = 60,
    MembraneStressComponent = 61,
    MembraneForceNorm       = 70,
    MembraneForceComponent  = 71,
    FlexuralStrainNorm      = 80,
    FlexuralStrainComponent = 81,
    FlexuralStressNorm      = 90,
    FlexuralStressComponent = 91,
    FlexuralMomentNorm      = 100,
    FlexuralMomentComponent = 101,
    Modal                   = 110,
    Buckling                = 120,
};

template <class T>
gsGoalFunctionBase<T> * getGoalFunction(enum GoalFunction GF, short_t component = 9)
{
    if      (GF==GoalFunction::DisplacementNorm)
        return new 
}

template<class T>
class gsGoalFunctionBase
{
public:
    gsGoalFunctionBase(gsExprAssembler<T> & assembler)
    :
    m_assembler(assembler)
    {}

protected:
    gsExprAssembler<T>  & m_assembler;
};



template<class T, GoalFunction _GoalFunction>
class gsGoalFunction
{
    auto function();
    auto derivative();
};

template<class T>
class gsGoalFunction<T,GoalFunction::DisplacementNorm> : public gsGoalFunctionBase<T>
{
    typedef gsGoalFunctionBase<T> Base;
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;

public:
    gsGoalFunction(       gsExprAssembler<T> & assembler,
                    const gsMultiPatch<T> & original,
                    const gsMultiPatch<T> & deformed)
    :
    Base(assembler),
    m_ori(original),
    m_def(deformed)
    {

    }
    
    // Check with C++11
    auto function() const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);
        auto expr = (Gdef-Gori).tr() * (Gdef-Gori) * meas(Gori);
        return expr;
    }

    auto derivative(typename gsExprAssembler<T>::space u) const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);
        auto expr = 2 * u * (Gdef-Gori) * meas(Gori);
        return expr;
    }

protected:
    using Base::m_assembler;
    const gsMultiPatch<T> & m_ori;
    const gsMultiPatch<T> & m_def;
};

template<class T>
class gsGoalFunction<T,GoalFunction::DisplacementComponent> : public gsGoalFunctionBase<T>
{
    typedef gsGoalFunctionBase<T> Base;
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;

public:
    gsGoalFunction(       gsExprAssembler<T> &  assembler,
                    const gsMultiPatch<T> &     original,
                    const gsMultiPatch<T> &     deformed,
                    const index_t               component)
    :
    Base(assembler),
    m_ori(original),
    m_def(deformed),
    m_com(component)
    {

    }

    auto function() const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);
        auto expr = (Gdef-Gori).tr() * gismo::expr::uv(m_com,3) * meas(Gori);
        return expr;
    }

    auto derivative(typename gsExprAssembler<T>::space u) const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        auto expr = u * gismo::expr::uv(m_com,3) * meas(Gori);
        return expr;
    }

protected:
    using Base::m_assembler;
    const gsMultiPatch<T> & m_ori;
    const gsMultiPatch<T> & m_def;
    index_t                 m_com;
};

template<class T>
class gsGoalFunction<T,GoalFunction::StretchNorm> : public gsGoalFunctionBase<T>
{
    typedef gsGoalFunctionBase<T> Base;
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;

public:
    gsGoalFunction(       gsExprAssembler<T> & assembler,
                    const gsMultiPatch<T> & original,
                    const gsMultiPatch<T> & deformed,
                    const gsMaterialMatrixContainer<T> & materialMatrices)
    :
    Base(assembler),
    m_ori(original),
    m_def(deformed),
    m_materialMatrices(materialMatrices)
    {
        m_mult12t = gsFunctionExpr<T>("1","0","0","0","1","0","0","0","0.5",2);
        m_Z.resize(1,1);
        m_Z.setZero();
        m_Tstretchf = gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>(m_materialMatrices,&m_def,m_Z);
    }

    auto function() const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);

        auto m12        = m_assembler.getCoeff(m_mult12t);
        auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
        auto Tstretch   = m_assembler.getCoeff(m_Tstretchf);
        auto lambda     = Cm     * reshape(Tstretch,3,3).tr();

        auto expr       = lambda.tr() * lambda * meas(Gori);
        return expr;
    }

    auto derivative(typename gsExprAssembler<T>::space u) const
    {

        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);

        auto m12 = m_assembler.getCoeff(m_mult12t);
        auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
        auto Cm_der = 2* flat( jac(Gdef).tr() * jac(u) ) * reshape(m12,3,3);

        auto Tstretch = m_assembler.getCoeff(m_Tstretchf);

        auto lambda     = Cm     * reshape(Tstretch,3,3).tr();
        auto lambda_der = Cm_der * reshape(Tstretch,3,3).tr();

        auto expr = 2 * lambda_der * lambda.tr() * meas(Gori);
        return expr;
    }

protected:
    using Base::m_assembler;
    const gsMultiPatch<T> & m_ori;
    const gsMultiPatch<T> & m_def;
    const gsMaterialMatrixContainer<T> & m_materialMatrices;

    gsFunctionExpr<T> m_mult12t;
    gsMatrix<T> m_Z;
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform> m_Tstretchf;

};

template<class T>
class gsGoalFunction<T,GoalFunction::StretchComponent> : public gsGoalFunctionBase<T>
{
    typedef gsGoalFunctionBase<T> Base;
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;

public:
    gsGoalFunction(       gsExprAssembler<T> & assembler,
                    const gsMultiPatch<T> & original,
                    const gsMultiPatch<T> & deformed,
                    const gsMaterialMatrixContainer<T> & materialMatrices,
                    const index_t               component)
    :
    Base(assembler),
    m_ori(original),
    m_def(deformed),
    m_materialMatrices(materialMatrices),
    m_com(component)
    {
        m_mult12t = gsFunctionExpr<T>("1","0","0","0","1","0","0","0","0.5",2);
        m_Z.resize(1,1);
        m_Z.setZero();
        m_Tstretchf = gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>(m_materialMatrices,&m_def,m_Z);
    }

    auto function() const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);

        auto m12        = m_assembler.getCoeff(m_mult12t);
        auto Cm         = flat( jac(Gdef).tr() * jac(Gdef) ) * reshape(m12,3,3);
        auto Tstretch   = m_assembler.getCoeff(m_Tstretchf);
        auto lambda     = Cm     * reshape(Tstretch,3,3).tr();

        auto expr       = lambda.tr() * gismo::expr::uv(m_com,3) * meas(Gori);
        return expr;
    }

    auto derivative(typename gsExprAssembler<T>::space u) const
    {

        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);

        auto m12 = m_assembler.getCoeff(m_mult12t);
        auto Cm_der = 2* flat( jac(Gdef).tr() * jac(u) ) * reshape(m12,3,3);
        auto Tstretch = m_assembler.getCoeff(m_Tstretchf);
        auto lambda_der = Cm_der * reshape(Tstretch,3,3).tr();

        auto expr = 2 * lambda_der * gismo::expr::uv(m_com,3) * meas(Gori);
        return expr;
    }

protected:
    using Base::m_assembler;
    const gsMultiPatch<T> & m_ori;
    const gsMultiPatch<T> & m_def;
    const gsMaterialMatrixContainer<T> & m_materialMatrices;
    index_t                 m_com;

    gsFunctionExpr<T> m_mult12t;
    gsMatrix<T> m_Z;
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform> m_Tstretchf;

};

template<class T>
class gsGoalFunction<T,GoalFunction::PStrainNorm> : public gsGoalFunctionBase<T>
{
    typedef gsGoalFunctionBase<T> Base;
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;

public:
    gsGoalFunction(       gsExprAssembler<T> & assembler,
                    const gsMultiPatch<T> & original,
                    const gsMultiPatch<T> & deformed,
                    const gsMaterialMatrixContainer<T> & materialMatrices)
    :
    Base(assembler),
    m_ori(original),
    m_def(deformed),
    m_materialMatrices(materialMatrices)
    {
        m_mult12t = gsFunctionExpr<T>("1","0","0","0","1","0","0","0","0.5",2);
        m_Z.resize(1,1);
        m_Z.setZero();
        m_Tstretchf = gsMaterialMatrixEval<T,MaterialOutput::StretchTransform>(m_materialMatrices,&m_def,m_Z);
    }

    auto function() const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);
        auto m12        = m_assembler.getCoeff(m_mult12t);
        auto Tstretch   = m_assembler.getCoeff(m_Tstretchf);
        auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
        auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();
        auto expr       = Emp.tr() * Emp * meas(Gori);
        return expr;
    }

    auto derivative(typename gsExprAssembler<T>::space u) const
    {
        geometryMap Gori = m_assembler.getMap(m_ori);
        geometryMap Gdef = m_assembler.getMap(m_def);
        auto m12        = m_assembler.getCoeff(m_mult12t);
        auto Tstretch   = m_assembler.getCoeff(m_Tstretchf);
        auto Em     = 0.5 * ( flat(jac(Gdef).tr()*jac(Gdef)) - flat(jac(Gori).tr()* jac(Gori)) ) ;
        auto Em_der = flat( jac(Gdef).tr() * jac(u) ) ;
        auto Emp    = (Em * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();
        auto Emp_der= (Em_der * reshape(m12,3,3)) * reshape(Tstretch,3,3).tr();
        auto expr = 2 * Emp_der * Emp.tr() * meas(Gori);
        return expr;
    }

protected:
    using Base::m_assembler;
    const gsMultiPatch<T> & m_ori;
    const gsMultiPatch<T> & m_def;
    const gsMaterialMatrixContainer<T> & m_materialMatrices;

    gsFunctionExpr<T> m_mult12t;
    gsMatrix<T> m_Z;
    gsMaterialMatrixEval<T,MaterialOutput::StretchTransform> m_Tstretchf;

};


/// TO DO: Others

/*
DO IT LIKE THIS
option 1
    gsGoalFunction<GoalFunction::Displacement> obj(xxx);
    GradientDescent<gsGoalFunction<GoalFunction::Displacement>> gd();
    gd.setGoalFunction(obj);

option 2 (e.g. gsExprAssembler)
    gsGoalFunction<GoalFunction::Displacement> obj(xxx);
    GradientDescent gd();
    gd.computeGoalFunction(obj); // automatically finds the type of obj (since only computeGoalFunction is templated)
*/

}
