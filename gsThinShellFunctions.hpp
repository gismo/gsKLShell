/** @file gsThinShellFunctions.hpp

    @brief Provides useful classes derived from gsFunction which can be used
    for visualization or coupling.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsThinShell2/gsThinShellFunctions.h>

namespace gismo
{
template <class T>
void gsShellStressFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),u.cols());
    gsExprEvaluator<T> ev(m_assembler);

    gsMaterialMatrix m_mm0 = m_materialMat;
    m_mm0.setMoment(0);
    gsMaterialMatrix m_mm2 = m_materialMat;
    m_mm2.setMoment(2);
    variable mm0 = m_assembler.getCoeff(m_mm0);
    variable mm2 = m_assembler.getCoeff(m_mm2);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    space       m_space = m_assembler.trialSpace(0);
    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();
    // variable m_thick = m_assembler.getCoeff(*m_thickFun, m_ori);

    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    auto Tmat   = cartcov(m_def);
    auto E_m     = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) * That;
    auto E_f     = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) * That;

    auto S_m    = (reshape(mm0,3,3) * E_m ) * Ttilde;
    auto S_f    = (reshape(mm2,3,3) * E_f ) * Ttilde;

    gsMatrix<T> tmp;
    tmp = ev.eval(S_m,u.col(0));
    gsDebugVar(tmp);

    // switch (m_stress_type)
    // {
    //     case stress_type::membrane :
    //         for (index_t k = 0; k != u.cols(); ++k)
    //             result.col(k) = ev.eval(S_m,u.col(k));
    //         break;

    //     case stress_type::flexural :
    //         for (index_t k = 0; k != u.cols(); ++k)
    //             result.col(k) = ev.eval(S_f,u.col(k));
    //         break;

    //     case stress_type::membrane_strain :
    //         for (index_t k = 0; k != u.cols(); ++k)
    //             result.col(k) = ev.eval(E_m,u.col(k));
    //         break;

    //     case stress_type::flexural_strain :
    //         for (index_t k = 0; k != u.cols(); ++k)
    //             result.col(k) = ev.eval(E_f,u.col(k));
    //         break;
    // }
}

} // namespace gismo ends
