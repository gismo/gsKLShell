/** @file gsThinShellFunctions.hpp

    @brief Provides evaluation function for stresses.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)

*/

#pragma once

#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsExprAssembler.h>

namespace gismo
{
template <class T>
void gsShellStressFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),u.cols());

    m_assembler.cleanUp();

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(*m_defpatches);

    // Initialize stystem
    // m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,m_defpatches);
    variable S0 = m_assembler.getCoeff(m_S0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,m_defpatches);
    variable S1 = m_assembler.getCoeff(m_S1);
    gsMaterialMatrixIntegrate<T,MaterialOutput::PStressN> m_Sp0(m_materialMatrices,m_defpatches);
    variable Sp0 = m_assembler.getCoeff(m_Sp0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::PStressM> m_Sp1(m_materialMatrices,m_defpatches);
    variable Sp1 = m_assembler.getCoeff(m_Sp1);
    gsMaterialMatrixIntegrate<T,MaterialOutput::Stretch> m_lambda(m_materialMatrices,m_defpatches);
    variable lambda = m_assembler.getCoeff(m_lambda);
    gsMaterialMatrixIntegrate<T,MaterialOutput::StretchDir> m_lambdadir(m_materialMatrices,m_defpatches);
    variable lambdadir = m_assembler.getCoeff(m_lambdadir);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    // auto Tmat   = cartcov(m_def);
    auto E_m    = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) * That;
    auto E_f    = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) * That;

    auto S_m    = S0.tr() * Ttilde;
    auto S_f    = S1.tr() * Ttilde;

    gsExprEvaluator<T> ev(m_assembler);
    gsMatrix<T> tmp;

    switch (m_stress_type)
    {
        case stress_type::membrane :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(S_m,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::flexural :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(S_f,u.col(k),m_patchID)).transpose();
            break;

        // TO BE IMPLEMENTED
        // -------------------------------------
        case stress_type::von_mises :
            GISMO_ERROR("stress_type::von_mises is not implemented");
            break;

        case stress_type::von_mises_membrane :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(S_m,u.col(k));
                result(0,k) = sqrt( tmp(0,0)*tmp(0,0) + tmp(0,1)*tmp(0,1) - tmp(0,0)*tmp(0,1) + 3*tmp(0,2)*tmp(0,2) );
            }
            break;

        case stress_type::von_mises_flexural :
            GISMO_ERROR("stress_type::von_mises_flexural is not implemented");
            break;

        case stress_type::total :
            GISMO_ERROR("stress_type::total is not implemented");
            break;
        // -------------------------------------

        case stress_type::membrane_strain :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(E_m.tr(),u.col(k),m_patchID);
            break;

        case stress_type::flexural_strain :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(E_f.tr(),u.col(k),m_patchID);
            break;
        case stress_type::principal_stretch :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(lambda,u.col(k),m_patchID);
            break;
        case stress_type::principal_stress_membrane :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(Sp0,u.col(k),m_patchID);
            break;
        case stress_type::principal_stress_flexural :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(Sp1,u.col(k),m_patchID);
            break;

        case stress_type::principal_stretch_dir1 :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(lambdadir,u.col(k),m_patchID);
                result.col(k) = tmp.reshape(3,3).col(0);
            }
            break;
        case stress_type::principal_stretch_dir2 :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(lambdadir,u.col(k),m_patchID);
                result.col(k) = tmp.reshape(3,3).col(1);
            }
            break;
        case stress_type::principal_stretch_dir3 :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(lambdadir,u.col(k),m_patchID);
                result.col(k) = tmp.reshape(3,3).col(2);
            }
            break;

    }
}

} // namespace gismo ends
