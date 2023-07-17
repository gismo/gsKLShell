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

namespace gismo
{
template <class T>
void gsShellStressFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),u.cols());

    gsExprEvaluator<T> ev;

    geometryMap m_ori   = ev.getMap(*m_patches);
    geometryMap m_def   = ev.getMap(*m_defpatches);

    gsMatrix<T> z(1,1);
    z.setZero();

    // gsMaterialMatrixEval<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,m_patches,m_defpatches,z);
    // variable S0 = ev.getVariable(m_S0);
    // gsMaterialMatrixEval<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,m_patches,m_defpatches,z);
    // variable S1 = ev.getVariable(m_S1);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_S0(m_materialMatrices,m_patches,m_defpatches);
    variable S0 = ev.getVariable(m_S0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_S1(m_materialMatrices,m_patches,m_defpatches);
    variable S1 = ev.getVariable(m_S1);
    gsMaterialMatrixIntegrate<T,MaterialOutput::PStressN> m_Sp0(m_materialMatrices,m_patches,m_defpatches);
    variable Sp0 = ev.getVariable(m_Sp0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::PStressM> m_Sp1(m_materialMatrices,m_patches,m_defpatches);
    variable Sp1 = ev.getVariable(m_Sp1);
    gsMaterialMatrixIntegrate<T,MaterialOutput::Stretch> m_lambda(m_materialMatrices,m_patches,m_defpatches);
    variable lambda = ev.getVariable(m_lambda);
    gsMaterialMatrixIntegrate<T,MaterialOutput::StretchDir> m_lambdadir(m_materialMatrices,m_patches,m_defpatches);
    variable lambdadir = ev.getVariable(m_lambdadir);

    gsFunctionExpr<> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    variable m_m12 = ev.getVariable(mult12t);

    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    auto E_m    = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) * reshape(m_m12,3,3) * That.tr();
    auto E_f    = ( expr::deriv2(m_ori,sn(m_ori).normalized().tr()) - expr::deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m12,3,3) * That.tr();

    auto S_m    = S0.tr() * Ttilde.tr();
    auto S_f    = S1.tr() * Ttilde.tr();

    gsMatrix<T> tmp;

    switch (m_stress_type)
    {
        case stress_type::displacement :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(m_def,u.col(k),m_patchID);
            break;

        case stress_type::membrane :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(S_m,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::flexural :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(S_f,u.col(k),m_patchID)).transpose();
            break;

        // -------------------------------------
        case stress_type::von_mises :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                gsMatrix<> S;
                gsMatrix<> Sm = (ev.eval(S_m,u.col(k),m_patchID)).transpose();
                gsMatrix<> Sf = (ev.eval(S_f,u.col(k),m_patchID)).transpose();
                S = Sm + Sf;
                result(0,k) = math::sqrt(S(0,0)*S(0,0)+S(1,0)*S(1,0)-S(0,0)*S(1,0)+3*S(2,0)*S(2,0)); // ASSUMES PLANE STRESS
                S = Sm - Sf;
                result(1,k) = math::sqrt(S(0,0)*S(0,0)+S(1,0)*S(1,0)-S(0,0)*S(1,0)+3*S(2,0)*S(2,0)); // ASSUMES PLANE STRESS
            }
            break;

        case stress_type::von_mises_membrane :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                gsMatrix<> S = (ev.eval(S_m,u.col(k),m_patchID)).transpose();
                result(0,k) = math::sqrt(S(0,0)*S(0,0)+S(1,0)*S(1,0)-S(0,0)*S(1,0)+3*S(2,0)*S(2,0)); // ASSUMES PLANE STRESS
            }
            break;

        case stress_type::von_mises_flexural :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                gsMatrix<> S = (ev.eval(S_f,u.col(k),m_patchID)).transpose();
                result(0,k) = math::sqrt(S(0,0)*S(0,0)+S(1,0)*S(1,0)-S(0,0)*S(1,0)+3*S(2,0)*S(2,0)); // ASSUMES PLANE STRESS
            }
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
