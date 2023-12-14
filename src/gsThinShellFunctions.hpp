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
    gsMatrix<T> Z(1,1);
    Z.setZero();

    result.setZero(targetDim(),u.cols());
    gsMatrix<T> z(1,1);
    z.setZero();

    gsExprEvaluator<T> ev;

    geometryMap m_ori   = ev.getMap(*m_patches);
    geometryMap m_def   = ev.getMap(*m_defpatches);

    // Initialize stystem
    // m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_N0(m_materialMatrices,m_defpatches);
    variable N0 = ev.getVariable(m_N0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_N1(m_materialMatrices,m_defpatches);
    variable N1 = ev.getVariable(m_N1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::CauchyVectorN> m_Nc0(m_materialMatrices,m_defpatches);
    variable Nc0 = ev.getVariable(m_Nc0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::CauchyVectorM> m_Nc1(m_materialMatrices,m_defpatches);
    variable Nc1 = ev.getVariable(m_Nc1);

    gsMaterialMatrixEval<T,MaterialOutput::StressN> m_S0(m_materialMatrices,m_defpatches,z);
    variable S0 = ev.getVariable(m_S0);
    gsMaterialMatrixEval<T,MaterialOutput::StressM> m_S1(m_materialMatrices,m_defpatches,z);
    variable S1 = ev.getVariable(m_S1);

    gsMaterialMatrixEval<T,MaterialOutput::CauchyStressN> m_Sc0(m_materialMatrices,m_defpatches,z);
    variable Sc0 = ev.getVariable(m_Sc0);
    gsMaterialMatrixEval<T,MaterialOutput::CauchyStressM> m_Sc1(m_materialMatrices,m_defpatches,z);
    variable Sc1 = ev.getVariable(m_Sc1);

    gsMaterialMatrixEval<T,MaterialOutput::PStress> m_Sp(m_materialMatrices,m_defpatches,z);
    variable Sp = ev.getVariable(m_Sp);
    gsMaterialMatrixEval<T,MaterialOutput::PStressDir> m_pstressdir(m_materialMatrices,m_defpatches,z);
    variable pstressdir = ev.getVariable(m_pstressdir);

    gsMaterialMatrixEval<T,MaterialOutput::PCauchyStressN> m_Sp0(m_materialMatrices,m_defpatches,z);
    variable Sp0 = ev.getVariable(m_Sp0);
    gsMaterialMatrixEval<T,MaterialOutput::PCauchyStressM> m_Sp1(m_materialMatrices,m_defpatches,z);
    variable Sp1 = ev.getVariable(m_Sp1);

    gsMaterialMatrixEval<T,MaterialOutput::PStrainN> m_Ep0(m_materialMatrices,m_defpatches,z);
    variable Ep0 = ev.getVariable(m_Ep0);
    gsMaterialMatrixEval<T,MaterialOutput::PStrainM> m_Ep1(m_materialMatrices,m_defpatches,z);
    variable Ep1 = ev.getVariable(m_Ep1);

    gsMaterialMatrixEval<real_t,MaterialOutput::TensionField> m_TF(m_materialMatrices,m_defpatches,z);
    variable TF  = ev.getVariable(m_TF);

    gsMaterialMatrixEval<T,MaterialOutput::Stretch> m_lambda(m_materialMatrices,m_defpatches,z);
    variable lambda = ev.getVariable(m_lambda);
    gsMaterialMatrixEval<T,MaterialOutput::StretchDir> m_lambdadir(m_materialMatrices,m_defpatches,z);
    variable lambdadir = ev.getVariable(m_lambdadir);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = ev.getVariable(mult2t);
    gsFunctionExpr<> mult12t("1","0","0","0","1","0","0","0","0.5",2);
    variable m_m12 = ev.getVariable(mult12t);

    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    auto Ttilde_def = cartcov(m_def);
    // auto Tmat   = cartcov(m_def);
    auto E_m    = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) * reshape(m_m12,3,3) * That.tr();
    auto E_f    = ( expr::deriv2(m_ori,sn(m_ori).normalized().tr()) - expr::deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) * That.tr();

    auto N_m    = N0.tr() * Ttilde.tr();
    auto N_f    = N1.tr() * Ttilde.tr();

    auto Nc_m   = Nc0.tr() * Ttilde_def.tr();
    auto Nc_f   = Nc1.tr() * Ttilde_def.tr();

    auto S_m    = S0.tr() * Ttilde.tr();
    auto S_f    = S1.tr() * Ttilde.tr();

    auto Sc_m   = Sc0.tr() * Ttilde_def.tr();
    auto Sc_f   = Sc1.tr() * Ttilde_def.tr();

    gsMatrix<T> tmp;

    switch (m_stress_type)
    {
        case stress_type::displacement :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(m_def-m_ori,u.col(k),m_patchID);
            break;

        case stress_type::membrane :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(Sc_m,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::flexural :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(Sc_f,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::membrane_PK2 :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(S_m,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::flexural_PK2 :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(S_f,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::membrane_force :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(Nc_m,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::flexural_moment :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(Nc_f,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::membrane_force_PK2 :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(N_m,u.col(k),m_patchID)).transpose();
            break;

        case stress_type::flexural_moment_PK2 :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = (ev.eval(N_f,u.col(k),m_patchID)).transpose();
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
        case stress_type::principal_membrane_strain:
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(Ep0,u.col(k),m_patchID);
            break;
        case stress_type::principal_flexural_strain :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(Ep1,u.col(k),m_patchID);
            break;

        case stress_type::principal_stretch :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(lambda,u.col(k),m_patchID);
            break;
        case stress_type::principal_stress :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(Sp,u.col(k));

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
        case stress_type::principal_stress_dir1 :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(pstressdir,u.col(k));
                result.col(k) = tmp.reshape(3,3).col(0);
            }
            break;
        case stress_type::principal_stress_dir2 :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(pstressdir,u.col(k));
                result.col(k) = tmp.reshape(3,3).col(1);
            }
            break;
        case stress_type::principal_stress_dir3 :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(pstressdir,u.col(k));
                result.col(k) = tmp.reshape(3,3).col(2);
            }
            break;
        case stress_type::tension_field :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                result.col(k) = ev.eval(TF,u.col(k),m_patchID);
            }
            break;
        case stress_type::total :
            gsWarn<<"Stress type 'total' not implemented\n";
            break;
    }
}

} // namespace gismo ends
