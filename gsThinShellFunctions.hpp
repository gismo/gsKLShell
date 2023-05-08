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
    gsMatrix<T> z(1,1);
    z.setZero();

    m_assembler.cleanUp();

    geometryMap m_ori   = m_assembler.getMap(m_patches);
    geometryMap m_def   = m_assembler.getMap(*m_defpatches);

    // Initialize stystem
    // m_assembler.initSystem(false);

    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> m_N0(m_materialMatrices,m_defpatches);
    variable N0 = m_assembler.getCoeff(m_N0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> m_N1(m_materialMatrices,m_defpatches);
    variable N1 = m_assembler.getCoeff(m_N1);

    gsMaterialMatrixIntegrate<T,MaterialOutput::CauchyVectorN> m_Nc0(m_materialMatrices,m_defpatches);
    variable Nc0 = m_assembler.getCoeff(m_Nc0);
    gsMaterialMatrixIntegrate<T,MaterialOutput::CauchyVectorM> m_Nc1(m_materialMatrices,m_defpatches);
    variable Nc1 = m_assembler.getCoeff(m_Nc1);

    gsMaterialMatrixEval<T,MaterialOutput::StressN> m_S0(m_materialMatrices,m_defpatches,z);
    variable S0 = m_assembler.getCoeff(m_S0);
    gsMaterialMatrixEval<T,MaterialOutput::StressM> m_S1(m_materialMatrices,m_defpatches,z);
    variable S1 = m_assembler.getCoeff(m_S1);

    gsMaterialMatrixEval<T,MaterialOutput::CauchyStressN> m_Sc0(m_materialMatrices,m_defpatches,z);
    variable Sc0 = m_assembler.getCoeff(m_Sc0);
    gsMaterialMatrixEval<T,MaterialOutput::CauchyStressM> m_Sc1(m_materialMatrices,m_defpatches,z);
    variable Sc1 = m_assembler.getCoeff(m_Sc1);

    gsMaterialMatrixEval<T,MaterialOutput::PCauchyStressN> m_Sp0(m_materialMatrices,m_defpatches,z);
    variable Sp0 = m_assembler.getCoeff(m_Sp0);
    gsMaterialMatrixEval<T,MaterialOutput::PCauchyStressM> m_Sp1(m_materialMatrices,m_defpatches,z);
    variable Sp1 = m_assembler.getCoeff(m_Sp1);

    gsMaterialMatrixEval<T,MaterialOutput::PStrainN> m_Ep0(m_materialMatrices,m_defpatches,z);
    variable Ep0 = m_assembler.getCoeff(m_Ep0);
    gsMaterialMatrixEval<T,MaterialOutput::PStrainM> m_Ep1(m_materialMatrices,m_defpatches,z);
    variable Ep1 = m_assembler.getCoeff(m_Ep1);

    gsMaterialMatrixEval<real_t,MaterialOutput::TensionField> m_TF(m_materialMatrices,m_defpatches,z);
    variable TF  = m_assembler.getCoeff(m_TF);

    gsMaterialMatrixEval<T,MaterialOutput::Stretch> m_lambda(m_materialMatrices,m_defpatches,z);
    variable lambda = m_assembler.getCoeff(m_lambda);
    gsMaterialMatrixEval<T,MaterialOutput::StretchDir> m_lambdadir(m_materialMatrices,m_defpatches,z);
    variable lambdadir = m_assembler.getCoeff(m_lambdadir);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    auto Ttilde_def = cartcov(m_def);
    // auto Tmat   = cartcov(m_def);
    auto E_m    = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) * That;
    auto E_f    = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) * That;

    auto N_m    = N0.tr() * Ttilde;
    auto N_f    = N1.tr() * Ttilde;

    auto Nc_m   = Nc0.tr() * Ttilde_def;
    auto Nc_f   = Nc1.tr() * Ttilde_def;

    auto S_m    = S0.tr() * Ttilde;
    auto S_f    = S1.tr() * Ttilde;

    auto Sc_m   = Sc0.tr() * Ttilde_def;
    auto Sc_f   = Sc1.tr() * Ttilde_def;

    gsExprEvaluator<T> ev(m_assembler);
    gsMatrix<T> tmp;

    switch (m_stress_type)
    {
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

        // TO BE IMPLEMENTED
        // -------------------------------------
        case stress_type::von_mises :
            GISMO_ERROR("stress_type::von_mises is not implemented");
            break;

        case stress_type::von_mises_membrane :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                tmp = ev.eval(Sc_m,u.col(k));
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
        case stress_type::tension_field :
            for (index_t k = 0; k != u.cols(); ++k)
            {
                result.col(k) = ev.eval(TF,u.col(k),m_patchID);
            }
            break;

    }
}

} // namespace gismo ends
