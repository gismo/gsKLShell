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

#include <gsThinShell2/gsThinShellFunctions.h>

namespace gismo
{
template <class T>
void gsShellStressFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),u.cols());

    m_assembler.cleanUp();

    m_assembler.getMap(m_patches);           // this map is used for integrals
    m_assembler.getMap(m_defpatches);

    // Initialize stystem
    m_assembler.initSystem();

    gsMaterialMatrix m_S0 = m_materialMat;
    m_S0.makeVector(0);
    gsMaterialMatrix m_S1 = m_materialMat;
    m_S1.makeVector(1);

    variable S0 = m_assembler.getCoeff(m_S0);
    variable S1 = m_assembler.getCoeff(m_S1);

    geometryMap m_ori   = m_assembler.exprData()->getMap();
    geometryMap m_def   = m_assembler.exprData()->getMap2();

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
    variable m_m2 = m_assembler.getCoeff(mult2t);

    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    auto Tmat   = cartcov(m_def);
    auto E_m    = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) ) * That;
    auto E_f    = ( deriv2(m_ori,sn(m_ori).normalized().tr()) - deriv2(m_def,sn(m_def).normalized().tr()) ) * reshape(m_m2,3,3) * That;

    auto S_m    = S0;
    auto S_f    = S1;

    auto Gdef   = jac(m_def);
    auto Gori   = jac(m_ori);
    auto normalDef = sn(m_def);
    auto normalOri = sn(m_ori);

    gsExprEvaluator ev(m_assembler);
    switch (m_stress_type)
    {
        case stress_type::membrane :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(S_m.tr(),u.col(k));
            break;

        case stress_type::flexural :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(S_f.tr(),u.col(k));
            break;

        case stress_type::membrane_strain :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(E_m.tr(),u.col(k));
            break;

        case stress_type::flexural_strain :
            for (index_t k = 0; k != u.cols(); ++k)
                result.col(k) = ev.eval(E_f.tr(),u.col(k));
            break;
        case stress_type::principal_stretch :
            gsMatrix<> Aori(3,3);
            gsMatrix<> Adef(3,3);
            gsMatrix<> tmp(2,2);
            gsMatrix<> evs;
            Eigen::SelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;

            T J0;
            // for (index_t k = 0; k != u.cols(); ++k)
            // {
            //     Aori.setZero();
            //     Adef.setZero();
            //     tmp.setZero();
            //     tmp = ev.eval(Gori,u.col(k)).block(0,0,2,2);
            //     tmp = tmp.transpose()*tmp;
            //     Aori.block(0,0,2,2) = tmp;

            //     tmp = ev.eval(Gdef,u.col(k)).block(0,0,2,2);
            //     tmp = tmp.transpose()*tmp;
            //     Adef.block(0,0,2,2) = tmp;

            //     J0 = math::sqrt( Adef.block(0,0,2,2).determinant() / Aori.block(0,0,2,2).determinant() );

            //     Adef(2,2) = math::pow(J0,-2.0);
            //     eigSolver.compute(Adef);
            //     evs = eigSolver.eigenvalues();
            //     for (index_t r=0; r!=3; r++)
            //         result(r,k) = math::sqrt(evs(r,0));
            // }

            for (index_t k = 0; k != u.cols(); ++k)
            {
                Aori.setZero();
                Adef.setZero();
                tmp.setZero();
                tmp = ev.eval(Gori,u.col(k));
                tmp = tmp.transpose()*tmp;
                Aori = tmp;

                tmp = ev.eval(Gdef,u.col(k));
                tmp = tmp.transpose()*tmp;
                Adef = tmp;

                J0 = math::sqrt( Adef.determinant() / Aori.determinant() );

                eigSolver.compute(Adef);
                evs = eigSolver.eigenvalues();
                for (index_t r=0; r!=2; r++)
                    result(r,k) = math::sqrt(evs(r,0));
                result(2,k) = 1./J0;

            }
            break;
    }
}

} // namespace gismo ends
