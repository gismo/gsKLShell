/** @file gsMaterialMatrixTFT.hpp

    @brief Provides linear material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/


#pragma once

#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>

using namespace gismo;

template <short_t d, class T>
class objectiveNL : public gsFunction<T>
{
public:

    objectiveNL()
    {}


    objectiveNL(gsMaterialMatrixBaseDim<d,T,true> * mm, index_t patch, const gsMatrix<T> & u, T z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> Ep, Sp;
        T theta, gamma, n1, n2, m1, m2;
        gsMatrix<T,2,2> n1_mat, n2_mat;
        gsVector<T,3> n1_vec, n2_vec;

        index_t colIdx;
        gsVector<T> gammas(1), thetas(1);

        for (index_t k=0; k!=u.cols(); k++)
        {
            thetas<<u(0,k);
            m_mm->setTheta(thetas);
            gammas<<u(1,k);
            m_mm->setGamma(gammas);

            gsMatrix<T> zMat(1,m_u.cols()); // init a matrix with 1 z-point per u-point
            zMat.setConstant(m_z);

            gsMatrix<T> Ss = m_mm->eval3D_vector(m_patch,m_u,zMat,MaterialOutput::Generic);

            n1 = math::cos(thetas.at(0));
            n2 = math::sin(thetas.at(0));
            m1 = -math::sin(thetas.at(0));
            m2 = math::cos(thetas.at(0));

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            result(0,k) = (n1_vec.transpose()*Ss)(0,0) * (n2_vec.transpose()*Ss)(0,0) + std::log(u(1,k));
            // result(1,k) = (n2_vec.transpose()*Ss)(0,0);

            // if (u(1,k) < 0)
            // {
            //     T tmp = math::exp(-u(1,k));
            //     for (index_t l = 0; l!=result.rows(); l++)
            //         result(l,k) = result(l,k)*(-tmp);
            //         // result(l,k) = std::abs(result(l,k))*(-tmp);
            //     // result.col(k).array() *= std::numeric_limits<T>::max();
            // }
        }
    }

    short_t domainDim() const
    {
        return 2;
    }

    short_t targetDim() const
    {
        return 1;
    }

private:
    gsMaterialMatrixBaseDim<d,T,true> * m_mm;
    index_t m_patch;
    const gsMatrix<T> m_u;
    T m_z;
};

namespace gismo
{

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixTFT<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    GISMO_ASSERT(out==MaterialOutput::MatrixA,"Tension Field Theory only works for membrane models, hence only outputs the A matrix");
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    gsMatrix<T> Cs = m_materialMat->eval3D_matrix(patch,u,z,out);

    gsMatrix<T> Es;
    if (out == MaterialOutput::MatrixA)
    {
        gsDebugVar("MatrixA");
        Es = m_materialMat->eval3D_strain(patch,u,z,MaterialOutput::VectorN); // returns flexural strain (although type is VectorN)
    }
    else if (out == MaterialOutput::MatrixD)
    {
        gsDebugVar("MatrixD");
        Es = m_materialMat->eval3D_strain(patch,u,z,MaterialOutput::VectorM); // returns flexural strain (although type is VectorM)
    }
    else
    {
        gsDebugVar("MatrixB/MatrixC");
        Es = gsMatrix<T>::Zero(3,Cs.cols()); // FOR LINEAR MATERIALS
    }

    gsMatrix<T> TF = m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);

    gsMatrix<T> result = Cs;

    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixTFT<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    GISMO_ASSERT(out==MaterialOutput::VectorN,"Tension Field Theory only works for membrane models, hence only outputs the N vector");

    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.rows())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    // WARNING: Check which values of out correspond to the output we want

    gsVector<T> gammas(1),thetas(1);
    gammas<<0;
    thetas<<0;
    m_materialMat->setTheta(thetas);
    m_materialMat->setGamma(gammas);

    gsMatrix<T> result = m_materialMat->eval3D_vector(patch,u,z,out);
    gsMatrix<T> TF = m_materialMat->eval3D_tensionfield(patch,u,z,MaterialOutput::TensionField);

    gsDebugVar(TF);

    index_t colIdx;
    T gamma, theta;
    for (index_t k=0; k!=u.cols(); k++)
    {
        for( index_t j=0; j < z.rows(); ++j ) // through-thickness points
        {
            colIdx = j*u.cols()+k;
            if (TF(0,colIdx) == 1) // taut
            {
                // do nothing
            }
            else if (TF(0,colIdx) == -1) // slack
            {
                // Set to zero
                result.col(colIdx).setZero();
            }
            else if (TF(0,colIdx) == 0) // wrinkled
            {
                index_t iter;
                gsVector<T> value(1);
                value<<0;
                gsVector<T> init(2);
                init<<0,0;

                objectiveNL<dim,T> obj(m_materialMat,patch,u.col(k),z(j,k));

                gsVector<T> x = gsVector<T>::LinSpaced(10,-0.5 * M_PI,0.5 * M_PI);
                std::vector<gsVector<T>> results; results.reserve(x.size());
                for (index_t k=0; k!=x.size(); k++)
                {
                    init.at(0) = x.at(k);
                    // init<< x.at(k),0;
                    iter = obj.newtonRaphson(value,init,false,1e-12,1000);

                    gsDebugVar(init);

                    if (init(0) > -M_PI && init(0) < M_PI && init(1) > 0)
                        results.push_back(init);
                }

                struct
                {
                    bool operator()(const gsVector<T> & a, const gsVector<T> & b) const
                    {
                        GISMO_ASSERT(a.size()==b.size(),"Sizes must agree!");
                        return std::lexicographical_compare(  a.begin(), a.end(),
                                                   b.begin(), b.end());
                    };
                }
                lexcomp;
                std::sort(results.begin(), results.end(),lexcomp);


                real_t comp_tol = 1e-5;
                auto comp = [&comp_tol](gsVector<T> & a, gsVector<T> & b)
                {
                    return (a-b).norm() < comp_tol;
                };

                auto last = std::unique(results.begin(), results.end(),comp);
                results.erase(last, results.end());

                GISMO_ASSERT(results.size()>0,"No suitable theta found");

                theta = results.at(0)(0);
                gamma = results.at(0)(1);



                // iter = obj.newtonRaphson(value,init,false,1e-12,1000);
                // gsDebugVar(iter);
                // GISMO_ASSERT(iter!=-1,"Newton-Raphson method did not converge!");

                // theta = init(0);
                // gamma = init(1);


                gsVector<T> gammas(1),thetas(1);
                gammas<<gamma;
                thetas<<theta;
                m_materialMat->setTheta(thetas);
                m_materialMat->setGamma(gammas);
                result.col(colIdx) = m_materialMat->eval3D_vector(patch,u,z,out);
            }
            else
                GISMO_ERROR("Tension field data not understood: "<<TF(0,colIdx));
        }
    }

    return result;
}

} // end namespace
