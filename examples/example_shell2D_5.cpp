/** @file example_shell2D.cpp

    @brief Simple 2D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixTFTLinear.h>
#include <gsKLShell/gsMaterialMatrixTFT.h>

#ifdef GISMO_WITH_IPOPT
#include <gsIpopt/gsOptProblem.h>
#endif
//#include <gsThinShell/gsNewtonIterator.h>


using namespace gismo;

template<typename T>
void toMatrix(gsMatrix<T> & voight)
{
    GISMO_ASSERT(voight.rows()==3 && voight.cols()==1, "Vector must be 3x1");
    voight.conservativeResize(4,1);
    std::swap(voight(1,0),voight(3,0));
    voight(2,0) *= 0.5;
    voight(1,0) = voight(2,0);
    voight.resize(2,2);
}

template <typename T>
class objective : public gsFunction<T>
{
public:
    objective(const gsMatrix<T> & C, const gsMatrix<T> & S)
    :
    m_C(C),
    m_S(S)
    {
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsVector<T,3> n1_vec, n2_vec;
        T theta, gamma;
        T n1, n2, m1, m2;

        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            gamma = - ( n1_vec.transpose() * m_S ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();

            T res = (n2_vec.transpose() * m_S).value() + gamma * (n2_vec.transpose() * m_C * n1_vec).value();

            result(0,k) = res;
        }
    }

    short_t domainDim() const
    {
        return 1;
    }

    short_t targetDim() const
    {
        return 1;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsVector<T,3> n1_vec, n2_vec, n3_vec, n4_vec;
        T theta, gamma, dgammadT;
        T n1, n2, m1, m2;

        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
            n4_vec<<m1*m1, m2*m2, 2*m1*m2;

            gamma = - ( n1_vec.transpose() * m_S ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();
            T tmp = ( n2_vec.transpose() * m_C * n1_vec ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();
            dgammadT = - 2 * gamma * tmp;
            // gsDebugVar(dgammadT);

            T res = (n3_vec.transpose() * m_S).value() + dgammadT * (n2_vec.transpose() * m_C * n1_vec).value()
                  + gamma * (n3_vec.transpose() * m_C * n1_vec).value() + 2*gamma * (n2_vec.transpose() * m_C * n2_vec).value();
            // T res = (n2_vec.transpose() * m_S).value() + gamma * (n2_vec.transpose() * m_C * n1_vec).value();

            result(0,k) = res;
        }
    }

private:
    const gsMatrix<T> m_C;
    const gsMatrix<T> m_S;
};

template <typename T>
class dummyCFunE : public gsFunction<T>
{
public:
    dummyCFunE() {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows()*u.rows()==targetDim(),"Something went wrong");
        result.resize(targetDim(),u.cols());
        for (index_t k=0; k!=u.cols(); k++)
        {
            gsAsMatrix<> C = result.reshapeCol(k,u.rows(),u.rows());
            for (index_t i=0; i!=u.rows(); i++)
                for (index_t j=0; j!=u.rows(); j++)
                    C(i,j) = u(i,k)*u(j,k);
        }

    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 9;
    }

    void deriv_mult_into(const gsMatrix<T> & E, const gsMatrix<T>& u, gsMatrix<T>&result)
    {
        GISMO_ASSERT(E.cols()==u.cols(),"Number of points is different");
        GISMO_ASSERT(E.rows()*E.rows()==targetDim(),"Dimensions are not correct");
        GISMO_ASSERT(E.rows()==domainDim(),"Dimensions are not correct");
        gsMatrix<T> ders;
        result.resize(targetDim(),u.cols());
        result.setZero();
        for (index_t k=0; k!=u.cols(); k++)
        {
            this->deriv_into(u.col(k),ders);
            ders = ders.transpose().blockTranspose(targetDim()).transpose();
            for (index_t d = 0; d!=ders.cols(); d++)
                result.col(k).segment(d*E.rows(),E.rows()) = ders.reshapeCol(d,E.rows(),E.rows()) * E;
        }
    }
};

template <typename T>
class CFunE : public gsFunction<T>
{
public:
    CFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T> value(1),init(1);
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        index_t iter;

        T theta;
        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            init.setZero();
            value.setZero();


            gsVector<T> x = gsVector<T>::LinSpaced(10,-0.5 * M_PI,0.5 * M_PI);
            std::vector<gsVector<T>> results; results.reserve(x.size());

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            for (index_t t=0; t!=x.size(); t++)
            {
                init.at(0) = x.at(t);

                objective<T> obj(C,S);
                iter = obj.newtonRaphson(value,init,false,1e-12,1000);
                GISMO_ASSERT(iter!=-1,"Newton iterations did not converge");

                theta = init.at(0);
                n1 = math::cos(theta);
                n2 = math::sin(theta);
                m1 = -math::sin(theta);
                m2 = math::cos(theta);

                n1_vec<<n1*n1, n2*n2, 2*n1*n2;
                n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

                gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

                gsVector<> res(2);
                res.at(0) = std::fmod(theta,M_PI);
                res.at(1) = gamma;

                if (res(1) > 0)
                {
                    results.push_back(res);
                }

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

            GISMO_ASSERT(results.size()>=1,"No suitable theta found");

            theta = results.at(0)(0);

            result(0,k) = theta;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 3;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T> value(1),init(1);
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        index_t iter;

        T theta;
        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            C = m_mm->C(E); // replace with C(E)
            result.reshapeCol(k,3,3) = C;
        }
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <typename T>
class thetaFunE : public gsFunction<T>
{
public:
    thetaFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T> value(1),init(1);
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        index_t iter;

        T theta;
        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            init.setZero();
            value.setZero();


            gsVector<T> x = gsVector<T>::LinSpaced(10,-0.5 * M_PI,0.5 * M_PI);
            std::vector<gsVector<T>> results; results.reserve(x.size());

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            for (index_t t=0; t!=x.size(); t++)
            {
                init.at(0) = x.at(t);

                objective<T> obj(C,S);
                iter = obj.newtonRaphson(value,init,false,1e-12,1000);
                GISMO_ASSERT(iter!=-1,"Newton iterations did not converge");

                theta = init.at(0);
                n1 = math::cos(theta);
                n2 = math::sin(theta);
                m1 = -math::sin(theta);
                m2 = math::cos(theta);

                n1_vec<<n1*n1, n2*n2, 2*n1*n2;
                n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

                gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

                gsVector<> res(2);
                res.at(0) = std::fmod(theta,M_PI);
                res.at(1) = gamma;

                if (res(1) > 0)
                {
                    results.push_back(res);
                }

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

            GISMO_ASSERT(results.size()>=1,"No suitable theta found");

            theta = results.at(0)(0);

            result(0,k) = theta;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 1;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T> value(1),init(1);
        gsVector<T,3> n1_vec, n2_vec, n3_vec;
        T n1, n2, m1, m2;

        index_t iter;

        T theta;
        T gamma;
        T dgammadT;
        T dfdT;
        gsMatrix<T> dgammadE;
        gsMatrix<T> dfdE;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);
            gsMatrix<T> Es = m_mm->eval3D_strain(m_patch,m_u,m_z,MaterialOutput::Generic);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            gsMatrix<T> res;
            this->eval_into(Es.col(0),res);

            theta = res(0,0);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);

            gsMatrix<T> I(3,3); I.setIdentity();

            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            T tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            dgammadT = - 2 * gamma * tmp;
            dgammadE = - ( n1_vec.transpose() * C * I) / ( n1_vec.transpose() * C * n1_vec ).value();

            // gsDebugVar(dgammadT);

            dfdT = (n3_vec.transpose() * S).value() + dgammadT * (n2_vec.transpose() * C * n1_vec).value()
                  + gamma * (n3_vec.transpose() * C * n1_vec).value() + 2*gamma * (n2_vec.transpose() * C * n2_vec).value();


            dfdE = n2_vec.transpose() * C * I + dgammadE * (n2_vec.transpose() * C * n1_vec).value();// + gamma * dCdE * n1_vec );

            gsDebugVar(dfdE);
            gsDebugVar(dfdT);
            // gsDebugVar(result.col(k));

            result.col(k) = -dfdE.transpose() / dfdT;

            ////// FOR LATER with dCdE

            // gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

            // T dCdE = 0;

            // // gsDebugVar(- ( n1_vec.transpose() * C ) / ( n1_vec.transpose() * C * n1_vec ).value());
            // // gsDebugVar(( n1_vec.transpose() * S * ( n1_vec.transpose() * dCdE * n1_vec ).value() ) / math::pow( (n1_vec.transpose() * C * n1_vec ).value(),2 ));

            // gsMatrix<T> dgammadE = - ( n1_vec.transpose() * C ) / ( n1_vec.transpose() * C * n1_vec ).value();
            //             //+ ( n1_vec.transpose() * S * ( n1_vec.transpose() * dCdE * n1_vec ).value() ) / math::pow( (n1_vec.transpose() * C * n1_vec ).value(),2 );

            // gsMatrix<T> I(3,3); I.setIdentity();
            // gsDebugVar(n2_vec.transpose() * ( C * I));
            // gsDebugVar(n2_vec.transpose() * (dgammadE * C * n1_vec).value());
            // gsMatrix<T> dTdE = n2_vec.transpose() * ( C * I) + n2_vec.transpose() * (dgammadE * C * n1_vec).value();//+ gamma * dCdE * n1_vec)
            // gsDebugVar(dTdE);
            // result.col(k) = dTdE.transpose();
        }
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <typename T>
class objFunE : public gsFunction<T>
{
public:
    objFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;
        T f;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);
            gsMatrix<T> Es = m_mm->eval3D_strain(m_patch,m_u,m_z,MaterialOutput::Generic);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
            gsMatrix<T> res;
            tfE.eval_into(Es.col(0),res);

            theta = res(0,0);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            gsDebugVar(gamma);
            f = (n2_vec.transpose() * C * n1_vec).value() + gamma * (n2_vec.transpose() * C * n1_vec).value();
            result(0,k) = f;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 1;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(this->targetDim()*this->domainDim(),u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;
        gsMatrix<T> dgammadE;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);
            gsMatrix<T> Es = m_mm->eval3D_strain(m_patch,m_u,m_z,MaterialOutput::Generic);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            T dCdE = 0;

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
            gsMatrix<T> res;
            tfE.eval_into(Es.col(0),res);

            theta = res(0,0);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

            dgammadE = - ( n1_vec.transpose() * C ) / ( n1_vec.transpose() * C * n1_vec ).value();
                    //+ ( n1_vec.transpose() * S * ( n1_vec.transpose() * dCdE * n1_vec ).value() ) / math::pow( (n1_vec.transpose() * C * n1_vec ).value(),2 );

            gsDebugVar(dgammadE * (n2_vec.transpose() * C * n1_vec).value());
            gsDebugVar(n2_vec.transpose() * C);

            gsMatrix<T> I(3,3); I.setIdentity();

            gsMatrix<T> dfdE = n2_vec.transpose() * C * I + dgammadE * (n2_vec.transpose() * C * n1_vec).value();// + gamma * dCdE * n1_vec );
            result.col(k) = dfdE.transpose();
        }
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <typename T>
class gammaFunE : public gsFunction<T>
{
public:
    gammaFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);
            gsMatrix<T> Es = m_mm->eval3D_strain(m_patch,m_u,m_z,MaterialOutput::Generic);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
            gsMatrix<T> res;
            tfE.eval_into(Es.col(0),res);

            theta = res(0,0);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            result(0,k) = gamma;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 1;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        gsMatrix<T> dgammadE;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);
            gsMatrix<T> Es = m_mm->eval3D_strain(m_patch,m_u,m_z,MaterialOutput::Generic);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            T dCdE = 0;

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
            gsMatrix<T> res;
            tfE.eval_into(Es.col(0),res);

            theta = res(0,0);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            dgammadE = - ( n1_vec.transpose() * C ) / ( n1_vec.transpose() * C * n1_vec ).value();
                    // + ( n1_vec.transpose() * S * ( n1_vec.transpose() * dCdE * n1_vec ).value() ) / math::pow( (n1_vec.transpose() * C * n1_vec ).value(),2 );

            gsDebugVar(result.col(k));
            gsDebugVar(dgammadE);

            result.col(k) = dgammadE.transpose();
        }
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <typename T>
class gammaFunT : public gsFunction<T>
{
public:
    gammaFunT(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, const gsMatrix<T> & E)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z),
    m_E(E)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E = m_E;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;

        gsMatrix<T> S, C;
        gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

        C = Cs.col(0); // replace with C(E)
        C.resize(3,3);
        toMatrix(E);
        S = m_mm->S(E);

        // Make voight notation
        S.resize(4,1);

        std::swap(S(3,0),S(1,0));
        S.conservativeResize(3,1);

        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            result(0,k) = gamma;
        }
    }

    short_t domainDim() const
    {
        return 1;
    }

    short_t targetDim() const
    {
        return 1;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E = m_E;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;
        T tmp;

        gsMatrix<T> S, C;
        gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

        C = Cs.col(0); // replace with C(E)
        C.resize(3,3);
        toMatrix(E);
        S = m_mm->S(E);

        // Make voight notation
        S.resize(4,1);

        std::swap(S(3,0),S(1,0));
        S.conservativeResize(3,1);

        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            result(0,k) = -2 * gamma * tmp;
        }
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
    const gsMatrix<T> m_E;
};


template <typename T>
class EwFunE : public gsFunction<T>
{
public:
    EwFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T,3> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;

        thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
        gsMatrix<T> res;
        tfE.eval_into(u,res);

        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);


            theta = res(0,k);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();


            result.col(k) = gamma * n1_vec;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 3;
    }

    void deriv_into2(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T> value(1),init(1);
        gsVector<T,3> n1_vec, n2_vec, n3_vec;
        T n1, n2, m1, m2;

        index_t iter;

        thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
        gsMatrix<T> res;
        tfE.eval_into(u,res);

        T theta;
        T gamma;
        T dgammadT;
        T dfdT;
        gsMatrix<T> dgammadE;
        gsMatrix<T> dfdE;
        gsMatrix<T> dTdE;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

            C = Cs.col(0); // replace with C(E)
            C.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            theta = res(0,k);

            gsDebugVar(theta);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);

            gsMatrix<T> I(3,3); I.setIdentity();

            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            gsDebugVar(gamma);
            T tmp = ( n2_vec.transpose() * C * n1_vec ).value() / ( n1_vec.transpose() * C * n1_vec ).value();
            dgammadT = - 2 * gamma * tmp;
            dgammadE = - ( n1_vec.transpose() * C * I) / ( n1_vec.transpose() * C * n1_vec ).value();
            dgammadE.transposeInPlace();
            gsDebugVar(dgammadE);

            dfdT = (n3_vec.transpose() * S).value() + dgammadT * (n2_vec.transpose() * C * n1_vec).value()
                  + gamma * (n3_vec.transpose() * C * n1_vec).value() + 2*gamma * (n2_vec.transpose() * C * n2_vec).value();


            dfdE = n2_vec.transpose() * C * I + dgammadE.transpose() * (n2_vec.transpose() * C * n1_vec).value();// + gamma * dCdE * n1_vec );

            dTdE = -dfdE.transpose() / dfdT;

            gsDebugVar(dfdE);
            gsDebugVar(dfdT);
            gsDebugVar(dTdE);
            // gsDebugVar(result.col(k));

            gsDebugVar(n1_vec.transpose()*dgammadE + ( n2_vec.transpose() * dgammadT + 2*gamma*n2_vec.transpose())*dTdE);

            result.col(k) = -dfdE.transpose() / dfdT;





            ////// FOR LATER with dCdE

            // gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();

            // T dCdE = 0;

            // // gsDebugVar(- ( n1_vec.transpose() * C ) / ( n1_vec.transpose() * C * n1_vec ).value());
            // // gsDebugVar(( n1_vec.transpose() * S * ( n1_vec.transpose() * dCdE * n1_vec ).value() ) / math::pow( (n1_vec.transpose() * C * n1_vec ).value(),2 ));

            // gsMatrix<T> dgammadE = - ( n1_vec.transpose() * C ) / ( n1_vec.transpose() * C * n1_vec ).value();
            //             //+ ( n1_vec.transpose() * S * ( n1_vec.transpose() * dCdE * n1_vec ).value() ) / math::pow( (n1_vec.transpose() * C * n1_vec ).value(),2 );

            // gsMatrix<T> I(3,3); I.setIdentity();
            // gsDebugVar(n2_vec.transpose() * ( C * I));
            // gsDebugVar(n2_vec.transpose() * (dgammadE * C * n1_vec).value());
            // gsMatrix<T> dTdE = n2_vec.transpose() * ( C * I) + n2_vec.transpose() * (dgammadE * C * n1_vec).value();//+ gamma * dCdE * n1_vec)
            // gsDebugVar(dTdE);
            // result.col(k) = dTdE.transpose();
        }
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <typename T>
class EwFunE4 : public gsFunction<T>
{
public:
    EwFunE4(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // voight notation of 4x1 tensors --> A11 A22 A12 A21

        result.resize(4,u.cols());
        gsMatrix<T> E, Evec;
        gsVector<T,4> n1_vec, n2_vec;
        T n1, n2, m1, m2;

        T theta;
        T gamma;

        thetaFunE<T> tfE(m_mm,m_patch,m_u,m_z);
        gsMatrix<T> res;
        tfE.eval_into(u,res);

        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = Evec = u.col(k);

            gsMatrix<T> S, C, Ctmp;
            gsMatrix<T> Cs = m_mm->eval3D_matrix(m_patch,m_u,m_z,MaterialOutput::MatrixA);

            Ctmp = Cs.col(0); // replace with C(E)
            Ctmp.resize(3,3);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);
            std::swap(S(3,0),S(1,0));
            gsDebugVar(S);

            E.resize(4,1);
            std::swap(E(3,0),E(1,0));
            gsDebugVar(E);

            C.resize(4,4);
            C.block(0,0,3,3) = Ctmp;
            C.row(3) = C.row(2);
            C.col(3) = C.col(2);


            theta = res(0,k);

            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, n1*n2, n2*n1;
            n2_vec<<m1*n1, m2*n2, m1*n2, m2*n1;
            gamma = - ( n1_vec.transpose() * S ).value() / ( n1_vec.transpose() * C * n1_vec ).value();


            result.col(k) = gamma * n1_vec;
        }
    }

    short_t domainDim() const
    {
        return 4;
    }

    short_t targetDim() const
    {
        return 4;
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <typename T>
class EfullFunE : public gsFunction<T>
{
public:
    EfullFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E;


        EwFunE<T> Ew(m_mm,m_patch,m_u,m_z);
        gsMatrix<T> res;
        Ew.eval_into(u,res);

        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);
            result.col(k) = E + res.col(k);
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 3;
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};


template <typename T>
class SfullFunE : public gsFunction<T>
{
public:
    SfullFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E,S;


        EfullFunE<T> Efull(m_mm,m_patch,m_u,m_z);
        gsMatrix<T> res;
        Efull.eval_into(u,res);

        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = res.col(k);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            result.col(k) = S;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 3;
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};


template <typename T>
class SFunE : public gsFunction<T>
{
public:
    SFunE(const gsMaterialMatrixBase<T> * mm, index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {
        GISMO_ASSERT(m_u.cols()==1 || m_z.cols()==1,"Currently only works for one point");
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E,S;


        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);
            toMatrix(E);
            S = m_mm->S(E);

            // Make voight notation
            S.resize(4,1);

            std::swap(S(3,0),S(1,0));
            S.conservativeResize(3,1);

            result.col(k) = S;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 3;
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
    const index_t m_patch;
    const gsMatrix<T> m_u;
    const gsMatrix<T> m_z;
};

template <short_t d, class T>
class objectiveNL2 : public gsFunction<T>
{
public:

    objectiveNL2()
    {}


    objectiveNL2(gsMaterialMatrixBaseDim<d,T,true> * mm, index_t patch, const gsMatrix<T> & u, T z)
    :
    m_mm(mm),
    m_patch(patch),
    m_u(u),
    m_z(z)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
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

            result(0,k) = (n1_vec.transpose()*Ss)(0,0);
            result(1,k) = (n2_vec.transpose()*Ss)(0,0);

            if (u(1,k) < 0)
            {
                T tmp = math::exp(-u(1,k));
                for (index_t l = 0; l!=result.rows(); l++)
                    result(l,k) = result(l,k)*(-tmp);
            }
        }
    }

    short_t domainDim() const
    {
        return 2;
    }

    short_t targetDim() const
    {
        return 2;
    }

private:
    gsMaterialMatrixBaseDim<d,T,true> * m_mm;
    index_t m_patch;
    const gsMatrix<T> m_u;
    T m_z;
};

template <typename T>
class gsTFTMat
{
public:
    gsTFTMat(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    { }

    void compute(T t)
    {
        theta = t;
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);
        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> denum = n1_vec.transpose() * m_C * n1_vec;

        C_I = m_C - 1 / (  n1_vec.transpose() * m_C * n1_vec ) * m_C * ( n1_vec * n1_vec.transpose() ) * m_C;

        T tmp = ( n2_vec.transpose() * m_C * n1_vec ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();

        gamma = - ( n1_vec.transpose() * m_C * m_e ).value() / ( n1_vec.transpose() * m_C * n1_vec ).value();

        // gsDebugVar((n1_vec.transpose() * m_C * n1_vec).value() + gamma * (n1_vec.transpose() * m_C * n1_vec).value());
        // gsDebugVar((n2_vec.transpose() * m_C * n1_vec).value() + gamma * (n2_vec.transpose() * m_C * n1_vec).value());

        gsMatrix<T> I(3,3); I.setIdentity();

        dgammadE = - ( n1_vec.transpose() * m_C * I ) / ( n1_vec.transpose() * m_C * n1_vec ).value();
        dgammadE.transposeInPlace();
        gsDebugVar(dgammadE);

        dgammadT = - 2 * gamma * tmp;
        gsDebugVar(dgammadT);

        dfdE = ( n2_vec.transpose() - tmp * n1_vec.transpose() ) * m_C * I;
        dfdE.transposeInPlace();

        dfdT = (n4_vec.transpose() * m_C * ( m_e + gamma * n1_vec )).value() +
                2 * gamma * ( n2_vec.transpose() * m_C * n2_vec - math::pow(( n1_vec.transpose() * m_C * n2_vec ).value(),2) / (n1_vec.transpose() * m_C * n1_vec).value() );

        dTdE = - dfdE / dfdT;

        gsMatrix<T,1,1> tmp2 = (n1_vec.transpose() * m_C * n2_vec);

        T df =  n4_vec.transpose() * m_C * (m_e + gamma * n1_vec)
                + 2 * gamma * ( n2_vec.transpose() * m_C * n2_vec
                - math::pow(tmp2(0,0),2) / (n1_vec.transpose() * m_C * n1_vec)
                )
                                ;

        GISMO_ASSERT(tmp2.rows()==1 && tmp2.cols()==1,"Must be scalar");



        gsMatrix<T,3,1> b = n2_vec - tmp * n1_vec;

        gsDebugVar(dfdE);
        gsDebugVar(b.transpose()*m_C);


        gsDebugVar(C_I);
        gsDebugVar(m_C + m_C * (n1_vec * dgammadE.transpose()));

        gsMatrix<T> C_I2 = m_C + m_C * (n1_vec * dgammadE.transpose());

        C_II = C_I + 2 * gamma / df * (m_C * b * b.transpose() * m_C);
        gsDebugVar(C_II);
        gsDebugVar(m_C * ( I + n1_vec * dgammadE.transpose() + dgammadT * n1_vec * dTdE.transpose() - 2*gamma * n2_vec * dTdE.transpose()));
        gsDebugVar(m_C * ( I + n1_vec * dgammadE.transpose() + dgammadT * n1_vec * dTdE.transpose() - 2*gamma * dTdE * n2_vec.transpose()));
        gsDebugVar(m_C * ( I + n1_vec * dgammadE.transpose() + dgammadT * n1_vec * dTdE.transpose() + 2*gamma * n2_vec * dTdE.transpose()));
        gsDebugVar(m_C * ( I + n1_vec * dgammadE.transpose() + dgammadT * n1_vec * dTdE.transpose() + 2*gamma * dTdE * n2_vec.transpose()));

        gsDebugVar(C_I2);
        gsDebugVar(m_C * dgammadT * n1_vec * dTdE.transpose());
        gsDebugVar(m_C * dgammadT * dTdE * n1_vec.transpose());
        gsDebug<<"\n";
        gsDebugVar(m_C * 2*gamma * n2_vec * dTdE.transpose());
        gsDebugVar(m_C * 2*gamma * dTdE * n2_vec.transpose());
        gsDebug<<"\n";
        gsDebugVar(m_C * 2*gamma * n2_vec * dTdE.transpose());

        Ew = gamma * n1_vec;
        E  = m_e + Ew;
        Sp = m_C * (m_e + Ew);

        dEw = n1_vec.transpose() * dgammadE
                + dgammadT * n1_vec.transpose() * dTdE
                + 2 * gamma * n2_vec.transpose() * dTdE;
    }

public:
    mutable gsMatrix<T> C_I, C_II;
    mutable T gamma, theta;
    mutable gsMatrix<T> Ew, E, Sp, dEw;
    mutable gsMatrix<T> dgammadE, dfdE, dTdE;
    mutable T dgammadT, dfdT;

private:
    const gsMatrix<T> m_C;
    const gsVector<T> m_e;
};

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool Compressibility = false;
    index_t material = 0;
    bool verbose = false;
    std::string fn;
    bool membrane = false;

    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("composite", "Composite material", composite);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    real_t L = 1;
    real_t B = 1;
    bool nonlinear = true;
    if (testCase == 0)
    {
        real_t mu = 1.5e6;
        thickness = 0.001;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);

        L = 1;
        B = 1;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= L;
        mp.patch(0).coefs().col(1) *= B;
        mp.addAutoBoundaries();
    }
    else if (testCase == 1)
    {
        E_modulus = 1;
        thickness = 1;
        // if (!Compressibility)
        //   PoissonRatio = 0.499;
        // else
          PoissonRatio = 0.3;

      L = 2;
      B = 1;

      mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
      mp.patch(0).coefs().col(0) *= L;
      mp.patch(0).coefs().col(1) *= B;
      mp.addAutoBoundaries();
    }
    else if (testCase == 2)
    {
        nonlinear = false;
        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0.3;

        L = 3;
        B = 1;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= L;
        mp.patch(0).coefs().col(1) *= B;
        mp.patch(0).coefs()(0,0) = B;

    }
    else if (testCase == 3)
    {
        nonlinear = true;
        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0.3;

        L = 3;
        B = 1;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= L;
        mp.patch(0).coefs().col(1) *= B;
    }
    else if (testCase == 4)
    {
        nonlinear = true;
        E_modulus = 100;
        thickness = 1;
        PoissonRatio = 0.3;

        L = 100;
        B = 100;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.patch(0).coefs().col(0) *= L;
        mp.patch(0).coefs().col(1) *= B;
    }


    gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";
    gsDebugVar(PoissonRatio);

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(2);
    tmp << 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    real_t pressure = 0.0;


    std::string fx = "0";
    std::string fy = "0";

    if (testCase == 0)
        fx = "2625";
    else if (testCase == 2)
    {
        real_t sigmax = 1e-1;
        char buffer[2000];
        sprintf(buffer,"%e ( 1 - y/%e)",sigmax,B);
        fx = buffer;
    }
    else if (testCase == 3)
        fx = "0.1";

    gsConstantFunction<> disply(0,2);
    gsFunctionExpr<> neuData(fx,fy,2);

    if (testCase == 0) // Uniaxial tension; use with hyperelastic material model!
    {
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
    }
    else if (testCase == 1)
    {
        for (index_t i=0; i!=2; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (2); load << 0.25, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 2)
    {
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
    }
    else if (testCase == 3)
    {
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

        bc.addCondition(boundary::north, condition_type::collapsed, 0, 0 ,false, 0);
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::north, condition_type::neumann, &neuData );
    }
    else if (testCase == 4)
    {
        disply.setValue(5,2);
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );

        bc.addCondition(boundary::east, condition_type::dirichlet, &disply, 0 ,false, 1);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false, 0);
        // bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false, 1);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (2); load << 10, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else
        GISMO_ERROR("Test case not known");

    //! [Refinement]

    // Linear isotropic material model and Neo-Hookean material
    gsConstantFunction<> force(tmp,2);
    gsConstantFunction<> pressFun(pressure,2);
    gsFunctionExpr<> t(std::to_string(thickness),2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
    gsFunctionExpr<> rho(std::to_string(Density),2);

    // Mooney-Rivlin material
    gsConstantFunction<> ratio(Ratio,2);

    // Ogden material
    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,2);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,2);
    gsConstantFunction<> alpha2(5.0,2);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,2);
    gsConstantFunction<> alpha3(-2.0,2);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,2);
    gsMaterialMatrixBase<real_t>* materialMatrix;

    // Linear anisotropic material model
    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,2);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,2);
    Ts[0] = &thicks;

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==1 || material==2) // NH & NH_ext
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==3) // MR
    {
      parameters.resize(3);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &ratio;
    }
    else if (material==4) // OG
    {
      parameters.resize(8);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
      parameters[4] = &mu2;
      parameters[5] = &alpha2;
      parameters[6] = &mu3;
      parameters[7] = &alpha3;
    }

    gsOptionList options;
    if      (material==0)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<2,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }

    gsMaterialMatrixLinear<2,real_t,true> materialMatrixTest(mp,t,parameters);

    gsMaterialMatrixTFTLinear<2,real_t> materialMatrixTFT(static_cast<gsMaterialMatrixLinear<2,real_t> *>(materialMatrix));
    gsMaterialMatrixTFT<2,real_t> materialMatrixTFT2(static_cast<gsMaterialMatrixLinear<2,real_t,true> *>(&materialMatrixTest));

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrix);

    assembler->setPointLoads(pLoads);

    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      time += stopwatch.stop();
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      time += stopwatch.stop();
      return assembler->rhs();
    };

    // Define Matrices
    stopwatch.restart();
    stopwatch2.restart();
    assembler->assemble();
    time += stopwatch.stop();

    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();

    // Solve linear problem
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);

    if (nonlinear)
    {
        real_t residual = vector.norm();
        real_t residual0 = residual;
        real_t residualOld = residual;
        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);
            residual = resVec.norm();

            gsInfo<<"Iteration: "<< it
               <<", residue: "<< residual
               <<", update norm: "<<updateVector.norm()
               <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
               <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
               <<"\n";

            residualOld = residual;

            if (updateVector.norm() < 1e-6)
                break;
            else if (it+1 == it)
                gsWarn<<"Maximum iterations reached!\n";
        }
    }

    totaltime += stopwatch2.stop();

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);
    }
    if (stress)
    {

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
        gsField<> Stretches(mp_def,stretches, true);

        gsPiecewiseFunction<> pstress_m;
        assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
        gsField<> pstressM(mp_def,pstress_m, true);

        gsPiecewiseFunction<> pstress_f;
        assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
        gsField<> pstressF(mp_def,pstress_f, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_def,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_def,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_def,stretch3, true);

        gsPiecewiseFunction<> VMStresses;
        assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
        gsField<> VMStress(mp_def,VMStresses, true);



        gsWriteParaview(membraneStress,"MembraneStress",5000);
        gsWriteParaview(VMStress,"MembraneStressVM",5000);
        gsWriteParaview(flexuralStress,"FlexuralStress",5000);
        gsWriteParaview(Stretches,"PrincipalStretch",5000);
        gsWriteParaview(pstressM,"PrincipalMembraneStress",5000);
        gsWriteParaview(pstressF,"PrincipalFlexuralStress",5000);
        gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
        gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
        gsWriteParaview(stretchDir2,"PrincipalDirection2",5000);
        gsWriteParaview(stretchDir3,"PrincipalDirection3",5000);


    }

    if (testCase==2)
    {
        gsPiecewiseFunction<> VMStresses;
        assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
        gsField<> VMStress(mp_def,VMStresses, true);
        gsVector<> pt(2);
        pt<<0,0;
        gsInfo<<"Stress in corner: "<<VMStresses.piece(0).eval(pt)<<"\n";
    }

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    gsMatrix<> z(1,1);
    z.setZero();
    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> mat(materialMatrix,&mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> vec(materialMatrix,&mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir> pdir(materialMatrix,&mp_def,z);

    gsMaterialMatrixEval<real_t,MaterialOutput::PStrainN> pstrain(materialMatrix,&mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::PStressN> pstress(materialMatrix,&mp_def,z);

    gsMaterialMatrixEval<real_t,MaterialOutput::TensionField> tensionfield(materialMatrix,&mp_def,z);

    gsMaterialMatrixEval<real_t,MaterialOutput::Strain > strain(materialMatrix,&mp_def,z);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> matTFT(&materialMatrixTFT,&mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> vecTFT(&materialMatrixTFT,&mp_def,z);

    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> matTFT2(&materialMatrixTFT2,&mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::VectorN> vecTFT2(&materialMatrixTFT2,&mp_def,z);


    gsVector<> Theta(1);
    Theta<<21.1939;
    gsVector<> Gamma(1);
    Gamma<<217.421;
    materialMatrixTest.setTheta(Theta);
    materialMatrixTest.setGamma(Gamma);
    gsMaterialMatrixEval<real_t,MaterialOutput::Strain> testE(&materialMatrixTest,&mp_def,z);


    gsField<> tensionField(mp_def, tensionfield, true);


    gsInfo<<"Plotting in Paraview...\n";
    gsWriteParaview<>( tensionField, "TensionField", 1000, true);

    gsVector<> pt(2);
    pt.setConstant(0.5);
    gsDebugVar(tensionfield.eval(pt));

    // gsDebugVar(matTFT.eval(pt).reshape(3,3));
    gsDebugVar(vecTFT.eval(pt));

    // // gsDebugVar(matTFT2.eval(pt).reshape(3,3));
    // gsDebugVar(vecTFT2.eval(pt));

    // gsDebugVar(testE.eval(pt));

    // return 0;

    gsMatrix<> result;
    mat.eval_into(pt,result);
    gsMatrix<real_t,3,3> C = result.reshape(3,3);
    gsDebugVar(C);

    pstrain.eval_into(pt,result);
    gsDebugVar(result);


    pstress.eval_into(pt,result);
    gsDebugVar(result);

    tensionfield.eval_into(pt,result);
    gsDebugVar(result);

    strain.eval_into(pt,result);
    gsVector<real_t> e = result;
    gsDebugVar(e);

    pdir.eval_into(pt,result);
    gsMatrix<> dir = result;
    gsDebugVar(dir.reshape(3,3));

    gsTFTMat<real_t> matCompute(C,e);

    gsMatrix<real_t,4,4> C4;
    C4.block(0,0,3,3) = C;
    C4.block(3,0,1,3) = C.row(2);
    C4.block(0,3,3,1) = C.col(2);
    C4(3,3) = C(2,2);
    gsDebugVar(C4);

    gsVector<> e4(4);
    e4.head(3) = e;
    e4.at(3) = e4.at(2) = 0.5 * e.at(2);
    gsDebugVar(e4);

    gsDebugVar(C4 * e4);
    gsDebugVar(C * e);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            THETA(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    thetaFunE<real_t> theta(materialMatrix,0,pt,z);
    theta.eval_into(e,result);
    gsDebugVar(result);

    gsMatrix<> THETA = result;

    matCompute.compute(THETA(0,0));


    gsDebugVar(THETA(0,0));
    gsDebugVar(matCompute.gamma);

    real_t n1, n2, m1, m2;
    n1 = math::cos(THETA(0,0));
    n2 = math::sin(THETA(0,0));
    m1 = -math::sin(THETA(0,0));
    m2 = math::cos(THETA(0,0));
    gsVector<> n1_vec(3), n2_vec(3);
    n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            f(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";
    objFunE<real_t> objE(materialMatrix,0,pt,z);
    objE.eval_into(e,result);
    gsDebugVar(result);
    objE.deriv_into(e,result);
    gsDebugVar(result);
    objE.deriv_into2(e,result);
    gsDebugVar(result);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            f(THETA)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsVector<> Stmp = C*e;
    objective<real_t> f(C,Stmp);
    f.deriv_into(THETA,result);
    gsDebugVar(result);
    f.deriv_into2(THETA,result);
    gsDebugVar(result);

    gsDebugVar(matCompute.dfdT);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            Dummy C(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";
    dummyCFunE<real_t> dummyC;
    gsMatrix<> dummyE(3,1); dummyE<<1,2,3;
    dummyC.eval_into(dummyE,result);
    gsDebugVar(result.reshape(3,3));

    // gsMatrix<> dCdE11, dCdE22, dCdE12;
    // for(index_t k=0; k!=results.rows(); k++)


    dummyC.deriv_into(dummyE,result);
gsDebugVar(result);
gsDebugVar(result.transpose().blockTranspose(3));
gsDebugVar(result.transpose().blockTranspose(9).transpose());
    result = result.transpose().blockTranspose(9).transpose();

    gsDebugVar(result.reshapeCol(0,3,3));
    gsDebugVar(result.reshapeCol(1,3,3));
    gsDebugVar(result.reshapeCol(2,3,3));

    dummyC.deriv_mult_into(dummyE,dummyE,result);

    // return 0;

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dTHETA(E) /dE                                               \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    // gsDebug<<"Numeric:\n";
    theta.deriv_into(e,result);
    gsDebugVar(result);
    theta.deriv_into2(e,result);
    gsDebugVar(result);

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.dTdE);

    gsInfo<<"Norm = "<<(result-matCompute.dTdE).norm()<<"\n";

    // gamma.deriv_into(e,result);
    // gsDebugVar(result);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            GAMMA(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    // gsDebug<<"Numeric:\n";
    gammaFunE<real_t> gamma(materialMatrix,0,pt,z);
    gamma.eval_into(e,result);
    gsDebugVar(result);

    // gsDebug<<"Analytical:\n";
    // gsDebugVar(matCompute.gamma);

    gsInfo<<"Norm = "<<std::abs(result.value()-matCompute.gamma)<<"\n";


    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dGAMMA(E) /dE                                               \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    gamma.deriv_into(e,result);
    gsDebugVar(result);
    gamma.deriv_into2(e,result);
    gsDebugVar(result);

    // gsDebug<<"Analytical:\n";
    // gsDebugVar(matCompute.dgammadE);

    gsInfo<<"Norm = "<<(result-matCompute.dgammadE).norm()<<"\n";

    // gamma.deriv_into(e,result);
    // gsDebugVar(result);

// return 0;


    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dGAMMA/dTHETA(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    // gsDebug<<"Numeric:\n";
    gammaFunT<real_t> gammaT(materialMatrix,0,pt,z,e);
    gammaT.deriv_into(THETA,result);
    gsDebugVar(result);

    gsMatrix<> result2;
    gammaT.deriv_into2(THETA,result2);
    gsDebugVar(result2);

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.dgammadT);


    // gamma.deriv_into(e,result);
    // gsDebugVar(result);

    gsInfo<<"Norm = "<<std::abs(result.value()-matCompute.dgammadT)<<"\n";


    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            Ew                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    EwFunE<real_t> Ew(materialMatrix,0,pt,z);
    Ew.eval_into(e,result);
    gsDebugVar(result);

    // ew = result;

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.Ew);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dEw/dE                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.dEw);

    gsDebug<<"Numeric:\n";
    Ew.deriv_into(e,result);
    gsDebugVar(result.reshape(3,3));
    gsDebugVar(result);

    gsDebug<<"Numeric:\n";
    result.resize(0,0);
    Ew.deriv_into2(e,result);
    gsDebugVar(result.reshape(3,3));
    gsDebugVar(result);


    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.dEw);
    return 0;


    gsDebugVar(C * result.reshape(3,3));
    gsDebugVar(C);
    // gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.C_II);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            E'                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsMatrix<> ep;

    gsDebug<<"Numeric:\n";
    EfullFunE<real_t> Efull(materialMatrix,0,pt,z);
    Efull.eval_into(e,result);
    gsDebugVar(result);

    ep = result;

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.E);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dE'/dE                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    Efull.deriv_into(e,result);
    gsDebugVar(result.reshape(3,3));

    gsDebugVar(C * result.reshape(3,3));

    // gsDebug<<"Analytical:\n";
    // gsDebugVar(matCompute.C_II);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            S'(E)                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    SfullFunE<real_t> Sfull(materialMatrix,0,pt,z);
    Sfull.eval_into(e,result);
    gsDebugVar(result);

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.Sp);

    // gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dS'/dE (E)                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    Sfull.deriv_into(e,result);
    gsDebugVar(result.reshape(3,3));

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.C_I );
    gsDebugVar(matCompute.C_II);


    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dS'(E)/dE *E                                                 \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    Sfull.deriv_into(e,result);
    gsDebugVar(result.reshape(3,3) * e);

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.C_I  * e);
    gsDebugVar(matCompute.C_II * e);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            S(E')                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    SFunE<real_t> S(materialMatrix,0,pt,z);
    S.eval_into(ep,result);
    gsDebugVar(result);

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.Sp);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dS'/dE                                                   \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    S.deriv_into(ep,result);
    gsDebugVar(result.reshape(3,3));

    gsDebug<<"Analytical:\n";
    gsDebugVar(C);

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dS'(E)/dE *E                                                 \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    S.deriv_into(ep,result);
    gsDebugVar(result.reshape(3,3) * ep);

    gsDebug<<"Analytical:\n";
    gsDebugVar(matCompute.C_I  * e);
    gsDebugVar(matCompute.C_II * e);





    gsMatrix<> dgammadE, dgammadT, dTdE;
    real_t gam;

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            GAMMA(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    gammaFunE<real_t> gammap(materialMatrix,0,pt,z);
    gammap.eval_into(e,result);

    gam = result.value();

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dGAMMA(E') /dE                                               \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    gamma.deriv_into(e,result);
    gsDebugVar(result);

    dgammadE = result;

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dGAMMA/dTHETA(E)                                                \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    gammaFunT<real_t> gammaTp(materialMatrix,0,pt,z,ep);
    gammaTp.deriv_into(THETA,result);
    gsDebugVar(result);

    dgammadT = result;

    gsDebug<<"---------------------------------------------------------------------------------\n";
    gsDebug<<"                            dTHETA(E)/dE                                               \n";
    gsDebug<<"---------------------------------------------------------------------------------\n";

    gsDebug<<"Numeric:\n";
    theta.deriv_into(e,result);
    gsDebugVar(result);

    dTdE = result;

    gsDebugVar(dgammadE * n1_vec.transpose());

    gsDebugVar(dgammadT.value() * dTdE * n1_vec.transpose());

    gsDebugVar(dTdE * n2_vec.transpose());

    gsDebugVar(dgammadE * n1_vec.transpose() + dgammadT.value() * dTdE * n1_vec.transpose() + dTdE * n2_vec.transpose());


    // gsDebug<<"---------------------------------------------------------------------------------\n";
    // gsDebug<<"                            S                                                    \n";
    // gsDebug<<"---------------------------------------------------------------------------------\n";

    // gsDebug<<"Numeric:\n";
    // SFun<real_t> S(materialMatrix);
    // S.eval_into(e4,result);
    // gsDebugVar(result);

    // gsDebug<<"Analytical:\n";

    // gsDebug<<"---------------------------------------------------------------------------------\n";
    // gsDebug<<"                            C                                                    \n";
    // gsDebug<<"---------------------------------------------------------------------------------\n";

    // gsDebug<<"Numeric:\n";
    // S.deriv_into(e4,result);
    // gsDebugVar(result.reshape(4,4));

    // gsDebug<<"Analytical:\n";
    // gsDebugVar(C );

    // // gsDebugVar(result.reshape(3,3) * e);

    // gsDebug<<"---------------------------------------------------------------------------------\n";
    // gsDebug<<"                            C*E                                                  \n";
    // gsDebug<<"---------------------------------------------------------------------------------\n";

    // gsDebug<<"Numeric:\n";
    // gsDebugVar(result.reshape(4,4) * e4);

    // gsDebug<<"Analytical:\n";
    // gsDebugVar(C * e );

    return 0;

    // value.resize(2);
    // value<<0,0;
    // init.resize(2);
    // init<<0,0;

    // objectiveNL<real_t> objNL(materialMatrix,e);
    // iter = objNL.newtonRaphson(value,init,false,1e-12,1000);
    // gsDebugVar(iter);
    // gsDebugVar(init);



    // gsMaterialMatrixBaseDim<2,real_t> baseDim(mp,mp_def);
    // baseDim._computeMetricUndeformed(0,pt,false);
    // baseDim._computeMetricDeformed(0,pt,false);

    // real_t lambda = E_modulus * PoissonRatio / ( (1. + PoissonRatio)*(1.-2.*PoissonRatio)) ;
    // real_t Cconstant = 2*lambda*mu/(lambda+2*mu);
    // gsDebugVar(Cconstant);
    // gsDebugVar(mu);

    // real_t theta = init.at(0);
    // real_t gamma = init.at(1);
    // gsMatrix<> Acon = baseDim._getAcon_ori(0,0.0);
    // gsDebugVar(Acon);
    // gsVector<> n(2);
    // n<<math::cos(theta),math::sin(theta);
    // gsMatrix<> Aconp= 2*gamma * n * n.transpose();
    // Acon += Aconp;
    // gsDebugVar(Acon);
    // gsDebugVar(Aconp);
    // gsDebugVar(n);

    // real_t n1 = math::cos(theta);
    // real_t n2 = math::sin(theta);

    // gsVector<> n1_vec(3);
    // n1_vec<<n1*n1, n2*n2, 2*n1*n2;
    // gsVector<> e_w = e + gamma*n1_vec;

    // gsDebugVar(e);
    // gsDebugVar(e_w);

    // // S_num.eval_into(e,result);
    // // gsDebugVar(result);

    // vectorS_e<real_t>  S_e (materialMatrix,gamma,theta,e);
    // vectorS_ew<real_t> S_ew(materialMatrix);

    // gsMatrix<> result1, result2, result3, result4;

    // S_e .eval_into(e,result1);
    // gsDebugVar(result1);
    // S_ew.eval_into(e_w,result2);
    // gsDebugVar(result2);

    // S_e .deriv_into(e,result3);
    // gsDebugVar(result3.reshape(3,3));
    // gsDebugVar(result3.reshape(3,3) * e);
    // gsDebugVar(result3.reshape(3,3) * e_w);

    // S_ew.deriv_into(e_w,result4);
    // gsDebugVar(result4.reshape(3,3));
    // gsDebugVar(result4.reshape(3,3) * e);
    // gsDebugVar(result4.reshape(3,3) * e_w);

    // gsDebugVar(C);


    // gsMatrix<> Cnew(3,3);
    // Cnew(0,0)             = Cconstant*Acon(0,0)*Acon(0,0) + mu*(Acon(0,0)*Acon(0,0) + Acon(0,0)*Acon(0,0)); // C1111
    // Cnew(1,1)             = Cconstant*Acon(1,1)*Acon(1,1) + mu*(Acon(1,1)*Acon(1,1) + Acon(1,1)*Acon(1,1)); // C2222
    // Cnew(2,2)             = Cconstant*Acon(0,1)*Acon(0,1) + mu*(Acon(0,0)*Acon(1,1) + Acon(0,1)*Acon(1,0)); // C1212
    // Cnew(1,0) = Cnew(0,1) = Cconstant*Acon(0,0)*Acon(1,1) + mu*(Acon(0,1)*Acon(0,1) + Acon(0,1)*Acon(0,1)); // C1122
    // Cnew(2,0) = Cnew(0,2) = Cconstant*Acon(0,0)*Acon(0,1) + mu*(Acon(0,0)*Acon(0,1) + Acon(0,1)*Acon(0,0)); // C1112
    // Cnew(2,1) = Cnew(1,2) = Cconstant*Acon(1,1)*Acon(0,1) + mu*(Acon(1,0)*Acon(1,1) + Acon(1,1)*Acon(1,0)); // C2212


    // gsMatrix<> Cnewp(3,3);
    // Cnewp(0,0)              = Cconstant*Aconp(0,0)*Aconp(0,0) + mu*(Aconp(0,0)*Aconp(0,0) + Aconp(0,0)*Aconp(0,0)); // C1111
    // Cnewp(1,1)              = Cconstant*Aconp(1,1)*Aconp(1,1) + mu*(Aconp(1,1)*Aconp(1,1) + Aconp(1,1)*Aconp(1,1)); // C2222
    // Cnewp(2,2)              = Cconstant*Aconp(0,1)*Aconp(0,1) + mu*(Aconp(0,0)*Aconp(1,1) + Aconp(0,1)*Aconp(1,0)); // C1212
    // Cnewp(1,0) = Cnewp(0,1) = Cconstant*Aconp(0,0)*Aconp(1,1) + mu*(Aconp(0,1)*Aconp(0,1) + Aconp(0,1)*Aconp(0,1)); // C1122
    // Cnewp(2,0) = Cnewp(0,2) = Cconstant*Aconp(0,0)*Aconp(0,1) + mu*(Aconp(0,0)*Aconp(0,1) + Aconp(0,1)*Aconp(0,0)); // C1112
    // Cnewp(2,1) = Cnewp(1,2) = Cconstant*Aconp(1,1)*Aconp(0,1) + mu*(Aconp(1,0)*Aconp(1,1) + Aconp(1,1)*Aconp(1,0)); // C2212

    // gsDebugVar(C);
    // gsDebugVar(Cnew);
    // gsDebugVar(Cnewp);

    // gsDebugVar(matCompute.C_I  * e);
    // gsDebugVar(matCompute.C_II * e);
    // gsDebugVar(Cnew * e);

    // objNL.eval_into(init,result);

    // gsVector<> x = gsVector<>::LinSpaced(100,-0.5 * M_PI,0.5 * M_PI);


    return 0;


// #ifdef GISMO_WITH_IPOPT
//     gsOptProblemExample<real_t> opt(C,e);

//     // for (index_t k = 0; k!=x.size(); k++)
//     // {
//     //     std::vector<real_t> X,dX,C,dC;
//     //     X.push_back(x[k]);
//     //     dX.resize(1);
//     //     gsAsVector<> der(dX);

//     //     C.resize(1);
//     //     gsAsVector<> con(C);

//     //     dC.resize(1);
//     //     gsAsVector<> dcon(dC);


//     //     opt.gradObj_into(gsAsConstVector<>(X),der);
//     //     opt.evalCon_into(gsAsConstVector<>(X),con);
//     //     opt.jacobCon_into(gsAsConstVector<>(X),dcon);
//     //     gsInfo<<x[k]<<","<<opt.evalObj(gsAsConstVector<>(X))<<","<<der<<","<<con<<","<<dcon<<"\n";
//     // }


//     // Run optimizer
//     opt.solve();

//     // Print some details in the output
//     gsDebugVar(opt.objective());
//     gsDebugVar(opt.iterations());
//     gsDebugVar(opt.currentDesign());
//     real_t theta = opt.currentDesign()(0,0);

//     matCompute.compute(theta);
//     gsDebugVar(matCompute.C_I);
//     gsDebugVar(matCompute.C_II);
//     gsDebugVar(matCompute.gamma);

// #else
//     gsInfo<<"GISMO_WITH_IPOPT is not enabled\n";
// #endif


    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////


    delete assembler;
    return EXIT_SUCCESS;

}// end main
