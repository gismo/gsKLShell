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
#include <gsKLShell/gsGLStrain.h>

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
    voight(1,0) = voight(2,0);
    voight.resize(2,2);
}

template <typename T>
class objective : public gsFunction<T>
{
public:
    objective(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        for (index_t k = 0; k!=u.cols(); k++)
        {
            T theta = u(0,k);
            T n1 = math::cos(theta);
            T n2 = math::sin(theta);
            T m1 = -math::sin(theta);
            T m2 = math::cos(theta);

            gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
            gsVector<T,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
            gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

            gsMatrix<T,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );

            gsMatrix<T,1,1> res = n2_vec.transpose() * m_C * m_e + gamma * n2_vec.transpose() * m_C * n1_vec;
            GISMO_ASSERT(res.rows()==1 && res.cols() ==1,"f is not scalar!");

            result(0,k) = res(0,0)*res(0,0);
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

private:
    const gsMatrix<T> m_C;
    const gsMatrix<T> m_e;
};


template <typename T>
class objectiveNL : public gsFunction<T>
{
public:

    objectiveNL()
    {}


    objectiveNL(const gsMaterialMatrixBase<T> * mm, const gsMatrix<T> & e)
    :
    m_mm(mm),
    m_e(e)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(2,u.cols());
        gsMatrix<T> Ep, Sp;
        T theta, gamma, n1, n2, m1, m2;
        gsVector<T,3> n1_vec, n2_vec;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            theta = u(0,k);
            gamma = u(1,k);
            n1 = math::cos(theta);
            n2 = math::sin(theta);
            m1 = -math::sin(theta);
            m2 = math::cos(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;
            n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;

            Ep = m_e + gamma * n1_vec;
            toMatrix(Ep);
            Sp = m_mm->S(Ep);

            // Make voight notation
            Sp.resize(4,1);

            std::swap(Sp(3,0),Sp(1,0));
            Sp.conservativeResize(3,1);

            result(0,k) = (n1_vec.transpose()*Sp)(0,0);
            result(1,k) = (n2_vec.transpose()*Sp)(0,0);
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
    const gsMaterialMatrixBase<T> * m_mm;
    const gsMatrix<T> m_e;
};

template <typename T>
class EfullFun : public gsFunction<T>
{
public:
    EfullFun(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Ep;
        gsVector<T> value(2),init(2);
        gsVector<T,3> n1_vec;
        T n1, n2;

        index_t iter;

        T theta;
        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);

            init.setZero();
            value.setZero();

            objectiveNL<T> obj(m_mm,E);
            iter = obj.newtonRaphson(value,init,false,1e-12,1000);

            theta = init.at(0);
            gamma = init.at(1);

            n1 = math::cos(theta);
            n2 = math::sin(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;

            Ep = E + gamma * n1_vec;
            result.col(k) = Ep;
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
};

template <typename T>
class EwFun : public gsFunction<T>
{
public:
    EwFun(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());
        gsMatrix<T> E, Ew;
        gsVector<T> value(2),init(2);
        gsVector<T,3> n1_vec;
        T n1, n2;

        index_t iter;

        T theta;
        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);

            init.setZero();
            value.setZero();

            objectiveNL<T> obj(m_mm,E);
            iter = obj.newtonRaphson(value,init,false,1e-12,1000);

            theta = init.at(0);
            gamma = init.at(1);

            gsDebugVar(theta);
            gsDebugVar(gamma);

            n1 = math::cos(theta);
            n2 = math::sin(theta);

            n1_vec<<n1*n1, n2*n2, 2*n1*n2;

            Ew = gamma * n1_vec;
            result.col(k) = Ew;
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
};


template <typename T>
class gammaFun : public gsFunction<T>
{
public:
    gammaFun(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E;
        gsVector<T> value(2),init(2);

        index_t iter;

        T gamma;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);

            init.setZero();
            value.setZero();

            objectiveNL<T> obj(m_mm,E);
            iter = obj.newtonRaphson(value,init,false,1e-12,1000);

            gamma = init.at(1);

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

private:
    const gsMaterialMatrixBase<T> * m_mm;
};

template <typename T>
class thetaFun : public gsFunction<T>
{
public:
    thetaFun(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1,u.cols());
        gsMatrix<T> E;
        gsVector<T> value(2),init(2);

        index_t iter;

        T theta;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);

            init.setZero();
            value.setZero();


            objectiveNL<T> obj(m_mm,E);
            iter = obj.newtonRaphson(value,init,false,1e-12,1000);

            theta = init.at(0);

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

private:
    const gsMaterialMatrixBase<T> * m_mm;
};

/// Takes e
template <typename T>
class SwFun : public gsFunction<T>
{
public:
    SwFun(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(3,u.cols());

        EfullFun<T> Efull_(m_mm);
        gsMatrix<T> E, Etmp, S, tmp;

        Efull_.eval_into(u,E);
        for (index_t k = 0; k!=u.cols(); k++)
        {
            Etmp = E.col(k);

            toMatrix(Etmp);

            S = m_mm->S(Etmp);

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
};

/// Takes e
template <typename T>
class SFun : public gsFunction<T>
{
public:
    SFun(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // result.resize(3,u.cols());
        result.resize(4,u.cols());

        EfullFun<T> Efull_(m_mm);
        gsMatrix<T> E, S;
        for (index_t k = 0; k!=u.cols(); k++)
        {
            E = u.col(k);

            S = m_mm->S(E.reshape(2,2));

            // Make voight notation
            S.resize(4,1);

            result.col(k) = S;
        }
    }

    short_t domainDim() const
    {
        return 3;
    }

    short_t targetDim() const
    {
        return 4;
    }

private:
    const gsMaterialMatrixBase<T> * m_mm;
};


/// Takes e
template <typename T>
class vectorS_e : public gsFunction<T>
{
public:
    vectorS_e(const gsMaterialMatrixBase<T> * mm, T gamma, T theta, const gsVector<T> & e)
    :
    m_mm(mm),
    m_gamma(gamma),
    m_theta(theta),
    m_e(e)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        T n1 = math::cos(m_theta);
        T n2 = math::sin(m_theta);

        gsVector<T> n1_vec(3);
        n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T> e_w = m_gamma*n1_vec;

        gsMatrix<T> E, S;
        result.resize(3,u.cols());

        for (index_t k = 0; k!=u.cols(); k++)
        {
            // E = u.col(k) + e_w;
            E = m_e+e_w;
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
    const T m_gamma, m_theta;
    const gsVector<T> m_e;
};

// takes e_w
template <typename T>
class vectorS_ew : public gsFunction<T>
{
public:
    vectorS_ew(const gsMaterialMatrixBase<T> * mm)
    :
    m_mm(mm)
    {

    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> E, S;
        result.resize(3,u.cols());
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
};



/**
 * @brief
 * Simple optimization example, to demonstrate the definition of an
 * optimization problem using the base class gsOptProblem.
 *
 *  This class implements the following NLP.
 *
 * min_x f(x) = (... )^2     (objetive function) --> square minimize it to zero
 *  s.t.
 *       beta > 0
 *       t sigma t > 0
 *       0.5 * pi <= theta <= 0.5 * pi
 *
 */
#ifdef GISMO_WITH_IPOPT

template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
{
public:

    gsOptProblemExample(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    {
        m_numDesignVars  = 1;
        m_numConstraints = 1;
        m_numConJacNonZero = 1;

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // theta has a lower bound of -1 and an upper bound of 1
    const T inf = std::numeric_limits<T>::infinity();
        m_desLowerBounds[0] = -0.5*M_PI;
        m_desUpperBounds[0] =  0.5*M_PI;

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // we have one equality constraint, so we set the bounds on
        // this constraint to be equal (and zero).
        m_conLowerBounds[0] = 0;
        m_conUpperBounds[0] = inf;

        // m_conLowerBounds[0] = -1;
        // m_conUpperBounds[0] = 1;

        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign(0,0) = 0.0;

        //
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);
    }


public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        T theta = u(0,0);
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);

        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );

        gsMatrix<T,1,1> result = n2_vec.transpose() * m_C * m_e + gamma * n2_vec.transpose() * m_C * n1_vec;
        GISMO_ASSERT(result.rows()==1 && result.cols() ==1,"f is not scalar!");

        return result(0,0) * result(0,0);
    }

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        T theta = u(0,0);
        result.resize(m_numConstraints,1);
        // return the value of the constraints: g(x)
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);

        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> gamma = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );
        result(0) = gamma(0,0);
    }

    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // at the moment only a full finite difference matrix is returned.

        gsVector<T> uu = u;
        gsVector<T> e1( this->m_numConstraints );
        gsVector<T> e2( this->m_numConstraints );
        gsAsVector<T> ee1( e1.data() , e1.rows() );
        gsAsVector<T> ee2( e2.data() , e2.rows() );

        index_t lastDesginVar = -1;

        // TODO: Replace by a better value or use AD...
        const T h = T(0.00001);

        for( index_t i = 0 ; i < this->m_numConJacNonZero; ++i )
        {
            index_t row = this->m_conJacRows[i];  // constrains
            index_t col = this->m_conJacCols[i];  // designVariables

            if( lastDesginVar != col )
            {
                gsAsConstVector<T> uuMap( uu.data() , uu.rows() );

                uu(col) -= h/2.;
                evalCon_into( uuMap, ee1);
                uu(col) += h;
                evalCon_into( uuMap, ee2 );
                uu(col) = u(col);

                lastDesginVar = col;
            }

            result(i) = (0.5*e1(row) - 0.5*e2(row)) / h;

        }

    }

private:

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;

    const gsMatrix<T> m_C;
    const gsMatrix<T> m_e;
};
#endif

template <typename T>
class gsTFTMat
{
public:
    gsTFTMat(const gsMatrix<T> & C, const gsMatrix<T> & e)
    :
    m_C(C),
    m_e(e)
    { }

    void compute(T theta)
    {
        T n1 = math::cos(theta);
        T n2 = math::sin(theta);
        T m1 = -math::sin(theta);
        T m2 = math::cos(theta);
        gsVector<T,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<T,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<T,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<T,1,1> denum = n1_vec.transpose() * m_C * n1_vec;

        C_I = m_C - 1 / (  n1_vec.transpose() * m_C * n1_vec ) * m_C * ( n1_vec * n1_vec.transpose() ) * m_C;


        gsMatrix<T> gamma_tmp = - ( n1_vec.transpose() * m_C * m_e ) / ( n1_vec.transpose() * m_C * n1_vec );
        gamma = gamma_tmp(0,0);

        gsMatrix<T,1,1> tmp2 = (n1_vec.transpose() * m_C * n2_vec);

        GISMO_ASSERT(tmp2.rows()==1 && tmp2.cols()==1,"Must be scalar");

        T df =  n4_vec.transpose() * m_C * (m_e + gamma * n1_vec)
                + 2 * gamma * ( n2_vec.transpose() * m_C * n2_vec
                - math::pow(tmp2(0,0),2) / (n1_vec.transpose() * m_C * n1_vec)
                )
                                ;


        gsMatrix<T,3,1> b = n2_vec - ( (n1_vec.transpose() * m_C * n2_vec)(0,0) / ( n1_vec.transpose() * m_C * n1_vec )(0,0)) * n1_vec;


        C_II = C_I + 2 * gamma / df * (m_C * b * b.transpose() * m_C);

        gsMatrix<T> Sp = m_C * (m_e + gamma * n1_vec);
        gsDebugVar(Sp);
        gsDebugVar(theta);
        gsDebugVar(gamma);
        gsDebugVar(n1_vec);
        gsDebugVar(n2_vec);
        gsDebugVar(n1_vec.transpose() * Sp);
        gsDebugVar(n2_vec.transpose() * Sp);

    }

public:
    mutable gsMatrix<T> C_I, C_II;
    mutable T gamma;

private:
    const gsMatrix<T> m_C;
    const gsVector<T> m_e;
};


template<class T>
T findMod(T a, T b)
{
    T mod;
    // Handling negative values
    if (a < 0)
        mod = -a;
    else
        mod =  a;
    if (b < 0)
        b = -b;

    // Finding mod by repeated subtraction

    while (mod >= b)
        mod = mod - b;

    // Sign of result typically depends
    // on sign of a.
    if (a < 0)
        return -mod;

    return mod;
}

template <class T>
class gsScalarRootFinder
{
protected:
    T m_x;
    T m_xmin;
    T m_xmax;
    T m_error;
    T m_tolerance;

    index_t m_maxIterations;
    index_t m_iteration;

    typedef std::function < T ( T const &) > function_t;
    function_t m_function;

public:
    gsScalarRootFinder(T xmin, T xmax, function_t fun)
    :
    m_xmin(xmin),
    m_xmax(xmax),
    m_function(fun),
    m_iteration(0),
    m_maxIterations(15),
    m_error(1.0),
    m_tolerance(1e-15)
    {

    }

    void compute()
    {
        // From: https://lemesurierb.github.io/elementary-numerical-analysis-python/notebooks/root-finding-without-derivatives-python.html
        T x_older = m_xmin;
        T x_more_recent = m_xmax;
        T x_new;
        T f_x_older = m_function(x_older);
        T f_x_more_recent = m_function(x_more_recent);
        T f_x_new;
        for (m_iteration = 0; m_iteration <= m_maxIterations; m_iteration++)
        {
            x_new = (x_older * f_x_more_recent - f_x_older * x_more_recent)/(f_x_more_recent - f_x_older);
            x_new = findMod(x_new , M_PI);

            f_x_new = m_function(x_new);

            gsDebugVar(x_new);
            gsDebugVar(x_more_recent);

            x_older = x_more_recent;
            x_more_recent = x_new;
            f_x_older = f_x_more_recent;
            f_x_more_recent = f_x_new;

            m_error = std::abs(x_older - x_more_recent);

            gsDebug<<m_iteration<<"\t"<<m_error<<"\n";

            if (m_error < m_tolerance)
                break;
        }
        if (m_iteration < m_maxIterations)
        {
            m_x = x_new;
        }
        else
        {
            gsWarn<<"Did not converge in "<<m_iteration<<"iterations\n";
            m_x = NAN;
        }


        // T a = m_xmin;
        // T b = m_xmax;
        // T fa = m_function(a);
        // T fb = m_function(b);
        // T c;
        // for (m_iteration = 0; m_iteration!= m_maxIterations; m_iteration++)
        // {
        //         c = (a * fb - fa * b)/(fb - fa);
        //         T fc = m_function(c);
        //         if (fa * fc < 0)
        //         {
        //             b = c;
        //             fb = fc;
        //         }
        //         else
        //         {
        //             a = c;
        //             fa = fc;
        //         }

        //         m_error = std::abs(b-a);

        //         gsDebug<<m_iteration<<"\t"<<m_error<<"\n";

        //     if (m_error < m_tolerance)
        //         break;
        // }
        // if (m_iteration < m_maxIterations)
        // {
        //     m_x = c;
        // }
        // else
        // {
        //     gsWarn<<"Did not converge in "<<m_iteration<<"iterations\n";
        //     m_x = NAN;
        // }

    }

    T result() { return m_x; }
    T error()  { return m_error; }
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

      L = 1;
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

    gsField<> tensionField(mp_def, tensionfield, true);
    gsInfo<<"Plotting in Paraview...\n";
    gsWriteParaview<>( tensionField, "TensionField", 1000, true);

    gsVector<> pt(2);

    pt.setConstant(0.25);


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

    gsWarn<<"This should not be pstrain.\n";
    pstrain.eval_into(pt,result);
    gsVector<real_t,3> e = result; // THIS IS ACTUALLY S!!
    gsDebugVar(e);

    pdir.eval_into(pt,result);
    gsMatrix<> dir = result;
    gsDebugVar(dir.reshape(3,3));

    typedef std::function < real_t ( real_t const &) > function_t;
    function_t f_fun = [&C,&e](real_t const & theta)
    {
        real_t n1 = math::cos(theta);
        real_t n2 = math::sin(theta);
        real_t m1 = -math::sin(theta);
        real_t m2 = math::cos(theta);

        gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<real_t,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * C * e ) / ( n1_vec.transpose() * C * n1_vec );
        gsMatrix<real_t,1,1> result = n2_vec.transpose() * C * e + gamma * n2_vec.transpose() * C * n1_vec;
        GISMO_ASSERT(result.rows()==1 && result.cols() ==1,"f is not scalar!");
        return result(0,0);
    };

    typedef std::function < real_t ( real_t const &) > function_t;
    function_t gamma_fun = [&C,&e](real_t const & theta)
    {
        real_t n1 = math::cos(theta);
        real_t n2 = math::sin(theta);
        real_t m1 = -math::sin(theta);
        real_t m2 = math::cos(theta);

        gsVector<real_t,3> n1_vec; n1_vec<<n1*n1, n2*n2, 2*n1*n2;
        gsVector<real_t,3> n2_vec; n2_vec<<m1*n1, m2*n2, m1*n2+m2*n1;
        gsVector<real_t,3> n3_vec; n3_vec<<m1*m1-n1*n1, m2*m2-n2*n2, 2*(m1*m2-n1*n2);
        gsVector<real_t,3> n4_vec; n4_vec<<m1*m1, m2*m2, 2*m1*m2;

        gsMatrix<real_t,1,1> gamma = - ( n1_vec.transpose() * C * e ) / ( n1_vec.transpose() * C * n1_vec );
        return gamma(0,0);
    };

    objective<real_t> obj(C,e);

    gsVector<> value(1);
    value<<0;
    gsVector<> init(1);
    init<<0;

    gsVector<> e4(4);
    e4<<e(0),0.5*e(2),0.5*e(2),e(1);

    gsTFTMat<real_t> matCompute(C,e);

    index_t iter;

    iter = obj.newtonRaphson(value,init,false,1e-12,1000);
    gsDebugVar(iter);

    matCompute.compute(init(0));
    gsDebugVar(matCompute.C_I);
    gsDebugVar(matCompute.C_II);
    gsDebugVar(matCompute.gamma);

    gsDebugVar(init);
    gsDebugVar(obj.eval(init));


    gammaFun<real_t> gamma(materialMatrix);
    gamma.eval_into(e,result);
    gsDebugVar(result);

    thetaFun<real_t> theta(materialMatrix);
    theta.eval_into(e,result);
    gsDebugVar(result);


    EwFun<real_t> Ew(materialMatrix);
    Ew.eval_into(e,result);
    gsDebugVar(result);


    EfullFun<real_t> Efull(materialMatrix);
    Efull.eval_into(e,result);
    gsDebugVar(result);

    SwFun<real_t> Sw(materialMatrix);
    Sw.eval_into(e,result);
    gsDebugVar(result);

    Sw.deriv_into(e,result);
    gsDebugVar(result.reshape(3,3));
    gsDebugVar(result.reshape(3,3) * e);

    gsDebugVar(matCompute.C_I );
    gsDebugVar(matCompute.C_II);

    gsDebugVar(matCompute.C_I  * e);
    gsDebugVar(matCompute.C_II * e);

    SFun<real_t> S(materialMatrix);
    S.eval_into(e4,result);
    gsDebugVar(result);

    S.deriv_into(e4,result);
    gsDebugVar(result.reshape(4,4));
    gsDebugVar(result.reshape(4,4) * e4);
    // gsDebugVar(result.reshape(3,3) * e);

    gsDebugVar(C );
    gsDebugVar(C * e );


    // gsDebugVar(Cnew * e);

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
