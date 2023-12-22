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
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsDofMapper.h>

#include <gsPde/gsBoundaryConditions.h>

#include <gsUtils/gsPointGrid.h>

#include <gsAssembler/gsAssemblerOptions.h>
#include <gsAssembler/gsExpressions.h>


#  define MatExprType  auto

namespace gismo{
namespace expr{

/**
 * @brief      Simple expression for the unit vector of length \a dim and with value 1 on \a index
 */
class unitVec_expr : public _expr<unitVec_expr >
{
public:
    typedef real_t Scalar;
private:
    index_t _index;
    index_t _dim;

public:
    unitVec_expr(const index_t index, const index_t dim) : _index(index), _dim(dim) { }

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    gsMatrix<Scalar> eval(const index_t) const
    {
        gsMatrix<Scalar> vec = gsMatrix<Scalar>::Zero(_dim,1);
        // vec.setZero();
        vec(_index,0) = 1;
        return vec;
    }

    index_t rows() const { return _dim; }
    index_t cols() const { return  1; }
    void parse(gsExprHelper<Scalar> & evList) const
    {  }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "uv("<<_dim <<")";}
};


/**
 * @brief      Expression for the first variation of the surface normal
 *
 * @tparam     E     Object type
 */
template<class E>
class var1_expr : public _expr<var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

    mutable gsMatrix<Scalar> res;
    mutable gsMatrix<Scalar> bGrads, cJac;
    mutable gsVector<Scalar,3> m_v, normal;
public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    var1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return cols_impl(_u); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        parse_impl<E>(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_ACTIVE | NEED_GRAD; // need actives for cardinality
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_DERIV;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar>>::value || util::is_same<U,gsFeVariable<Scalar>>::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_DERIV;

        grad(_u).parse(evList); //

        _u.parse(evList);
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        const index_t A = _u.cardinality()/_u.dim(); // _u.data().actives.rows()
        res.resize(_u.cardinality(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();

        bGrads = _u.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure = (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // Jac(u) ~ Jac(G) with alternating signs ?..
                m_v.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col3d(1) )
                              - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col3d(0) )) / measure;

                // ---------------  First variation of the normal
                res.row(s+j).noalias() = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
            }
        }
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        grad_expr<U> vGrad = grad_expr<U>(u);
        bGrads = vGrad.eval(k);

        return make_vector(k);
    }

    // Specialisation for a solution
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        grad_expr<gsFeSolution<Scalar>> sGrad =  grad_expr<gsFeSolution<Scalar>>(_u);
        bGrads = sGrad.eval(k);

        return make_vector(k);
    }

    const gsMatrix<Scalar> & make_vector(const index_t k) const
    {
        res.resize(rows(), cols()); // rows()*
        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();

        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure = (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        m_v.noalias() = ( ( bGrads.col3d(0) ).cross( cJac.col3d(1) )
                      -   ( bGrads.col3d(1) ).cross( cJac.col3d(0) ) ) / measure;

        // ---------------  First variation of the normal
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, index_t >::type
    cols_impl(const U & u)  const
    {
        return _u.data().dim.second;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value || util::is_same<U,gsFeSolution<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.dim();
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeVariable<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.source().targetDim();
    }

};

/**
 * @brief      Second variation of the surface normal times a vector.
 *
 * @tparam     E1    Type of u
 * @tparam     E2    Type of v
 * @tparam     E3    Type of the vector
 */
template<class E1, class E2, class E3>
class var2dot_expr : public _expr<var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    typename E3::Nested_t _Ef;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0 };

    var2dot_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& Ef) : _u(u),_v(v), _G(G), _Ef(Ef) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2, evEf, result;
    mutable gsVector<Scalar> m_u, m_v, normal, m_uv, m_u_der, n_der, n_der2, tmp; // memomry leaks when gsVector<T,3>, i.e. fixed dimension
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        res.resize(_u.cardinality(), _u.cardinality());

        normal = _G.data().normal(k);
        normal.normalize();
        uGrads = _u.data().values[1].col(k);
        vGrads = _v.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();

        const index_t cardU = _u.data().values[0].rows(); // number of actives per component of u
        const index_t cardV = _v.data().values[0].rows(); // number of actives per component of v
        const Scalar measure = (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        evEf = _Ef.eval(k);

        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    m_u.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col3d(1) )
                                     -vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col3d(0) ))
                                    / measure;

                    const short_t s = d*cardU;

                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;
                        m_v.noalias() = ( vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col3d(1) )
                                         -vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col3d(0) ))
                                        / measure;

                        // n_der.noalias() = (m_v - ( normal.dot(m_v) ) * normal);
                        n_der.noalias() = (m_v - ( normal*m_v.transpose() ) * normal); // outer-product version

                        m_uv.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                          +vecFun(c, vGrads.at(2*i  ) ).cross( vecFun(d, uGrads.at(2*j+1) ) ))
                                          / measure; //check

                        m_u_der.noalias() = (m_uv - ( normal.dot(m_v) ) * m_u);
                        // m_u_der.noalias() = (m_uv - ( normal*m_v.transpose() ) * m_u); // outer-product version TODO

                        // ---------------  Second variation of the normal
                        tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;
                        // tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;

                        // Evaluate the product
                        result = evEf * tmp;

                        res(s + j, r + i ) = result(0,0);
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const
    {
        return 1; // because the resulting matrix has scalar entries for every combination of active basis functions
    }

    index_t cols() const
    {
        return 1; // because the resulting matrix has scalar entries for every combination of active basis functions
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD;
        evList.add(_v);
        _v.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD;
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_DERIV;
        _Ef.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const { os << "var2("; _u.print(os); os <<")"; }
};

/**
 * @brief      Second variation of the surface normal times the second derivative of the geometry map times a vector.
 *
 * @tparam     E1    Type of u
 * @tparam     E2    Type of v
 * @tparam     E3    Type of the vector
 */
template<class E1, class E2, class E3>
class var2deriv2dot_expr : public _expr<var2deriv2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    typename E3::Nested_t _Ef;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0 };

    var2deriv2dot_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& Ef) : _u(u),_v(v), _G(G), _Ef(Ef) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2, evEf, result;
    mutable gsVector<Scalar> m_u, m_v, normal, m_uv, m_u_der, n_der, n_der2, tmp; // memomry leaks when gsVector<T,3>, i.e. fixed dimension
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        res.resize(_u.cardinality(), _u.cardinality());

        normal = _G.data().normal(k);
        normal.normalize();
        uGrads = _u.data().values[1].col(k);
        vGrads = _v.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        cDer2 = _G.data().values[2].reshapeCol(k, _G.data().dim.second, _G.data().dim.second);

        const index_t cardU = _u.data().values[0].rows(); // number of actives per component of u
        const index_t cardV = _v.data().values[0].rows(); // number of actives per component of v
        const Scalar measure = (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        evEf = _Ef.eval(k);

        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    m_u.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col3d(1) )
                                     -vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col3d(0) ))
                                    / measure;

                    const short_t s = d*cardU;

                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;
                        m_v.noalias() = ( vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col3d(1) )
                                         -vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col3d(0) ))
                                        / measure;

                        // n_der.noalias() = (m_v - ( normal.dot(m_v) ) * normal);
                        n_der.noalias() = (m_v - ( normal*m_v.transpose() ) * normal); // outer-product version

                        m_uv.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                          +vecFun(c, vGrads.at(2*i  ) ).cross( vecFun(d, uGrads.at(2*j+1) ) ))
                                          / measure; //check

                        m_u_der.noalias() = (m_uv - ( normal.dot(m_v) ) * m_u);
                        // m_u_der.noalias() = (m_uv - ( normal*m_v.transpose() ) * m_u); // outer-product version TODO

                        // ---------------  Second variation of the normal
                        tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;
                        // tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;

                        // Evaluate the product
                        tmp = cDer2 * tmp; // E_f_der2, last component
                        tmp.row(2) *= 2.0;
                        result = evEf * tmp;

                        res(s + j, r + i ) = result(0,0);
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const
    {
        return 1; // because the resulting matrix has scalar entries for every combination of active basis functions
    }

    index_t cols() const
    {
        return 1; // because the resulting matrix has scalar entries for every combination of active basis functions
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD;
        evList.add(_v);
        _v.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD;
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_2ND_DER;
        _Ef.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const { os << "var2dot("; _u.print(os), _v.print(os), _G.print(os), _Ef.print(os); os <<")"; }
};

/// Second variation of the normal
template<class E1, class E2>
class var2_expr : public _expr<var2_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E1::Space, ScalarValued= 0, ColBlocks= 0 };

    var2_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G) : _u(u),_v(v), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2, result;
    mutable gsVector<Scalar> m_u, m_v, normal, m_uv, m_u_der, n_der, n_der2, tmp; // memomry leaks when gsVector<T,3>, i.e. fixed dimension

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,_v,k); }

    index_t rows() const
    {
        return 3;
    }

    index_t cols() const
    {
        return 1;
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        grad(_u).parse(evList); //
        grad(_v).parse(evList); //

        _u.parse(evList);
        evList.add(_u);
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;

        _v.parse(evList);
        evList.add(_v);
        _v.data().flags |= NEED_GRAD | NEED_ACTIVE;

        _G.parse(evList);
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_DERIV2 ;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    enum{rowSpan = 0, colSpan = 0};

    void print(std::ostream &os) const { os << "var2("; _u.print(os), _v.print(os), _G.print(os); os <<")"; }

private:
// Specialisation for a space
    template<class U, class V> inline
    typename util::enable_if< !(util::is_same<U,gsFeVariable<Scalar> >::value && util::is_same<V,gsFeVariable<Scalar> >::value), const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const V & v, const index_t k)  const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<class U, class V> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value && util::is_same<V,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const V & v, const index_t k)  const
    {
        grad_expr<U> uGrad = grad_expr<U>(u);
        uGrads = uGrad.eval(k);
        grad_expr<V> vGrad = grad_expr<V>(v);
        vGrads = vGrad.eval(k);

        res.resize(rows(), cols()); // rows()*

        normal = _G.data().normal(k);
        normal.normalize();

        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        cDer2 = _G.data().values[2].reshapeCol(k, _G.data().dim.second, _G.data().dim.second);

        const Scalar measure = (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        m_u.noalias() = ( ( uGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                        - ( uGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) )
                        / measure;

        m_v.noalias() = ( ( vGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                        - ( vGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) )
                        / measure;

        // n_der.noalias() = (m_v - ( normal.dot(m_v) ) * normal);
        n_der.noalias() = (m_v - ( normal*m_v.transpose() ) * normal); // outer-product version

        m_uv.noalias() = ( ( uGrads.col(0).template head<3>() ).cross( vGrads.col(1).template head<3>() )
                         - ( vGrads.col(1).template head<3>() ).cross( uGrads.col(0).template head<3>() ) )
                        / measure;

        m_u_der.noalias() = (m_uv - ( normal.dot(m_v) ) * m_u);
        // m_u_der.noalias() = (m_uv - ( normal*m_v.transpose() ) * m_u); // outer-product version TODO

        // ---------------  Second variation of the normal
        tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;
        // tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;

        // Evaluate the product
        tmp = cDer2 * tmp; // E_f_der2, last component
        tmp.row(2) *= 2.0;
        res = tmp;
        return res;
    }

    // // Specialisation for a solution
    // template<class U> inline
    // typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    // eval_impl(const U & u, const index_t k)  const
    // {
    //     GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
    //     solGrad_expr<Scalar> sGrad =  solGrad_expr<Scalar>(_u);
    //     res.resize(rows(), cols()); // rows()*

    //     normal = _G.data().normal(k);// not normalized to unit length
    //     normal.normalize();
    //     bGrads = sGrad.eval(k);
    //     cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
    //     const Scalar measure =  _G.data().measures.at(k);

    //     m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
    //                   -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

    //     // ---------------  First variation of the normal
    //     // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
    //     res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
    //     return res;
    // }

};


/**
 * @brief      Expression for the first variation of the outer tangent
 *
 * @tparam     E     Object type
 */
template<class E>
class tvar1_expr : public _expr<tvar1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

    mutable gsVector<Scalar,3> onormal, tangent, utangent, dtan;
    mutable gsVector<Scalar> tmp;
    mutable gsMatrix<Scalar> bGrads, cJac, res;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    tvar1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

#   define Eigen gsEigen
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }
    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return _u.dim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        parse_impl<E>(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "tvar("; _G.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV;

        grad(_u).parse(evList); //

        _u.parse(evList);
    }

    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);

        const index_t A = _u.cardinality()/_u.dim();
        res.resize(_u.cardinality(), cols()); // rows()*
        res.setZero();
        cJac = _G.data().jacobian(k);
        cJac.colwise().normalize();

        onormal = _G.data().outNormal(k);
        onormal.normalize();
        tmp = cJac.transpose() * onormal;
        tmp.colwise().normalize(); //normalize the inner product for fair comparison

        Scalar tol = 1e-8;
        /*
            We can check which column of the Jacobian corresponds to the outer normal vector or to the tangent.
            The tangent is a covariant vector and hence the column of the Jacobian should be equal to the tangent.
            The normal is a contravariant vector and hence the corresponding column of the Jacobian times the outward normal should give 1. We use this property.
        */
        index_t colIndex = -1;
        if      ( (math::abs(tmp.at(0)) < tol) && (math::abs(tmp.at(1)) > 1-tol ) )         // then the normal is vector 2 and the tangent vector 1
            colIndex = 0;
        else if ( (math::abs(tmp.at(1)) < tol) && (math::abs(tmp.at(0)) > 1-tol ) )     // then the normal is vector 1 and the tangent vector 2
            colIndex = 1;
        else                    // then the normal is unknown??
            gsInfo<<"warning: choice unknown\n";

        // tangent = cJac.col(colIndex);
        tangent_expr<Scalar> tan_expr = tangent_expr<Scalar>(_G);
        tangent = tan_expr.eval(k);
        utangent = tangent.normalized();

        index_t sign = 0;
        if (colIndex!=-1)
        {
            Scalar dot = tangent.dot(cJac.col(colIndex));
            sign = (Scalar(0) < dot) - (dot < Scalar(0));
        }
        else
        {
            gsWarn<<"No suitable tangent and outer normal vector found for point "<<_G.data().values[0].transpose()<<"\n";
            return res;
        }

        // Now we will compute the derivatives of the basis functions
        bGrads = _u.data().values[1].col(k);
        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // The tangent vector is in column colIndex in cJac and thus in 2*j+colIndex in bGrads.
                // Furthermore, as basis function for dimension d, it has a nonzero in entry d, and zeros elsewhere
                dtan = sign*vecFun(d, bGrads.at(2*j+colIndex));
                res.row(s+j).noalias() = (1 / tangent.norm() * ( dtan - ( utangent.transpose() * dtan ) * utangent ) ).transpose();
            }
        }
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        res.resize(rows(), cols());
        res.setZero();

        cJac = _G.data().jacobian(k);
        cJac.colwise().normalize();

        onormal = _G.data().outNormal(k);
        onormal.normalize();
        tmp = cJac.transpose() * onormal;
        tmp.colwise().normalize(); //normalize the inner product for fair comparison

        Scalar tol = 1e-8;

        /*
            We can check which column of the Jacobian corresponds to the outer normal vector or to the tangent.
            The tangent is a covariant vector and hence the column of the Jacobian should be equal to the tangent.
            The normal is a contravariant vector and hence the corresponding column of the Jacobian times the outward normal should give 1. We use this property.
        */
        index_t colIndex = -1;
        if      ( (math::abs(tmp.at(0)) < tol) && (math::abs(tmp.at(1)) > 1-tol ) )         // then the normal is vector 2 and the tangent vector 1
            colIndex = 0;
        else if ( (math::abs(tmp.at(1)) < tol) && (math::abs(tmp.at(0)) > 1-tol ) )     // then the normal is vector 1 and the tangent vector 2
            colIndex = 1;
        else                    // then the normal is unknown??
            gsInfo<<"warning: choice unknown\n";

        // tangent = cJac.col(colIndex);
        tangent_expr<Scalar> tan_expr = tangent_expr<Scalar>(_G);
        tangent = tan_expr.eval(k);
        // utangent = tangent.normalized();

        index_t sign = 0;
        if (colIndex!=-1)
        {
            Scalar dot = tangent.dot(cJac.col(colIndex));
            sign = (Scalar(0) < dot) - (dot < Scalar(0));
        }
        else
        {
            gsWarn<<"No suitable tangent and outer normal vector found for point "<<_G.data().values[0].transpose()<<"\n";
            return res;
        }

        bGrads = _u.data().values[1].col(k);
        dtan = sign*bGrads.col(colIndex);
        res.noalias() = (1 / tangent.norm() * ( dtan - ( tangent * dtan ) * tangent / (tangent.norm() * tangent.norm()) )).transpose();
        return res;
    }

};

/**
 * @brief      Expression for the first variation of the outer normal
 *
 * @tparam     E     Object type
 */
template<class E>
class ovar1_expr : public _expr<ovar1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

    mutable gsVector<Scalar,3> tangent, normal, tvar, snvar;

    mutable gsMatrix<Scalar> tvarMat, snvarMat, res;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    ovar1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const { return eval_impl(_u,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return _u.dim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        parse_impl<E>(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "ovar("; _G.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV; // all needed?

        tv(_G).parse(evList);
        tvar1(_u,_G).parse(evList);
        var1(_u,_G).parse(evList);
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV; // all needed?

        grad(_u).parse(evList); //
        tv(_G).parse(evList);

        _u.parse(evList);
    }

    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);
        const index_t A = _u.cardinality()/_u.dim(); // _u.data().actives.rows()
        res.resize(_u.cardinality(), cols()); // rows()*

        tangent_expr<Scalar> tan_expr = tangent_expr<Scalar>(_G);
        tangent = tan_expr.eval(k);
        tangent.normalize();

        // For the normal vector variation
        normal  =  _G.data().normal(k);
        normal.normalize();

        tvar1_expr<E> tvar_expr = tvar1_expr<E>(_u,_G);
        tvarMat = tvar_expr.eval(k);

        var1_expr<E> snvar_expr = var1_expr<E>(_u,_G);
        snvarMat = snvar_expr.eval(k);

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                tvar = tvarMat.row(s+j);
                snvar= snvarMat.row(s+j);

                // VARIATION OF THE OUTER NORMAL
                res.row(s+j).noalias() = tvar.cross(normal) + tangent.cross(snvar);
            }
        }
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);
        res.resize(rows(), cols());

        tangent_expr<Scalar> tan_expr = tangent_expr<Scalar>(_G);
        tangent = tan_expr.eval(k);
        tangent.normalize();

        // For the normal vector variation
        normal  =  _G.data().normal(k);
        normal.normalize();

        tvar1_expr<E> tvar_expr = tvar1_expr<E>(_u,_G);
        tvar = tvar_expr.eval(k);

        var1_expr<E> snvar_expr = var1_expr<E>(_u,_G);
        snvar = snvar_expr.eval(k);


        // VARIATION OF THE OUTER NORMAL
        res.noalias() = tvar.cross(normal) + tangent.cross(snvar);

        return res;
    }
};

/**
 * @brief      Expression for the second variation of the outer normal times a vector
 *
 * @tparam     E1    Type of u
 * @tparam     E2    Type of v
 * @tparam     E3    Type of the vector
 */
template<class E1, class E2, class E3>
class ovar2dot_expr : public _expr<ovar2dot_expr<E1,E2,E3> >
    {
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    typename E3::Nested_t _C;

    mutable gsVector<Scalar,3> normal, onormal, tangent, utangent, dtanu, dtanv, tvaru, tvarv, tvar2,
                            mu, mv, muv, mu_der, snvaru, snvarv, snvar2, nvar2;
    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, res, eC;

    mutable gsMatrix<Scalar> tvaruMat, tvarvMat, tmp;
public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0 };

    ovar2dot_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& C) : _u(u), _v(v), _G(G), _C(C) { }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);
        GISMO_ASSERT(_C.cols()*_C.rows()==3, "Size of vector is incorrect");
        res.resize(_u.cardinality(), _u.cardinality());
        res.setZero();

        eC = _C.eval(k);
        eC.resize(1,3);

        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        cJac.colwise().normalize();

        onormal = _G.data().outNormal(k);
        onormal.normalize();
        tmp = cJac.transpose() * onormal;
        tmp.colwise().normalize(); //normalize the inner product for fair comparison

        Scalar tol = 1e-8;

        /*
            We can check which column of the Jacobian corresponds to the outer normal vector or to the tangent.
            The tangent is a covariant vector and hence the column of the Jacobian should be equal to the tangent.
            The normal is a contravariant vector and hence the corresponding column of the Jacobian times the outward normal should give 1. We use this property.
        */
        index_t colIndex = -1;
        if ( (math::abs(tmp.at(0)) < tol) && (math::abs(tmp.at(1)) > 1-tol ) )         // then the normal is vector 2 and the tangent vector 1
            colIndex = 0;
        else if ( (math::abs(tmp.at(1)) < tol) && (math::abs(tmp.at(0)) > 1-tol ) )     // then the normal is vector 1 and the tangent vector 2
            colIndex = 1;
        else                    // then the normal is unknown??
            gsInfo<<"warning: choice unknown\n";

        // tangent = cJac.col(colIndex);
        // utangent = tangent / tangent.norm();

        // Required for the normal vector variation
        normal = _G.data().normal(k);
        normal.normalize();
        uGrads = _u.data().values[1].col(k);
        vGrads = _v.data().values[1].col(k);

        const index_t cardU = _u.data().values[0].rows(); // number of actives per component of u
        const index_t cardV = _v.data().values[0].rows(); // number of actives per component of v
        const Scalar measure = (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        tangent_expr<Scalar> tan_expr = tangent_expr<Scalar>(_G);
        tangent = tan_expr.eval(k);
        utangent = tangent.normalized();

        index_t sign = 0;
        if (colIndex!=-1)
        {
            Scalar dot = tangent.dot(cJac.col(colIndex));
            sign = (Scalar(0) < dot) - (dot < Scalar(0));
        }
        else
        {
            gsWarn<<"No suitable tangent and outer normal vector found for point "<<_G.data().values[0].transpose()<<"\n";
            return res;
        }

        // // For the normal vector variation
        // normal  =  _G.data().normal(k);
        // normal.normalize();

        // tvar1_expr<E1> tvaru_expr = tvar1_expr<E1>(_u,_G);
        // tvaruMat = tvaru_expr.eval(k);
        // tvar1_expr<E2> tvarv_expr = tvar1_expr<E2>(_v,_G);
        // tvarvMat = tvarv_expr.eval(k);

        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    const short_t s = d*cardU;

                    // // first variation of the tangent
                    // tvaru = tvaruMat.row(s+j);
                    // first variation of the tangent (colvector)
                    dtanu = sign*vecFun(d, uGrads.at(2*j+colIndex));
                    tvaru = 1 / tangent.norm() * ( dtanu - ( utangent.dot(dtanu) ) * utangent );

                    // first variation of the surface normal (colvector)
                    mu.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col3d(1) )
                                    -vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col3d(0) ))
                                    / measure;
                    snvaru.noalias() = (mu - ( normal.dot(mu) ) * normal);

                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;

                        // // first variation of the tangent
                        // tvarv = tvarvMat.row(r+i);
                        // first variation of the tangent (colvector)
                        dtanv = sign*vecFun(c, vGrads.at(2*i+colIndex));
                        tvarv = 1 / tangent.norm() * ( dtanv - ( utangent.dot(dtanv) ) * utangent );

                        // first variation of the surface normal (colvector)
                        mv.noalias() = ( vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col3d(1) )
                                        -vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col3d(0) ))
                                        / measure;
                        snvarv.noalias() = (mv - ( normal.dot(mv) ) * normal);


                        // See point 2 of
                        // Herrema, Austin J. et al. 2021. “Corrigendum to ‘Penalty Coupling of Non-Matching Isogeometric Kirchhoff–Love Shell Patches with Application to Composite Wind Turbine Blades’ [Comput. Methods Appl. Mech. Engrg. 346 (2019) 810–840].” Computer Methods in Applied Mechanics and Engineering 373: 113488.
                        // Second variation of the tangent (colvector)
                        tvar2 = -1 / tangent.norm() * (
                                                        ( tvarv.dot(dtanu) * utangent )
                                                      + ( utangent.dot(dtanv) * tvaru )
                                                      + ( utangent.dot(dtanu) * tvarv )
                                                                            );

                        // Second variation of the surface normal (colvector)
                        muv.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                         +vecFun(c, vGrads.at(2*i  ) ).cross( vecFun(d, uGrads.at(2*j+1) ) ))
                                          / measure;

                        mu_der.noalias() = (muv - ( normal.dot(mv) ) * mu);

                        snvar2 = mu_der - (mu.dot(snvarv) + normal.dot(mu_der) ) * normal - (normal.dot(mu) ) * snvarv;

                        // Second variation of the outer normal (colvector)
                        nvar2 = tvar2.cross(normal) + tvaru.cross(snvarv) + tvarv.cross(snvaru) + utangent.cross(snvar2);

                        res(s + j, r + i ) = (eC * nvar2).value();
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD;
        evList.add(_v);
        _v.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD;
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV;
        _C.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const { os << "nvar2("; _u.print(os); os <<")"; }
};

/**
 * @brief      Expression that takes the second derivative of an expression and multiplies it with a row vector
 *
 * @tparam     E1    Expression
 * @tparam     E2    Row vector
 */
template<class E1, class E2>
class deriv2dot_expr : public _expr<deriv2dot_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum{   Space = (E1::Space == 1 || E2::Space == 1) ? 1 : 0,
            ScalarValued= 0,
            ColBlocks= 0
        };

    typedef typename E1::Scalar Scalar;

    deriv2dot_expr(const E1 & u, const E2 & v) : _u(u), _v(v) { }

    mutable gsMatrix<Scalar> res,tmp, vEv;

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const
    {
        return 1; //since the product with another vector is computed
    }

    index_t cols() const
    {
        return cols_impl(_u);
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        parse_impl<E1>(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const
    {
        if      (E1::Space == 1 && E2::Space == 0)
            return _u.rowVar();
        else if (E1::Space == 0 && E2::Space == 1)
            return _v.rowVar();
        else
            return gsNullExpr<Scalar>::get();
    }

    const gsFeSpace<Scalar> & colVar() const
    {
        if      (E1::Space == 1 && E2::Space == 0)
            return _v.colVar();
        else if (E1::Space == 0 && E2::Space == 1)
            return _u.colVar();
        else
            return gsNullExpr<Scalar>::get();
    }

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); _v.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);
        evList.add(_u);   // We manage the flags of _u "manually" here (sets data)
        _u.data().flags |= NEED_DERIV2; // define flags

        _v.parse(evList); // We need to evaluate _v (_v.eval(.) is called)

        // Note: evList.parse(.) is called only in exprAssembler for the global expression
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);
        evList.add(_u);   // We manage the flags of _u "manually" here (sets data)
        _u.data().flags |= NEED_DERIV2; // define flags

        _v.parse(evList); // We need to evaluate _v (_v.eval(.) is called)

        // Note: evList.parse(.) is called only in exprAssembler for the global expression
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);
        evList.add(_u);   // We manage the flags of _u "manually" here (sets data)
        _u.data().flags |= NEED_DERIV2; // define flags

        _v.parse(evList); // We need to evaluate _v (_v.eval(.) is called)

        // Note: evList.parse(.) is called only in exprAssembler for the global expression
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); //
        hess(_u).parse(evList); //

        // evList.add(_u);   // We manage the flags of _u "manually" here (sets data)
        _v.parse(evList); // We need to evaluate _v (_v.eval(.) is called)
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        /*
            Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
            The hessian of the geometry map c has the form: hess(c)
            [d11 c1, d11 c2, d11 c3]
            [d22 c1, d22 c2, d22 c3]
            [d12 c1, d12 c2, d12 c3]
            And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
            So we simply evaluate for every active basis function v_k the product hess(c).v_k
        */

        // evaluate the geometry map of U
        tmp =_u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second );
        vEv = _v.eval(k);
        res = vEv * tmp.transpose();
        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        /*
            We assume that the basis has the form v*e_i where e_i is the unit vector with 1 on index i and 0 elsewhere
            This implies that hess(v) = [hess(v_1), hess(v_2), hess(v_3)] only has nonzero entries in column i. Hence,
            hess(v) . normal = hess(v_i) * n_i (vector-scalar multiplication. The result is then of the form
            [hess(v_1)*n_1 .., hess(v_2)*n_2 .., hess(v_3)*n_3 ..]. Here, the dots .. represent the active basis functions.
        */
        const index_t numAct = u.data().values[0].rows();   // number of actives of a basis function
        const index_t cardinality = u.cardinality();        // total number of actives (=3*numAct)
        res.resize(rows()*cardinality, cols() );
        tmp.transpose() =_u.data().values[2].reshapeCol(k, cols(), numAct );
        vEv = _v.eval(k);

        for (index_t i = 0; i!=_u.dim(); i++)
            res.block(i*numAct, 0, numAct, cols() ) = tmp * vEv.at(i);
        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        /*
            Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
            The laplacian of the variable has the form: hess(v)
            [d11 c1, d22 c1, d12 c1]
            [d11 c2, d22 c2, d12 c2]
            [d11 c3, d22 c3, d12 c3]
            And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
            So we simply evaluate for every active basis function v_k the product hess(c).v_k
        */

        gsMatrix<> tmp2;
        tmp = u.data().values[2].col(k);
        index_t nDers = _u.source().domainDim() * (_u.source().domainDim() + 1) / 2;
        index_t dim = _u.source().targetDim();
        tmp2.resize(nDers, dim);
        for (index_t comp = 0; comp != u.source().targetDim(); comp++)
            tmp2.col(comp) = tmp.block(comp * nDers, 0, nDers, 1); //star,length

        vEv = _v.eval(k);
        res = vEv * tmp2.transpose();
        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        /*
            Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
            The hessian of the geometry map c has the form: hess(c)
            [d11 c1, d22 c1, d12 c1]
            [d11 c2, d22 c2, d12 c2]
            [d11 c3, d22 c3, d12 c3]
            And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
            So we simply evaluate for every active basis function v_k the product hess(c).v_k
        */

        hess_expr<gsFeSolution<Scalar>> sHess = hess_expr<gsFeSolution<Scalar>>(_u); // NOTE: This does not parse automatically!
        tmp = sHess.eval(k);  //.transpose();
        vEv = _v.eval(k);
        res = vEv * tmp;
        return res;
    }

    template<class U> inline
    typename util::enable_if<   util::is_same<U,gsGeometryMap<Scalar> >::value ||
                                util::is_same<U,gsFeVariable<Scalar>  >::value, index_t >::type
    cols_impl(const U & u)  const
    {
        return _u.targetDim();
    }

    template<class U> inline
    typename util::enable_if<   util::is_same<U,gsFeSpace<Scalar>    >::value ||
                                util::is_same<U,gsFeSolution<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.dim();
    }
};



/**
 * @brief      Computes the second derivative of an expression
 *
 *  It assumes that the vector of basis functions is of the form v = u*e_i where u
 *  is the scalar basis function u: [0,1]^3 -> R^1 and e_i is the unit vector with a 1 on index i and a 0 elsewhere.
 *  Let us define the following blocks
 *  hess1(u) =              hess2(u) =              hess3(u) =
 *  [d11 u , 0 , 0 ]    |   [0 , d11 u , 0 ]     |  [0 , 0 , d11 u ]
 *  [d22 u , 0 , 0 ]    |   [0 , d22 u , 0 ]     |  [0 , 0 , d22 u ]
 *  [d12 u , 0 , 0 ]    |   [0 , d12 u , 0 ]     |  [0 , 0 , d12 u ]
 *
 *  Then the deriv2(u) is defined as follows (for k number of actives)
 *  [hess1(u)_1]
 *  ...
 *  [hess1(u)_k]
 *  [hess2(u)_1]
 *  ...
 *  [hess2(u)_k]
 *  [hess3(u)_1]
 *  ...
 *  [hess3(u)_k]
 *
 *
 * @tparam     E     Expression type
 */
template<class E>
class deriv2_expr : public _expr<deriv2_expr<E> >
{
    typename E::Nested_t _u;

public:
    // enum {ColBlocks = E::rowSpan }; // ????
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks = (1==E::Space?1:0) };

    typedef typename E::Scalar Scalar;

    deriv2_expr(const E & u) : _u(u) { }

    mutable gsMatrix<Scalar> res, tmp;

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const //(components)
    {
        return rows_impl(_u); // _u.dim() for space or targetDim() for geometry
    }

    index_t cols() const
    {
        return cols_impl(_u);
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); os <<")"; }

    private:
        template<class U> inline
        typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k)  const
        {
            /*
                Here, we compute the hessian of the geometry map.
                The hessian of the geometry map c has the form: hess(c)
                [d11 c1, d11 c2, d11 c3]
                [d22 c1, d22 c2, d22 c3]
                [d12 c1, d12 c2, d12 c3]

                The geometry map has components c=[c1,c2,c3]
            */
            // evaluate the geometry map of U
            res = _u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second );
            return res;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k)  const
        {
            /*
                Here, we compute the hessian of the geometry map.
                The hessian of the geometry map c has the form: hess(c)
                [d11 c1, d11 c2, d11 c3]
                [d22 c1, d22 c2, d22 c3]
                [d12 c1, d12 c2, d12 c3]

                The geometry map has components c=[c1,c2,c3]
            */
            // evaluate the geometry map of U
            tmp =  _u.data().values[2];
            res.resize(rows(),cols());
            for (index_t comp = 0; comp != _u.source().targetDim(); comp++)
                res.col(comp) = tmp.block(comp*rows(),0,rows(),1); //star,length
            return res;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k)  const
        {
            /*
                Here, we compute the hessian of the geometry map.
                The hessian of the geometry map c has the form: hess(c)
                [d11 c1, d11 c2, d11 c3]
                [d22 c1, d22 c2, d22 c3]
                [d12 c1, d12 c2, d12 c3]

                The geometry map has components c=[c1,c2,c3]
            */
            // evaluate the geometry map of U
            hess_expr<gsFeSolution<Scalar>> sHess = hess_expr<gsFeSolution<Scalar>>(_u);
            res = sHess.eval(k).transpose();
            return res;
        }

        /// Spexialization for a space
        template<class U> inline
        typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k) const
        {
            /*
                Here, we compute the hessian of the basis with n actives.
                The hessian of the basis u has the form: hess(u)
                    active 1                active 2                        active n = cardinality
                [d11 u1, d11 u2, d11 u3] [d11 u1, d11 u2, d11 u3] ... [d11 u1, d11 u2, d11 u3]
                [d22 u1, d22 u2, d22 u3] [d22 u1, d22 u2, d22 u3] ... [d22 u1, d22 u2, d22 u3]
                [d12 u1, d12 u2, d12 u3] [d12 u1, d12 u2, d12 u3] ... [d12 u1, d12 u2, d12 u3]

                Here, the basis function has components u = [u1,u2,u3]. Since they are evaluated for scalars
                we use blockDiag to make copies for all components ui

                    active 1     active 2     active k = cardinality/dim   active 1           active 2k       active 1           active 2k
                [d11 u, 0, 0] [d11 u, 0, 0] ... [d11 u, 0, 0]            [0, d11 u, 0]  ... [0, d11 u, 0]  [0, d11 u, 0]  ... [0, d11 u, 0]
                [d22 u, 0, 0] [d22 u, 0, 0] ... [d22 u, 0, 0]            [0, d22 u, 0]  ... [0, d22 u, 0]  [0, d22 u, 0]  ... [0, d22 u, 0]
                [d12 u, 0, 0] [d12 u, 0, 0] ... [d12 u, 0, 0]            [0, d12 u, 0]  ... [0, d12 u, 0]  [0, d12 u, 0]  ... [0, d12 u, 0]

            */
            const index_t numAct = u.data().values[0].rows();   // number of actives of a basis function
            const index_t cardinality = u.cardinality();        // total number of actives (=3*numAct)

            res.resize(rows(), _u.dim() *_u.cardinality()); // (3 x 3*cardinality)
            res.setZero();

            tmp = _u.data().values[2].reshapeCol(k, cols(), numAct );
            for (index_t d = 0; d != cols(); ++d)
            {
                const index_t s = d*(cardinality + 1);
                for (index_t i = 0; i != numAct; ++i)
                    res.col(s+i*_u.cols()) = tmp.col(i);
            }

            return res;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value || util::is_same<U,gsGeometryMap<Scalar> >::value || util::is_same<U,gsFeSpace<Scalar> >::value, index_t >::type
        rows_impl(const U & u)  const
        {
            return _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, index_t >::type
        rows_impl(const U & u) const
        {
            return _u.parDim() * ( _u.parDim() + 1 ) / 2;
        }

        ////////////
        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value || util::is_same<U,gsGeometryMap<Scalar> >::value, index_t >::type
        cols_impl(const U & u)  const
        {
            return u.source().targetDim();
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value || util::is_same<U,gsFeSpace<Scalar> >::value, index_t >::type
        cols_impl(const U & u) const
        {
            return _u.dim();
        }

};

/**
 * @brief      Computes the product of expressions \a E1 and \a E2 and multiplies with a vector \a E3 in voight notation
 *
 * @tparam     E1    Type of expression 1
 * @tparam     E2    Type of expression 2
 * @tparam     E3    Type of the vector expression
 */
template<class E1, class E2, class E3>
class flatdot_expr  : public _expr<flatdot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _A;
    typename E2::Nested_t _B;
    typename E3::Nested_t _C;
    mutable gsMatrix<Scalar> eA, eB, eC, tmp, res;

public:
    enum {Space = 3, ScalarValued = 0, ColBlocks = 0};

public:

    flatdot_expr(_expr<E1> const& A, _expr<E2> const& B, _expr<E3> const& C) : _A(A),_B(B),_C(C)
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t Ac = _A.cols();
        const index_t An = _A.cardinality();
        const index_t Bc = _B.cols();
        const index_t Bn = _B.cardinality();

        eA = _A.eval(k);
        eB = _B.eval(k);
        eC = _C.eval(k);

        GISMO_ASSERT(Bc==_A.rows(), "Dimensions: "<<Bc<<","<< _A.rows()<< "do not match");
        GISMO_STATIC_ASSERT(E1::Space==1, "First entry should be rowSpan"  );
        GISMO_STATIC_ASSERT(E2::Space==2, "Second entry should be colSpan.");

        res.resize(An, Bn);

        // to do: evaluate the jacobians internally for the basis functions and loop over dimensions as well, since everything is same for all dimensions.
        for (index_t i = 0; i!=An; ++i) // for all actives u
            for (index_t j = 0; j!=Bn; ++j) // for all actives v
            {
                tmp.noalias() = eB.middleCols(i*Bc,Bc) * eA.middleCols(j*Ac,Ac);

                tmp(0,0) *= eC.at(0);
                tmp(0,1) *= eC.at(2);
                tmp(1,0) *= eC.at(2);
                tmp(1,1) *= eC.at(1);
                res(i,j) = tmp.sum();
            }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _A.parse(evList);_B.parse(evList);_C.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _A.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _B.colVar(); }
    index_t cardinality_impl() const { return _A.cardinality_impl(); }

    void print(std::ostream &os) const { os << "flatdot("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/**
 * @brief      Computes the product of expressions \a E1 and \a E2 and multiplies with a vector \a E3 in voight notation
 *
 *              NOTE: SPECIALIZED FOR THE SECOND VARIATION OF THE CURVATURE!
 *
 * @tparam     E1    Type of expression 1
 * @tparam     E2    Type of expression 2
 * @tparam     E3    Type of the vector expression
 */
template<class E1, class E2, class E3>
class flatdot2_expr  : public _expr<flatdot2_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _A;
    typename E2::Nested_t _B;
    typename E3::Nested_t _C;
    mutable gsMatrix<Scalar> eA, eB, eC, res, tmp;

public:
    enum {Space = E1::Space, ScalarValued = 0, ColBlocks = 0};

    flatdot2_expr(_expr<E1> const& A, _expr<E2> const& B, _expr<E3> const& C) : _A(A),_B(B),_C(C)
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t Ac = _A.cols();
        const index_t An = _A.cardinality();
        const index_t Bn = _B.cardinality();

        eA = _A.eval(k);
        eB = _B.eval(k);
        eC = _C.eval(k);

        GISMO_ASSERT(_B.rows()==_A.cols(), "Dimensions: "<<_B.rows()<<","<< _A.cols()<< "do not match");
        GISMO_STATIC_ASSERT(E1::Space==1, "First entry should be rowSpan");
        GISMO_STATIC_ASSERT(E2::Space==2, "Second entry should be colSpan.");
        GISMO_ASSERT(_C.cols()==_B.rows(), "Dimensions: "<<_C.rows()<<","<< _B.rows()<< "do not match");

        res.resize(An, Bn);
        for (index_t i = 0; i!=An; ++i)
            for (index_t j = 0; j!=Bn; ++j)
            {
                tmp = eA.middleCols(i*Ac,Ac) * eB.col(j);   // E_f_der2
                tmp.row(2) *= 2.0;                          // multiply the third row of E_f_der2 by 2 for voight notation
                res(i,j) = (eC.row(0) * tmp.col(0)).value();          // E_f^T * mm * E_f_der2
            }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _A.parse(evList);_B.parse(evList);_C.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _A.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _B.colVar(); }
    index_t cardinality_impl() const { return _A.cardinality_impl(); }

    void print(std::ostream &os) const { os << "flatdot2("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/// Expression for the transformation matrix FROM local covariant TO local cartesian bases, based on a geometry map
/// Use of this expression:
/// Let E be a tensor in local covariant coordinates. Then E' = cartcov(G) * E = E.tr() * cartcov(G).tr() is the tensor in local Cartesian basis
///
template<class T> class cartcovinv_expr ;

template<class T>
class cartcov_expr : public _expr<cartcov_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    cartcov_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar> covBasis;
    mutable gsMatrix<Scalar,3,3> result;
    mutable gsVector<Scalar> normal, tmp;
    mutable gsVector<Scalar> e1, e2, a1, a2;

    gsMatrix<Scalar> eval(const index_t k) const
    {
        if (_G.targetDim()==3)
        {
            // Compute covariant bases in deformed and undeformed configuration
            normal = _G.data().normals.col(k);
            normal.normalize();
            covBasis.resize(3,3);
            covBasis.leftCols(2) = _G.data().jacobian(k);
            covBasis.col(2)      = normal;
        }
        else if (_G.targetDim()==2)
        {
            // Compute covariant bases in deformed and undeformed configuration
            covBasis.resize(2,2);
            covBasis = _G.data().jacobian(k);
        }
        else
            GISMO_ERROR("Not implemented");

        a1 = covBasis.col(0);
        a2 = covBasis.col(1);

        e1 = a1.normalized();
        e2 = (a2-(a2.dot(e1)*e1)).normalized();

        result(0,0) = (e1.dot(a1))*(a1.dot(e1));    // 1111
        result(0,1) = (e1.dot(a2))*(a2.dot(e1));    // 1122
        result(0,2) = 2*(e1.dot(a1))*(a2.dot(e1));  // 1112
        // Row 1
        result(1,0) = (e2.dot(a1))*(a1.dot(e2));    // 2211
        result(1,1) = (e2.dot(a2))*(a2.dot(e2));    // 2222
        result(1,2) = 2*(e2.dot(a1))*(a2.dot(e2));  // 2212
        // Row 2
        result(2,0) = (e1.dot(a1))*(a1.dot(e2));    // 1211
        result(2,1) = (e1.dot(a2))*(a2.dot(e2));    // 1222
        result(2,2) = (e1.dot(a1))*(a2.dot(e2)) + (e1.dot(a2))*(a1.dot(e2));    // 1212
        return result;
    }

    cartcovinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return cartcovinv_expr<T>(_G);
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcov("; _G.print(os); os <<")"; }
};

template<class T>
class cartcovinv_expr : public _expr<cartcovinv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    cartcovinv_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<T> temp;

    gsMatrix<T> eval(const index_t k) const
    {
        cartcov_expr<Scalar> cartcov =  cartcov_expr<Scalar>(_G);
        temp = (cartcov.eval(k)).reshape(3,3).inverse();
        return temp;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        cartcov(_G).parse(evList); //

        evList.add(_G);
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcovinv("; _G.print(os); os <<")"; }
};

/// Expression for the transformation matrix FROM local contravariant TO local cartesian bases, based on a geometry map
/// Use of this expression:
/// Let S be a tensor in local covariant coordinates. Then S' = cartcon(G) * S = S.tr() * cartcon(G).tr() is the tensor in local Cartesian basis
template<class T> class cartconinv_expr ;

template<class T>
class cartcon_expr : public _expr<cartcon_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    cartcon_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar> covBasis, conBasis, covMetric, conMetric, cartBasis;
    mutable gsMatrix<Scalar,3,3> result;
    mutable gsVector<Scalar> normal, tmp;
    mutable gsVector<Scalar> e1, e2, ac1, ac2, a1, a2;

    gsMatrix<Scalar> eval(const index_t k) const
    {
        if (_G.targetDim()==3)
        {
            // Compute covariant bases in deformed and undeformed configuration
            normal = _G.data().normals.col(k);
            normal.normalize();

            covBasis.resize(3,3);
            conBasis.resize(3,3);
            covBasis.leftCols(2) = _G.data().jacobian(k);
            covBasis.col(2)      = normal;
            covMetric = covBasis.transpose() * covBasis;

            conMetric = covMetric.inverse();

            conBasis.col(0) = conMetric(0,0)*covBasis.col(0)+conMetric(0,1)*covBasis.col(1)+conMetric(0,2)*covBasis.col(2);
            conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1)+conMetric(1,2)*covBasis.col(2);
        }
        else if (_G.targetDim()==2)
        {
            // Compute covariant bases in deformed and undeformed configuration
            normal.resize(3);
            normal<<0,0,1;
            normal.normalize();
            covBasis.resize(2,2);
            conBasis.resize(2,2);
            // Compute covariant bases in deformed and undeformed configuration
            covBasis = _G.data().jacobian(k);
            covMetric = covBasis.transpose() * covBasis;

            conMetric = covMetric.inverse();

            conBasis.col(0) = conMetric(0,0)*covBasis.col(0)+conMetric(0,1)*covBasis.col(1);
            conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1);
        }
        else
            GISMO_ERROR("Not implemented");

        a1 = covBasis.col(0);
        a2 = covBasis.col(1);

        e1 = a1.normalized();
        e2 = (a2-(a2.dot(e1)*e1)).normalized();

        ac1 = conBasis.col(0);
        ac2 = conBasis.col(1);

        result(0,0) = (e1.dot(ac1))*(ac1.dot(e1));      // 1111
        result(0,1) = (e1.dot(ac2))*(ac2.dot(e1));      // 1122
        result(0,2) = 2*(e1.dot(ac1))*(ac2.dot(e1));    // 1112
        // Row 1
        result(1,0) = (e2.dot(ac1))*(ac1.dot(e2));      // 2211
        result(1,1) = (e2.dot(ac2))*(ac2.dot(e2));      // 2222
        result(1,2) = 2*(e2.dot(ac1))*(ac2.dot(e2));    // 2212
        // Row 2
        result(2,0) = (e1.dot(ac1))*(ac1.dot(e2));      // 1211
        result(2,1) = (e1.dot(ac2))*(ac2.dot(e2));      // 1222
        result(2,2) = (e1.dot(ac1))*(ac2.dot(e2)) + (e1.dot(ac2))*(ac1.dot(e2)); // 1212
        return result;

    }

    cartconinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return cartconinv_expr<T>(_G);
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcon("; _G.print(os); os <<")"; }
};

template<class T>
class cartconinv_expr : public _expr<cartconinv_expr<T> >
{
private:
    typename gsGeometryMap<T>::Nested_t _G;
    mutable gsMatrix<T> temp;

public:
    typedef T Scalar;

    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    cartconinv_expr(const gsGeometryMap<T> & G) : _G(G) { }

    gsMatrix<T> eval(const index_t k) const
    {
        cartcon_expr<Scalar> cartcon =  cartcon_expr<Scalar>(_G);
        temp = (cartcon.eval(k)).reshape(3,3).inverse();
        return temp;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        cartcon(_G).parse(evList); //

        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartconinv("; _G.print(os); os <<")"; }
};

EIGEN_STRONG_INLINE
unitVec_expr uv(const index_t index, const index_t dim) { return unitVec_expr(index,dim); }

template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

template<class E1, class E2> EIGEN_STRONG_INLINE
var2_expr<E1,E2> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G)
{ return var2_expr<E1,E2>(u,v, G); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
var2dot_expr<E1,E2,E3> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & Ef)
{ return var2dot_expr<E1,E2,E3>(u,v, G, Ef); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
var2deriv2dot_expr<E1,E2,E3> var2deriv2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & Ef)
{ return var2deriv2dot_expr<E1,E2,E3>(u,v, G, Ef); }

template<class E> EIGEN_STRONG_INLINE
tvar1_expr<E> tvar1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return tvar1_expr<E>(u, G); }

template<class E> EIGEN_STRONG_INLINE
ovar1_expr<E> ovar1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return ovar1_expr<E>(u, G); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
ovar2dot_expr<E1,E2,E3> ovar2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & C)
{ return ovar2dot_expr<E1,E2,E3>(u,v, G, C); }

// template<class E1, class E2> EIGEN_STRONG_INLINE
// hessdot_expr<E1,E2> hessdot(const E1 & u, const E2 & v) { return hessdot_expr<E1,E2>(u, v); }

template<class E> EIGEN_STRONG_INLINE
deriv2_expr<E> deriv2(const E & u) { return deriv2_expr<E>(u); }

template<class E1, class E2> EIGEN_STRONG_INLINE
deriv2dot_expr<E1, E2> deriv2(const E1 & u, const E2 & v) { return deriv2dot_expr<E1, E2>(u,v); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
flatdot_expr<E1,E2,E3> flatdot(const E1 & u, const E2 & v, const E3 & w)
{ return flatdot_expr<E1,E2,E3>(u, v, w); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
flatdot2_expr<E1,E2,E3> flatdot2(const E1 & u, const E2 & v, const E3 & w)
{ return flatdot2_expr<E1,E2,E3>(u, v, w); }

template<class E> EIGEN_STRONG_INLINE
cartcov_expr<E> cartcov(const gsGeometryMap<E> & G) { return cartcov_expr<E>(G); }

template<class E> EIGEN_STRONG_INLINE
cartcon_expr<E> cartcon(const gsGeometryMap<E> & G) { return cartcon_expr<E>(G); }

}
}
