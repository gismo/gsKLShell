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
#include <gismo.h>

#  define MatExprType  auto

namespace gismo{
namespace expr{

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
    enum{ Space = 0 };

    gsMatrix<Scalar> eval(const index_t) const
    {
        gsMatrix<Scalar> vec = gsMatrix<Scalar>::Zero(_dim,1);
        // vec.setZero();
        vec(_index,0) = 1;
        return vec;
    }

    index_t rows() const { return _dim; }
    index_t cols() const { return  1; }
    void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & ) const {  }

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "uv("<<_dim <<")";}
};

/// Expression for the outer tangent (3D)
template<class T>
class otangent_expr : public _expr<otangent_expr<T> >
{
public:
    typedef T Scalar;

private:

    typename gsGeometryMap<Scalar>::Nested_t _G;

public:

    otangent_expr(const gsGeometryMap<Scalar> & G) : _G(G) { }

    mutable gsVector<Scalar,3> onormal, normal;
    mutable gsMatrix<Scalar> res;

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);

        gsMatrix<Scalar> Jacobian = _G.data().jacobian(k);

        onormal = _G.data().outNormal(k);
        normal =  _G.data().normal(k);
        res = normal.cross(onormal);
        return res;

        // gsMatrix<T> result(_G.data().dim.second,1);
        // result.setZero();
        // return result;
        // if (_G.data().dim.second!=3)
        //     return normal.head(_G.data().dim.second);
        // else
        //     return normal;
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    enum{rowSpan = 0, colSpan = 0};

    void setFlag() const
    {
        _G.data().flags |= NEED_OUTER_NORMAL;
        _G.data().flags |= NEED_NORMAL;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_OUTER_NORMAL;
        _G.data().flags |= NEED_NORMAL;
    }

    // Normalized to unit length
    normalized_expr<otangent_expr<Scalar> > normalized()
    { return normalized_expr<otangent_expr<Scalar> >(*this); }

    void print(std::ostream &os) const { os << "tv("; _G.print(os); os <<")"; }
};

/// Expression for the first variation of the normal
template<class E>
class var1_expr : public _expr<var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:

    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E::Space };

    var1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> bGrads, cJac;
    mutable gsVector<Scalar,3> m_v, normal;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::rowSpan, colSpan = 0};

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

private:
// Specialisation for a space
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        const index_t A = _u.cardinality()/_u.targetDim();
        res.resize(_u.cardinality(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = _u.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // Jac(u) ~ Jac(G) with alternating signs ?..
                m_v.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                              - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() )) / measure;

                // ---------------  First variation of the normal
                // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
                res.row(s+j).noalias() = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
            }
        }
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        res.resize(rows(), cols()); // rows()*
        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        grad_expr<U> vGrad = grad_expr<U>(_u);

        bGrads = vGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }

    // Specialisation for a solution
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        solGrad_expr<Scalar> sGrad =  solGrad_expr<Scalar>(_u);
        res.resize(rows(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = sGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }
};

template<class E1, class E2>
class var1dif_expr : public _expr<var1dif_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E1::Space };

    var1dif_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G) : _u(u), _v(v), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> bGrads, cJac;
    mutable gsVector<Scalar,3> m_v, normal;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,_v,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _v.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _v.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E1::rowSpan, colSpan = 0};

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

private:
    template<class U, class V> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value && util::is_same<V,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const V & v, const index_t k)  const
    {
        res.resize(rows(), cols()); // rows()*
        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        grad_expr<U> uGrad = grad_expr<U>(_u);
        grad_expr<V> vGrad = grad_expr<V>(_v);

        bGrads = uGrad.eval(k) - vGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }

    template<class U, class V> inline
     typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value && util::is_same<V,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const V & v, const index_t k)  const
    {
        GISMO_ASSERT(1==_v.data().actives.cols(), "Single actives expected");
        grad_expr<U> uGrad = grad_expr<U>(_u);
        solGrad_expr<Scalar> vGrad =  solGrad_expr<V>(_v);
        res.resize(rows(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = uGrad.eval(k) - vGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }
};

/// Second variation of the normal times a vector.
template<class E1, class E2, class E3>
class var2_expr : public _expr<var2_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    typename E3::Nested_t _Ef;

public:
    enum{ Space = E1::Space };

    var2_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& Ef) : _u(u),_v(v), _G(G), _Ef(Ef) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2, evEf, result;
    mutable gsVector<Scalar> m_u, m_v, normal, m_uv, m_u_der, n_der, n_der2, tmp; // memomry leaks when gsVector<T,3>, i.e. fixed dimension
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
        const Scalar measure =  _G.data().measures.at(k);

        evEf = _Ef.eval(k);

        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    m_u.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                                     -vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() ))
                                    / measure;

                    const short_t s = d*cardU;

                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;
                        m_v.noalias() = ( vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col(1).template head<3>() )
                                         -vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col(0).template head<3>() ))
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

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_2ND_DER | NEED_MEASURE;
        _Ef.setFlag();
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |=  NEED_NORMAL | NEED_DERIV | NEED_2ND_DER | NEED_MEASURE;
        _Ef.setFlag();
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    enum{rowSpan = 1, colSpan = 1};

    void print(std::ostream &os) const { os << "var2("; _u.print(os); os <<")"; }
};


/// Product of the second derivative of a space or a map and a vector
template<class E1, class E2>
class deriv2dot_expr : public _expr<deriv2dot_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum{ Space = E1::Space };

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

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV2;
        _v.setFlag();
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    enum{rowSpan = E1::rowSpan, colSpan = E2::rowSpan};

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); _v.print(os); os <<")"; }

private:
    /// Specialization for a geometry map
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
        tmp = _u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second);
        vEv = _v.eval(k);
        res = vEv * tmp.transpose();
        return res;
    }

    /// Specialization for a space
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
        const index_t numAct = u.data().values[0].rows(); // number of actives of a basis function
        const index_t cardinality = u.cardinality();      // total number of actives (=3*numAct)
        res.resize(rows() * cardinality, cols());
        tmp.transpose() = _u.data().values[2].reshapeCol(k, cols(), numAct);
        vEv = _v.eval(k);

        for (index_t i = 0; i != _u.dim(); i++)
            res.block(i * numAct, 0, numAct, cols()) = tmp * vEv.at(i);

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

    /// Specialization for a solution
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

        solHess_expr<Scalar> sHess = solHess_expr<Scalar>(_u);
        tmp = sHess.eval(k).transpose();
        vEv = _v.eval(k);
        res = vEv * tmp.transpose();
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
 * @brief   Compute the hessian of a basis

    It assumes that the vector of basis functions is of the form v = u*e_i where u
    is the scalar basis function u: [0,1]^3 -> R^1 and e_i is the unit vector with a 1 on index i and a 0 elsewhere.
    Let us define the following blocks
    hess1(u) =              hess2(u) =              hess3(u) =
    [d11 u , 0 , 0 ]    |   [0 , d11 u , 0 ]     |  [0 , 0 , d11 u ]
    [d22 u , 0 , 0 ]    |   [0 , d22 u , 0 ]     |  [0 , 0 , d22 u ]
    [d12 u , 0 , 0 ]    |   [0 , d12 u , 0 ]     |  [0 , 0 , d12 u ]

    Then the deriv2(u) is defined as follows (for k number of actives)
    [hess1(u)_1]
    ...
    [hess1(u)_k]
    [hess2(u)_1]
    ...
    [hess2(u)_k]
    [hess3(u)_1]
    ...
    [hess3(u)_k]
*/
template<class E>
class deriv2_expr : public _expr<deriv2_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:

    typename E::Nested_t _u;

public:
    enum {ColBlocks = E::rowSpan };
    enum{ Space = E::Space };

    deriv2_expr(const E & u) : _u(u) { }

    mutable gsMatrix<Scalar> res, tmp;

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const //(components)
    {
        // return _u.data().values[2].rows() / _u.data().values[0].rows(); // numHessian dimensions
        // return _u.source().targetDim(); // no. dimensions should be 3
        return rows_impl(_u); // _u.dim() for space or targetDim() for geometry
    }

    index_t cols() const
    {
        return _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV2|NEED_ACTIVE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); os <<")"; }

    private:
        /// Spexialization for a geometry map
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
            solHess_expr<Scalar> sHess = solHess_expr<U>(_u);
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

            // res = _u.data().values[2].col(k).transpose().blockDiag(_u.targetDim()); // analoguous to jacobian..
            // res = res.transpose();
            // gsDebugVar(res);
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
            return _u.dim();
        }

};


// template<class E1, class E2>
// class hessdot_expr : public _expr<hessdot_expr<E1,E2> >
// {
//     typename E1::Nested_t _u;
//     typename E2::Nested_t _v;

// public:
//     enum{ Space = E1::Space };

//     typedef typename E1::Scalar Scalar;

//     hessdot_expr(const E1 & u, const E2 & v) : _u(u), _v(v) {}

//     mutable gsMatrix<Scalar> res, hess, tmp;
//     mutable gsMatrix<Scalar> normalMat;

//     MatExprType eval(const index_t k) const
//     {
//         const gsFuncData<Scalar> & udata = _u.data(); // udata.values[2].col(k)
//         const index_t numAct = udata.values[0].rows();
//         const gsAsConstMatrix<Scalar> ders(udata.values[2].col(k).data(), 3, numAct );

//         tmp = _v.eval(k);

//         res.resize(rows(), cols() );


//             for (index_t i = 0; i!=tmp.rows(); ++i)
//             {
//                 res.block(i*numAct, 0, numAct, 3).noalias() = ders.transpose() * tmp.at(i);
//             }

//         return res;
//     }

//     index_t rows() const
//     {
//         return _u.dim() * _u.data().values[0].rows();
//     }

//     index_t cols() const
//     {
//         return // ( E2::rowSpan() ? rows() : 1 ) *
//             _u.data().values[2].rows() / _u.data().values[0].rows();//=3
//     }

//     void setFlag() const
//     {
//         _u.data().flags |= NEED_2ND_DER;
//         _v.setFlag();
//     }

//     void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
//     {
//         //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
//         evList.push_sorted_unique(&_u.source());
//         _u.data().flags |= NEED_2ND_DER;
//     }

//     const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
//     const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

//     static constexpr bool rowSpan() {return E1::rowSpan(); }
//     static constexpr bool colSpan() {return E2::rowSpan(); }

//     void print(std::ostream &os) const { os << "hessdot("; _u.print(os); os <<")"; }
// };


/// ??
template<class E1, class E2, class E3>
class flatdot_expr  : public _expr<flatdot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, Space = E1::Space};
private:
    typename E1::Nested_t _A;
    typename E2::Nested_t _B;
    typename E3::Nested_t _C;
    mutable gsMatrix<Scalar> eA, eB, eC, tmp, res;

public:

    flatdot_expr(_expr<E1> const& A, _expr<E2> const& B, _expr<E3> const& C) : _A(A),_B(B),_C(C)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
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
        GISMO_STATIC_ASSERT(E1::rowSpan, "First entry should be rowSpan"  );
        GISMO_STATIC_ASSERT(E2::colSpan, "Second entry should be colSpan.");

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
    void setFlag() const { _A.setFlag();_B.setFlag();_C.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _A.parse(evList);_B.parse(evList);_C.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _A.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _B.colVar(); }
    index_t cardinality_impl() const { return _A.cardinality_impl(); }

    enum{rowSpan = E1::rowSpan, colSpan = E2::colSpan};

    void print(std::ostream &os) const { os << "flatdot("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/// ???
template<class E1, class E2, class E3>
class flatdot2_expr  : public _expr<flatdot2_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, Space = E1::Space};
    enum {ColBlocks = 0 };

private:
    typename E1::Nested_t _A;
    typename E2::Nested_t _B;
    typename E3::Nested_t _C;
    mutable gsMatrix<Scalar> eA, eB, eC, res, tmp;

public:

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
        GISMO_STATIC_ASSERT(E1::rowSpan, "First entry should be rowSpan");
        GISMO_STATIC_ASSERT(E2::colSpan, "Second entry should be colSpan.");
        GISMO_ASSERT(_C.cols()==_B.rows(), "Dimensions: "<<_C.rows()<<","<< _B.rows()<< "do not match");

        res.resize(An, Bn);
        for (index_t i = 0; i!=An; ++i)
            for (index_t j = 0; j!=Bn; ++j)
            {
                tmp = eA.middleCols(i*Ac,Ac) * eB.col(j);   // E_f_der2
                tmp.row(2) *= 2.0;                          // multiply the third row of E_f_der2 by 2 for voight notation
                res(i,j) = eC.row(0) * tmp.col(0);          // E_f^T * mm * E_f_der2
            }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }
    void setFlag() const { _A.setFlag();_B.setFlag();_C.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _A.parse(evList);_B.parse(evList);_C.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _A.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _B.colVar(); }
    index_t cardinality_impl() const { return _A.cardinality_impl(); }

    enum{rowSpan = E1::rowSpan, colSpan = E2::colSpan};

    void print(std::ostream &os) const { os << "flatdot2("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/// Expression for the transformation matrix FROM local covariant TO local cartesian bases, based on a geometry map
template<class T> class cartcovinv_expr ;

template<class T>
class cartcov_expr : public _expr<cartcov_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartcov_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar> covBasis, conBasis, covMetric, conMetric, cartBasis;
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
            conBasis.resize(3,3);
            covBasis.leftCols(2) = _G.data().jacobian(k);
            covBasis.col(2)      = normal;
            covMetric = covBasis.transpose() * covBasis;

            conMetric = covMetric.inverse();

            conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1)+conMetric(1,2)*covBasis.col(2);

            e1 = covBasis.col(0); e1.normalize();
            e2 = conBasis.col(1); e2.normalize();

            a1 = covBasis.col(0);
            a2 = covBasis.col(1);

            result(0,0) = (e1.dot(a1))*(a1.dot(e1));
            result(0,1) = (e1.dot(a2))*(a2.dot(e2));
            result(0,2) = 2*(e1.dot(a1))*(a2.dot(e1));
            // Row 1
            result(1,0) = (e2.dot(a1))*(a1.dot(e2));
            result(1,1) = (e2.dot(a2))*(a2.dot(e2));
            result(1,2) = 2*(e2.dot(a1))*(a2.dot(e2));
            // Row 2
            result(2,0) = (e1.dot(a1))*(a1.dot(e2));
            result(2,1) = (e1.dot(a2))*(a2.dot(e2));
            result(2,2) = (e1.dot(a1))*(a2.dot(e2)) + (e1.dot(a2))*(a1.dot(e2));

            // return result.inverse(); // !!!!
            return result;
        }
        else if (_G.targetDim()==2)
        {
            // Compute covariant bases in deformed and undeformed configuration
            covBasis.resize(2,2);
            conBasis.resize(2,2);
            covBasis = _G.data().jacobian(k);
            covMetric = covBasis.transpose() * covBasis;
            conMetric = covMetric.inverse();
            conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1);

            e1 = covBasis.col(0); e1.normalize();
            e2 = conBasis.col(1); e2.normalize();
            // e3 = normal;

            a1 = covBasis.col(0);
            a2 = covBasis.col(1);

            result(0,0) = (e1.dot(a1))*(a1.dot(e1));
            result(0,1) = (e1.dot(a2))*(a2.dot(e2));
            result(0,2) = 2*(e1.dot(a1))*(a2.dot(e1));
            // Row 1
            result(1,0) = (e2.dot(a1))*(a1.dot(e2));
            result(1,1) = (e2.dot(a2))*(a2.dot(e2));
            result(1,2) = 2*(e2.dot(a1))*(a2.dot(e2));
            // Row 2
            result(2,0) = (e1.dot(a1))*(a1.dot(e2));
            result(2,1) = (e1.dot(a2))*(a2.dot(e2));
            result(2,2) = (e1.dot(a1))*(a2.dot(e2)) + (e1.dot(a2))*(a1.dot(e2));

            return result;
        }
        else
            GISMO_ERROR("Not implemented");
    }

    cartcovinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return cartcovinv_expr<T>(_G);
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    enum{rowSpan = 0, colSpan = 0};

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
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

    cartcovinv_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<T> temp;

    gsMatrix<T> eval(const index_t k) const
    {
        temp = (cartcov_expr<gsGeometryMap<T> >(_G).eval(k)).reshape(3,3).inverse();
        return temp;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    enum{rowSpan = 0, colSpan = 0};

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcovinv("; _G.print(os); os <<")"; }
};


/// Expression for the transformation matrix FROM local contravariant TO local cartesian bases, based on a geometry map
template<class T> class cartconinv_expr ;

template<class T>
class cartcon_expr : public _expr<cartcon_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartcon_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar> covBasis, conBasis, covMetric, conMetric, cartBasis;
    mutable gsMatrix<Scalar,3,3> result;
    mutable gsVector<Scalar> normal, tmp;
    mutable gsVector<Scalar> e1, e2, ac1, ac2;

    gsMatrix<Scalar> eval(const index_t k) const
    {
        if (_G.targetDim()==3)
        {
            // Compute covariant bases in deformed and undeformed configuration
            normal = _G.data().normals.col(k);
            normal.normalize();
            covBasis.resize(3,3);
            conBasis.resize(3,3);
            // Compute covariant bases in deformed and undeformed configuration
            normal = _G.data().normals.col(k);
            normal.normalize();
            covBasis.leftCols(2) = _G.data().jacobian(k);
            covBasis.col(2)      = normal;
            covMetric = covBasis.transpose() * covBasis;

            conMetric = covMetric.inverse();

            conBasis.col(0) = conMetric(0,0)*covBasis.col(0)+conMetric(0,1)*covBasis.col(1)+conMetric(0,2)*covBasis.col(2);
            conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1)+conMetric(1,2)*covBasis.col(2);

            e1 = covBasis.col(0); e1.normalize();
            e2 = conBasis.col(1); e2.normalize();
            // e3 = normal;

            ac1 = conBasis.col(0);
            ac2 = conBasis.col(1);

            result(0,0) = (e1.dot(ac1))*(ac1.dot(e1));
            result(0,1) = (e1.dot(ac2))*(ac2.dot(e2));
            result(0,2) = 2*(e1.dot(ac1))*(ac2.dot(e1));
            // Row 1
            result(1,0) = (e2.dot(ac1))*(ac1.dot(e2));
            result(1,1) = (e2.dot(ac2))*(ac2.dot(e2));
            result(1,2) = 2*(e2.dot(ac1))*(ac2.dot(e2));
            // Row 2
            result(2,0) = (e1.dot(ac1))*(ac1.dot(e2));
            result(2,1) = (e1.dot(ac2))*(ac2.dot(e2));
            result(2,2) = (e1.dot(ac1))*(ac2.dot(e2)) + (e1.dot(ac2))*(ac1.dot(e2));

            return result;
        }
        else if (_G.targetDim()==2)
        {
            // Compute covariant bases in deformed and undeformed configuration
            normal = _G.data().normals.col(k);
            normal.normalize();
            covBasis.resize(2,2);
            conBasis.resize(2,2);
            // Compute covariant bases in deformed and undeformed configuration
            covBasis = _G.data().jacobian(k);
            covMetric = covBasis.transpose() * covBasis;

            conMetric = covMetric.inverse();

            conBasis.col(0) = conMetric(0,0)*covBasis.col(0)+conMetric(0,1)*covBasis.col(1);
            conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1);

            e1 = covBasis.col(0); e1.normalize();
            e2 = conBasis.col(1); e2.normalize();

            ac1 = conBasis.col(0);
            ac2 = conBasis.col(1);

            result(0,0) = (e1.dot(ac1))*(ac1.dot(e1));
            result(0,1) = (e1.dot(ac2))*(ac2.dot(e2));
            result(0,2) = 2*(e1.dot(ac1))*(ac2.dot(e1));
            // Row 1
            result(1,0) = (e2.dot(ac1))*(ac1.dot(e2));
            result(1,1) = (e2.dot(ac2))*(ac2.dot(e2));
            result(1,2) = 2*(e2.dot(ac1))*(ac2.dot(e2));
            // Row 2
            result(2,0) = (e1.dot(ac1))*(ac1.dot(e2));
            result(2,1) = (e1.dot(ac2))*(ac2.dot(e2));
            result(2,2) = (e1.dot(ac1))*(ac2.dot(e2)) + (e1.dot(ac2))*(ac1.dot(e2));

            return result;
        }
        else
            GISMO_ERROR("Not implemented");

    }

    cartconinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return cartconinv_expr<T>(_G);
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    enum{rowSpan = 0, colSpan = 0};

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcon("; _G.print(os); os <<")"; }
};

template<class T>
class cartconinv_expr : public _expr<cartconinv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartconinv_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<T> temp;

    gsMatrix<T> eval(const index_t k) const
    {
        temp = (cartcon_expr<gsGeometryMap<T> >(_G).eval(k)).reshape(3,3).inverse();
        return temp;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    enum{rowSpan = 0, colSpan = 0};

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartconinv("; _G.print(os); os <<")"; }
};

EIGEN_STRONG_INLINE
unitVec_expr uv(const index_t index, const index_t dim) { return unitVec_expr(index,dim); }

template<class T> EIGEN_STRONG_INLINE
otangent_expr<T> otangent(const gsGeometryMap<T> & u) { return otangent_expr<T>(u); }

template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

template<class E1, class E2> EIGEN_STRONG_INLINE
var1dif_expr<E1,E2> var1dif(const E1 & u,const E2 & v, const gsGeometryMap<typename E1::Scalar> & G) { return var1dif_expr<E1,E2>(u,v, G); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
var2_expr<E1,E2,E3> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & Ef)
{ return var2_expr<E1,E2,E3>(u,v, G, Ef); }

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
