/** @file gsEmbeddingUtils.h

    @brief Utilities for gsEmbedding. Mainly expressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        C. Chianese     (2025 UniNa)
        H.M. Verhelst   (2025 TU Delft)
*/

#pragma once

//! [Include namespace] 
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsDofMapper.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsUtils/gsPointGrid.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsAssembler/gsExpressions.h>

namespace gismo{
namespace expr{

template<class T>
class curve_tangent_expr : public _expr<curve_tangent_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S),_C(C),
                        Spatch(_S.source().piece(0)) //<<--- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, dtheta, sJac, res;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].reshapeCol(k,1,_C.data().dim.second).transpose();
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();
        res.resize(3,1);

        if (_S.source().targetDim() == 2)
        {
            res.col(0).head(2) = sJac*dtheta;
            res(2,0)           = 0;
        }
        else if (_S.source().targetDim()==3)
        {
            res = sJac * dtheta;
        }
        else
        {
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
        }

        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_normal_expr : public _expr<curve_normal_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_normal_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C) : _S(S), _C(C),
                        Spatch(S.source().piece(0)) //<<--- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2, but is "<<_S.source().domainDim());
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1, but is "<<_C.source().domainDim());
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, sJac, res;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // Make this work for a 2D planar surface AND for a 3D surface:
    // - 2D: normal vector is always (0,0,1)
    // - 3D: normal vector is the cross product of the tangent vectors
    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d) to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        if (_S.source().targetDim() == 2)
        {
            res.resize(3,1);
            res << 0, 0, 1;
            return res;
        }
        else if (_S.source().targetDim()==3)
        {
            theta = _C.data().values[0].col(k);
            sJac  = Spatch.deriv(theta);
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            sJac.transposeInPlace();
            res = sJac.col3d(0).cross(sJac.col3d(1));
            // res.normalize();
            return res;
        }
        else
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_binormal_expr : public _expr<curve_binormal_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S),_C(C) 
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2, but is "<<_S.source().domainDim());
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1, but is "<<_C.source().domainDim());
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> cnvec, ctvec, res;
    //const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d) to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        curve_normal_expr<Scalar>  cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        res = ctvec.col3d(0).cross(cnvec.col3d(0));
        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C); // WHY IS THIS NEEDED?
        _C.data().flags |= NEED_VALUE;

        cnv(_S,_C).parse(evList);
        ctv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

/**
 * @brief First variation of basis vectors wrt DOFs.
 */

template<class E>
class curve_tangent_var1_expr : public _expr<curve_tangent_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_var1_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _u(u), _S(S), _C(C)
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, res;
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
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();
        bGrads = mb.basis(0).deriv(theta);

        res.resize(A*_u.dim(),cols());

        for (index_t d = 0; d!= _u.dim(); ++d) // for all DOF components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                res.row(s+j) = vecFun(d, bGrads.at(2*j  ))* dtheta(0,0)
                                + vecFun(d, bGrads.at(2*j+1))* dtheta(1,0);
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return _u.rowVar();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv_var1("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class E>
class curve_normal_var1_expr : public _expr<curve_normal_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_normal_var1_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _u(u), _S(S), _C(C),
                            Spatch(_S.source().piece(0)) //<<--- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, bGrads, sJac, cnvec, res, var;
    const gsFunctionSet<Scalar> & Spatch;
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
        theta  = _C.data().values[0].col(k);
        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();

        res.resize(A*_u.dim(), cols());

        if (_S.source().targetDim() == 2)
        {
            res.setZero();
            return res;
        }
        else if (_S.source().targetDim()==3)
        {
            bGrads = mb.basis(0).deriv(theta);
            sJac  = Spatch.deriv(theta);
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            sJac.transposeInPlace();
            curve_normal_expr<Scalar> cnv(_S,_C);
            cnvec = cnv.eval(k);
            //const Scalar measure = cnvec.norm();
            cnvec.normalize();

            for (index_t d = 0; d!= cols(); ++d) // for all DOF components
            {
                const short_t s = d*A;
                for (index_t j = 0; j!= A; ++j) // for all actives
                {
                    //first variation of non-unit normal vector
                    res.row(s+j) =  (vecFun(d, bGrads.at(2*j  )).cross(sJac.col3d(1))
                                    - vecFun(d, bGrads.at(2*j+1)).cross(sJac.col3d(0))).transpose();

                    // // first variation of non-unit normal vector (divided by measure)
                    // auto var =  (vecFun(d, bGrads.at(2*j  )).cross(sJac.col3d(1))
                    //       - vecFun(d, bGrads.at(2*j+1)).cross(sJac.col3d(0)))/measure;

                    // // first variation of unit normal vector
                    // res.row(s+j) = (var - ((cnvec.col3d(0)*var.transpose()) * cnvec.col3d(0))).transpose();
                }
            }
            return res;
        }
        else
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        cnv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "cnv_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_var1_expr : public _expr<curve_binormal_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_var1_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _u(u), _S(S), _C(C)
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> ctvec, cnvec, cbvec, ctvec_var1, cnvec_var1, theta, var, res;
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
        theta = _C.data().values[0].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();

        var.resize(A*_u.dim(), cols());
        res.resize(A*_u.dim(), cols());

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctvec_var1 = ctv_var1.eval(k).transpose();
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_binormal_expr<Scalar> cbv(_S,_C);
        cbvec = cbv.eval(k);
        const Scalar measure = cbvec.norm();
        cbvec.normalize();

        for (index_t i = 0; i!= A*_u.dim(); ++i)
        {
            //first variation of non-unit binormal vector
            res.row(i) = ctvec_var1.col3d(i).cross(cnvec.col3d(0)) +
                            ctvec.col3d(0).cross(cnvec_var1.col3d(i));

            //first variation of non-unit binormal vector (divided by measure)
            // var.row(i) = (ctvec_var1.col3d(i).cross(cnvec.col3d(0)) +
                            // ctvec.col3d(0).cross(cnvec_var1.col3d(i)))/measure;

            //first variation of unit binormal vector
            // res.row(i) = (var.row(i).transpose() - (cbvec.col3d(0)*var.row(i)) *cbvec.col3d(0)).transpose();
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cbv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "cbv_var1("; _u.print(os); os <<")"; }

};

/**
 * @brief Derivative of basis vectors along the curve.
 */

template<class T>
class curve_tangent_varalong_expr : public _expr<curve_tangent_varalong_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_varalong_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S), _C(C),
                                Spatch(S.source().piece(0)) //<<--- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {}

    mutable gsMatrix<Scalar> theta, ddtheta, sJac, sHess, res;
    mutable gsVector<Scalar,3> vec;
    mutable gsEigen::Array<Scalar,2,1> dtheta;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        ddtheta = _C.data().values[2].col(k);
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();
        sHess = Spatch.deriv2(theta);
        sHess.resize(numRows,_S.source().targetDim());
        sHess.transposeInPlace();
        vec.head(2) = dtheta.pow(2);
        vec(2) = 2*dtheta(0,0)*dtheta(1,0);
        
        res = sJac*ddtheta + sHess*vec;
        //res.normalize();
        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV | NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv_vara("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_normal_varalong_expr : public _expr<curve_normal_varalong_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S), _C(C),
                                Spatch(S.source().piece(0)) //<<--- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {}

    mutable gsMatrix<Scalar> theta, dtheta, sJac, sHess, res, vecA, vecB;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();
        sHess = Spatch.deriv2(theta);
        sHess.resize(numRows,_S.source().targetDim());
        sHess.transposeInPlace();
        sHess.col(2).swap(sHess.col(1));
        vecA = sHess.leftCols(2)*dtheta;
        vecB = sHess.rightCols(2)*dtheta;

        res  = vecA.col3d(0).cross(sJac.col3d(1)) + sJac.col3d(0).cross(vecB.col3d(0));
        //res.normalize();
        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv_varalong("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_normal_varalong_normalized_expr : public _expr<curve_normal_varalong_normalized_expr<T>>
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_normalized_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S), _C(C) {}

    mutable gsMatrix<Scalar> theta, cnvec, cnvec_vara, res;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta = _C.data().values[0].col(k);

        curve_normal_expr<T> cnv(_S,_C);
        cnvec = cnv.eval(k);
        const Scalar measure = cnvec.norm();
        curve_normal_varalong_expr<T> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        const Scalar measure1 = (cnvec.col3d(0).dot(cnvec_vara.col3d(0)))/measure;
        
        res  = cnvec_vara.col3d(0)/measure - (cnvec.col3d(0)*measure1)/std::pow(measure,2);
        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv_varalong_normalized("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_binormal_varalong_expr : public _expr<curve_binormal_varalong_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S), _C(C) {}

    mutable gsMatrix<Scalar> ctvec, ctvec_vara, cnvec, cnvec_vara, res;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        curve_tangent_varalong_expr<T> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_normal_varalong_expr<T> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_tangent_expr<T> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_normal_expr<T> cnv(_S,_C);
        cnvec = cnv.eval(k);

        res = ctvec_vara.col3d(0).cross(cnvec.col3d(0)) +
              ctvec.col3d(0).cross(cnvec_vara.col3d(0));

        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // WHY IS THIS NEEDED?
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;

        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cbv_varalong("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_binormal_varalong_normalized_expr : public _expr<curve_binormal_varalong_normalized_expr<T>>
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_normalized_expr(const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _S(S), _C(C) {}

    mutable gsMatrix<Scalar> theta, cbvec, cbvec_vara, res;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta = _C.data().values[0].col(k);

        curve_binormal_expr<T> cbv(_S,_C);
        cbvec = cbv.eval(k);
        const Scalar measure = cbvec.norm();
        curve_binormal_varalong_expr<T> cbv_vara(_S,_C);
        cbvec_vara = cbv_vara.eval(k);
        const Scalar measure1 = (cbvec.col3d(0).dot(cbvec_vara.col3d(0)))/measure;
        
        res  = cbvec_vara.col3d(0)/measure - (cbvec.col3d(0)*measure1)/std::pow(measure,2);
        return res;
    }

    index_t rows() const { return 3; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        cbv(_S,_C).parse(evList);
        cbv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cbv_varalong_normalized("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

/**
 * @brief First variation of the derivative of basis vectors along the curve wrt DOFs.
 */

template<class E>
class curve_tangent_varalong_var1_expr : public _expr<curve_tangent_varalong_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_varalong_var1_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _u(u), _S(S), _C(C) {}

    mutable gsMatrix<Scalar> theta, res;
    mutable gsEigen::Array<Scalar,2,1> dtheta;
    mutable gsVector<Scalar> bGrads, bHess, ddtheta;
    mutable gsVector<Scalar,3> vec;
    mutable Scalar dotA, dotB;
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
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        ddtheta = _C.data().values[2].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); 
        bGrads = mb.basis(0).deriv(theta);
        bHess =  mb.basis(0).deriv2(theta);
        vec.head(2) = dtheta.pow(2);
        vec(2) = 2*dtheta(0,0)*dtheta(1,0);

        res.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d) // for all DOF components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j)  // for all actives
            {
                dotA = bHess.middleRows(3*j,3).dot(vec);
                dotB = bGrads.middleRows(2*j,2).dot(ddtheta);
                res.row(s+j) = vecFun(d, dotA) + vecFun(d, dotB);
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV | NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "ctv_vara_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_normal_varalong_var1_expr : public _expr<curve_normal_varalong_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_var1_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _u(u), _S(S), _C(C),
                                    Spatch(S.source().piece(0)) //<<--- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
                                    {}

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, bHess, localbHess, sJac, sHess, res, vecA, vecB;
    mutable gsVector<Scalar,2> dot;
    const gsFunctionSet<Scalar> & Spatch;
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
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();
        bGrads = mb.basis(0).deriv(theta);
        bHess =  mb.basis(0).deriv2(theta);

        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();
        sHess = Spatch.deriv2(theta);
        sHess.resize(numRows,_S.source().targetDim());
        sHess.transposeInPlace();
        sHess.col(2).swap(sHess.col(1));
        vecA = sHess.leftCols(2) * dtheta;
        vecB = sHess.rightCols(2) * dtheta;

        res.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d) // for all DOF components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                secDerToHessian(bHess.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dot = localbHess * dtheta;

                //not-normalized result
                res.row(s+j) =  vecFun(d,dot(0)).cross(sJac.col3d(1)) +
                                vecA.col3d(0).cross(vecFun(d,bGrads.at(2*j+1))) +
                                vecFun(d,bGrads.at(2*j  )).cross(vecB.col3d(0)) +
                                sJac.col3d(0).cross(vecFun(d,dot(1)));
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "cnv_vara_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_normal_varalong_var1_normalized_expr : public _expr<curve_normal_varalong_var1_normalized_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_var1_normalized_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C):
                                                _u(u), _S(S), _C(C) {}

    mutable gsMatrix<Scalar> theta, cnvec, cnvec_vara, cnvec_var1, cnvec_vara_var1, res;
    mutable Scalar measure, measure1, measurer, measure1r;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();

        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_normal_varalong_var1_expr<E> cnv_vara_var1(_u,_S,_C);
        cnvec_vara_var1 = cnv_vara_var1.eval(k).transpose();
        measure = cnvec.norm();
        measure1 = (cnvec.col3d(0).dot(cnvec_vara.col3d(0)))/measure;

        res.resize(A*_u.dim(), cols());

        for (index_t r = 0; r!= A*_u.dim(); ++r) // loop over cardinality
        {
            measurer =  (cnvec.col3d(0).dot(cnvec_var1.col3d(r)))/measure;
            measure1r = cnvec_var1.col3d(r).dot(cnvec_vara.col3d(0))/measure +
                        cnvec.col3d(0).dot(cnvec_vara_var1.col3d(r))/measure -
                        (cnvec.col3d(0).dot(cnvec_vara.col3d(0))*measurer)/std::pow(measure,2);

            res.row(r) = cnvec_vara_var1.col3d(r)/measure - cnvec_vara.col3d(0)*measurer/std::pow(measure,2) -
                         cnvec_var1.col3d(r)*measure1/std::pow(measure,2) - cnvec.col3d(0)*measure1r/std::pow(measure,2) +
                         2*measurer*measure1*cnvec.col3d(0)/std::pow(measure,3);
        }

        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cnv_vara_var1(_u,_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "cnv_vara_var1_normalized("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_varalong_var1_expr : public _expr<curve_binormal_varalong_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_var1_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C): _u(u), _S(S), _C(C) {}

    mutable gsMatrix<Scalar> ctvec, ctvec_vara, ctvec_var1, ctvec_vara_var1, theta;
    mutable gsMatrix<Scalar> cnvec, cnvec_vara, cnvec_var1, cnvec_vara_var1, res; 
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {   
        theta  = _C.data().values[0].col(k);
        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();
        
        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_varalong_expr<Scalar> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctvec_var1 = ctv_var1.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E> ctv_vara_var1(_u,_S,_C);
        ctvec_vara_var1 = ctv_vara_var1.eval(k).transpose();

        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_normal_varalong_var1_expr<E> cnv_vara_var1(_u,_S,_C);
        cnvec_vara_var1 = cnv_vara_var1.eval(k).transpose();

        res.resize(A*_u.dim(),cols());

        for (index_t i = 0; i!= A*_u.dim(); ++i)
        {
            //not-normalized result
            res.row(i) =  ctvec_vara_var1.col3d(i).cross(cnvec.col3d(0)) +
                            ctvec_vara.col3d(0).cross(cnvec_var1.col3d(i)) +
                            ctvec_var1.col3d(i).cross(cnvec_vara.col3d(0)) +
                            ctvec.col3d(0).cross(cnvec_vara_var1.col3d(i));
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        ctv_vara_var1(_u,_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cnv_vara_var1(_u,_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "cbv_vara_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_varalong_var1_normalized_expr : public _expr<curve_binormal_varalong_var1_normalized_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_var1_normalized_expr(const E &u, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C):
                                                _u(u), _S(S), _C(C) {}

    mutable gsMatrix<Scalar> theta, cbvec, cbvec_vara, cbvec_var1, cbvec_vara_var1, res;
    mutable Scalar measure, measure1, measurer, measure1r;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();

        curve_binormal_expr<Scalar> cbv(_S,_C);
        cbvec = cbv.eval(k);
        curve_binormal_varalong_expr<Scalar> cbv_vara(_S,_C);
        cbvec_vara = cbv_vara.eval(k);
        curve_binormal_var1_expr<E> cbv_var1(_u,_S,_C);
        cbvec_var1 = cbv_var1.eval(k).transpose();
        curve_binormal_varalong_var1_expr<E> cbv_vara_var1(_u,_S,_C);
        cbvec_vara_var1 = cbv_vara_var1.eval(k).transpose();
        measure = cbvec.norm();
        measure1 = (cbvec.col3d(0).dot(cbvec_vara.col3d(0)))/measure;

        res.resize(A*_u.dim(), cols());

        for (index_t r = 0; r!= A*_u.dim(); ++r) // loop over cardinality
        {
            measurer =  (cbvec.col3d(0).dot(cbvec_var1.col3d(r)))/measure;
            measure1r =  cbvec_var1.col3d(r).dot(cbvec_vara.col3d(0))/measure +
                         cbvec.col3d(0).dot(cbvec_vara_var1.col3d(r))/measure -
                        (cbvec.col3d(0).dot(cbvec_vara.col3d(0))*measurer)/std::pow(measure,2);

            res.row(r) = cbvec_vara_var1.col3d(r)/measure - cbvec_vara.col3d(0)*measurer/std::pow(measure,2) -
                         cbvec_var1.col3d(r)*measure1/std::pow(measure,2) - cbvec.col3d(0)*measure1r/std::pow(measure,2) +
                         2*measurer*measure1*cbvec.col3d(0)/std::pow(measure,3);
        }

        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        cbv(_S,_C).parse(evList);
        cbv_vara(_S,_C).parse(evList);
        cbv_var1(_u,_S,_C).parse(evList);
        cbv_vara_var1(_u,_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "cbv_vara_var1_normalized("; _u.print(os); os <<")"; }

};

/**
 * @brief Second variation of derivative of basis vector wrt DOFs times a vector.
 */

template<class E1, class E2, class E3>
class curve_normal_var2dot_expr : public _expr<curve_normal_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_normal_var2dot_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &C, _expr<E3> const& vec):
                              _u(u), _v(v), _C(C), _vec(vec) {}

    mutable gsMatrix<Scalar> theta, bGradsu, bGradsv, evec, res;

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
        theta  = _C.data().values[0].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows();
        const index_t Av = mbv.basis(0).active(theta).rows();

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);

        res.resize(Au*_u.dim(), Av*_v.dim());
        
        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {   
                        res(s+j,r+k) = (vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                        vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1)))).dot(evec.col(0));
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
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_normal_var2dot("; _u.print(os); os <<")"; }

};

template<class E1, class E2, class E3>
class curve_binormal_var2dot_expr : public _expr<curve_binormal_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_var2dot_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &C,
                                _expr<E3> const& vec): _u(u), _v(v), _S(S), _C(C), _vec(vec) {}

    mutable gsMatrix<Scalar> theta, bGradsu, bGradsv, evec, res;
    mutable gsMatrix<Scalar> ctvec, ctvec_var1u, ctvec_var1v, cnvec_var1u, cnvec_var1v;
    mutable gsVector<Scalar,3> cnvec_var2;

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
        theta  = _C.data().values[0].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows();
        const index_t Av = mbv.basis(0).active(theta).rows();

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_var1_expr<E1> ctv_var1u(_u,_S,_C);
        ctvec_var1u = ctv_var1u.eval(k).transpose();
        curve_tangent_var1_expr<E2> ctv_var1v(_v,_S,_C);
        ctvec_var1v = ctv_var1v.eval(k).transpose();
        curve_normal_var1_expr<E1> cnv_var1u(_u,_S,_C);
        cnvec_var1u = cnv_var1u.eval(k).transpose();
        curve_normal_var1_expr<E2> cnv_var1v(_v,_S,_C);
        cnvec_var1v = cnv_var1v.eval(k).transpose();
                
        res.resize(Au*_u.dim(), Av*_v.dim());
        
        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {   
                        cnvec_var2 = vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                     vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1)));

                        res(s+j,r+k) = (ctvec_var1u.col3d(s+j).cross(cnvec_var1v.col3d(r+k)) +
                                        ctvec_var1v.col3d(r+k).cross(cnvec_var1u.col3d(s+j)) + 
                                        ctvec.col3d(0).cross(cnvec_var2)).dot(evec.col3d(0));
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
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        ctv_var1(_v,_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cnv_var1(_v,_S,_C).parse(evList);
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_binormal_var2dot("; _u.print(os); os <<")"; }
};

/**
 * @brief Second variation of derivative of basis vector along the curve wrt DOFs times a vector.
 */

 template<class E1, class E2, class E3>
 class curve_normal_varalong_var2dot_expr : public _expr<curve_normal_varalong_var2dot_expr<E1,E2,E3> >
 {
 public:
     typedef typename E1::Scalar Scalar;
 
 private:
     typename E1::Nested_t _u;
     typename E2::Nested_t _v;
     typename gsGeometryMap<Scalar>::Nested_t _S;
     typename gsGeometryMap<Scalar>::Nested_t _C;
     typename E3::Nested_t _vec;
 
 public:
     enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};
 
     curve_normal_varalong_var2dot_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &S,
                                        const gsGeometryMap<Scalar> &C, _expr<E3> const &vec): _u(u), _v(v), _S(S), _C(C), _vec(vec)
                                        {}
 
     mutable gsMatrix<Scalar> theta, dtheta, bGradsu, bGradsv, bHessu, bHessv, localbHess, dotu, dotv, evec, res;
     mutable gsMatrix<Scalar> cnvec, cnvec_vara, cnvec_var1u, cnvec_var1v, cnvec_vara_var1u, cnvec_vara_var1v;
     mutable gsVector<Scalar,3> cnvec_var2, cnvec_vara_var2;
     mutable gsMatrix<gsVector<Scalar,3>> matrix;
     mutable Scalar measurer, measures, measure1, measure1r, measure1s, measurers, measure1rs;
 
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
         theta  = _C.data().values[0].col(k);
         dtheta = _C.data().values[1].col(k);
         evec = _vec.eval(k);
 
         const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
         const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
         const index_t Au = mbu.basis(0).active(theta).rows(); 
         const index_t Av = mbv.basis(0).active(theta).rows();
 
         bGradsu = mbu.basis(0).deriv(theta);
         bGradsv = mbv.basis(0).deriv(theta);
         bHessu =  mbu.basis(0).deriv2(theta);
         bHessv =  mbv.basis(0).deriv2(theta);
 
         curve_normal_expr<Scalar> cnv(_S,_C);
         cnvec = cnv.eval(k);
         const Scalar measure = cnvec.norm();
         curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
         cnvec_vara = cnv_vara.eval(k);
         curve_normal_var1_expr<E1> cnv_var1u(_u,_S,_C);
         cnvec_var1u = cnv_var1u.eval(k).transpose();
         curve_normal_var1_expr<E2> cnv_var1v(_v,_S,_C);
         cnvec_var1v = cnv_var1v.eval(k).transpose();
         curve_normal_varalong_var1_expr<E1> cnv_vara_var1u(_u,_S,_C);
         cnvec_vara_var1u = cnv_vara_var1u.eval(k).transpose();
         curve_normal_varalong_var1_expr<E2> cnv_vara_var1v(_v,_S,_C);    
         cnvec_vara_var1v = cnv_vara_var1v.eval(k).transpose();
 
         matrix.resize(Au*_u.dim(), Av*_v.dim());
         res.resize(Au*_u.dim(), Av*_v.dim());
 
         for (index_t d = 0; d!= _u.dim(); ++d) // for all components
         {
             const short_t s = d*Au;
             for (index_t j = 0; j!= Au; ++j) // for all actives
             {
                 secDerToHessian(bHessu.col(0).segment(3*j,3),2,localbHess);
                 localbHess.resize(2,2);
                 dotu = localbHess * dtheta;
 
                 for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                 {
                     const short_t r = e*Av;
                     for (index_t k = 0; k!= Av; ++k) // for all actives
                     {
                          secDerToHessian(bHessv.col(0).segment(3*k,3),2,localbHess);
                          localbHess.resize(2,2);
                          dotv = localbHess * dtheta;

                          cnvec_var2 = vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                       vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1))); 

                          cnvec_vara_var2 = vecFun(d,dotu(0)).cross(vecFun(e,bGradsv.at(2*k+1))) + 
                                            vecFun(e,dotv(0)).cross(vecFun(d,bGradsu.at(2*j+1))) +
                                            vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,dotv(1))) +
                                            vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,dotu(1)));

                          measure1 = (cnvec.col3d(0).dot(cnvec_vara.col3d(0)))/measure; 
                          measurer = (cnvec.col3d(0).dot(cnvec_var1u.col3d(s+j)))/measure; 
                          measures = (cnvec.col3d(0).dot(cnvec_var1v.col3d(r+k)))/measure; 

                          measure1r = cnvec_var1u.col3d(s+j).dot(cnvec_vara.col3d(0))/measure +
                                      cnvec.col3d(0).dot(cnvec_vara_var1u.col3d(s+j))/measure -
                                      cnvec.col3d(0).dot(cnvec_vara.col3d(0))*measurer/std::pow(measure,2); 
                          
                          measure1s = cnvec_var1v.col3d(r+k).dot(cnvec_vara.col3d(0))/measure +
                                      cnvec.col3d(0).dot(cnvec_vara_var1v.col3d(r+k))/measure -
                                      cnvec.col3d(0).dot(cnvec_vara.col3d(0))*measures/std::pow(measure,2); 

                          measurers = cnvec_var1v.col3d(r+k).dot(cnvec_var1u.col3d(s+j))/measure +
                                      cnvec.col3d(0).dot(cnvec_var2)/measure - 
                                      cnvec.col3d(0).dot(cnvec_var1u.col3d(s+j))*measures/std::pow(measure,2);

                          measure1rs = cnvec_var2.dot(cnvec_vara.col3d(0))/measure + cnvec_var1u.col3d(s+j).dot(cnvec_vara_var1v.col3d(r+k))/
                                       measure - cnvec_var1u.col3d(s+j).dot(cnvec_vara.col3d(0))*measures/std::pow(measure,2) +
                                       cnvec_var1v.col3d(r+k).dot(cnvec_vara_var1u.col3d(s+j))/measure + cnvec.col3d(0).dot(cnvec_vara_var2)/
                                       measure - cnvec.col3d(0).dot(cnvec_vara_var1u.col3d(s+j))*measures/std::pow(measure,2) -
                                       cnvec_var1v.col3d(r+k).dot(cnvec_vara.col3d(0))*measurer/std::pow(measure,2) -
                                       cnvec.col3d(0).dot(cnvec_vara_var1v.col3d(r+k))*measurer/std::pow(measure,2) -
                                       cnvec.col3d(0).dot(cnvec_vara.col3d(0))*measurers/std::pow(measure,2) +
                                       2*measurer*measures*cnvec.col3d(0).dot(cnvec_vara.col3d(0))/std::pow(measure,3); 

                          matrix(s+j,r+k) = cnvec_vara_var2/measure - cnvec_vara_var1u.col3d(s+j)*measures/std::pow(measure,2) -
                                            cnvec_vara_var1v.col3d(r+k)*measurer/std::pow(measure,2) - cnvec_vara*measurers/std::pow(measure,2) +
                                            2*measurer*measures*cnvec_vara/std::pow(measure,3) - cnvec_var2*measure1/std::pow(measure,2) -
                                            cnvec_var1u.col3d(s+j)*measure1s/std::pow(measure,2) + 2*measure1*measures* cnvec_var1u.col3d(s+j)/
                                            std::pow(measure,3) - cnvec_var1v.col3d(r+k)*measure1r/std::pow(measure,2) - cnvec.col3d(0)*measure1rs/
                                            std::pow(measure,2) + 2*measure1r*measures*cnvec.col3d(0)/std::pow(measure,3) + 2*measure1*measurer*
                                            cnvec_var1v.col3d(r+k)/std::pow(measure,3) + 2*measure1s*measurer*cnvec.col3d(0)/std::pow(measure,3) +
                                            2*measure1*measurers*cnvec.col3d(0)/std::pow(measure,3) - 6*measure1*measurer*measures*cnvec.col3d(0)/
                                            std::pow(measure,5);

                          res(s+j,r+k) = matrix(s+j,r+k).dot(evec.col(0));
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
         evList.add(_C);
         _C.data().flags |= NEED_VALUE | NEED_DERIV;
         cnv(_S,_C).parse(evList);
         cnv_vara(_S,_C).parse(evList);     
         cnv_var1(_u,_S,_C).parse(evList);   
         cnv_var1(_v,_S,_C).parse(evList);   
         cnv_vara_var1(_u,_S,_C).parse(evList); 
         cnv_vara_var1(_v,_S,_C).parse(evList);  
         _vec.parse(evList);
     }
 
     const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
     const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
     index_t cardinality_impl() const { return _u.cardinality_impl(); }
 
     void print(std::ostream &os) const { os << "curve_normal_varalong_var2dot("; _u.print(os); os <<")"; }
 };

template<class E1, class E2, class E3>
class curve_binormal_varalong_var2dot_expr : public _expr<curve_binormal_varalong_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_var2dot_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &S,
                                         const gsGeometryMap<Scalar> &C, _expr<E3> const &vec): _u(u), _v(v), _S(S), _C(C), _vec(vec)
                                         {}

    mutable gsMatrix<Scalar> theta, dtheta, bGradsu, bGradsv, bHessu, bHessv, localbHess, dotu, dotv, evec, res;
    mutable gsMatrix<Scalar> ctvec, ctvec_vara, ctvec_var1u, ctvec_var1v, ctvec_vara_var1u, ctvec_vara_var1v;
    mutable gsMatrix<Scalar> cnvec_var1u, cnvec_var1v, cnvec_vara_var1u, cnvec_vara_var1v;
    mutable gsMatrix<Scalar> cbvec, cbvec_vara, cbvec_var1u, cbvec_var1v, cbvec_vara_var1u, cbvec_vara_var1v;
    mutable gsVector<Scalar,3> cnvec_var2, cbvec_var2, cnvec_vara_var2, cbvec_vara_var2;
    mutable gsMatrix<gsVector<Scalar,3>> matrix;
    mutable Scalar measurer, measures, measure1, measure1r, measure1s, measurers, measure1rs;

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
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows(); 
        const index_t Av = mbv.basis(0).active(theta).rows();

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);
        bHessu =  mbu.basis(0).deriv2(theta);
        bHessv =  mbv.basis(0).deriv2(theta);

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_varalong_expr<Scalar> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_tangent_var1_expr<E1> ctv_var1u(_u,_S,_C);    
        ctvec_var1u = ctv_var1u.eval(k).transpose();
        curve_tangent_var1_expr<E2> ctv_var1v(_v,_S,_C);    
        ctvec_var1v = ctv_var1v.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E1> ctv_vara_var1u(_u,_S,_C);
        ctvec_vara_var1u = ctv_vara_var1u.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E2> ctv_vara_var1v(_v,_S,_C);
        ctvec_vara_var1v = ctv_vara_var1v.eval(k).transpose();
        curve_normal_var1_expr<E1> cnv_var1u(_u,_S,_C);
        cnvec_var1u = cnv_var1u.eval(k).transpose();
        curve_normal_var1_expr<E2> cnv_var1v(_v,_S,_C);
        cnvec_var1v = cnv_var1v.eval(k).transpose();
        curve_normal_varalong_var1_expr<E1> cnv_vara_var1u(_u,_S,_C);
        cnvec_vara_var1u = cnv_vara_var1u.eval(k).transpose();
        curve_normal_varalong_var1_expr<E2> cnv_vara_var1v(_v,_S,_C);    
        cnvec_vara_var1v = cnv_vara_var1v.eval(k).transpose();
        curve_binormal_expr<Scalar> cbv(_S,_C);
        cbvec = cbv.eval(k);
        const Scalar measure = cbvec.norm();
        curve_binormal_varalong_expr<Scalar> cbv_vara(_S,_C);
        cbvec_vara = cbv_vara.eval(k);
        curve_binormal_var1_expr<E1> cbv_var1u(_u,_S,_C);
        cbvec_var1u = cbv_var1u.eval(k).transpose();
        curve_binormal_var1_expr<E2> cbv_var1v(_v,_S,_C);
        cbvec_var1v = cbv_var1v.eval(k).transpose();
        curve_binormal_varalong_var1_expr<E1> cbv_vara_var1u(_u,_S,_C);
        cbvec_vara_var1u = cbv_vara_var1u.eval(k).transpose();
        curve_binormal_varalong_var1_expr<E2> cbv_vara_var1v(_u,_S,_C);
        cbvec_vara_var1v = cbv_vara_var1v.eval(k).transpose();

        matrix.resize(Au*_u.dim(), Av*_v.dim());
        res.resize(Au*_u.dim(), Av*_v.dim());

        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                secDerToHessian(bHessu.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dotu = localbHess * dtheta;

                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {
                        secDerToHessian(bHessv.col(0).segment(3*k,3),2,localbHess);
                        localbHess.resize(2,2);
                        dotv = localbHess * dtheta;

                        cnvec_var2 = vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                     vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1)));
                        
                        cnvec_vara_var2 = vecFun(d,dotu(0)).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                          vecFun(e,dotv(0)).cross(vecFun(d,bGradsu.at(2*j+1))) +
                                          vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,dotv(1))) +
                                          vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,dotu(1)));

                        cbvec_var2 = ctvec_var1u.col3d(s+j).cross(cnvec_var1v.col3d(r+k)) +
                                     ctvec_var1v.col3d(r+k).cross(cnvec_var1u.col3d(s+j)) + 
                                     ctvec.col3d(0).cross(cnvec_var2);

                        cbvec_vara_var2 = ctvec_vara_var1u.col3d(s+j).cross(cnvec_var1v.col3d(r+k)) +
                                          ctvec_vara_var1v.col3d(r+k).cross(cnvec_var1u.col3d(s+j)) +
                                          ctvec_vara.col3d(0).cross(cnvec_var2) + ctvec.col3d(0).cross(cnvec_vara_var2) +
                                          ctvec_var1u.col3d(s+j).cross(cnvec_vara_var1v.col3d(r+k)) +
                                          ctvec_var1v.col3d(r+k).cross(cnvec_vara_var1u.col3d(s+j));

                        measure1 = (cbvec.col3d(0).dot(cbvec_vara.col3d(0)))/measure;
                        measurer = (cbvec.col3d(0).dot(cbvec_var1u.col3d(s+j)))/measure;
                        measures = (cbvec.col3d(0).dot(cbvec_var1v.col3d(r+k)))/measure;

                        measure1r = cbvec_var1u.col3d(s+j).dot(cbvec_vara.col3d(0))/measure +
                                    cbvec.col3d(0).dot(cbvec_vara_var1u.col3d(s+j))/measure -
                                    cbvec.col3d(0).dot(cbvec_vara.col3d(0))*measurer/std::pow(measure,2); 

                        measure1s = cbvec_var1v.col3d(r+k).dot(cbvec_vara.col3d(0))/measure +
                                    cbvec.col3d(0).dot(cbvec_vara_var1v.col3d(r+k))/measure -
                                    cbvec.col3d(0).dot(cbvec_vara.col3d(0))*measures/std::pow(measure,2);

                        measurers = cbvec_var1v.col3d(r+k).dot(cbvec_var1u.col3d(s+j))/measure +
                                    cbvec.col3d(0).dot(cbvec_var2)/measure - 
                                    cbvec.col3d(0).dot(cbvec_var1u.col3d(s+j))*measures/std::pow(measure,2);

                        measure1rs = cbvec_var2.dot(cbvec_vara.col3d(0))/measure + cbvec_var1u.col3d(s+j).dot(cbvec_vara_var1v.col3d(r+k))/
                                     measure - cbvec_var1u.col3d(s+j).dot(cbvec_vara.col3d(0))*measures/std::pow(measure,2) +
                                     cbvec_var1v.col3d(r+k).dot(cbvec_vara_var1u.col3d(s+j))/measure + cbvec.col3d(0).dot(cbvec_vara_var2)/
                                     measure - cbvec.col3d(0).dot(cbvec_vara_var1u.col3d(s+j))*measures/std::pow(measure,2) -
                                     cbvec_var1v.col3d(r+k).dot(cbvec_vara.col3d(0))*measurer/std::pow(measure,2) -
                                     cbvec.col3d(0).dot(cbvec_vara_var1v.col3d(r+k))*measurer/std::pow(measure,2) -
                                     cbvec.col3d(0).dot(cbvec_vara.col3d(0))*measurers/std::pow(measure,2) +
                                     2*measurer*measures*cbvec.col3d(0).dot(cbvec_vara.col3d(0))/std::pow(measure,3); 

                        matrix(s+j,r+k) = cbvec_vara_var2/measure - cbvec_vara_var1u.col3d(s+j)*measures/std::pow(measure,2) -
                                          cbvec_vara_var1v.col3d(r+k)*measurer/std::pow(measure,2) - cbvec_vara*measurers/std::pow(measure,2) +
                                          2*measurer*measures*cbvec_vara/std::pow(measure,3) - cbvec_var2*measure1/std::pow(measure,2) -
                                          cbvec_var1u.col3d(s+j)*measure1s/std::pow(measure,2) + 2*measure1*measures* cbvec_var1u.col3d(s+j)/
                                          std::pow(measure,3) - cbvec_var1v.col3d(r+k)*measure1r/std::pow(measure,2) - cbvec.col3d(0)*measure1rs/
                                          std::pow(measure,2) + 2*measure1r*measures*cbvec.col3d(0)/std::pow(measure,3) + 2*measure1*measurer*
                                          cbvec_var1v.col3d(r+k)/std::pow(measure,3) + 2*measure1s*measurer*cbvec.col3d(0)/std::pow(measure,3) +
                                          2*measure1*measurers*cbvec.col3d(0)/std::pow(measure,3) - 6*measure1*measurer*measures*cbvec.col3d(0)/
                                          std::pow(measure,5);

                        res(s+j,r+k) = matrix(s+j,r+k).dot(evec.col(0));
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
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList); 
        ctv_var1(_v,_S,_C).parse(evList);  
        ctv_vara_var1(_u,_S,_C).parse(evList);  
        ctv_vara_var1(_v,_S,_C).parse(evList);     
        cnv_var1(_u,_S,_C).parse(evList);   
        cnv_var1(_v,_S,_C).parse(evList);   
        cnv_vara_var1(_u,_S,_C).parse(evList); 
        cnv_vara_var1(_v,_S,_C).parse(evList);  
        cbv(_S,_C).parse(evList);
        cbv_vara(_S,_C).parse(evList);
        cbv_var1(_u,_S,_C).parse(evList);
        cbv_var1(_v,_S,_C).parse(evList);
        cbv_vara_var1(_u,_S,_C).parse(evList);
        cbv_vara_var1(_v,_S,_C).parse(evList);
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_binormal_varalong_var2dot("; _u.print(os); os <<")"; }
};

template<class E1, class E2, class E3, class E4>
class var1_dot_othervar1_expr : public _expr<var1_dot_othervar1_expr<E1,E2,E3,E4>>
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _that;
    typename E4::Nested_t _other;
    
public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    var1_dot_othervar1_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &C, _expr<E3> const &that, _expr<E4> const &other): 
                            _u(u), _v(v), _C(C), _that(that), _other(other){}

    mutable gsMatrix<Scalar> ethat, eother, theta, res;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        ethat = _that.eval(k).transpose();
        eother = _other.eval(k).transpose();

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows();
        const index_t Av = mbv.basis(0).active(theta).rows();

        res.resize(Au*_u.dim(), Av*_v.dim());

        for (index_t i = 0; i!= Au*_u.dim(); ++i)
        {
            for (index_t j = 0; j!= Av*_v.dim(); ++j)
            {
                res(i,j) = ethat.col3d(i).dot(eother.col3d(j));
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
            evList.add(_C);
            _C.data().flags |= NEED_VALUE;
            _that.parse(evList);
            _other.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _u.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "var1_dot_othervar1("; _u.print(os); os <<")"; }
};

/**
 * @brief Mass matrix.
 */

 template<class E1, class E2>
 class massmatrix_expr : public _expr< massmatrix_expr<E1,E2> >
 {
 public:
     typedef typename E1::Scalar Scalar;
 
 private:
     typename E1::Nested_t _u;
     typename E2::Nested_t _v;
     typename gsGeometryMap<Scalar>::Nested_t _C;
     typename E1::Scalar _density;
 
 public:
     enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};
 
     massmatrix_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &C, const Scalar &density):
                     _u(u), _v(v), _C(C), _density(density) {}
 
     mutable gsMatrix<Scalar> theta, bu, bv, res;
 
 #   define Eigen gsEigen
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 #   undef Eigen
 
     const gsMatrix<Scalar> & eval(const index_t k) const
     {
         theta  = _C.data().values[0].col(k);
 
         const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
         const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
         const index_t Au = mbu.basis(0).active(theta).rows();
         const index_t Av = mbv.basis(0).active(theta).rows();
         bu = mbu.basis(0).eval(theta);
         bv = mbv.basis(0).eval(theta);
         gsDebugVar(bu);
 
         res.resize(Au*_u.dim(), Av*_v.dim());
         res.setConstant(_density);
 
         for (index_t d = 0; d!= _u.dim(); ++d) // for all components
         {
             const short_t s = d*Au;
             for (index_t j = 0; j!= Au; ++j) // for all actives
             {
                 for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                 {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {
                        res(s+j,r+k) *= bu.at(j) * bv.at(k);
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
         evList.add(_C);
         _C.data().flags |= NEED_VALUE | NEED_DERIV;
     }
 
     const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
     const gsFeSpace<Scalar> & colVar() const {return _u.rowVar();}
     index_t cardinality_impl() const { return _u.cardinality_impl(); }
 
     void print(std::ostream &os) const { os << "massmatrix("; _u.print(os); os <<")"; }
 };

 template<class T> EIGEN_STRONG_INLINE
 curve_tangent_expr<T> ctv(const gsGeometryMap<T> &S,
                           const gsGeometryMap<T> &C) { return curve_tangent_expr<T>(S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_tangent_var1_expr<E> ctv_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                            const gsGeometryMap<typename E::Scalar> &C) { return curve_tangent_var1_expr<E>(u,S,C); }
 
 template<class T> EIGEN_STRONG_INLINE
 curve_tangent_varalong_expr<T> ctv_vara(const gsGeometryMap<T> &S,
                                const gsGeometryMap<T> &C) { return curve_tangent_varalong_expr<T>(S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_tangent_varalong_var1_expr<E> ctv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                     const gsGeometryMap<typename E::Scalar> &C) { return curve_tangent_varalong_var1_expr<E>(u,S,C); }
 
 template<class T> EIGEN_STRONG_INLINE
 curve_normal_expr<T>  cnv(const gsGeometryMap<T> &S,
                           const gsGeometryMap<T> &C) { return curve_normal_expr<T>(S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_normal_var1_expr<E> cnv_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                           const gsGeometryMap<typename E::Scalar> &C) { return curve_normal_var1_expr<E>(u,S,C); }
 
 template<class T> EIGEN_STRONG_INLINE
 curve_normal_varalong_expr<T>  cnv_vara(const gsGeometryMap<T> &S, const gsGeometryMap<T> &C) 
                                        { return curve_normal_varalong_expr<T>(S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_varalong_normalized_expr<T> cnv_vara_normalized(const gsGeometryMap<T> &S, const gsGeometryMap<T> &C) 
                                        { return curve_normal_varalong_normalized_expr<T>(S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_normal_varalong_var1_expr<E> cnv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                                  const gsGeometryMap<typename E::Scalar> &C)
                                                { return curve_normal_varalong_var1_expr<E>(u,S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_normal_varalong_var1_normalized_expr<E> cnv_vara_var1_normalized(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                               const gsGeometryMap<typename E::Scalar> &C)
                                               { return curve_normal_varalong_var1_normalized_expr<E>(u,S,C); }
 
 template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
 curve_normal_var2dot_expr<E1,E2,E3> cnv_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &C, const E3 &vec)
                                                { return curve_normal_var2dot_expr<E1,E2,E3>(u,v,C,vec); }
 
 template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
 curve_normal_varalong_var2dot_expr<E1,E2,E3> cnv_vara_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &S,
                                              const gsGeometryMap<typename E1::Scalar> &C, const E3 &vec) 
                                              { return curve_normal_varalong_var2dot_expr<E1,E2,E3>(u,v,S,C,vec); }
 
 template<class T> EIGEN_STRONG_INLINE
 curve_binormal_expr<T>  cbv(const gsGeometryMap<T> &S,
                             const gsGeometryMap<T> &C) { return curve_binormal_expr<T>(S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_binormal_var1_expr<E> cbv_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                             const gsGeometryMap<typename E::Scalar> &C) { return curve_binormal_var1_expr<E>(u,S,C); }
 
 template<class T> EIGEN_STRONG_INLINE
 curve_binormal_varalong_expr<T> cbv_vara(const gsGeometryMap<T> &S,
                                 const gsGeometryMap<T> &C) { return curve_binormal_varalong_expr<T>(S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_varalong_normalized_expr<T> cbv_vara_normalized(const gsGeometryMap<T> &S, const gsGeometryMap<T> &C)
                                                            { return curve_binormal_varalong_normalized_expr<T>(S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_binormal_varalong_var1_expr<E> cbv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                                    const gsGeometryMap<typename E::Scalar> &C)
                                                   { return curve_binormal_varalong_var1_expr<E>(u,S,C); }
 
 template<class E> EIGEN_STRONG_INLINE
 curve_binormal_varalong_var1_normalized_expr<E> cbv_vara_var1_normalized(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                                 const gsGeometryMap<typename E::Scalar> &C)
                                                 { return curve_binormal_varalong_var1_normalized_expr<E>(u,S,C); }
 
 template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
 curve_binormal_var2dot_expr<E1,E2,E3> cbv_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &S,
                                                   const gsGeometryMap<typename E1::Scalar> &C, const E3 &vec) 
                                                  { return curve_binormal_var2dot_expr<E1,E2,E3>(u,v,S,C,vec); }
 
 template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
 curve_binormal_varalong_var2dot_expr<E1,E2,E3> cbv_vara_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &S,
                                                const gsGeometryMap<typename E1::Scalar> &C,
                                                const E3 &vec) { return curve_binormal_varalong_var2dot_expr<E1,E2,E3>(u,v,S,C,vec); }
 
 template<class E1, class E2, class E3, class E4> EIGEN_STRONG_INLINE
 var1_dot_othervar1_expr<E1,E2,E3,E4> var1_dot_var1(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &C,
                                                    const E3 &that, const E4 &other)
                                                   {return var1_dot_othervar1_expr<E1,E2,E3,E4>(u,v,C,that,other);}
 
 template<class E1, class E2> EIGEN_STRONG_INLINE
 massmatrix_expr<E1,E2> massmat(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &C, const typename E1::Scalar &density)
                                                 { return massmatrix_expr<E1,E2>(u,v,C,density); }
}
}