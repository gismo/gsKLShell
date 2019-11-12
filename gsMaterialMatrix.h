/** @file gsMaterialMatrix.h

    @brief Provides material matrices for the thin shell class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

namespace gismo
{

template <class T>
class gsMaterialMatrix : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _YoungsModulus;
    const gsFunction<T> * _PoissonRatio;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<real_t,3,3> F0;
    mutable gsMatrix<T> Emat,Nmat;
    mutable real_t lambda, mu, E, nu, C_constant;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrix(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
                   const gsFunction<T> & PoissonRatio) :
    _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrix() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrix)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrix<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrix(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
        return *_mm_piece;
    }

    //class .. matMatrix_z
    // should contain eval_into(thickness variable)

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points = u;

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _YoungsModulus->eval_into(_tmp.values[0], Emat);
        _PoissonRatio->eval_into(_tmp.values[0], Nmat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

            F0.leftCols(2) = _tmp.jacobian(i);
            F0.col(2)      = _tmp.normal(i).normalized();
            F0 = F0.inverse();
            F0 = F0 * F0.transpose(); //3x3

            // Evaluate material properties on the quadrature point
            E = Emat(0,i);
            nu = Nmat(0,i);
            lambda = E * nu / ( (1. + nu)*(1.-2.*nu)) ;
            mu     = E / (2.*(1. + nu)) ;

            C_constant = 4*lambda*mu/(lambda+2*mu);

            C(0,0) = C_constant*F0(0,0)*F0(0,0) + 2*mu*(2*F0(0,0)*F0(0,0));
            C(1,1) = C_constant*F0(1,1)*F0(1,1) + 2*mu*(2*F0(1,1)*F0(1,1));
            C(2,2) = C_constant*F0(0,1)*F0(0,1) + 2*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
            C(1,0) =
            C(0,1) = C_constant*F0(0,0)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(0,1));
            C(2,0) =
            C(0,2) = C_constant*F0(0,0)*F0(0,1) + 2*mu*(2*F0(0,0)*F0(0,1));
            C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(1,1));

            //gsDebugVar(C);
        }
    }

    // piece(k) --> for patch k

};

} //namespace gismo
