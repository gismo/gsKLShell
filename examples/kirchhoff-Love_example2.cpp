/** @file kirchhoff-Love_example.cpp

    @brief Solver for kirchhoff-Love shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst & A. Mantzaflaris
*/

//! [Include namespace]

#include <gismo.h>
#include <gsKLShell/src/gsThinShellUtils.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>

namespace gismo{
namespace expr{

/**
 * @brief      Second variation of the surface normal times a vector.
 *
 * @tparam     E1    Type of u
 * @tparam     E2    Type of v
 * @tparam     E3    Type of the vector
 */
template<class E1, class E2, class E3>
class var2mod_expr : public _expr<var2mod_expr<E1,E2,E3> >
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

    var2mod_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& Ef) : _u(u),_v(v), _G(G), _Ef(Ef) { }

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
        const Scalar measure =  _G.data().measures.at(k);

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
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
        _Ef.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const { os << "var2("; _u.print(os); os <<")"; }
};

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
var2mod_expr<E1,E2,E3> var2mod(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & Ef)
{ return var2mod_expr<E1,E2,E3>(u,v, G, Ef); }

} // namespace expr
} // namespace gismo

//! [Include namespace]
using namespace gismo;
//! [Include namespace]

// Input is parametric coordinates of the surface \a mp
template <class T>
class gsMaterialMatrixD : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _YoungsModulus;
    const gsFunction<T> * _PoissonRatio;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<T,3,3> F0;
    mutable gsMatrix<T> Emat,Nmat;
    mutable T lambda, mu, E, nu, C_constant;

    mutable std::vector<gsMaterialMatrixD> m_pieces;
public:
    /// Shared pointer for gsMaterialMatrixD
    typedef memory::shared_ptr< gsMaterialMatrixD > Ptr;

    /// Unique pointer for gsMaterialMatrixD
    typedef memory::unique_ptr< gsMaterialMatrixD > uPtr;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    gsMaterialMatrixD() { }

    gsMaterialMatrixD(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
                   const gsFunction<T> & PoissonRatio) :
    _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixD() { }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixD)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    const gsFunction<T> & piece(const index_t k) const
    {
#       pragma omp critical
        if (m_pieces.empty())
        {
            m_pieces.resize(_mp->nPieces());
            for (index_t k = 0; k!=_mp->nPieces(); ++k)
                m_pieces[k] = gsMaterialMatrixD(_mp->piece(k),
                                               _YoungsModulus->piece(k),
                                               _PoissonRatio->piece(k) );
        }
        return m_pieces[k];
    }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        #pragma omp critical
        {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points = u;

        static_cast<const gsFunction<T>&>(_mp->piece(0)).computeMap(_tmp); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

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

            C_constant = 2*lambda*mu/(lambda+2*mu);

            C(0,0) = C_constant*F0(0,0)*F0(0,0) + 1*mu*(2*F0(0,0)*F0(0,0));
            C(1,1) = C_constant*F0(1,1)*F0(1,1) + 1*mu*(2*F0(1,1)*F0(1,1));
            C(2,2) = C_constant*F0(0,1)*F0(0,1) + 1*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
            C(1,0) =
            C(0,1) = C_constant*F0(0,0)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(0,1));
            C(2,0) =
            C(0,2) = C_constant*F0(0,0)*F0(0,1) + 1*mu*(2*F0(0,0)*F0(0,1));
            C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(1,1));
        }
        }//end omp critical
    }

};

template <class T>
class Test : public gismo::gsFunction<T>
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
public:
    Test( gsMultiPatch<T> & mp,
            gsMultiBasis<T> & mb,
            gsFunction<T> & thickness,
            gsFunction<T> & force,
            gsMaterialMatrixD<T> & materialmatrix,
            gsBoundaryConditions<T> & bc,
            T alpha_d,
            T alpha_r,
            index_t continuity = -1,
            index_t verbose = 0                    )
    :
    m_mp(mp),
    m_mb(mb),
    m_thick(&thickness),
    m_force(&force),
    m_mm(materialmatrix),
    m_bcs(bc),
    m_alpha_d(alpha_d),
    m_alpha_r(alpha_r),
    m_continuity(continuity),
    m_verbose(verbose)
    {
        //! [Assembler setup]
        m_assembler = gsExprAssembler<T>(1,1);

        // Elements used for numerical integration
        m_assembler.setIntegrationElements(m_mb);
        m_evaluator = gsExprEvaluator<T>(m_assembler);

        // Set the discretization space
        space u = m_assembler.getSpace(m_mb, 3);
        gsBoundaryConditions<T> m_bcs_empty;
        u.setup(m_bcs, dirichlet::interpolation, -1);
        // u.setup();
        m_assembler.initSystem();
        gsDebugVar(m_assembler.numDofs());
    }

    short_t domainDim() const
    { return m_assembler.numDofs(); }

    short_t targetDim() const
    { return 3; }

    /// Evaluates the non-zero spline functions at value u.
    void eval_into(const gsMatrix<T> & solVec, gsMatrix<T>& result) const
    {
        result.resize(this->targetDim(),solVec.cols());
        space u = m_assembler.trialSpace(0);
        gsMultiPatch<T> m_def;
        // Set the geometry map
        // geometryMap G = m_assembler.getMap(m_mp);
        geometryMap defG = m_assembler.getMap(m_def);

        for (index_t p = 0; p!=solVec.cols(); p++)
        {
            gsMatrix<T> sol = solVec.col(p);
            solution u_sol = m_assembler.getSolution(u,sol);

            gsMatrix<T> cc;
            m_def = m_mp;
            for ( size_t k =0; k!=m_mp.nPatches(); ++k) // Deform the geometry
            {
                u_sol.extract(cc, k);
                m_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
            }


            gsVector<T> pt(2);
            pt<<0.5,1.0;
            result.col(p) = m_evaluator.eval(unv(defG),pt,0);
        }
    }

    // void deriv_into(const gsMatrix<T> & solVec, gsMatrix<T>& result) const
    // {

    // }

protected:
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsFunction<T> * m_thick;
    gsFunction<T> * m_force;
    gsMaterialMatrixD<T> m_mm;
    gsBoundaryConditions<T> m_bcs;
    T m_alpha_d, m_alpha_r;
    index_t m_continuity, m_verbose;


    mutable gsExprAssembler<T> m_assembler;
    mutable gsExprEvaluator<T> m_evaluator;
};


template <class T>
class Residual : public gismo::gsFunction<T>
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
public:
    Residual(   gsMultiPatch<T> & mp,
                gsMultiBasis<T> & mb,
                gsFunction<T> & thickness,
                gsFunction<T> & force,
                gsMaterialMatrixD<T> & materialmatrix,
                gsBoundaryConditions<T> & bc,
                T alpha_d,
                T alpha_r,
                index_t continuity = -1,
                index_t verbose = 0                    )
    :
    m_mp(mp),
    m_mb(mb),
    m_thick(&thickness),
    m_force(&force),
    m_mm(materialmatrix),
    m_bcs(bc),
    m_alpha_d(alpha_d),
    m_alpha_r(alpha_r),
    m_continuity(continuity),
    m_verbose(verbose)
    {
        //! [Assembler setup]
        m_assembler = gsExprAssembler<T>(1,1);

        // Elements used for numerical integration
        m_assembler.setIntegrationElements(m_mb);
        m_evaluator = gsExprEvaluator<T>(m_assembler);

        // Set the discretization space
        space u = m_assembler.getSpace(m_mb, 3);
        u.setup(m_bcs, dirichlet::interpolation, -1);
        m_assembler.initSystem();
    }

    short_t domainDim() const
    { return m_assembler.numDofs(); }

    short_t targetDim() const
    { return m_assembler.numDofs(); }

    /// Evaluates the non-zero spline functions at value u.
    void eval_into(const gsMatrix<T> & solVec, gsMatrix<T>& result) const
    {
        result.resize(this->targetDim(),solVec.cols());
        space u = m_assembler.trialSpace(0);
        gsMultiPatch<T> m_def;
        // Set the geometry map
        geometryMap G = m_assembler.getMap(m_mp); // the last map counts
        geometryMap defG = m_assembler.getMap(m_def);

        gsDebugVar(m_mp);

        auto mm = m_assembler.getCoeff(m_mm); // evaluates in the parametric domain, but the class transforms E and nu to physical

        auto tt = m_assembler.getCoeff(*m_thick, G); // evaluates in the physical domain

        gsFunctionExpr<T> mult2t("1","0","0","0","1","0","0","0","2",3);
        auto m2 = m_assembler.getCoeff(mult2t, G); // evaluates in the physical domain

        auto ff = m_assembler.getCoeff(*m_force,G); // evaluates in the physical domain

        // Membrane components
        auto E_m = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) ;
        auto S_m = E_m * reshape(mm,3,3);
        auto N   = tt.val() * S_m;

        auto E_m_der = flat( jac(defG).tr() * jac(u) ) ;
        auto S_m_der = E_m_der * reshape(mm,3,3);
        auto N_der   = tt.val() * S_m_der;

        auto E_m_der2 = flatdot( jac(u),jac(u).tr(), N );

        // Flexural components
        auto E_f = ( deriv2(G,usn(G).tr()) - deriv2(defG,usn(defG).tr()) ) * reshape(m2,3,3) ;
        auto S_f = E_f * reshape(mm,3,3);
        auto M   = tt.val() * tt.val() * tt.val() / 12.0 * S_f;

        auto E_f_der = -( deriv2(u,usn(defG).tr() ) + deriv2(defG,var1(u,defG) ) ) * reshape(m2,3,3);
        auto S_f_der = E_f_der * reshape(mm,3,3);
        auto M_der   = tt.val() * tt.val() * tt.val() / 12.0 * S_f_der;

        auto E_f_der2 = - (flatdot2( deriv2(u), var1(u,defG).tr(), M  ).symmetrize() + var2(u,u,defG, M ));

        auto F        = ff;

        auto That   = cartcon(G);
        auto Ttilde = cartcov(G);
        auto E_m_plot = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) * That;
        auto S_m_plot = E_m_plot * reshape(mm,3,3) * Ttilde;

        // // For Neumann (same for Dirichlet/Nitsche) conditions
        auto g_N = m_assembler.getBdrFunction(G);
        // auto g_N = ff;

        auto du = ((defG.left()-G.left()) - (defG.right()-G.right()));

        auto dN_lr1= (usn(defG.left()).tr()*usn(defG.right())).val();
        auto dN_lr2= (usn(   G.left()).tr()*usn(   G.right())).val();
        auto dN_lr = (usn(defG.left()).tr()*usn(defG.right())
                        - usn(G.left()).tr()*usn(G.right())).val();

        auto dN_rl = (usn(defG.right()).tr()*usn(defG.left())
                        - usn(G.right()).tr()*usn(G.left())).val();

        auto dnN_lr1= (unv(defG.left()).tr()*usn(defG.right())).val();
        auto dnN_lr2= (unv(   G.left()).tr()*usn(   G.right())).val();
        auto dnN_lr= (unv(defG.left()).tr()*usn(defG.right())
                        - unv(G.left()).tr()*usn(G.right())).val();

        auto dnN_rl= (unv(defG.right()).tr()*usn(defG.left())
                        - unv(G.right()).tr()*usn(G.left())).val();

        for (index_t p = 0; p!=solVec.cols(); p++)
        {
            gsMatrix<T> sol = solVec.col(p);
            solution u_sol = m_assembler.getSolution(u,sol);

            gsMatrix<T> cc;
            m_def = m_mp;
            for ( size_t k =0; k!=m_mp.nPatches(); ++k) // Deform the geometry
            {
                u_sol.extract(cc, k);
                m_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
            }

            m_assembler.initSystem();

            // RHS is minussed!
            if (m_continuity > -1)
            m_assembler.assembleIfc(m_mb.topology().interfaces(),
                            m_alpha_d * u.left() * du
                            ,
                            -m_alpha_d * u.right()* du
                             );

            // Penalty of out-of-plane coupling
            // dW^pr / du_r --> first line
            if (m_continuity > 0)
            m_assembler.assembleIfc(m_mb.topology().interfaces(),
                            m_alpha_r * dN_lr * var1(u.left(),defG.left())   * usn(defG.right())
                            ,
                            m_alpha_r * dN_lr * var1(u.right(),defG.right()) * usn(defG.left() )
                            ,
                            // Symmetry
                            m_alpha_r * dN_rl * var1(u.right(),defG.right()) * usn(defG.left())
                            ,
                            m_alpha_r * dN_rl * var1(u.left(),defG.left()) * usn(defG.right() )
                             );

            // Penalty of in-plane coupling
            // dW^pr / du_r --> second line
            //
            // This term seems problematic!!
            if (m_continuity > 0)
            m_assembler.assembleIfc(m_mb.topology().interfaces(),
                              m_alpha_r * dnN_lr * ovar1(u.left() ,defG.left() ) * usn(defG.right()) //nonzero [0, 0 ; BA, BB] -- term 1
                            ,
                              m_alpha_r * dnN_lr *  var1(u.left() ,defG.left() ) * unv(defG.right()) //nonzero [0, 0 ; BA, BB] -- term 2
                            ,
                            // Symmetry
                              m_alpha_r * dnN_rl * ovar1(u.right(),defG.right()) * usn(defG.left()) //nonzero [AA, AB ; 0, 0] -- term 3
                            ,
                              m_alpha_r * dnN_rl *  var1(u.right(),defG.right()) * unv(defG.left() ) //nonzero [-AA, -AB ; 0, 0] -- term 4
                             );

            // Weak Dirichlet term
            // RHS is minussed!
            m_assembler.assembleBdr(
                m_bcs.get("Weak Dirichlet")
                ,
                m_alpha_d * (u * (defG - G) - u * (g_N) )
            );

            // Weak Clamping term
            // RHS is minussed!
            m_assembler.assembleBdr(
                m_bcs.get("Weak Clamped")
                ,
                m_alpha_r * (usn(defG).tr() * nv(G) - usn(G).tr() * nv(G)).val() * var1(u,defG) * nv(G) * tv(G).norm()
            );

            // Neumann term
            m_assembler.assembleBdr(m_bcs.get("Neumann"), u * g_N * tv(G).norm() );

            m_assembler.assemble(
                -u * F * meas(G) + ( ( N * E_m_der.tr() + M * E_f_der.tr() ) * meas(G) ).tr() //+ pressure * u * usn(defG) * meas(G)
                );

            result.col(p) = m_assembler.rhs();
        }
    }
protected:
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsFunction<T> * m_thick;
    gsFunction<T> * m_force;
    gsMaterialMatrixD<T> m_mm;
    gsBoundaryConditions<T> m_bcs;
    T m_alpha_d, m_alpha_r;
    index_t m_continuity, m_verbose;


    mutable gsExprAssembler<T> m_assembler;
    mutable gsExprEvaluator<T> m_evaluator;
};

//! [Forward declarations]
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);
//! [Forward declarations]

template <class T>
gsMultiPatch<T> Panel(T Lp, T Wp, T Hw, T Wf, T x = 0, T y = 0, T z = 0);

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool weak = false;
    bool nonlinear = false;
    index_t verbose = 0;
    index_t continuity = -1;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1.0;

    real_t ALPHA_D = 1e3;
    real_t ALPHA_R = 1e0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addReal("D", "alphaD","alphaD",ALPHA_D);
    cmd.addReal("R", "alphaR","alphaR",ALPHA_R);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "v", "verbose", "Verbosity: 0 = nothing, 1 = iterations, 2 = matrix and vector", verbose);
    cmd.addInt( "c", "continuity", "C^x continuity", continuity);
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("weak", "Weak BCs", weak);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [set test case data]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    if (testCase < 4)
    {
        if (testCase==1)
            fn = "planar/two_squares.xml";
        else if (testCase==2)
            fn = "multipatches/T-beam.xml";
        else if (testCase==3)
            fn = "multipatches/I-beam.xml";
        gsReadFile<>(fn, mp);

    }
    else if (testCase==4)
        mp = Panel(10.,5.,1.,1.);

    gsWrite(mp,"mp_xml");

    PoissonRatio = 0.0;
    if (testCase==1)
        mp.embed(3);

    //! [Refinement]
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview(mp,"mp",1000,true);
    gsWriteParaview(mp.basis(0),"basis",1000,true);

    gsMultiBasis<> dbasis(mp);

    if (verbose>0) gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    if (verbose>0) gsInfo << dbasis.basis(0)<<"\n";
    //! [Refinement]

    //! [pre-define boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    real_t pressure = 0.0;
    //! [pre-define boundary conditions]

    gsPiecewiseFunction<> force(mp.nPatches());
    gsPiecewiseFunction<> t(mp.nPatches());
    //! [Boundary condition case 1]
    if (testCase == 1)
    {
        mp.clearTopology();
        mp.computeTopology();
        // for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        //     for (index_t d = 0; d!=3; d++)
        //         bc.addCondition(bit->patch, bit->side(), condition_type::dirichlet, 0, 0, false, d);

        if (weak)
        {
            bc.addCondition(0, boundary::east, condition_type::weak_dirichlet, 0, 0, false, -1);
            bc.addCondition(1, boundary::west, condition_type::weak_dirichlet, 0, 0, false, -1);
        }
        else
        {
            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1, boundary::west, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0, boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1, boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            }
        }


        // if (weak)
        // {
        //     bc.addCondition(0, boundary::east, condition_type::weak_dirichlet, 0, 0, false, -1);
        //     bc.addCondition(1, boundary::west, condition_type::weak_dirichlet, 0, 0, false, -1);
        // }
        // else
        // {
        // for (index_t d = 0; d!=3; d++)
        //     {
        //         bc.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);
        //         bc.addCondition(1, boundary::west, condition_type::dirichlet, 0, 0, false, d);
        //     }
        // }
        // if (weak)
        //     bc.addCondition(1, boundary::west, condition_type::weak_clamped, 0, 0, false, 2);
        // else
        //     bc.addCondition(1, boundary::west, condition_type::clamped, 0, 0, false, 2);

        tmp << 0,0,-1;
        gsConstantFunction<> piece0(tmp,3);
        force.addPiece(piece0);
        tmp << 0,0,-1;
        gsConstantFunction<> piece1(tmp,3);
        force.addPiece(piece1);

        gsConstantFunction<> tpiece0(1,3);
        t.addPiece(tpiece0);
        gsConstantFunction<> tpiece1(1,3);
        t.addPiece(tpiece1);
    }
    else if (testCase == 2)
    {
        mp.clearTopology();
        mp.computeTopology();
        for (index_t d = 0; d!=3; d++)
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::dirichlet, 0, 0, false, d);
        if (weak)
        {
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::weak_clamped, 0, 0, false, 2);
        }
        else
        {
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::clamped, 0, 0, false, 2);
        }

        PoissonRatio = 0.3;

        tmp << 0,0,0;
        gsConstantFunction<> piece0(tmp,3);
        force.addPiece(piece0);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> piece1(tmp,3);
        force.addPiece(piece1);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> piece2(tmp,3);
        force.addPiece(piece2);

        gsConstantFunction<> tpiece0(1,3);
        t.addPiece(tpiece0);
        gsConstantFunction<> tpiece1(1,3);
        t.addPiece(tpiece1);
        gsConstantFunction<> tpiece2(1,3);
        t.addPiece(tpiece2);
    }
    else if (testCase == 3)
    {
        mp.clearTopology();
        mp.computeTopology();
        gsVector<> neu(3);
        neu<<0,0,1;
        gsConstantFunction<> neuData(neu,3);
        for (index_t d = 0; d!=3; d++)
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::dirichlet, 0, 0, false, d);
        if (weak)
        {
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::weak_clamped, 0, 0, false, 2);
        }
        else
        {
            for (size_t p=0; p!=mp.nPatches(); ++p)
                bc.addCondition(p, boundary::east, condition_type::clamped, 0, 0, false, 2);
        }

        /// DO NOT WORK!!!
        // bc.addCondition(0, boundary::west, condition_type::neumann, &neuData, 0, 0);
        // bc.addCondition(1, boundary::west, condition_type::neumann, &neuData);
        // bc.addCondition(2, boundary::west, condition_type::neumann, &neuData);


        // if (weak)
        //     bc.addCondition(1, boundary::west, condition_type::weak_clamped, 0, 0, false, 2);
        // else
        //     bc.addCondition(1, boundary::west, condition_type::clamped, 0, 0, false, 2);

        tmp << 0,0,0;
        gsConstantFunction<> piece0(tmp,3);
        force.addPiece(piece0);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> piece1(tmp,3);
        force.addPiece(piece1);
        tmp << 0,0,-1e-3;
        gsConstantFunction<> piece2(tmp,3);
        force.addPiece(piece2);
        tmp << 0,0,0;
        gsConstantFunction<> piece3(tmp,3);
        force.addPiece(piece3);
        tmp << 0,0,0;
        gsConstantFunction<> piece4(tmp,3);
        force.addPiece(piece4);

        gsConstantFunction<> tpiece0(1,3);
        t.addPiece(tpiece0);
        gsConstantFunction<> tpiece1(1,3);
        t.addPiece(tpiece1);
        gsConstantFunction<> tpiece2(1,3);
        t.addPiece(tpiece2);
        gsConstantFunction<> tpiece3(1,3);
        t.addPiece(tpiece3);
        gsConstantFunction<> tpiece4(1,3);
        t.addPiece(tpiece4);
    }
    else if (testCase == 4)
    {
        mp.clearTopology();
        mp.computeTopology();
        for (index_t d = 0; d!=3; d++)
        {
            bc.addCondition(0, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(1, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::south, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(2, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(2, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(3, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(3, boundary::south, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(4, boundary::north, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(4, boundary::south, condition_type::dirichlet, 0, 0, false, d);

        }
        // if (weak)
        // {
        //     for (size_t p=0; p!=mp.nPatches(); ++p)
        //         bc.addCondition(p, boundary::east, condition_type::weak_clamped, 0, 0, false, 2);
        // }
        // else
        // {
        //     for (size_t p=0; p!=mp.nPatches(); ++p)
        //         bc.addCondition(p, boundary::east, condition_type::clamped, 0, 0, false, 2);
        // }

        tmp << 0,0,1e-3;
        gsConstantFunction<> piece0(tmp,3);
        force.addPiece(piece0);
        tmp << 0,0,1e-3;
        gsConstantFunction<> piece1(tmp,3);
        force.addPiece(piece1);
        tmp << 0,0,0;
        gsConstantFunction<> piece2(tmp,3);
        force.addPiece(piece2);
        tmp << 0,0,0;
        gsConstantFunction<> piece3(tmp,3);
        force.addPiece(piece3);
        tmp << 0,0,0;
        gsConstantFunction<> piece4(tmp,3);
        force.addPiece(piece4);


        gsConstantFunction<> tpiece0(1,3);
        t.addPiece(tpiece0);
        gsConstantFunction<> tpiece1(1,3);
        t.addPiece(tpiece1);
        gsConstantFunction<> tpiece2(0.5,3);
        t.addPiece(tpiece2);
        gsConstantFunction<> tpiece3(0.1,3);
        t.addPiece(tpiece3);
        gsConstantFunction<> tpiece4(0.1,3);
        t.addPiece(tpiece4);
    }

    //! [Assembler setup]
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    typedef gsExprAssembler<>::element     element;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp); // the last map counts
    geometryMap defG = A.getMap(mp_def);


    // Set the discretization space
    space u = A.getSpace(dbasis, 3);

    // Solution vector and solution variable
    gsMatrix<> random;
    solution u_sol = A.getSolution(u,random);

    // gsFunctionExpr<> materialMat("1","0","0","0","1","0","0","0","1",3);
    gsFunctionExpr<> E(util::to_string(E_modulus),3);
    gsFunctionExpr<> nu(util::to_string(PoissonRatio),3);
    gsMaterialMatrixD<real_t> materialMat(mp, E, nu);
    auto mm = A.getCoeff(materialMat); // evaluates in the parametric domain, but the class transforms E and nu to physical

    auto tt = A.getCoeff(t, G); // evaluates in the physical domain

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    auto m2 = A.getCoeff(mult2t, G); // evaluates in the physical domain

    // gsFunctionExpr<> force("0","0","x",3);
    auto ff = A.getCoeff(force,G); // evaluates in the physical domain

    gsSparseSolver<>::CGDiagonal solver;

    //! [Assembler setup]

    //! [System assembly]

    // Initialize the system
    u.setup(bc, dirichlet::interpolation, -1);

    A.initSystem();

    if (verbose>0) gsInfo<<"Number of degrees of freedom: "<< A.numDofs() <<"\n"<<std::flush;

    /*
        We provide the following functions:
            E_m         membrane strain tensor.
            E_m_der     first variation of E_m.
            E_m_der2    second variation of E_m MULTIPLIED BY S_m.
            E_f         flexural strain tensor.
            E_f_der     second variation of E_f.
            E_f_der2    second variation of E_f MULTIPLIED BY S_f.

        Where:
            G       the undeformed geometry,
            defG    the deformed geometry,
            mm      the material matrix,
            m2      an auxillary matrix to multiply the last row of a tensor with 2
    **/

    // Membrane components
    auto E_m = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) ;
    auto S_m = E_m * reshape(mm,3,3);
    auto N   = tt.val() * S_m;

    auto E_m_der = flat( jac(defG).tr() * jac(u) ) ;
    auto S_m_der = E_m_der * reshape(mm,3,3);
    auto N_der   = tt.val() * S_m_der;

    auto E_m_der2 = flatdot( jac(u),jac(u).tr(), N );

    // Flexural components
    auto E_f = ( deriv2(G,usn(G).tr()) - deriv2(defG,usn(defG).tr()) ) * reshape(m2,3,3) ;
    auto S_f = E_f * reshape(mm,3,3);
    auto M   = tt.val() * tt.val() * tt.val() / 12.0 * S_f;

    auto E_f_der = -( deriv2(u,usn(defG).tr() ) + deriv2(defG,var1(u,defG) ) ) * reshape(m2,3,3);
    auto S_f_der = E_f_der * reshape(mm,3,3);
    auto M_der   = tt.val() * tt.val() * tt.val() / 12.0 * S_f_der;

    auto E_f_der2 = - (flatdot2( deriv2(u), var1(u,defG).tr(), M  ).symmetrize() + var2(u,u,defG, M ));

    auto F        = ff;

    auto That   = cartcon(G);
    auto Ttilde = cartcov(G);
    auto E_m_plot = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) * That;
    auto S_m_plot = E_m_plot * reshape(mm,3,3) * Ttilde;

    // // For Neumann (same for Dirichlet/Nitsche) conditions
    auto g_N = A.getBdrFunction(G);
    // auto g_N = ff;



    auto mmA = tt.val()*mm;
    auto mmD = tt.val() * tt.val() * tt.val() / 12.0 * mm;

    auto cart2cov = gismo::expr::cartcov(G);
    auto con2cartI = gismo::expr::cartcon(G).inv();

    auto mmAcart = (con2cartI * reshape(mmA,3,3) * cart2cov);
    auto mmDcart = (con2cartI * reshape(mmD,3,3) * cart2cov);


    element el   = A.getElement();
    auto h       = (el.area(G.left()) + el.area(G.right())) / 2;
    auto alpha_d = ALPHA_D * reshape(mmAcart,9,1).max().val() / h;
    auto alpha_r = ALPHA_R * reshape(mmDcart,9,1).max().val() / h;

    auto du = ((defG.left()-G.left()) - (defG.right()-G.right()));

    if (continuity > -1)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_d * u.left() * du
                    ,
                    -alpha_d * u.right()* du
                    ,
                     alpha_d * u.left() * u.left().tr()
                    ,
                    -alpha_d * u.right()* u.left() .tr()
                    ,
                    -alpha_d * u.left() * u.right().tr()
                    ,
                     alpha_d * u.right()* u.right().tr()
                     );

    if (verbose>1) gsDebugVar(A.matrix().toDense());
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    auto dN_lr1= (usn(defG.left()).tr()*usn(defG.right())).val();
    auto dN_lr2= (usn(   G.left()).tr()*usn(   G.right())).val();
    auto dN_lr = (usn(defG.left()).tr()*usn(defG.right())
                    - usn(G.left()).tr()*usn(G.right())).val();

    auto dN_rl = (usn(defG.right()).tr()*usn(defG.left())
                    - usn(G.right()).tr()*usn(G.left())).val();

    auto dnN_lr1= (unv(defG.left()).tr()*usn(defG.right())).val();
    auto dnN_lr2= (unv(   G.left()).tr()*usn(   G.right())).val();
    auto dnN_lr= (unv(defG.left()).tr()*usn(defG.right())
                    - unv(G.left()).tr()*usn(G.right())).val();

    auto dnN_rl= (unv(defG.right()).tr()*usn(defG.left())
                    - unv(G.right()).tr()*usn(G.left())).val();

    // Penalty of out-of-plane coupling
    // dW^pr / du_r --> first line
    if (continuity > 0)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_r * dN_lr * var1(u.left(),defG.left())   * usn(defG.right())
                    ,
                     alpha_r * dN_lr * var1(u.right(),defG.right()) * usn(defG.left() )
                    ,
                    // Symmetry
                     alpha_r * dN_rl * var1(u.right(),defG.right())   * usn(defG.left())
                    ,
                     alpha_r * dN_rl * var1(u.left(),defG.left()) * usn(defG.right() )
                     );

    if (verbose>1) gsDebugVar(A.matrix().toDense());
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    // Penalty of in-plane coupling
    // dW^pr / du_r --> second line
    if (continuity > 0)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_r * dnN_lr* ovar1(u.left(),defG.left())  * usn(defG.right())
                    ,
                     alpha_r * dnN_lr* var1(u.right(),defG.right()) * unv(defG.left() )
                    ,
                    // Symmetry
                     alpha_r * dnN_rl* ovar1(u.right(),defG.right())  * usn(defG.left())
                    ,
                     alpha_r * dnN_rl* var1(u.left(),defG.left()) * unv(defG.right() )
                     );

    if (verbose>1) gsDebugVar(A.matrix().toDense());
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    // Penalty of out-of-plane coupling
    // dW^pr / du_r --> first line
    if (continuity > 0)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_r * dN_lr * var2(u.left() ,u.left() ,defG.left() ,usn(defG.right()).tr() )      // left left
                    ,
                     alpha_r * dN_lr * ( var1(u.left() ,defG.left() ) * var1(u.right(),defG.right()).tr() )// left right
                    ,
                     alpha_r * dN_lr * ( var1(u.right(),defG.right()) * var1(u.left() ,defG.left() ).tr() )// right left
                    ,
                     alpha_r * dN_lr * var2( u.right(),u.right(),defG.right(),usn(defG.left() ).tr() )     // right right
                    ,
                    // Symmetry
                     alpha_r * dN_rl * var2(u.right() ,u.right() ,defG.right() ,usn(defG.left()).tr() )      // right right
                    ,
                     alpha_r * dN_rl * ( var1(u.right() ,defG.right() ) * var1(u.left(),defG.left()).tr() )// right left
                    ,
                     alpha_r * dN_rl * ( var1(u.left(),defG.left()) * var1(u.right() ,defG.right() ).tr() )// left right
                    ,
                     alpha_r * dN_rl * var2( u.left(),u.left(),defG.left(),usn(defG.right() ).tr() )     // left left
                     );

    if (verbose>1) gsDebugVar(A.matrix().toDense());
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    // Penalty of out-of-plane coupling
    // dW^pr / du_r --> second line
    if (continuity > 0)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()    // left left
                    ,
                     alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()   // left right
                    ,
                     alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()   // right left
                    ,
                     alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()  // right right
                    ,
                    // Symmetry
                     alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()    // right right
                    ,
                     alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()   // right left
                    ,
                     alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()   // left right
                    ,
                     alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()  // left left
                     );

    if (verbose>1) gsDebugVar(A.matrix().toDense());
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    // Penalty of in-plane coupling
    // dW^pr / du_r --> third line
    if (continuity > 0)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_r * dnN_lr * ovar2(u.left(),u.left(),defG.left(),usn(defG.right()).tr()) // left left
                    + // Symmetry
                     alpha_r * dnN_rl * ovar2(u.left(),u.left(),defG.left(),usn(defG.right()).tr()) // left left
                    ,
                     alpha_r * dnN_lr * ( ovar1(u.left() ,defG.left() ) * var1(u.right(),defG.right()).tr() ) // left right
                    + // Symmetry
                     alpha_r * dnN_rl * ( ovar1(u.left() ,defG.left() ) * var1(u.right(),defG.right()).tr() ) // right left
                    ,
                     alpha_r * dnN_lr * ( ovar1(u.right(),defG.right()) * var1(u.left() ,defG.left() ).tr() ) // right left
                    + // Symmetry
                     alpha_r * dnN_rl * ( ovar1(u.right(),defG.right()) * var1(u.left() ,defG.left() ).tr() ) // right left
                    ,
                     alpha_r * dnN_lr * ovar2(u.right(),u.right(),defG.right(),usn(defG.left()).tr()) // right right
                    + // Symmetry
                     alpha_r * dnN_rl * ovar2(u.right(),u.right(),defG.right(),usn(defG.left()).tr()) // right right
                     );

    if (verbose>1) gsDebugVar(A.matrix().toDense());
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    // Penalty of in-plane coupling
    // dW^pr / du_r --> fourth line
    if (continuity > 0)
    A.assembleIfc(dbasis.topology().interfaces(),
                     alpha_r * ( ovar1(u.left() ,defG.left() ) * usn(defG.right()) ) * ( ovar1(u.left() ,defG.left() ) * usn(defG.right()) ).tr() // left left
                    + // Symmetry
                     alpha_r * (  var1(u.left() ,defG.left() ) * unv(defG.right()) ) * (  var1(u.left() ,defG.left() ) * unv(defG.right()) ).tr() // left left
                    ,
                     alpha_r * ( ovar1(u.left() ,defG.left() ) * usn(defG.right()) ) * (  var1(u.right(),defG.right()) * unv(defG.left() ) ).tr() // left right
                    + // Symmetry
                     alpha_r * (  var1(u.left() ,defG.left() ) * unv(defG.right()) ) * ( ovar1(u.right(),defG.right()) * usn(defG.left() ) ).tr() // left right
                    ,
                     alpha_r * (  var1(u.right(),defG.right()) * unv(defG.left() ) ) * ( ovar1(u.left() ,defG.left() ) * usn(defG.right()) ).tr() // right left
                    + // Symmetry
                     alpha_r * ( ovar1(u.right(),defG.right()) * usn(defG.left() ) ) * (  var1(u.left() ,defG.left() ) * unv(defG.right()) ).tr() // right left
                    ,
                     alpha_r * (  var1(u.right(),defG.right()) * unv(defG.left() ) ) * (  var1(u.right(),defG.right()) * unv(defG.left() ) ).tr() // right right
                    + // Symmetry
                     alpha_r * ( ovar1(u.right(),defG.right()) * usn(defG.left() ) ) * ( ovar1(u.right(),defG.right()) * usn(defG.left() ) ).tr() // right right
                     );

    // gsDebugVar(u.mapper());
    // gsDebugVar(u.mapper().asVector(0).transpose());
    // gsDebugVar(u.mapper().asVector(1).transpose());
    // gsDebugVar(u.mapper().asVector(2).transpose());

    // gsExprEvaluator<> ev2(A);
    // gsVector<> pt1(2);
    // pt1<<0,0.5;
    // gsVector<> pt2(2);
    // pt2<<1,0.5;

    // gsDebugVar(ev2.eval(defG,pt1,0));
    // gsDebugVar(ev2.eval(var1(u,defG),pt1,0));
    // gsDebugVar(ev2.eval(unv(defG),pt1,0));
    // gsDebugVar(ev2.eval(var1(u,defG) * unv(defG),pt1,0));

    // gsDebugVar(ev2.eval(defG,pt2,1));
    // gsDebugVar(ev2.eval(var1(u,defG),pt2,1));
    // gsDebugVar(ev2.eval(unv(defG),pt2,1));
    // gsDebugVar(ev2.eval(var1(u,defG) * unv(defG),pt2,1));

    if (verbose>1)
    {
        gsDebugVar(A.matrix().rows());
        gsDebugVar(A.matrix().cols());
        gsDebugVar(A.matrix().toDense());
    }
    if (verbose>1) gsDebugVar(A.rhs().transpose());

    A.assembleBdr(
        bc.get("Weak Dirichlet")
        ,
        alpha_d * u * u.tr()
        ,
        alpha_d * (u * (defG - G) - u * (g_N) )
    );

    // RHS is minussed, why?
    A.assembleBdr(
        bc.get("Weak Clamped")
        ,
        alpha_r * (usn(defG).tr() * nv(G) - usn(G).tr() * nv(G)).val() * var2(u,u,defG,nv(G).tr()) * tv(G).norm()
        +
        alpha_r * (var1(u,defG) * nv(G)) * (var1(u,defG) * nv(G)).tr() * meas(G)
        ,
        alpha_r * (usn(defG).tr() * nv(G) - usn(G).tr() * nv(G)).val() * var1(u,defG) * nv(G) * tv(G).norm()
    );

    // For Neumann conditions
    A.assembleBdr(bc.get("Neumann"), u * g_N * tv(G).norm() );

    A.assemble(
        (N_der * (E_m_der).tr() + M_der * (E_f_der).tr() ) * meas(G)
        ,
        u * F  * meas(G) + pressure * u * usn(defG) * meas(G)
        );

    //! [System assembly]

    // Residual<real_t> ResidualFun(mp,dbasis,t,force,materialMat,bc,alpha_d,alpha_r,continuity,verbose);

    // Test<real_t> TestFun(mp,dbasis,t,force,materialMat,bc,alpha_d,alpha_r,continuity,verbose);


    gsMatrix<> result;

    gsMatrix<> zeros(A.numDofs(),1);
    zeros.setZero();

    // ResidualFun.jacobian_into(zeros,result);
    // gsDebugVar(result);
    // gsDebugVar(A.matrix().toDense());

    //! [Linear solve]
    // solve system
    solver.compute( A.matrix() );
    gsMatrix<> solVector = solver.solve(A.rhs());

    // update deformed patch
    gsMatrix<> cc;

    u_sol.setSolutionVector(solVector);
    for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
    {
        // // extract deformed geometry
        u_sol.extract(cc, k);
        mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    // gsVector<> pt(2);
    // pt.setConstant(0.25);
    // ev.options().setInt("plot.npts", 5000);
    // ev.writeParaview(S_m_plot,defG,"stress");
    // if (verbose>0) gsInfo<<"Stresses\n"<<ev.eval(E_m_plot,pt)<<"\n";
    //! [Linear solve]

    //! [Nonlinear solve]
    real_t residual = A.rhs().norm();
    real_t residual0 = residual;
    real_t residualOld = residual;
    gsMatrix<> updateVector = solVector;

    if (nonlinear)
    {
        index_t itMax = 20;
        real_t tol = 1e-8;
        for (index_t it = 0; it != itMax; ++it)
        {
            A.initSystem();
            // assemble system

            if (verbose>1) gsDebugVar(solVector.transpose());

            // RHS is minussed!
            if (continuity > -1)
            A.assembleIfc(dbasis.topology().interfaces(),
                            -alpha_d * u.left() * du
                            ,
                             alpha_d * u.right()* du
                            ,
                             alpha_d * u.left() * u.left().tr()
                            ,
                            -alpha_d * u.right()* u.left() .tr()
                            ,
                            -alpha_d * u.left() * u.right().tr()
                            ,
                             alpha_d * u.right()* u.right().tr()
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Penalty of out-of-plane coupling
            // dW^pr / du_r --> first line
            if (continuity > 0)
            A.assembleIfc(dbasis.topology().interfaces(),
                            -alpha_r * dN_lr * var1(u.left(),defG.left())   * usn(defG.right())
                            ,
                            -alpha_r * dN_lr * var1(u.right(),defG.right()) * usn(defG.left() )
                            ,
                            // Symmetry
                            -alpha_r * dN_rl * var1(u.right(),defG.right())   * usn(defG.left())
                            ,
                            -alpha_r * dN_rl * var1(u.left(),defG.left()) * usn(defG.right() )
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Penalty of in-plane coupling
            // dW^pr / du_r --> second line
            if (continuity > 0)
            A.assembleIfc(dbasis.topology().interfaces(),
                            -alpha_r * dnN_lr* ovar1(u.left(),defG.left())  * usn(defG.right())
                            ,
                            -alpha_r * dnN_lr* var1(u.right(),defG.right()) * unv(defG.left() )
                            ,
                            // Symmetry
                            -alpha_r * dnN_rl* ovar1(u.right(),defG.right())  * usn(defG.left())
                            ,
                            -alpha_r * dnN_rl* var1(u.left(),defG.left()) * unv(defG.right() )
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Penalty of out-of-plane coupling
            // dW^pr / du_r --> first line
            if (continuity > 0)
            A.assembleIfc(dbasis.topology().interfaces(),
                             alpha_r * dN_lr * var2mod(u.left() ,u.left() ,defG.left() ,usn(defG.right()).tr() )      // left left
                            ,
                             alpha_r * dN_lr * ( var1(u.left() ,defG.left() ) * var1(u.right(),defG.right()).tr() )// left right
                            ,
                             alpha_r * dN_lr * ( var1(u.right(),defG.right()) * var1(u.left() ,defG.left() ).tr() )// right left
                            ,
                             alpha_r * dN_lr * var2mod( u.right(),u.right(),defG.right(),usn(defG.right() ).tr() )     // right right
                            ,
                            // Symmetry
                             alpha_r * dN_rl * var2mod(u.right() ,u.right() ,defG.right() ,usn(defG.left()).tr() )      // right right
                            ,
                             alpha_r * dN_rl * ( var1(u.right(),defG.right() ) * var1(u.left(),defG.left()).tr() )// right left
                            ,
                             alpha_r * dN_rl * ( var1(u.left() ,defG.left() ) * var1(u.right() ,defG.right() ).tr() )// left right
                            ,
                             alpha_r * dN_rl * var2mod( u.left() ,u.left() ,defG.left() ,usn(defG.right() ).tr() )     // left left
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Penalty of out-of-plane coupling
            // dW^pr / du_r --> second line
            if (continuity > 0)
            A.assembleIfc(dbasis.topology().interfaces(),
                             alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()    // left left
                            ,
                             alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()   // left right
                            ,
                             alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()   // right left
                            ,
                             alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()  // right right
                            ,
                            // Symmetry
                             alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()    // right right
                            ,
                             alpha_r * ( var1(u.right(),defG.right()) * usn(defG.left()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()   // right left
                            ,
                             alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.right(),defG.right()) * usn(defG.left()) ).tr()   // left right
                            ,
                             alpha_r * ( var1(u.left(),defG.left()) * usn(defG.right()) ) * ( var1(u.left(),defG.left()) * usn(defG.right()) ).tr()  // left left
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Penalty of in-plane coupling
            // dW^pr / du_r --> third line
            if (continuity > 0)
            A.assembleIfc(dbasis.topology().interfaces(),
                             alpha_r * dnN_lr * ovar2(u.left(),u.left(),defG.left(),usn(defG.right()).tr())            // left left
                            ,
                             alpha_r * dnN_lr * ( ovar1(u.left() ,defG.left() ) * var1(u.right(),defG.right()).tr() ) // left right
                            ,
                             alpha_r * dnN_lr * ( ovar1(u.right(),defG.right()) * var1(u.left() ,defG.left() ).tr() ) // left right
                            ,
                             alpha_r * dnN_lr * ovar2(u.right(),u.right(),defG.right(),usn(defG.left()).tr())          // right right
                            ,
                            // Symmetry
                             alpha_r * dnN_rl * ovar2(u.right(),u.right(),defG.right(),usn(defG.left()).tr())            // right right -- dterm3
                            ,
                             alpha_r * dnN_rl * ( ovar1(u.right(),defG.right()) * var1(u.left() ,defG.left() ).tr() ) // right left
                            ,
                             alpha_r * dnN_rl * ( ovar1(u.left() ,defG.left() ) * var1(u.right(),defG.right()).tr() ) // right left -- dterm3       --->>> gives -0.00163163 in top-right
                            ,
                             alpha_r * dnN_rl * ovar2(u.left(),u.left(),defG.left(),usn(defG.right()).tr())          // left left
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Penalty of in-plane coupling
            // dW^pr / du_r --> fourth line
            if (continuity > 0)
            A.assembleIfc(dbasis.topology().interfaces(),
                             alpha_r * ( ovar1(u.left() ,defG.left()) * usn(defG.right()) ) * ( ovar1(u.left(),defG.left() ) * usn(defG.right()) ).tr()    // left left -- dterm1
                            ,
                             alpha_r * ( ovar1(u.left() ,defG.left()) * usn(defG.right()) ) * ( var1(u.right(),defG.right()) * unv(defG.left() ) ).tr()   // left right
                            ,
                             alpha_r * ( var1(u.right(),defG.right()) * unv(defG.left())  ) * ( ovar1(u.left(),defG.left() ) * usn(defG.right()) ).tr()   // right left
                            ,
                             alpha_r * ( var1(u.right(),defG.right()) * unv(defG.left())  ) * ( var1(u.right(),defG.right()) * unv(defG.left() ) ).tr()  // right right -- dterm1
                            ,
                            // Symmetry
                             alpha_r * ( ovar1(u.right() ,defG.right()) * usn(defG.left()) ) * ( ovar1(u.right(),defG.right() ) * usn(defG.left()) ).tr()    // right right -- dterm3
                            ,
                             alpha_r * ( ovar1(u.right() ,defG.right()) * usn(defG.left()) ) * ( var1(u.left(),defG.left()) * unv(defG.right() ) ).tr()   // right left
                            ,
                             alpha_r * ( var1(u.left(),defG.left()) * unv(defG.right())  ) * ( ovar1(u.right(),defG.right() ) * usn(defG.left()) ).tr()   // left right -- dterm3
                            ,
                             alpha_r * ( var1(u.left(),defG.left()) * unv(defG.right())  ) * ( var1(u.left(),defG.left()) * unv(defG.right() ) ).tr()  // left left
                             );

            if (verbose>1) gsDebugVar(A.matrix().toDense());
            if (verbose>1) gsDebugVar(A.rhs().transpose());

            // Weak Dirichlet term
            // RHS is minussed!
            A.assembleBdr(
                bc.get("Weak Dirichlet")
                ,
                alpha_d * u * u.tr()
                ,
                -alpha_d * (u * (defG - G) - u * (g_N) )
            );

            // Weak Clamping term
            // RHS is minussed!
            A.assembleBdr(
                bc.get("Weak Clamped")
                ,
                alpha_r * (usn(defG).tr() * nv(G) - usn(G).tr() * nv(G)).val() * var2(u,u,defG,nv(G).tr()) * tv(G).norm()
                +
                alpha_r * (var1(u,defG) * nv(G)) * (var1(u,defG) * nv(G)).tr() * meas(G)
                ,
                -alpha_r * (usn(defG).tr() * nv(G) - usn(G).tr() * nv(G)).val() * var1(u,defG) * nv(G) * tv(G).norm()
            );

            // Neumann term
            A.assembleBdr(bc.get("Neumann"), u * g_N * tv(G).norm() );

            A.assemble(
                (
                  N_der * E_m_der.tr()
                    +
                  E_m_der2
                    +
                  M_der * E_f_der.tr()
                    +
                  E_f_der2
                  ) * meas(G)
                    -
                  pressure * u * var1(u,defG) .tr() * meas(G)
                , u * F * meas(G) + pressure * u * usn(defG) * meas(G) - ( ( N * E_m_der.tr() + M * E_f_der.tr() ) * meas(G) ).tr()
                );



            // NormalFunction w(...);
            // w.deriv_into
            // w.gsFunction<real_t>::deriv_into

            // ResidualFun.deriv_into(solVector,result);
            // result.resize(ResidualFun.targetDim(),ResidualFun.targetDim());
            // gsDebugVar(result);

            // return 0;
            // solve system
            gsSparseMatrix<> mat = result.sparseView();
            solver.compute( A.matrix() );
            // solver.compute( mat );
            // gsDebugVar(A.matrix().toDense());
            updateVector = solver.solve(A.rhs()); // this is the UPDATE

            // gsDebugVar((A.matrix().toDense() - result).lpNorm<gsEigen::Infinity>());


            // gsDebugVar(A.rhs().transpose());
            // ResidualFun.eval_into(solVector,result);
            // gsDebugVar(result.transpose());



            // return 0;

            solVector += updateVector;
            residual = A.rhs().norm();

            if (verbose>0) gsInfo <<"Iteration: "<< it
                                  <<", residue: "<< residual
                                  <<", update norm: "<<updateVector.norm()
                                  <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
                                  <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
                                  <<"\n";

            residualOld = residual;

            // update deformed patch
            u_sol.setSolutionVector(updateVector);
            for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
            {
                // // extract deformed geometry
                u_sol.extract(cc, k);
                mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
            }

            if (residual < tol)
                break;
        }
    }
    //! [Nonlinear solve]

    //! [Construct solution]
    u_sol.setSolutionVector(solVector);
    mp_def = mp;
    for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
    {
        u_sol.extract(cc, k);
        mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    if (verbose>0) gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    if (verbose>0) gsInfo <<"Minimum deformation coef: "
                          << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    if (verbose>0) gsInfo <<"Area (undeformed) = "<<ev.integral(meas(G))<<"\tArea (undeformed) = "<<ev.integral(meas(defG))<<"\n";
    //! [Construct solution]

    //! [Evaluate solution]
    gsMatrix<> coords(2,1);
    if (testCase==1)
      coords<<0,0;
    else if (testCase==2)
      coords<<0.5,0;
    else if (testCase==3)
      coords<<0,0;
    else if (testCase==4)
      coords<<1,1;
    else
      coords<<0,0;

    gsMatrix<> res(3,1);
    mp.patch(0).eval_into(coords,res);
    real_t x=res.at(0);
    real_t y=res.at(1);
    real_t z=res.at(2);

    deformation.patch(0).eval_into(coords,res);
    real_t ux=res.at(0);
    real_t uy=res.at(1);
    real_t uz=res.at(2);

    if (verbose>0) gsInfo<<"Deformation on point ("<<x<<","<<y<<","<<z<<"):\n";
    if (verbose>0) gsInfo<<std::setprecision(20)<<"("<<ux<<","<<uy<<","<<uz<<")"<<"\n";
    //! [Evaluate solution]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp, deformation);
        if (verbose>0) gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // gsFileManager::open("solution.pvd");
    }
    //! [Export visualization in ParaView]
    return EXIT_SUCCESS;

}// end main

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B)
{
  int q = p;
  int m = n;
  gsMultiPatch<T> mp = RectangularDomain(n, m, p, q, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> Panel(T Lp, T Wp, T Hw, T Wf, T x, T y, T z)

{
    gsMultiPatch<T> result, tmp;

    // Base plate, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< Wp/2,0,0;
    result.patch(0).coefs().row(2)<< 0,Lp,0;
    result.patch(0).coefs().row(3)<< Wp/2,Lp,0;

    // Base plate, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< -Wp/2,0,0;
    result.patch(1).coefs().row(1)<< 0,0,0;
    result.patch(1).coefs().row(2)<< -Wp/2,Lp,0;
    result.patch(1).coefs().row(3)<< 0,Lp,0;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(2).embed(3);
    result.patch(2).coefs().row(0)<< 0,0,0;
    result.patch(2).coefs().row(1)<< 0,Lp,0;
    result.patch(2).coefs().row(2)<< 0,0,Hw;
    result.patch(2).coefs().row(3)<< 0,Lp,Hw;

    // Flange, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(3).embed(3);
    result.patch(3).coefs().row(0)<< 0,0,Hw;
    result.patch(3).coefs().row(1)<< Wf/2,0,Hw;
    result.patch(3).coefs().row(2)<< 0,Lp,Hw;
    result.patch(3).coefs().row(3)<< Wf/2,Lp,Hw;

    // Flange, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(4).embed(3);
    result.patch(4).coefs().row(0)<< -Wf/2,0,Hw;
    result.patch(4).coefs().row(1)<< 0,0,Hw;
    result.patch(4).coefs().row(2)<< -Wf/2,Lp,Hw;
    result.patch(4).coefs().row(3)<< 0,Lp,Hw;

    for (index_t p = 0; p!=5; p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.addInterface(&result.patch(1),2,&result.patch(0),1);
    result.addInterface(&result.patch(2),3,&result.patch(0),1);
    result.addInterface(&result.patch(2),4,&result.patch(3),1);
    result.addInterface(&result.patch(2),4,&result.patch(4),2);

    result.addAutoBoundaries();

    return result;
}
