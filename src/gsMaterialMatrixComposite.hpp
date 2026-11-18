/** @file gsMaterialMatrixComposite.hpp

    @brief Provides a material matrix for laminates

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

/*
    To Do [updated 16-06-2020]:(rh)
    - Make beta (compressible materials) and material parameters universal for all integration points over the thickness. So get them out of the dPsi functions etc and move them into the integration loops as global variables.

*/



#pragma once

#include <gsKLShell/src/gsMaterialMatrixComposite.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsCore/gsFunction.h>

namespace gismo
{

template <short_t dim, class T >
gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
                        const gsFunctionSet<T>                  & mp,
                        const std::vector< gsFunctionSet<T> *>  & thickness,
                        const std::vector< gsFunctionSet<T> *>  & G,
                        const std::vector< gsFunctionSet<T> *>  & alpha,
                        const std::vector< gsFunctionSet<T> *>  & rho           )
                        :
                        Base(&mp)
{
    index_t nLayers = thickness.size();
    GISMO_ASSERT(nLayers==(index_t)G.size(),    "Size error in layer input");
    GISMO_ASSERT(nLayers==(index_t)alpha.size(),"Size error in layer input");
    GISMO_ASSERT(nLayers==(index_t)rho.size(),  "Size error in layer input");

    for (size_t k=0; k!=thickness.size(); k++)
    {
        GISMO_ENSURE(thickness[k]!=NULL,"Function "<<k<<" of the thickness vector is not defined!");
        GISMO_ENSURE(G[k]!=NULL,"Function "<<k<<" of the shear moduli vector is not defined!");
        GISMO_ENSURE(alpha[k]!=NULL,"Function "<<k<<" of the alpha vector is not defined!");
    }
    for (size_t k=0; k!=rho.size(); k++)
        GISMO_ENSURE(rho[k]!=NULL,"Function "<<k<<" of the density vector is not defined!");

    _initialize(nLayers);

    for (index_t k=0; k!=nLayers; k++)
    {
        m_Ts[k] = memory::make_shared_not_owned(thickness[k]);
        m_Gs[k] = memory::make_shared_not_owned(G[k]);
        m_As[k] = memory::make_shared_not_owned(alpha[k]);
        m_Rs[k] = memory::make_shared_not_owned(rho[k]);
    }
}

template <short_t dim, class T >
gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
                        const gsFunctionSet<T>                  & mp,
                        const std::vector< gsFunctionSet<T> *>  & thickness,
                        const std::vector< gsFunctionSet<T> *>  & G,
                        const std::vector< gsFunctionSet<T> *>  & alpha         )
                        :
                        Base(&mp)
{
    index_t nLayers = thickness.size();
    GISMO_ASSERT(nLayers==(index_t)G.size(),    "Size error in layer input");
    GISMO_ASSERT(nLayers==(index_t)alpha.size(),"Size error in layer input");
    
    for (size_t k=0; k!=thickness.size(); k++)
    {
        GISMO_ENSURE(thickness[k]!=NULL,"Function "<<k<<" of the thickness vector is not defined!");
        GISMO_ENSURE(G[k]!=NULL,"Function "<<k<<" of the shear moduli vector is not defined!");
        GISMO_ENSURE(alpha[k]!=NULL,"Function "<<k<<" of the alpha vector is not defined!");
    }

    _initialize(nLayers);
    for (index_t k=0; k!=nLayers; k++)
    {
        m_Ts[k] = memory::make_shared_not_owned(thickness[k]);
        m_Gs[k] = memory::make_shared_not_owned(G[k]);
        m_As[k] = memory::make_shared_not_owned(alpha[k]);
    }
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_defaultOptions()
{
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T >
std::ostream & gsMaterialMatrixComposite<dim,T>::print(std::ostream &os) const
{
    os  <<"---------------------------------------------------------------------\n"
        <<"----------------------Composite Material Info------------------------\n"
        <<"---------------------------------------------------------------------\n\n";

    os  <<"Material model: \t";
    os  <<"Saint-Venant Kirchhoff";
    os  <<"\n";

    os  <<"---------------------------------------------------------------------\n\n";
    return os;
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_initialize(const index_t nLayers)
{
    GISMO_ENSURE(nLayers > 0, "gsMaterialMatrixComposite: a laminate needs at least one ply");

    // Set default options
    this->_defaultOptions();

    m_nLayers = nLayers;
    m_Ts.resize(m_nLayers);
    m_Gs.resize(m_nLayers);
    m_As.resize(m_nLayers);
    m_Rs.resize(m_nLayers);
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u, bool /*basis*/) const
{
    this->_computeMetricUndeformed(patch,u);

    if (Base::m_defpatches->nPieces()!=0)
        this->_computeMetricDeformed(patch,u);
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT(m_Ts.size()==m_Rs.size(),"Size of vectors of thickness and densities is not equal: " << m_Ts.size()<<" & "<<m_Rs.size());

    gsMapData<T> map;
    map.points = u;
    map.flags = NEED_VALUE;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    result.resize(1, u.cols());
    result.setZero();
    GISMO_ENSURE(m_Rs.front() != nullptr, "gsMaterialMatrixComposite::density_into: no density was given; construct with the rho argument to use density_into");
    std::vector<gsMatrix<T>> Tcontainer(m_nLayers);
    std::vector<gsMatrix<T>> Rcontainer(m_nLayers);
    for (size_t k=0; k!= m_Rs.size(); k++)
    {
        Tcontainer[k] = m_Ts[k]->piece(patch).eval(map.values[0]);
        Rcontainer[k] = m_Rs[k]->piece(patch).eval(map.values[0]);

        GISMO_ASSERT(Tcontainer[k].rows()==1,"Thickness has the wrong size, must be scalar");
        GISMO_ASSERT(Rcontainer[k].rows()==1,"Densities has the wrong size, must be scalar");
    }


    for (index_t i = 0; i != m_nLayers; ++i) // layers
        for (index_t k = 0; k != u.cols(); ++k) // points
            result(0,k) += Tcontainer[i](0,k)*Rcontainer[i](0,k);

}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    gsMapData<T> map;
    map.points = u;
    map.flags = NEED_VALUE;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    result.resize(1, u.cols());
    result.setZero();
    gsMatrix<T> Tply;
    for (size_t k = 0; k != m_Ts.size(); ++k)
    {
        m_Ts[k]->piece(patch).eval_into(map.values[0], Tply);
        GISMO_ASSERT(Tply.rows()==1,"Thickness has the wrong size, must be scalar");
        result += Tply;
    }
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & /*z*/, enum MaterialOutput out) const
{
    // z is ignored: the laminate is integrated through the thickness analytically, ply by ply.
    // Output: 9 x u.cols(); column k is the 3x3 Voigt matrix (column-major, reshapeCol(k,3,3))
    // in contravariant components at point k, with
    //   MatrixA:              A = sum_i t_i C_i
    //   MatrixB or MatrixC:   B = sum_i t_i z_i C_i
    //   MatrixD:              D = sum_i (t_i z_i^2 + t_i^3/12) C_i
    // Any other output throws.
    GISMO_ENSURE(out==MaterialOutput::MatrixA || out==MaterialOutput::MatrixB ||
                 out==MaterialOutput::MatrixC || out==MaterialOutput::MatrixD,
                 "gsMaterialMatrixComposite::eval3D_matrix: MaterialOutput must be MatrixA, MatrixB, MatrixC or MatrixD");
    gsMatrix<T> result;
    this->_eval3D_laminate(patch,u,out,result,nullptr);
    return result;
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_eval3D_laminate(const index_t patch, const gsMatrix<T> & u, enum MaterialOutput out, gsMatrix<T> & primary, gsMatrix<T> * second) const
{
    // One pass over points and plies for both matrices; the optional second one is always B.
    const bool isA  = (out==MaterialOutput::MatrixA);
    const bool isBC = (out==MaterialOutput::MatrixB || out==MaterialOutput::MatrixC);
    const bool withSecond = (second != nullptr);

    this->_computePoints(patch,u);
    primary.resize(9, u.cols());
    if (withSecond)
        second->resize(9, u.cols());
    gsMatrix<T,3,3> Dmat, Cmat, Bmat;

    T t, t_tot, t_temp, z;

    std::vector<gsMatrix<T>> Tcontainer(m_nLayers);
    std::vector<gsMatrix<T>> Rcontainer(m_nLayers);
    std::vector<gsMatrix<T>> Acontainer(m_nLayers);
    std::vector<gsMatrix<T>> Gcontainer(m_nLayers);

    gsMapData<T> map;
    map.points = u;
    map.flags = NEED_VALUE;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    gsMatrix<T> angles;
    for (size_t k=0; k!= m_Gs.size(); k++)
    {
        Gcontainer[k] = m_Gs[k]->piece(patch).eval(map.values[0]);
        Tcontainer[k] = m_Ts[k]->piece(patch).eval(map.values[0]);
        angles = m_As[k]->piece(patch).eval(map.values[0]);

        GISMO_ASSERT(Gcontainer[k].rows()==9,"G has the wrong size, must be 3x3");
        GISMO_ASSERT(Tcontainer[k].rows()==1,"Thickness has the wrong size, must be scalar");

        Acontainer[k] = _transformationMatrix(angles,map.values[0]);
    }


    for (index_t k = 0; k != u.cols(); ++k)
    {
        Cmat.setZero();
        Bmat.setZero();
        this->_getMetric(k,0.0); // on point i, with height 0.0

        // Compute total thickness (sum of entries)
        t_tot = 0;
        for (size_t i=0; i!=Tcontainer.size(); i++)
            t_tot += Tcontainer[i](0,k);

        // now we use t_temp to add the thickness of all plies iteratively
        t_temp = 0.0;
        t = 0.0;

        for (index_t i = 0; i != m_nLayers; ++i) // loop over laminates
        {
            // Lookup all quantities
            t = Tcontainer[i](0,k);

            // Transform the matrix: T^T D T
            Dmat = Acontainer[i].reshapeCol(k,3,3).transpose() * Gcontainer[i].reshapeCol(k,3,3) * Acontainer[i].reshapeCol(k,3,3);

            // Signed height of the ply mid-surface above the laminate mid-plane, along the
            // normal n; ply 0 is at the most negative z. B = +sum_i t_i z_i C_i = int z C dz.
            z = (t_temp + t/2.0) - t_tot/2.0;

            if      (isA)
                    Cmat += Dmat * t; // A
            else if (isBC)
                    Cmat += Dmat * t*z; // B
            else
                    Cmat += Dmat * ( t*math::pow(z,2) + math::pow(t,3)/12.0 ); // D
            if (withSecond)
                Bmat += Dmat * t*z; // B
            t_temp += t;
        }

        GISMO_ASSERT(std::abs(t_tot-t_temp)/t_tot < 1e-12,"Total thickness after loop is wrong. t_temp = "<<t_temp<<", sum(thickness) = "<<t_tot<<", difference = "<<t_tot-t_temp);

        // Basis change from the local Cartesian frame {e1,e2} to the contravariant basis.
        // The contravariant stiffness tensor is
        //   C^{abcd} = C_{ijkl} (e_i.a^a)(e_j.a^b)(e_k.a^c)(e_l.a^d),
        // with C_{ijkl} the ply stiffness in {e1,e2}, e1 = a1/|a1|, e2 = a^2/|a^2|.
        // With the Voigt map eps_cart = Q eps_cov (engineering shear), invariance of the strain
        // energy gives the contravariant Voigt matrix C^con = Q^T C_cart Q. Q carries the factor 2 on its shear row, so
        // _transformation(E,gcon) = Q^T carries it on its shear column.
        // The congruence of a symmetric matrix is symmetric and needs no matrix inversion.
        const gsMatrix<T> & gcov = m_data.mine().m_gcov_ori;
        const gsMatrix<T> & gcon = m_data.mine().m_gcon_ori;
        gsMatrix<T> E(gcov.rows(),2);
        E.col(0) = gcov.col(0).normalized();
        E.col(1) = gcon.col(1).normalized();
        const gsMatrix<T> Tc = this->_transformation(E,gcon);
        Cmat = Tc * Cmat * Tc.transpose();
        primary.reshapeCol(k,3,3) = Cmat;

        if (withSecond)
        {
            Bmat = Tc * Bmat * Tc.transpose();
            second->reshapeCol(k,3,3) = Bmat;
        }
    }
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & /*z*/, enum MaterialOutput out) const
{
    // Resultants of the membrane strain eps = 0.5*(a_def - a_ori) and the curvature change
    // kappa = b_ori - b_def:
    //   N = A eps + B kappa,   M = B^T eps + D kappa.
    // Input: u in-plane points; z is ignored (the laminate is integrated through the thickness analytically).
    // Output: one column per point.

    if (out!=MaterialOutput::VectorN && out!=MaterialOutput::VectorM)
        GISMO_ERROR("Output type is not VectorN or VectorM!");

    gsMatrix<T> result(3, u.cols());

    // N needs A and B, M needs B and D. The laminate pass computes the points, which sets
    // the metric state used by _getMetric below.
    const bool isN = (out == MaterialOutput::VectorN);
    gsMatrix<T> Cmat, Bmat;
    this->_eval3D_laminate(patch,u,isN ? MaterialOutput::MatrixA : MaterialOutput::MatrixD,Cmat,&Bmat);
    gsMatrix<T,2,2> Eij, Kij;
    gsMatrix<T,3,1> eps, kap;
    gsMatrix<T,3,3> Ck, Bk;

    for (index_t k = 0; k != u.cols(); ++k)
    {
        this->_getMetric(k,0.0); // on point k, with height 0.0

        Eij = 0.5*(m_data.mine().m_Acov_def - m_data.mine().m_Acov_ori);
        Kij = (m_data.mine().m_Bcov_ori - m_data.mine().m_Bcov_def);

        // Voigt fold with engineering shear: [E11, E22, 2E12]
        eps << Eij(0,0), Eij(1,1), Eij(0,1) + Eij(1,0);
        kap << Kij(0,0), Kij(1,1), Kij(0,1) + Kij(1,0);

        Ck = Cmat.reshapeCol(k,3,3);
        Bk = Bmat.reshapeCol(k,3,3);
        if (isN)
            result.col(k) = Ck * eps + Bk * kap;
        else
            result.col(k) = Bk.transpose() * eps + Ck * kap;
    }
    return result;
}

//-----------------------------------------------------------------------------------------------------------------------
template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::_transformationMatrix(const gsMatrix<T> & phi, const gsMatrix<T> & u) const
{
    gsMatrix<T> result(9,u.cols());
    GISMO_ASSERT(phi.rows()==1,"Angles has the wrong size, must be scalar");

    for (index_t k=0; k!=u.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> Tmat = result.reshapeCol(k,3,3);
        // Make transformation matrix
        Tmat(0,0) = Tmat(1,1) = math::pow(math::cos(phi(0,k)),2);
        Tmat(0,1) = Tmat(1,0) = math::pow(math::sin(phi(0,k)),2);
        Tmat(2,0) = Tmat(0,2) = Tmat(2,1) = Tmat(1,2) = math::sin(phi(0,k)) * math::cos(phi(0,k));
        Tmat(2,0) *= -2.0;
        Tmat(2,1) *= 2.0;
        Tmat(1,2) *= -1.0;
        Tmat(2,2) = math::pow(math::cos(phi(0,k)),2) - math::pow(math::sin(phi(0,k)),2);
    }

    // Compute laminate stiffness matrix
    return result;
}


} // end namespace
