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
    // Set default options
    this->_defaultOptions();

    m_nLayers = nLayers;
    m_Ts.resize(m_nLayers);
    m_Gs.resize(m_nLayers);
    m_As.resize(m_nLayers);
    m_Rs.resize(m_nLayers);
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u, bool basis) const
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

    std::vector<gsMatrix<T>> Tcontainer(m_nLayers);
    for (size_t k=0; k!= m_Rs.size(); k++)
    {
        Tcontainer[k] = m_Ts[k]->eval(map.values[0]);
        GISMO_ASSERT(Tcontainer[k].rows()==1,"Thickness has the wrong size, must be scalar");
    }

    for (index_t i = 0; i != m_nLayers; ++i) // layers
        for (index_t k = 0; k != u.cols(); ++k) // points
            result.row(i) = Tcontainer[i].row(0);
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u);
    // Initialize and result
    gsMatrix<T> result(9, u.cols());
    gsMatrix<T,3,3> Dmat, Cmat, Tmat;

    T t, t_tot, z_mid, t_temp, dz;

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
        this->_getMetric(k,0.0); // on point i, with height 0.0

        // Compute total thickness (sum of entries)
        t_tot = 0;
        for (size_t i=0; i!=Tcontainer.size(); i++)
            t_tot += Tcontainer[i](0,k);

        // compute mid-plane height of total plate
        z_mid = t_tot / 2.0;

        // now we use t_temp to add the thickness of all plies iteratively
        t_temp = 0.0;
        t = 0.0;

        for (index_t i = 0; i != m_nLayers; ++i) // loop over laminates
        {
            // Lookup all quantities
            t = Tcontainer[i](0,k);

            // Transform the matrix: T^T D T
            Dmat = Acontainer[i].reshapeCol(k,3,3).transpose() * Gcontainer[i].reshapeCol(k,3,3) * Acontainer[i].reshapeCol(k,3,3);

            // distance from mid of the ply to mid of the plate
            // dz = math::abs(z_mid - (t/2.0 + t_temp) ); // distance from mid-plane of plate
            dz = z_mid - (t/2.0 + t_temp); // distance from mid-plane of plate

            if      (out==MaterialOutput::MatrixA)
                    Cmat += Dmat * t; // A
            else if (out==MaterialOutput::MatrixB || out==MaterialOutput::MatrixC)
                    Cmat += Dmat * t*dz; // B
            else if (out==MaterialOutput::MatrixD)
                    Cmat += Dmat * ( t*math::pow(dz,2) + math::pow(t,3)/12.0 ); // D
            else
                GISMO_ERROR("MaterialOutput unknown");
            t_temp += t;
        }

        GISMO_ASSERT(std::abs(t_tot-t_temp)/t_tot < 1e-12,"Total thickness after loop is wrong. t_temp = "<<t_temp<<", sum(thickness) = "<<t_tot<<", difference = "<<t_tot-t_temp);

        gsVector<T,3> a1,ac1,e1,a2,ac2,e2;
        a1 = m_data.mine().m_gcov_ori.col(0);               a2 = m_data.mine().m_gcov_ori.col(1);
        ac1 = m_data.mine().m_gcon_ori.col(0);              ac2 = m_data.mine().m_gcon_ori.col(1);
        e1 = m_data.mine().m_gcov_ori.col(0).normalized();  e2 = m_data.mine().m_gcon_ori.col(1).normalized();

        Cmat = _cart2cov(a1,a2,e1,e2) * Cmat;
        Cmat = Cmat * _con2cart(ac1,ac2,e1,e2);

        result.reshapeCol(k,3,3) = Cmat;

    }
    return result;
}

template <short_t dim, class T >
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const
{
    // gsInfo<<"TO DO: evaluate moments using thickness";
    // Input: u in-plane points
    //        z matrix with, per point, a column with z integration points
    // Output: (n=u.cols(), m=z.cols())
    //          [(u1,z1) (u2,z1) ..  (un,z1), (u1,z2) ..  (un,z2), ..,  (u1,zm) .. (un,zm)]

    this->_computePoints(patch,u);
    gsMatrix<T> result(3, u.cols());
    enum MaterialOutput _out;
    if      (out==MaterialOutput::VectorN)
        _out=MaterialOutput::MatrixA;
    else if (out==MaterialOutput::VectorM)
        _out=MaterialOutput::MatrixD;
    else
        GISMO_ERROR("Output unknown");

    gsMatrix<T> Cmat = this->eval3D_matrix(patch,u,z,_out);
    gsMatrix<T> Eij(2,2);

    for (index_t k = 0; k != u.cols(); ++k)
    {
        this->_getMetric(k,0.0); // on point i, with height 0.0

        if      (out == MaterialOutput::VectorN) // To be used with multiplyZ_into
            Eij = 0.5*(m_data.mine().m_Acov_def - m_data.mine().m_Acov_ori);
        else if (out == MaterialOutput::VectorM) // To be used with multiplyZ_into
            Eij = (m_data.mine().m_Bcov_ori - m_data.mine().m_Bcov_def);
        else if (out == MaterialOutput::Generic) // To be used with multiplyLinZ_into or integrateZ_into
            Eij = 0.5*(m_data.mine().m_Acov_def - m_data.mine().m_Acov_ori) + z*(m_data.mine().m_Bcov_ori - m_data.mine().m_Bcov_def);
        else
            GISMO_ERROR("Output type is not VectorN, VectorM or Generic!");

        Eij(0,1) += Eij(1,0);
        std::swap(Eij(1,0), Eij(1,1));
        Eij.resize(4,1);
        Eij.conservativeResize(3,1);
        result.col(k) = Cmat.reshapeCol(k,3,3) * Eij;
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

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::_cart2cov(const gsVector<T> & a1, const gsVector<T> & a2, const gsVector<T> & e1, const gsVector<T> & e2) const
{
    gsMatrix<T,3,3> Tmat;
    Tmat.setZero();
    // Covariant to local cartesian
    Tmat(0,0) = (e1.dot(a1))*(a1.dot(e1));
    Tmat(0,1) = (e1.dot(a2))*(a2.dot(e2));
    Tmat(0,2) = 2*(e1.dot(a1))*(a2.dot(e1));
    // Row 1
    Tmat(1,0) = (e2.dot(a1))*(a1.dot(e2));
    Tmat(1,1) = (e2.dot(a2))*(a2.dot(e2));
    Tmat(1,2) = 2*(e2.dot(a1))*(a2.dot(e2));
    // Row 2
    Tmat(2,0) = (e1.dot(a1))*(a1.dot(e2));
    Tmat(2,1) = (e1.dot(a2))*(a2.dot(e2));
    Tmat(2,2) = (e1.dot(a1))*(a2.dot(e2)) + (e1.dot(a2))*(a1.dot(e2));

    Tmat = Tmat.template block<3,3>(0,0).inverse(); // !!!!
    return Tmat;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::_con2cart(const gsVector<T> & ac1, const gsVector<T> & ac2, const gsVector<T> & e1, const gsVector<T> & e2) const
{
    gsMatrix<T,3,3> Tmat;
    Tmat.setZero();
    // Contravariant to local cartesian
    Tmat(0,0) = (e1.dot(ac1))*(ac1.dot(e1));
    Tmat(0,1) = (e1.dot(ac2))*(ac2.dot(e2));
    Tmat(0,2) = 2*(e1.dot(ac1))*(ac2.dot(e1));
    // Row 1
    Tmat(1,0) = (e2.dot(ac1))*(ac1.dot(e2));
    Tmat(1,1) = (e2.dot(ac2))*(ac2.dot(e2));
    Tmat(1,2) = 2*(e2.dot(ac1))*(ac2.dot(e2));
    // Row 2
    Tmat(2,0) = (e1.dot(ac1))*(ac1.dot(e2));
    Tmat(2,1) = (e1.dot(ac2))*(ac2.dot(e2));
    Tmat(2,2) = (e1.dot(ac1))*(ac2.dot(e2)) + (e1.dot(ac2))*(ac1.dot(e2));
    return Tmat;
}


} // end namespace
