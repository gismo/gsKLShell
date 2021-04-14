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
    To Do [updated 16-06-2020]:
    - Make beta (compressible materials) and material parameters universal for all integration points over the thickness. So get them out of the dPsi functions etc and move them into the integration loops as global variables.

*/



#pragma once

#include <gsKLShell/gsMaterialMatrixComposite.h>
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsCore/gsFunction.h>

namespace gismo
{

template <short_t dim, class T>
gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const gsFunction<T>                                 & thickness,
                            const gsFunction<T>                                 & E11,
                            const gsFunction<T>                                 & E22,
                            const gsFunction<T>                                 & G12,
                            const gsFunction<T>                                 & nu12,
                            const gsFunction<T>                                 & nu21,
                            const gsFunction<T>                                 & phi           )
                            :
                            m_patches(&mp),
                            m_thickness(&thickness),
                            m_E11(&E11),
                            m_E22(&E22),
                            m_G12(&G12),
                            m_nu12(&nu12),
                            m_nu21(&nu21),
                            m_phi(&phi)
{
    _initialize();
}

template <short_t dim, class T>
gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const gsFunction<T>                                 & thickness,
                            const gsFunction<T>                                 & E11,
                            const gsFunction<T>                                 & E22,
                            const gsFunction<T>                                 & G12,
                            const gsFunction<T>                                 & nu12,
                            const gsFunction<T>                                 & nu21,
                            const gsFunction<T>                                 & phi,
                            const gsFunction<T>                                 & rho           )
                            :
                            m_patches(&mp),
                            m_thickness(&thickness),
                            m_E11(&E11),
                            m_E22(&E22),
                            m_G12(&G12),
                            m_nu12(&nu12),
                            m_nu21(&nu21),
                            m_phi(&phi),
                            m_rho(&rho)
{
    _initialize();
}

template <short_t dim, class T>
gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const gsFunction<T>                                 & thickness,
                            const std::vector<gsFunction<T>*>                   & pars          )
                            :
                            m_patches(&mp),
                            m_thickness(&thickness),
                            m_pars(pars),
                            m_rho(nullptr)
{
    _initializeParameters();
    _initialize();
}

template <short_t dim, class T>
gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const gsFunction<T>                                 & thickness,
                            const std::vector<gsFunction<T>*>                   & pars,
                            const gsFunction<T>                                 & rho           )
                            :
                            m_patches(&mp),
                            m_thickness(&thickness),
                            m_pars(pars),
                            m_rho(&rho)
{
    _initializeParameters();
    _initialize();
}

// // ---------------------------------------------------------------------------------------------------
// template <short_t dim, class T>
// gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
//                             const gsFunctionSet<T>                              & mp,
//                             const gsVector<T>                                   & E11,
//                             const gsVector<T>                                   & E22,
//                             const gsVector<T>                                   & G12,
//                             const gsVector<T>                                   & nu12,
//                             const gsVector<T>                                   & nu21,
//                             const gsVector<T>                                   & thickness,
//                             const gsVector<T>                                   & phi           )
//                             :
//                             m_patches(&mp)
// {
//     m_E11       = new gsConstantFunction<T>(E11,3);
//     m_E22       = new gsConstantFunction<T>(E22,3);
//     m_G12       = new gsConstantFunction<T>(G12,3);
//     m_nu12      = new gsConstantFunction<T>(nu12,3);
//     m_nu21      = new gsConstantFunction<T>(nu21,3);
//     m_thickness = new gsConstantFunction<T>(thickness,3);
//     m_phi       = new gsConstantFunction<T>(phi,3);
//     _initialize();
// }

// template <short_t dim, class T>
// gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
//                             const gsFunctionSet<T>                              & mp,
//                             const gsVector<T>                                   & E11,
//                             const gsVector<T>                                   & E22,
//                             const gsVector<T>                                   & G12,
//                             const gsVector<T>                                   & nu12,
//                             const gsVector<T>                                   & nu21,
//                             const gsVector<T>                                   & thickness,
//                             const gsVector<T>                                   & phi,
//                             const gsVector<T>                                   & rho           )
//                             :
//                             m_patches(&mp)
// {
//     m_E11       = new gsConstantFunction<T>(E11,3);
//     m_E22       = new gsConstantFunction<T>(E22,3);
//     m_G12       = new gsConstantFunction<T>(G12,3);
//     m_nu12      = new gsConstantFunction<T>(nu12,3);
//     m_nu21      = new gsConstantFunction<T>(nu21,3);
//     m_thickness = new gsConstantFunction<T>(thickness,3);
//     m_phi       = new gsConstantFunction<T>(phi,3);
//     m_rho       = new gsConstantFunction<T>(rho,3);
//     _initialize();
// }
// template <short_t dim, class T>
// gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
//                             const gsFunctionSet<T>                              & mp,
//                             const gsFunctionSet<T>                              & mp_def,
//                             const gsVector<T>                                   & E11,
//                             const gsVector<T>                                   & E22,
//                             const gsVector<T>                                   & G12,
//                             const gsVector<T>                                   & nu12,
//                             const gsVector<T>                                   & nu21,
//                             const gsVector<T>                                   & thickness,
//                             const gsVector<T>                                   & phi           )
//                             :
//                             m_patches(&mp),
//                             Base::m_defpatches(&mp_def)
// {
//     m_E11       = new gsConstantFunction<T>(E11,3);
//     m_E22       = new gsConstantFunction<T>(E22,3);
//     m_G12       = new gsConstantFunction<T>(G12,3);
//     m_nu12      = new gsConstantFunction<T>(nu12,3);
//     m_nu21      = new gsConstantFunction<T>(nu21,3);
//     m_thickness = new gsConstantFunction<T>(thickness,3);
//     m_phi       = new gsConstantFunction<T>(phi,3);
//     _initialize();
// }

// template <short_t dim, class T>
// gsMaterialMatrixComposite<dim,T>::gsMaterialMatrixComposite(
//                             const gsFunctionSet<T>                              & mp,
//                             const gsFunctionSet<T>                              & mp_def,
//                             const gsVector<T>                                   & E11,
//                             const gsVector<T>                                   & E22,
//                             const gsVector<T>                                   & G12,
//                             const gsVector<T>                                   & nu12,
//                             const gsVector<T>                                   & nu21,
//                             const gsVector<T>                                   & thickness,
//                             const gsVector<T>                                   & phi,
//                             const gsVector<T>                                   & rho           )
//                             :
//                             m_patches(&mp),
//                             Base::m_defpatches(&mp_def)
// {
//     m_E11       = new gsConstantFunction<T>(E11,3);
//     m_E22       = new gsConstantFunction<T>(E22,3);
//     m_G12       = new gsConstantFunction<T>(G12,3);
//     m_nu12      = new gsConstantFunction<T>(nu12,3);
//     m_nu21      = new gsConstantFunction<T>(nu21,3);
//     m_thickness = new gsConstantFunction<T>(thickness,3);
//     m_phi       = new gsConstantFunction<T>(phi,3);
//     m_rho       = new gsConstantFunction<T>(rho,3);
//     _initialize();
// }
// ---------------------------------------------------------------------------------------------------

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_defaultOptions()
{
    m_options.addInt("NumGauss","Number of Gaussian points through thickness",4);
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::info() const
{
    gsInfo  <<"---------------------------------------------------------------------\n"
            <<"---------------------Hyperelastic Material Info----------------------\n"
            <<"---------------------------------------------------------------------\n\n";

    gsInfo  <<"Material model: \t";
    gsInfo<<" Composite Saint-Venant Kirchhoff";
    gsInfo<<"\n";

    gsInfo  <<"---------------------------------------------------------------------\n\n";

}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_initializeParameters()
{
    GISMO_ASSERT(m_pars.size()==6,"parameters should contain 6 functions, not "<<m_pars.size());
    m_E11   = m_pars[0];
    m_E22   = m_pars[1];
    m_G12   = m_pars[2];
    m_nu12  = m_pars[3];
    m_nu21  = m_pars[4];
    m_phi   = m_pars[5];
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_initialize()
{
    // Set default options
    this->_defaultOptions();

    // set flags
    m_map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;

    m_nLayers = m_thickness->targetDim();
}


template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::_computePoints(const index_t patch, const gsMatrix<T> & u) const
{
    m_map.flags = NEED_JACOBIAN | NEED_DERIV | NEED_NORMAL | NEED_VALUE | NEED_DERIV2;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(patch)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
    this->_computeMetricUndeformed();

    if (Base::m_defpatches->nPieces()!=0)
    {
        m_map_def.flags = m_map.flags;
        m_map_def.points = u;
        static_cast<const gsFunction<T>&>(Base::m_defpatches->piece(patch)).computeMap(m_map_def); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()
        this->_computeMetricDeformed();
    }

    GISMO_ASSERT(m_E11->targetDim()==m_E22->targetDim(),"Size of vectors of Youngs Moduli not equal: " << m_E11->targetDim()<<" & "<<m_E22->targetDim());
    GISMO_ASSERT(m_nu12->targetDim()==m_nu21->targetDim(),"Size of vectors of Poisson Ratios is not equal: " << m_nu12->targetDim()<<" & "<<m_nu21->targetDim());
    GISMO_ASSERT(m_E11->targetDim()==m_nu12->targetDim(),"Size of vectors of Youngs Moduli and Poisson Ratios is not equal: " << m_E11->targetDim()<<" & "<<m_nu12->targetDim());
    GISMO_ASSERT(m_E11->targetDim()==m_G12->targetDim(),"Size of vectors of Youngs Moduli and Shear Moduli is not equal: " << m_E11->targetDim()<<" & "<<m_G12->targetDim());
    GISMO_ASSERT(m_thickness->targetDim()==m_phi->targetDim(),"Size of vectors of thickness and angles is not equal: " << m_thickness->targetDim()<<" & "<<m_phi->targetDim());
    GISMO_ASSERT(m_E11->targetDim()==m_thickness->targetDim(),"Size of vectors of material properties and laminate properties is not equal: " << m_E11->targetDim()<<" & "<<m_thickness->targetDim());
    GISMO_ASSERT(m_E11->targetDim()!=0,"No laminates defined");


    // Compute properties per ply
    m_Tmat.resize(m_nLayers,u.cols());
    m_E1mat.resize(m_nLayers,u.cols());
    m_E2mat.resize(m_nLayers,u.cols());
    m_G12mat.resize(m_nLayers,u.cols());
    m_nu12mat.resize(m_nLayers,u.cols());
    m_nu21mat.resize(m_nLayers,u.cols());
    m_phiMat.resize(m_nLayers,u.cols());
    m_rhoMat.resize(m_nLayers,u.cols());

    m_E1mat = m_E11 ->eval(m_map.values[0]);
    m_E2mat = m_E22 ->eval(m_map.values[0]);
    m_G12mat = m_G12 ->eval(m_map.values[0]);
    m_nu12mat = m_nu12->eval(m_map.values[0]);
    m_nu21mat = m_nu21->eval(m_map.values[0]);
    m_Tmat = m_thickness->eval(m_map.values[0]);
    m_phiMat = m_phi->eval(m_map.values[0]);
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT(m_thickness->targetDim()==m_rho->targetDim(),"Size of vectors of thickness and densities is not equal: " << m_thickness->targetDim()<<" & "<<m_rho->targetDim());

    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    result.resize(1, u.cols());
    m_thickness->eval_into(m_map.values[0], m_Tmat);
    m_rho->eval_into(m_map.values[0], m_rhoMat);

    for (index_t i = 0; i != m_nLayers; ++i) // layers
        for (index_t k = 0; k != u.cols(); ++k) // points
            result(i,k) = m_Tmat(i,k)*m_rhoMat(i,k);
}

template <short_t dim, class T >
void gsMaterialMatrixComposite<dim,T>::thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    m_map.flags = NEED_VALUE;
    m_map.points = u;
    static_cast<const gsFunction<T>&>(m_patches->piece(0)   ).computeMap(m_map); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

    m_thickness->eval_into(m_map.values[0], result);
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
    gsMatrix<T,3,3> Dmat, Cmat;

    T t_tot, z_mid, t_temp, dz;
    T E1, E2, G12, nu12, nu21, t, phi;

    for (index_t k = 0; k != u.cols(); ++k)
    {
        Cmat.setZero();
        Dmat.setZero();
        this->_getMetric(k,0.0); // on point i, with height 0.0

        // Compute total thickness (sum of entries)
        t_tot = m_Tmat.col(k).sum();

        // compute mid-plane height of total plate
        z_mid = t_tot / 2.0;

        // now we use t_temp to add the thickness of all plies iteratively
        t_temp = 0.0;
        t = 0.0;

        for (index_t i = 0; i != m_nLayers; ++i) // loop over laminates
        {
            // Lookup all quantities
            E1 = m_E1mat(i,k);
            E2 = m_E2mat(i,k);
            G12 = m_G12mat(i,k);
            nu12 = m_nu12mat(i,k);
            nu21 = m_nu21mat(i,k);
            t = m_Tmat(i,k);
            phi = m_phiMat(i,k);

            Dmat = _computeMatrix(E1,E2,G12,nu12,nu21,phi);

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
        a1 = m_gcov_ori.col(0);               a2 = m_gcov_ori.col(1);
        ac1 = m_gcon_ori.col(0);              ac2 = m_gcon_ori.col(1);
        e1 = m_gcov_ori.col(0).normalized();  e2 = m_gcon_ori.col(1).normalized();

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
            Eij = 0.5*(m_Acov_def - m_Acov_ori);
        else if (out == MaterialOutput::VectorM) // To be used with multiplyZ_into
            Eij = (m_Bcov_ori - m_Bcov_def);
        else if (out == MaterialOutput::Generic) // To be used with multiplyLinZ_into or integrateZ_into
            Eij = 0.5*(m_Acov_def - m_Acov_ori) + z*(m_Bcov_ori - m_Bcov_def);
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
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::_computeMatrix(const T E11, const T E22, const T G12, const T nu12, const T nu21, const T phi) const
{
    GISMO_ASSERT(nu21*E11 == nu12*E22, "No symmetry in material properties for ply! (nu12*E2!=nu21*E1):\n"<<
            "\tnu12 = "<<nu12<<"\t E22 = "<<E22<<"\t nu12*E22 = "<<nu12*E22<<"\n"
          <<"\tnu21 = "<<nu21<<"\t E11 = "<<E11<<"\t nu21*E11 = "<<nu21*E11);

    gsMatrix<T,3,3> Dmat, Tmat;
    // Fill material matrix
    Dmat(0,0) = E11 / (1-nu12*nu21);
    Dmat(1,1) = E22 / (1-nu12*nu21);
    Dmat(2,2) = G12;
    Dmat(0,1) = nu21*E11 / (1-nu12*nu21);
    Dmat(1,0) = nu12*E22 / (1-nu12*nu21);
    Dmat(2,0) = Dmat(0,2) = Dmat(2,1) = Dmat(1,2) = 0.0;

    // Make transformation matrix
    Tmat(0,0) = Tmat(1,1) = math::pow(math::cos(phi),2);
    Tmat(0,1) = Tmat(1,0) = math::pow(math::sin(phi),2);
    Tmat(2,0) = Tmat(0,2) = Tmat(2,1) = Tmat(1,2) = math::sin(phi) * math::cos(phi);
    Tmat(2,0) *= -2.0;
    Tmat(2,1) *= 2.0;
    Tmat(1,2) *= -1.0;
    Tmat(2,2) = math::pow(math::cos(phi),2) - math::pow(math::sin(phi),2);
    // Compute laminate stiffness matrix
    return Tmat.transpose() * Dmat * Tmat;
}

template <short_t dim, class T>
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::_cart2cov(const gsVector<T> a1, const gsVector<T> a2, const gsVector<T> e1, const gsVector<T> e2) const
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
gsMatrix<T> gsMaterialMatrixComposite<dim,T>::_con2cart(const gsVector<T> ac1, const gsVector<T> ac2, const gsVector<T> e1, const gsVector<T> e2) const
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
//-----------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixComposite<dim,T>::_computeMetricDeformed() const
{
    _computeMetricDeformed_impl<dim>();
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixComposite<dim,T>::_computeMetricDeformed_impl() const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();
    m_Bcov_def_mat.resize(4,m_map_def.points.cols());    m_Bcov_def_mat.setZero();

    m_acov_def_mat.resize(2*3,m_map_def.points.cols());    m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*3,m_map_def.points.cols());    m_acon_def_mat.setZero();
    m_ncov_def_mat.resize(2*3,m_map_def.points.cols());    m_ncov_def_mat.setZero();

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        m_acov_def_mat.reshapeCol(k,3,2)   = m_map_def.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,3,2);
        acov = m_map_def.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map_def.deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map_def.normal(k).normalized();

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_Bcov_def_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_def_mat.reshapeCol(k,2,2);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

        gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_def_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
    }
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixComposite<dim,T>::_computeMetricDeformed_impl() const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_def_mat.resize(4,m_map_def.points.cols());    m_Acov_def_mat.setZero();
    m_Acon_def_mat.resize(4,m_map_def.points.cols());    m_Acon_def_mat.setZero();

    m_acov_def_mat.resize(2*2,m_map_def.points.cols());    m_acov_def_mat.setZero();
    m_acon_def_mat.resize(2*2,m_map_def.points.cols());    m_acon_def_mat.setZero();

    for (index_t k=0; k!= m_map_def.points.cols(); k++)
    {
        m_acov_def_mat.reshapeCol(k,2,2)   = m_map_def.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_def_mat.reshapeCol(k,2,2);
        acov = m_map_def.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_def_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_def_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_def_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_def_mat.reshapeCol(k,2,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixComposite<dim,T>::_computeMetricUndeformed() const
{
    _computeMetricUndeformed_impl<dim>();
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixComposite<dim,T>::_computeMetricUndeformed_impl() const
{
    gsMatrix<T> deriv2;
    gsMatrix<T,3,1> normal;
    gsMatrix<T,2,2> mixedB;
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();
    m_Bcov_ori_mat.resize(4,m_map.points.cols());    m_Bcov_ori_mat.setZero();

    m_acov_ori_mat.resize(2*3,m_map.points.cols());    m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*3,m_map.points.cols());    m_acon_ori_mat.setZero();
    m_ncov_ori_mat.resize(2*3,m_map.points.cols());    m_ncov_ori_mat.setZero();

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        m_acov_ori_mat.reshapeCol(k,3,2)   = m_map.jacobian(k);
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,3,2);
        acov = m_map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);

        // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
        deriv2    = m_map.deriv2(k);
        deriv2.resize(3,3);
        normal    = m_map.normal(k).normalized();

        tmp(0,0) = deriv2.row(0).dot(normal);
        tmp(1,1) = deriv2.row(1).dot(normal);
        tmp(0,1) = tmp(1,0) = deriv2.row(2).dot(normal);

        m_Bcov_ori_mat.reshapeCol(k,2,2) = tmp;
        gsAsMatrix<T,Dynamic,Dynamic> metricBcov = m_Bcov_ori_mat.reshapeCol(k,2,2);

        // Mixed tensor
        for (index_t i=0; i < 2; i++)
            for (index_t j=0; j < 2; j++)
                mixedB(i,j) = metricAcon(i,0)*metricBcov(0,j) + metricAcon(i,1)*metricBcov(1,j);

        gsAsMatrix<T,Dynamic,Dynamic> ncov = m_ncov_ori_mat.reshapeCol(k,3,2);
        for (index_t i=0; i < 2; i++)
            ncov.col(i)     = -mixedB(0,i)*acov.col(0) -mixedB(1,i)*acov.col(1);
    }
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixComposite<dim,T>::_computeMetricUndeformed_impl() const
{
    gsMatrix<T,2,2> tmp;

    m_Acov_ori_mat.resize(4,m_map.points.cols());    m_Acov_ori_mat.setZero();
    m_Acon_ori_mat.resize(4,m_map.points.cols());    m_Acon_ori_mat.setZero();

    m_acov_ori_mat.resize(2*2,m_map.points.cols());    m_acov_ori_mat.setZero();
    m_acon_ori_mat.resize(2*2,m_map.points.cols());    m_acon_ori_mat.setZero();

    for (index_t k=0; k!= m_map.points.cols(); k++)
    {
        gsAsMatrix<T,Dynamic,Dynamic> acov = m_acov_ori_mat.reshapeCol(k,2,2);
        acov = m_map.jacobian(k);

        tmp = acov.transpose() * acov;
        m_Acov_ori_mat.reshapeCol(k,2,2) = tmp;
        m_Acon_ori_mat.reshapeCol(k,2,2) = tmp.inverse();

        gsAsMatrix<T,Dynamic,Dynamic> metricAcon = m_Acon_ori_mat.reshapeCol(k,2,2);

        gsAsMatrix<T,Dynamic,Dynamic> acon = m_acon_ori_mat.reshapeCol(k,2,2);
        for (index_t i=0; i < 2; i++)
            acon.col(i)     = metricAcon(i,0)*acov.col(0) + metricAcon(i,1)*acov.col(1);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixComposite<dim,T>::_getMetric(index_t k, T z) const
{
    this->_getMetricDeformed(k,z);
    this->_getMetricUndeformed(k,z);

    T ratio = m_Gcov_def.determinant() / m_Gcov_ori.determinant();
    GISMO_ENSURE(ratio > 0, "Jacobian determinant is negative! det(Gcov_def) = "<<m_Gcov_def.determinant()<<"; det(Gcov_ori) = "<<m_Gcov_ori.determinant());
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixComposite<dim,T>::_getMetricDeformed(index_t k, T z) const
{
    _getMetricDeformed_impl<dim>(k,z);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixComposite<dim,T>::_getMetricDeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_ncov_def_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    m_Bcov_def = m_Bcov_def_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,3,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,3,2);
    m_ncov_def = m_ncov_def_mat.reshapeCol(k,3,2);
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def - 2.0 * z * m_Bcov_def + z*z * m_ncov_def.transpose()*m_ncov_def;
    m_Gcov_def(2,2) = 1.0;
    m_Gcon_def = m_Gcov_def.inverse();
    // Compute full basis
    gsMatrix<T,3,1> normal = m_map_def.normal(k).normalized();
    m_gcov_def.leftCols(2) = m_acov_def + z * m_ncov_def;
    m_gcov_def.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_def.col(c) = m_Gcon_def(c,0) * m_gcov_def.col(0)
                            + m_Gcon_def(c,1) * m_gcov_def.col(1)
                            + m_Gcon_def(c,2) * m_gcov_def.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_def_L = m_Gcov_def;
    // m_Gcov_def.block(0,0,2,2) -= z*z * m_ncov_def.transpose()*m_ncov_def;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixComposite<dim,T>::_getMetricDeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_def_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_def_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_def_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_def = m_Acov_def_mat.reshapeCol(k,2,2);
    m_Acon_def = m_Acon_def_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_def = m_acov_def_mat.reshapeCol(k,2,2);
    m_acon_def = m_acon_def_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_Gcov_def.setZero();
    m_Gcov_def.block(0,0,2,2)= m_Acov_def;
    m_Gcov_def(2,2) = 1.0;
    m_Gcon_def = m_Gcov_def.inverse();
    // Compute full basis
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_gcov_def.setZero();
    m_gcov_def.block(0,0,2,2) = m_acov_def;
    m_gcov_def.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_def.col(c) = m_Gcon_def(c,0) * m_gcov_def.col(0)
                            + m_Gcon_def(c,1) * m_gcov_def.col(1)
                            + m_Gcon_def(c,2) * m_gcov_def.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_def_L = m_Gcov_def;
    // m_Gcov_def.block(0,0,2,2) -= z*z * m_ncov_def.transpose()*m_ncov_def;
}

//--------------------------------------------------------------------------------------------------------------------------------------

template <short_t dim, class T>
void gsMaterialMatrixComposite<dim,T>::_getMetricUndeformed(index_t k, T z) const
{
    _getMetricUndeformed_impl<dim>(k,z);
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==3, void>::type
gsMaterialMatrixComposite<dim,T>::_getMetricUndeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Bcov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_ncov_ori_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    m_Bcov_ori = m_Bcov_ori_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,3,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,3,2);
    m_ncov_ori = m_ncov_ori_mat.reshapeCol(k,3,2);
    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori - 2.0 * z * m_Bcov_ori + z*z * m_ncov_ori.transpose()*m_ncov_ori;
    m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori = m_Gcov_ori.inverse();
    // Compute full basis
    gsMatrix<T,3,1> normal = m_map.normal(k).normalized();
    m_gcov_ori.block(0,0,3,2) = m_acov_ori + z * m_ncov_ori;
    m_gcov_ori.col(2) = normal;
    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_ori.col(c) = m_Gcon_ori(c,0) * m_gcov_ori.col(0)
                            + m_Gcon_ori(c,1) * m_gcov_ori.col(1)
                            + m_Gcon_ori(c,2) * m_gcov_ori.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_ori_L = m_Gcov_ori;
    // m_Gcov_ori.block(0,0,2,2) -= z*z * m_ncov_ori.transpose()*m_ncov_ori;
}

template <short_t dim, class T>
template <short_t _dim>
typename std::enable_if<_dim==2, void>::type
gsMaterialMatrixComposite<dim,T>::_getMetricUndeformed_impl(index_t k, T z) const
{
    GISMO_ENSURE(m_Acov_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_Acon_ori_mat.cols()!=0,"Is the metric initialized?");
    GISMO_ENSURE(m_acov_ori_mat.cols()!=0,"Is the basis initialized?");
    GISMO_ENSURE(m_acon_ori_mat.cols()!=0,"Is the basis initialized?");

    // metrics
    m_Acov_ori = m_Acov_ori_mat.reshapeCol(k,2,2);
    m_Acon_ori = m_Acon_ori_mat.reshapeCol(k,2,2);
    // basis vectors
    m_acov_ori = m_acov_ori_mat.reshapeCol(k,2,2);
    m_acon_ori = m_acon_ori_mat.reshapeCol(k,2,2);
    // Compute full metric
    m_Gcov_ori.setZero();
    m_Gcov_ori.block(0,0,2,2)= m_Acov_ori;
    m_Gcov_ori(2,2) = 1.0;
    m_Gcon_ori = m_Gcov_ori.inverse();
    // Compute full basis
    gsMatrix<T,3,1> normal;
    normal << 0,0,1;
    m_gcov_ori.setZero();
    m_gcov_ori.block(0,0,2,2) = m_acov_ori;
    m_gcov_ori.col(2) = normal;

    for (index_t c = 0; c!=3; c++)
    {
        m_gcon_ori.col(c) = m_Gcon_ori(c,0) * m_gcov_ori.col(0)
                            + m_Gcon_ori(c,1) * m_gcov_ori.col(1)
                            + m_Gcon_ori(c,2) * m_gcov_ori.col(2);
    }

    // // create a Linearised covariant tensor for SvK models
    // m_Gcov_ori_L = m_Gcov_ori;
    // m_Gcov_ori.block(0,0,2,2) -= z*z * m_ncov_ori.transpose()*m_ncov_ori;
}

} // end namespace
