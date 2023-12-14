/** @file gsMaterialMatrixLinear.h

    @brief Provides linear material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/src/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{


/**
 * @brief      This class defines a linear material
 *
 * @tparam     dim   The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T     Real type
 *
 * @ingroup    MaterialMatrix
 *
 */
template <  short_t dim,
            class T     >
class gsMaterialMatrixLinear : public gsMaterialMatrixBaseDim<dim,T>
{
public:

    typedef T Scalar_t;

    GISMO_CLONE_FUNCTION(gsMaterialMatrixLinear)

    using Base = gsMaterialMatrixBaseDim<dim,T>;

    typedef typename Base::function_ptr function_ptr;

    enum {Linear=1};

    /**
     * @brief      Empty constructor
     */
    gsMaterialMatrixLinear();

    /**
     * @brief      Constructor without material parameters
     *
     * @param[in]  mp             Original geometry
     * @param[in]  thickness      Thickness function
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness);

    /**
     * @brief      Constructor without deformed multipatch and density
     *
     * @param[in]  mp             Original geometry
     * @param[in]  thickness      Thickness function
     * @param[in]  YoungsModulus  The youngs modulus
     * @param[in]  PoissonRatio   The poisson ratio
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness,
                        const gsFunctionSet<T> & YoungsModulus,
                        const gsFunctionSet<T> & PoissonRatio);

    /**
     * @brief      Full constructor
     *
     * @param[in]  mp             Original geometry
     * @param[in]  thickness      Thickness function
     * @param[in]  YoungsModulus  The youngs modulus
     * @param[in]  PoissonRatio   The poisson ratio
     * @param[in]  Density        The density
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness,
                        const gsFunctionSet<T> & YoungsModulus,
                        const gsFunctionSet<T> & PoissonRatio,
                        const gsFunctionSet<T> & Density);

protected:
    /**
     * @brief      Full constructor
     *
     * @param[in]  mp             Original geometry
     * @param[in]  thickness      Thickness function
     * @param[in]  YoungsModulus  The youngs modulus
     * @param[in]  PoissonRatio   The poisson ratio
     * @param[in]  Density        The density
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> * mp,
                        const gsFunctionSet<T> * thickness,
                        const gsFunctionSet<T> & YoungsModulus,
                        const gsFunctionSet<T> & PoissonRatio,
                        const gsFunctionSet<T> * Density);
public:
    /**
     * @brief      Constructor without density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars);

    /**
     * @brief      Full constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars,
                        const gsFunctionSet<T> & Density);

    /**
     * @brief      Constructor without density and multipatch
     *
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixLinear(    const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars);

    /**
     * @brief      Constructor without multipatch
     *
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> & thickness,
                        const std::vector<gsFunctionSet<T> *> &pars,
                        const gsFunctionSet<T> & Density);

protected:
    /**
     * @brief      Full constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrixLinear(   const gsFunctionSet<T> * mp,
                        const gsFunctionSet<T> * thickness,
                        const std::vector<gsFunctionSet<T> *> &pars,
                        const gsFunctionSet<T> * Density);
public:

    gsMaterialMatrixLinear( const gsMaterialMatrixLinear<dim,T> & other);

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isMatIntegrated() const override {return MatIntegration::Constant; }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isVecIntegrated() const override {return MatIntegration::Constant; }

    /// See \ref gsMaterialMatrixBase for details
    void defaultOptions();

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix (const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_dmatrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector (const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyVector (const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T>& u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const override;
    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector_C(const gsMatrix<T> & Cmat, const index_t patch, const gsVector<T> & u, const T z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstressDir(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyPStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_stress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_CauchyStress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const override;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_detF (const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const override;

    /// Sets the YoungsModulus
    void setYoungsModulus(const gsFunctionSet<T> & YoungsModulus) { Base::setParameter(0,YoungsModulus); }

    /// Gets the YoungsModulus
    const function_ptr getYoungsModulus() const { return Base::getParameter(0); }

    /// Sets the Poisson's Ratio
    void setPoissonsRatio(const gsFunctionSet<T> & PoissonsRatio) { Base::setParameter(1,PoissonsRatio); }

    /// Gets the Poisson's Ratio
    const function_ptr getPoissonsRatio() const { return Base::getParameter(1); }

    /// See \ref gsMaterialMatrixBase for details
    std::ostream &print(std::ostream &os) const override;

    gsMatrix<T> S(const gsMatrix<T> & strain) const;

    gsMatrix<T> C(const gsMatrix<T> & strain) const;

    /// Computes the vector S as function of the deformation tensor C=FTF
    gsMatrix<T> S(const gsMatrix<T> & C, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;

    /// Computes the matrix C as function of the deformation tensor C=FTF
    gsMatrix<T> C(const gsMatrix<T> & C, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;

    /// Computes the derivative of the matrix C as function of the deformation tensor C=FTF
    gsMatrix<T> dC(const gsMatrix<T> & C, const index_t patch, const gsMatrix<T> & u, const gsMatrix<T> & z, enum MaterialOutput out) const;

public:
    /// Shared pointer for gsMaterialMatrixLinear
    typedef memory::shared_ptr< gsMaterialMatrixLinear > Ptr;

    /// Unique pointer for gsMaterialMatrixLinear
    typedef memory::unique_ptr< gsMaterialMatrixLinear > uPtr;

protected:
    /**
     * @brief      Initializes the object.
     *
     * Initializes options, flags and defines the number of parameters
     *
     */
    void _initialize();

protected:
    /**
     * @brief      Computes the linear material matrix entry with indices \a i \a j \a k \a l
     *
     * The entry is computed by \f$ \mathcal{C}^{ijkl} = \frac{2\lambda\mu}{\lambda+2\mu}a^{ij}a^{kl} + \mu (a^{ik}a^{jl} + a^{il}a^{jk})\f$ where \f$\lambda\f$ and \f$\mu\f$
     * are the Lamé parameters and \f$a^{ij}=\mathbf{a}^i\cdot \mathbf{a}^j$ with \f$\mathbf{a}^i\f$ is the contravariant vector \f$ i \f$
     *
     * @param[in]  i,j,k,l  Indices
     *
     * @return     Cijkl
     */
    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;

    /**
     * @brief      Computes the linear force/moment entry with indices \a i \a j at height z
     *
     * Computes the thickness-integrated stress tensor, i.e. the normal force (0th thickness-moment) or the bending moment (1st thickness-moment).
     * Sij is computed as \f$ \mathcal{C}^{ijkl} : \mathbf{E}_{ij} \f$ where \f$\mathbf{E}_{ij} = a_{ij} z b_{ij}\f$ with \f$ a_{ij}\f$ the in-plane metric and \f$b_{ij}\f$ the curvature.
     * According to this definition, MaterialOutput out==VectorN would return the 0th moment, hence the \f$a_{ij}\f$-part would be relevant and for VectorM the \f$b_{ij}$ part is relevant.
     * The term could also be integrated (out==Generalized).
     *
     * @param[in]  i,j   Indices
     * @param[in]  z     Through-thickness coordinate
     * @param[in]  out   Output specification
     *
     * @return     Sij
     */
    T _Sij    (const index_t i, const index_t j, const gsMatrix<T> & z) const;


    using Base::_getMetric;

    /**
     * @brief      Computes the strain tensor
     *
     * @param[in]  z     Through-thickness coordinate
     * @param[in]  out   Output specification
     *
     * @return     E
     */
    gsMatrix<T> _E      (const T z, enum MaterialOutput out) const;

protected:
    // constructor
    using Base::m_thickness;
    using Base::m_pars;
    using Base::m_density;

    // Geometric data
    using Base::m_data;

    using Base::m_options;

};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMaterialMatrixLinear
   */
  void pybind11_init_gsMaterialMatrixLinear2(pybind11::module &m);
  void pybind11_init_gsMaterialMatrixLinear3(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixLinear.hpp)
#endif
