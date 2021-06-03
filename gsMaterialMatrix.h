/** @file gsMaterialMatrix.h

    @brief Provides hyperelastic material matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{


/**
 * @brief      This class defines hyperelastic material matrices
 *
 * @tparam     dim    The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T      Real type
 * @tparam     matId  Encoded ID for material (see \ref encodeMat_id)
 * @tparam     comp   Compressibility flag
 * @tparam     mat    Material flag (see \ref Material)
 * @tparam     imp    Implementation flag (see \ref Implementation)
 *
 * @ingroup    MaterialMatrix
 */
template <  short_t dim,
            class T,
            short_t matId,
            bool comp,
            enum Material mat = decodeMat_id<matId>::material,
            enum Implementation imp  = decodeMat_id<matId>::implementation
         >
class gsMaterialMatrix :    public gsMaterialMatrixBaseDim<dim,T>
{
public:
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /**
     * @brief      Constructor without density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    /**
     * @brief      Full constructor
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     * @param[in]  Density    Density function
     */
    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars,
                        const gsFunction<T> & Density);

    /// Destructor
    gsMaterialMatrix() { }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isMatIntegrated() const {return MatIntegration::NotIntegrated; }

    /// See \ref gsMaterialMatrixBase for details
    inline enum MatIntegration isVecIntegrated() const {return MatIntegration::NotIntegrated; }

    /// See \ref gsMaterialMatrixBase for details
    gsOptionList & options() {return m_options;}

    /// See \ref gsMaterialMatrixBase for details
    void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
        m_numPars = m_pars.size();
    }

    /// See \ref gsMaterialMatrixBase for details
    void info() const;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

protected:
    /**
     * @brief      Initializes the object.
     *
     * Initializes options, flags and defines the number of parameters
     *
     */
    void _initialize();

    /**
     * @brief      Sets default options
     */
    void _defaultOptions();

private:
    /**
     * @brief      Implementation of stretch_into, specialization for compressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<_com, void>::type _stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /**
     * @brief      Implementation of stretch_into, specialization for incompressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<!_com, void>::type _stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /**
     * @brief      Implementation of stretchDir_into, specialization for compressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<_com, void>::type _stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /**
     * @brief      Implementation of stretchDir_into, specialization for incompressible materials
     *
     * @param[in]  u       The in-plane coordinates to be evaluated on
     * @param      result  The result
     *
     * @tparam     _com    Compressibility parameter (true: compressible, false: incompressible)
     *
     */
    template<bool _com>
    typename std::enable_if<!_com, void>::type _stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

protected:
    /**
     * @brief      Evalluates the incompressible material matrix
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The matrices (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval_Incompressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evalluates the incompressible stress vector
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The vectors (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval_Incompressible_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evalluates the incompressible principal stresses
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The principal stresses (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval_Incompressible_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evalluates the compressible material matrix
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The matrices (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval_Compressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evalluates the compressible stress vector
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The vectors (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval_Compressible_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Evalluates the compressible principal stresses
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @return     The principal stresses (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    gsMatrix<T> _eval_Compressible_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

private:

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the material matrix, specialization for SvK material (incompressible)
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The material matrix (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the material matrix, specialization compressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The material matrix (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the material matrix, specialization for incompressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The material matrix (9x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for SvK material (incompressible)
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for incompressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for compressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The stress vector (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for SvK material (incompressible)
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The principal stresses (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for incompressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The principal stresses (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    /**
     * @brief      Implementation of the 3D (in-plane+thickness) evaluator of the stress vector, specialization for compressible models
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane coordinates to be evaluated on
     * @param[in]  z      The through-thickness coordinate
     *
     * @tparam     _mat   The material model
     * @tparam     _com   Compressibility parameter (true: compressible, false: incompressible)
     *
     * @return     The principal stresses (3x1 per column), stored per thickness per in-plane point [(u1,t1) (u2,t1),...,(u1,t2),(u2,t2)]
     */
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;


protected:

    /**
     * @brief      Returns an entry of the material tensor C for incompressible materials
     *
     * @param[in]  i,j,k,l     The indices
     *
     * @return     C^{ijkl}
     */
    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;

    /**
     * @brief      Returns an entry of the material tensor C for compressible materials without static condensation
     *
     * @param[in]  i,j,k,l     The indices
     * @param[in]  c           The deformation tensor
     * @param[in]  cinv        The inverse of the deformation tensor
     *
     * @return     C^{ijkl}
     */
    T _Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Returns an entry of the material tensor C for compressible materials, with static condensation included
     *
     * @param[in]  i,j,k,l     The indices
     * @param[in]  c           The deformation tensor
     * @param[in]  cinv        The inverse of the deformation tensor
     *
     * @return     C^{ijkl}
     */
    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Returns an entry of the stress tensor S for incompressible materials
     *
     * @param[in]  i,j     The indices
     *
     * @return     S^{ij}
     */
    T _Sij    (const index_t i, const index_t j) const;

    /**
     * @brief      Returns an entry of the stress tensor S for compressible materials
     *
     * @param[in]  i,j     The indices
     * @param[in]  c       The deformation tensor
     * @param[in]  cinv    The inverse of the deformation tensor
     *
     * @return     S^{ij}
     */
    T _Sij    (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Returns an entry of the diagonal of the stress tensor S for incompressible materials
     *
     * @param[in]  i       The index
     *
     * @return     S^{ii}
     */
    T _Sii    (const index_t i) const;

    /**
     * @brief      Returns an entry of the diagonal of the stress tensor S for incompressible materials
     *
     * @param[in]  i       The index
     * @param[in]  c       The deformation tensor
     *
     * @return     S^{ii}
     */
    T _Sii    (const index_t i, const gsMatrix<T> & c) const;

private:
    ///////////////////////////////////////////////////////
    // Incompressible formulations                       //
    ///////////////////////////////////////////////////////
    /// Specialization for incompressible Cijkl(i,j,k,l) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH  && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR  && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Cijkl(i,j,k,l) for Extended NH materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Cijkl(i,j,k,l) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Specialization for incompressible Cijkl(i,j,k,l) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;


    ///////////////////////////////////////////////////////
    // Compressible formulations                         //
    ///////////////////////////////////////////////////////
    /**
     * @brief      Specialization for compressible Cijkl(i,j,k,l,c,cinv) for all materials implemented generally
     *
     * @note        The specialization is required because static condensation for spectrally implemented material models is performed before the transformation. For all other implementation, static condensation is performed using the Cijkl3D function.
     *
     */
    template<enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral,   T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl(i,j,k,l,c,cinv) for all materials implemented not spectrally
    template<enum Implementation _imp>
    typename std::enable_if<!(_imp==Implementation::Spectral),T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for Extended NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Cijkl3D(i,j,k,l,c,cinv) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;


    /// Specialization for incompressible Sij(i,j) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Sij(i,j) for Extended NH materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for incompressible Sij(i,j) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    /// Specialization for incompressible Sij(i,j) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Sij_impl(const index_t i, const index_t j) const;


    /// Specialization for compressible Sij(i,j,c,cinv) for SvK materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for compressible Sij(i,j,c,cinv) for NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for MR materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for OG materials implemented analytically (not implemented)
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Specialization for compressible Sij(i,j,c,cinv) for Extended NH materials implemented analytically
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for all materials implemented spectrally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Specialization for compressible Sij(i,j,c,cinv) for all materials implemented generally
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

protected:
    /**
     * @brief      Provides the derivative of the incompressible strain energy density function w.r.t. component C_{ij} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C
     *
     * @return     The derivative of psi w.r.t. C_{ij}
     */
    T _dPsi   (const index_t i, const index_t j) const;

    /**
     * @brief      Provides the derivative of the compressible strain energy density function w.r.t. component C_{ij} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The derivative of psi w.r.t. C_{ij}
     */
    T _dPsi   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the derivative of the volumetric part of the compressible strain energy density function w.r.t. component C_{ij} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The derivative of psi w.r.t. C_{ij}
     */
    T _dPsi_vol(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the second (mixed) derivative of the incompressible strain energy density function w.r.t. components C_{ij} and C_{kl} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C_{ij}
     * @param[in]  k,l   The indices of C_{kl}
     *
     * @return     The second (mixed) derivative of psi w.r.t. C_{ij} and C_{kl}
     */
    T _d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l) const;

    /**
     * @brief      Provides the second (mixed) derivative of the compressible strain energy density function w.r.t. components C_{ij} and C_{kl} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C_{ij}
     * @param[in]  k,l   The indices of C_{kl}
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The second (mixed) derivative of psi w.r.t. C_{ij} and C_{kl}
     */
    T _d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the second (mixed) derivative of the volumetric part of the compressible strain energy density function w.r.t. components C_{ij} and C_{kl} of the deformation tensor
     *
     * @param[in]  i,j   The indices of C_{ij}
     * @param[in]  k,l   The indices of C_{kl}
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The second (mixed) derivative of psi w.r.t. C_{ij} and C_{kl}
     */
    T _d2Psi_vol(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /**
     * @brief      Provides the derivative of the first invariant compressible materials w.r.t. component C_{ij} of the deformation tensor
     *
     * @note       \f$I_1=\trace{\mathbf{C}}\f$
     *
     * @param[in]  i,j   The indices of C
     *
     * @return     The derivative the first invariant w.r.t. C_{ij}
     */
    T _dI_1   (const index_t i, const index_t j) const;

    /**
     * @brief      Provides the derivative of the second invariant for compressible materials w.r.t. component C_{ij} of the deformation tensor
     *
     * @note       \f$I_1=\frac{1}{2}\left( \trace{\mathbf{C}^2} + (\trace{\mathbf{C}})^2 \right)\f$
     *
     * @param[in]  i,j   The indices of C
     * @param[in]  c     The deformation tensor
     * @param[in]  cinv  The inverse of the deformation tensor
     *
     * @return     The derivative the second invariant w.r.t. C_{ij}
     */
    T _dI_2   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;



private:
    /*
        Note: when adding a new implementation, make sure to add the condition in the 'other' implementation : !(_mat==xx || ...)
    */

    /// Implementation of _dPsi(i,j) for NH materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _dPsi_impl(const index_t i, const index_t j) const;

    /// Implementation of _dPsi(i,j) for MR materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _dPsi_impl(const index_t i, const index_t j) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR), T>::type _dPsi_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};
    // add other

    /// Implementation of _dPsi(i,j,c,cinv) for NH materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _dPsi(i,j,c,cinv) for MR materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _dPsi(i,j,c,cinv) for Extended NH materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH_ext, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR || _mat==Material::NH_ext), T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    /// Implementation of _d2Psi(i,j) for NH materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    /// Implementation of _d2Psi(i,j) for MR materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR), T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    /// Implementation of _d2Psi(i,j,c,cinv) for NH materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _d2Psi(i,j,c,cinv) for MR materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    /// Implementation of _d2Psi(i,j,c,cinv) for Extended NH materials
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH_ext, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR || _mat==Material::NH_ext), T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};
protected:
    ///////////////////////////////////////////////////////
    // Stretch based formulation                         //
    ///////////////////////////////////////////////////////

    /**
     * @brief      First derivative of a strain energy density function \f$\Psi\f$ w.r.t. the stretch \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     *
     */
    T _dPsi_da   (const index_t a) const;

    /**
     * @brief      Second derivative of a strain energy density function \f$\Psi\f$ w.r.t. the stretches \f$\lambda_a\f$ and \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     *
     */
    T _d2Psi_dab (const index_t a, const index_t b) const;

    /**
     * @brief      First derivative of the volumetric part of a strain energy density function \f$\Psi_{vol}\f$ w.r.t. the stretch \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     */
    T _dPsi_da_vol(const index_t a) const;

    /**
     * @brief      Second derivative of the volumetric part of a strain energy density function \f$\Psi_{vol}\f$ w.r.t. the stretches \f$\lambda_a\f$ and \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     */
    T _d2Psi_dab_vol(const index_t a, const index_t b) const;

    /**
     * @brief      First derivative of the compressibilty function \f$\J\f$ w.r.t. the stretche \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     */
    T _dJ_da     (const index_t a) const;

    /**
     * @brief      First derivative of the compressibilty function \f$\J\f$ w.r.t. the stretches \f$\lambda_a\f$ and \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     */
    T _d2J_dab   (const index_t a, const index_t b) const;

    /**
     * @brief      Lagrange multiplier for incompressible materials
     */
    T _p()                                          const;

    /**
     * @brief      First derivative of the Lagrange multiplier for incompressible materials w.r.t. the stretch \f$\lambda_a\f$
     *
     * @param[in]  a     The index a
     */
    T _dp_da     (const index_t a) const;

    /**
     * @brief     Component \f$a\f$ of the stress
     *
     * @param[in]  a     The index a
     */
    T _Sa        (const index_t a) const;

    /**
     * @brief     First derivative of the \f$a^{\text{th}\f$ component of the stress w.r.t. the stretch \f$\lambda_b\f$
     *
     * @param[in]  a,b     The indices a,b
     */
    T _dSa_db    (const index_t a, const index_t b) const;

    /**
     * @brief     The material matrix for stretch-based implementations
     *
     * @param[in]  a,b,c,d     The indices a,b,c,d
     */
    T _Cabcd     (const index_t a, const index_t b, const index_t c, const index_t d) const;

private:
    // ----------------------------------------------------------------------------------

    /// Specialization of _dPsi_da(a) for compressible NH materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible NH materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::NH), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for compressible MR materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::MR), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible MR materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::MR), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for compressible OG materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::OG), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible OG materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::OG), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for compressible Extended NH materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH_ext), T>::type _dPsi_da_impl(const index_t a) const;

    /// Specialization of _dPsi_da(a) for incompressible Extended NH materials (not implemented)
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::NH_ext), T>::type _dPsi_da_impl(const index_t a) const
    {GISMO_NO_IMPLEMENTATION};

    // other
    template<enum Material _mat, bool _com>
    typename std::enable_if<
                            !(     _mat==Material::NH
                                || _mat==Material::MR
                                || _mat==Material::OG
                                || _mat==Material::NH_ext
                              )
                                                                , T>::type _dPsi_da_impl(const index_t a) const
    {GISMO_NO_IMPLEMENTATION}

    // ----------------------------------------------------------------------------------

    /// Specialization of _d2Psi_dab(a,b) for compressible NH materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible NH materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::NH), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for compressible MR materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::MR), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible MR materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::MR), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for compressible OG materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::OG), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible OG materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::OG), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for compressible Extended NH materials
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH_ext), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    /// Specialization of _d2Psi_dab(a,b) for incompressible Extended NH materials (not implemented)
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::NH_ext), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const
    {GISMO_NO_IMPLEMENTATION};

    // other
    template<enum Material _mat, bool _com>
    typename std::enable_if<
                            !(     _mat==Material::NH
                                || _mat==Material::MR
                                || _mat==Material::OG
                                || _mat==Material::NH_ext
                              )
                                                                , T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const
    {GISMO_NO_IMPLEMENTATION}

    // ----------------------------------------------------------------------------------

    /// Specialization of _Sa(a) for compressible materials
    template<bool _com>
    typename std::enable_if<_com , T>::type _Sa_impl(const index_t a) const;

    /// Specialization of _Sa(a) for incompressible materials
    template<bool _com>
    typename std::enable_if<!_com, T>::type _Sa_impl(const index_t a) const;

    // ----------------------------------------------------------------------------------

    /// Specialization of _dSa_db(a,b) for compressible materials
    template<bool _com>
    typename std::enable_if<_com , T>::type _dSa_db_impl(const index_t a, const index_t b) const;

    /// Specialization of _dSa_db(a,b) for incompressible materials
    template<bool _com>
    typename std::enable_if<!_com, T>::type _dSa_db_impl(const index_t a, const index_t b) const;

    // ----------------------------------------------------------------------------------

    /// Specialization of _Cabcd(a,b,c,d) for compressible materials
    template<bool _com>
    typename std::enable_if<_com , T>::type _Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;

    /// Specialization of _Cabcd(a,b,c,d) for incompressible materials
    template<bool _com>
    typename std::enable_if<!_com, T>::type _Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;

    // ----------------------------------------------------------------------------------

protected:

    /**
     * @brief      Computes the map, the metric quantities and the parameters on
     *             specified points.
     *
     * @param[in]  patch  The patch index
     * @param[in]  u      The in-plane point coordinates
     */
    void _computePoints(const index_t patch, const gsMatrix<T> & u) const;


private:
    /// Implementation of \ref _computePoints for planar geometries
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::OG, void>::type _computePoints_impl(const gsMatrix<T> & u) const;

    /// Implementation of \ref _computePoints for surface geometries
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::OG, void>::type _computePoints_impl(const gsMatrix<T> & u) const;

protected:
    // general
    index_t m_numPars; // how many parameters for the material model?

    // constructor
    using Base::m_patches;
    using Base::m_defpatches;
    const gsFunction<T> * m_thickness;
    std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_density;

    // // Geometric data point
    // mutable gsMapData<T> m_map, m_map_def;
    // mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    // mutable gsMatrix<T,dim,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def;
    // mutable gsMatrix<T,3,2> m_ncov_ori, m_ncov_def;
    // mutable gsMatrix<T,3,3> m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    // mutable gsMatrix<T,3,3> m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    // mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    // mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;
    // mutable T           m_J0_sq;


    // Material parameters and kinematics
    mutable gsMatrix<T> m_parmat;
    mutable gsVector<T> m_parvals;
    mutable gsMatrix<T> m_Tmat,m_rhomat;

    // Geometric data point
    using Base::m_map;
    using Base::m_map_def;

    using Base::m_Acov_ori;
    using Base::m_Acon_ori;
    using Base::m_Acov_def;
    using Base::m_Acon_def;
    using Base::m_Bcov_ori;
    using Base::m_Bcon_ori;
    using Base::m_Bcov_def;
    using Base::m_Bcon_def;
    using Base::m_acov_ori;
    using Base::m_acon_ori;
    using Base::m_acov_def;
    using Base::m_acon_def;
    using Base::m_ncov_ori;
    using Base::m_ncov_def;
    using Base::m_Gcov_ori;
    using Base::m_Gcon_ori;
    using Base::m_Gcov_def;
    using Base::m_Gcon_def;
    using Base::m_Gcov_ori_L;
    using Base::m_Gcov_def_L;
    using Base::m_gcov_ori;
    using Base::m_gcov_def;
    using Base::m_gcon_ori;
    using Base::m_gcon_def;
    using Base::m_Acov_ori_mat;
    using Base::m_Acon_ori_mat;
    using Base::m_Acov_def_mat;
    using Base::m_Acon_def_mat;
    using Base::m_Bcov_ori_mat;
    using Base::m_Bcov_def_mat;
    using Base::m_acov_ori_mat;
    using Base::m_acon_ori_mat;
    using Base::m_acov_def_mat;
    using Base::m_acon_def_mat;
    using Base::m_ncov_ori_mat;
    using Base::m_ncov_def_mat;

    using Base::m_stretches;
    using Base::m_stretchvec;

    using Base::m_J0_sq;
    using Base::m_J_sq;


    gsOptionList m_options;

private:
    static int delta(const int a, const int b)
    {
        return (a==b) ? 1 : 0;
    }

    static int idelta(const int a, const int b)
    {
        return (a!=b) ? 1 : 0;
    }
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrix.hpp)
#endif
