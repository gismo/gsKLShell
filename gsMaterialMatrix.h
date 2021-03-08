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

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <  short_t dim,
            class T,
            short_t matId,
            bool comp,
            enum Material mat = decodeMat_id<matId>::material,
            enum Implementation imp  = decodeMat_id<matId>::implementation
         >
class gsMaterialMatrix : public gsMaterialMatrixBase<T>
{
public:
    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    /// Default empty constructor
    gsMaterialMatrix() { }

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars,
                        const gsFunction<T> & Density);

    inline enum MatIntegration isMatIntegrated() const {return MatIntegration::NotIntegrated; }
    inline enum MatIntegration isVecIntegrated() const {return MatIntegration::NotIntegrated; }

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt

    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
        m_numPars = m_pars.size();
    }

    void info() const;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

protected:
    void _initialize();
    void _defaultOptions();

private:
    template<bool _com>
    typename std::enable_if<_com, void>::type _stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<bool _com>
    typename std::enable_if<!_com, void>::type _stretch_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    template<bool _com>
    typename std::enable_if<_com, void>::type _stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<bool _com>
    typename std::enable_if<!_com, void>::type _stretchDir_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

protected:
    gsMatrix<T> _eval_Incompressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval_Incompressible_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval_Incompressible_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    gsMatrix<T> _eval_Compressible_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval_Compressible_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    gsMatrix<T> _eval_Compressible_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

private:

    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_matrix_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_vector_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_mat==Material::SvK, gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && !(_mat==Material::SvK), gsMatrix<T>>::type _eval3D_pstress_impl(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z) const;


protected:

    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    T _Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _Sij    (const index_t i, const index_t j) const;
    T _Sij    (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _Sii    (const index_t i) const;
    T _Sii    (const index_t i, const gsMatrix<T> & c) const;

private:

    // Incompressible
    // template<enum Material _mat>
    // typename std::enable_if<_mat==0, T>::type _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH  && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR  && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

    // Compressible
    template<enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral,   T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Implementation _imp>
    typename std::enable_if<!(_imp==Implementation::Spectral),T>::type
    _Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Sij_impl(const index_t i, const index_t j) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Sij_impl(const index_t i, const index_t j) const;

    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::SvK && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::MR && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::OG && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_mat==Material::NH_ext && _imp==Implementation::Analytical, T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Spectral , T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat, enum Implementation _imp>
    typename std::enable_if<_imp==Implementation::Generalized , T>::type
    _Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

protected:

    T _dPsi   (const index_t i, const index_t j) const;
    T _dPsi   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _dPsi_vol(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    T _d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _d2Psi_vol(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    T _dI_1   (const index_t i, const index_t j) const;
    T _dI_2   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

private:
    /*
        Note: when adding a new implementation, make sure to add the condition in the 'other' implementation : !(_mat==xx || ...)
    */
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _dPsi_impl(const index_t i, const index_t j) const;
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _dPsi_impl(const index_t i, const index_t j) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR), T>::type _dPsi_impl(const index_t i, const index_t j) const
    {GISMO_NO_IMPLEMENTATION};
    // add other

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH_ext, T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR || _mat==Material::NH_ext), T>::type _dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR), T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
    {GISMO_NO_IMPLEMENTATION};

    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::MR, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::NH_ext, T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // other
    template<enum Material _mat>
    typename std::enable_if<!(_mat==Material::NH || _mat==Material::MR || _mat==Material::NH_ext), T>::type _d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
    {GISMO_NO_IMPLEMENTATION};
protected:

    // Stretch based formulation
    T _dPsi_da   (const index_t a) const;
    T _d2Psi_dab (const index_t a, const index_t b) const;
    T _dPsi_da_vol(const index_t a) const;
    T _d2Psi_dab_vol(const index_t a, const index_t b) const;
    T _dJ_da     (const index_t a) const;
    T _d2J_dab   (const index_t a, const index_t b) const;
    T _p()                                          const;
    T _dp_da     (const index_t a) const;
    T _Sa        (const index_t a) const;
    T _dSa_db    (const index_t a, const index_t b) const;
    T _Cabcd     (const index_t a, const index_t b, const index_t c, const index_t d) const;

private:
    // ----------------------------------------------------------------------------------

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH), T>::type _dPsi_da_impl(const index_t a) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::NH), T>::type _dPsi_da_impl(const index_t a) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::MR), T>::type _dPsi_da_impl(const index_t a) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::MR), T>::type _dPsi_da_impl(const index_t a) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::OG), T>::type _dPsi_da_impl(const index_t a) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::OG), T>::type _dPsi_da_impl(const index_t a) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH_ext), T>::type _dPsi_da_impl(const index_t a) const;
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

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::NH), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::MR), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::MR), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::OG), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;
    template<enum Material _mat, bool _com>
    typename std::enable_if<!_com && (_mat==Material::OG), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;

    template<enum Material _mat, bool _com>
    typename std::enable_if<_com && (_mat==Material::NH_ext), T>::type _d2Psi_dab_impl(const index_t a, const index_t b) const;
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

    template<bool _com>
    typename std::enable_if<_com , T>::type _Sa_impl(const index_t a) const;
    template<bool _com>
    typename std::enable_if<!_com, T>::type _Sa_impl(const index_t a) const;

    // ----------------------------------------------------------------------------------

    template<bool _com>
    typename std::enable_if<_com , T>::type _dSa_db_impl(const index_t a, const index_t b) const;
    template<bool _com>
    typename std::enable_if<!_com, T>::type _dSa_db_impl(const index_t a, const index_t b) const;

    // ----------------------------------------------------------------------------------

    template<bool _com>
    typename std::enable_if<_com , T>::type _Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;
    template<bool _com>
    typename std::enable_if<!_com, T>::type _Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;

    // ----------------------------------------------------------------------------------

protected:

    void _computePoints(const index_t patch, const gsMatrix<T> & u) const;
    void _computeMetricDeformed() const;
    void _computeMetricUndeformed() const;
    void _getMetric(index_t k, T z) const;
    void _getMetricDeformed(index_t k, T z) const;
    void _getMetricUndeformed(index_t k, T z) const;
    void _computeStretch(const gsMatrix<T> & C ) const;
    std::pair<gsVector<T>,gsMatrix<T>> _evalStretch(const gsMatrix<T> & C ) const;

private:
    template<enum Material _mat>
    typename std::enable_if<_mat==Material::OG, void>::type _computePoints_impl(const gsMatrix<T> & u) const;
    template<enum Material _mat>
    typename std::enable_if<_mat!=Material::OG, void>::type _computePoints_impl(const gsMatrix<T> & u) const;

    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricDeformed_impl() const;
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricDeformed_impl() const;

    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _computeMetricUndeformed_impl() const;
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _computeMetricUndeformed_impl() const;

    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetric_impl(index_t k, T z) const;
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetric_impl(index_t k, T z) const;

    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricDeformed_impl(index_t k, T z) const;
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricDeformed_impl(index_t k, T z) const;

    template<short_t _dim>
    typename std::enable_if<_dim==2, void>::type _getMetricUndeformed_impl(index_t k, T z) const;
    template<short_t _dim>
    typename std::enable_if<_dim==3, void>::type _getMetricUndeformed_impl(index_t k, T z) const;

    // template<short_t _dim>
    // typename std::enable_if<_dim==2, void>::type _computeStretch_impl(const gsMatrix<T> & C ) const;
    // template<short_t _dim>
    // typename std::enable_if<_dim==3, void>::type _computeStretch_impl(const gsMatrix<T> & C ) const;

    // template<short_t _dim>
    // typename std::enable_if<_dim==2, std::pair<gsVector<T>,gsMatrix<T>>>::type _evalStretch_impl(const gsMatrix<T> & C ) const;
    // template<short_t _dim>
    // typename std::enable_if<_dim==3, std::pair<gsVector<T>,gsMatrix<T>>>::type _evalStretch_impl(const gsMatrix<T> & C ) const;


protected:
    // general
    index_t m_numPars; // how many parameters for the material model?

    // constructor
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
    const gsFunction<T> * m_thickness;
    std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_density;

    // Geometric data point
    mutable gsMapData<T> m_map, m_map_def;
    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,dim,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def;
    mutable gsMatrix<T,3,2> m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T,3,3> m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3> m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;


    // Material parameters and kinematics
    mutable gsMatrix<T> m_parmat;
    mutable gsVector<T> m_parvals;
    mutable gsMatrix<T> m_Tmat,m_rhomat;
    mutable T           m_J0, m_J0_sq, m_J, m_J_sq;
    mutable gsMatrix<T> m_stretches, m_stretchvec;

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
