/** @file gsMaterialMatrixComposite.h

    @brief Provides a material matrix for laminates

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
 * @brief      This class defines a linear material laminate
 *
 * @tparam     dim   The dimension of the problem (2 = planar, 3 = surface)
 * @tparam     T     Real type
 *
 * @ingroup    MaterialMatrix
 *
 */
template <  short_t dim,
            class T
         >
class gsMaterialMatrixComposite : public gsMaterialMatrixBaseDim<dim,T>
{
public:
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    /// Default empty constructor
    gsMaterialMatrixComposite() { }

    gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const std::vector< gsFunctionSet<T> *>              & thickness,
                            const std::vector< gsFunctionSet<T> *>              & G,
                            const std::vector< gsFunctionSet<T> *>              & alpha,
                            const std::vector< gsFunctionSet<T> *>              & rho           );

    gsMaterialMatrixComposite(
                            const gsFunctionSet<T>                              & mp,
                            const std::vector< gsFunctionSet<T> *>              & thickness,
                            const std::vector< gsFunctionSet<T> *>              & G,
                            const std::vector< gsFunctionSet<T> *>              & alpha         );


    enum MatIntegration isMatIntegrated() const {return MatIntegration::Integrated; }
    enum MatIntegration isVecIntegrated() const {return MatIntegration::Integrated; }

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    // template COM
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COM
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// See \ref gsMaterialMatrixBase for details
    void transform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {GISMO_NO_IMPLEMENTATION;}

    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const
    { GISMO_NO_IMPLEMENTATION; }

    void setParameters(const std::vector<gsFunction<T>*> &pars)
    { GISMO_NO_IMPLEMENTATION; }
    void info() const;

public:
    /// Shared pointer for gsMaterialMatrixComposite
    typedef memory::shared_ptr< gsMaterialMatrixComposite > Ptr;

    /// Unique pointer for gsMaterialMatrixComposite
    typedef memory::unique_ptr< gsMaterialMatrixComposite > uPtr;

protected:
    void _initialize();
    void _initializeParameters();
    void _defaultOptions();

    gsMatrix<T> _transformationMatrix(const gsMatrix<T> & phi, const gsMatrix<T> & u) const;
    gsMatrix<T> _cart2cov(const gsVector<T> a1, const gsVector<T> a2, const gsVector<T> e1, const gsVector<T> e2) const;
    gsMatrix<T> _con2cart(const gsVector<T> ac1, const gsVector<T> ac2, const gsVector<T> e1, const gsVector<T> e2) const;


protected:

    // template MAT
    void _computePoints(const index_t patch, const gsMatrix<T> & u, bool basis = true) const;

protected:
    // constructor
    using Base::m_patches;
    using Base::m_defpatches;
    const gsFunction<T> * m_thickness;
    const gsFunction<T> * m_E11;
    const gsFunction<T> * m_E22;
    const gsFunction<T> * m_G12;
    const gsFunction<T> * m_nu12;
    const gsFunction<T> * m_nu21;
    const gsFunction<T> * m_phi;
    std::vector<gsFunction<T>* > m_pars; // TO DO: change to uPtr
    const gsFunction<T> * m_rho;

    // Composite
    index_t m_nLayers;
    mutable gsMatrix<T>                 m_Tmat, m_E1mat, m_E2mat, m_G12mat, m_nu12mat, m_nu21mat, m_phiMat, m_rhoMat;

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

    const std::vector<gsFunctionSet<T> * > m_Ts;
    const std::vector<gsFunctionSet<T> * > m_Gs;
    const std::vector<gsFunctionSet<T> * > m_As;
    const std::vector<gsFunctionSet<T> * > m_Rs;
    mutable std::vector< gsMatrix<T> > m_Gcontainer, m_Tcontainer, m_Acontainer, m_Rcontainer;
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrixComposite.hpp)
#endif
