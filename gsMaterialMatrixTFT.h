/** @file gsMaterialMatrixTFT.h

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

#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixBaseDim.h>
#include <gsKLShell/gsMaterialMatrixUtils.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
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
            class T
         >
class gsMaterialMatrixTFT : public gsMaterialMatrixBaseDim<dim,T>
{
public:
    using Base = gsMaterialMatrixBaseDim<dim,T>;

    /**
     * @brief      Constructor without deformed multipatch and density
     *
     * @param[in]  mp         Original geometry
     * @param[in]  thickness  Thickness function
     * @param[in]  pars       Vector with parameters (E, nu)
     */
    gsMaterialMatrixTFT(gsMaterialMatrixLinear<dim,T> * materialmatrix)
    :
    m_materialMat(materialmatrix)
    {
        // if (!dynamic_cast<const gsMaterialMatrixLinear<dim,T> *>(m_materialMat))
            // GISMO_ERROR("Material matrix must be linear");
    }

    /// Destructor
    gsMaterialMatrixTFT() { }

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isMatIntegrated() const {return m_materialMat->isMatIntegrated(); }

    /// See \ref gsMaterialMatrixBase for details
    /// Here, we use the integrated matrix and vector to make our computations!
    inline enum MatIntegration isVecIntegrated() const {return m_materialMat->isVecIntegrated(); }

    /// See \ref gsMaterialMatrixBase for details
    gsOptionList & options()
    {return m_materialMat->options(); }

    /// See \ref gsMaterialMatrixBase for details
    void setOptions(gsOptionList opt)
    {m_materialMat->setOptions(opt); }

    /// See \ref gsMaterialMatrixBase for details
    void density_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->density_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void stretch_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->stretch_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void stretchDir_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { m_materialMat->stretchDir_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void thickness_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->thickness_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void covtransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->covtransform_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    void contransform_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result) const
    { m_materialMat->contransform_into( patch,u,result ); }

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_matrix(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_vector(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const;

    /// See \ref gsMaterialMatrixBase for details
    gsMatrix<T> eval3D_pstress(const index_t patch, const gsMatrix<T> & u, const gsMatrix<T>& z, enum MaterialOutput out = MaterialOutput::Generic) const
    { return m_materialMat->eval3D_pstress( patch,u,z,out); }

    /// See \ref gsMaterialMatrixBase for details
    void setParameters(const std::vector<gsFunction<T>*> &pars)
    { m_materialMat->setParameters(pars); }

    /// See \ref gsMaterialMatrixBase for details
    void info() const
    { m_materialMat->info(); }

    /// See \ref gsMaterialMatrixBase for details
    void setDeformed(const gsFunctionSet<T> * deformed)
    {
        gsDebugVar(deformed);
        m_materialMat->setDeformed(deformed);
    }


public:
    /// Shared pointer for gsMaterialMatrixTFT
    typedef memory::shared_ptr< gsMaterialMatrixTFT > Ptr;

    /// Unique pointer for gsMaterialMatrixTFT
    typedef memory::unique_ptr< gsMaterialMatrixTFT > uPtr;

protected:
    mutable gsMaterialMatrixLinear<dim,T> * m_materialMat;

    // // general
    // index_t m_numPars; // how many parameters for the material model?

    // // constructor
    // using Base::m_patches;
    // using Base::m_defpatches;
    // const gsFunction<T> * m_thickness;
    // std::vector<gsFunction<T>* > m_pars;
    // const gsFunction<T> * m_density;

    // mutable gsMatrix<T> m_Emat,m_Nmat,m_Tmat,m_rhomat;
    // mutable real_t m_lambda, m_mu;

    // mutable gsMatrix<T>                 m_parmat;
    // mutable gsVector<T>                 m_parvals;

    // mutable gsMatrix<T> m_pstress, m_pstressvec;

    // // Geometric data point
    // using Base::m_map;
    // using Base::m_map_def;

    // using Base::m_Acov_ori;
    // using Base::m_Acon_ori;
    // using Base::m_Acov_def;
    // using Base::m_Acon_def;
    // using Base::m_Bcov_ori;
    // using Base::m_Bcon_ori;
    // using Base::m_Bcov_def;
    // using Base::m_Bcon_def;
    // using Base::m_acov_ori;
    // using Base::m_acon_ori;
    // using Base::m_acov_def;
    // using Base::m_acon_def;
    // using Base::m_ncov_ori;
    // using Base::m_ncov_def;
    // using Base::m_Gcov_ori;
    // using Base::m_Gcon_ori;
    // using Base::m_Gcov_def;
    // using Base::m_Gcon_def;
    // using Base::m_Gcov_ori_L;
    // using Base::m_Gcov_def_L;
    // using Base::m_gcov_ori;
    // using Base::m_gcov_def;
    // using Base::m_gcon_ori;
    // using Base::m_gcon_def;
    // using Base::m_Acov_ori_mat;
    // using Base::m_Acon_ori_mat;
    // using Base::m_Acov_def_mat;
    // using Base::m_Acon_def_mat;
    // using Base::m_Bcov_ori_mat;
    // using Base::m_Bcov_def_mat;
    // using Base::m_acov_ori_mat;
    // using Base::m_acon_ori_mat;
    // using Base::m_acov_def_mat;
    // using Base::m_acon_def_mat;
    // using Base::m_ncov_ori_mat;
    // using Base::m_ncov_def_mat;

    // using Base::m_stretches;
    // using Base::m_stretchvec;

    // using Base::m_J0_sq;
    // using Base::m_J_sq;

};

template <class T>
class gsScalarRootFinder
{
protected:
    T m_x;
    T m_xmin;
    T m_xmax;
    T m_error;
    T m_tolerance;

    index_t m_maxIterations;
    index_t m_iteration;

    typedef std::function < T ( T const &) > function_t;
    function_t m_function;

public:
    gsScalarRootFinder(T xmin, T xmax, function_t fun)
    :
    m_xmin(xmin),
    m_xmax(xmax),
    m_function(fun),
    m_iteration(0),
    m_maxIterations(15),
    m_error(1.0),
    m_tolerance(1e-15)
    {
        gsDebugVar(m_function(m_xmin));
        gsDebugVar(m_function(m_xmax));
    }

    void compute()
    {
        // From: https://lemesurierb.github.io/elementary-numerical-analysis-python/notebooks/root-finding-without-derivatives-python.html
        T x_older = m_xmin;
        T x_more_recent = m_xmax;
        T x_new;
        T f_x_older = m_function(x_older);
        T f_x_more_recent = m_function(x_more_recent);
        T f_x_new;
        for (m_iteration = 0; m_iteration <= m_maxIterations; m_iteration++)
        {
            gsDebugVar(m_iteration);
            x_new = (x_older * f_x_more_recent - f_x_older * x_more_recent)/(f_x_more_recent - f_x_older);
            f_x_new = m_function(x_new);

            x_older = x_more_recent;
            x_more_recent = x_new;
            f_x_older = f_x_more_recent;
            f_x_more_recent = f_x_new;

            m_error = std::abs(x_older - x_more_recent);

            gsDebugVar(m_error);

            if (m_error < m_tolerance)
                break;
        }
        if (m_iteration < m_maxIterations)
        {
            gsDebugVar("converged");
            m_x = x_new;
        }
        else
            gsDebugVar("not converged");

    }

    T result() { return m_x; }
    T error()  { return m_error; }
};

} // namespace


// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsMaterialMatrixTFT.hpp)
// #endif
