/** @file gsThinShellFunctions.h

    @brief Provides evaluation function for stresses.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)

*/

#pragma once

#include <gsThinShell2/gsMaterialMatrix.h>

namespace gismo
{

/** @brief Specifies the type of stresses to compute.
 *
 *         Currently, gsWriteParaview can only plot vector-valued functions with an output dimension up to three.
 *         Therefore it not possible to plot all stress components as components of a single vector-valued function.
*/
struct stress_type
{
    enum type
    {
        von_mises          = 0,  /// compute only von Mises stress
        von_mises_membrane = 1,  /// compute only von Mises stress - membrane stresses
        von_mises_flexural = 2,  /// compute only von Mises stress - flexural stresses
        membrane           = 3,  /// compute normal and shear stresses due to membrane component
        flexural           = 4,  /// compute normal and shear stresses due to membrane component
        total              = 5,  /// compute normal and shear stresses due to both components
        membrane_strain    = 6,  /// compute normal and shear stresses due to both components
        flexural_strain    = 7,  /// compute normal and shear stresses due to both components
        principal_stretch  = 8,  /// principal stretches
    };
};

/** @brief Compute Cauchy stresses for a previously computed/defined displacement field.
 *         Can be pushed into gsPiecewiseFunction to construct gsField for visualization in Paraview.
*/
template <class T>
class gsShellStressFunction : public gsFunction<T>
{
public:

    gsShellStressFunction(const gsMultiPatch<T> & geometry,
                           const gsMultiPatch<T> & deformed,
                           const gsMaterialMatrix<T> & mm,
                           index_t patch,
                           stress_type::type type,
                           const gsExprAssembler<T> & assembler
                           )
        : m_patches(geometry),
          m_defpatches(deformed),
          m_materialMat(mm),
          m_patchID(patch), ///WHAT DO WE DO WITH THIS?
          m_stress_type(type),
          m_assembler(assembler)
    {}

    virtual short_t domainDim() const
    {
        return 2;
    }

    virtual short_t targetDim() const
    {
        return 3;
    }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Each column of the output matrix (result) corresponds to a set of stress components:
     *         all_2D:    s_11 s_22 s_12
     *         normal_3D: s_11 s_22 s_33
     *         shear_3D:  s_12 s_13 s_23
     *         or is one number in case if stress_type::von_mises is chosen.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    const gsMultiPatch<T>& m_patches;
    const gsMultiPatch<T>& m_defpatches;
    index_t m_patchID;
    stress_type::type m_stress_type;
    gsMaterialMatrix<T> m_materialMat;
    mutable gsExprAssembler<> m_assembler;

}; // class definition ends


} // namespace ends


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellFunctions.hpp)
#endif
