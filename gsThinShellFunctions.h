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

#include <gsKLShell/gsMaterialMatrix.h>

namespace gismo
{

/** @brief Specifies the type of stresses to compute.
 *
 *         Currently, gsWriteParaview can only plot vector-valued functions with an output dimension up to three.
 *         Therefore it not possible to plot all stress components as components of a single vector-valued function.
 *
 *  \ingroup KLShell
 *
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
        principal_stress_membrane  = 9,  /// principal stress membrane
        principal_stress_flexural  = 10,  /// principal stress bending
        principal_stretch_dir1  = 11,  /// principal stretch directions
        principal_stretch_dir2  = 12,  /// principal stretch directions
        principal_stretch_dir3  = 13,  /// principal stretch directions
    };
};

/** @brief Compute Cauchy stresses for a previously computed/defined displacement field.
 *         Can be pushed into gsPiecewiseFunction to construct gsField for visualization in Paraview.
 *
 *  \ingroup KLShell
 *
*/
template <class T>
class gsShellStressFunction : public gsFunction<T>
{
public:

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  geometry   The undeformed geometry
     * @param[in]  deformed   The deformed geometry
     * @param      mm         The material matrix
     * @param[in]  patch      The patch index
     * @param[in]  type       The stress type
     * @param[in]  assembler  The shell assembler
     */
    gsShellStressFunction(const gsMultiPatch<T> & geometry,
                           const gsMultiPatch<T> & deformed,
                           gsMaterialMatrixBase<T> * mm,
                           index_t patch,
                           stress_type::type type,
                           const gsExprAssembler<T> & assembler
                           )
        : m_patches(geometry),
          m_defpatches(deformed),
          m_materialMat(mm),
          m_patchID(patch),
          m_stress_type(type),
          m_assembler(assembler)
    {}

    virtual short_t domainDim() const
    {
        return 2;
    }

    virtual short_t targetDim() const
    {
        switch (m_stress_type)
        {
            default:
                return 0;
                break;
            case stress_type::membrane :
                return 3;
                break;

            case stress_type::flexural :
                return 3;
                break;

            // TO BE IMPLEMENTED
            // -------------------------------------
            case stress_type::von_mises :
                return 1;
                break;

            case stress_type::von_mises_membrane :
                return 1;
                break;

            case stress_type::von_mises_flexural :
                return 1;
                break;

            case stress_type::total :
                return 1;
                break;
            // -------------------------------------

            case stress_type::membrane_strain :
                return 3;
                break;

            case stress_type::flexural_strain :
                return 3;
                break;
            case stress_type::principal_stretch :
                return 3;
                break;
            case stress_type::principal_stress_membrane :
                return 2;
                break;
            case stress_type::principal_stress_flexural :
                return 2;
                break;
            case stress_type::principal_stretch_dir1 :
                return 3;
                break;
            case stress_type::principal_stretch_dir2 :
                return 3;
                break;
            case stress_type::principal_stretch_dir3 :
                return 3;
                break;
            /*
                DEFAULT includes:
                stress_type::membrane;
                stress_type::flexural;
                stress_type::membrane_strain;
                stress_type::flexural_strain;
                stress_type::principal_stretch;
                stress_type::principal_stretch_dir1;
                stress_type::principal_stretch_dir2;
                stress_type::principal_stretch_dir3
            */
        }
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
    mutable gsMaterialMatrixBase<T> * m_materialMat;
    index_t m_patchID;
    stress_type::type m_stress_type;
    mutable gsExprAssembler<> m_assembler;

}; // class definition ends


} // namespace ends


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellFunctions.hpp)
#endif
