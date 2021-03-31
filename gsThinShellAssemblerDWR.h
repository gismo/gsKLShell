/** @file gsThinShellAssemblerDWR.h

    @brief Provides DWR assembly routines for the Kirchhoff-Love shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsThinShellFunctions.h>

#include <gsPde/gsPointLoads.h>


namespace gismo
{

enum class GoalFunction : short_t
{
    DisplacementNorm    = 0,
    Displacement        = 1,
    MembraneStrain      = 2,
    MembraneStress      = 3,
};

template<class T> class gsThinShellAssemblerDWRBase;


template <short_t d, class T, bool bending, enum GoalFunction GF>
class gsThinShellAssemblerDWR : public gsThinShellAssemblerDWRBase<T>,
                                public gsThinShellAssembler<d,T,bending>
{
public:
    typedef gsThinShellAssembler<d,T,bending> Base;

    gsThinShellAssemblerDWR(
                                const gsMultiPatch<T> & patches,
                                const gsMultiBasis<T> & basisL,
                                const gsMultiBasis<T> & basisH,
                                const gsBoundaryConditions<T> & bconditions,
                                const gsFunction<T> & surface_force,
                                gsMaterialMatrixBase<T> * materialmatrix
                            );

    void assembleMatrixL();
    void assembleMatrixH();

    void assemblePrimalL();
    void assemblePrimalH();

    void assembleDualL() { return assembleDualL_impl<GF>(); }
    void assembleDualH() { return assembleDualH_impl<GF>(); }

    const gsSparseMatrix<T> & matrixL() const { return m_matrixL; }
    const gsSparseMatrix<T> & matrixH() const { return m_matrixH; }

    const gsMatrix<T> primalL() const { return m_pL; }
    const gsMatrix<T> primalH() const { return m_pH; }

    const gsMatrix<T> dualL();
    const gsMatrix<T> dualH();


private:
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::DisplacementNorm, void>::type
    assembleDualL_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::Displacement, void>::type
    assembleDualL_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStrain, void>::type
    assembleDualL_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStress, void>::type
    assembleDualL_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF!=GoalFunction::DisplacementNorm &&
                            _GF!=GoalFunction::Displacement &&
                            _GF!=GoalFunction::MembraneStrain &&
                            _GF!=GoalFunction::MembraneStress
                            , void>::type
    assembleDualL_impl()
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::DisplacementNorm, void>::type
    assembleDualH_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::Displacement, void>::type
    assembleDualH_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStrain, void>::type
    assembleDualH_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStress, void>::type
    assembleDualH_impl();

    template<enum GoalFunction _GF>
    typename std::enable_if<_GF!=GoalFunction::DisplacementNorm &&
                            _GF!=GoalFunction::Displacement &&
                            _GF!=GoalFunction::MembraneStrain &&
                            _GF!=GoalFunction::MembraneStress
                            , void>::type
    assembleDualH_impl()
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    mutable gsThinShellAssemblerBase<T> * m_assemblerL;
    mutable gsThinShellAssemblerBase<T> * m_assemblerH;
    mutable gsMaterialMatrixBase<T> * m_materialMatL;
    mutable gsMaterialMatrixBase<T> * m_materialMatH;

    gsVector<T> m_pL, m_pH, m_dL, m_dH;
    gsSparseMatrix<T> m_matrixL;
    gsSparseMatrix<T> m_matrixH;

    gsMultiBasis<T> m_basisL, m_basisH;

    // using Base::m_patches;
    // using Base::m_basis;
    // using Base::m_bcs;
    // using Base::m_forceFun;
    // using Base::m_materialMat;

    // using Base::m_assembler;

};

/**
 * @brief      Base class for the gsThinShellAssembler
 *
 * @tparam     T     Real type
 *
 * @ingroup    KLShell
 */
template <class T>
class gsThinShellAssemblerDWRBase
{
public:
    /// Default empty constructor
    virtual ~gsThinShellAssemblerDWRBase() {};

    virtual void assembleMatrixL() = 0;
    virtual void assembleMatrixH() = 0;

    virtual void assembleDualL() = 0;
    virtual void assembleDualH() = 0;
};

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssemblerDWRs.hpp)
#endif
