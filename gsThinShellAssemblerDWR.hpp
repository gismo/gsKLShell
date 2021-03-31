/** @file gsThinShellAssemblerDWR.hpp

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

#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixBase.h>
#include <gsKLShell/gsMaterialMatrixEval.h>

namespace gismo
{

template <short_t d, class T, bool bending, enum GoalFunction GF>
gsThinShellAssemblerDWR<d, T, bending, GF>::gsThinShellAssemblerDWR(
                                                                const gsMultiPatch<T> & patches,
                                                                const gsMultiBasis<T> & basisL,
                                                                const gsMultiBasis<T> & basisH,
                                                                const gsBoundaryConditions<T> & bconditions,
                                                                const gsFunction<T> & surface_force,
                                                                gsMaterialMatrixBase<T> * materialmatrix
                                                            )
                                                            :
                                                            Base(patches,basisL,bconditions,surface_force,materialmatrix),
                                                            m_basisL(basisL),
                                                            m_basisH(basisH)
{
    gsDebugVar(Base::m_patches);
    gsDebugVar(Base::m_basis.basis(0).size());
    Base::m_basis = m_basisH;
    gsDebugVar(Base::m_basis.basis(0).size());
}


template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::assembleMatrixL()
{
    Base::m_basis = m_basisL;
    Base::_initialize();
    Base::assemble();
    m_matrixL = Base::matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::assembleMatrixH()
{
    Base::m_basis = m_basisH;
    Base::_initialize();
    Base::assemble();
    m_matrixH = Base::matrix();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::assemblePrimalL()
{
    Base::m_basis = m_basisL;
    Base::_initialize();
    Base::assemble();
    m_pL = Base::rhs();
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
void gsThinShellAssemblerDWR<d, T, bending, GF>::assemblePrimalH()
{
    Base::m_basis = m_basisH;
    Base::_initialize();
    Base::assemble();
    m_pH = Base::rhs();
}


template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::DisplacementNorm, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualL_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::Displacement, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualL_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStrain, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualL_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStress, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualL_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::DisplacementNorm, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualH_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::Displacement, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualH_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStrain, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualH_impl()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, class T, bool bending, enum GoalFunction GF>
template<enum GoalFunction _GF>
typename std::enable_if<_GF==GoalFunction::MembraneStress, void>::type
gsThinShellAssemblerDWR<d, T, bending, GF>::assembleDualH_impl()
{
    GISMO_NO_IMPLEMENTATION;
}



}// namespace gismo
