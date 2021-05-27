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
    Displacement        = 0,
    Stretch             = 1,
    MembraneStrain      = 10,
    MembranePStrain     = 11,
    MembraneStress      = 12,
    MembranePStress     = 13,
    MembraneForce       = 14,
    FlexuralStrain      = 20,
    FlexuralStress      = 22,
    FlexuralMoment      = 24,
    Modal               = 100,
    Buckling            = 200,
};

template<class T> class gsThinShellAssemblerDWRBase;


template <short_t d, class T, bool bending>
class gsThinShellAssemblerDWR : public gsThinShellAssemblerDWRBase<T>,
                                public gsThinShellAssembler<d,T,bending>
{
public:
    typedef gsThinShellAssembler<d,T,bending> Base;

    virtual ~gsThinShellAssemblerDWR()
    {
        delete m_assemblerL;
        delete m_assemblerH;
    };

    // empty constructor
    gsThinShellAssemblerDWR() {};


    gsThinShellAssemblerDWR(
                                const gsMultiPatch<T> & patches,
                                const gsMultiBasis<T> & basisL,
                                const gsMultiBasis<T> & basisH,
                                const gsBoundaryConditions<T> & bconditions,
                                const gsFunction<T> & surface_force,
                                gsMaterialMatrixBase<T> * materialmatrix
                            );


    void assembleMassL(bool lumped = false)
    { m_massL = _assembleMass(m_assemblerL, lumped); }
    void assembleMassH(bool lumped = false)
    { m_massH = _assembleMass(m_assemblerH, lumped); }
    void assembleMatrixL()
    { m_matrixL = _assembleMatrix(m_assemblerL); }
    void assembleMatrixH()
    { m_matrixH = _assembleMatrix(m_assemblerH); }
    void assembleMatrixL(const gsMultiPatch<T> & deformed)
    { m_matrixL = _assembleMatrix(m_assemblerL,deformed); }
    void assembleMatrixH(const gsMultiPatch<T> & deformed)
    { m_matrixH = _assembleMatrix(m_assemblerH,deformed); }

    void assemblePrimalL()
    { m_pL = _assemblePrimal(m_assemblerL); }
    void assemblePrimalH()
    { m_pH = _assemblePrimal(m_assemblerH); }
    void assemblePrimalL(const gsMultiPatch<T> & deformed)
    { m_pL = _assemblePrimal(m_assemblerL,deformed); }
    void assemblePrimalH(const gsMultiPatch<T> & deformed)
    { m_pH = _assemblePrimal(m_assemblerH,deformed); }

    void assembleDualL(const gsMultiPatch<T> & primal)
    { m_dL += _assembleDual(m_assemblerL,primal); }
    void assembleDualH(const gsMultiPatch<T> & primal)
    { m_dH += _assembleDual(m_assemblerH,primal); }

    void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { m_dL += _assembleDual(points,m_assemblerL,primal); }
    void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { m_dH += _assembleDual(points,m_assemblerH,primal); }

    const gsSparseMatrix<T> & matrixL() const { return m_matrixL; }
    const gsSparseMatrix<T> & matrixH() const { return m_matrixH; }

    const gsSparseMatrix<T> & massL() const { return m_massL; }
    const gsSparseMatrix<T> & massH() const { return m_massH; }

    const gsMatrix<T> primalL() const { return m_pL; }
    const gsMatrix<T> primalH() const { return m_pH; }

    const gsMatrix<T> dualL() const { return m_dL; }
    const gsMatrix<T> dualH() const { return m_dH; }

    void constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result);
    void constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result);

    void constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed);
    void constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed);
    gsMultiPatch<T> constructSolutionL(const gsMatrix<T> & solVector);
    gsMultiPatch<T> constructSolutionH(const gsMatrix<T> & solVector);
    gsMultiPatch<T> constructDisplacementL(const gsMatrix<T> & solVector);
    gsMultiPatch<T> constructDisplacementH(const gsMatrix<T> & solVector);

    // Linear elasticity
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH);

    // Nonlinear elasticity
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);

    // Eigenvalues
    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal);

    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);

    T computeGoal(const gsMultiPatch<T> & deformed);
    T computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH) const;
    T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, const gsMultiPatch<T> &deformed) const;

    void setGoal(enum GoalFunction GF, short_t component = 9) { m_goalFunction = GF; m_component = component; }

protected:

    void _setBasis(const gsMultiBasis<T> & basis);

    gsSparseMatrix<T>   _assembleMass(gsThinShellAssemblerBase<T> * assembler, bool lumped = false);

    gsSparseMatrix<T>   _assembleMatrix(gsThinShellAssemblerBase<T> * assembler);
    gsSparseMatrix<T>   _assembleMatrix(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed);
    gsVector<T>         _assemblePrimal(gsThinShellAssemblerBase<T> * assembler);
    gsVector<T>         _assemblePrimal(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed);

    /**
     * @brief      Assembles the dual as a domain integral
     *
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed);

    /**
     * @brief      Assembles the dual as a domain integral
     *
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed);

    template<int _d, bool _bending>
    typename std::enable_if<(_d==3 && _bending), T>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);


    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), T>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);


protected:
    mutable gsThinShellAssemblerBase<T> * m_assemblerL;
    mutable gsThinShellAssemblerBase<T> * m_assemblerH;
    mutable gsMaterialMatrixBase<T> * m_materialMatL;
    mutable gsMaterialMatrixBase<T> * m_materialMatH;

    gsVector<T> m_pL, m_pH, m_dL, m_dH;
    gsSparseMatrix<T> m_matrixL;
    gsSparseMatrix<T> m_matrixH;
    gsSparseMatrix<T> m_massL;
    gsSparseMatrix<T> m_massH;

    gsMultiBasis<T> m_basisL, m_basisH;

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // using Base::m_patches;
    // using Base::m_basis;

    using Base::m_patches;
    using Base::m_defpatches;
    using Base::m_materialMat;
    using Base::m_forceFun;
    using Base::m_options;
    using Base::m_foundFun;
    using Base::m_foundInd;
    using Base::m_pressFun;
    using Base::m_pressInd;

    enum GoalFunction m_goalFunction;
    short_t m_component;
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

    virtual void assembleMassL(bool lumped = false) =0;
    virtual void assembleMassH(bool lumped = false) =0;
    virtual void assembleMatrixL() =0;
    virtual void assembleMatrixH() =0;
    virtual void assembleMatrixL(const gsMultiPatch<T> & deformed) =0;
    virtual void assembleMatrixH(const gsMultiPatch<T> & deformed) =0;

    virtual void assemblePrimalL() =0;
    virtual void assemblePrimalH() =0;
    virtual void assemblePrimalL(const gsMultiPatch<T> & deformed) =0;
    virtual void assemblePrimalH(const gsMultiPatch<T> & deformed) =0;


    virtual void assembleDualL(const gsMultiPatch<T> & primal) =0;
    virtual void assembleDualH(const gsMultiPatch<T> & primal) =0;

    virtual void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal) =0;
    virtual void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal) =0;

    virtual const gsSparseMatrix<T> & matrixL() const =0;
    virtual const gsSparseMatrix<T> & matrixH() const =0;

    virtual const gsSparseMatrix<T> & massL() const =0;
    virtual const gsSparseMatrix<T> & massH() const =0;

    virtual const gsMatrix<T> primalL() const =0;
    virtual const gsMatrix<T> primalH() const =0;

    virtual const gsMatrix<T> dualL() const =0;
    virtual const gsMatrix<T> dualH() const =0;

    virtual void constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) =0;
    virtual void constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) =0;

    virtual void constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) =0;
    virtual void constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) =0;
    virtual gsMultiPatch<T> constructSolutionL(const gsMatrix<T> & solVector) =0;
    virtual gsMultiPatch<T> constructSolutionH(const gsMatrix<T> & solVector) =0;
    virtual gsMultiPatch<T> constructDisplacementL(const gsMatrix<T> & solVector) =0;
    virtual gsMultiPatch<T> constructDisplacementH(const gsMatrix<T> & solVector) =0;

    // Linear elasticity ;
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH) =0;

    // Nonlinear elasticity
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed) =0;

    // Eigenvalues
    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal) =0;
    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

    virtual T computeGoal(const gsMultiPatch<T> & deformed) =0;

    virtual T computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed) =0;

    virtual T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH) const = 0;
    virtual T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, const gsMultiPatch<T> &deformed) const = 0;

    virtual void setGoal(enum GoalFunction GF, short_t component = 9) = 0;
};

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssemblerDWR.hpp)
#endif
