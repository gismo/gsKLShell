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
    MembranePStress     = 4,
    MembraneForce       = 5,
    Modal               = 10,
    Buckling            = 20,
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


    void assembleMassL(bool lumped = false)
    { m_massL = _assembleMass(m_basisL, lumped); }
    void assembleMassH(bool lumped = false)
    { m_massH = _assembleMass(m_basisH, lumped); }
    void assembleMatrixL()
    { m_matrixL = _assembleMatrix(m_basisL); }
    void assembleMatrixH()
    { m_matrixH = _assembleMatrix(m_basisH); }
    void assembleMatrixL(const gsMultiPatch<T> & deformed)
    { m_matrixL = _assembleMatrix(m_basisL,deformed); }
    void assembleMatrixH(const gsMultiPatch<T> & deformed)
    { m_matrixH = _assembleMatrix(m_basisH,deformed); }

    void assemblePrimalL()
    { m_pL = _assemblePrimal(m_basisL); }
    void assemblePrimalH()
    { m_pH = _assemblePrimal(m_basisH); }
    void assemblePrimalL(const gsMultiPatch<T> & deformed)
    { m_pL = _assemblePrimal(m_basisL,deformed); }
    void assemblePrimalH(const gsMultiPatch<T> & deformed)
    { m_pH = _assemblePrimal(m_basisH,deformed); }


    void assembleDualL(const gsMultiPatch<T> & primal)
    { m_dL += _assembleDual(m_basisL,primal); }
    void assembleDualH(const gsMultiPatch<T> & primal)
    { m_dH += _assembleDual(m_basisH,primal); }

    void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { m_dL += _assembleDual(points,m_basisL,primal); }
    void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { m_dH += _assembleDual(points,m_basisH,primal); }

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
                      const gsMultiPatch<T> & deformed);
    T computeGoalEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & deformed)
    {GISMO_NO_IMPLEMENTATION;}

    T computeGoal(const gsMultiPatch<T> & deformed)
    { return computeGoal_impl<GF>(deformed); }

    T computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
    { return computeGoal_impl<GF>(points, deformed); }

protected:

    void _setBasis(const gsMultiBasis<T> & basis);

    gsSparseMatrix<T>   _assembleMass(const gsMultiBasis<T> & basis, bool lumped = false);

    gsSparseMatrix<T>   _assembleMatrix(const gsMultiBasis<T> & basis);
    gsSparseMatrix<T>   _assembleMatrix(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & deformed);
    gsVector<T>         _assemblePrimal(const gsMultiBasis<T> & basis);
    gsVector<T>         _assemblePrimal(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & deformed);

    /**
     * @brief      Assembles the dual as a domain integral
     *
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & deformed)
    {
       assembleDual_impl<GF>(basis,deformed);
       return Base::m_assembler.rhs();
    }

    /**
     * @brief      Assembles the dual as a domain integral
     *
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(const gsMatrix<T> & points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & deformed)
    {
       return assembleDual_impl<GF>(points,basis,deformed);
    }


    void _assembleDual_expr(auto expr, const gsMultiPatch<T> & deformed, const gsMultiBasis<T> & basis);

private:
    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::DisplacementNorm, void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::Displacement, void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStrain, void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStress, void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembranePStress, void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);


    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneForce, void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF!=GoalFunction::DisplacementNorm &&
                            _GF!=GoalFunction::Displacement &&
                            _GF!=GoalFunction::MembraneStrain &&
                            _GF!=GoalFunction::MembraneStress &&
                            _GF!=GoalFunction::MembranePStress &&
                            _GF!=GoalFunction::MembraneForce
                            , void>::type
    assembleDual_impl(const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
    {
        GISMO_NO_IMPLEMENTATION;
    }

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::DisplacementNorm, gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::Displacement, gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStrain, gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStress, gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembranePStress, gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneForce, gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF!=GoalFunction::DisplacementNorm &&
                            _GF!=GoalFunction::Displacement &&
                            _GF!=GoalFunction::MembraneStrain &&
                            _GF!=GoalFunction::MembraneStress &&
                            _GF!=GoalFunction::MembranePStress &&
                            _GF!=GoalFunction::MembraneForce
                            , gsVector<T>>::type
    assembleDual_impl(const gsMatrix<T>& points, const gsMultiBasis<T> & basis, const gsMultiPatch<T> & primal)
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<int _d, bool _bending>
    typename std::enable_if<(_d==3 && _bending), T>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);


    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), T>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);

    template<enum GoalFunction _GF>
    typename std::enable_if<!(_GF==GoalFunction::Buckling || _GF==GoalFunction::Modal), T>::type
    computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed);

    template<enum GoalFunction _GF>
    typename std::enable_if<(_GF==GoalFunction::Buckling), T>::type
    computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed);

    template<enum GoalFunction _GF>
    typename std::enable_if<(_GF==GoalFunction::Modal), T>::type
    computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::DisplacementNorm, T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::Displacement, T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStrain, T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStress, T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembranePStress, T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneForce, T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF!=GoalFunction::DisplacementNorm &&
                            _GF!=GoalFunction::Displacement &&
                            _GF!=GoalFunction::MembraneStrain &&
                            _GF!=GoalFunction::MembraneStress &&
                            _GF!=GoalFunction::MembranePStress &&
                            _GF!=GoalFunction::MembraneForce
                            , T>::type
    computeGoal_impl(const gsMultiPatch<T> & deformed)
    {
        GISMO_NO_IMPLEMENTATION;
    }

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::DisplacementNorm, T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::Displacement, T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStrain, T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneStress, T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembranePStress, T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF==GoalFunction::MembraneForce, T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    // template<short_t _d, bool _bending, enum GoalFunction _GF>
    template<enum GoalFunction _GF>
    typename std::enable_if<_GF!=GoalFunction::DisplacementNorm &&
                            _GF!=GoalFunction::Displacement &&
                            _GF!=GoalFunction::MembraneStrain &&
                            _GF!=GoalFunction::MembraneStress &&
                            _GF!=GoalFunction::MembranePStress &&
                            _GF!=GoalFunction::MembraneForce
                            , T>::type
    computeGoal_impl(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed)
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
    gsSparseMatrix<T> m_massL;
    gsSparseMatrix<T> m_massH;

    gsMultiBasis<T> m_basisL, m_basisH;

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // using Base::m_patches;
    // using Base::m_basis;
    // using Base::m_bcs;

    using Base::m_patches;
    using Base::m_defpatches;
    using Base::m_assembler;
    using Base::m_materialMat;
    using Base::m_forceFun;
    using Base::m_options;
    using Base::m_foundFun;
    using Base::m_foundInd;
    using Base::m_pressFun;
    using Base::m_pressInd;

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
                      const gsMultiPatch<T> & deformed) =0;
    virtual T computeGoalEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & deformed) =0;

    virtual T computeGoal(const gsMultiPatch<T> & deformed) =0;

    virtual T computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed) =0;
};

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssemblerDWRs.hpp)
#endif
