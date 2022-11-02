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

namespace gismo
{

enum class GoalFunction : short_t
{
    Displacement        = 0,
    Stretch             = 1,
    PStrain             = 2,
    PStress             = 3,
    MembraneStrain      = 10,
    MembraneStress      = 11,
    MembraneForce       = 12,
    FlexuralStrain      = 20,
    FlexuralStress      = 21,
    FlexuralMoment      = 22,
    Modal               = 100,
    Buckling            = 200,
};

template<class T> class gsThinShellAssemblerDWRBase;
//template <class T> class gsThinShellDWRFunction;


template <short_t d, class T, bool bending>
class gsThinShellAssemblerDWR : public gsThinShellAssembler<d,T,bending>,
                                public gsThinShellAssemblerDWRBase<T>

{
public:
    typedef gsThinShellAssembler<d,T,bending> Base;
    typedef typename gsThinShellAssemblerDWRBase<T>::bContainer  bContainer;

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

    /// See \ref gsThinShellAssemblerBase for details
    gsOptionList & optionsL() {return m_assemblerL->options();}
    /// See \ref gsThinShellAssemblerBase for details
    gsOptionList & optionsH() {return m_assemblerH->options();}

    /// See \ref gsThinShellAssemblerBase for details
    gsExprAssembler<T> assemblerL() {return m_assemblerL->assembler(); }
    /// See \ref gsThinShellAssemblerBase for details
    gsExprAssembler<T> assemblerH() {return m_assemblerH->assembler(); }

    /// See \ref gsThinShellAssemblerBase for details
    void setOptions(gsOptionList & options) { m_assemblerL->setOptions(options);  m_assemblerH->setOptions(options); }

    /// See \ref gsThinShellAssemblerBase for details
    void setPointLoads(const gsPointLoads<T> & pLoads) { m_pLoads = pLoads; m_assemblerL->setPointLoads(pLoads);  m_assemblerH->setPointLoads(pLoads); }

    /// See \ref gsThinShellAssemblerBase for details
    void setFoundation(const gsFunction<T> & foundation) { m_assemblerL->setFoundation(foundation);  m_assemblerH->setFoundation(foundation); }

    /// See \ref gsThinShellAssemblerBase for details
    void setPressure(const gsFunction<T> & pressure) { m_assemblerL->setPressure(pressure);  m_assemblerH->setPressure(pressure); }

    /// See \ref gsThinShellAssemblerBase for details
    void updateBCs(const gsBoundaryConditions<T> & bconditions) { m_assemblerL->updateBCs(bconditions);  m_assemblerH->updateBCs(bconditions); }

    /// See \ref gsThinShellAssemblerBase for details
    void setBasisL(const gsMultiBasis<T> & basis) { m_assemblerL->setBasis(basis); }
    void setBasisH(const gsMultiBasis<T> & basis) { m_assemblerH->setBasis(basis); }

    /// See \ref gsThinShellAssemblerBase for details
    void setUndeformed(const gsMultiPatch<T> & patches) { m_patches = patches; m_assemblerL->setUndeformed(patches);  m_assemblerH->setUndeformed(patches); }

    /// See \ref gsThinShellAssemblerBase for details
    void homogenizeDirichlet() { m_assemblerL->homogenizeDirichlet();  m_assemblerH->homogenizeDirichlet(); }

    /// See \ref gsThinShellAssemblerBase for details
    index_t numDofsL() const { return m_assemblerL->numDofs(); }
    index_t numDofsH() const { return m_assemblerH->numDofs(); }

    void setSpaceBasisL(const gsFunctionSet<T> & spaceBasis) { m_assemblerL->setSpaceBasis(spaceBasis); }
    void setSpaceBasisH(const gsFunctionSet<T> & spaceBasis) { m_assemblerH->setSpaceBasis(spaceBasis); }


    ////////////////////////////////////////////////////////////////////////////
    //--------------------- SYSTEM ASSEMBLY ----------------------------------//
    ////////////////////////////////////////////////////////////////////////////
    void assembleL();
    void assembleH();

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
    { m_dL = _assembleDual(m_assemblerL,primal); }
    void assembleDualH(const gsMultiPatch<T> & primal)
    { m_dH = _assembleDual(m_assemblerH,primal); }
    \
    void assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal)
    { m_dL = _assembleDual(bnds,m_assemblerL,primal); }
    void assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal)
    { m_dH = _assembleDual(bnds,m_assemblerH,primal); }

    void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { m_dL = _assembleDual(points,m_assemblerL,primal); }
    void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { m_dH = _assembleDual(points,m_assemblerH,primal); }

    void assembleDualL(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { m_dL = _assembleDual(m_assemblerL,primal,deformed); }
    void assembleDualH(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { m_dH = _assembleDual(m_assemblerH,primal,deformed); }

    void assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { m_dL = _assembleDual(bnds,m_assemblerL,primal,deformed); }
    void assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { m_dH = _assembleDual(bnds,m_assemblerH,primal,deformed); }

    void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { m_dL = _assembleDual(points,m_assemblerL,primal,deformed); }
    void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { m_dH = _assembleDual(points,m_assemblerH,primal,deformed); }

    const gsSparseMatrix<T> & matrixL() const { return m_matrixL; }
    const gsSparseMatrix<T> & matrixH() const { return m_matrixH; }

    const gsSparseMatrix<T> & massL() const { return m_massL; }
    const gsSparseMatrix<T> & massH() const { return m_massH; }

    const gsMatrix<T> primalL() const { return m_pL; }
    const gsMatrix<T> primalH() const { return m_pH; }

    const gsMatrix<T> dualL() const { return m_dL; }
    const gsMatrix<T> dualH() const { return m_dH; }

    void updateMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result);
    void updateMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result);

    void constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result);
    void constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result);

    void constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed);
    void constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed);
    gsMultiPatch<T> constructSolutionL(const gsMatrix<T> & solVector);
    gsMultiPatch<T> constructSolutionH(const gsMatrix<T> & solVector);
    gsMultiPatch<T> constructDisplacementL(const gsMatrix<T> & solVector);
    gsMultiPatch<T> constructDisplacementH(const gsMatrix<T> & solVector);
    gsVector<T> constructSolutionVectorL(const gsMultiPatch<T> & deformed);
    gsVector<T> constructSolutionVectorH(const gsMultiPatch<T> & deformed);

    gsMatrix<T> projectL2_L(const gsFunction<T> &fun);
    gsMatrix<T> projectL2_H(const gsFunction<T> &fun);

    // Linear elasticity
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH);
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH);
    std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH);
    std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    // Nonlinear elasticity
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);
    std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed);
    std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    // Eigenvalues
    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal);
    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal);
    std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal);
    std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);
    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);
    std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);
    std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    T computeGoal(const gsMultiPatch<T> & deformed);
    T computeGoal(const bContainer & bnds, const gsMultiPatch<T> & deformed);
    T computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed);

    T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH) const;
    T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, const gsMultiPatch<T> &deformed) const;

    void setGoal(enum GoalFunction GF, short_t component = 9) { m_goalFunction = GF; m_component = component; }

    /// See \ref gsThinShellAssemblerBase for details
    void constructStress(const gsMultiPatch<T> & deformed, gsPiecewiseFunction<T> & result, stress_type::type type)
    {
        m_assemblerL->constructStress(deformed,result,type);
    }

    T error() const { return m_error; }
    std::vector<T> errors() const { return m_errors; }
    std::vector<T> absErrors() const
    {
        std::vector<T> result = m_errors;
        for (typename std::vector<T>::iterator it = result.begin(); it!=result.end(); it++)
            *it = std::abs(*it);
        return result;
    }

    gsDofMapper getMapperL() { return m_assemblerL->getMapper(); };
    gsDofMapper getMapperH() { return m_assemblerH->getMapper(); };

protected:

    gsSparseMatrix<T>   _assembleMass(gsThinShellAssemblerBase<T> * assembler, bool lumped = false);

    std::pair<gsSparseMatrix<T>,gsVector<T>>   _assemble(gsThinShellAssemblerBase<T> * assembler);

    gsSparseMatrix<T>   _assembleMatrix(gsThinShellAssemblerBase<T> * assembler);
    gsSparseMatrix<T>   _assembleMatrix(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed);
    gsVector<T>         _assemblePrimal(gsThinShellAssemblerBase<T> * assembler);
    gsVector<T>         _assemblePrimal(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed);

    /**
     * @brief      Assembles the dual as a domain integral
     *
     * @param      assembler  The assembler
     * @param[in]  primal     The primal
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
    {
        gsMultiPatch<T> deformed = m_patches;
        for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
            deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

        return _assembleDual(assembler,primal,deformed);
    }
    gsVector<T>         _assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);

    /**
     * @brief      Assembles the dual on boundaries
     *
     * @param[in]  points     The points
     * @param      assembler  The assembler
     * @param[in]  primal     The primal
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(const bContainer & bnds, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
    {
        gsMultiPatch<T> deformed = m_patches;
        for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
            deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

        return _assembleDual(bnds,assembler,primal,deformed);
    }
    gsVector<T>         _assembleDual(const bContainer & bnds, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);

    /**
     * @brief      Assembles the dual on points
     *
     * @param[in]  points     The points
     * @param      assembler  The assembler
     * @param[in]  primal     The primal
     * @param[in]  basis     The basis
     * @param[in]  deformed  The deformed geometry
     *
     * @return     RHS vector
     */
    gsVector<T>         _assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal)
    {
        gsMultiPatch<T> deformed = m_patches;
        for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
            deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

        return _assembleDual(points,assembler,primal,deformed);
    }
    gsVector<T>         _assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed);

    template<int _elWise>
    void computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                            std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _d, bool _bending, int _elWise>
    typename std::enable_if<(_d==3 && _bending), void>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,
                        std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _d, bool _bending, int _elWise>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,
                        std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _elWise>
    void computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,
                      std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _elWise>
    void computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,
                      std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false);

    void _applyLoadsToElWiseError(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, std::vector<T> & errors) const;
    void _applyLoadsToError(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, T & error) const;


protected:
    mutable gsThinShellAssemblerBase<T> * m_assemblerL;
    mutable gsThinShellAssemblerBase<T> * m_assemblerH;

    gsVector<T> m_pL, m_pH, m_dL, m_dH;
    gsSparseMatrix<T> m_matrixL;
    gsSparseMatrix<T> m_matrixH;
    gsSparseMatrix<T> m_massL;
    gsSparseMatrix<T> m_massH;

    gsPointLoads<T> m_pLoads;

    T m_error;
    std::vector<T> m_errors;

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsBoundaryConditions<T> m_bcs;
    // using Base::m_bcs;
    // using Base::m_basis;

    using Base::m_patches;
    // using Base::m_defpatches;
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
class gsThinShellAssemblerDWRBase //: public virtual gsThinShellAssemblerBase<T>
{
public:
    typedef gsBoxTopology::bContainer  bContainer;

    /// Default empty constructor
    virtual ~gsThinShellAssemblerDWRBase() {};

    virtual gsOptionList & optionsL() =0;
    virtual gsOptionList & optionsH() =0;

    virtual gsExprAssembler<T> assemblerL() =0;
    virtual gsExprAssembler<T> assemblerH() =0;

    virtual void setOptions(gsOptionList & options) =0;
    virtual void setPointLoads(const gsPointLoads<T> & pLoads) =0;

    virtual void setFoundation(const gsFunction<T> & foundation) =0;
    virtual void setPressure(const gsFunction<T> & pressure) =0;

    virtual void updateBCs(const gsBoundaryConditions<T> & bconditions) =0;

    virtual void setBasisL(const gsMultiBasis<T> & basis) =0;
    virtual void setBasisH(const gsMultiBasis<T> & basis) =0;

    virtual void setUndeformed(const gsMultiPatch<T> & patches) =0;

    virtual void homogenizeDirichlet() =0;

    virtual index_t numDofsL() const =0;
    virtual index_t numDofsH() const =0;

    virtual void setSpaceBasisL(const gsFunctionSet<T> & spaceBasis) =0;
    virtual void setSpaceBasisH(const gsFunctionSet<T> & spaceBasis) =0;

    virtual void assembleL() =0;
    virtual void assembleH() =0;

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

    virtual void assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal) =0;
    virtual void assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal) =0;

    virtual void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal) =0;
    virtual void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal) =0;

    virtual void assembleDualL(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;
    virtual void assembleDualH(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

    virtual void assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;
    virtual void assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

    virtual void assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;
    virtual void assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

    virtual const gsSparseMatrix<T> & matrixL() const =0;
    virtual const gsSparseMatrix<T> & matrixH() const =0;

    virtual const gsSparseMatrix<T> & massL() const =0;
    virtual const gsSparseMatrix<T> & massH() const =0;

    virtual const gsMatrix<T> primalL() const =0;
    virtual const gsMatrix<T> primalH() const =0;

    virtual const gsMatrix<T> dualL() const =0;
    virtual const gsMatrix<T> dualH() const =0;

    virtual void updateMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) =0;
    virtual void updateMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) =0;

    virtual void constructMultiPatchL(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) =0;
    virtual void constructMultiPatchH(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) =0;

    virtual void constructSolutionL(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) =0;
    virtual void constructSolutionH(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) =0;
    virtual gsMultiPatch<T> constructSolutionL(const gsMatrix<T> & solVector) =0;
    virtual gsMultiPatch<T> constructSolutionH(const gsMatrix<T> & solVector) =0;
    virtual gsMultiPatch<T> constructDisplacementL(const gsMatrix<T> & solVector) =0;
    virtual gsMultiPatch<T> constructDisplacementH(const gsMatrix<T> & solVector) =0;
    virtual gsVector<T> constructSolutionVectorL(const gsMultiPatch<T> & deformed) =0;
    virtual gsVector<T> constructSolutionVectorH(const gsMultiPatch<T> & deformed) =0;

    virtual gsMatrix<T> projectL2_L(const gsFunction<T> &fun) =0;
    virtual gsMatrix<T> projectL2_H(const gsFunction<T> &fun) =0;

    // Linear elasticity ;
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH) = 0;
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH) = 0;
    virtual std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH) = 0;
    virtual std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    // Nonlinear elasticity
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed) = 0;
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed) = 0;
    virtual std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed) = 0;
    virtual std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    // Eigenvalues
    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal) = 0;
    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal) = 0;
    virtual std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal) = 0;
    virtual std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) = 0;
    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) = 0;
    virtual std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) = 0;
    virtual std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename, unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    virtual T computeGoal(const gsMultiPatch<T> & deformed) =0;
    virtual T computeGoal(const bContainer & bnds, const gsMultiPatch<T> & deformed) =0;
    virtual T computeGoal(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed) =0;

    virtual T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH) const = 0;
    virtual T matrixNorm(const gsMultiPatch<T> &dualL, const gsMultiPatch<T> &dualH, const gsMultiPatch<T> &deformed) const = 0;

    virtual void setGoal(enum GoalFunction GF, short_t component = 9) = 0;

    /// See \ref gsThinShellAssemblerBase for details
    virtual void constructStress(const gsMultiPatch<T> & deformed,
                               gsPiecewiseFunction<T> & result,
                               stress_type::type type) = 0;

    virtual T error() const =0;
    virtual std::vector<T> errors() const =0;
    virtual std::vector<T> absErrors() const =0;

    virtual gsDofMapper getMapperL()=0;
    virtual gsDofMapper getMapperH()=0;


};

// template<short_t d, typename T, bool bending>
// class gsThinShellDWRFunction : public gsFunctionSet<T>
// {
//     typedef gsExprAssembler<>::geometryMap geometryMap;
//     typedef gsExprAssembler<>::space       space;
//     typedef gsExprAssembler<>::solution    solution;

// public:
//     gsThinShellDWRFunction( const gsFunctionSet<T> & mp,
//                             const gsFunctionSet<T> & mp_def,
//                             const gsFunctionSet<T> & dualL,
//                             const gsFunctionSet<T> & dualH,
//                             gsMaterialMatrixBase<T> * materialmatrix)
//     :
//     m_patches(&mp),
//     m_deformed(&mp_def),
//     m_dualL(&dualL),
//     m_dualH(&dualH),
//     m_materialMat(materialmatrix)
//     {

//     }

//     gsThinShellDWRFunction( const gsFunctionSet<T> & mp,
//                             const gsFunctionSet<T> & dualL,
//                             const gsFunctionSet<T> & dualH,
//                             gsMaterialMatrixBase<T> * materialmatrix)
//     :
//     m_patches(&mp),
//     m_deformed(nullptr),
//     m_dualL(&dualL),
//     m_dualH(&dualH),
//     m_materialMat(materialmatrix)
//     {

//     }

//     void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
//     {
//         this->eval_into_impl<d,bending>(u,res);
//     }

//     short_t domainDim() const { return m_patches.domainDim();}

//     short_t targetDim() const { return 1;}

// private:
//     template<int _d, bool _bending>
//     typename std::enable_if<!(_d==3 && _bending), void>::type
//     eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& res) const
//     {
//         gsWarn<<"Neumann boundary conditions left out\n";
//         gsExprEvaluator<T> ev;
//         // Geometries
//         geometryMap Gori = ev.getMap(*m_patches);           // this map is used for integrals
//         geometryMap Gdef = ev.getMap(*m_deformed);

//         gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_materialMat,m_deformed);
//         auto S0  = ev.getCoeff(S0f);

//         auto F      = ev.getCoeff(*m_forceFun, Gori);

//         auto zsolL  = ev.getCoeff(*m_dualL);
//         auto zsolH  = ev.getCoeff(*m_dualH);
//             // auto m_thick = exprAssembler.getCoeff(*m_thickFun, m_ori);

//         auto N        = S0.tr();
//         // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
//         auto Em_der   = flat( jac(Gdef).tr() * (jac(zsolH) - jac(zsolL)) ) ;

//         auto Fext = (zsolH-zsolL).tr() * F;
//         auto Fint = ( N * Em_der.tr() );

//         auto expr = ( Fext - Fint ) * meas(Gori);
//     }

//     template<int _d, bool _bending>
//     typename std::enable_if< (_d==3 && _bending), void>::type
//     eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& res) const
//     {
//         gsWarn<<"Neumann boundary conditions left out\n";
//         gsExprEvaluator<T> ev;
//         // Geometries
//         geometryMap Gori = ev.getMap(m_patches);           // this map is used for integrals
//         geometryMap Gdef = ev.getMap(deformed);

//         gsMaterialMatrixIntegrate<T,MaterialOutput::VectorN> S0f(m_materialMat,m_deformed);
//         gsMaterialMatrixIntegrate<T,MaterialOutput::VectorM> S1f(m_materialMat,m_deformed);
//         auto S0  = ev.getCoeff(S0f);
//         auto S1  = ev.getCoeff(S1f);

//         gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
//         auto m2 = ev.getCoeff(mult2t);

//         auto F      = ev.getCoeff(*m_forceFun, Gori);
//         auto zsolL  = ev.getCoeff(dualL);
//         auto zsolH  = ev.getCoeff(dualH);
//             // auto m_thick = ev.getCoeff(*m_thickFun, m_ori);

//         Base::homogenizeDirichlet();

//         auto N        = S0.tr();
//         // auto Em_der   = flat( jac(Gdef).tr() * jac(m_space) ) ;
//         auto Em_der   = flat( jac(Gdef).tr() * (jac(zsolH) - jac(zsolL)) ) ;

//         auto M        = S1.tr(); // output is a column
//         // auto Ef_der   = -( deriv2(space,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(space,Gdef) ) ) * reshape(m2,3,3);
//         auto Ef_der   = -( deriv2(zsolH,sn(Gdef).normalized().tr() )
//                         -  deriv2(zsolL,sn(Gdef).normalized().tr() ) + deriv2(Gdef,var1(zsolH,Gdef) )
//                                                                      - deriv2(Gdef,var1(zsolL,Gdef) ) ) * reshape(m2,3,3);

//         auto Fext = (zsolH-zsolL).tr() * F;
//         auto Fint = ( N * Em_der.tr() + M * Ef_der.tr() );

//         auto expr = ( Fext - Fint ) * meas(Gori);
//     }

// protected:
//     gsThinShellAssemblerDWRBase<T> * m_assembler;
//     const gsFunctionSet<T> * m_patches;
//     const gsFunctionSet<T> * m_deformed;
//     const gsFunctionSet<T> * m_dualL;
//     const gsFunctionSet<T> * m_dualH;
//     mutable gsMaterialMatrixBase<T> * m_materialMat;

// };


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssemblerDWR.hpp)
#endif
