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

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsThinShellFunctions.h>

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

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the enum: GoalFunction
   */
  void pybind11_enum_GoalFunction(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11


template<class T> class gsThinShellAssemblerDWRBase;

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

    void setSpaceBasisL(const gsMultiPatch<T> & spaceBasis) { m_assemblerL->setSpaceBasis(spaceBasis); }
    void setSpaceBasisH(const gsMultiPatch<T> & spaceBasis) { m_assemblerH->setSpaceBasis(spaceBasis); }


    ////////////////////////////////////////////////////////////////////////////
    //--------------------- SYSTEM ASSEMBLY ----------------------------------//
    ////////////////////////////////////////////////////////////////////////////
    ThinShellAssemblerStatus assembleL();
    ThinShellAssemblerStatus assembleH();

    ThinShellAssemblerStatus assembleMassL(bool lumped = false)
    { return _assembleMass(m_assemblerL, m_massL, lumped); }
    ThinShellAssemblerStatus assembleMassH(bool lumped = false)
    { return _assembleMass(m_assemblerH, m_massH, lumped); }
    ThinShellAssemblerStatus assembleMatrixL()
    { return _assembleMatrix(m_assemblerL, m_matrixL); }
    ThinShellAssemblerStatus assembleMatrixH()
    { return _assembleMatrix(m_assemblerH, m_matrixH); }
    ThinShellAssemblerStatus assembleMatrixL(const gsMultiPatch<T> & deformed)
    { return _assembleMatrix(m_assemblerL,deformed, m_matrixL); }
    ThinShellAssemblerStatus assembleMatrixH(const gsMultiPatch<T> & deformed)
    { return _assembleMatrix(m_assemblerH,deformed, m_matrixH); }

    ThinShellAssemblerStatus assemblePrimalL()
    { return _assemblePrimal(m_assemblerL, m_pL); }
    ThinShellAssemblerStatus assemblePrimalH()
    { return _assemblePrimal(m_assemblerH, m_pH); }
    ThinShellAssemblerStatus assemblePrimalL(const gsMultiPatch<T> & deformed)
    { return _assemblePrimal(m_assemblerL,deformed, m_pL); }
    ThinShellAssemblerStatus assemblePrimalH(const gsMultiPatch<T> & deformed)
    { return _assemblePrimal(m_assemblerH,deformed, m_pH); }

    ThinShellAssemblerStatus assembleDualL(const gsMultiPatch<T> & primal)
    { return _assembleDual(m_assemblerL,primal, m_dL); }
    ThinShellAssemblerStatus assembleDualH(const gsMultiPatch<T> & primal)
    { return _assembleDual(m_assemblerH,primal, m_dH); }

    ThinShellAssemblerStatus assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal)
    { return _assembleDual(bnds,m_assemblerL,primal, m_dL); }
    ThinShellAssemblerStatus assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal)
    { return _assembleDual(bnds,m_assemblerH,primal, m_dH); }

    ThinShellAssemblerStatus assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { return _assembleDual(points,m_assemblerL,primal, m_dL); }
    ThinShellAssemblerStatus assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal)
    { return _assembleDual(points,m_assemblerH,primal, m_dH); }

    ThinShellAssemblerStatus assembleDualL(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { return _assembleDual(m_assemblerL,primal,deformed, m_dL); }
    ThinShellAssemblerStatus assembleDualH(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { return _assembleDual(m_assemblerH,primal,deformed, m_dH); }

    ThinShellAssemblerStatus assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { return _assembleDual(bnds,m_assemblerL,primal,deformed, m_dL); }
    ThinShellAssemblerStatus assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { return _assembleDual(bnds,m_assemblerH,primal,deformed, m_dH); }

    ThinShellAssemblerStatus assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { return _assembleDual(points,m_assemblerL,primal,deformed, m_dL); }
    ThinShellAssemblerStatus assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed)
    { return _assembleDual(points,m_assemblerH,primal,deformed, m_dH); }

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

    T deformationNorm(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> & original)
    { return m_assemblerL->deformationNorm(deformed,original); }

    // Linear elasticity
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    // Nonlinear elasticity
    T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    // Linear elasticity
    T computeSquaredError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeSquaredErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeSquaredErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    // Nonlinear elasticity
    T computeSquaredError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeSquaredErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeSquaredErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);


    // Eigenvalues
    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);
    std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

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

    ThinShellAssemblerStatus _assembleMass(gsThinShellAssemblerBase<T> * assembler, gsSparseMatrix<T> & result, bool lumped = false);

    ThinShellAssemblerStatus _assemble(gsThinShellAssemblerBase<T> * assembler, std::pair<gsSparseMatrix<T>,gsVector<T>> & result);

    ThinShellAssemblerStatus _assembleMatrix(gsThinShellAssemblerBase<T> * assembler, gsSparseMatrix<T> & result);
    ThinShellAssemblerStatus _assembleMatrix(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed, gsSparseMatrix<T> & result);
    ThinShellAssemblerStatus _assemblePrimal(gsThinShellAssemblerBase<T> * assembler, gsVector<T> & result);
    ThinShellAssemblerStatus _assemblePrimal(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & deformed, gsVector<T> & result);

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
    ThinShellAssemblerStatus _assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, gsVector<T> & result)
    {
        gsMultiPatch<T> deformed = m_patches;
        for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
            deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

        return _assembleDual(assembler,primal,deformed,result);
    }
    ThinShellAssemblerStatus _assembleDual(gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed, gsVector<T> & result);

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
    ThinShellAssemblerStatus _assembleDual(const bContainer & bnds, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, gsVector<T> & result)
    {
        gsMultiPatch<T> deformed = m_patches;
        for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
            deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

        return _assembleDual(bnds,assembler,primal,deformed,result);
    }
    ThinShellAssemblerStatus _assembleDual(const bContainer & bnds, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed, gsVector<T> & result);

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
    ThinShellAssemblerStatus _assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, gsVector<T> & result)
    {
        gsMultiPatch<T> deformed = m_patches;
        for ( size_t k =0; k!=primal.nPatches(); ++k) // Deform the geometry
            deformed.patch(k).coefs() += primal.patch(k).coefs();  // Gdef points to mp_def, therefore updated

        return _assembleDual(points,assembler,primal,deformed,result);
    }
    ThinShellAssemblerStatus _assembleDual(const gsMatrix<T> & points, gsThinShellAssemblerBase<T> * assembler, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed, gsVector<T> & result);

    template<int _elWise>
    void computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                            std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _d, bool _bending, int _elWise>
    typename std::enable_if<(_d==3 && _bending), void>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                        std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _d, bool _bending, int _elWise>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    computeError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                        std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _elWise>
    void computeSquaredError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads,
                            std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _d, bool _bending, int _elWise>
    typename std::enable_if<(_d==3 && _bending), void>::type
    computeSquaredError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                        std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _d, bool _bending, int _elWise>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    computeSquaredError_impl(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads,
                        std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);


    template<int _elWise>
    void computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,
                      std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

    template<int _elWise>
    void computeErrorEig_impl(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,
                      std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false);

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
    using Base::m_materialMatrices;
    using Base::m_forceFun;
    using Base::m_options;
    using Base::m_foundFun;
    using Base::m_foundInd;
    using Base::m_pressFun;
    using Base::m_pressInd;

    enum GoalFunction m_goalFunction;
    short_t m_component;
};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsThinShellAssemblerDWR
   */
  void pybind11_init_gsThinShellAssemblerDWR2(pybind11::module &m);
  void pybind11_init_gsThinShellAssemblerDWR3(pybind11::module &m);
  void pybind11_init_gsThinShellAssemblerDWR3nb(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

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

    virtual void setSpaceBasisL(const gsMultiPatch<T> & spaceBasis) =0;
    virtual void setSpaceBasisH(const gsMultiPatch<T> & spaceBasis) =0;

    virtual ThinShellAssemblerStatus assembleL() =0;
    virtual ThinShellAssemblerStatus assembleH() =0;

    virtual ThinShellAssemblerStatus assembleMassL(bool lumped = false) =0;
    virtual ThinShellAssemblerStatus assembleMassH(bool lumped = false) =0;
    virtual ThinShellAssemblerStatus assembleMatrixL() =0;
    virtual ThinShellAssemblerStatus assembleMatrixH() =0;
    virtual ThinShellAssemblerStatus assembleMatrixL(const gsMultiPatch<T> & deformed) =0;
    virtual ThinShellAssemblerStatus assembleMatrixH(const gsMultiPatch<T> & deformed) =0;

    virtual ThinShellAssemblerStatus assemblePrimalL() =0;
    virtual ThinShellAssemblerStatus assemblePrimalH() =0;
    virtual ThinShellAssemblerStatus assemblePrimalL(const gsMultiPatch<T> & deformed) =0;
    virtual ThinShellAssemblerStatus assemblePrimalH(const gsMultiPatch<T> & deformed) =0;


    virtual ThinShellAssemblerStatus assembleDualL(const gsMultiPatch<T> & primal) =0;
    virtual ThinShellAssemblerStatus assembleDualH(const gsMultiPatch<T> & primal) =0;

    virtual ThinShellAssemblerStatus assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal) =0;
    virtual ThinShellAssemblerStatus assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal) =0;

    virtual ThinShellAssemblerStatus assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal) =0;
    virtual ThinShellAssemblerStatus assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal) =0;

    virtual ThinShellAssemblerStatus assembleDualL(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;
    virtual ThinShellAssemblerStatus assembleDualH(const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

    virtual ThinShellAssemblerStatus assembleDualL(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;
    virtual ThinShellAssemblerStatus assembleDualH(const bContainer & bnds, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

    virtual ThinShellAssemblerStatus assembleDualL(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;
    virtual ThinShellAssemblerStatus assembleDualH(const gsMatrix<T> & points, const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed) =0;

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

    virtual T deformationNorm(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> & original) =0;

    // Linear elasticity ;
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    // Nonlinear elasticity
    virtual T computeError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    // Linear elasticity ;
    virtual T computeSquaredError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeSquaredErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeSquaredErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    // Nonlinear elasticity
    virtual T computeSquaredError(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeSquaredErrorElements(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeSquaredErrorDofs(const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH, const gsMultiPatch<T> & deformed, bool withLoads=false, std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;


    // Eigenvalues
    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

    virtual T computeErrorEig(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigElements(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;
    virtual std::vector<T> computeErrorEigDofs(const T evPrimalL, const T evDualL, const T evDualH,
                      const gsMultiPatch<T> & dualL, const gsMultiPatch<T> & dualH,
                      const gsMultiPatch<T> & primal, const gsMultiPatch<T> & deformed,std::string filename = std::string(), unsigned np=1000, bool parametric=false, bool mesh=false) = 0;

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

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsThinShellAssemblerDWRBase
   */
  void pybind11_init_gsThinShellAssemblerDWRBase(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssemblerDWR.hpp)
#endif
