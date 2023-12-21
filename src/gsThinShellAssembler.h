/** @file gsThinShellAssembler.h

    @brief Provides linear and nonlinear assemblers for thin shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsThinShellFunctions.h>

#include <gsPde/gsPointLoads.h>

#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

namespace gismo
{

/**
 * @brief      Defines the coupling type over interfaces
 */
struct shell_coupling
{
    enum type
    {
        automatic           = -1, // (applies weak_penalty when needed, strong_Dpatch elsewhere)
        weak_penalty        = 0,
        weak_nitsche        = 1,
        strong_C0           = 2,
        strong_Dpatch       = 3,
        strong_ASG1         = 4,
    };

    enum continuity
    {
        Cminus1             = -1,
        C0                  = 0,
        C1                  = 1,
    };
};

enum class ThinShellAssemblerStatus
{
    Success,         ///< Assembly is successful
    AssemblyError,   ///< Assembly failed due to an error in the expression (e.g. overflow)
    DimensionError   ///< Assembly failed due to a dimension error
};


template<class T> class gsThinShellAssemblerBase;
/**
 * @brief      Assembles the system matrix and vectors for 2D and 3D shell
 *             problems, including geometric nonlinearities and loading
 *             nonlinearities. The material nonlinearities are handled by the
 *             @ref gsMaterialMatrixIntegrate class.
 *
 * @tparam     d        The dimension (2 = planar, 3 = surface)
 * @tparam     T        Real type
 * @tparam     bending  True: Assemble bending terms; False: Do not assemble bending terms
 *
 *
 * @ingroup    KLShell
 */
template <short_t d, class T, bool bending>
class gsThinShellAssembler : public gsThinShellAssemblerBase<T>
{
public:
    typedef gsBoxTopology::ifContainer ifContainer;

public:

    /**
     * @brief      Constructor for the shell assembler
     *
     * @param[in]  patches         The geometry
     * @param[in]  basis           The basis
     * @param[in]  bconditions     The boundary condition
     * @param[in]  surface_force   The surface force
     * @param      materialmatrix  The material matrix container
     */
    gsThinShellAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & basis,
                        const gsBoundaryConditions<T> & bconditions,
                        const gsFunction<T> & surface_force,
                        const gsMaterialMatrixContainer<T> & materialmatrices);

    /**
     * @brief      Constructor for the shell assembler
     *
     * @param[in]  patches         The geometry
     * @param[in]  basis           The basis
     * @param[in]  bconditions     The boundary condition
     * @param[in]  surface_force   The surface force
     * @param      materialmatrix  The material matrix class
     */
    gsThinShellAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & basis,
                        const gsBoundaryConditions<T> & bconditions,
                        const gsFunction<T> & surface_force,
                        gsMaterialMatrixBase<T> * materialmatrix);

    /**
     * @brief      Constructor for te shell assembler
     *
     * @param[in]  patches         The geometry
     * @param[in]  basis           The basis
     * @param[in]  bconditions     The boundary condition
     * @param[in]  surface_force   The surface force
     * @param      materialmatrix  The material matrix class
     */
    // gsThinShellAssembler(const gsMultiPatch<T> & patches,
    //                     const gsMultiBasis<T> & basis,
    //                     const gsBoundaryConditions<T> & bconditions,
    //                     const gsFunction<T> & surface_force,
    //                     const gsPiecewiseFunction<T> & materialmatrices);

    /// Default empty constructor
    gsThinShellAssembler() { }

    /// Copy constructor (makes deep copy)
    gsThinShellAssembler( const gsThinShellAssembler& other )
    {
        operator=(other);
    }

    /// Move constructor
    gsThinShellAssembler( gsThinShellAssembler&& other )
    {
        operator=(give(other));
    }

    /// Assignment operator
    gsThinShellAssembler& operator= ( const gsThinShellAssembler& other );

    /// Move assignment operator
    gsThinShellAssembler& operator= ( gsThinShellAssembler&& other );

    /// See \ref gsThinShellAssemblerBase for details
    gsOptionList & options() {return m_options;}

    /// See \ref gsThinShellAssemblerBase for details
    gsExprAssembler<T> assembler() {return m_assembler; }

    /// See \ref gsThinShellAssemblerBase for details
    void setOptions(gsOptionList & options);

    //--------------------- PROBLEM FORMULATION-------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    void setPointLoads(const gsPointLoads<T> & pLoads){ m_pLoads = pLoads; }
    const gsPointLoads<T> & getPointLoads(){ return m_pLoads; }
    void setPointMass(const gsPointLoads<T> & pMass){ m_pMass = pMass; }
    const gsPointLoads<T> & getPointMass(){ return m_pMass; }

    /// See \ref gsThinShellAssemblerBase for details
    void setFoundation(const gsFunction<T> & foundation) { m_foundFun = &foundation; m_foundInd = true; }

    /// See \ref gsThinShellAssemblerBase for details
    void setPressure(const gsFunction<T> & pressure) { m_pressFun = &pressure; m_pressInd = true; }

    /// See \ref gsThinShellAssemblerBase for details
    void updateBCs(const gsBoundaryConditions<T> & bconditions)
    {
        m_bcs = bconditions;
        space m_space = m_assembler.trialSpace(0);
        this->_assembleDirichlet();

        m_ddofs = m_space.fixedPart();
        m_mapper = m_space.mapper();
    }

    /// See \ref gsThinShellAssemblerBase for details
    void homogenizeDirichlet();

    /// See \ref gsThinShellAssemblerBase for details
    index_t numDofs() const {return m_assembler.numDofs();}

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus status() const { return m_status; }

    ////////////////////////////////////////////////////////////////////////////
    //--------------------- SYSTEM ASSEMBLY ----------------------------------//
    ////////////////////////////////////////////////////////////////////////////

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemble();

    /// See \ref gsThinShellAssemblerBase for details
    const gsFunctionSet<T> & getSpaceBasis() { return *m_spaceBasis; }

    /// See \ref gsThinShellAssemblerBase for details
    gsDofMapper getMapper() const { return m_mapper; } ;

private:
    /// Specialisation of assemble() for surfaces (3D)
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type assemble_impl();

    /// Specialisation of assemble() for planar geometries (2D)
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type assemble_impl();

public:
    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleMass(const bool lumped = false);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleFoundation();

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemble(const gsFunctionSet<T> & deformed, const bool matrix = true, const bool homogenize = true);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemble(const gsMatrix<T>     & solVector, const bool matrix = true, const bool homogenize = true);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleMatrix(const gsFunctionSet<T>   & deformed  );

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleMatrix(const gsMatrix<T>       & solVector );

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleMatrix(const gsFunctionSet<T> & deformed, const gsFunctionSet<T> & previous, gsMatrix<T> & update);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleMatrix(const gsMatrix<T> & solVector, const gsMatrix<T> & prevVector);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureMatrix(const gsFunction<T>   & pressFun  );

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureMatrix(const T pressure);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureMatrix(const gsFunction<T>   & pressFun, const gsFunctionSet<T> & deformed);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureMatrix(const T pressure, const gsFunctionSet<T> & deformed);

private:
    /// Implementation of assembleMatrix for surfaces (3D)
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
    assembleMatrix_impl(const gsFunctionSet<T>   & deformed  );

    /// Implementation of assembleMatrix for planar geometries (2D)
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
    assembleMatrix_impl(const gsFunctionSet<T>   & deformed  );

    /// Implementation of assembleMatrix for surfaces (3D)
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
    assembleMatrix_impl(const gsFunctionSet<T> & deformed, const gsFunctionSet<T> & previous, gsMatrix<T> & update);

    /// Implementation of assembleMatrix for planar geometries (2D)
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
    assembleMatrix_impl(const gsFunctionSet<T> & deformed, const gsFunctionSet<T> & previous, gsMatrix<T> & update)
    { GISMO_NO_IMPLEMENTATION; }

public:
    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleVector(const gsFunctionSet<T>   & deformed, const bool homogenize = true);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleVector(const gsMatrix<T>       & solVector, const bool homogenize = true);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureVector(const gsFunction<T>   & pressFun  );

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureVector(const T pressure);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureVector(const gsFunction<T>   & pressFun, const gsFunctionSet<T> & deformed);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureVector(const T pressure, const gsFunctionSet<T> & deformed);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assemblePressureVector(const gsFunctionSet<T> & deformed);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleFoundationVector(const gsFunctionSet<T> & deformed);

    /// See \ref gsThinShellAssemblerBase for details
    ThinShellAssemblerStatus assembleFoundationVector(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed);

private:
    /// Implementation of assembleVector for surfaces (3D)
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, ThinShellAssemblerStatus>::type
    assembleVector_impl(const gsFunctionSet<T>   & deformed, const bool homogenize);

    /// Implementation of assembleVector for planar geometries (2D)
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), ThinShellAssemblerStatus>::type
    assembleVector_impl(const gsFunctionSet<T>   & deformed, const bool homogenize);

public:
    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> boundaryForce(const gsFunctionSet<T> & deformed, const std::vector<patchSide> & patchSides) const;
    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> boundaryForce(const gsFunctionSet<T> & deformed, const patchSide & ps) const
    {
        std::vector<patchSide> patchSides(1);
        patchSides[0] = ps;
        return boundaryForce(deformed,patchSides);
    }


private:
    /// Implementation of the boundary force vector for surfaces (3D)
    template<short_t _d, bool _bending>
    typename std::enable_if<(_d==3) && _bending, gsMatrix<T> >::type
    boundaryForce_impl(const gsFunctionSet<T> & deformed, const std::vector<patchSide> & patchSides) const;

    /// Implementation of the boundary force vector for planar geometries (2D)
    template<short_t _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), gsMatrix<T> >::type
    boundaryForce_impl(const gsFunctionSet<T> & deformed, const std::vector<patchSide> & patchSides) const;

public:

    //--------------------- GEOMETRY ACCESS --------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    const gsMultiPatch<T> & geometry()      const  {return m_patches;}
    gsMultiPatch<T> & geometry()  {return m_patches;}
    void setGeometry(const gsMultiPatch<T> & patches)
    {
        this->_getOptions();
        m_patches = patches;
        this->_initialize();
    }

    /// See \ref gsThinShellAssemblerBase for details
    void setUndeformed(const gsMultiPatch<T> & patches)
    {
        this->setGeometry(patches);
    }

    //--------------------- BASIS ACCESS --------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    const gsMultiBasis<T> & basis()      const  {return m_basis;}
    const gsMultiBasis<T> & getBasis()   const  { return this->basis(); }

    gsMultiBasis<T> & basis()                   {return m_basis;}
    gsMultiBasis<T> & getBasis()                { return this->basis(); }

    void setBasis(const gsMultiBasis<T> & basis)
    {
        m_basis = basis;
        this->_getOptions();
        this->_initialize();
    }

    /// See \ref gsThinShellAssemblerBase for details
    void setSpaceBasis(const gsFunctionSet<T> & spaceBasis)
    {
        m_spaceBasis = &spaceBasis;
        this->_getOptions();
        this->_initialize();
    }

    // / See \ref gsThinShellAssemblerBase for details
    // const gsFunctionSet<T> & defGeometry()   const  {return *m_defpatches;}

    /// See \ref gsThinShellAssemblerBase for details
    T getArea(const gsFunctionSet<T> & geometry);

    /// See \ref gsThinShellAssemblerBase for details
    T getDisplacementNorm(const gsFunctionSet<T> & deformed);

    /// See \ref gsThinShellAssemblerBase for details
    T getElasticEnergy(const gsFunctionSet<T> & deformed);

    //--------------------- MATERIAL ACCESS --------------------------------//
    gsMaterialMatrixContainer<T> materials()    const  {return m_materialMatrices;}
    gsMaterialMatrixBase<T> * material(const index_t p)    const  {return m_materialMatrices.piece(p);}

    //--------------------- SYSTEM ACCESS ----------------------------------//
    const gsSparseMatrix<T> & matrix()      const   {return m_assembler.matrix();}
    gsSparseMatrix<T> & massMatrix() {return m_mass;}

    const gsMatrix<T>       & rhs()         const {return m_rhs.size()==0 ? m_assembler.rhs() : m_rhs;}
    // const gsMatrix<T>       & rhs()     const {return m_assembler.rhs();}

    //--------------------- INTERFACE HANDLING -----------------------------//
    void addStrongC0(const gsBoxTopology::ifContainer & interfaces);
    void addStrongC1(const gsBoxTopology::ifContainer & interfaces);
    void addWeakC0(const gsBoxTopology::ifContainer & interfaces);
    void addWeakC1(const gsBoxTopology::ifContainer & interfaces);
    void addUncoupled(const gsBoxTopology::ifContainer & interfaces);
    void initInterfaces();

    //--------------------- SOLUTION CONSTRUCTION --------------------------//
    gsMultiPatch<T> constructMultiPatch(const gsMatrix<T> & solVector) const;
    void updateMultiPatch(const gsMatrix<T> & solVector, gsMultiPatch<T> & mp) const;

    /// See \ref gsThinShellAssemblerBase for details
protected:
    gsMultiPatch<T> _constructSolution(const gsMatrix<T> & solVector, const gsMultiPatch<T> & undeformed) const;
public:
    // gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector, const gsMultiPatch<T> & undeformed) const;
    gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const;

    /// See \ref gsThinShellAssemblerBase for details
    void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    /// See \ref gsThinShellAssemblerBase for details
    gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const;

    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> fullSolutionVector(const gsMatrix<T> & vector) const;

    /// See \ref gsThinShellAssemblerBase for details
    void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    /// See \ref gsThinShellAssemblerBase for details
    gsVector<T> constructSolutionVector(const gsMultiPatch<T> & deformed) const;

    //--------------------- SPECIALS ----------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    void constructStress(   const gsFunctionSet<T> & deformed,
                            gsPiecewiseFunction<T> & result,
                            stress_type::type type);

    /// See \ref gsThinShellAssemblerBase for details
    void constructStress(   const gsFunctionSet<T> & original,
                            const gsFunctionSet<T> & deformed,
                            gsPiecewiseFunction<T> & result,
                            stress_type::type type);

    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> computePrincipalStretches(const gsMatrix<T> & points, const gsFunctionSet<T> & deformed, const T z=0);

    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> computePrincipalStresses(const gsMatrix<T> & points, const gsFunctionSet<T> & deformed, const T z=0);

    /// See \ref gsThinShellAssemblerBase for details
    void projectL2_into(const gsFunction<T> &fun, gsMatrix<T> & result);

    /// See \ref gsThinShellAssemblerBase for details
    void projectL2_into(const gsFunction<T> &fun, gsMultiPatch<T> & result);

    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> projectL2(const gsFunction<T> &fun);

    /// See \ref gsThinShellAssemblerBase for details
    void plotSolution(std::string string, const gsMatrix<T> & solVector);

    /// See \ref gsThinShellAssemblerBase for details
    T deformationNorm(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> & original);

    /// See \ref gsThinShellAssemblerBase for details
    gsDofMapper getMapper() { return m_mapper; };

    /// See \ref gsThinShellAssemblerBase for details
    T interfaceErrorC0(const gsFunctionSet<T> & deformed)
    { return interfaceErrorC0(deformed,m_patches.interfaces()); }
    T interfaceErrorC0(const gsFunctionSet<T> & deformed, const ifContainer & iFaces);

    /// See \ref gsThinShellAssemblerBase for details
    T interfaceErrorG1(const gsFunctionSet<T> & deformed)
    { return interfaceErrorG1(deformed,m_patches.interfaces()); }
    T interfaceErrorG1(const gsFunctionSet<T> & deformed, const ifContainer & iFaces);

    /// See \ref gsThinShellAssemblerBase for details
    T interfaceErrorNormal(const gsFunctionSet<T> & deformed)
    { return interfaceErrorNormal(deformed,m_patches.interfaces()); }
    T interfaceErrorNormal(const gsFunctionSet<T> & deformed, const ifContainer & iFaces);

    /// See \ref gsThinShellAssemblerBase for details
    T interfaceErrorGaussCurvature(const gsFunctionSet<T> & deformed)
    { return interfaceErrorGaussCurvature(deformed,m_patches.interfaces()); }
    T interfaceErrorGaussCurvature(const gsFunctionSet<T> & deformed, const ifContainer & iFaces);

    /// See \ref gsThinShellAssemblerBase for details
    T interfaceErrorMeanCurvature(const gsFunctionSet<T> & deformed)
    { return interfaceErrorMeanCurvature(deformed,m_patches.interfaces()); }
    T interfaceErrorMeanCurvature(const gsFunctionSet<T> & deformed, const ifContainer & iFaces);

protected:
    /// Initializes the method
    void _initialize();
    void _defaultOptions();
    void _getOptions();

    void _assembleNeumann();

    template <bool _matrix>
    void _assemblePressure(const gsFunction<T> & pressFun);
    template <bool _matrix>
    void _assemblePressure(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    template <bool _matrix>
    void _assembleFoundation(const gsFunction<T> & foundFun);
    template <bool _matrix>
    void _assembleFoundation(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed);

    template <bool _matrix>
    void _assembleWeakBCs();
    template <bool _matrix>
    void _assembleWeakBCs(const gsFunctionSet<T> & deformed);

    template <bool _matrix>
    void _assembleWeakIfc();
    template <bool _matrix>
    void _assembleWeakIfc(const gsFunctionSet<T> & deformed);

    void _assembleDirichlet();

    void _applyLoads();
    void _applyMass();

    void _ifcTest(const T tol = 1e-2);
    bool _isInPlane(const boundaryInterface & ifc, const T tol = 1e-2);

private:
    template<short_t _d>
    typename std::enable_if<(_d==3), void>::type
    _assembleNeumann_impl();

    template<short_t _d>
    typename std::enable_if<!(_d==3), void>::type
    _assembleNeumann_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3), void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3), void>::type
    _assemblePressure_impl(const gsFunction<T> & pressFun, const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assembleFoundation_impl(const gsFunction<T> & foundFun);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assembleFoundation_impl(const gsFunction<T> & foundFun);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3), void>::type
    _assembleFoundation_impl(const gsFunction<T> & foundFun);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assembleFoundation_impl(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assembleFoundation_impl(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3), void>::type
    _assembleFoundation_impl(const gsFunction<T> & foundFun, const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assembleWeakBCs_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assembleWeakBCs_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && _matrix, void>::type
    _assembleWeakBCs_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && !_matrix, void>::type
    _assembleWeakBCs_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assembleWeakBCs_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assembleWeakBCs_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && _matrix, void>::type
    _assembleWeakBCs_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && !_matrix, void>::type
    _assembleWeakBCs_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assembleWeakIfc_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assembleWeakIfc_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && _matrix, void>::type
    _assembleWeakIfc_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && !_matrix, void>::type
    _assembleWeakIfc_impl();

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && _matrix, void>::type
    _assembleWeakIfc_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<(_d==3) && !_matrix, void>::type
    _assembleWeakIfc_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && _matrix, void>::type
    _assembleWeakIfc_impl(const gsFunctionSet<T> & deformed);

    template<short_t _d, bool _matrix>
    typename std::enable_if<!(_d==3) && !_matrix, void>::type
    _assembleWeakIfc_impl(const gsFunctionSet<T> & deformed);

protected:
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    typedef gsExprAssembler<>::element     element;

    gsDofMapper m_mapper;

    gsExprAssembler<> m_assembler;
    gsExprEvaluator<> m_evaluator;

    gsMultiPatch<T> m_patches;
    // const gsFunctionSet<T> * m_defpatches;
    const gsFunctionSet<T> * m_itpatches;
    mutable gsMultiBasis<T> m_basis;
    const gsFunctionSet<T> *  m_spaceBasis;
    gsBoundaryConditions<T> m_bcs;

    mutable gsMatrix<T> m_ddofs;

    gsSparseMatrix<T> m_mass;

    const gsFunction<T> * m_forceFun;
    const gsFunction<T> * m_foundFun;
    const gsFunction<T> * m_pressFun;

    gsMaterialMatrixContainer<T> m_materialMatrices;

    gsPointLoads<T>  m_pLoads, m_pMass;

    mutable gsMatrix<T> m_solvector;

    gsMatrix<T> m_rhs;

    mutable gsOptionList m_options;

    mutable bool m_foundInd;
    mutable bool m_pressInd;

    mutable index_t m_continuity;

    mutable T m_alpha_d_bc,m_alpha_r_bc,m_alpha_d_ifc,m_alpha_r_ifc;
    mutable index_t m_IfcDefault;

    mutable ifContainer m_inPlane, m_outPlane, m_uncoupled, m_strongC0, m_weakC0, m_strongC1, m_weakC1, m_unassigned;

    mutable ThinShellAssemblerStatus m_status;
};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsThinShellAssembler
   */
  void pybind11_init_gsThinShellAssembler2(pybind11::module &m);
  void pybind11_init_gsThinShellAssembler3(pybind11::module &m);
  void pybind11_init_gsThinShellAssembler3nb(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

/**
 * @brief      Base class for the gsThinShellAssembler
 *
 * @tparam     T     Real type
 *
 * @ingroup    KLShell
 */
template <class T>
class gsThinShellAssemblerBase
{
public:
    /// Default deconstructor
    typedef gsBoxTopology::ifContainer ifContainer;

public:

    /// Default deconstructor
    gsThinShellAssemblerBase() {};

    /// Default empty constructor
    virtual ~gsThinShellAssemblerBase() {};

    /// Returns the options of the assembler
    virtual gsOptionList & options() = 0;

    /// Returns the internal expression assembler
    virtual gsExprAssembler<T> assembler() =0;

    /// Sets the options of the assembler
    virtual void setOptions(gsOptionList & options) = 0;

    /// Registers a \ref gsPointLoads object for point loads acting on the shell
    virtual void setPointLoads(const gsPointLoads<T> & pLoads) = 0;
    /// Gets the registered point loads
    virtual const gsPointLoads<T> & getPointLoads() = 0;
    /// Registers a \ref gsPointLoads object for a point mass acting on the shell. The point masss must be 1-dimensional
    virtual void setPointMass(const gsPointLoads<T> & pLoads) = 0;
    /// Gets the registered point masses
    virtual const gsPointLoads<T> & getPointMass() = 0;

    /**
     * @brief      Registers a stiffness function to be used for handling an elastic foundation, only relevant for 3D shells, with out-of-plane deformations
     *
     * @param[in]  foundation  The foundation stiffness function (3D vector)
     */
    virtual void setFoundation(const gsFunction<T> & foundation) = 0;

    /**
     * @brief      Registers a scalar function acting as pressure (in normal direction) in normal direction
     *
     * @note       Since the pressure acts in normal direction, and since the normals of the shell change under deformation, this load is nonlinear!
     *
     * @param[in]  pressure  The pressure
     */
    virtual void setPressure(const gsFunction<T> & pressure) = 0;

    /**
     * @brief      Overwrites the boundary conditions
     *
     * @param[in]  bconditions  The boundary conditions
     */
    virtual void updateBCs(const gsBoundaryConditions<T> & bconditions) = 0;

    /**
     * @brief      Overwrites the basis
     *
     * @param[in]  bconditions  The basis
     */
    virtual void setBasis(const gsMultiBasis<T> & basis) = 0;

    /**
     * @brief      Gets the basis
     *
     */
    virtual const gsMultiBasis<T> & getBasis() = 0;

    /**
     * @brief      Overwrites the undeformed geometry
     *
     * @param[in]  bconditions  The undeformed geometry
     */
    virtual void setUndeformed(const gsMultiPatch<T> & patches) = 0;

    /// Sets the Dirichlet BCs to zero
    virtual void homogenizeDirichlet() = 0;

    /// Returns the number of degrees of freedom in the assembler.
    virtual index_t numDofs() const  = 0;

    /// Returns the assembler status
    virtual ThinShellAssemblerStatus status() const = 0;

    /// Assembles the linear system and corresponding right-hand side
    virtual ThinShellAssemblerStatus assemble() = 0;

    /// Set the basis that is used for assembly (but not for quadrature!)
    virtual void setSpaceBasis(const gsFunctionSet<T> & spaceBasis) = 0;

    /// Get the basis that is used for assembly (but not for quadrature!)
    virtual const gsFunctionSet<T> & getSpaceBasis()  = 0;
    // virtual const gsMultiBasis<T> & getSpaceBasis() = 0;

    /// Assembles the mass matrix (including density and thickness!); if lumped=true, a lumped mass matrix will be constructed,
    virtual ThinShellAssemblerStatus assembleMass(const bool lumped = false) = 0;

    /// Assembles the elastic foundation matrix
    virtual ThinShellAssemblerStatus assembleFoundation() = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix and the residual for an iteration of Newton's method
     *
     * @note       .rhs() returns the negative residual!
     *
     * @param[in]  deformed  The deformed multipatch
     * @param[in]  Matrix    True if the matrix should be assembled
     */
    virtual ThinShellAssemblerStatus assemble(const gsFunctionSet<T> & deformed, const bool matrix = true, const bool homogenize = true) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix and the residual for an iteration of Newton's method
     *
     * @note       .rhs() returns the negative residual!
     *
     * @param[in]  deformed  The solution vector
     * @param[in]  Matrix    True if the matrix should be assembled
     */
    virtual ThinShellAssemblerStatus assemble(const gsMatrix<T>     & solVector, const bool matrix = true, const bool homogenize = true) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix (nonlinear)
     *
     * @param[in]  deformed  The deformed geometry
     */
    virtual ThinShellAssemblerStatus assembleMatrix(const gsFunctionSet<T>   & deformed  ) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix (nonlinear)
     *
     * @param[in]  deformed  The solution vector
     */
    virtual ThinShellAssemblerStatus assembleMatrix(const gsMatrix<T>       & solVector ) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix (nonlinear) using the Mixed Integration Point (MIP) method
     *
     * For more details, see
     *  Leonetti, L., Magisano, D., Madeo, A., Garcea, G., Kiendl, J., & Reali, A. (2019).
     *  A simplified Kirchhoff–Love large deformation model for elastic shells and its effective isogeometric formulation.
     *  Computer Methods in Applied Mechanics and Engineering, 354, 369–396.
     *  https://doi.org/10.1016/j.cma.2019.05.025
     *
     * @param[in]  deformed  The deformed geometry
     * @param[in]  previous  The previous geometry
     * @param      update    The update vector
     */
    virtual ThinShellAssemblerStatus assembleMatrix(const gsFunctionSet<T> & deformed, const gsFunctionSet<T> & previous, gsMatrix<T> & update) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix (nonlinear) using the Mixed Integration Point (MIP) method
     *
     * For more details, see
     *  Leonetti, L., Magisano, D., Madeo, A., Garcea, G., Kiendl, J., & Reali, A. (2019).
     *  A simplified Kirchhoff–Love large deformation model for elastic shells and its effective isogeometric formulation.
     *  Computer Methods in Applied Mechanics and Engineering, 354, 369–396.
     *  https://doi.org/10.1016/j.cma.2019.05.025
     *
     * @param[in]  solVector   The current  solution vector
     * @param[in]  prevVector  The previous solution vector
     */
    virtual ThinShellAssemblerStatus assembleMatrix(const gsMatrix<T> & solVector, const gsMatrix<T> & prevVector) = 0;

    /**
     * @brief      Assembles the residual vector
     *
     * @param[in]  deformed  The deformed geometry
     */
    virtual ThinShellAssemblerStatus assembleVector(const gsFunctionSet<T>   & deformed, const bool homogenize = true) = 0;

    /**
     * @brief      Assembles the residual vector
     *
     * @param[in]  deformed  The solution vector
     */
    virtual ThinShellAssemblerStatus assembleVector(const gsMatrix<T>       & solVector, const bool homogenize = true) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system matrix (linear)
     *
     * @param[in]  pressFun  The pressure function
     */
    virtual ThinShellAssemblerStatus assemblePressureMatrix(const gsFunction<T>   & pressFun  ) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system matrix (linear)
     *
     * @param[in]  pressure  The pressure value
     */
    virtual ThinShellAssemblerStatus assemblePressureMatrix(const T pressure) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system matrix (non-linear)
     *
     * @param[in]  pressFun  The pressure function
     * @param[in]  deformed  The deformed shape
     */
    virtual ThinShellAssemblerStatus assemblePressureMatrix(const gsFunction<T>   & pressFun, const gsFunctionSet<T> & deformed) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system matrix (non-linear)
     *
     * @param[in]  pressure  The pressure value
     * @param[in]  deformed  The deformed shape
     */
    virtual ThinShellAssemblerStatus assemblePressureMatrix(const T pressure, const gsFunctionSet<T> & deformed) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system vector (linear)
     *
     * @param[in]  pressFun  The pressure function
     */
    virtual ThinShellAssemblerStatus assemblePressureVector(const gsFunction<T>   & pressFun  ) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system vector (linear)
     *
     * @param[in]  pressure  The pressure value
     */
    virtual ThinShellAssemblerStatus assemblePressureVector(const T pressure) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system vector (non-linear)
     *
     * @param[in]  pressFun  The pressure value
     * @param[in]  deformed  The deformed shape
     */
    virtual ThinShellAssemblerStatus assemblePressureVector(const gsFunction<T>   & pressFun, const gsFunctionSet<T> & deformed) = 0;

    /**
     * @brief      Assembles the pressure contribution in the system vector (non-linear)
     *
     * @param[in]  pressure  The pressure value
     * @param[in]  deformed  The deformed shape
     */
    virtual ThinShellAssemblerStatus assemblePressureVector(const T pressure, const gsFunctionSet<T> & deformed) = 0;

    /**
     * @brief      Computes the force on a set of boundaries
     *
     *             This function is typically used when you want to know the
     *             load on a boundary on which a displacement is applied
     *
     * @param[in]  deformed    The deformed geometry
     * @param      patchSides  patch sides
     *
     * @return     The sim of the loads on the control points.
     */
    virtual gsMatrix<T> boundaryForce(const gsFunctionSet<T>   & deformed , const std::vector<patchSide>& patchSides) const = 0;
    virtual gsMatrix<T> boundaryForce(const gsFunctionSet<T>   & deformed , const patchSide & ps) const = 0;

    /// Returns the undeformed geometry
    virtual const gsMultiPatch<T> & geometry()    const = 0;
    virtual gsMultiPatch<T> & geometry() = 0;

    ///
    virtual const gsMultiBasis<T> & basis()      const  = 0;
    virtual gsMultiBasis<T> & basis() = 0;

    virtual void setGeometry(const gsMultiPatch<T> & patches) = 0;

    // /// Returns the deformed geometry
    // virtual const gsFunctionSet<T> & defGeometry() const = 0;

    /// Returns the material matrices used in the class
    virtual gsMaterialMatrixContainer<T> materials()          const = 0;
    /// Returns the material matrix on patch p used in the class
    virtual gsMaterialMatrixBase<T> * material(const index_t p)          const = 0;

    /// Returns the area of \a geometry
    virtual T getArea(const gsFunctionSet<T> & geometry) = 0;

    /// Returns the displacement norm, i.e. norm = sqrt(u'*u/area)
    virtual T getDisplacementNorm(const gsFunctionSet<T> & deformed) = 0;

    /// Returns the elastic energy norm, i.e. norm = 0.5 * u'*F_int
    virtual T getElasticEnergy(const gsFunctionSet<T> & deformed) = 0;

    /// Returns a reference to the system matrix that is assembled
    virtual const gsSparseMatrix<T> & matrix()  const  = 0;

    /// Returns a reference to the mass matrix that is assembled
    virtual gsSparseMatrix<T> & massMatrix() = 0;

    /// Returns a reference to the right-hand side vector that is assembled
    virtual const gsMatrix<T>       & rhs()     const  = 0;

    //--------------------- INTERFACE HANDLING -----------------------------//
    virtual void addStrongC0(const gsBoxTopology::ifContainer & interfaces) = 0;
    virtual void addStrongC1(const gsBoxTopology::ifContainer & interfaces) = 0;
    virtual void addWeakC0(const gsBoxTopology::ifContainer & interfaces) = 0;
    virtual void addWeakC1(const gsBoxTopology::ifContainer & interfaces) = 0;
    virtual void addUncoupled(const gsBoxTopology::ifContainer & interfaces) = 0;
    virtual void initInterfaces() = 0;

    /// Construct solution field from computed solution vector \a solVector and returns a multipatch
    virtual gsMultiPatch<T> constructMultiPatch(const gsMatrix<T> & solVector) const = 0;
    virtual void updateMultiPatch(const gsMatrix<T> & solVector, gsMultiPatch<T> & mp) const =0;

    /// Construct deformed shell geometry from computed solution vector \a solVector and returns a multipatch
    virtual gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const  = 0;

    /// Construct deformed shell geometry from computed solution vector \a solVector and returns the result in \a deformed
    virtual void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const = 0;

    /// Construct displacement field from computed solution vector \a solVector and returns a multipatch
    virtual gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const = 0;

    /// Construct displacement field from computed solution vector \a solVector and returns the result in \a deformed
    virtual void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const = 0;

    /// Reconstruct the solution vector based on the currently stored boundary conditions (thus the mapper).
    virtual gsMatrix<T> fullSolutionVector(const gsMatrix<T> & vector) const = 0;

    /// Reconstruct the solution vector based on the currently stored boundary conditions (thus the mapper).
    virtual gsVector<T> constructSolutionVector(const gsMultiPatch<T> & deformed) const = 0;

    /// Construct Cauchy stress tensor for visualization
    virtual void constructStress(const gsFunctionSet<T> & deformed,gsPiecewiseFunction<T> & result,stress_type::type type) = 0;

    /// Compute the principal stretches in \a points given a \a deformed geometry. Optionally, the stretches can be computed on through-thickness coordinate \a z
    virtual gsMatrix<T> computePrincipalStretches(const gsMatrix<T> & points, const gsFunctionSet<T> & deformed, const T z=0) = 0;

    /// Compute the principal stretches in \a points given a \a deformed geometry. Optionally, the stretches can be computed on through-thickness coordinate \a z
    virtual gsMatrix<T> computePrincipalStresses(const gsMatrix<T> & points, const gsFunctionSet<T> & deformed, const T z=0) = 0;

    /// Returns the gsDofMapper
    virtual   gsDofMapper getMapper() = 0;

    /// Returns the C1 error over the interface
    virtual T interfaceErrorC0(const gsFunctionSet<T> & deformed) = 0;
    virtual T interfaceErrorC0(const gsFunctionSet<T> & deformed, const ifContainer & iFaces) = 0;

    /// Returns the G1 error over the interface
    virtual T interfaceErrorG1(const gsFunctionSet<T> & deformed) = 0;
    virtual T interfaceErrorG1(const gsFunctionSet<T> & deformed, const ifContainer & iFaces) = 0;

    /// Returns the normal vector error over the interface
    virtual T interfaceErrorNormal(const gsFunctionSet<T> & deformed) = 0;
    virtual T interfaceErrorNormal(const gsFunctionSet<T> & deformed, const ifContainer & iFaces) = 0;

    /// Returns the Gaussian curvature error over the interface
    virtual T interfaceErrorGaussCurvature(const gsFunctionSet<T> & deformed) = 0;
    virtual T interfaceErrorGaussCurvature(const gsFunctionSet<T> & deformed, const ifContainer & iFaces) = 0;

    /// Returns the mean curvature error over the interface
    virtual T interfaceErrorMeanCurvature(const gsFunctionSet<T> & deformed) = 0;
    virtual T interfaceErrorMeanCurvature(const gsFunctionSet<T> & deformed, const ifContainer & iFaces) = 0;

    /// Projects function \a fun on the basis and geometry stored in the class and returns the coefficients in \a result
    virtual void projectL2_into(const gsFunction<T> &fun, gsMatrix<T> & result) = 0;

    /// Projects function \a fun on the basis and geometry stored in the class and returns a multipatch in \a result
    virtual void projectL2_into(const gsFunction<T> &fun, gsMultiPatch<T> & result) = 0;

    /// Projects function \a fun on the basis and geometry stored in the class and returns the coefficients as a matrix
    virtual gsMatrix<T> projectL2(const gsFunction<T> &fun) = 0;

    virtual void plotSolution(std::string string, const gsMatrix<T> & solVector) = 0;

    virtual T deformationNorm(const gsMultiPatch<T> & deformed, const gsMultiPatch<T> & original) = 0;

    virtual gsDofMapper getMapper() const = 0;


};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsThinShellAssemblerBase
   */
  void pybind11_init_gsThinShellAssemblerBase(pybind11::module &m);
  void pybind11_enum_gsThinShellAssemblerStatus(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssembler.hpp)
#endif
