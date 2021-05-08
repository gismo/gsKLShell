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

#include <gsKLShell/gsMaterialMatrix.h>
#include <gsKLShell/gsThinShellFunctions.h>

#include <gsPde/gsPointLoads.h>


namespace gismo
{

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

    /**
     * @brief      Constructor for te shell assembler
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

    /// Default empty constructor
    gsThinShellAssembler() { }

    /// See \ref gsThinShellAssemblerBase for details
    gsOptionList & options() {return m_options;}

    /// See \ref gsThinShellAssemblerBase for details
    gsExprAssembler<T> assembler() {return m_assembler; }

    /// See \ref gsThinShellAssemblerBase for details
    void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

    //--------------------- PROBLEM FORMULATION-------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    void setPointLoads(const gsPointLoads<T> & pLoads){ m_pLoads = pLoads; }

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

    ////////////////////////////////////////////////////////////////////////////
    //--------------------- SYSTEM ASSEMBLY ----------------------------------//
    ////////////////////////////////////////////////////////////////////////////

    void assemble();

private:
    /// Specialisation of assemble() for surfaces (3D)
    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, void>::type assemble_impl();

    /// Specialisation of assemble() for planar geometries (2D)
    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), void>::type assemble_impl();

public:
    /// See \ref gsThinShellAssemblerBase for details
    void assembleMass(bool lumped = false);

    /// See \ref gsThinShellAssemblerBase for details
    void assembleFoundation();

    /// See \ref gsThinShellAssemblerBase for details
    void assemble(const gsMultiPatch<T> & deformed,     bool Matrix = true);

    /// See \ref gsThinShellAssemblerBase for details
    void assemble(const gsMatrix<T>     & solVector,    bool Matrix = true);

    /// See \ref gsThinShellAssemblerBase for details
    void assembleMatrix(const gsMultiPatch<T>   & deformed  );

    /// See \ref gsThinShellAssemblerBase for details
    void assembleMatrix(const gsMatrix<T>       & solVector );

private:
    /// Implementation of assembleMatrix for surfaces (3D)
    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, void>::type
    assembleMatrix_impl(const gsMultiPatch<T>   & deformed  );

    /// Implementation of assembleMatrix for planar geometries (2D)
    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    assembleMatrix_impl(const gsMultiPatch<T>   & deformed  );

public:
    /// See \ref gsThinShellAssemblerBase for details
    void assembleVector(const gsMultiPatch<T>   & deformed  );

    /// See \ref gsThinShellAssemblerBase for details
    void assembleVector(const gsMatrix<T>       & solVector );

private:
    /// Implementation of assembleVector for surfaces (3D)
    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, void>::type
    assembleVector_impl(const gsMultiPatch<T>   & deformed  );

    /// Implementation of assembleVector for planar geometries (2D)
    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    assembleVector_impl(const gsMultiPatch<T>   & deformed  );

public:
    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> boundaryForceVector(const gsMultiPatch<T>   & deformed , patchSide& ps, index_t com );

    gsMatrix<T> boundaryForce(const gsMultiPatch<T>   & deformed , patchSide& ps);


private:
    /// Implementation of the boundary force vector for surfaces (3D)
    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, gsMatrix<T> >::type
    boundaryForceVector_impl(const gsMultiPatch<T>   & deformed , patchSide& ps, index_t com );

    /// Implementation of the boundary force vector for planar geometries (2D)
    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), gsMatrix<T> >::type
    boundaryForceVector_impl(const gsMultiPatch<T>   & deformed , patchSide& ps, index_t com );

    /// Implementation of the boundary force vector for surfaces (3D)
    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, gsMatrix<T> >::type
    boundaryForce_impl(const gsMultiPatch<T>   & deformed , patchSide& ps);

    /// Implementation of the boundary force vector for planar geometries (2D)
    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), gsMatrix<T> >::type
    boundaryForce_impl(const gsMultiPatch<T>   & deformed , patchSide& ps);

public:

    //--------------------- GEOMETRY ACCESS --------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    const gsMultiPatch<T> & geometry()      const  {return m_patches;}

    /// See \ref gsThinShellAssemblerBase for details
    const gsMultiPatch<T> & defGeometry()   const  {return m_defpatches;}

    /// See \ref gsThinShellAssemblerBase for details
    T getArea(const gsMultiPatch<T> & geometry);

    //--------------------- MATERIAL ACCESS --------------------------------//
    gsMaterialMatrixBase<T> * material()    const  {return m_materialMat;}

    //--------------------- SYSTEM ACCESS ----------------------------------//
    const gsSparseMatrix<T> & matrix()      const   {return m_assembler.matrix();}
    // gsSparseMatrix<T> & matrix() {return const_cast <gsSparseMatrix<T> &>(m_assembler.matrix());}

    const gsMatrix<T>       & rhs()         const {return m_rhs.size()==0 ? m_assembler.rhs() : m_rhs;}
    // const gsMatrix<T>       & rhs()     const {return m_assembler.rhs();}

    //--------------------- SOLUTION CONSTRUCTION ----------------------------------//
    gsMultiPatch<T> constructMultiPatch(const gsMatrix<T> & solVector) const;

    /// See \ref gsThinShellAssemblerBase for details
    gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const;

    /// See \ref gsThinShellAssemblerBase for details
    void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    /// See \ref gsThinShellAssemblerBase for details
    gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const;

    /// See \ref gsThinShellAssemblerBase for details
    void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    /// See \ref gsThinShellAssemblerBase for details
    gsVector<T> constructSolutionVector(const gsMultiPatch<T> & deformed) const;

    //--------------------- SPECIALS ----------------------------------//
    /// See \ref gsThinShellAssemblerBase for details
    void constructStress(const gsMultiPatch<T> & deformed,
                               gsPiecewiseFunction<T> & result,
                               stress_type::type type);

    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> computePrincipalStretches(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed, const T z=0);

    /// See \ref gsThinShellAssemblerBase for details
    void projectL2_into(const gsFunction<T> &fun, gsMatrix<T> & result);

    /// See \ref gsThinShellAssemblerBase for details
    void projectL2_into(const gsFunction<T> &fun, gsMultiPatch<T> & result);

    /// See \ref gsThinShellAssemblerBase for details
    gsMatrix<T> projectL2(const gsFunction<T> &fun);

protected:
    /// Initializes the method
    void _initialize();
    void _defaultOptions();
    void _getOptions() const;

    void _assembleNeumann();
    void _assembleWeakBCs();
    void _assembleWeakBCs(const gsMultiPatch<T> & deformed);
    void _assembleDirichlet();

    void _applyLoads();

private:
    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, void>::type
    _assembleNeumann_impl();

    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    _assembleNeumann_impl();

    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, void>::type
    _assembleWeakBCs_impl();

    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    _assembleWeakBCs_impl();

    template<int _d, bool _bending>
    typename std::enable_if<_d==3 && _bending, void>::type
    _assembleWeakBCs_impl(const gsMultiPatch<T> & deformed);

    template<int _d, bool _bending>
    typename std::enable_if<!(_d==3 && _bending), void>::type
    _assembleWeakBCs_impl(const gsMultiPatch<T> & deformed);

protected:
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsSparseSolver<>::CGDiagonal m_solver;

    std::vector<gsDofMapper>  m_dofMappers;
    gsDofMapper m_mapper;

    gsExprAssembler<> m_assembler;
    gsExprEvaluator<> m_evaluator;

    gsMultiPatch<T> m_patches;
    gsMultiPatch<T> m_defpatches;
    mutable gsMultiBasis<T> m_basis;
    gsBoundaryConditions<T> m_bcs;

    mutable gsMatrix<T> m_ddofs;

    const gsFunction<T> * m_forceFun;
    const gsFunction<T> * m_thickFun;
    const gsFunction<T> * m_foundFun;
    const gsFunction<T> * m_pressFun;
    typename gsFunction<T>::Ptr m_YoungsModulus;
    typename gsFunction<T>::Ptr m_PoissonsRatio;

    mutable gsMaterialMatrixBase<T> * m_materialMat;

    gsPointLoads<T>  m_pLoads;

    mutable gsMatrix<T> m_solvector;

    gsMatrix<T> m_rhs;

    mutable gsOptionList m_options;

    mutable bool m_foundInd;
    mutable bool m_pressInd;

    mutable index_t m_type; // shell_type

    mutable T m_alpha_d,m_alpha_r; // shell_type

};

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

    /// Sets the Dirichlet BCs to zero
    virtual void homogenizeDirichlet() = 0;

    /// Returns the number of degrees of freedom in the assembler.
    virtual index_t numDofs() const  = 0;

    /// Assembles the linear system and corresponding right-hand side
    virtual void assemble() = 0;

    /// Assembles the mass matrix (including density and thickness!); if lumped=true, a lumped mass matrix will be constructed,
    virtual void assembleMass(bool lumped = false) = 0;

    /// Assembles the elastic foundation matrix
    virtual void assembleFoundation() = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix and the residual for an iteration of Newton's method
     *
     * @note       .rhs() returns the negative residual!
     *
     * @param[in]  deformed  The deformed multipatch
     * @param[in]  Matrix    True if the matrix should be assembled
     */
    virtual void assemble(const gsMultiPatch<T> & deformed,     bool Matrix = true) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix and the residual for an iteration of Newton's method
     *
     * @note       .rhs() returns the negative residual!
     *
     * @param[in]  deformed  The solution vector
     * @param[in]  Matrix    True if the matrix should be assembled
     */
    virtual void assemble(const gsMatrix<T>     & solVector,    bool Matrix = true) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix (nonlinear)
     *
     * @param[in]  deformed  The deformed geometry
     */
    virtual void assembleMatrix(const gsMultiPatch<T>   & deformed  ) = 0;

    /**
     * @brief      Assembles the tangential stiffness matrix (nonlinear)
     *
     * @param[in]  deformed  The solution vector
     */
    virtual void assembleMatrix(const gsMatrix<T>       & solVector ) = 0;

    /**
     * @brief      Assembles the residual vector
     *
     * @param[in]  deformed  The deformed geometry
     */
    virtual void assembleVector(const gsMultiPatch<T>   & deformed  ) = 0;

    /**
     * @brief      Assembles the residual vector
     *
     * @param[in]  deformed  The solution vector
     */
    virtual void assembleVector(const gsMatrix<T>       & solVector ) = 0;

    /**
     * @brief      Computes the force on a boundary
     *
     * This function is typically used when you want to know the load on aboundary on which a displacement is applied
     *
     * @param[in]  deformed  The deformed geometry
     * @param      ps        The patch side
     * @param[in]  com       The component
     *
     * @return     The loads on the control points. The sum is the total load on the boundary
     */
    virtual gsMatrix<T> boundaryForceVector(const gsMultiPatch<T>   & deformed , patchSide& ps, int com ) = 0;

    virtual gsMatrix<T> boundaryForce(const gsMultiPatch<T>   & deformed , patchSide& ps) = 0;

    /// Returns the undeformed geometry
    virtual const gsMultiPatch<T> & geometry()    const = 0;

    /// Returns the deformed geometry
    virtual const gsMultiPatch<T> & defGeometry() const = 0;

    /// Returns the material matrix used in the class
    virtual gsMaterialMatrixBase<T> * material()          const = 0;

    /// Returns the area of \a geometry
    virtual T getArea(const gsMultiPatch<T> & geometry) = 0;

    /// Returns a reference to the system matrix that is assembled
    virtual const gsSparseMatrix<T> & matrix()  const  = 0;

    /// Returns a reference to the right-hand side vector that is assembled
    virtual const gsMatrix<T>       & rhs()     const  = 0;

    /// Construct solution field from computed solution vector \a solVector and returns a multipatch
    virtual gsMultiPatch<T> constructMultiPatch(const gsMatrix<T> & solVector) const = 0;

    /// Construct deformed shell geometry from computed solution vector \a solVector and returns a multipatch
    virtual gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const  = 0;

    /// Construct deformed shell geometry from computed solution vector \a solVector and returns the result in \a deformed
    virtual void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const = 0;

    /// Construct displacement field from computed solution vector \a solVector and returns a multipatch
    virtual gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const = 0;

    /// Construct displacement field from computed solution vector \a solVector and returns the result in \a deformed
    virtual void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const = 0;

    /// Reconstruct the solution vector based on the currently stored boundary conditions (thus the mapper).
    virtual gsVector<T> constructSolutionVector(const gsMultiPatch<T> & deformed) const = 0;

    /// Construct Cauchy stress tensor for visualization
    virtual void constructStress(const gsMultiPatch<T> & deformed,gsPiecewiseFunction<T> & result,stress_type::type type) = 0;

    /// Compute the principal stretches in \a points given a \a deformed geometry. Optionally, the stretches can be computed on through-thickness coordinate \a z
    virtual gsMatrix<T> computePrincipalStretches(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed, const T z=0) = 0;


    /// Projects function \a fun on the basis and geometry stored in the class and returns the coefficients in \a result
    virtual void projectL2_into(const gsFunction<T> &fun, gsMatrix<T> & result) = 0;

    /// Projects function \a fun on the basis and geometry stored in the class and returns a multipatch in \a result
    virtual void projectL2_into(const gsFunction<T> &fun, gsMultiPatch<T> & result) = 0;

    /// Projects function \a fun on the basis and geometry stored in the class and returns the coefficients as a matrix
    virtual gsMatrix<T> projectL2(const gsFunction<T> &fun) = 0;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssembler.hpp)
#endif
