/** @file gsThinShellAssembler.h

    @brief Provides system matrices for the elasticity problem on thin shells.

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

/** @brief Assembles system matrices for thin shell linear and nonlinear elasticity problems.

    \tparam T coefficient type

    \ingroup gsKLShell
*/

// Desired template parameters


template <short_t d, class T, bool bending>
class gsThinShellAssembler
{
public:
    // typedef gsExprAssembler<T> Base;
    // typedef std::vector<std::pair<patchSide,int> > clamped_t;
public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] basis
    \param[in] bconditions is a gsBoundaryConditions object describing the boundary conditions
    \param[in] surface_force is a gsFunction object describing the 3D surface force
                and/or 2*\em thickness *body_force

    \ingroup gsKLShell
*/

    /*
        Templates arguments:
        - geometric dimension       int     DIM
        - bending term              bool    BEN
        - foundation                bool    FOU --> should be an option?
        - follower pressure         bool    PRE --> should be an option?

    */
    /// Default empty constructor
    gsThinShellAssembler() { }

    gsThinShellAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & basis,
                        const gsBoundaryConditions<T> & bconditions,
                        const gsFunction<T> & surface_force,
                        const gsMaterialMatrix<T> & materialmatrix);

    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}

    //--------------------- PROBLEM FORMULATION-------------------------------//
    void setPointLoads(const gsPointLoads<T> & pLoads){ m_pLoads = pLoads; }
    void setFoundation(const gsFunction<T> & foundation) { m_foundFun = &foundation; m_foundInd = true; }
    void setPressure(const gsFunction<T> & pressure) { m_pressFun = &pressure; m_pressInd = true; }

    void setPlate()
    {
        m_type = 2;
        m_dim = 2;
        this->initialize();
        // GISMO_ASSERT(m_patches.geoDim()==static_cast<index_t>(m_dim),"Domain dimension does not match element dimension? Domain dimension: "<<m_patches.geoDim()<<"; element dimension: "<<m_dim);
    }
    void setMembrane()
    {
        m_type = 1;
    }

    void updateBCs(const gsBoundaryConditions<T> & bconditions)
    {
        m_bcs = bconditions;
        space m_space = m_assembler.trialSpace(0);
        this->assembleDirichlet();

        m_ddofs = m_space.fixedPart();
        m_mapper = m_space.mapper();
    }
    void homogenizeDirichlet();

    index_t numDofs() const {return m_assembler.numDofs();}

    // void setFoundation(const T & foundation) { m_foundFun = gsConstantFunction<T>(foundation,2); m_foundInd = true; }
    // NOTE: can we improve the m_foundInd indicator? e.g. by checking if m_foundation exists?


    //--------------------- SYSTEM ASSEMBLY ----------------------------------//

    /// @brief Assembles the stiffness matrix and the RHS for the LINEAR ELASTICITY
    /// set *assembleMatrix* to false to only assemble the RHS;

    // template DIM BEN
    void assemble();
    // template DIM BEN FOU PRE
    void assembleShell();
    // template DIM BEN FOU PRE
    void assembleMembrane();

    void assembleMass();
    void assembleFoundation();


    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for displacement formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    void assemble(const gsMultiPatch<T> & deformed,     bool Matrix = true);
    void assemble(const gsMatrix<T>     & solVector,    bool Matrix = true);


    // template DIM BEN FOU PRE
    void assembleMatrix(const gsMultiPatch<T>   & deformed  );
    // template DIM BEN FOU PRE
    void assembleMatrixShell(const gsMultiPatch<T>   & deformed  );
    // template DIM BEN FOU PRE
    void assembleMatrixMembrane(const gsMultiPatch<T>   & deformed  );
    // template DIM BEN FOU PRE
    void assembleMatrix(const gsMatrix<T>       & solVector );

    // template DIM BEN FOU PRE
    void assembleVector(const gsMultiPatch<T>   & deformed  );
    // template DIM BEN FOU PRE
    void assembleVectorShell(const gsMultiPatch<T>   & deformed  );
    // template DIM BEN FOU PRE
    void assembleVectorMembrane(const gsMultiPatch<T>   & deformed  );
    // template DIM BEN FOU PRE
    void assembleVector(const gsMatrix<T>       & solVector );

    // template DIM BEN
    gsMatrix<T> boundaryForceVector(const gsMultiPatch<T>   & deformed , patchSide& ps, int com );
    // template DIM BEN
    gsMatrix<T> boundaryForceVectorShell(const gsMultiPatch<T>   & deformed , patchSide& ps, int com );
    // template DIM BEN
    gsMatrix<T> boundaryForceVectorMembrane(const gsMultiPatch<T>   & deformed , patchSide& ps, int com );

    //--------------------- GEOMETRY ACCESS --------------------------------//
    const gsMultiPatch<T> & geometry()    const  {return m_patches;}
    const gsMultiPatch<T> & defGeometry() const  {return m_defpatches;}
    const gsMultiPatch<T> & strips()      const  {return m_strips;}
    const gsMultiPatch<T> & defStrips()   const  {return m_defstrips;}

    T getArea(const gsMultiPatch<T> & geometry);


    //--------------------- SYSTEM ACCESS ----------------------------------//
    const gsSparseMatrix<T> & matrix()  const   {return m_assembler.matrix();}
    // gsSparseMatrix<T> & matrix() {return const_cast <gsSparseMatrix<T> &>(m_assembler.matrix());}

    // const gsMatrix<T>       & rhs()     const {return m_assembler.rhs();}
    const gsMatrix<T>       & rhs()     const {return m_rhs.size()==0 ? m_assembler.rhs() : m_rhs;}

    //--------------------- SOLUTION CONSTRUCTION ----------------------------------//

    /// @brief Construct deformed shell geometry from computed solution vector
    virtual gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const;
    virtual void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const;
    void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;



    //--------------------- SPECIALS ----------------------------------//

    /// @brief Construct Cauchy stress tensor for visualization (only valid for linear elasticity)
    void constructStress(const gsMultiPatch<T> & deformed,
                               gsPiecewiseFunction<T> & result,
                               stress_type::type type);

    // gsField<T> constructStress(const gsMultiPatch<T> & deformed,
                               // stress_type::type type);

    gsMatrix<T> computePrincipalStretches(const gsMatrix<T> & points, const gsMultiPatch<T> & deformed, const T z=0);

    void projectL2_into(const gsFunction<T> &fun, gsMatrix<T> & result);
    void projectL2_into(const gsFunction<T> &fun, gsMultiPatch<T> & result);
    gsMatrix<T> projectL2(const gsFunction<T> &fun);

protected:
    void initialize();
    void defaultOptions();
    void getOptions() const;

    void defineComponents();

    void assembleNeumann();
    void assembleDirichlet();
    void assembleClamped();

    void applyLoads();


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
    gsMultiBasis<T> m_basis;
    gsBoundaryConditions<T> m_bcs;

    mutable gsMatrix<T> m_ddofs;

    gsMaterialMatrix<T> m_materialMat;

    const gsFunction<T> * m_forceFun;
    const gsFunction<T> * m_thickFun;
    const gsFunction<T> * m_foundFun;
    const gsFunction<T> * m_pressFun;
    typename gsFunction<T>::Ptr m_YoungsModulus;
    typename gsFunction<T>::Ptr m_PoissonsRatio;

    gsPointLoads<T>  m_pLoads;

    mutable gsMatrix<T> m_solvector;

    mutable size_t m_dim;


    gsMatrix<T> m_rhs;

    mutable gsOptionList m_options;

    mutable bool m_nl_loads;
    mutable bool m_foundInd;
    mutable bool m_pressInd;

    mutable index_t m_type; // shell_type

    /// @brief Specifies the material law to use
    struct nl_loads
    {
        enum type
        {
            off = false,
            on  = true,
        };
    };
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssembler.hpp)
#endif
