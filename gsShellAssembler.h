/** @file gsThinShellAssembler.h

    @brief Provides system matrices for the elasticity problem on thin shells.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris
*/

#pragma once

#include <gsThinShell2/gsThinShellUtils.h>

namespace gismo
{

/** @brief Assembles system matrices for thin shell linear and nonlinear elasticity problems.

    \tparam T coefficient type

    \ingroup gsThinShell
*/
template <class T>
class gsThinShellAssembler
{
public:
    // typedef gsExprAssembler<T> Base;
    // typedef std::vector<std::pair<patchSide,int> > clamped_t;
public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] thickness is half of the shell thickness
    \param[in] E_modulus is the E-modulus of the chosen material
    \param[in] poissons_ratio is the poisson ratio of the chosen material
    \param[in] bconditions is a gsBoundaryConditions object describing the boundary conditions
    \param[in] surface_force is a gsFunction object describing the 3D surface force
                and/or 2*\em thickness *body_force
    \param[in] clamped contains pairs of patchside objects and x,y or z direction aimed at getting
                zero derivative for that direction
    \param[in] pLoads is a gsPointLoads object describing the points and corresponding loads

    \ingroup Assembler
*/
    gsThinShellAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & basis,
                        const gsBoundaryConditions<T> & bconditions,
                        const gsFunction<T> & surface_force,
                        const gsFunction<T> & thickness,
                        T YoungsModulus,
                        T PoissonsRatio);

    gsThinShellAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & basis,
                        const gsBoundaryConditions<T> & bconditions,
                        const gsFunction<T> & surface_force,
                        const gsFunction<T> & thickness,
                        const gsFunction<T> & YoungsModulus,
                        const gsFunction<T> & PoissonsRatio);



    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    //--------------------- SYSTEM ASSEMBLY ----------------------------------//

    /// @brief Assembles the stiffness matrix and the RHS for the LINEAR ELASTICITY
    /// set *assembleMatrix* to false to only assemble the RHS;
    virtual void assemble();

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for displacement formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    virtual void assemble(const gsMultiPatch<T> & deformed,     bool assembleMatrix = true);
    virtual void assemble(const gsMatrix<T>     & solVector,    bool assembleMatrix = true);

    virtual void assembleMatrix(const gsMultiPatch<T>   & deformed  );
    virtual void assembleMatrix(const gsMatrix<T>       & solVector );

    virtual void assembleVector(const gsMultiPatch<T>   & deformed  );
    virtual void assembleVector(const gsMatrix<T>       & solVector );

    //--------------------- SYSTEM ACCESS ----------------------------------//
    virtual gsSparseMatrix<T>   matrix()   {return m_assembler.matrix();}
    virtual gsVector<T>         rhs()      {return m_assembler.rhs();}

    //--------------------- SOLUTION CONSTRUCTION ----------------------------------//

    /// @brief Construct deformed shell geometry from computed solution vector
    virtual gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const;
    virtual void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    virtual gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const;
    virtual void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;



    //--------------------- SPECIALS ----------------------------------//

    // /// @brief Construct Cauchy stress tensor for visualization (only valid for linear elasticity)
    // void constructCauchyStresses(const gsMultiPatch<T> & displacement,
    //                              gsPiecewiseFunction<T> & result,
    //                              stress_type::type type = stress_type::von_mises) const;

    /// @brief Check whether the displacement field is valid, i.e. J = det(F) > 0;
    /// return -1 if yes or a number of the first invalid patch
    virtual index_t checkSolution(const gsMultiPatch<T> & solution) const;

    /// @brief Return minJ/maxJ
    virtual T solutionJacRatio(const gsMultiPatch<T> & solution) const;

protected:



protected:
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> m_assembler;
    gsExprEvaluator<> m_evaluator;

    geometryMap m_ori;
    geometryMap m_def;

    space m_space;

    gsMultiPatch<T> m_patches;
    gsMultiPatch<T> m_defpatches;
    gsMultiBasis<T> m_basis;
    gsBoundaryConditions<T> m_bcs;

    const gsFunction<T> * m_YoungsModulus;
    const gsFunction<T> * m_PoissonsRatio;
    const gsFunction<T> * m_forceFun;
    const gsFunction<T> * m_thickFun;
    variable m_force;
    variable m_thick;

    gsMatrix<T> m_solvector;
    solution m_solution;

    // material matrix
    variable m_materialmat;

    // matrix for multiplication of last entries of components.
    variable m_m2;

    /*
        Make type aliasses for function expressions
    */
        template <typename T1, typename T2, typename T3 > using var2_t = gismo::expr::var2_expr<T1,T2,T3>;
        template <typename T1, typename T2, typename T3 > using flatdot_t = gismo::expr::flatdot_expr<T1,T2,T3 >;
        template <typename T1, typename T2, typename T3 > using flatdot2_t= gismo::expr::flatdot2_expr<T1,T2,T3 >;
        template <typename T1, typename T2  > using mult_t = gismo::expr::mult_expr<T1,T2,false >;
        template <typename T1, typename T2  > using add_t  = gismo::expr::add_expr<T1,T2>;
        template <typename T1, typename T2  > using sub_t  = gismo::expr::sub_expr<T1,T2>;
        template <typename T1, typename T2  > using der2d_t= gismo::expr::deriv2dot_expr<T1,T2>;

        template <typename T1> using jacG_t      = gismo::expr::jacG_expr<T1>;
        template <typename T1> using jac_t       = gismo::expr::jac_expr<T1>;
        template <typename T1> using sn_t        = gismo::expr::normal_expr<T1>;
        template <typename T1> using var1_t      = gismo::expr::var1_expr<T1>;
        template <typename T1> using der2_t      = gismo::expr::deriv2_expr<T1>;
        template <typename T1> using normalized_t= gismo::expr::normalized_expr<T1>;
        template <typename T1> using symmetrize_t= gismo::expr::symmetrize_expr<T1>;
        template <typename T1> using flat_t      = gismo::expr::flat_expr<T1>;
        template <typename T1> using tr_t        = gismo::expr::tr_expr<T1>;
        template <typename T1> using u_t         = gismo::expr::gsFeSpace<T1>;
        template <typename T1> using G_t         = gismo::expr::gsGeometryMap<T1>;
        template <typename T1> using var_t       = gismo::expr::gsFeVariable<T1>;
        template <typename T1> using reshape_t   = gismo::expr::reshape_expr<T1>;

        template <typename T1> using Em_t  =
        mult_t
        < T1,
            sub_t
            <
                flat_t
                <
                    mult_t< tr_t< jacG_t<T1> >, jacG_t<T1> >
                >
                ,
                flat_t
                <
                    mult_t< tr_t< jacG_t<T1> >, jacG_t<T1> >
                >
            >
        >;

        template <typename T1> using Em_der_t  =
        flat_t
        <
            mult_t< tr_t< jacG_t<T1> >, jac_t<u_t<T1>> >
        >;

        template <typename T1> using Em_der2_t  =
        flatdot_t
        <
            jac_t< u_t<T1> >,
            tr_t< jac_t< u_t<T1> > >,
            mult_t< Em_t<T1>,reshape_t< var_t<T1> > >
        >;

        template <typename T1> using Ef_t  =
        mult_t
        <
            sub_t
            <
                der2d_t<G_t<T1>, tr_t< normalized_t< sn_t<T1> > > >
                ,
                der2d_t<G_t<T1>, tr_t< normalized_t< sn_t<T1> > > >
            >
            ,
            reshape_t< var_t<T1> >
        >;

        template <typename T1> using Ef_der_t  =
        mult_t
        <
            add_t
            <
                der2d_t<u_t<T1>, tr_t< normalized_t< sn_t<T1> > > >
                ,
                der2d_t<G_t<T1>, var1_t< u_t<T1> > >
            >
            ,
            reshape_t< var_t<T1> >
        >;

        template <typename T1> using Ef_der2_t  =
        add_t
        <
            symmetrize_t
            <
                flatdot2_t
                <
                    der2_t< u_t<T1> >,
                    tr_t< var1_t<u_t<T1> > >,
                    mult_t < Ef_t<T1>, reshape_t< var_t <T1> > >
                >
            >
            ,
            var2_t< u_t<T1>, u_t<T1>, mult_t < Ef_t<T1>, reshape_t< var_t <T1> > > >
        >;

    Em_t<T>         m_Em;
    Em_der_t<T>     m_Em_der;
    Em_der2_t<T>    m_Em_der2;
    Ef_t<T>         m_Ef;
    Ef_der_t<T>     m_Ef_der;
    Ef_der2_t<T>    m_Ef_der2;

};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssembler.hpp)
#endif
