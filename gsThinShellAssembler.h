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

#include <gsThinShell2/gsMaterialMatrix.h>
#include <gsThinShell2/gsThinShellFunctions.h>

#include <gsPde/gsPointLoads.h>


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
                        const gsMaterialMatrix<T> & materialmatrix);

    // gsThinShellAssembler(const gsMultiPatch<T> & patches,
    //                     const gsMultiBasis<T> & basis,
    //                     const gsBoundaryConditions<T> & bconditions,
    //                     const gsFunction<T> & surface_force,
    //                     const gsFunction<T> & thickness,
    //                     const gsFunction<T> & YoungsModulus,
    //                     const gsFunction<T> & PoissonsRatio);


    /// @brief Returns the list of default options for assembly
    gsOptionList defaultOptions();

    //--------------------- PROBLEM FORMULATION-------------------------------//
    void setPointLoads(const gsPointLoads<T> & pLoads){ m_pLoads = pLoads; }

    //--------------------- SYSTEM ASSEMBLY ----------------------------------//

    /// @brief Assembles the stiffness matrix and the RHS for the LINEAR ELASTICITY
    /// set *assembleMatrix* to false to only assemble the RHS;
    void assemble();

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for displacement formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    void assemble(const gsMultiPatch<T> & deformed,     bool Matrix = true);
    void assemble(const gsMatrix<T>     & solVector,    bool Matrix = true);

    void assemble(int i);

    void assembleMatrix(const gsMultiPatch<T>   & deformed  );
    void assembleMatrix(const gsMatrix<T>       & solVector );

    void assembleVector(const gsMultiPatch<T>   & deformed  );
    void assembleVector(const gsMatrix<T>       & solVector );

    //--------------------- GEOMETRY ACCESS --------------------------------//
    const gsMultiPatch<T> & geometry()    const  {return m_patches;}
    const gsMultiPatch<T> & defGeometry() const  {return m_defpatches;}
    const gsMultiPatch<T> & strips()      const  {return m_strips;}
    const gsMultiPatch<T> & defStrips()   const  {return m_defstrips;}


    //--------------------- SYSTEM ACCESS ----------------------------------//
    const gsSparseMatrix<T> & matrix()  const   {return m_assembler.matrix();}
    // gsSparseMatrix<T> & matrix() {return const_cast <gsSparseMatrix<T> &>(m_assembler.matrix());}

    const gsMatrix<T>       & rhs()     const {return m_assembler.rhs();}

    //--------------------- SOLUTION CONSTRUCTION ----------------------------------//

    /// @brief Construct deformed shell geometry from computed solution vector
    gsMultiPatch<T> constructSolution(const gsMatrix<T> & solVector) const;
    void constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;

    gsMultiPatch<T> constructDisplacement(const gsMatrix<T> & solVector) const;
    void constructDisplacement(const gsMatrix<T> & solVector, gsMultiPatch<T> & deformed) const;



    //--------------------- SPECIALS ----------------------------------//

    /// @brief Construct Cauchy stress tensor for visualization (only valid for linear elasticity)
    void constructStress(const gsMultiPatch<T> & deformed,
                               gsPiecewiseFunction<T> & result,
                               stress_type::type type);


    /// @brief Check whether the displacement field is valid, i.e. J = det(F) > 0;
    /// return -1 if yes or a number of the first invalid patch
    index_t checkSolution(const gsMultiPatch<T> & solution) const;

    /// @brief Return minJ/maxJ
    T solutionJacRatio(const gsMultiPatch<T> & solution) const;

protected:
    typedef typename std::vector< gsMatrix<unsigned> > gsStripIndices; // index_t instead of unsigned


    void initialize();
    void defineComponents();

    void assembleNeumann();
    void assembleDirichlet();
    void assembleClamped();

    void createBendingStrips(gsMultiPatch<T> & mp, gsMultiPatch<T> & strips, gsStripIndices & indices);
    void assembleBendingStrips();


    void applyLoads();


protected:
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    std::vector<gsDofMapper>  m_dofMappers;

    gsExprAssembler<> m_assembler;
    gsExprEvaluator<> m_evaluator;

    gsStripIndices m_stripIndices, m_defstripIndices;
    gsMultiPatch<T> m_strips, m_defstrips;
    gsExprAssembler<T> m_stripAssembler;

    // geometryMap m_ori;
    // geometryMap m_def;
    // space m_space;
    //mutable solution m_solution;

    //expr::gsFeVariable<T> * m_force;
    // variable m_thick;
    // variable m_materialMat; // material matrix
    // variable m_m2; // matrix for multiplication of last entries of components.

    gsMultiPatch<T> m_patches;
    gsMultiPatch<T> m_defpatches;
    gsMultiBasis<T> m_basis;
    gsBoundaryConditions<T> m_bcs;

    gsMaterialMatrix<T> m_materialMat;

    const gsFunction<T> * m_forceFun;
    const gsFunction<T> * m_thickFun;
    typename gsFunction<T>::Ptr m_YoungsModulus;
    typename gsFunction<T>::Ptr m_PoissonsRatio;

    gsPointLoads<T>  m_pLoads;

    mutable gsMatrix<T> m_solvector;

    // /*
    //     Make type aliasses for function expressions
    // */
    // template <typename T1, typename T2, typename T3 > using var2_t = gismo::expr::var2_expr<T1,T2,T3>;
    // template <typename T1, typename T2, typename T3 > using flatdot_t = gismo::expr::flatdot_expr<T1,T2,T3 >;
    // template <typename T1, typename T2, typename T3 > using flatdot2_t= gismo::expr::flatdot2_expr<T1,T2,T3 >;
    // template <typename T1, typename T2  > using mult_t = gismo::expr::mult_expr<T1,T2,false >;
    // template <typename T1, typename T2  > using divide_t= gismo::expr::divide_expr<T1,T2 >;
    // template <typename T1, typename T2  > using add_t  = gismo::expr::add_expr<T1,T2>;
    // template <typename T1, typename T2  > using sub_t  = gismo::expr::sub_expr<T1,T2>;
    // template <typename T1, typename T2  > using der2d_t= gismo::expr::deriv2dot_expr<T1,T2>;

    // template <typename T1> using jacG_t      = gismo::expr::jacG_expr<T1>;
    // template <typename T1> using jac_t       = gismo::expr::jac_expr<T1>;
    // template <typename T1> using sn_t        = gismo::expr::normal_expr<T1>;
    // template <typename T1> using var1_t      = gismo::expr::var1_expr<T1>;
    // template <typename T1> using der2_t      = gismo::expr::deriv2_expr<T1>;
    // template <typename T1> using normalized_t= gismo::expr::normalized_expr<T1>;
    // template <typename T1> using symmetrize_t= gismo::expr::symmetrize_expr<T1>;
    // template <typename T1> using flat_t      = gismo::expr::flat_expr<T1>;
    // template <typename T1> using tr_t        = gismo::expr::tr_expr<T1>;
    // template <typename T1> using u_t         = gismo::expr::gsFeSpace<T1>;
    // template <typename T1> using G_t         = gismo::expr::gsGeometryMap<T1>;
    // template <typename T1> using var_t       = gismo::expr::gsFeVariable<T1>;
    // template <typename T1> using val_t       = gismo::expr::value_expr<T1>;
    // template <typename T1> using reshape_t   = gismo::expr::reshape_expr<T1>;

    // template <typename T1> using Em_t  =
    // mult_t
    // < T,
    //     sub_t
    //     <
    //         flat_t
    //         <
    //             mult_t< tr_t< jacG_t<T1> >, jacG_t<T1> >
    //         >
    //         ,
    //         flat_t
    //         <
    //             mult_t< tr_t< jacG_t<T1> >, jacG_t<T1> >
    //         >
    //     >
    // >;

    // template <typename T1> using Sm_t  =
    // mult_t
    // <
    //     Em_t<T1>
    //     ,
    //     reshape_t< var_t <T1> >
    // >;

    // template <typename T1> using N_t  =
    // mult_t< val_t<var_t<T1>>, Sm_t<T1> >;

    // template <typename T1> using Em_der_t  =
    // flat_t
    // <
    //     mult_t< tr_t< jacG_t<T1> >, jac_t<u_t<T1>> >
    // >;

    // template <typename T1> using Sm_der_t  =
    // mult_t
    // <
    //     Em_der_t<T1>
    //     ,
    //     reshape_t< var_t <T1> >
    // >;

    // template <typename T1> using N_der_t  =
    // mult_t< val_t<var_t<T1>>, Sm_der_t<T1> >;

    // template <typename T1> using E_m_der2_t  =
    // flatdot_t
    // <
    //     jac_t< u_t<T1> >,
    //     tr_t< jac_t< u_t<T1> > >,
    //     N_t<T1>
    // >;

    // template <typename T1> using Ef_t  =
    // mult_t
    // <
    //     sub_t
    //     <
    //         der2d_t<G_t<T1>, tr_t< normalized_t< sn_t<T1> > > >
    //         ,
    //         der2d_t<G_t<T1>, tr_t< normalized_t< sn_t<T1> > > >
    //     >
    //     ,
    //     reshape_t< var_t<T1> >
    // >;

    // template <typename T1> using Sf_t  =
    // mult_t
    // <
    //     Ef_t<T1>
    //     ,
    //     reshape_t< var_t <T1> >
    // >;

    // template <typename T1> using M_t  =
    // mult_t
    // <
    //     divide_t
    //     <
    //         mult_t
    //         <
    //             mult_t
    //             <
    //                 val_t<var_t<T1>>
    //                 ,
    //                 val_t<var_t<T1>>
    //             >
    //             ,
    //             val_t<var_t<T1>>
    //         >
    //     ,
    //     T1
    //     >
    //     ,
    //     Sf_t<T1>
    // >;

    // template <typename T1> using Ef_der_t  =
    // mult_t
    // <
    //     add_t
    //     <
    //         der2d_t<u_t<T1>, tr_t< normalized_t< sn_t<T1> > > >
    //         ,
    //         der2d_t<G_t<T1>, var1_t< u_t<T1> > >
    //     >
    //     ,
    //     reshape_t< var_t<T1> >
    // >;

    // template <typename T1> using Sf_der_t  =
    // mult_t
    // <
    //     Ef_der_t<T1>
    //     ,
    //     reshape_t< var_t <T1> >
    // >;

    // template <typename T1> using M_der_t  =
    // mult_t
    // <
    //     divide_t
    //     <
    //         mult_t
    //         <
    //             mult_t
    //             <
    //                 val_t<var_t<T1>>
    //                 ,
    //                 val_t<var_t<T1>>
    //             >
    //             ,
    //             val_t<var_t<T1>>
    //         >
    //     ,
    //     T1
    //     >
    //     ,
    //     Sf_der_t<T1>
    // >;


    // template <typename T1> using E_f_der2_t  =
    // add_t
    // <
    //     symmetrize_t
    //     <
    //         flatdot2_t
    //         <
    //             der2_t< u_t<T1> >,
    //             tr_t< var1_t<u_t<T1> > >,
    //             M_t<T1>
    //         >
    //     >
    //     ,
    //     var2_t< u_t<T1>, u_t<T1>, M_t<T1> >
    // >;

    // template <typename T1> using force_t = var_t<T1>;

    // Em_t<T>         m_Em;
    // // Sm_t<T>         m_Sm;
    // N_t<T>          m_N;
    // Em_der_t<T>     m_Em_der;
    // // Sm_der_t<T>     m_Sm_der;
    // N_der_t<T>      m_N_der;
    // E_m_der2_t<T>    m_Em_der2;

    // Ef_t<T>         m_Ef;
    // // Sf_t<T>         m_Sf;
    // M_t<T>          m_M;
    // Ef_der_t<T>     m_Ef_der;
    // Sf_der_t<T>     m_Sf_der;
    // M_der_t<T>      m_M_der;
    // E_f_der2_t<T>    m_Ef_der2;

    // force_t<T>      m_ff;

};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThinShellAssembler.hpp)
#endif
