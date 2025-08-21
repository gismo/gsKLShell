/** @file shapeOpt_strainenergyC.cpp

    @brief Structural optimization problem of
           a geometrically nonlinear Kirchhoff-Love shell
           with embedded ribs (design variables: curve control points)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Authors: C. Chianese, H.M. Verhelst
*/

#include <gismo.h>
#include <gsKLShell/gsKLShell.h>
#include <gsKLShell/src/gsEmbeddingUtils.h>
#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsGradientDescent.h>

#ifdef gsOptim_ENABLED
#include <gsOptim/gsOptim.h>
#endif

template<class T>
class gsShapeOptProblem : public gsOptProblem<T>
{
    using Base = gsOptProblem<T>;
public:
    gsShapeOptProblem(  gsThinShellAssemblerBase<T> *assembler,
                        const gsDofMapper           &mapper,
                        const gsMultiPatch<T>       &geom,
                        const gsMultiPatch<T>       &curve,
                        const gsSparseMatrix<T>     &K_shell,
                        const gsVector<T>           &F_s,
                        const gsVector<T>           &materialParam)
    :
    m_assembler(assembler),
    m_mapper(mapper),
    m_geom(geom),
    m_curve(curve),
    m_K_shell(K_shell),
    m_F_s(F_s),
    m_materialParam(materialParam)
    m_delta_s(0.00001)
    {
        m_numDesignVars = m_mapper.freeSize();   // number of design variables
        m_numDofs = m_assembler->numDofs();      // number of unconstrained control points * dimensionality
        m_curDesign = this->vectorUpdate(m_curve, m_mapper);
        gsDebug << "Current design variables: " << m_curDesign.transpose() << "\n";
        gsDebugVar(m_numDesignVars);

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);
        for (index_t p = 0; p != m_curve.nPatches(); ++p)
        {
            for (index_t i = 0; i != m_curve.patch(p).coefs().rows(); ++i)
            {
                const index_t glx = mapper.index(i,p,0);
                if (mapper.is_free_index(glu))
                {
                    // m_desLowerBounds[glu] = m_curDesign(glu,0) - 0.5;  // u-coordinate
                    // m_desUpperBounds[glu] = m_curDesign(glu,0) + 0.5;
                    m_desLowerBounds[glu] = 0.0;  // u-coordinate
                    m_desUpperBounds[glu] = 1.0;
                }
                
                const index_t glv = mapper.index(i,p,1);
                if (mapper.is_free_index(glv))
                {
                    // m_desLowerBounds[glv] = m_curDesign[glv] - 0.5;  // v-coordinate
                    // m_desUpperBounds[glv] = m_curDesign[glv] + 0.5;
                    m_desLowerBounds[glv] = 0.0;  // v-coordinate
                    m_desUpperBounds[glv] = 1.0;
                }
            }
        }

        // m_numConstraints = 0;
        // m_conJacRows.resize(m_numConstraints);
        // m_conJacCols.resize(m_numConstraints);
        // m_conLowerBounds.resize(m_numConstraints);
        // m_conUpperBounds.resize(m_numConstraints);
    }

    // void updateBounds(const gsAsConstVector<T> &m_curDesign)
    // {
    //     for (index_t i = 0; i < m_numDesignVars; ++i)
    //     {
    //         m_desLowerBounds[i] = std::max(m_curDesign(i, 0) - 0.15, T(-0.2));
    //         m_desUpperBounds[i] = std::min(m_curDesign(i, 0) + 0.15, T(1.2));
    //     }
    // }

    static gsVector<T> vectorUpdate(const gsMultiPatch<T> &curve, const gsDofMapper &mapper)
    {
        gsVector<T> result(mapper.freeSize());
        for (index_t p = 0; p != curve.nPatches(); ++p)
        {
            for (index_t i = 0; i != curve.patch(p).coefs().rows(); ++i)
            {
                for (index_t j = 0; j != curve.patch(p).coefs().cols(); ++j)
                {
                    const index_t gl = mapper.index(i,p,j);
                    if (mapper.is_free_index(gl))
                        result[gl] = curve.patch(p).coefs()(i,j);
                }
            }
        }
        return result;
    }

    static void geomUpdate(const gsAsConstVector<T> &u, gsMultiPatch<T> &curve, const gsDofMapper &mapper)
    {
        for (index_t p = 0; p != curve.nPatches(); ++p)
        {
            for (index_t i = 0; i != curve.patch(p).coefs().rows(); ++i)
            {
                for (index_t j = 0; j != curve.patch(p).coefs().cols(); ++j)
                {
                    const index_t gl = mapper.index(i, p, j);
                    if (mapper.is_free_index(gl))
                        curve.patch(p).coefs()(i,j) = u[gl];
                }
            }
        }
    }

    static void geomUpdate(const gsVector<T> &u, gsMultiPatch<T> &curve, const gsDofMapper &mapper)
    {
        gsAsConstVector<T> tmp(u.data(),u.size());
        gsShapeOptProblem<T>::geomUpdate(tmp,curve,mapper);
    }

    gsVector<T> solveStateEquation(gsMultiPatch<T> &curveAn, gsVector<T> &F_s, gsSparseMatrix<T> &K_shell, gsSparseMatrix<T> &K_s,
                                   const std::vector<gsMatrix<T>> &allQuPointsCurve, const std::vector<gsVector<T>> &allQuWeights)
    {
        T EA = m_materialParam[0];        T EI_min = m_materialParam[1];
        T EI_max = m_materialParam[2];    T GI_p   = m_materialParam[3];

        ThinShellAssemblerStatus status_embedded = m_assembler->
            assembleLinearEmbeddedCurve(curveAn,EA,EI_min,EI_max,GI_p,allQuPointsCurve,allQuWeights);
        GISMO_ENSURE(status_embedded == ThinShellAssemblerStatus::Success, "Embedded beam assembly failed");
        K_s = K_shell +  m_assembler->matrix();

        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(K_s);
        gsVector<> u_s = solver.solve(F_s);
        return u_s;
    }

    T evalObj(const gsAsConstVector<T> &u) const override
    {
        gsDebug << "Computing objective at point " << u.transpose() << "\n";

        // Update embedded rib geometry from current design u
        gsMultiPatch<> tmpCurve = m_curve;
        this->geomUpdate(u,tmpCurve,m_mapper);

        // h-refine embedded rib for conforming quadrature
        gsMultiPatch<> tmpCurveAn = tmpCurve;
        std::vector<gsMatrix<T>> allquPointsCurve;
        std::vector<gsVector<T>> allquWeights;

        for (index_t p = 0; p < tmpCurveAn.nPatches(); ++p)
        {
            gsMatrix<T> quPoints;         gsVector<T> quWeights;
            embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn.patch(p),quPoints,quWeights,m_assembler);
            allquPointsCurve.push_back(quPoints);
            allquWeights.push_back(quWeights);
        }

        // Assemble the linear system for the current design of the rib-enforced shell
        gsSparseMatrix<T> K_s;
        gsVector<T> u_s = solveStateEquation(tmpCurveAn, m_F_s, m_K_shell, K_s, allquPointsCurve, allquWeights);

        // Return the objective function, i.e.strain energy, at current design
        T obj = 0.5 * u_s.transpose() * K_s * u_s;
        gsDebug << "Objective: " << obj << " at point " << u.transpose() << "\n";
        return obj;
    }

    //void gradObj_analytical_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const
    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const override
    {
        gsDebug << "Computing gradient at point " << u.transpose() << "\n";
        result.resize(m_numDesignVars);

        // Update embedded rib geometry from current design u
        gsMultiPatch<> tmpCurve = m_curve;
        this->geomUpdate(u,tmpCurve,m_mapper);

        // h-refine embedded rib for conforming quadrature
        gsMultiPatch<> tmpCurveAn = tmpCurve;
        std::vector<gsMatrix<T>> allquPointsCurve;
        std::vector<gsVector<T>> allquWeights;

        for (index_t p = 0; p < tmpCurve.nPatches(); ++p)
        {
            gsMatrix<T> quPoints;         gsVector<T> quWeights;
            embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn.patch(p),quPoints,quWeights,m_assembler);
            allquPointsCurve.push_back(quPoints);
            allquWeights.push_back(quWeights);
        }

        // Assemble the current design of the rib-enforced shell
        gsSparseMatrix<T> K_s;
        gsVector<T> u_s = solveStateEquation(tmpCurveAn, m_F_s, m_K_shell, K_s, allquPointsCurve, allquWeights);
        
        // Compute pseudo load matrix R*
        gsMatrix<> R_star(m_numDofs,m_numDesignVars);
        for (index_t p = 0; p != tmpCurve.nPatches(); ++p)
        {
            for (index_t i = 0; i != tmpCurve.patch(p).coefs().rows(); ++i)
            {
                for (index_t j = 0; j != tmpCurve.patch(p).coefs().cols(); ++j)
                {
                    index_t gl = m_mapper.index(i,p,j);
                    if (!m_mapper.is_free_index(gl)) continue;

                    // Perturb optimization rib topology
                    gsMultiPatch<> tmpCurve_splusds = tmpCurve;
                    tmpCurve_splusds.patch(p).coefs()(i,j) += m_delta_s;
                    gsMultiPatch<> tmpCurveAn_splusds = tmpCurve_splusds;

                    // Refine to match analysis resolution
                    allquPointsCurve.clear();
                    allquWeights.clear();

                    for (index_t p = 0; p < tmpCurveAn_splusds.nPatches(); ++p)
                    {
                        gsMatrix<T> quPoints;         gsVector<T> quWeights;
                        embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn_splusds.patch(p),quPoints,quWeights,m_assembler);
                        allquPointsCurve.push_back(quPoints);
                        allquWeights.push_back(quWeights);
                    }

                    // Assemble perturbed system
                    T EA = m_materialParam[0];        T EI_min = m_materialParam[1];
                    T EI_max = m_materialParam[2];    T GI_p   = m_materialParam[3];

                    ThinShellAssemblerStatus status_embedded = m_assembler->
                        assembleLinearEmbeddedCurve(tmpCurveAn_splusds,EA,EI_min,EI_max,GI_p,allQuPointsCurve,allQuWeights);
                    GISMO_ENSURE(status_embedded == ThinShellAssemblerStatus::Success, "Embedded beam assembly failed");
                    K_splusds = m_K_shell +  m_assembler->matrix();

                    R_star.col(gl) =  (-0.5 * (K_splusds - K_s) * u_s)/ m_delta_s; //since linear rhs does not depend on embedded ribs
                }
            }
        }
        // Return sensitivity vector df/ds
        result = u_s.transpose() * R_star;
        gsDebug << "Sensitivity vector: " << result.transpose() << "at point " << u.transpose() << "\n";
        //gsDebugVar(result);
    }

    void gradObj_FDM_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const
    {
        this->gradObj_into(u, result);
    }

protected:
    gsThinShellAssemblerBase<T>    *m_assembler;
    const gsDofMapper              &m_mapper;
    const gsMultiPatch<T>          &m_geom;
    const gsMultiPatch<T>          &m_curve;
    const gsSparseMatrix<T>        &m_K_shell;
    const gsVector<T>              &m_F_shell;
    index_t                         m_numDofs;
    T                               m_EA;
    T                               m_EI_min;
    T                               m_EI_max;
    T                               m_GI_p;
    T                               m_delta_s;
    using Base::m_numDesignVars;
    using Base::m_curDesign;
    using Base::m_desLowerBounds;
    using Base::m_desUpperBounds;
    // using Base::m_numConstraints;
    // using Base::m_conLowerBounds;
    // using Base::m_conUpperBounds;
    // using Base::m_conJacRows;
    // using Base::m_conJacCols;
};

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    index_t numRefine  = 0;
    index_t numRefineOpt = 0;
    gsCmdLine cmd("Rib topology optimization based on strain energy.");
    cmd.addInt( "r", "numRefine", "Number of uniform h-refinement steps to perform on shell",  numRefine );
    cmd.addInt( "R", "numRefineOpt", "Number of uniform h-refinement steps to perform on ribs",  numRefineOpt );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Shell reference geometry for analysis]
    gsMultiPatch<> mp_surf;
    gsKnotVector<> kv_u(0,1,2,3);   //start,end, number of interior knots, mult of end knots
    gsKnotVector<> kv_v(0,1,2,3);
    gsTensorBSplineBasis<2, real_t> basis_s(kv_u, kv_v);
    gsMatrix<> coefs_s (basis_s.size(), 3);
    coefs_s << 0, 0, 0,
               0.166, 0, 0.130,
               0.5, 0, 0.215,
               0.833, 0, 0.130,
               1, 0, 0,
               0, 0.166, 0.130,
               0.166, 0.166, 0.260,
               0.5, 0.166, 0.347,
               0.833, 0.166, 0.260,
               1, 0.166, 0.130,
               0, 0.5, 0.215,
               0.166, 0.5, 0.347,
               0.5, 0.5, 0.433,
               0.833, 0.5, 0.347,
               1, 0.5, 0.215,
               0, 0.833, 0.130,
               0.166, 0.833, 0.260,
               0.5, 0.833, 0.347,
               0.833, 0.833, 0.260,
               1, 0.833, 0.130,
               0, 1, 0,
               0.166, 1, 0.130,
               0.5, 1, 0.215,
               0.833, 1, 0.130,
               1, 1, 0;
    gsTensorBSpline<2, real_t>  surf(basis_s, coefs_s); //original surface
    mp_surf.addPatch(surf);
    mp_surf.addAutoBoundaries();
    mp_surf.degreeElevate(1,-1); //set degree to 2 in both directions
    for (int r = 0; r < numRefine; ++r)
         mp_surf.uniformRefine();
    gsMultiBasis<> mbasis_surf(mp_surf);

    // surf is already a reference, not a pointer
    // gsGeometry<real_t> &surfgeo = mp_surf.patch(0);
    // gsTensorBSpline<2, real_t>* surf = dynamic_cast< gsTensorBSpline<2, real_t>* >(&surfgeo);

    gsInfo << "\nShell reference geometry\n";
    gsInfo << "Patches: "<< mp_surf.nPatches() <<", degree: "<< mbasis_surf.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surf.basis(0)<<"\n";
    gsWriteParaview(mp_surf, "initialDesign", 1000, true, true);
    //! [Shell reference geometry for analysis]

    //! [Embedded ribs]
    gsKnotVector<real_t> kv_c = surf.knots(0);
    gsBSplineBasis<> basis_c(kv_c);

    gsEigen::ArrayXXd cpvec (surf.knots(0).size() - surf.degree(0) - 1, 1);
    cpvec = (surf.coefs().block(0,0,cpvec.rows(),1))/r;
    auto cpvec_flipped = cpvec.reverse();
    gsMatrix<real_t> coef_c1(basis_c.size(), surf.parDim());
    coef_c1.col(0) = cpvec;
    coef_c1.col(1) = cpvec_flipped;
    gsMatrix<real_t> coef_c2(basis_c.size(), surf.parDim());
    coef_c2.col(0) = cpvec;
    coef_c2.col(1) = cpvec;

    gsBSpline<> ribA(basis_c, coef_c1);
    gsBSpline<> ribB(basis_c, coef_c2);

    gsMultiPatch<> mp_curve;
    mp_curve.addPatch(ribA);
    mp_curve.addPatch(ribB);

    gsMatrix<real_t> coef_LLDPE(basis_c.size(), surf.parDim());
    coef_LLDPE.col(0) = cpvec;
    coef_LLDPE.col(1).setOnes();
    gsBSpline<> LLDPE_top(basis_c, coef_LLDPE);    // top buoyant breakwater
    mp_curve.addPatch(LLDPE_top);

    coef_LLDPE.col(1).setZero();
    gsBSpline<> LLDPE_bottom(basis_c, coef_LLDPE); // bottom buoyant breakwater
    mp_curve.addPatch(LLDPE_bottom);

    coef_LLDPE.col(0).setZero();
    coef_LLDPE.col(1) = cpvec;
    gsBSpline<> LLDPE_left(basis_c, coef_LLDPE);   // left buoyant breakwater
    mp_curve.addPatch(LLDPE_left);

    coef_LLDPE.col(0).setOnes();
    coef_LLDPE.col(1) = cpvec;
    gsBSpline<> LLDPE_right(basis_c, coef_LLDPE);  // right buoyant breakwater
    mp_curve.addPatch(LLDPE_right);

    for (int r = 0; r < numRefineOpt; ++r)
    {
        for (index_t i = 0; i < mp_curve.nPatches(); ++i)
             mp_curve.patch(i).uniformRefine();
    }

    gsMultiBasis<> mbasis_curve(mp_curve);
    gsInfo << "\nEmbedded rib for optimization\n";
    gsInfo << "Patches: "<< mp_curve.nPatches() <<", degree: "<< mbasis_curve.minCwiseDegree() <<"\n";
    gsInfo << mbasis_curve.basis(0)<<"\n";
    gsWriteParaview(mp_curve,  "ribs",  1000, true,  true);
    //! [Embedded ribs]

    //! [Material properties of shell and embedded ribs]
    real_t E_modulus = 1.5e9; // [Pa] HDPE
    real_t PoissonRatio = 0.45;
    real_t density = 950; // [kg/m^3]
    real_t thickness = 5.0e-3; // [m]

    real_t E_modulus_b = 1.5e9; // [Pa] HDPE
    real_t PoissonRatio_b = 0.45;
    real_t thickness = 0.10; // [m] 3e-3
    real_t height = 0.15;   // [m] 30e-3
    real_t G_modulus_b = 0.5 * E_modulus_b / (1 + PoissonRatio_b);
    real_t EA = E_modulus_b * (height * thickness);               //axial rigidity 
    real_t EI_min = E_modulus_b * (height * pow(thickness,3))/12; //minimum flexural rigidity
    real_t EI_max = E_modulus_rib * (thickness_rib * pow(height_rib,3))/12; //maximum flexural rigidity
    real_t GI_p = G_modulus_b/E_modulus_b * (EI_min + EI_max);  //torsional rigidity

    gsVector<real_t> materialParam(4);
    materialParam << EA, EI_min, EI_max, GI_p;
    //! [Material properties of shell and embedded ribs]

    //! [Make material functions: linear isotropic model]
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    //! [Make material functions: linear isotropic model]

    //! [Set boundary conditions and loads]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp_surf);
    bc.addCornerValue(boundary::southwest, 0.0, 0, 0, -1);
    bc.addCornerValue(boundary::southeast, 0.0, 0, 0, -1);
    bc.addCornerValue(boundary::northwest, 0.0, 0, 0, -1);
    bc.addCornerValue(boundary::northeast, 0.0, 0, 0, -1);

    // Pressure resulting from shell gravity loading, weight of PV layers and buoyancy on shell underside
    gsVector<> tmp(3);
    tmp << 0,0,-26.8794; //[N/m^2]
    gsConstantFunction<> force(tmp,3);

    //Buoyant line loads on shell edges
    gsVector<> buoyancy(3);
    buoyancy << 0,0,78.933; // [N/m] 
    gsConstantFunction<> neuData(buoyancy,3);
    bc.addCondition(0,boundary::west, condition_type::neumann,  &neuData);
    bc.addCondition(0,boundary::east, condition_type::neumann,  &neuData);
    bc.addCondition(0,boundary::north, condition_type::neumann, &neuData);
    bc.addCondition(0,boundary::south, condition_type::neumann, &neuData);
    //! [Set boundary conditions and loads]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    materialMatrix = getMaterialMatrix<3,real_t>(mp_surf,t,parameters,rho,options);
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t,true>(mp_surf,mbasis_surf,bc,force,materialMatrix);
    //! [Make assembler]

     //! [Assemble shell linear part]
    ThinShellAssemblerStatus status = assembler->assemble();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Shell assembly failed");
    gsSparseMatrix<> K_shell = assembler->matrix();
    gsVector<> F_s = assembler->rhs();
    //! [Assemble shell linear part]

    //! [Collect control points on rib ends]
    gsDofMapper mapper(mbasis_curve,mp_curve.geoDim());
    for ( size_t p = 0; p < mbasis_curve.nBases(); ++p)
    {
        for (short_t c = 1; c!=3; c++)
        {
            boxCorner corner(c);
            index_t idx = mbasis_curve.basis(p).functionAtCorner(corner);
            for (short_t d = 0; d!=mp_curve.geoDim(); ++d)
                mapper.eliminateDof(idx,p,d);
        }
    } 
    mapper.finalize();
    //! [Collect control points on rib ends]

    //! [Optimizer setup]
    gsShapeOptProblem<real_t> problem(assembler,mapper,mp_surf,mp_curve,K_shell,F_s,materialParam);

    gsOptimizer<real_t> *optimizer;
#ifdef gsOptim_ENABLED
    optimizer = new gsOptim<real_t>::LBFGS(&problem);
#else
    optimizer = new gsGradientDescent<>(&problem);
    optimizer->options().setReal("MinGradientLength",1e-9);
    optimizer->options().setReal("MinStepLength",1e-9);
    optimizer->options().setReal("MaxStepLength", 1e-2);
#endif
    optimizer->options().setInt("MaxIterations",100);
    optimizer->options().setInt("Verbose",1);
    //optimizer->options().setReal("GradErrTol",1e-8);
    //! [Optimizer setup]

    gsVector<> reshaped = gsShapeOptProblem<real_t>::vectorUpdate(mp_curve,mapper);
    gsAsConstVector<> initialDesign(reshaped.data(), reshaped.size());
    // gsDebugVar(problem.evalObj(initialDesign));
    // gsDebugVar(initialDesign.transpose());

    //! [Solve]
    // Start optimization
    optimizer->solve(initialDesign);
    //! [Solve]

    // Get the optimized design
    gsVector<> optimizedDesign = optimizer->currentDesign();
    gsShapeOptProblem<real_t>::geomUpdate(optimizedDesign,mp_curve,mapper);

    // Plot optimized design
    gsWrite(mp_surf, "ShellShape"); //.xml file of shell geometry
    gsWrite(mp_curve, "OptimalTopology"); //.xml file of rib geometry
    gsWriteParaview(mp_curve, "OptimizedDesign", 1000, true, false);

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;
}