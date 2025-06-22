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
                        const gsVector<T>           &F_shell,
                        const T                     &EA,
                        const T                     &EI_min,
                        const T                     &EI_max,
                        const T                     &GI_p)
    :
    m_assembler(assembler),
    m_mapper(mapper),
    m_geom(geom),
    m_curve(curve),
    m_K_shell(K_shell),
    m_F_shell(F_shell),
    m_EA(EA),
    m_EI_min(EI_min),
    m_EI_max(EI_max),
    m_GI_p(GI_p),
    m_delta_s(0.00001)
    {
        m_numDesignVars = m_mapper.freeSize();   // number of design variables
        m_numDofs = m_assembler->numDofs();      // number of unconstrained control points * dimensionality
        m_curDesign = this->vectorUpdate(m_curve, m_mapper);
        gsDebug << "Current design variables: " << m_curDesign.transpose() << "\n";
        gsDebugVar(m_numDesignVars);

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);
        for (index_t i = 0; i != m_numDesignVars; ++i)
        {
        //     // m_desLowerBounds[i] = std::max(m_curDesign(i,0) - 0.2, T(0.0));
        //     // m_desUpperBounds[i] = std::min(m_curDesign(i,0) + 0.2, T(1.0 - m_delta_s));
        //     // m_desLowerBounds[i] = T(0.0 + m_delta_s); // lower bound for design variable
        //     // m_desUpperBounds[i] = T(1.0 - m_delta_s); // upper bound for design variable
            m_desLowerBounds[i] = std::max(m_curDesign(i, 0) - 0.15, T(-0.2));
            m_desUpperBounds[i] = std::min(m_curDesign(i, 0) + 0.15, T(1.2));
        }

        // m_numConstraints = 0;
        // m_conJacRows.resize(m_numConstraints);
        // m_conJacCols.resize(m_numConstraints);
        // m_conLowerBounds.resize(m_numConstraints);
        // m_conUpperBounds.resize(m_numConstraints);
    }

    void updateBounds(const gsAsConstVector<T> &m_curDesign)
    {
        for (index_t i = 0; i < m_numDesignVars; ++i)
        {
            m_desLowerBounds[i] = std::max(m_curDesign(i, 0) - 0.15, T(-0.2));
            m_desUpperBounds[i] = std::min(m_curDesign(i, 0) + 0.15, T(1.2));
        }
    }

    static gsVector<T> vectorUpdate(const gsMultiPatch<T> &curve, const gsDofMapper &mapper)
    {
        gsVector<T> result(mapper.freeSize());
        for (index_t i = 0; i != curve.patch(0).coefs().rows(); ++i)
        {
            for (index_t j = 0; j != curve.patch(0).coefs().cols(); ++j)
            {
                const index_t gl = mapper.index(i,0,j);
                if (mapper.is_free_index(gl))
                result[gl] = curve.patch(0).coefs()(i,j);
            }
        }
        return result;
    }

    static void geomUpdate(const gsAsConstVector<T> &u, gsMultiPatch<T> &curve, const gsDofMapper &mapper)
    {
        for (index_t i=0; i!=curve.patch(0).coefs().rows(); ++i)
        {
            for (index_t j=0; j!=curve.patch(0).coefs().cols(); ++j)
            {
                const index_t gl = mapper.index(i,0,j);
                if (mapper.is_free_index(gl))
                curve.patch(0).coefs()(i,j) = u[gl];
            }
        }
    }

    static void geomUpdate(const gsVector<T> &u, gsMultiPatch<T> &curve, const gsDofMapper &mapper)
    {
        gsAsConstVector<T> tmp(u.data(),u.size());
        gsShapeOptProblem<T>::geomUpdate(tmp,curve,mapper);
    }

    T evalObj(const gsAsConstVector<T> &u) const override
    {
        //Update bounds at the beginning of each objective evaluation
        const_cast<gsShapeOptProblem<T>*>(this)->updateBounds(u);
        //gsDebug << "\nLower bounds: " << m_desLowerBounds.transpose() << "\n";
        //gsDebug << "\nUpper bounds: " << m_desUpperBounds.transpose() << "\n";

        gsDebug << "Computing objective at point " << u.transpose() << "\n";
        // Update embedded rib geometry from current design u
        gsMultiPatch<> tmpCurveOpt = m_curve;
        this->geomUpdate(u,tmpCurveOpt,m_mapper);

        // h-refine embedded rib for conforming quadrature
        index_t numPatches = tmpCurveOpt.nPatches();
        gsMultiPatch<> tmpCurveAn = tmpCurveOpt;
        gsMatrix<T> quPoints, allQuPoints, allQuWeights;
        gsVector<T> quWeights;

        embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn.patch(0),quPoints,quWeights,m_assembler);
        index_t numQuadPoints = quPoints.cols();

        allQuPoints.resize(numPatches, numQuadPoints);
        allQuWeights.resize(numPatches, numQuadPoints);
        allQuPoints.topRows(1) = quPoints;               //store quadrature points for the first patch
        allQuWeights.topRows(1) = quWeights.transpose(); //store quadrature weights for the first patch
        //gsDebug << "quPoints from obj:" << quPoints << "\n";
        //gsDebug << quPoints.size() <<"\n";

        for (index_t p = 1; p < numPatches; ++p)
        {
            embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn.patch(p),quPoints,quWeights,m_assembler);
            allQuPoints.middleRows(p,1) = quPoints;
            allQuWeights.middleRows(p,1) = quWeights.transpose();
        }

        // Re-assemble the system after updating the rib geometry for analysis
        ThinShellAssemblerStatus status_embedded = m_assembler->assembleLinearEmbeddedCurve(tmpCurveAn,
                                               m_EA, m_EI_min, m_EI_max, m_GI_p, allQuPoints, allQuWeights);
        GISMO_ENSURE(status_embedded == ThinShellAssemblerStatus::Success, "Embedded beam assembly failed");
        gsSparseMatrix<> K_s = m_K_shell + m_assembler->matrix();
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(K_s);
        gsVector<> u_s = solver.solve(m_F_shell);
        T obj = 0.5 * u_s.transpose() * K_s * u_s;

        // Return the objective function, i.e.strain energy, at current design
        gsDebug << "Objective: " << obj << " at point " << u.transpose() << "\n";
        return obj;
    }

    //void gradObj_analytical_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const
    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const override
    {
        // Update bounds at the beginning of each gradient evaluation
        const_cast<gsShapeOptProblem<T>*>(this)->updateBounds(u);
        //gsDebug << "\nLower bounds: " << m_desLowerBounds.transpose() << "\n";
        //gsDebug << "\nUpper bounds: " << m_desUpperBounds.transpose() << "\n";

        gsDebug << "Computing gradient at point " << u.transpose() << "\n";
        result.resize(m_numDesignVars);

        // Update embedded rib geometry from current design u
        gsMultiPatch<> tmpCurveOpt = m_curve;
        this->geomUpdate(u,tmpCurveOpt,m_mapper);

        // h-refine embedded rib for conforming quadrature
        index_t numPatches = tmpCurveOpt.nPatches();
        gsMultiPatch<> tmpCurveAn = tmpCurveOpt;
        gsMatrix<T> quPoints, allQuPoints, allQuWeights;
        gsVector<T> quWeights;

        embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn.patch(0),quPoints,quWeights,m_assembler);
        index_t numQuadPoints = quPoints.cols();

        allQuPoints.resize(numPatches, numQuadPoints);
        allQuWeights.resize(numPatches, numQuadPoints);
        allQuPoints.topRows(1) = quPoints;               //store quadrature points for the first patch
        allQuWeights.topRows(1) = quWeights.transpose(); //store quadrature weights for the first patch

        for (index_t p = 1; p < numPatches; ++p)
        {
            embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn.patch(p),quPoints,quWeights,m_assembler);
            allQuPoints.middleRows(p, 1) = quPoints;
            allQuWeights.middleRows(p, 1) = quWeights.transpose();
        }

        //gsDebug << "quPoints from grad for unperturbed:" << quPoints << "\n";
        //gsDebug << quPoints.size() <<"\n";
        
        // Assemble base system with current design
        ThinShellAssemblerStatus status = m_assembler->assembleLinearEmbeddedCurve(tmpCurveAn,
                                 m_EA, m_EI_min, m_EI_max, m_GI_p, allQuPoints, allQuWeights);
        GISMO_ENSURE(status == ThinShellAssemblerStatus::Success, "Embedded assembly failed");
        gsSparseMatrix<> K_embedded = m_assembler->matrix();
        gsSparseMatrix<> K_s = m_K_shell + K_embedded;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(K_s);
        gsVector<> u_s = solver.solve(m_F_shell);
        //gsDebugVar(u_s.transpose());
        
        // Compute pseudo load matrix R*
        gsMatrix<> R_star(m_numDofs,m_numDesignVars);
        for (index_t i = 0; i != tmpCurveOpt.patch(0).coefs().rows(); ++i)
        {
            for (index_t j = 0; j != tmpCurveOpt.patch(0).coefs().cols(); ++j)
            { 
                index_t gl = m_mapper.index(i,0,j);
                if (!m_mapper.is_free_index(gl)) continue;

                gsMultiPatch<> tmpCurveOpt_splusds = tmpCurveOpt;
                tmpCurveOpt_splusds.patch(0).coefs()(i,j) += m_delta_s;
                gsMultiPatch<> tmpCurveAn_splusds = tmpCurveOpt_splusds;
                embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn_splusds.patch(0),quPoints,quWeights,m_assembler);
                index_t numQuadPoints = quPoints.cols();
                allQuPoints.resize(numPatches, numQuadPoints);
                allQuWeights.resize(numPatches, numQuadPoints);
                allQuPoints.topRows(1) = quPoints;
                allQuWeights.topRows(1) = quWeights.transpose();
                //gsDebug << "quPoints from grad for perturbed:" << quPoints<< "\n";
                //gsDebug << quPoints.size() <<"\n";
                
                for (index_t p = 1; p < numPatches; ++p)
                {
                    embeddedQuadraturePoints(m_geom.patch(0),tmpCurveAn_splusds.patch(p),quPoints,quWeights,m_assembler);
                    allQuPoints.middleRows(p, 1) = quPoints;
                    allQuWeights.middleRows(p, 1) = quWeights.transpose();
                }

                ThinShellAssemblerStatus status = m_assembler->assembleLinearEmbeddedCurve(tmpCurveAn_splusds,
                                         m_EA, m_EI_min, m_EI_max, m_GI_p, allQuPoints, allQuWeights);
                GISMO_ENSURE(status == ThinShellAssemblerStatus::Success, "Embedded assembly failed");
                gsSparseMatrix<> K_embedded_splusds = m_assembler->matrix();
                gsSparseMatrix<> K_splusds = m_K_shell + K_embedded_splusds;
                
                R_star.col(gl) =  (-0.5 * (K_splusds - K_s) * u_s)/ m_delta_s; //since linear rhs does not depend on embedded ribs
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
    index_t numRefine  = 0;  //2
    index_t numRefineRib = 0;
    gsCmdLine cmd("Shell shape optimization based on strain energy.");
    cmd.addInt( "r", "numRefine", "Number of uniform h-refinement steps to perform on shell",  numRefine );
    cmd.addInt( "R", "numRefineRib", "Number of uniform h-refinement steps to perform on ribs",  numRefineRib );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Shell reference geometry for analysis]
    gsMultiPatch<> mp_surf;
    mp_surf.addPatch(gsNurbsCreator<>::BSplineSquare(10));
    mp_surf.embed(3);
    mp_surf.addAutoBoundaries();
    mp_surf.degreeElevate(1,-1); //set degree to 2 in both directions
    for (int r = 0; r < numRefine; ++r)
         mp_surf.uniformRefine();
    gsMultiBasis<> mbasis_surf(mp_surf);
    //gsGeometry<real_t> &surfgeo = mp_surf.patch(0);
    //gsTensorBSpline<2, real_t>* surf = dynamic_cast< gsTensorBSpline<2, real_t>* >(&surfgeo);
    gsInfo << "\nShell reference geometry\n";
    gsInfo << "Patches: "<< mp_surf.nPatches() <<", degree: "<< mbasis_surf.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surf.basis(0)<<"\n";
    gsWriteParaview(mp_surf, "initialDesign", 1000, true, true);
    //! [Shell reference geometry for analysis]

    //! [Embedded ribs]
    std::vector<real_t> vect_c = {0,0,0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1,1,1};
    gsKnotVector<real_t> kv_c(vect_c,2); 
    gsBSplineBasis<> basis_c(kv_c);
    gsMatrix<real_t> coefs_c(basis_c.size(), 2); 
    coefs_c << 0, 0.875,  //(u,v)
               0.14, 0.897,
               0.41, 0.979,
               0.668, 1.11,
               0.95, 1.098,
               1.17, 0.869,
               1.07, 0.59,
               0.85, 0.407,
               0.733, 0.144,
               0.736, 0;
    gsBSpline<> curveOpt(basis_c, give(coefs_c));
    gsMultiPatch<> mp_curveOpt;
    mp_curveOpt.addPatch(curveOpt);
    for (int r = 0; r < numRefineRib; ++r)
         mp_curveOpt.uniformRefine();
    gsMultiBasis<> mbasis_curveOpt(mp_curveOpt);
    gsInfo << "\nEmbedded rib for optimization\n";
    gsInfo << "Patches: "<< mp_curveOpt.nPatches() <<", degree: "<< mbasis_curveOpt.minCwiseDegree() <<"\n";
    gsInfo << mbasis_curveOpt.basis(0)<<"\n";
    gsWriteParaview(mp_curveOpt,  "ribNet",  1000, true,  false);
    //! [Embedded ribs]

    //! [Material properties of shell and embedded ribs]
    real_t E_modulus = 210e9; // [Pa]
    real_t PoissonRatio = 0.30;
    real_t Density = 7.85e3; // [kg/m^3]
    real_t thickness = 3.0e-3; // [m]

    real_t E_modulus_b = E_modulus; // [Pa]
    real_t PoissonRatio_b = 0.30;
    real_t thickness_b = 3e-3; // [m]
    real_t height_b = 30e-3; // [m]
    real_t G_modulus_b = 0.5 * E_modulus_b / (1+PoissonRatio_b);
    real_t EA = E_modulus_b * (height_b * thickness_b);                      //cross-sectional area
    real_t EI_min = E_modulus_b * (height_b * pow(thickness_b,3))/12;        //minimum principal inertia moment of beam cross-section
    real_t EI_max = E_modulus_b * (thickness_b * pow(height_b,3))/12;        //maximum principal inertia moment of beam cross-section
    real_t GI_p = G_modulus_b * (EI_min/E_modulus_b + EI_max/E_modulus_b);   //polar inertia moment of beam cross-section
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

    gsVector<> tmp(3);
    tmp << 0,0,-10; //[N/m^2]
    gsConstantFunction<> force(tmp,3); //pressure over shell

    bc.addCornerValue(boundary::southwest, 0.0, 0, 0, -1);
    bc.addCornerValue(boundary::southeast, 0.0, 0, 0, -1);
    bc.addCornerValue(boundary::northwest, 0.0, 0, 0, -1);
    bc.addCornerValue(boundary::northeast, 0.0, 0, 0, -1);
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
    gsVector<> F_shell = assembler->rhs();
    //! [Assemble shell linear part]

    //! [Collect control points on rib ends]
    gsDofMapper mapper(mbasis_curveOpt,mp_curveOpt.geoDim());
    for (short_t c = 1; c!=3; c++)
    {
        boxCorner corner(c);
        index_t idx = mbasis_curveOpt.basis(0).functionAtCorner(corner);
        for (short_t d = 0; d!=mp_curveOpt.geoDim(); ++d)
            mapper.eliminateDof(idx,0,d);
    }
    mapper.finalize();
    //! [Collect control points on rib ends]

    //! [Optimizer setup]
    gsShapeOptProblem<real_t> problem(assembler,mapper,mp_surf,mp_curveOpt,K_shell,F_shell,EA,EI_min,EI_max,GI_p);

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

    gsVector<> reshaped = gsShapeOptProblem<real_t>::vectorUpdate(mp_curveOpt,mapper);
    gsAsConstVector<> initialDesign(reshaped.data(), reshaped.size());
    gsDebugVar(problem.evalObj(initialDesign));
    gsDebugVar(initialDesign.transpose());

    //! [Solve]
    // Start optimization
    optimizer->solve(initialDesign);
    //! [Solve]

    // Get the optimized design
    gsVector<> optimizedDesign = optimizer->currentDesign();
    gsShapeOptProblem<real_t>::geomUpdate(optimizedDesign,mp_curveOpt,mapper);

    // Plot optimized design
    gsWriteParaview(mp_curveOpt, "OptimizedDesign", 1000, true, false);

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;
}