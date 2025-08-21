/** @file shapeOpt_strainenergy.cpp

    @brief Structural optimization problem of
           a geometrically nonlinear Kirchhoff-Love shell
           with embedded ribs

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
    gsShapeOptProblem(  gsThinShellAssemblerBase<T>     *assembler,
                        const gsDofMapper               &mapper,
                        const gsMultiPatch<T>           &geom,
                        const gsMultiPatch<T>           &rib,
                        const std::vector<gsMatrix<T>>  &allquPointsCurve_rib,
                        const std::vector<gsVector<T>>  &allquWeights_rib,
                        const gsVector<T>               &materialParameters,
                        index_t                         &numRefineAn,
                        index_t                         &numRefineOpt)
    :
    m_assembler(assembler),
    m_mapper(mapper),
    m_geom(geom),
    m_rib(rib),
    m_allquPointsCurve_rib(allquPointsCurve_rib),
    m_allquWeights_rib(allquWeights_rib),
    m_materialParameters(materialParameters),
    m_numRefineAn(numRefineAn),
    m_numRefineOpt(numRefineOpt),
    m_delta_s(0.00001)
    {
        m_numDesignVars = m_mapper.freeSize();   // number of design variables
        m_numDofs = m_assembler->numDofs();      // number of unconstrained control points * dimensionality
        m_curDesign = this->vectorUpdate(m_geom, m_mapper);

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);
        for (index_t i = 0; i != geom.patch(0).coefs().rows(); ++i)
        {
            const index_t glx = mapper.index(i,0,0);
            if (mapper.is_free_index(glx))
            {
                m_desLowerBounds[glx] = m_curDesign(glx,0) - 0.5;  // x-coordinate
                m_desUpperBounds[glx] = m_curDesign(glx,0) + 0.5;
            }

            const index_t gly = mapper.index(i,0,1);
            if (mapper.is_free_index(gly))
            {
                m_desLowerBounds[gly] = m_curDesign(gly,0) - 0.5;  // y-coordinate
                m_desUpperBounds[gly] = m_curDesign(gly,0) + 0.5;
            }

            const index_t glz = mapper.index(i,0,2);
            if (mapper.is_free_index(glz))
            {
                m_desLowerBounds[glz] = m_curDesign(glz,0) - 0.2;  // z-coordinate
                m_desUpperBounds[glz] = m_curDesign(glz,0) + 0.2;
            }
        }

        //for (index_t i = 0; i != m_numDesignVars; ++i)
        // {
        //     m_desLowerBounds[i] = m_curDesign(i,0) - 0.7;
        //     m_desUpperBounds[i] = m_curDesign(i,0) + 0.7;
        // }

        // m_numConstraints = 0;
        // m_conJacRows.resize(m_numConstraints);
        // m_conJacCols.resize(m_numConstraints);
        // m_conLowerBounds.resize(m_numConstraints);
        // m_conUpperBounds.resize(m_numConstraints);
    }

    static gsVector<T> vectorUpdate(const gsMultiPatch<T> & geom, const gsDofMapper & mapper)
    {
        gsVector<T> result(mapper.freeSize());
        for (index_t i = 0; i != geom.patch(0).coefs().rows(); ++i)
        {
            for (index_t j = 0; j != geom.patch(0).coefs().cols(); ++j)
            {
                const index_t gl = mapper.index(i,0,j);
                if (mapper.is_free_index(gl))
                result[gl] = geom.patch(0).coefs()(i,j);
            }
        }
        return result;
    }

    static void geomUpdate(const gsAsConstVector<T> &u, gsMultiPatch<T> & geom, const gsDofMapper & mapper)
    {
        for (index_t i=0; i!=geom.patch(0).coefs().rows(); ++i)
        {
            for (index_t j=0; j!=geom.patch(0).coefs().cols(); ++j)
            {
                const index_t gl = mapper.index(i,0,j);
                if (mapper.is_free_index(gl))
                geom.patch(0).coefs()(i,j) = u[gl];
            }
        }
    }

    static void geomUpdate(const gsVector<T> &u, gsMultiPatch<T> & geom, const gsDofMapper & mapper)
    {
        gsAsConstVector<T> tmp(u.data(),u.size());
        gsShapeOptProblem<T>::geomUpdate(tmp,geom,mapper);
    }

    T evalObj(const gsAsConstVector<T> &u) const override
    {
        // Update geometry coefficients from current design u
        gsDebug << "Computing objective at point " << u.transpose() << "\n";
        gsMultiPatch<> tmpGeom = m_geom;
        this->geomUpdate(u,tmpGeom,m_mapper);

        // Re-assemble the rib-enforced system after updating the analysis geometry
        gsMultiPatch<> anGeom = tmpGeom;
        index_t m_numRefineDiff = m_numRefineAn - m_numRefineOpt;
        for (int r = 0; r < m_numRefineDiff; ++r)
             anGeom.uniformRefine();
        m_assembler->setGeometry(anGeom);
        ThinShellAssemblerStatus status = m_assembler->assemble();
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        gsInfo<<"Setting up shell assembly\n";
        gsSparseMatrix<> K_s = m_assembler->matrix();
        gsVector<> F_s = m_assembler->rhs();
        gsInfo<<"Shell assembly done\n";

        gsInfo <<"Setting up rib assembly\n";
        T EA_rib = m_materialParameters[0];
        T EI_min_rib = m_materialParameters[1];
        T EI_max_rib = m_materialParameters[2];
        T GI_p_rib = m_materialParameters[3];
        ThinShellAssemblerStatus status_rib;
        status_rib = m_assembler->assembleLinearEmbeddedCurve(m_rib,EA_rib,EI_min_rib,EI_max_rib,GI_p_rib,m_allquPointsCurve_rib,m_allquWeights_rib);
        GISMO_ENSURE(status_rib==ThinShellAssemblerStatus::Success,"Assembly for rib embedding failed");
        gsSparseMatrix<> K_embedded = m_assembler->matrix();
        K_s += K_embedded;

        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(K_s);
        gsVector<> u_s = solver.solve(F_s);

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

        // Update geometry coefficients from current design u
        gsMultiPatch<> tmpGeom = m_geom;
        this->geomUpdate(u,tmpGeom,m_mapper);

        // Re-assemble the rib-enforced system after updating the geometry
        gsMultiPatch<> anGeom = tmpGeom;
        index_t m_numRefineDiff = m_numRefineAn - m_numRefineOpt;
        for (int r = 0; r < m_numRefineDiff; ++r)
             anGeom.uniformRefine();
        m_assembler->setGeometry(anGeom);
        ThinShellAssemblerStatus status = m_assembler->assemble();
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        //gsInfo<<"Setting up shell assembly\n";
        gsSparseMatrix<> K_s = m_assembler->matrix();
        gsVector<> F_s = m_assembler->rhs();
        //gsInfo<<"Shell assembly done\n";

        //gsInfo <<"Setting up rib assembly\n";
        T EA_rib = m_materialParameters[0];
        T EI_min_rib = m_materialParameters[1];
        T EI_max_rib = m_materialParameters[2];
        T GI_p_rib = m_materialParameters[3];
        ThinShellAssemblerStatus status_rib;
        status_rib = m_assembler->assembleLinearEmbeddedCurve(m_rib,EA_rib,EI_min_rib,EI_max_rib,GI_p_rib,m_allquPointsCurve_rib,m_allquWeights_rib);
        GISMO_ENSURE(status_rib==ThinShellAssemblerStatus::Success,"Assembly for rib embedding failed");
        gsSparseMatrix<> K_embedded = m_assembler->matrix();
        K_s += K_embedded;

        gsSparseSolver<>::CGDiagonal m_solver;
        m_solver.compute(K_s);
        gsVector<> u_s = m_solver.solve(F_s);

        // Compute pseudo load matrix R*
        gsMatrix<> R_star(m_numDofs,m_numDesignVars);
        for (index_t i = 0; i != tmpGeom.patch(0).coefs().rows(); ++i)
        {
            for (index_t j = 0; j != tmpGeom.patch(0).coefs().cols(); ++j)
            {
                index_t gl = m_mapper.index(i,0,j);
                if (!m_mapper.is_free_index(gl)) continue;
                //gsMultiPatch<> anGeom_splusds = anGeom;
                //anGeom_splusds.patch(0).coefs()(i,j) += m_delta_s;
                gsMultiPatch<> tmpGeom_splusds = tmpGeom;
                tmpGeom_splusds.patch(0).coefs()(i,j) += m_delta_s;
                gsMultiPatch<> anGeom_splusds = tmpGeom_splusds;
                for (int r =0; r < m_numRefineDiff; ++r)
                        anGeom_splusds.uniformRefine();
                m_assembler->setGeometry(anGeom_splusds);
                m_assembler->assemble();
                gsSparseMatrix<> K_splusds = m_assembler->matrix();
                gsVector<> F_splusds = m_assembler->rhs();

                ThinShellAssemblerStatus status_rib;
                status_rib = m_assembler->assembleLinearEmbeddedCurve(m_rib,EA_rib,EI_min_rib,EI_max_rib,GI_p_rib,m_allquPointsCurve_rib,m_allquWeights_rib);
                GISMO_ENSURE(status_rib==ThinShellAssemblerStatus::Success,"Assembly for rib embedding failed");
                gsSparseMatrix<> K_embedded = m_assembler->matrix();
                K_splusds += K_embedded;

                R_star.col(gl) = ((F_splusds - F_s) - 0.5 * (K_splusds - K_s)* u_s)/m_delta_s;
            }
        }

        // Return sensitivity vector df/ds
        result = u_s.transpose() * R_star;
        //gsDebug << "Gradient: " << result.transpose() << " at point " << u.transpose() << "\n";
        //gsDebug << "Gradient norm: " << result.norm() << "\n";
    }

    void gradObj_FDM_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const
    {
        this->gradObj_into(u, result);
    }

protected:
    gsThinShellAssemblerBase<T>    *m_assembler;
    const gsDofMapper              &m_mapper;
    const gsMultiPatch<T>          &m_geom;
    const gsMultiPatch<T>          &m_rib;
    const std::vector<gsMatrix<T>> &m_allquPointsCurve_rib;
    const std::vector<gsVector<T>> &m_allquWeights_rib;
    const gsVector<T>              &m_materialParameters;
    index_t                         m_numRefineAn;
    index_t                         m_numRefineOpt;
    index_t                         m_numDofs;
    using Base::m_numDesignVars;
    using Base::m_curDesign;
    // using Base::m_numConstraints;
    using Base::m_desLowerBounds;
    using Base::m_desUpperBounds;
    // using Base::m_conLowerBounds;
    // using Base::m_conUpperBounds;
    // using Base::m_conJacRows;
    // using Base::m_conJacCols;
    T                               m_delta_s;
};

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    index_t numRefineAn  = 0;
    index_t numRefineOpt  = 0;
    std::string outputDir = "./output";

    gsCmdLine cmd("Strain-energy based linear optimization of rib-enforced shells by adjustment of shell geometry.");
    cmd.addInt( "A", "rAn", "Number of uniform h-refinement steps to perform before analysis",  numRefineAn );
    cmd.addInt( "O", "rOpt", "Number of uniform h-refinement steps to perform before optimization",  numRefineOpt );
    cmd.addString("o", "output", "Output directory", outputDir);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ENSURE(numRefineAn >= numRefineOpt,"Mesh refinement for analysis not coarser than for optimization");
    outputDir += gsFileManager::getNativePathSeparator();
    if (!gsFileManager::fileExists(outputDir))
        gsFileManager::mkdir(outputDir);

     //! [Shell reference geometry for analysis and optimization]
    gsMultiPatch<> mp_surfOpt;
    real_t r = 10; // [m]
    mp_surfOpt.addPatch(gsNurbsCreator<>::BSplineSquare(r));
    mp_surfOpt.embed(3);
    mp_surfOpt.addAutoBoundaries();
    mp_surfOpt.degreeElevate(1,-1); //set degree to 2 in both directions
    gsMultiPatch<> mp_surfAn = mp_surfOpt;

    for (int r = 0; r < numRefineOpt; ++r)
        mp_surfOpt.uniformRefine(); //h-refine;

    gsMultiBasis<> mbasis_surfOpt(mp_surfOpt);
    gsInfo << "\nShell reference geometry for optimization\n";
    gsInfo << "Patches: "<< mp_surfOpt.nPatches() <<", degree: "<< mbasis_surfOpt.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surfOpt.basis(0)<<"\n";
    gsDebug << mp_surfOpt.patch(0).coefs() << "\n";
    gsWriteParaview(mp_surfOpt, outputDir + "initialDesignOpt", 1000, true, true);

    for (int r = 0; r < numRefineAn; ++r)
         mp_surfAn.uniformRefine();

    gsMultiBasis<> mbasis_surfAn(mp_surfAn);
    gsInfo << "\nShell reference geometry for analysis\n";
    gsInfo << "Patches: "<< mp_surfAn.nPatches() <<", degree: "<< mbasis_surfAn.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surfAn.basis(0)<<"\n";
    gsDebug << mp_surfAn.patch(0).coefs() << "\n";
    gsWriteParaview(mp_surfAn, outputDir + "initialShellAn", 1000, true, true);

    gsGeometry<real_t> &surfgeo = mp_surfAn.patch(0);
    gsTensorBSpline<2, real_t>* surf = dynamic_cast< gsTensorBSpline<2, real_t>* >(&surfgeo);
    //! [Shell reference geometry for analysis and optimization]

    //! [Embedded beam features for analysis]
    gsKnotVector<real_t> kv_c = surf->knots(0);
    gsBSplineBasis<> basis_c(kv_c);

    gsEigen::ArrayXXd cpvec (surf->knots(0).size() - surf->degree(0) - 1, 1);
    cpvec = (surf->coefs().block(0,0,cpvec.rows(),1))/r;

    auto cpvec_flipped = cpvec.reverse();
    gsMatrix<real_t> coef_c1(basis_c.size(), surf->parDim());
    coef_c1.col(0) = cpvec;
    coef_c1.col(1) = cpvec_flipped;
    gsMatrix<real_t> coef_c2(basis_c.size(), surf->parDim());
    coef_c2.col(0) = cpvec;
    coef_c2.col(1) = cpvec;

    gsBSpline<> ribA(basis_c, coef_c1);
    gsBSpline<> ribB(basis_c, coef_c2);

    gsMultiPatch<> mp_rib;
    mp_rib.addPatch(ribA);
    mp_rib.addPatch(ribB);

    gsMatrix<real_t> coef_LLDPE(basis_c.size(), surf->parDim());
    coef_LLDPE.col(0) = cpvec;
    coef_LLDPE.col(1).setOnes();
    gsBSpline<> LLDPE_top(basis_c, coef_LLDPE);    // top buoyant breakwater
    mp_rib.addPatch(LLDPE_top);

    coef_LLDPE.col(1).setZero();
    gsBSpline<> LLDPE_bottom(basis_c, coef_LLDPE); // bottom buoyant breakwater
    mp_rib.addPatch(LLDPE_bottom);

    coef_LLDPE.col(0).setZero();
    coef_LLDPE.col(1) = cpvec;
    gsBSpline<> LLDPE_left(basis_c, coef_LLDPE);   // left buoyant breakwater
    mp_rib.addPatch(LLDPE_left);

    coef_LLDPE.col(0).setOnes();
    coef_LLDPE.col(1) = cpvec;
    gsBSpline<> LLDPE_right(basis_c, coef_LLDPE);  // right buoyant breakwater
    mp_rib.addPatch(LLDPE_right);

    gsWriteParaview(mp_rib,  outputDir + "ribs",  1000, true,  true);
    gsMultiBasis<> mbasis_rib(mp_rib);
    //! [Embedded beam features for analysis]

    //! [Mechanical properties of shell and embedded entities]
    real_t E_modulus = 1.5e9; // [Pa] HDPE
    real_t PoissonRatio = 0.45;
    real_t density = 950; // [kg/m^3]
    real_t thickness = 5.0e-3; // [m]

    real_t E_modulus_rib = 1.5e9; // [Pa] HDPE
    real_t PoissonRatio_rib = 0.45;
    real_t thickness_rib = 0.10; // [m] 3e-3
    real_t height_rib = 0.15;   // [m] 30e-3
    real_t G_modulus_rib = 0.5 * E_modulus_rib / (1 + PoissonRatio_rib);
    real_t EA_rib = E_modulus_rib * (height_rib * thickness_rib);               //axial rigidity
    real_t EI_min_rib = E_modulus_rib * (height_rib * pow(thickness_rib,3))/12; //minimum flexural rigidity
    real_t EI_max_rib = E_modulus_rib * (thickness_rib * pow(height_rib,3))/12; //maximum flexural rigidity
    real_t GI_p_rib = G_modulus_rib/E_modulus_rib * (EI_min_rib + EI_max_rib);  //torsional rigidity

    gsVector<real_t> materialParameters(4);
    materialParameters << EA_rib, EI_min_rib, EI_max_rib, GI_p_rib;
    //! [Mechanical properties of shell and embedded entities]

    //! [Make material functions]
    // Linear isotropic material model
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(density),3);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    //! [Make material functions]

    //! [Set boundary conditions and loads]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp_surfAn);
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
    materialMatrix = getMaterialMatrix<3,real_t>(mp_surfAn,t,parameters,rho,options);
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t,true>(mp_surfAn,mbasis_surfAn,bc,force,materialMatrix);
    //assembler->setPointLoads(pLoads);
    //! [Make assembler]

    //! [h-refine embedded curves based on mp_surfAn for conforming quadrature]
    index_t numPatches_rib = mp_rib.nPatches();

    std::vector<gsMatrix<>> allquPointsCurve_rib;
    std::vector<gsVector<>> allquWeights_rib;

    for (index_t p = 0; p < numPatches_rib; ++p)
    {
        gsMatrix<> quPointsCurve;         gsVector<> quWeights;
        embeddedQuadraturePoints(mp_surfAn.patch(0),mp_rib.patch(p),quPointsCurve,quWeights);
        allquPointsCurve_rib.push_back(quPointsCurve);
        allquWeights_rib.push_back(quWeights);
    }
    //! [h-refine embedded curves based on mp_surfAn for conforming quadrature]

    //! [Freeze z-dof on shell boundaries]
    gsDofMapper mapper(mbasis_surfOpt, mp_surfOpt.geoDim());
    for (index_t side = 0; side < 4; ++side)
    {
        gsMatrix<index_t> boundaryIndices = mbasis_surfOpt.basis(0).boundary(side);

        for (index_t i = 0; i < boundaryIndices.rows(); ++i)
        {
            index_t globalIndex = boundaryIndices(i, 0); // only one column
            mapper.eliminateDof(globalIndex,0,2);
        }
    }
    mapper.finalize();
    //! [Freeze z-dof on shell boundaries]

    //! [Optimizer setup]
    gsShapeOptProblem<real_t> problem(assembler,mapper,mp_surfOpt,mp_rib,
                                      allquPointsCurve_rib,allquWeights_rib,
                                      materialParameters,numRefineAn,numRefineOpt);

    gsOptimizer<real_t> *optimizer;
#ifdef gsOptim_ENABLED
    optimizer = new gsOptim<real_t>::LBFGS(&problem);
#else
    optimizer = new gsGradientDescent<>(&problem);
    optimizer->options().setReal("MinGradientLength",1e-9);
    optimizer->options().setReal("MinStepLength",1e-9);
#endif
    optimizer->options().setInt("MaxIterations",50);
    optimizer->options().setInt("Verbose",1);
    //optimizer->options().setReal("GradErrTol",1e-8);
    //! [Optimizer setup]

    gsVector<> reshaped = gsShapeOptProblem<real_t>::vectorUpdate(mp_surfOpt, mapper);
    gsAsConstVector<> initialDesign(reshaped.data(), reshaped.size());

    // gsMatrix<> mat_FDM(mapper.freeSize(),1);
    // gsAsVector<> sensitivities_FDM(mat_FDM.data(),mat_FDM.rows());
    // problem.gradObj_FDM_into(initialDesign,sensitivities_FDM);
    // gsInfo<<"\nNumerical sensitivity vector:\n";
    // gsDebugVar(sensitivities_FDM.transpose());
    // gsDebugVar(sensitivities_FDM.norm());
    // return EXIT_SUCCESS;

    //! [Solve]
    // Start optimization
    optimizer->solve(initialDesign);
    //! [Solve]

    // Get the optimized design
    gsVector<> optimizedDesign = optimizer->currentDesign();
    gsShapeOptProblem<real_t>::geomUpdate(optimizedDesign,mp_surfOpt,mapper);

    // Plot optimized design
    gsWrite(mp_surfOpt, outputDir + "OptimalShell"); //.xml file of optimal shell geometry
    gsWrite(mp_rib, outputDir + "RibTopology"); //.xml file of rib geometry
    gsWriteParaview(mp_surfOpt, outputDir + "OptimizedDesign", 1000, true, false);

    // ****** VALIDATION ****** //

    // //Evaluate the objective function
    // gsShapeOptProblem<real_t> SOP(assembler,mapper,mp_surfOpt,numRefineAn,numRefineOpt);
    // gsDebug<<SOP.evalObj(initialDesign)<<"\n";

    // //Evaluate the sensitivity vector
    // gsMatrix<> mat(mapper.freeSize(),1);
    // gsAsVector<> sensitivities(mat.data(),mat.rows());
    // SOP.gradObj_analytical_into(initialDesign,sensitivities);
    // gsInfo<<"\nAnalytical sensitivity vector:\n";
    // gsDebugVar(sensitivities.transpose());

    // gsMatrix<> mat_FDM(mapper.freeSize(),1);
    // gsAsVector<> sensitivities_FDM(mat_FDM.data(),mat_FDM.rows());
    // SOP.gradObj_FDM_into(initialDesign,sensitivities_FDM);
    // gsInfo<<"\nNumerical sensitivity vector:\n";
    // gsDebugVar(sensitivities_FDM.transpose());

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;
}
