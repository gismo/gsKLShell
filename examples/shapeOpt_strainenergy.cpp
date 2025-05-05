/** @file assembly_pseudoR.cpp

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
#include <gsAssembler/gsEmbeddingUtils.h>
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
                        index_t                     &numRefineAn,
                        index_t                     &numRefineOpt)
    :
    m_assembler(assembler),
    m_mapper(mapper),
    m_geom(geom),
    m_numRefineAn(numRefineAn),
    m_numRefineOpt(numRefineOpt),
    m_delta_s(0.00001)
    {
        m_numDesignVars = m_mapper.freeSize();   // number of design variables
        m_numDofs = m_assembler->numDofs();      // number of unconstrained control points * dimensionality
        m_curDesign = this->vectorUpdate(m_geom, m_mapper);

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);
        for (index_t i = 0; i != m_numDesignVars; ++i)
        {
            m_desLowerBounds[i] = m_curDesign(i,0) - 0.7;
            m_desUpperBounds[i] = m_curDesign(i,0) + 0.7;
        }

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
        gsMultiPatch<> tmpGeom = m_geom;
        this->geomUpdate(u,tmpGeom,m_mapper);

        // Re-assemble the system after updating the analysis geometry
        gsMultiPatch<> anGeom = tmpGeom;
        index_t m_numRefineDiff = m_numRefineAn - m_numRefineOpt;
        for (int r =0; r < m_numRefineDiff; ++r)
             anGeom.uniformRefine();
        m_assembler->setGeometry(anGeom);
        ThinShellAssemblerStatus status = m_assembler->assemble();
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(m_assembler->matrix());
        gsVector<> u_s = solver.solve(m_assembler->rhs());

        // Return the objective function, i.e.strain energy, at current design
        return 0.5 * u_s.transpose() * m_assembler->matrix() * u_s;
    }

    //void gradObj_analytical_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const
    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const override
    {
        result.resize(m_numDesignVars);

        // Update geometry coefficients from current design u
        gsMultiPatch<> tmpGeom = m_geom;
        this->geomUpdate(u,tmpGeom,m_mapper);

        // Re-assemble the system after updating the geometry
        gsMultiPatch<> anGeom = tmpGeom;
        index_t m_numRefineDiff = m_numRefineAn - m_numRefineOpt;
        for (int r =0; r < m_numRefineDiff; ++r)
             anGeom.uniformRefine();
        m_assembler->setGeometry(anGeom);
        ThinShellAssemblerStatus status = m_assembler->assemble();
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        gsSparseMatrix<> K_s = m_assembler->matrix();
        gsVector<> F_s = m_assembler->rhs();
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
                if (m_mapper.is_free_index(gl))
                {
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
                    R_star.col(gl) = ((F_splusds - F_s) - 0.5 * (K_splusds - K_s)* u_s)/m_delta_s;
                }
            }
        }

        // Return sensitivity vector df/ds
        result = u_s.transpose() * R_star;
    }

    void gradObj_FDM_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const
    {
        this->gradObj_into(u, result);
    }

protected:
    gsThinShellAssemblerBase<T>    *m_assembler;
    const gsDofMapper              &m_mapper;
    const gsMultiPatch<T>          &m_geom;
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
    index_t numRefineAn  = 2;
    index_t numRefineOpt  = 1;
    GISMO_ASSERT(numRefineAn >= numRefineOpt,"Mesh refinement for analysis not coarser than for optimization");

    gsCmdLine cmd("Shell shape optimization based on strain energy.");
    cmd.addInt( "A", "rAn", "Number of uniform h-refinement steps to perform before analysis",  numRefineAn );
    cmd.addInt( "O", "rOpt", "Number of uniform h-refinement steps to perform before optimization",  numRefineOpt );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Shell material properties]
    real_t E_modulus = 210e9; // [Pa] 210e9
    real_t PoissonRatio = 0.30;
    real_t Density = 7.85e3; // [kg/m^3]
    real_t thickness = 3.0e-3; // [m]
    //! [Shell material properties]

    //! [Shell reference geometry for analysis and optimization]
    gsKnotVector<> kv1(0, 1, 2, 3);
    gsKnotVector<> kv2(0, 1, 2, 3);
    gsTensorBSplineBasis<2, real_t> basis_s(kv1, kv2);
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
    gsTensorBSpline<2, real_t>  surf(basis_s, coefs_s);

    gsMultiPatch<> mp_surfOpt;
    mp_surfOpt.addPatch(surf);
    for (int r =0; r < numRefineOpt; ++r)
         mp_surfOpt.uniformRefine();

    gsMultiBasis<> mbasis_surfOpt(mp_surfOpt);
    gsInfo << "\nShell reference geometry for optimization\n";
    gsInfo << "Patches: "<< mp_surfOpt.nPatches() <<", degree: "<< mbasis_surfOpt.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surfOpt.basis(0)<<"\n";
    gsWriteParaview(mp_surfOpt, "initialDesignOpt", 1000, true, true);

    gsMultiPatch<> mp_surfAn;
    mp_surfAn.addPatch(surf);
    mp_surfAn.addAutoBoundaries();
    for (int r =0; r < numRefineAn; ++r)
         mp_surfAn.uniformRefine();

    gsMultiBasis<> mbasis_surfAn(mp_surfAn);
    gsInfo << "\nShell reference geometry for analysis\n";
    gsInfo << "Patches: "<< mp_surfAn.nPatches() <<", degree: "<< mbasis_surfAn.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surfAn.basis(0)<<"\n";
    gsWriteParaview(mp_surfAn, "initialDesignAn", 1000, true, true);
    //! [Shell reference geometry for analysis and optimization]

    //! [Set boundary conditions and loads]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp_surfAn);

    // Pressure under shell
    gsVector<> tmp(3);
    tmp << 0,0,0; //[N/m^2]
    gsConstantFunction<> force(tmp,3);

    gsVector<> point(2);
    point<<0.5,0.5;
    gsVector<> load(3);
    load << 0,0,-5; //[N/m^2]
    gsPointLoads<real_t> pLoads;
    pLoads.addLoad(point,load,0,true);

    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, i);
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false,  i);
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false,  i);
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, i);
        // bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)
        // bc.addCornerValue(boundary::southeast, 0.0, 0, 0); // (corner,value, patch, unknown)
        // bc.addCornerValue(boundary::northwest, 0.0, 0, 0); // (corner,value, patch, unknown)
        // bc.addCornerValue(boundary::northeast, 0.0, 0, 0); // (corner,value, patch, unknown)

    }
    //! [Set boundary conditions and loads]

    //! [Make material functions]
    // Linear isotropic material model
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    //! [Make material functions]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    materialMatrix = getMaterialMatrix<3,real_t>(mp_surfAn,t,parameters,rho,options);
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t,true>(mp_surfAn,mbasis_surfAn,bc,force,materialMatrix);
    assembler->setPointLoads(pLoads);
    //! [Make assembler]

    //! [Collect control points on shell sides]
    gsDofMapper mapper(mbasis_surfOpt,mp_surfOpt.geoDim());
    for (short_t c = 1; c!=5; c++)
    {
        boxCorner corner(c);
        index_t idx = mbasis_surfOpt.basis(0).functionAtCorner(corner);
        for (short_t d = 0; d!=mp_surfOpt.geoDim(); ++d)
            mapper.eliminateDof(idx,0,d);
    }
    mapper.finalize();
    //! [Collect control points on shell sides]

    //! [Optimizer setup]
    gsShapeOptProblem<real_t> problem(assembler, mapper, mp_surfOpt , numRefineAn, numRefineOpt);

    gsOptimizer<real_t> *optimizer;
#ifdef gsOptim_ENABLED
    optimizer = new gsOptim<real_t>::LBFGS(&problem);
#else
    optimizer = new gsGradientDescent<>(&problem);
    optimizer->options().setReal("MinGradientLength",1e-9);
    optimizer->options().setReal("MinStepLength",1e-9);
#endif
    optimizer->options().setInt("MaxIterations",100);
    optimizer->options().setInt("Verbose",1);
    //optimizer->options().setReal("GradErrTol",1e-8);
    //! [Optimizer setup]

    gsVector<> reshaped = gsShapeOptProblem<real_t>::vectorUpdate(mp_surfOpt, mapper);
    gsAsConstVector<> initialDesign(reshaped.data(), reshaped.size());

    //! [Solve]
    // Start optimization
    optimizer->solve(initialDesign);
    //! [Solve]

    // Get the optimized design
    gsVector<> optimizedDesign = optimizer->currentDesign();
    gsShapeOptProblem<real_t>::geomUpdate(optimizedDesign,mp_surfOpt,mapper);

    // Plot optimized design
    gsWriteParaview(mp_surfOpt, "OptimizedDesign", 1000, true, false);

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