/** @file example_shell3D.cpp

    @brief Examples for the surface thin shells including the shell obstacle course

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;

    bool membrane = false;
    bool composite = false;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 0 = square plate with pressure; 1 = Scordelis Lo Roof; 2 = quarter hemisphere; 3 = pinched cylinder",  testCase );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("composite", "Composite material", composite);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    std::string fn;
    bool nonlinear = true;

    real_t aDim = 10.0;
    real_t bDim = 1.0;
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.patch(0).coefs().col(0) *= aDim;
    mp.patch(0).coefs().col(1) *= bDim;
    mp.computeTopology();
    mp.embed(3);
    E_modulus = 210e9;
    thickness = 1.0;
    PoissonRatio = 0.0;

    //! [Read input file]
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);


    //! [Refinement]
    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    // Basis for the exact solution
    gsMultiBasis<> exBasis(mp);
    exBasis.setDegree(4);

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    // Pinned-Pinned beam with Manufactured Solution
    real_t A  = bDim * thickness;
    real_t EA = E_modulus * A;
    real_t EI = E_modulus * bDim * math::pow(thickness,3) / 12.0;
    real_t q = 1e6;

    gsDebugVar(EI);
    gsDebugVar(EA);

    std::string fx,fy,fz;

    char buffer_u_inix[200];
    sprintf(buffer_u_inix,"%e*%e^2*(%e^3 - 6*%e*x^2 + 4*x^3)*x*(%e - x)/(48*%e^2)",EA,q,aDim,aDim,aDim,EI);
    fx = buffer_u_inix;

    fy = "0";

    char buffer_u_iniz[200];
    sprintf(buffer_u_iniz,"((x*%e*(%e - x)*(%e - 2*x)^2*(%e^2 + 2*%e*x - 2*x^2)^2*%e^2 + 768*%e^3)*%e)/(768*%e^3)",EA,aDim,aDim,aDim,aDim,q,EI,q,EI);
    fz = buffer_u_iniz;

    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    char buffer_ex[200];
    std::string ex_x, ex_y, ex_z;
    ex_x = "0";
    ex_y = "0";
    sprintf(buffer_ex,"%e * x * (%e^3 - 2*%e*x^2 + x^3) / (24 * %e)",q,aDim,aDim,EI);
    ex_z = buffer_ex;

    gsFunctionExpr<> ex(ex_x,ex_y,ex_z,3);
    gsField<> exf(mp,ex);
    gsWriteParaview<>(exf,"exact",1000);
    //! [Refinement]

    // Linear isotropic material model
    gsFunctionExpr<> force(fx,fy,fz,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    gsMaterialMatrixBase<real_t>* materialMatrix;
    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    // Set assembler
    gsThinShellAssemblerBase<real_t>* assembler;

    // Set stopwatch
    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      time += stopwatch.stop();
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      time += stopwatch.stop();
      return assembler->rhs();
    };

    gsField<> solField;
    gsMultiPatch<> deformation;
    std::vector<real_t> errors(numRefine+1);
    for (index_t r = 0; r!=numRefine+1; r++)
    {
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);
        // Define Matrices
        stopwatch.restart();
        stopwatch2.restart();
        assembler->assemble();
        time += stopwatch.stop();

        gsSparseMatrix<> matrix = assembler->matrix();
        gsVector<> vector = assembler->rhs();

        // Solve linear problem
        gsVector<> solVector;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute( matrix );
        solVector = solver.solve(vector);

        if (nonlinear)
        {
            real_t residual = vector.norm();
            real_t residual0 = residual;
            real_t residualOld = residual;
            gsVector<real_t> updateVector = solVector;
            gsVector<real_t> resVec = Residual(solVector);
            gsSparseMatrix<real_t> jacMat;
            for (index_t it = 0; it != 100; ++it)
            {
                jacMat = Jacobian(solVector);
                solver.compute(jacMat);
                updateVector = solver.solve(resVec); // this is the UPDATE
                solVector += updateVector;

                resVec = Residual(solVector);
                residual = resVec.norm();

                gsInfo<<"Iteration: "<< it
                   <<", residue: "<< residual
                   <<", update norm: "<<updateVector.norm()
                   <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
                   <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
                   <<"\n";

                residualOld = residual;

                if (updateVector.norm() < 1e-6)
                    break;
                else if (it+1 == it)
                    gsWarn<<"Maximum iterations reached!\n";
            }
        }

        totaltime += stopwatch2.stop();

        mp_def = assembler->constructSolution(solVector);
        deformation = assembler->constructDisplacement(solVector);

        solField = gsField(mp_def, deformation);

        errors[r] = solField.distanceL2(exf);

        gsDebugVar(solField.distanceL2(exf));

        mp.uniformRefine();
        mp_def.uniformRefine();
        dbasis = gsMultiBasis(mp);
    }


    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> exA(1,1);
    exA.setIntegrationElements(exBasis);
    space u = exA.getSpace(exBasis,3);
    geometryMap G = exA.getMap(mp);
    variable ff = exA.getCoeff(ex, G);

    gsMatrix<> projection;
    solution u_sol = exA.getSolution(u, projection);

    exA.initSystem();

    exA.assemble(u * u.tr(), u * ff);
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute(exA.matrix());
    projection = solver.solve(exA.rhs());

    gsMultiPatch<> exact;
    u_sol.extract(exact);

    gsField<> exactfield(mp,exact);
    gsDebugVar(exactfield.distanceL2(exf));

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);

        gsWriteParaview<>( exactfield, "Exact", 1000, false);

    }
    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsWriteParaview(membraneStress,"MembraneStress");
        gsWriteParaview(flexuralStress,"FlexuralStress");
    }
    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    return EXIT_SUCCESS;

}// end main
