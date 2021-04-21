/** @file example_shell2D.cpp

    @brief Simple 2D examples for the shell class

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
    bool Compressibility = false;
    index_t material = 0;
    bool verbose = false;
    std::string fn;
    bool membrane = false;

    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("composite", "Composite material", composite);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase == 0)
    {
        real_t mu = 1.5e6;
        thickness = 0.001;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);
    }
    else if (testCase == 1)
    {
        E_modulus = 1;
        thickness = 1;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;
    }
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();

    gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";
    gsDebugVar(PoissonRatio);

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(2);
    tmp << 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    real_t pressure = 0.0;
    if (testCase == 0) // Uniaxial tension; use with hyperelastic material model!
    {
        gsVector<> neu(2);
        neu << 2625, 0;
        gsConstantFunction<> neuData(neu,2);

        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
    }
    else if (testCase == 1)
    {
        for (index_t i=0; i!=2; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (2); load << 0.25, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else
        GISMO_ERROR("Test case not known");

    //! [Refinement]

    // Linear isotropic material model and Neo-Hookean material
    gsConstantFunction<> force(tmp,2);
    gsConstantFunction<> pressFun(pressure,2);
    gsFunctionExpr<> t(std::to_string(thickness),2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
    gsFunctionExpr<> rho(std::to_string(Density),2);

    // Mooney-Rivlin material
    gsConstantFunction<> ratio(Ratio,2);

    // Ogden material
    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,2);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,2);
    gsConstantFunction<> alpha2(5.0,2);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,2);
    gsConstantFunction<> alpha3(-2.0,2);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,2);
    gsMaterialMatrixBase<real_t>* materialMatrix;

    // Linear anisotropic material model
    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,2);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,2);
    Ts[0] = &thicks;

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==1 || material==2) // NH & NH_ext
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==3) // MR
    {
      parameters.resize(3);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &ratio;
    }
    else if (material==4) // OG
    {
      parameters.resize(8);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
      parameters[4] = &mu2;
      parameters[5] = &alpha2;
      parameters[6] = &mu3;
      parameters[7] = &alpha3;
    }

    gsOptionList options;
    if      (material==0)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<2,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrix);

    assembler->setPointLoads(pLoads);

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

    totaltime += stopwatch2.stop();

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);
    }
    if (stress)
    {

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
        gsField<> Stretches(mp_def,stretches, true);

        gsPiecewiseFunction<> pstress_m;
        assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
        gsField<> pstressM(mp_def,pstress_m, true);

        gsPiecewiseFunction<> pstress_f;
        assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
        gsField<> pstressF(mp_def,pstress_f, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_def,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_def,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_def,stretch3, true);


        gsWriteParaview(membraneStress,"MembraneStress");
        gsWriteParaview(flexuralStress,"FlexuralStress");
        gsWriteParaview(Stretches,"PrincipalStretch");
        gsWriteParaview(pstressM,"PrincipalMembraneStress");
        gsWriteParaview(pstressF,"PrincipalFlexuralStress");
        gsWriteParaview(stretchDir1,"PrincipalDirection1");
        gsWriteParaview(stretchDir1,"PrincipalDirection1");
        gsWriteParaview(stretchDir2,"PrincipalDirection2");
        gsWriteParaview(stretchDir3,"PrincipalDirection3");
    }

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    return EXIT_SUCCESS;

}// end main
