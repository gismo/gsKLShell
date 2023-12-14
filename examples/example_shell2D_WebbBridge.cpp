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
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>

#include <gsKLShell/gsMaterialMatrixTFT.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    real_t pi = 4*math::atan(1);
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool Compressibility = false;
    index_t material = 0;
    bool verbose = false;
    bool nonlinear = false;
    std::string fn;

    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 1.0e9;
    real_t PoissonRatio = 1.0/3.0;
    real_t Density = 1.0;
    real_t thickness = 1.0e-3;

    real_t Rbar = 1.2;

    gsCmdLine cmd("2D shell example.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addReal( "R", "Rbar", "Rbar",  Rbar );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("nl", "Nonlinear analysis", nonlinear);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    real_t H = 25e-3;
    real_t W  = 75e-3;
    real_t W1 = 2./5.*W;
    real_t W2 = 1./5.*W;

    //! [Read input file]
    gsMultiPatch<> mp, mp_mid;
    gsMultiPatch<> mp_def;

    gsTensorBSpline<2,real_t> patch = *(gsNurbsCreator<>::BSplineSquare());
    patch.coefs().col(0) *=W;
    patch.coefs().col(1) *=H;
    patch.insertKnot(0.4,0);
    patch.insertKnot(0.6,0);
    mp.addPatch(patch);

    // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,W1,H));
    // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(W1,0,W1+W2,H));
    // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(W1+W2,0,W,H));

    // mp.degreeReduce(1);
    // mp.computeTopology();
    // gsDebugVar(mp.patch(0));

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    // Set the deformed configuration
    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true,true);
    gsWrite(mp_def,"mp");

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";
    //! [Make geometry and refine/elevate]

    //! [Make material functions]
    // Linear isotropic material model and Neo-Hookean material
    gsVector<> tmp(2);
    tmp << 0, 0;
    gsConstantFunction<> force(tmp,2);
    gsFunctionExpr<> t(std::to_string(thickness),2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
    gsFunctionExpr<> rho(std::to_string(Density),2);

    gsConstantFunction<> mu1(749.18,2);
    gsConstantFunction<> alpha1(17.14,2);

    gsWarn<<"Obtain E and nu\n";
    // Define parameters vector depending on material law
    std::vector<gsFunctionSet<>*> parameters;
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
    else if (material==4) // OG
    {
      parameters.resize(4);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
    }
    else
        GISMO_ERROR("No material law known");
    //! [Make material functions]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    // Make gsMaterialMatrix depending on the user-defined choices

    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
    materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);

    //! [Define jacobian and residual]
    gsThinShellAssemblerBase<real_t>* assembler;
    gsSparseMatrix<> matrix;
    gsVector<> vector;
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;

    gsInfo<<"---------------------------------------------------------------------\n";
    gsInfo<<"-------------------------Stage 1-------------------------------------\n";
    gsInfo<<"---------------------------------------------------------------------\n";
    //! [Set boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsVector<> dzvec(2);
    dzvec<<0,-6e-3;
    // gsConstantFunction<> DZ(-6e-3,2);
    // gsConstantFunction<> DZ(dzvec,2);
    // bc.addCondition(0,boundary::west , condition_type::dirichlet, 0, 0 ,false, -1);
    // bc.addCondition(1,boundary::south, condition_type::weak_dirichlet,&DZ);
    // bc.addCondition(2,boundary::east , condition_type::dirichlet, 0, 0 ,false, -1);

    gsFunctionExpr<> DZ("0","if ((x>=0.4) and (x<=0.6),-6e-3,0)",2);
    bc.addCondition(0,boundary::west , condition_type::dirichlet, 0, 0 ,false, -1);
    bc.addCondition(0,boundary::south, condition_type::weak_dirichlet,&DZ);
    bc.addCondition(0,boundary::east , condition_type::dirichlet, 0, 0 ,false, -1);

    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrix);
    assembler->options().setInt("Continuity",0);


    assembler->assemble();
    matrix = assembler->matrix();
    vector = assembler->rhs();

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      return assembler->rhs();
    };

    //! [Solve linear problem]
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
    //! [Solve non-linear problem]

    gsMultiPatch<> deformation2 = assembler->constructDisplacement(solVector);
    //! [Define jacobian and residual]
    mp_def = assembler->constructSolution(solVector);
    gsWriteParaview(mp_def,"mp_def",100,true);

    gsInfo <<"Maximum deformation coef: "
           << deformation2.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation2.patch(0).coefs().colwise().minCoeff() <<".\n";
    //! [Construct and evaluate solution]
    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation2);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "FinalDeformation", 10000, true);
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

        gsPiecewiseFunction<> pstrain_m;
        assembler->constructStress(mp_def,pstrain_m,stress_type::principal_membrane_strain);
        gsField<> pstrainM(mp_def,pstrain_m, true);

        gsPiecewiseFunction<> pstrain_f;
        assembler->constructStress(mp_def,pstrain_f,stress_type::principal_flexural_strain);
        gsField<> pstrainF(mp_def,pstrain_f, true);

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

        gsPiecewiseFunction<> VMStresses;
        assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
        gsField<> VMStress(mp_def,VMStresses, true);


        gsPiecewiseFunction<> TFes;
        assembler->constructStress(mp_def,TFes,stress_type::tension_field);
        gsField<> TF(mp_def,TFes, true);


        gsWriteParaview(membraneStress,"MembraneStress",5000);
        // gsWriteParaview(VMStress,"MembraneStressVM",5000);
        // gsWriteParaview(Stretches,"PrincipalStretch",5000);
        gsWriteParaview(pstrainM,"PrincipalMembraneStrain",5000);
        // gsWriteParaview(pstrainF,"PrincipalFlexuralStrain",5000);
        gsWriteParaview(pstressM,"PrincipalMembraneStress",5000);
        // gsWriteParaview(pstressF,"PrincipalFlexuralStress",5000);
        // gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
        // gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
        // gsWriteParaview(stretchDir2,"PrincipalDirection2",5000);
        // gsWriteParaview(stretchDir3,"PrincipalDirection3",5000);
        gsWriteParaview(TF,"tensionfield",5000);


    }

    if (testCase==2)
    {
        gsPiecewiseFunction<> VMStresses;
        assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
        gsField<> VMStress(mp_def,VMStresses, true);
        gsVector<> pt(2);
        pt<<0,0;
        gsInfo<<"Stress in corner: "<<VMStresses.piece(0).eval(pt)<<"\n";
    }
    // ! [Export visualization in ParaView]

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    return 0;

    gsMatrix<> z(1,1);
    z.setZero();

    gsVector<> pt(2);
    pt<<0.5,0.5;

    gsMaterialMatrixLinear<2,real_t> * materialMatrixLinear = new gsMaterialMatrixLinear<2,real_t>(mp,t,parameters,rho);
    gsMaterialMatrixTFT<2,real_t> * materialMatrixTFT = new gsMaterialMatrixTFT<2,real_t>(materialMatrixLinear);

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> matA(materialMatrixTFT,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> matB(materialMatrixTFT,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> matC(materialMatrixTFT,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> matD(materialMatrixTFT,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> vecN(materialMatrixTFT,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> vecM(materialMatrixTFT,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::PStrainN> pstrain(materialMatrix,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::PStressN> pstress(materialMatrix,&mp_def);

    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixA> matAL(materialMatrixLinear,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixB> matBL(materialMatrixLinear,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixC> matCL(materialMatrixLinear,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::MatrixD> matDL(materialMatrixLinear,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> vecNL(materialMatrixLinear,&mp_def);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> vecML(materialMatrixLinear,&mp_def);


    gsMaterialMatrixEval<real_t,MaterialOutput::TensionField> tensionfield(materialMatrixLinear,&mp_def,z);


    gsMaterialMatrixEval<real_t,MaterialOutput::MatrixA> mat(materialMatrixLinear,&mp_def,z);
    gsMaterialMatrixEval<real_t,MaterialOutput::StretchDir> pdir(materialMatrixLinear,&mp_def,z);

    gsDebugVar(mat.piece(0).eval(pt));
    gsDebugVar(matA.piece(0).eval(pt));
    gsDebugVar(matAL.piece(0).eval(pt));

    gsDebugVar(matB.piece(0).eval(pt));
    gsDebugVar(matBL.piece(0).eval(pt));

    gsDebugVar(matC.piece(0).eval(pt));
    gsDebugVar(matCL.piece(0).eval(pt));

    gsDebugVar(matD.piece(0).eval(pt));
    gsDebugVar(matDL.piece(0).eval(pt));

    gsDebugVar(matD.piece(0).eval(pt));
    gsDebugVar(matDL.piece(0).eval(pt));

    gsDebugVar(vecN.piece(0).eval(pt));
    gsDebugVar(vecNL.piece(0).eval(pt));

    gsDebugVar(vecM.piece(0).eval(pt));
    gsDebugVar(vecML.piece(0).eval(pt));


    gsField<> tensionField(mp_def, tensionfield, true);
    gsInfo<<"Plotting in Paraview...\n";
    gsWriteParaview<>( tensionField, "TensionField", 10000, true);

    // gsThinShellAssemblerBase<real_t> * assembler2;
    // gsThinShellAssemblerBase<real_t> * assembler3;
    // assembler2 = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrixLinear);
    // assembler2->assemble();
    // assembler3 = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrixTFT);
    // assembler3->assemble();

    // assembler->assemble();
    // gsDebugVar(assembler->matrix().toDense());
    // gsDebugVar(assembler2->matrix().toDense());
    // gsDebugVar(assembler3->matrix().toDense());

    delete assembler, assembler;
    return EXIT_SUCCESS;
}// end main
