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

#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>

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
  index_t Compressibility = 0;
  index_t material = 0;
  bool verbose = false;
  std::string fn;

  bool SingularPoint= false;
  bool quasiNewton  = false;
  int quasiNewtonInt= -1;
  bool adaptive     = false;
  int step          = 10;
  int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)

  index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral
  bool TFT = false;

  // Arc length method options
  real_t dL         = 0; // General arc length
  real_t dLb        = 0.5; // Ard length to find bifurcation
  real_t tol        = 1e-3;
  real_t tolU       = 1e-3;
  real_t tolF       = 1e-3;
  real_t tau        = 1e4;
  real_t relax      = 1.0;
  index_t maxit     = 20;

  std::string assemblerOptionsFile = "options/assembler_options.xml";

  gsCmdLine cmd("2D shell example.");
  cmd.addInt( "e", "degreeElevation",
    "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
  cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
  cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
  cmd.addInt( "M", "Material", "Material law",  material );
  cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
  cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
  cmd.addSwitch("TFT", "Use Tension Field Theory", TFT);

  cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
  cmd.addReal("L","dLb", "arc length", dLb);
  cmd.addReal("l","dL", "arc length after bifurcation", dL);
  cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

  cmd.addReal("F","factor", "factor for bifurcation perturbation", tau);
  cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
  cmd.addInt("N", "maxsteps", "Maximum number of steps", step);


  cmd.addString( "f", "file", "Input XML file", fn );
  cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
  cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
  cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  //! [Parse command line]

  std::string output = "solution";
  std::string dirname = "ArcLengthResults";

  gsOptionList assemblerOptions;
  gsReadFile<>(assemblerOptionsFile,assemblerOptions);

  real_t aDim,bDim;
  real_t thickness = 0.14e-3;
  real_t E_modulus     = 1;
  real_t PoissonRatio = 0;
  real_t Density = 1e0;
  real_t Ratio = 7.0;

  if ((!Compressibility) && (material!=0))
    PoissonRatio = 0.5;
  else
    PoissonRatio = 0.499;

  real_t mu, C01,C10;
  if (material==3)
  {
    C10 = 6.21485502e4; // c1/2
    C01 = 15.8114570e4; // c2/2
    Ratio = C10/C01;
    mu = 2*(C01+C10);
  }
  else
  {
    C10 = 19.1010178e4;
    mu = 2*C10;
  }
  E_modulus = 2*mu*(1+PoissonRatio);
  gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"; mu = "<<mu<<"; ratio = "<<Ratio<<"\n";

  //! [Read input file]
  gsMultiPatch<> mp;
  gsMultiPatch<> mp_def;

  bool nonlinear = true;
  real_t length = 2;
  real_t width = 1;

  mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
  mp.patch(0).coefs().col(0) *= length;
  mp.patch(0).coefs().col(1) *= width;
  mp.addAutoBoundaries();
  mp.computeTopology();


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
  gsConstantFunction<> t(thickness,2);
  gsConstantFunction<> E(E_modulus,2);
  gsConstantFunction<> nu(PoissonRatio,2);
  gsConstantFunction<> rho(Density,2);
  gsConstantFunction<> ratio(Ratio,2);

  std::vector<gsFunction<>*> parameters;
  if (material==0) // SvK
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

  gsMaterialMatrixBase<real_t>* materialMatrix;
  gsOptionList options;
  if      (material==0 && impl==1)
  {
    parameters.resize(2);
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
  }
  else
  {
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
    materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
  }

  //! [Define jacobian and residual]
  gsStopwatch stopwatch,stopwatch2;
  real_t time = 0.0;
  real_t totaltime = 0.0;

  gsInfo<<"---------------------------------------------------------------------\n";
  gsInfo<<"-------------------------Stage 1-------------------------------------\n";
  gsInfo<<"---------------------------------------------------------------------\n";
  //! [Set boundary conditions]
  gsBoundaryConditions<> bc;
  bc.setGeoMap(mp);
  gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

  bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,0);

  bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
  bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);

  bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.

  real_t Load = 2*1e0;
  gsVector<> point(2); point<< 1.0, 0.5 ;
  gsVector<> load (2); load << Load,0.0;
  pLoads.addLoad(point, load, 0 );

  // Construct the gsThinShellAssembler
  gsThinShellAssemblerBase<real_t>* assembler;
  gsMaterialMatrixTFT<2,real_t,false> materialMatrixTFT(mp,t,materialMatrix);
  materialMatrixTFT.options().setReal("SlackMultiplier",1e-2);
  materialMatrixTFT.options().setSwitch("Explicit",false);
  materialMatrixTFT.updateDeformed(&mp_def);

  if (TFT)
    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,&materialMatrixTFT);
  else
    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,bc,force,materialMatrix);

  assembler->options().setInt("Continuity",0);
  assembler->setPointLoads(pLoads);
  assembler->setOptions(assemblerOptions);
  //! [Make assembler]

  // Function for the Jacobian
  typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
  typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
  Jacobian_t Jacobian = [&time,&stopwatch,&assembler](gsVector<real_t> const &x)
  {
    gsMultiPatch<> tmp;
    stopwatch.restart();
    assembler->constructSolution(x,tmp);
    assembler->assembleMatrix(tmp);
    time += stopwatch.stop();
    gsSparseMatrix<real_t> m = assembler->matrix();
    return m;
  };
  // Function for the Residual
  ALResidual_t ALResidual = [&time,&stopwatch,&assembler](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
  {
    gsMultiPatch<> tmp;
    stopwatch.restart();
    assembler->constructSolution(x,tmp);
    assembler->assembleVector(tmp);
    gsVector<real_t> Fint = -(assembler->rhs() - force);
    gsVector<real_t> result = Fint - lam * force;
    time += stopwatch.stop();
    return result; // - lam * force;
  };
  //! [Define jacobian and residual]

  stopwatch.restart();
  stopwatch2.restart();
  assembler->assemble();
  //! [Assemble linear part]
  gsVector<> Force = assembler->rhs();
  //! [Assemble linear part]


  gsALMBase<real_t> * arcLength;
  if (method==0)
    arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
  else if (method==1)
    arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
  else if (method==2)
    arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
  else
    GISMO_ERROR("Method "<<method<<" unknown");

  arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
  arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
  arcLength->options().setReal("Length",dLb);
  arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
  arcLength->options().setSwitch("AdaptiveLength",adaptive);
  arcLength->options().setInt("AdaptiveIterations",5);
  arcLength->options().setReal("Perturbation",tau);
  arcLength->options().setReal("Scaling",0.0);
  arcLength->options().setReal("Tol",tol);
  arcLength->options().setReal("TolU",tolU);
  arcLength->options().setReal("TolF",tolF);
  arcLength->options().setInt("MaxIter",maxit);
  arcLength->options().setSwitch("Verbose",true);
  arcLength->options().setReal("Relaxation",relax);
  if (quasiNewtonInt>0)
  {
    quasiNewton = true;
    arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
  }
  arcLength->options().setSwitch("Quasi",quasiNewton);


  gsDebug<<arcLength->options();
  arcLength->applyOptions();
  arcLength->initialize();

  gsParaviewCollection collection(dirname + "/" + output);
  gsParaviewCollection Smembrane(dirname + "/" + "membrane");
  gsParaviewCollection Sflexural(dirname + "/" + "flexural");
  gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
  gsParaviewCollection Stensionfield(dirname + "/" + "tensionfield");
  gsMultiPatch<> deformation = mp;

  // Make objects for previous solutions
  real_t Lold = 0;
  gsMatrix<> Uold = Force;
  Uold.setZero();

  gsMatrix<> solVector;
  real_t indicator = 0.0;
  arcLength->setIndicator(indicator); // RESET INDICATOR
  bool bisected = false;
  real_t dLb0 = dLb;
  for (index_t k=0; k<step; k++)
  {
    gsInfo<<"Load step "<< k<<"\n";
  // assembler->constructSolution(solVector,solution);
    materialMatrixTFT.updateDeformed(&mp_def);
    arcLength->step();

  // gsInfo<<"m_U = "<<arcLength->solutionU()<<"\n";
    if (!(arcLength->converged()))
    {
      gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
      dLb = dLb / 2.;
      arcLength->setLength(dLb);
      arcLength->setSolution(Uold,Lold);
      bisected = true;
      k -= 1;
      continue;
  // if (plot)
  // {
  //   solVector = arcLength->solutionU();
  //   Uold = solVector;
  //   Lold = arcLength->solutionL();
  //   assembler->constructSolution(solVector,mp_def);

  //   deformation = mp_def;
  //   deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

  //   gsField<> solField(mp,deformation);
  //   std::string fileName = dirname + "/" + output + util::to_string(k);
  //   gsWriteParaview<>(solField, fileName, 5000);
  //   fileName = output + util::to_string(k) + "0";
  //   collection.addTimestep(fileName,k,".vts");
  // }
  // break;
    }

    if (SingularPoint)
    {
      arcLength->computeStability(arcLength->solutionU(),quasiNewton);
      if (arcLength->stabilityChange())
      {
        gsInfo<<"Bifurcation spotted!"<<"\n";
        arcLength->computeSingularPoint(1e-4, 5, Uold, Lold, 1e-10, 0, false);
        arcLength->switchBranch();
        dLb0 = dLb = dL;
        arcLength->setLength(dLb);
      }
    }
    indicator = arcLength->indicator();

    solVector = arcLength->solutionU();
    Uold = solVector;
    Lold = arcLength->solutionL();
    assembler->constructSolution(solVector,mp_def);

    if (testCase==4 || testCase==8 || testCase==9)
    {
      gsMatrix<> pts(2,1);
      pts<<0.5,0.5;
      if (testCase==8 || testCase==9)
      {
        pts.resize(2,3);
        pts.col(0)<<0.0,1.0;
        pts.col(1)<<0.5,1.0;
        pts.col(2)<<1.0,1.0;
      }

      gsMatrix<> lambdas = assembler->computePrincipalStretches(pts,mp_def,0);
      std::streamsize ss = std::cout.precision();
      std::cout <<std::setprecision(20)
      <<"lambdas = \n"<<lambdas<<"\n";
      std::cout<<std::setprecision(ss);

      real_t S = Lold / 1e-3 / lambdas(0) / lambdas(2);
      real_t San = mu * (math::pow(lambdas(1),2)-1/lambdas(1));
      gsDebugVar(S);
      gsDebugVar(San);
      gsDebugVar(abs(S-San));
    }

    deformation = mp_def;
  deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

  // gsDebugVar(mp_def.patch(0).coefs());

  if (plot)
  {
    gsField<> solField(mp,deformation);

    std::string fileName = dirname + "/" + output + util::to_string(k);
    gsWriteParaview<>(solField, fileName, 1000,true);
    fileName = output + util::to_string(k) + "0";
    collection.addPart(fileName + ".vts",k);
    collection.addPart(fileName + "_mesh.vtp",k);
  }
  if (stress)
  {
    std::string fileName;

    gsPiecewiseFunction<> membraneStresses, flexuralStresses, membraneStresses_p, TensionFields;
    gsField<> membraneStress, flexuralStress, membraneStress_p, TensionField;

    ////////////////////////////////////////////////////////////////////////
    assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
    membraneStress = gsField<>(mp,membraneStresses,true);

    fileName = dirname + "/" + "membrane" + util::to_string(k);
    gsWriteParaview( membraneStress, fileName, 1000);
    fileName = "membrane" + util::to_string(k) + "0";
    Smembrane.addPart(fileName + ".vts",k);

    ////////////////////////////////////////////////////////////////////////
    assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
    flexuralStress = gsField<>(mp,flexuralStresses, true);

    fileName = dirname + "/" + "flexural" + util::to_string(k);
    gsWriteParaview( flexuralStress, fileName, 1000);
    fileName = "flexural" + util::to_string(k) + "0";
    Sflexural.addPart(fileName + ".vts",k);

    // ////////////////////////////////////////////////////////////////////////
    // assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
    // membraneStress_p = gsField<>(mp,membraneStresses_p, true);

    // fileName = dirname + "/" + "membrane_p" + util::to_string(k);
    // gsWriteParaview( membraneStress_p, fileName, 1000);
    // fileName = "membrane_p" + util::to_string(k) + "0";
    // Smembrane_p.addPart(fileName + ".vts",k);

    ////////////////////////////////////////////////////////////////////////
    assembler->constructStress(mp_def,TensionFields,stress_type::tension_field);
    TensionField = gsField<>(mp,TensionFields, true);

    fileName = dirname + "/" + "tensionfield" + util::to_string(k);
    gsWriteParaview( TensionField, fileName, 1000);
    fileName = "tensionfield" + util::to_string(k) + "0";
    Stensionfield.addPart(fileName + ".vts",k);
  }

  if (!bisected)
  {
    dLb = dLb0;
    arcLength->setLength(dLb);
  }
  bisected = false;

  }

  if (plot)
  {
    collection.save();
    Smembrane.save();
    Sflexural.save();
    Smembrane_p.save();
    Stensionfield.save();
  }

  delete materialMatrix;
  delete assembler;
  delete arcLength;

  return EXIT_SUCCESS;

}// end main
