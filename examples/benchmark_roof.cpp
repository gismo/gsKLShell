/** @file example_shell3D.cpp

    @brief Simple 3D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixContainer.h>
#include <gsKLShell/gsMaterialMatrixEval.h>
#include <gsKLShell/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsCore/gsPiecewiseFunction.h>

using namespace gismo;

void writeCSV( const gsMatrix<> & mat, std::string const & filename )
{
    std::string tmp = gsFileManager::getExtension(filename);
    if (tmp != "csv" )
        tmp = filename + ".csv";
    else
        tmp = filename;

    std::ofstream file(tmp.c_str());

    for ( index_t i = 0 ; i != mat.rows(); ++i )
    {
        file << mat(i,0);
        for ( index_t j = 1 ; j != mat.cols(); ++j )
            file << ", " << mat(i,j) ;
        file<< "\n";
    }
    file.close();
}

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool nonlinear = false;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;

    bool NURBS = true;
    bool GAUSS = true;

    gsCmdLine cmd("2D shell example.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 0 = square plate with pressure; 1 = Scordelis Lo Roof; 2 = quarter hemisphere; 3 = pinched cylinder",  testCase );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch( "nl", "Print information", nonlinear );
    cmd.addSwitch( "noNURBS", "Don't use NURBS", NURBS );
    cmd.addSwitch( "noGAUSS", "Don't use GAUSS", GAUSS );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Define material parameters and geometry per example]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    std::string fn;
    if (testCase == 1 )
    {
        thickness = 0.25;
        E_modulus = 4.32E8;
        fn = "surface/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
        PoissonRatio = 0.0;

        if (!NURBS)
        {
            gsTensorNurbsBasis<2,real_t> * basis;
            if (basis = dynamic_cast<gsTensorNurbsBasis<2,real_t> * >(&mp.basis(0)));
            basis->weights().setOnes();
        }
    }
    else
        GISMO_ERROR("TESTCASE UNKNOWN");
    //! [Define material parameters and geometry per example]

    //! [Refine and elevate]
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);
    //! [Refine and elevate]

    gsMultiBasis<> dbasis(mp,!NURBS);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";


    //! [Set boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsPiecewiseFunction<> force(mp.nPatches());
    gsPiecewiseFunction<> t(mp.nPatches());
    // gsPiecewiseFunction<> nu(mp.nPatches());

    gsVector<> tmp(3);
    tmp.setZero();
    // gsMatrix<> refPoint(2,5);
    gsMatrix<> refPoint(2,2);
    real_t refPatch = 0;
    if (testCase == 1)
    {
        tmp << 0, 0, -90;

        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 1 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 2 );

        bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 1 );
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 2 );

        // Surface forces
        gsConstantFunction<> force0(tmp,3);
        force.addPiece(force0);

        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);

        // gsFunctionExpr<> nu0(std::to_string(PoissonRatio), 3);
        // nu.addPiece(nu0);

        // refPoint.col(0)<<0.5,0.0;
        // refPoint.col(1)<<0.5,0.5;
        // refPoint.col(2)<<0.5,1.0;
        // refPoint.col(3)<<0.0,0.5;
        // refPoint.col(4)<<1.0,0.5;
        refPoint.col(0)<<0.5,0.0;
        refPoint.col(1)<<0.5,1.0;
    }
    else
        GISMO_ERROR("TESTCASE UNKNOWN");

    //! [Make material functions]
    // Linear isotropic material model
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    // gsConstantFunction<> t(thickness,3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> nu(PoissonRatio,3);
    // gsFunctionExpr<> force("0","0","0",3);

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,t,E,nu);

    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    // Set the penalty parameter for the interface C1 continuity
    assembler->addWeakC0(mp.topology().interfaces());
    assembler->addWeakC1(mp.topology().interfaces());
    assembler->initInterfaces();
    //! [Make assembler]

    if (!GAUSS)
    {
        assembler->options().setInt ("quRule",3);
        assembler->options().setReal("quA",2*mp.basis(0).maxDegree());
        assembler->options().setInt ("quB",0);
    }

    // // PatchRule
    // gsMatrix<> points;
    // gsVector<> weights;
    // gsMatrix<> TensorPatch(mp.parDim(),0);
    // gsMatrix<> allweights(1,0);
    // index_t start;
    // gsBasis<real_t>::domainIter domIt = mp.basis(0).makeDomainIterator();
    // gsQuadRule<real_t>::uPtr patchRule = gsQuadrature::getPtr(mp.basis(0), assembler->options());
    // for (; domIt->good(); domIt->next() )
    // {
    //     //  Patch-rule
    //     patchRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
    //                     points, weights);
    //     gsInfo  <<"* \t PatchRule \n"
    //             <<"- points:\n"<<points<<"\n"
    //             <<"- weights:\n"<<weights.transpose()<<"\n"
    //             <<"- weights*4:\n"<<weights.transpose().array()*4<<"\n";

    //     start = TensorPatch.cols();
    //     TensorPatch.conservativeResize(Eigen::NoChange,TensorPatch.cols()+points.cols());
    //     TensorPatch.block(0,start,TensorPatch.rows(),points.cols()) = points;

    //     start = allweights.cols();
    //     allweights.conservativeResize(Eigen::NoChange,allweights.cols()+weights.transpose().cols());
    //     allweights.block(0,start,allweights.rows(),points.cols()) = weights.transpose();
    // }

    // gsDebugVar(allweights);
    // gsDebugVar(allweights.sum());

    // gsDebugVar(TensorPatch);

    // Set stopwatch
    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    //! [Define jacobian and residual]
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
    //! [Define jacobian and residual]

    stopwatch.restart();
    stopwatch2.restart();
    assembler->assemble();
    time += stopwatch.stop();

    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    //! [Assemble linear part]

    // writeCSV(matrix.toDense(),"Matrix");
    // gsDebugVar(matrix.toDense());

    // writeCSV(vector,"Vector");
    // gsDebugVar(vector.transpose());

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler->numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    // writeCSV(solVector,"Solution");
    // gsDebugVar(solVector.transpose());

    //! [Solve non-linear problem]
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

    totaltime += stopwatch2.stop();

    //! [Construct and evaluate solution]
    mp_def = assembler->constructSolution(solVector);
    gsMultiPatch<> deformation = assembler->constructDisplacement(solVector);
    //! [Construct and evaluate solution]


    //! [Construct and evaluate solution]
    gsMatrix<> refVals = deformation.patch(refPatch).eval(refPoint);
    // real_t numVal;
    // if      (testCase == 0 || testCase == 1 || testCase == 3)
    //     numVal = refVals.at(2);
    // else
    //     numVal = refVals.at(1);

    // gsInfo << "Displacement at reference point: "<<numVal<<"\n";
    gsInfo << "Displacement at reference point: \n"<<std::setprecision(12)<<refVals<<"\n";
    //! [Construct and evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        // gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField, "Deformation", 1000, true);
        gsWriteParaview<>( solField, "Deformation", 1000, false);

        if (testCase==3)
            gsWarn<<"Deformations are possibly zero in Paraview, due to the default precision (1e-5).\n";
    }
    if (stress)
    {
        if (testCase==2)
        {
            gsWarn<<"Stresses cannot be plotted for this case due to the singularity at the top of the geometry.\n";
        }
        else
        {
            gsPiecewiseFunction<> membraneStresses;
            assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
            gsField<> membraneStress(mp_def,membraneStresses, true);

            gsPiecewiseFunction<> flexuralStresses;
            assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
            gsField<> flexuralStress(mp_def,flexuralStresses, true);

            gsWriteParaview(membraneStress,"MembraneStress",1000);
            gsWriteParaview(flexuralStress,"FlexuralStress",1000);
        }
    }
    // ! [Export visualization in ParaView]

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;

}// end main
