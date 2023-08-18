/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsUtils/gsQuasiInterpolate.h>

#include <gsAssembler/gsExprAssembler.h>

using namespace gismo;

void writeToFile(const std::string & bufferstring, std::ofstream & file, const std::string & name)
{
    file.open(name);
    file<<bufferstring;
    file.close();
    gsInfo<<"Data written to "<<name<<"\n";
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool mesh       = false;
    bool stress     = false;
    bool write      = false;
    bool last       = false;
    bool nonlinear  = false;
    bool homogeneous= false;

    index_t numRefine  = 0;
    index_t numRefine0 = 1;
    index_t degree = 3;
    index_t smoothness = 2;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    std::string fn1,fn2,fn3;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    fn3 = "options/solver_options.xml";

    real_t ifcDirichlet = 1.0;
    real_t ifcClamped = 1.0;

    std::string dirname = ".";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "DirBc", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "ClaBc", "Clamped BC penalty scalar",  bcClamped );

    cmd.addReal( "d", "DirIfc", "Dirichlet penalty scalar",  ifcDirichlet );
    cmd.addReal( "c", "ClaIfc", "Clamped penalty scalar",  ifcClamped );

    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addString( "o", "out", "Dir name of the output",  dirname );

    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "R", "preRefine", "Refinement before the loop.",  numRefine0);

    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("last", "last case only",last);
    cmd.addSwitch( "nl", "Print information", nonlinear );
    cmd.addSwitch("homogeneous", "homogeneous dirichlet BCs",homogeneous);

    //! [Parse command line]
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd;

    // Geometry 
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    gsInfo<<"Reading geometry from "<<fn1<<"...";
    gsReadFile<>(fn1, mp);
    if (mp.nInterfaces()==0 && mp.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it...";
        mp.computeTopology();
    }
    gsInfo<<"Finished\n";
    if (mp.geoDim()==2)
        mp.embed(3);

    // Boundary conditions
    gsBoundaryConditions<> bc;
    if (homogeneous)
    {
        for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet, 0, false, 0, -1);
    }
    else
    {
        fd.read(fn2);
        index_t num = 0;
        gsInfo<<"Reading BCs from "<<fn2<<"...";
        num = fd.template count<gsBoundaryConditions<>>();
        GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
        fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
        gsInfo<<"Finished\n";

    }
    bc.setGeoMap(mp);

    // Distributed load
    gsFunctionExpr<> force;
    gsInfo<<"Reading force function from "<<fn2<<" (ID=21) ...";
    fd.getId(21, force);
    gsInfo<<"Finished\n";

    // Point loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;
    gsInfo<<"Reading point load point locations from "<<fn2<<" (ID=30) ...";
    if ( fd.hasId(30) ) fd.getId(30,points);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load point vectors from "<<fn2<<" (ID=31) ...";
    if ( fd.hasId(31) ) fd.getId(31,loads);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load point patch indices from "<<fn2<<" (ID=32) ...";
    if ( fd.hasId(32) ) fd.getId(32,pid_ploads);
    gsInfo<<"Finished\n";

    if ( !fd.hasId(30) || !fd.hasId(31) || !fd.hasId(32) )
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refPars, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<fn2<<" (ID=50) ...";
    if ( fd.hasId(50) ) fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches from "<<fn2<<" (ID=51) ...";
    if ( fd.hasId(51) ) fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values from "<<fn2<<" (ID=52) ...";
    if ( fd.hasId(52) ) fd.getId(52,refValue);

    if ( !fd.hasId(50) || !fd.hasId(51) || !fd.hasId(52) )
        refValue = gsMatrix<>::Zero(mp.geoDim(),refPoints.cols());
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");

    if (refPoints.rows()==2)
    {
        refPars = refPoints;
        gsInfo<<"Reference points are provided in parametric coordinates.\n";
    }
    else if (refPoints.rows()==3)
        gsInfo<<"Reference points are provided in physical coordinates.\n";
    else
        gsInfo<<"No reference points are provided.\n";

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    gsInfo<<"Reading thickness from "<<fn2<<" (ID=10) ...";
    fd.getId(10,t);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Young's Modulus from "<<fn2<<" (ID=11) ...";
    fd.getId(11,E);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Poisson ratio from "<<fn2<<" (ID=12) ...";
    fd.getId(12,nu);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading density from "<<fn2<<" (ID=13) ...";
    fd.getId(13,rho);
    gsInfo<<"Finished\n";
    
    // Solver options
    gsInfo<<"Reading solver options from "<<fn1<<"...";
    gsOptionList solverOptions;
    fd.read(fn3);
    fd.template getFirst<gsOptionList>(solverOptions);
    gsInfo<<"Finished\n";
    //! [Read input file]

    //! [Create output directory]
    if ((plot || write) && !dirname.empty())
        gsFileManager::mkdir(dirname);
    //! [Create output directory]

    //! [Refine and elevate]
    // p-refine
    GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine each basis
    for (int r =0; r < numRefine0; ++r)
        mp.uniformRefine(1,degree-smoothness);
    numRefine -= numRefine0;

    if (last)
    {
        // h-refine each basis
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine(1,degree-smoothness);
        numRefine = 0;
    }

    mp_def = mp;
    if (plot) gsWriteParaview<>( mp_def    , dirname + "/mp", 1000, mesh);

    // Make multi-basis object
    gsMultiBasis<> dbasis(mp);

    //! [Make assembler]
    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters,rho);

    // Construct the gsThinShellAssembler
    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), linferr(numRefine+1),
        b2err(numRefine+1), b1err(numRefine+1), binferr(numRefine+1);

    gsVector<> numDofs(numRefine+1);
    gsVector<> DisplacementNorm(numRefine+1);
    gsVector<> EnergyNorm(numRefine+1);
    gsMatrix<> refs(numRefine+1,3*refPoints.cols());

    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> solVector;

    gsField<> solField;
    gsMultiPatch<> mp_def, deformation;
    for( index_t r = 0; r<=numRefine; ++r)
    {
        // Construct the gsThinShellAssembler
        assembler = gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,&materialMatrix);

        assembler.options().setReal("WeakDirichlet",bcDirichlet);
        assembler.options().setReal("WeakClamped",bcClamped);
        // Set the penalty parameter for the interface C1 continuity
        assembler.options().setInt("Continuity",-1);
        assembler.options().setReal("IfcDirichlet",ifcDirichlet);
        assembler.options().setReal("IfcClamped",ifcClamped);
        assembler.addWeakC0(mp.topology().interfaces());
        assembler.addWeakC1(mp.topology().interfaces());
        assembler.initInterfaces();

        assembler.setPointLoads(pLoads);

        // Nonlinear
        //! [Define jacobian and residual]
        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
        Jacobian_t Jacobian = [&assembler](gsVector<real_t> const &x)
        {
            gsMultiPatch<> mp_def;
            assembler.constructSolution(x,mp_def);
            assembler.assembleMatrix(mp_def);
            gsSparseMatrix<real_t> m = assembler.matrix();
            return m;
        };
        // Function for the Residual
        Residual_t Residual = [&assembler](gsVector<real_t> const &x)
        {
            gsMultiPatch<> mp_def;
            assembler.constructSolution(x,mp_def);
            assembler.assembleVector(mp_def);
            return assembler.rhs();
        };
        //! [Define jacobian and residual]

        //! [Assemble linear part]
        assembler.assemble();
        gsSparseMatrix<> matrix = assembler.matrix();
        gsVector<> vector = assembler.rhs();
        //! [Assemble linear part]

        //! [Solve linear problem]
        gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
        gsVector<> solVector;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute( matrix );
        solVector = solver.solve(vector);
        //! [Solve linear problem]


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

        //! [Construct and evaluate solution]
        mp_def = assembler.constructSolution(solVector);
        gsMultiPatch<> deformation = assembler.constructDisplacement(solVector);
        //! [Construct and evaluate solution]

        //! [Export reference point data]
        solField = gsField<>(mp_def, deformation);
        if (refPoints.cols()!=0)
        {
            // Reference points are provided in the physical domain and should be mapped to the parametric domain
            if (refPoints.rows()==3) 
            {
                refPars.resize(2,refPoints.cols());
                for (index_t p = 0; p!=refPoints.cols(); p++)
                {
                    mp.patch(refPatches(0,p)).invertPoints(refPoints.col(p),result,1e-10);
                    if (result.at(0)==std::numeric_limits<real_t>::infinity()) // if failed
                        gsWarn<<"Point inversion failed\n";
                    refPars.col(p) = result;
                }
            }

            // Store reference data
            for (index_t p=0; p!=refPars.cols(); p++)
                refs.block(r,p*mp.geoDim(),1,mp.geoDim()) = solField.value(refPars.col(p),refPatches(0,p)).transpose();
        }

        numDofs[r] = assembler.numDofs();
        DisplacementNorm[r] = solVector.transpose() * solVector;
        if (!nonlinear)
            EnergyNorm[r] = solVector.transpose() * matrix * solVector;
        else
            EnergyNorm[r] = solVector.transpose() * Jacobian(solVector) * solVector;

        mp.uniformRefine(1,degree-smoothness);
        dbasis = gsMultiBasis<>(mp);
    }
    //! [Solver loop]

    //! [Export reference point data]
    std::stringstream buffer;
    std::ofstream file;
    buffer<<"numDoFs";
    for (index_t p=0; p!=refPars.cols(); ++p)
        buffer<<",x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p);
    buffer<<",DisplacementNorm,"<<"Energynorm\n";

    for (index_t k=0; k<=numRefine; ++k)
    {
        buffer<<numDofs(k);
        for (index_t p=0; p!=refPars.cols(); ++p)
        {
            buffer<<std::setprecision(12)<<","<<refs(k,3*p)<<","<<refs(k,3*p+1)<<","<<refs(k,3*p+2);
        }
        buffer<<"\n";
    }
    gsInfo<<buffer.str();
    if (write) writeToFile(buffer.str(),file,dirname + "/solutions.csv");
    //! [Export reference point data]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, false);

    }
    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        assembler.constructStress(mp,mp_def,membraneStresses,stress_type::membrane);
        gsWriteParaview(mp,membraneStresses,dirname + "/" + "MembraneStress",5000);

        gsPiecewiseFunction<> membraneStressesVM;
        assembler.constructStress(mp,mp_def,membraneStressesVM,stress_type::von_mises_membrane);
        gsWriteParaview(mp,membraneStressesVM,dirname + "/" + "MembraneStressVM",5000);

        gsPiecewiseFunction<> flexuralStresses;
        assembler.constructStress(mp,mp_def,flexuralStresses,stress_type::flexural);
        gsWriteParaview(mp,flexuralStresses,dirname + "/" + "FlexuralStress",5000);
    }
    // ! [Export visualization in ParaView]
    return EXIT_SUCCESS;
}
