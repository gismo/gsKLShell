/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/gsMPBESBasis.h>
#include <gsUnstructuredSplines/gsMPBESSpline.h>
#include <gsUnstructuredSplines/gsDPatch.h>
#include <gsUnstructuredSplines/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/gsC1SurfSpline.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

using namespace gismo;


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool mesh       = false;
    bool stress     = false;
    bool write      = false;
    bool last       = false;
    bool nonlinear  = false;

    index_t numRefine  = 0;
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

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "DirBc", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "ClaBc", "Clamped BC penalty scalar",  bcClamped );

    cmd.addReal( "d", "DirIfc", "Dirichlet penalty scalar",  ifcDirichlet );
    cmd.addReal( "c", "ClaIfc", "Clamped penalty scalar",  ifcClamped );

    cmd.addString( "G", "geom","File containing the geometry",  fn1 );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  fn2 );
    cmd.addString( "O", "opt", "File containing solver options",  fn3 );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch("last", "last case only",last);
    cmd.addSwitch( "nl", "Print information", nonlinear );

    //! [Parse command line]
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Define material parameters and geometry per example]
    gsMultiPatch<> mp;
    gsBoundaryConditions<> bc;

    gsFileData<> fd;
    gsInfo<<"Reading geometry from "<<fn1<<"...";
    gsReadFile<>(fn1, mp);
    if (mp.nInterfaces()==0 && mp.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it...";
        mp.computeTopology();
    }
    gsInfo<<"Finished\n";

    fd.read(fn2);
    index_t num = 0;
    gsInfo<<"Reading BCs from "<<fn2<<"...";
    num = fd.template count<gsBoundaryConditions<>>();
    GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
    gsInfo<<"Finished\n";

    bc.setGeoMap(mp);

    // Loads
    gsFunctionExpr<> force, pressure;
    gsInfo<<"Reading force function from "<<fn2<<" (ID=21) ...";
    fd.getId(21, force); // id=1: source function
    gsInfo<<"Finished\n";
    // fd.getId(22, pressure); // id=1: source function ------- TO DO!
    // gsInfo<<"Pressure function "<< force << "\n";

    // Loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;
    if ( fd.hasId(30) )
        fd.getId(30,points);
    if ( fd.hasId(31) )
        fd.getId(31,loads);

    if ( fd.hasId(32) )
        fd.getId(32,pid_ploads);
    else
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    gsInfo<<pLoads;

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<fn2<<" (ID=50) ...";
    if ( fd.hasId(50) )
        fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches from "<<fn2<<" (ID=51) ...";
    if ( fd.hasId(51) )
        fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values from "<<fn2<<" (ID=52) ...";
    if ( fd.hasId(52) )
        fd.getId(52,refValue);
    else
        refValue = gsMatrix<>::Zero(mp.geoDim(),refPoints.cols());
    gsInfo<<"Finished\n";
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");

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

    if (mp.geoDim()==2)
        mp.embed(3);

    // p-refine
    GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    mp.degreeElevate(degree-mp.patch(0).degree(0));

        // h-refine each basis
        for (int r =0; r < 2; ++r)
        {
            mp.uniformRefine(1,degree-smoothness);
        }
        numRefine -= 2;
    if (last)
    {
        // h-refine each basis
        for (int r =0; r < numRefine; ++r)
        {
            mp.uniformRefine(1,degree-smoothness);
        }
        numRefine = 0;
    }

    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    for (size_t p=0; p!=mp.nPatches(); p++)
        gsInfo <<"Basis "<<p<<": "<< dbasis.basis(0)<<"\n";

    if (plot) gsWriteParaview(mp,"mp",1000,mesh,false);
    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    //     gsDebugVar(mp.patch(p));

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
    gsMatrix<> refs(numRefine+1,3*refPoints.cols());

    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> solVector;

    gsField<> solField;
    gsMultiPatch<> mp_def, deformation;
    for( index_t r = 0; r<=numRefine; ++r)
    {

        gsInfo<<"--------------------------------------------------------------\n";

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

        assembler.assemble();
        //! [Assemble linear part]
        gsSparseMatrix<> matrix = assembler.matrix();
        gsVector<> vector = assembler.rhs();
        //! [Assemble linear part]

        // Linear solve
        gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
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


        //! [Construct and evaluate solution]
        mp_def = assembler.constructSolution(solVector);
        deformation = assembler.constructDisplacement(solVector);
        //! [Construct and evaluate solution]

        solField = gsField<>(mp_def, deformation);
        if (refPoints.cols()!=0)
        {
            for (index_t p=0; p!=refPoints.cols(); p++)
                // refs.block(r,p*mp.geoDim(),1,mp.geoDim()) = def.piece(refPatches(0,p)).eval(refPoints.col(p)).transpose();
                refs.block(r,p*mp.geoDim(),1,mp.geoDim()) = solField.value(refPoints.col(p),refPatches(0,p)).transpose();
        }

        numDofs[r] = assembler.numDofs();

        mp.uniformRefine();
        dbasis = gsMultiBasis<>(mp);
    }
    //! [Solver loop]

    gsInfo<<"numDoFs";
    for (index_t p=0; p!=refPoints.cols(); ++p)
        gsInfo<<",x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p);
    gsInfo<<"\n";

    for (index_t k=0; k<=numRefine; ++k)
    {
        gsInfo<<numDofs(k);
        for (index_t p=0; p!=refPoints.cols(); ++p)
        {
            gsInfo<<","<<refs(k,3*p)<<","<<refs(k,3*p+1)<<","<<refs(k,3*p+2)<<",";
        }
        gsInfo<<"\n";
    }

        // ! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField, "Deformation", 1000, true);
        gsWriteParaview<>( solField, "Deformation", 1000, false);

    }
    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        assembler.constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler.constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsWriteParaview(membraneStress,"MembraneStress",1000);
        gsWriteParaview(flexuralStress,"FlexuralStress",1000);
    }
    return EXIT_SUCCESS;
}
