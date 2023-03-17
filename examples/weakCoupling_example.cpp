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
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsThinShellUtils.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool stress     = false;
    bool nonlinear  = false;
    bool homogeneous = false;

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
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch( "nl", "Print information", nonlinear );
    cmd.addSwitch("homogeneous", "homogeneous dirichlet BCs",homogeneous);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Define material parameters and geometry per example]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
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
    if (mp.geoDim()==2)
        mp.embed(3);

    if (homogeneous)
    {
        for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet, 0, false, 0, -1);
    }
    else
    {
        index_t num = 0;
        gsInfo<<"Reading BCs from "<<fn2<<"...";
        num = fd.template count<gsBoundaryConditions<>>();
        GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
        fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
        gsInfo<<"Finished\n";

    }


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
    gsMatrix<> refPoints, refPars, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<fn2<<" (ID=50) ...";
    if ( fd.hasId(50) )
        fd.getId(50,refPoints);
    if (refPoints.rows()==2)
    {
        refPars = refPoints;
        gsInfo<<"Reference points are provided in parametric coordinates.\n";
    }
    else if (refPoints.rows()==3)
        gsInfo<<"Reference points are provided in physical coordinates.\n";
    else
        gsInfo<<"No reference points are provided.\n";
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

    //! [Refine and elevate]
    // p-refine
    mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
    {
        mp.uniformRefine(1,degree-smoothness);
    }

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);
    //! [Refine and elevate]

    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    for (size_t p=0; p!=mp.nPatches(); p++)
        gsInfo <<"Basis "<<p<<": "<< dbasis.basis(0)<<"\n";

    //! [Make assembler]
    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters,rho);

    // Construct the gsThinShellAssembler
    gsThinShellAssembler<3, real_t, true> assembler(mp,dbasis,bc,force,&materialMatrix);

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
    //! [Make assembler]

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
      assembler.constructSolution(x,mp_def);
      assembler.assembleMatrix(mp_def);
      time += stopwatch.stop();
      gsSparseMatrix<real_t> m = assembler.matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler.constructSolution(x,mp_def);
      assembler.assembleVector(mp_def);
      time += stopwatch.stop();
      return assembler.rhs();
    };
    //! [Define jacobian and residual]

    stopwatch.restart();
    stopwatch2.restart();
    assembler.assemble();
    time += stopwatch.stop();

    //! [Assemble linear part]
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

    totaltime += stopwatch2.stop();

    //! [Construct and evaluate solution]
    mp_def = assembler.constructSolution(solVector);
    gsMultiPatch<> deformation = assembler.constructDisplacement(solVector);
    //! [Construct and evaluate solution]

    gsExprEvaluator<> ev;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    geometryMap m_ori   = ev.getMap(mp);
    geometryMap m_def   = ev.getMap(mp_def);

    gsMatrix<> pts(2,4);
    pts.col(0)<<0.1,0.1;
    pts.col(1)<<0.1,0.9;
    pts.col(2)<<0.9,0.1;
    pts.col(3)<<0.9,0.9;
    auto test = (flat(jac(m_def).tr()*jac(m_def))).tr();
    auto That   = cartcon(m_ori);
    auto Ttilde = cartcov(m_ori);
    auto cartJac_ori = That * jac(m_ori);
    auto cartJac_def = That * jac(m_def);
    auto E_m    = 0.5 * ( flat(cartJac_def.tr()*cartJac_def) - flat(cartJac_ori.tr()*cartJac_ori) );
    // auto E_m    = 0.5 * ( flat(jac(m_def).tr()*jac(m_def)) - flat(jac(m_ori).tr()* jac(m_ori)) );

    ev.setIntegrationElements(dbasis);
    ev.writeParaview(E_m,m_ori,"test1");
    ev.writeParaview(flat(cartJac_ori.tr()*cartJac_ori) * meas(m_ori),m_ori,"test2");
    ev.writeParaview(flat(cartJac_def.tr()*cartJac_def) * meas(m_ori),m_ori,"test3");

    ev.writeParaview(flat(jac(m_ori).tr()*jac(m_ori)) * meas(m_ori),m_ori,"test4");
    ev.writeParaview(flat(jac(m_def).tr()*jac(m_def)) * meas(m_ori),m_ori,"test5");


    //! [Construct and evaluate solution]
    gsField<> solField(mp_def, deformation);
    if (refPoints.cols()!=0)
    {
        gsMatrix<> result;
        if (refPoints.rows()==3) // then they are provided in the physical domain and should be mapped to the parametric domain
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

        gsInfo<<"Physical coordinates of points\n";
        for (index_t p=0; p!=refPars.cols(); p++)
        {
            gsInfo<<",x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p);
        }
        gsInfo<<"\n";

        for (index_t p=0; p!=refPars.cols(); ++p)
        {
            mp.patch(refPatches(0,p)).eval_into(refPars.col(p),result);
            gsInfo<<result.row(0)<<","<<result.row(1)<<","<<result.row(2)<<",";
        }
        gsInfo<<"\n";
        gsMatrix<> refs(1,mp.geoDim()*refPoints.cols());
        for (index_t p=0; p!=refPars.cols(); p++)
            refs.block(0,p*mp.geoDim(),1,mp.geoDim()) = mp.piece(refPatches(0,p)).eval(refPars.col(p)).transpose();
        gsInfo<<"Reference point coordinates\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<refs(0,mp.geoDim()*p)<<"\t"<<refs(0,mp.geoDim()*p+1)<<"\t"<<refs(0,mp.geoDim()*p+2)<<"\t";
        gsInfo<<"\n";

        for (index_t p=0; p!=refPars.cols(); p++)
            refs.block(0,p*mp.geoDim(),1,mp.geoDim()) = solField.value(refPars.col(p),refPatches(0,p)).transpose();
        gsInfo<<"Computed values\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<refs(0,mp.geoDim()*p)<<"\t"<<refs(0,mp.geoDim()*p+1)<<"\t"<<refs(0,mp.geoDim()*p+2)<<"\t";
        gsInfo<<"\n";

        gsInfo<<"Reference values\n"; // provided as mp.geoDim() x points.cols() matrix
        for (index_t p=0; p!=refValue.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refValue.cols(); ++p)
            for (index_t d=0; d!=mp.geoDim(); d++)
                gsInfo<<refValue(d,p)<<"\t";
        gsInfo<<"\n";
    }
    // ! [Export visualization in ParaView]
    if (plot)
    {
        // gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField, "Deformation", 1000, true);
        gsWriteParaview<>( solField, "Deformation", 1000, false);

    }
    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        gsDebugVar("MembraneStress construction");
        assembler.constructStress(mp,mp_def,membraneStresses,stress_type::membrane);
        gsWriteParaview(mp,membraneStresses,"MembraneStress",5000);

        gsPiecewiseFunction<> membraneStressesVM;
        gsDebugVar("MembraneStress (VM) construction");
        assembler.constructStress(mp,mp_def,membraneStressesVM,stress_type::von_mises_membrane);
        gsWriteParaview(mp,membraneStressesVM,"MembraneStressVM",5000);

        gsPiecewiseFunction<> flexuralStresses;
        gsDebugVar("FlexuralStress construction");
        assembler.constructStress(mp,mp_def,flexuralStresses,stress_type::flexural);
        gsWriteParaview(mp,flexuralStresses,"FlexuralStress",5000);
    }
    // ! [Export visualization in ParaView]

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    return EXIT_SUCCESS;

}// end main
