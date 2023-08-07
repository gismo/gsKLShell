/** @file example_pointload.cpp

    @brief Asymmetric loading applied with an offset angle (pi/50)
    from the crown of a semicircular arch.

    Example 13 from Liu et al 2006

    Liu, B., Xing, Y., Wang, Z., Lu, X., & Sun, H. (2017) Non-uniform rational Lagrange functions and its applications to isogeometric analysis of in-plane and flexural vibration of thin plates. Computer Methods in Applied Mechanics and Engineering, 321, 173-208. http://dx.doi.org/10.1016/j.cma.2017.04.007

    Author(s): J. Li
 **/

//![Include namespace]
#include <gismo.h>
#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

using namespace gismo;

int main(int argc, char *argv[]){

    //! [Parse command line]
    bool plot  = true;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool Compressibility = false;
    index_t material = 0;
    bool verbose = false;
    std::string fn;
    bool membrane = false;
    int steps = 100;

    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

//    real_t Ratio = 7.0;

    gsCmdLine cmd("2D open sectiorial membrane.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
//    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
//    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("comp", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("composite", "Composite material", composite);


    std::string assemberOptionsFile("options/solver_options.xml");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);
    //! [Parse command line]

    //! [Material properties]
    real_t E_modulus = 7.1e10; //Pa
    real_t PoissonRatio = 0.3;
    real_t Density = 2700.0; //kg/m^3
    real_t thickness = 0.005; //mm
    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
//    real_t omega = 6.8138; // omega2 = 8.2666 omega3 = 12.855 omega4 = 13.742
//    gsInfo<<"Eigen frequency: "<< omega;
    //! [Material properties]

    //! [Make geometry, refine and evaluate]
    gsMultiPatch<> mp, mp_vib;
    real_t r0 = 0.5; //m
    real_t r1 = 1.0;
//    gsTensorNurbs<2> geo = *gsNurbsCreator<>::NurbsQuarterAnnulus(r0,r1);
    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(r0,r1));
//    mp.embed(3);
    mp.addAutoBoundaries();
    real_t PI = 3.1415926535;
    real_t Area = 0.25*PI*(pow(r1,2) - pow(r0,2));

    real_t EA = E_modulus*Area;
//    real_t EI = 1.0/12.0*(width*math::pow(thickness*2,3))*E_modulus;



    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();
    // Set the deformed configuration
    mp_vib = mp;
    if (plot) gsWriteParaview<>( mp_vib  , "mp", 1000, true);
    gsMultiBasis<> dbasis(mp);


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";
    //! [Make geometry, refine and evaluate]

    //! [Set boundary conditions]
    std::vector< std::pair<patchSide,int> > clamped;
    real_t Pressure_outer = 0.00001;
    real_t Pressure_inner = 0.01;
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    gsFunctionExpr<> neu_outer, neu_inner;
    neu_outer = gsFunctionExpr<>(std::to_string(Pressure_outer/r1) + "*(x + y)",2);
    neu_inner = gsFunctionExpr<>(std::to_string(Pressure_inner/r0) + "*(x + y)",2);


    BCs.addCondition(0,boundary::north,condition_type::dirichlet, 0, 0, false, 0 );
    BCs.addCondition(0,boundary::east, condition_type::dirichlet, &neu_outer );
    BCs.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 0 );
    BCs.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 1 );
    BCs.addCondition(0,boundary::west, condition_type::dirichlet, &neu_inner );
    BCs.addCondition(0,boundary::west, condition_type::clamped, 0, 0, false, 0 );
    BCs.addCondition(0,boundary::west, condition_type::clamped, 0, 0, false, 1 );
    BCs.addCondition(0,boundary::south,condition_type::dirichlet, 0, 0, false, 1 );

    //! [Set boundary conditions]

    //! [Make material functions]
    gsVector<> tmp(2);
    tmp << 0,0;
    gsFunctionExpr<> force("0","0",2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> rho(std::to_string(Density),2);
    gsConstantFunction<> nu(PoissonRatio,2);
    gsFunctionExpr<> t(std::to_string(thickness), 2);




    index_t kmax = 1; // number of layers
    std::vector<gsFunctionSet<> * > Gs(kmax); // Material matrices
    std::vector<gsFunctionSet<> * > Ts(kmax); // Thickness per layer
    std::vector<gsFunctionSet<> * > Phis(kmax); // Fiber angle per layer

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;

    // Define fiber angle
    gsConstantFunction<> phi;
    phi.setValue(0,2);
    Phis[0] = &phi;

    // Define thickness
    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsConstantFunction<> thicks(thickness/kmax,2);
    Ts[0] = &thicks;

    //! [Make material functions]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);



//    gsMaterialMatrixContainer<real_t> materialMats(mp.nPatches());
//    for (size_t p = 0; p!=mp.nPatches(); p++)
//        materialMats.add(materialMatrix);

    //Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2, real_t,false>(mp,dbasis,BCs,force,materialMatrix);
    assembler->setOptions(opts);

    assembler->assemble();
    //! [Make assembler]

    //! [Define jacobian and residual]
    gsStopwatch stopwatch,stopwatch1;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    //Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_vib](gsVector<real_t> const &x)
    {
        ThinShellAssemblerStatus status;
        stopwatch.restart();
        assembler->constructSolution(x,mp_vib);
        status = assembler->assembleMatrix(mp_vib);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        time += stopwatch.stop();
        gsSparseMatrix<real_t> m = assembler->matrix();
        return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_vib](gsVector<real_t> const &x)
    {
        ThinShellAssemblerStatus status;
        stopwatch.restart();
        assembler->constructSolution(x,mp_vib);
        status = assembler->assembleVector(mp_vib);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        time += stopwatch.stop();
        return assembler->rhs();
    };

    //! [Define jacobian and residual]
    ThinShellAssemblerStatus status;
    stopwatch.restart();
    stopwatch1.restart();
    status = assembler->assemble();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    time += stopwatch.stop();

    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    //! [Solve non-linear problem]
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
    //! [Solve non-linear problem]

    totaltime += stopwatch1.stop();

    //! [Construct and evaluate solution]
    mp_vib = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_vib;
    for (size_t k = 0; k != mp_vib.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    //! [Construct and evaluate solution]



// ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_vib, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);
    }
    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_vib,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_vib,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_vib,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_vib,flexuralStresses, true);

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_vib,stretches,stress_type::principal_stretch);
        gsField<> Stretches(mp_vib,stretches, true);

        gsPiecewiseFunction<> pstress_m;
        assembler->constructStress(mp_vib,pstress_m,stress_type::principal_stress_membrane);
        gsField<> pstressM(mp_vib,pstress_m, true);

        gsPiecewiseFunction<> pstress_f;
        assembler->constructStress(mp_vib,pstress_f,stress_type::principal_stress_flexural);
        gsField<> pstressF(mp_vib,pstress_f, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_vib,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_vib,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_vib,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_vib,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_vib,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_vib,stretch3, true);


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
    // ! [Export visualization in ParaView]

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    return EXIT_SUCCESS;

}
