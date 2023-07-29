#include <gismo.h>

#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsThinShellAssembler.h>


#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>

using namespace gismo;


template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsNurbs<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);


int main(int argc, char* argv[])
{
    // Input options
    bool plot = true; // If set to true, paraview file is generated and launched on exit.
    bool xml = false;
    bool stress = true;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 0;
//    bool Compressibility = false;
    bool verbose = false;
    int step          = 100;
    
    // Material properties
    index_t material = 2;
    bool composite = true;
    bool SingularPoint = false;
    bool quasiNewton = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral
    
    real_t E_modulus = 200;
    real_t thickness = 3.464101615;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0; //Couldn't find

    
    
    

    // Read solver option file
    std::string assemberOptionsFile("options/solver_options.xml");

    // Arc length method options
    real_t dL         = -1; // General arc length
    real_t dLb = 0.1; // Arc length to find bifurcation
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;
    
    bool write = true;
    std::string wn("data.csv");
    
    
    // Arc length method options
    gsCmdLine cmd("Arc-length analysis of a semicircular arch with pinned-pinned boundary subjected to an asymmetrical force.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt("t", "testcase",
        "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free",
       testCase);
    cmd.addInt( "m", "Material", "Material law",  material );
//    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);
//
//    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
//    cmd.addSwitch("weak", "Impose boundary conditions weakly", weak);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("xml", "Write geometry into XML files",xml);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    
    //! [Read input file]
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    
    //! [Set material parameters]
    if(testCase ==0)
    {
        thickness = 3.464101615;
        PoissonRatio = 0.499;
    }
    else
        gsInfo<<"Not available";
    //![Make geometry and refine/elevate]
    real_t aDim = 100;
    real_t bDim = 0.2886751345;

// Definition for NurbsHalfAnnulus
//    template<class T> typename gsNurbsCreator<T>::TensorNurbs2Ptr
//    gsNurbsCreator<T>::NurbsHalfAnnulus(T const & r0, T const & r1, T const & x, T const & y)
//    {
//        gsKnotVector<T> KVy (0,1,0,2) ;
//        gsKnotVector<T> KVx (0,1,0,4) ;
//        gsMatrix<T> C(8,2);
//        C<< r0,0,
//            r0,2*r0,
//            -r0,2*r0,
//            -r0,0,
//            r1,0,
//            r1,2*r1,
//            -r1,2*r1,
//            -r1,0;
//
//        float PI = 3.1415926535;
//        C.col(0).array() *= PI;
//        C.col(1).array() *= PI;
//        C.col(0).array() += x;
//        C.col(1).array() += y;
//
//        gsMatrix<T> ww(8,1);
//        ww.setOnes();
//        ww.at(1)= 0.333333333333333;
//        ww.at(2)= 0.333333333333333 ;
//        ww.at(5)= 0.333333333333333 ;
//        ww.at(6)= 0.333333333333333 ;
//        return TensorNurbs2Ptr(new gsTensorNurbs<2,T>(KVx,KVy, give(C), give(ww)));
//    }
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    gsTensorNurbs<2> geo = *gsNurbsCreator<>::NurbsHalfAnnulus((aDim-bDim)/2, aDim/2);
    mp.addPatch(geo);
    
//    template<class T> typename gsNurbsCreator<T>::TensorNurbs3Ptr
//    gsNurbsCreator<T>::Nurbs3DHalfAnnulus(T const & r0, T const r1, T const &x, T const &y, T const &z, T const & thickness)
//    {
//        gsKnotVector<T> KVx (0,1,0,4);
//        gsKnotVector<T> KVy (0,1,0,2);
//        gsKnotVector<T> KVz (0,1,0,2);
//        gsMatrix<T> C(16,3);
//        
//        C<< r0,0, 0,
//            r0,2*r0,0,
//            -r0,2*r0,0,
//            -r0,0,0,
//            r1,0,0,
//            r1,2*r1,0,
//            -r1,2*r1,0,
//            -r1,0,0,
//            r0,0, thickness,
//            r0,2*r0,thickness,
//            -r0,2*r0,thickness,
//            -r0,0,thickness,
//            r1,0,thickness,
//            r1,2*r1,thickness,
//            -r1,2*r1,thickness,
//        -r1,0,thickness;
//        float PI = math::atan(1)*4;
//        C.col(0).array() *= PI;
//        C.col(1).array() *= PI;
//        C.col(0).array() += x;
//        C.col(1).array() += y;
//        C.col(2).array() += z;
//        
//        gsMatrix<T> ww(16,1);
//        ww.setOnes();
//        ww.at(1)= 1./3.;
//        ww.at(2)= 1./3. ;
//        ww.at(5)= 1./3. ;
//        ww.at(6)= 1./3. ;
//        ww.at(9)= 1./3. ;
//        ww.at(10)= 1./3. ;
//        ww.at(13)= 1./3. ;
//        ww.at(14)= 1./3. ;
//        return TensorNurbs3Ptr(new gsTensorNurbs<3, T>(KVx, KVy, KVz, give(C), give(ww)));
//    }
//    gsTensorNurbs<3> geo_3d = *gsNurbsCreator<>::Nurbs3DHalfAnnulus((aDim-bDim)/2, aDim/2,0,0,0,1);
    
    mp.addAutoBoundaries(); //Make all patch sides which are not yet declared as interface or boundary to a boundary.
    std::string output("single_patch_arc");
    
    // p-refine
    if(numElevate!=0)
        mp.degreeElevate(numElevate);
    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();
    
    // Set the deformed configuration
    mp_def = mp;
    if (plot) gsWriteParaview<>( mp_def    , "mp", 1000, true);
    gsMultiBasis<> dbasis(mp);
    
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";
    
    //! [Set boundary conditions]
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    gsVector<> tmp(2);
    tmp << 0,0;
    
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    
    real_t pressure = 0.0;
    
    if (testCase == 0)
    {
        
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,-1);
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, -1);
        
        //Point load
        gsVector<> point(2);
        point<< 0.54, 1.0;
        gsVector<> load(2);
        load << 0.0, -1.0;
        pLoads.addLoad(point, load, 0); //Point, Load vector, Patch id
    }
    else
        GISMO_ERROR("Test case not known");
    //! [Set boundary conditions]
    
    //! [Make material functions]
    gsConstantFunction<> force(tmp,2);
    gsFunctionExpr<> t(std::to_string(thickness),2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
    gsFunctionExpr<> rho(std::to_string(Density),2);
    
    // Linear anisotropic material model (only one layer for example purposes)
    index_t kmax = 1; // number of layers
    
    std::vector<gsFunctionSet<> * > Gs(kmax); // Material matrices
    std::vector<gsFunctionSet<> * > Ts(kmax); // Thickness per layer
    std::vector<gsFunctionSet<> * > Phis(kmax); // Fiber angle per layer
    
    //Make material matrix
    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;
    
    // Define thickness
    gsConstantFunction<> thicks(thickness/kmax,2);
    Ts[0] = &thicks;
    
    
    // Define fiber angle
    gsConstantFunction<> phi;
    phi.setValue(0,2);
    Phis[0] = &phi;
    
    
    // Define parameters vector depending on material law
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

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    
    if      (material==0) //Linear
    {
        if (composite) // Composite
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
    else if (material==1 || material==2)
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }
    else
        GISMO_ERROR("Material "<<material<<" not supported");
    
    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2, real_t, false>(mp,dbasis,BCs,force,materialMatrix);
    
    // Add point loads to the assembler
    assembler->setPointLoads(pLoads);
    
    //! [Make assembler]

    //! [Define jacobian and residual]

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
//      GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };

    // Function for the Residual
    Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleVector(mp_def);
      GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
      return assembler->rhs();
    };

    //! [Define jacobian and residual]
    ThinShellAssemblerStatus status;
    status = assembler->assemble();

    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");


    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();

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

//
//    //! [Construct and evaluate solution]
//    mp_def = assembler->constructSolution(solVector);
//
//    gsMultiPatch<> deformation = mp_def;
//    for (size_t k = 0; k != mp_def.nPatches(); ++k)
//        deformation.patch(k).coefs() -= mp.patch(k).coefs();
//
//    gsInfo <<"Maximum deformation coef: "
//           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
//    gsInfo <<"Minimum deformation coef: "
//           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    //! [Construct and evaluate solution]
    // ! [Export visualization in ParaView]
//    if (plot)
//    {
//        gsField<> solField(mp_def, deformation);
//        gsInfo<<"Plotting in Paraview...\n";
//        gsWriteParaview<>( solField, "Deformation", 1000, true);
//    }
//    if (stress)
//    {
//        gsPiecewiseFunction<> membraneStresses;
//        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
//        gsField<> membraneStress(mp_def,membraneStresses, true);
//
//        gsPiecewiseFunction<> flexuralStresses;
//        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
//        gsField<> flexuralStress(mp_def,flexuralStresses, true);
//
//        gsPiecewiseFunction<> stretches;
//        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
//        gsField<> Stretches(mp_def,stretches, true);
//
//        gsPiecewiseFunction<> pstress_m;
//        assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
//        gsField<> pstressM(mp_def,pstress_m, true);
//
//        gsPiecewiseFunction<> pstress_f;
//        assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
//        gsField<> pstressF(mp_def,pstress_f, true);
//
//        gsPiecewiseFunction<> stretch1;
//        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
//        gsField<> stretchDir1(mp_def,stretch1, true);
//
//        gsPiecewiseFunction<> stretch2;
//        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
//        gsField<> stretchDir2(mp_def,stretch2, true);
//
//        gsPiecewiseFunction<> stretch3;
//        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
//        gsField<> stretchDir3(mp_def,stretch3, true);
//
//
//        gsWriteParaview(membraneStress,"MembraneStress");
//        gsWriteParaview(flexuralStress,"FlexuralStress");
//        gsWriteParaview(Stretches,"PrincipalStretch");
//        gsWriteParaview(pstressM,"PrincipalMembraneStress");
//        gsWriteParaview(pstressF,"PrincipalFlexuralStress");
//        gsWriteParaview(stretchDir1,"PrincipalDirection1");
//        gsWriteParaview(stretchDir1,"PrincipalDirection1");
//        gsWriteParaview(stretchDir2,"PrincipalDirection2");
//        gsWriteParaview(stretchDir3,"PrincipalDirection3");
//    }
    // ! [Export visualization in ParaView]

//    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
//    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    return EXIT_SUCCESS;

}

