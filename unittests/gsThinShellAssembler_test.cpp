/** @file gsThinShellAssembler_test.cpp

    @brief Provides unittests for the gsThinShellAssembler class

    * Balloon: unit-test based on a hyperelastic balloon inflated by a follower pressure.
               This test allows to test the follower pressure as well as the hyperelastic material models (incompressible)

    * UAT:     unit-test based on a uni-axial tension test
               This test allows to test the hyperelastic material models (incompressible and compressible)

    * Modal:   unit-test based on a modal analysis
               This test allows to test the mass matrix and the compressible material model


    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsKLShell/gsKLShell.h>

SUITE(gsThinShellAssembler_test)                 // The suite should have the same name as the file
{

    std::pair<real_t,real_t> balloon_numerical(const index_t material, const index_t impl);
    real_t balloon_analytical(const index_t material, const index_t impl, const real_t r);
    void balloon_CHECK(const index_t material, const index_t impl);

    std::pair<real_t,real_t> UAT_analytical(const index_t material, const index_t impl, const bool Compressibility);
    std::pair<real_t,real_t> UAT_numerical(const index_t material, const index_t impl, const bool Compressibility);
    void UAT_CHECK(const index_t material, const index_t impl, const bool Compressibility);

    gsVector<real_t> Modal_analytical();
    gsVector<real_t> Modal_numerical(const bool composite);
    void Modal_CHECK(const bool composite);

    /// The imposed axial stretch of the UAT fixture. UAT_numerical applies it as a
    /// Dirichlet displacement (lambda-1 on the east edge) and UAT_analytical builds
    /// its closed forms at it. The two used to carry INDEPENDENT copies of the
    /// literal 2.0 -- a second way for the oracle to go stale silently, alongside
    /// the hardcoded Jacobians D2 removed. One source now. (task 61)
    const real_t UAT_lambda = 2.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(Balloon_NH_Analytical)
    {
     index_t mat = 1, impl = 1;
     balloon_CHECK(mat, impl);
    }
    TEST(Balloon_NH_Generic)
    {
     index_t mat = 1, impl = 2;
     balloon_CHECK(mat, impl);
    }
    TEST(Balloon_NH_Spectral)
    {
     index_t mat = 1, impl = 3;
     balloon_CHECK(mat, impl);
    }

    TEST(Balloon_MR_Analytical)
    {
     index_t mat = 3, impl = 1;
     balloon_CHECK(mat, impl);
    }
    TEST(Balloon_MR_Generic)
    {
     index_t mat = 3, impl = 2;
     balloon_CHECK(mat, impl);
    }
    TEST(Balloon_MR_Spectral)
    {
     index_t mat = 3, impl = 3;
     balloon_CHECK(mat, impl);
    }

    TEST(Balloon_OG_Spectral)
    {
     index_t mat = 4, impl = 3;
     balloon_CHECK(mat, impl);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(UAT_NH_Incomp_Analytical)
    {
     index_t mat = 1, impl = 1;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_NH_Comp_Analytical)
    {
     index_t mat = 1, impl = 1;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_NH_Incomp_Generic)
    {
     index_t mat = 1, impl = 2;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_NH_Comp_Generic)
    {
     index_t mat = 1, impl = 2;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_NH_Incomp_Spectral)
    {
     index_t mat = 1, impl = 3;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_NH_Comp_Spectral)
    {
     index_t mat = 1, impl = 3;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }

    TEST(UAT_MR_Incomp_Analytical)
    {
     index_t mat = 3, impl = 1;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_MR_Comp_Analytical)
    {
     index_t mat = 3, impl = 1;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_MR_Incomp_Generic)
    {
     index_t mat = 3, impl = 2;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_MR_Comp_Generic)
    {
     index_t mat = 3, impl = 2;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_MR_Incomp_Spectral)
    {
     index_t mat = 3, impl = 3;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_MR_Comp_Spectral)
    {
     index_t mat = 3, impl = 3;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }

    TEST(UAT_OG_Incomp_Spectral)
    {
     index_t mat = 4, impl = 3;
     bool comp = false;
     UAT_CHECK(mat, impl, comp);
    }
    TEST(UAT_OG_Comp_Spectral)
    {
     index_t mat = 4, impl = 3;
     bool comp = true;
     UAT_CHECK(mat, impl, comp);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(Modal_isotropic)
    {
     bool comp = false;
     Modal_CHECK(comp);
    }
    TEST(Modal_composite)
    {
     bool comp = true;
     Modal_CHECK(comp);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::pair<real_t,real_t> balloon_numerical(const index_t material, const index_t impl)
    {
        //! [Parse command line]
        index_t numRefine  = 1;
        index_t numElevate = 1;

        real_t E_modulus = 1.0;
        real_t Density = 1.0;
        real_t Ratio = 7.0;

        real_t thickness = 0.1;
        real_t mu = 4.225e5;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        real_t PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);

        //! [Read input file]
        gsMultiPatch<> mp, mp_def;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.embed(3);

        gsReadFile<>("surfaces/eighth_sphere.xml", mp);

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numRefine; ++i)
          mp.patch(0).uniformRefine();

        mp_def = mp;

        //! [Refinement]
        gsMultiBasis<> dbasis(mp);

        gsBoundaryConditions<> bc;
        bc.setGeoMap(mp);

        GISMO_ENSURE(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Pressure
        real_t pressure = 10e3;

        //! [Refinement]

        // Linear isotropic material model
        gsVector<> tmp(3);
        tmp.setZero();
        gsConstantFunction<> force(tmp,3);
        gsConstantFunction<> pressFun(pressure,3);
        gsFunctionExpr<> t(std::to_string(thickness), 3);
        gsFunctionExpr<> E(std::to_string(E_modulus),3);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
        gsFunctionExpr<> rho(std::to_string(Density),3);
        gsConstantFunction<> ratio(Ratio,3);

        gsConstantFunction<> alpha1fun(alpha1,3);
        gsConstantFunction<> mu1fun(mu1,3);
        gsConstantFunction<> alpha2fun(alpha2,3);
        gsConstantFunction<> mu2fun(mu2,3);
        gsConstantFunction<> alpha3fun(alpha3,3);
        gsConstantFunction<> mu3fun(mu3,3);

        std::vector<gsFunctionSet<>*> parameters(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
        gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

        if (material==4)
        {
            parameters.resize(8);
            parameters[0] = &E;
            parameters[1] = &nu;
            parameters[2] = &mu1fun;
            parameters[3] = &alpha1fun;
            parameters[4] = &mu2fun;
            parameters[5] = &alpha2fun;
            parameters[6] = &mu3fun;
            parameters[7] = &alpha3fun;
        }

        gsOptionList options;
        if      (material==0)
        {
            GISMO_ERROR("This test is not available for SvK models");
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

        assembler->setPressure(pressFun);

        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
        Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          ThinShellAssemblerStatus status = assembler->assembleMatrix(mp_def);
          GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
          gsSparseMatrix<real_t> m = assembler->matrix();
          return m;
        };
        // Function for the Residual
        Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          ThinShellAssemblerStatus status = assembler->assembleVector(mp_def);
          GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
          return assembler->rhs();
        };

        // Define Matrices
        assembler->assemble();

        gsSparseMatrix<> matrix = assembler->matrix();
        gsVector<> vector = assembler->rhs();

        // Solve linear problem
        gsVector<> solVector;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute( matrix );
        solVector = solver.solve(vector);

        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        bool   converged  = false;
        index_t nIt       = 0;
        real_t updateNorm = updateVector.norm();
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);

            nIt        = it+1;
            updateNorm = updateVector.norm();
            if (updateNorm < 1e-6)
            {
                converged = true;
                break;
            }
        }
        // (task 68, F11) TWIN of the defect task 65 fixed in UAT_numerical (see the
        // note at the Newton loop of UAT_numerical in this same file). The guard
        // removed just above read
        //     else if (it+1 == it) gsWarn<<"Maximum iterations reached!\n";
        // which is ALWAYS FALSE, so 100 fruitless Newton steps exited quietly with a
        // garbage solution and the balloon tests could report an unconverged result
        // as a pass. Non-convergence is now asserted and reported.
        gsInfo << "[BALLOON_NEWTON] mat "<<material<<" impl "<<impl
               << " : converged "<<converged<<" in "<<nIt<<" its, |du| = "<<updateNorm
               << " , |R| = "<<resVec.norm()<<"\n";
        CHECK(converged);

        mp_def = assembler->constructSolution(solVector);

        gsMultiPatch<> deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        gsMatrix<> pt(2,2);
        pt.col(0)<<0.0,1.0;
        pt.col(1)<<1.0,1.0;


        gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

        real_t tol = 10e-10;
        GISMO_ENSURE((lambdas.col(0)-lambdas.col(1)).norm() < tol, "Stretches must be equal over the balloon");

        real_t r = (mp_def.patch(0).eval(pt).col(0)).norm();

        // Get the total force on the tension boundary
        real_t P = pressure * assembler->getArea(mp) / assembler->getArea(mp_def);

        std::pair<real_t,real_t> result;
        result.first = P;
        result.second = r;

        delete assembler;

        return result;
    }

    real_t balloon_analytical(const index_t material, const index_t impl, const real_t r)
    {
        GISMO_UNUSED(impl);

        real_t Pan;
        real_t R = 10.;

        real_t Ratio = 7.0;

        real_t thickness = 0.1;
        real_t mu = 4.225e5;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        real_t lambda = r/R;

        if      (material==1)
        {
            // Pan = 2*(thickness/R)*(mu*(1.0/lambdas(0)-lambdas(0)));
            Pan = 2*(thickness/R)*(mu*(math::pow(lambda,2-3)-math::pow(lambda,-2*2-3)));
        }
        else if (material==3)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            real_t m1 = c1*mu;
            real_t m2 = -c2*mu;
            real_t a1 = 2;
            real_t a2 = -2;
            Pan = 2*(thickness/R)*(m1*(math::pow(lambda,a1-3)-math::pow(lambda,-2*a1-3))+m2*(math::pow(lambda,a2-3)-math::pow(lambda,-2*a2-3)));
        }
        else if (material==4)
        {
            Pan=2*(thickness/R)*(
                mu1*(math::pow(lambda,alpha1-3)-math::pow(lambda,-2*alpha1-3))
                +mu2*(math::pow(lambda,alpha2-3)-math::pow(lambda,-2*alpha2-3))
                +mu3*(math::pow(lambda,alpha3-3)-math::pow(lambda,-2*alpha3-3)) );
        }
        else
            GISMO_ERROR("Material not treated");

        return Pan;
    }

    void balloon_CHECK(const index_t material, const index_t impl)
    {
     // Ogden (material 4) is only implemented for the Spectral implementation
     // (impl 3): getMaterialMatrix.h:235-253 raises GISMO_ERROR for every other
     // one. This is a SKIP, so it must RETURN -- without the return the vacuous
     // CHECK(true) was recorded and the body ran on anyway, straight into that
     // GISMO_ERROR as an unhandled exception (measured, task 61). No OG/non-
     // Spectral combination is registered today, so the guard is dead defensive
     // code; it becomes live the moment one is added.
     if (material==4 && impl!=3)
     {
          CHECK(true);
          return;
     }

     std::pair<real_t,real_t> num = balloon_numerical(material,impl);
     real_t P = num.first;
     real_t rnum = num.second;

     real_t Pan = balloon_analytical(material,impl,rnum);

     CHECK_CLOSE(std::abs(P-Pan)/Pan,0,1e-3);
    }

    /*  (task 65) RUNTIME identity gate.

        This whole task exists because twelve of the fourteen UAT_* tests were named
        for materials they did not run: every body passed mat = 1, impl = 1. That was
        invisible in a source diff and visible only at runtime. So the arguments are
        NOT trusted here -- gsMaterialMatrixNonlinear::print() reports the object's own
        TEMPLATE parameters <matId, comp> (gsMaterialMatrixNonlinear.hpp:153-192), i.e.
        what getMaterialMatrix ACTUALLY constructed. This asserts the constructed type,
        not the argument that was typed.
    */
    void CHECK_material_identity(const gsMaterialMatrixBase<real_t> & mm,
                                 const index_t material, const index_t impl,
                                 const bool Compressibility)
    {
        /*  THE EXPECTATION COMES FROM THE TEST NAME, NOT FROM THE ARGUMENTS.

            An earlier version of this helper compared the constructed object against
            the `material`/`impl` ARGUMENTS -- and it was MEASURED not to fire when two
            tests were reverted to the original fictional `mat = 1, impl = 1` (task 65,
            poison round 1: both suites stayed green). Of course: reverting the argument
            moves BOTH sides of that comparison. The defect this task exists to prevent
            is a mismatch between the test NAME and what the test RUNS, so the name is
            the only admissible anchor. With the name on one side and the constructed
            template instantiation on the other, the argument appears nowhere in the
            chain and the original defect becomes unrepresentable.
        */
        const std::string tn = UnitTest::CurrentTest::Details()->testName;

        std::string matName, implName, compName;
        if      (tn.find("_NH_")!=std::string::npos) matName = "Neo-Hookean\n";
        else if (tn.find("_MR_")!=std::string::npos) matName = "Mooney-Rivlin\n";
        else if (tn.find("_OG_")!=std::string::npos) matName = "Ogden\n";
        // the trailing newline matters: "Neo-Hookean" is a PREFIX of the NH_ext name

        if      (tn.find("_Analytical") !=std::string::npos) implName = "Analytical implementation";
        else if (tn.find("_Generic")    !=std::string::npos ||
                 tn.find("_Generalized")!=std::string::npos) implName = "Generalized implementation";
        else if (tn.find("_Spectral")   !=std::string::npos) implName = "Spectral implementation";

        // "_Incomp" must be tested BEFORE "_Comp"
        if      (tn.find("_Incomp")!=std::string::npos) compName = "\tIncompressible ";
        else if (tn.find("_Comp")  !=std::string::npos) compName = "\tCompressible ";

        if (matName.empty() || implName.empty() || compName.empty())
        {
            // A test that reaches the material path but whose name does not say which
            // material it runs is exactly the condition this gate exists to forbid.
            gsInfo << "[MATERIAL] FAIL: test name '"<<tn<<"' does not encode "
                      "material / implementation / compressibility\n";
            CHECK(false);
            return;
        }

        // gsMaterialMatrixNonlinear::print() reports the object's own TEMPLATE
        // parameters <matId, comp> (gsMaterialMatrixNonlinear.hpp:153-192), i.e. what
        // getMaterialMatrix ACTUALLY constructed -- not what was requested.
        std::ostringstream oss;
        mm.print(oss);
        const std::string s = oss.str();

        const bool okMat  = (s.find(matName)  != std::string::npos);
        const bool okImpl = (s.find(implName) != std::string::npos);
        const bool okComp = (s.find(compName) != std::string::npos);
        gsInfo << "[MATERIAL] "<<tn<<" : name wants "
               << compName.substr(1) << matName.substr(0,matName.size()-1) << " / " << implName
               << " ; args (mat "<<material<<", impl "<<impl<<", comp "<<Compressibility
               << ") CONSTRUCTED "
               << (okMat&&okImpl&&okComp ? "MATCH" : "MISMATCH") << "\n";
        if (!(okMat&&okImpl&&okComp))
            gsInfo << "[MATERIAL] constructed object reports:\n"<<s;
        CHECK(okMat);
        CHECK(okImpl);
        CHECK(okComp);
    }

    std::pair<real_t,real_t> UAT_numerical(const index_t material, const index_t impl, const bool Compressibility)
    {
        //! [Parse command line]
        index_t numRefine  = 1;
        index_t numElevate = 1;

        real_t E_modulus = 1.0;
        real_t PoissonRatio;
        real_t Density = 1.0;
        real_t Ratio = 7.0;

        real_t mu = 1.5e6;
        real_t thickness = 0.001;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        E_modulus = 2*mu*(1+PoissonRatio);

        //! [Parse command line]

        //! [Read input file]
        gsMultiPatch<> mp, mp_def;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree

        if (numElevate!=0)
            mp.degreeElevate(numElevate);

        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();

        mp_def = mp;

        //! [Refinement]
        gsMultiBasis<> dbasis(mp);

        gsBoundaryConditions<> bc;
        bc.setGeoMap(mp);

        gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

        const real_t lambda = UAT_lambda;   // (task 61) was an independent literal 2.0
        gsConstantFunction<> displx(lambda-1.0,2);

        GISMO_ENSURE(mp.targetDim()==2,"Geometry must be planar (targetDim=2)!");
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );

        bc.addCondition(boundary::east, condition_type::dirichlet, &displx, 0, false, 0 );

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

        //! [Refinement]

        // Linear isotropic material model
        gsVector<> tmp(2);
        tmp.setZero();
        gsConstantFunction<> force(tmp,2);
        gsFunctionExpr<> t(std::to_string(thickness),2);
        gsFunctionExpr<> E(std::to_string(E_modulus),2);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
        gsFunctionExpr<> rho(std::to_string(Density),2);
        gsConstantFunction<> ratio(Ratio,2);

        gsConstantFunction<> alpha1fun(alpha1,2);
        gsConstantFunction<> mu1fun(mu1,2);
        gsConstantFunction<> alpha2fun(alpha2,2);
        gsConstantFunction<> mu2fun(mu2,2);
        gsConstantFunction<> alpha3fun(alpha3,2);
        gsConstantFunction<> mu3fun(mu3,2);

        std::vector<gsFunctionSet<>*> parameters(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
        gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

        if (material==4)
        {
            parameters.resize(8);
            parameters[0] = &E;
            parameters[1] = &nu;
            parameters[2] = &mu1fun;
            parameters[3] = &alpha1fun;
            parameters[4] = &mu2fun;
            parameters[5] = &alpha2fun;
            parameters[6] = &mu3fun;
            parameters[7] = &alpha3fun;
        }

        gsOptionList options;
        if      (material==0)
        {
            GISMO_ERROR("This test is not available for SvK models");
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }

        // (task 65) what was CONSTRUCTED, not what was asked for.
        CHECK_material_identity(*materialMatrix,material,impl,Compressibility);

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<2, real_t, false >(mp,dbasis,bc,force,materialMatrix);

        assembler->setPointLoads(pLoads);

        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
        Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          ThinShellAssemblerStatus status = assembler->assembleMatrix(mp_def);
          GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
          gsSparseMatrix<real_t> m = assembler->matrix();
          return m;
        };
        // Function for the Residual
        Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          ThinShellAssemblerStatus status = assembler->assembleVector(mp_def);
          GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
          return assembler->rhs();
        };

        // Define Matrices
        assembler->assemble();

        gsSparseMatrix<> matrix = assembler->matrix();
        gsVector<> vector = assembler->rhs();

        // Solve linear problem
        gsVector<> solVector;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute( matrix );
        solVector = solver.solve(vector);

        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        bool   converged  = false;
        index_t nIt       = 0;
        real_t updateNorm = updateVector.norm();
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);

            nIt        = it+1;
            updateNorm = updateVector.norm();
            if (updateNorm < 1e-6)
            {
                converged = true;
                break;
            }
        }
        // (task 65) NON-CONVERGENCE USED TO BE SILENT. The guard removed just above read
        //     else if (it+1 == it) gsWarn<<"Maximum iterations reached!\n";
        // which is ALWAYS FALSE, so 100 fruitless Newton steps exited quietly with a
        // garbage solution. That matters now that materials other than Neo-Hookean are
        // actually run here: a non-converged solve would surface downstream as a large
        // L/S deviation and be misread as a material (or oracle) defect. It is a THIRD,
        // distinct failure mode, and it is now asserted and reported.
        gsInfo << "[UAT_NEWTON] mat "<<material<<" impl "<<impl
               << (Compressibility ? " compressible" : " incompressible")
               << " : converged "<<converged<<" in "<<nIt<<" its, |du| = "<<updateNorm
               << " , |R| = "<<resVec.norm()<<"\n";
        CHECK(converged);

        mp_def = assembler->constructSolution(solVector);

        gsMultiPatch<> deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Check solutions
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // NOTE: the imposed stretch is UAT_lambda, shared with UAT_analytical. The
        // compressible closed forms there no longer hardcode a Jacobian for it: J is
        // solved at the actual lambda (task 61, D2).

        // Compute stretches (should be the same everywhere)
        // Ordering: lambda(0) < lambda(1); lambda(2) is ALWAYS the through-thickness stretch
        gsVector<> pt(2);
        pt<<1,0;
        gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

        // Get the total force on the tension boundary
        patchSide ps(0,boundary::east);
        gsMatrix<> forceVector = assembler->boundaryForce(mp_def,ps);
        real_t sideForce = forceVector.sum();
        real_t S   = -sideForce / (thickness*lambdas(0)*lambdas(2));
        real_t L   = lambdas(0);

        std::pair<real_t,real_t> result;
        result.first = L;
        result.second = S;

        delete assembler;

        return result;
    }

    /*  TASK 61 / D2 -- the compressible UAT oracles used to HARDCODE the Jacobian:
            J = 1.088778638;// specific for lambda==2!!
        and the same, twice more, for the NH and MR branches. The comments were
        honest, but nothing enforced them: changing the imposed stretch left the
        analytic San/Lan built on the Jacobian of a DIFFERENT stretch. J is now
        SOLVED at the actual lambda.

        All three compressible branches below are Ogden-form,
            Psi = sum_p (mu_p/alpha_p) (b1^a_p + b2^a_p + b3^a_p - 3)
                  + K/4 (J^2 - 1 - 2 ln J),        b_i = lambda_i * J^(-1/3),
        which is exactly what their San expressions differentiate -- verified term
        by term against the code as written, and NOT assumed: NH is
        {mu_p} = {mu}, {alpha_p} = {2}; MR is {c1*mu, -c2*mu}, {2, -2}; Ogden is
        {mu1,mu2,mu3}, {alpha1,alpha2,alpha3}. In uniaxial tension along 1 with
        lambda_2 = lambda_3 = sqrt(J/lambda),
            J * sigma_22 = sum_p (mu_p/3)(b2^a_p - b1^a_p) + K/4 (2 J^2 - 2),
        and the compressible uniaxial condition is sigma_22 = sigma_33 = 0.
        UAT_J_residual is that right-hand side; UAT_solveJ bisects it.

        THE OLD CONSTANTS WERE CORRECT. Recomputed independently at 30 digits
        (mpmath) with this file's mu / mu_p / alpha_p / PoissonRatio = 0.45:
            NH  1.10559856482753778   vs the file's 1.105598565
            MR  1.09990584204437305   vs the file's 1.099905842
            OG  1.08877863796876059   vs the file's 1.088778638
        every digit they carried agrees. UAT_pinLegacyJ keeps that agreement
        ASSERTED at lambda == 2, so a broken solver cannot pass quietly either.
    */
    real_t UAT_J_residual(const real_t J, const real_t lambda, const real_t K,
                          const std::vector<real_t> & mu_p,
                          const std::vector<real_t> & alpha_p)
    {
        const real_t Jm13 = math::pow(J,-1./3.);
        const real_t b1   = lambda * Jm13;                     // lambda_bar_1
        const real_t b2   = math::pow(J/lambda,0.5) * Jm13;    // lambda_bar_2 = _3
        real_t r = 0.25*K*(2*J*J - 2);
        for (size_t p = 0; p != mu_p.size(); ++p)
            r += mu_p[p]/3. * (math::pow(b2,alpha_p[p]) - math::pow(b1,alpha_p[p]));
        return r;
    }

    real_t UAT_solveJ(const real_t lambda, const real_t K,
                      const std::vector<real_t> & mu_p,
                      const std::vector<real_t> & alpha_p)
    {
        // Bracket generously and CHECK the bracket rather than assuming it: the
        // residual -> -inf as J -> 0 (the isochoric terms blow up with the wrong
        // sign there) and ~ K/2 * J^2 -> +inf for large J.
        real_t a = 1e-3, b = 1e3;
        GISMO_ENSURE(UAT_J_residual(a,lambda,K,mu_p,alpha_p) < 0 &&
                     UAT_J_residual(b,lambda,K,mu_p,alpha_p) > 0,
                     "UAT_solveJ: sigma_22(J) is not bracketed on ["<<a<<","<<b<<"] "
                     "at lambda = "<<lambda);
        // Bisection. The width test is in units of eps so that it terminates at the
        // last representable bit for float, double and multiprecision alike; the
        // 200-iteration cap is only a guard against a non-terminating tolerance.
        for (index_t it = 0;
             it != 200 && (b-a) > 4*std::numeric_limits<real_t>::epsilon()*b; ++it)
        {
            const real_t m = 0.5*(a+b);
            if (UAT_J_residual(m,lambda,K,mu_p,alpha_p) < 0) a = m; else b = m;
        }
        return 0.5*(a+b);
    }

    /// Pin the solved @a J against the constant this file used to hardcode, while
    /// the stretch is still the one that constant was computed for. Silent at any
    /// other lambda -- which is the whole point of D2.
    void UAT_pinLegacyJ(const real_t J, const real_t lambda, const real_t legacy)
    {
        // Only at the stretch the constants were computed for -- and only when the
        // arithmetic can tell the difference. NO single absolute tolerance works for
        // every real_t here, so this is a precision GATE rather than a rescale:
        //   - the legacy constants are ROUNDED, sitting 1.7e-10 from the true root,
        //     so any tolerance must EXCEED 1.7e-10;
        //   - at real_t = float the bisection in UAT_solveJ stops at a width of
        //     4*eps*J ~ 5.2e-7, so a tolerance float could satisfy must exceed
        //     2.6e-7 -- which is VACUOUS in double, where a 3.5e-8 perturbation of
        //     the constant is exactly what was used to falsify this check.
        // Pinned in double and above; skipped below, deliberately and on the record.
        if (lambda == 2.0 && std::numeric_limits<real_t>::epsilon() < 1e-12)
            CHECK_CLOSE(legacy,J,1e-8);   // ABSOLUTE (CHECK_CLOSE always is); J = O(1)
    }

    std::pair<real_t,real_t> UAT_analytical(const index_t material, const index_t impl, const bool Compressibility)
    {
        GISMO_UNUSED(impl);

        real_t PoissonRatio;
        real_t Ratio = 7.0;

        real_t mu = 1.5e6;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        const real_t lambda = UAT_lambda;   // (task 61) was an independent literal 2.0

        real_t San,J,K,Lan;
        if      (material==1 && Compressibility)
        {
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = UAT_solveJ(lambda,K,std::vector<real_t>{mu},std::vector<real_t>{2.0});
            UAT_pinLegacyJ(J,lambda,1.105598565);   // was hardcoded, "specific for lambda==2!!"
            San = lambda*(0.5*mu*(-(2*(math::pow(lambda,2)+2*J/lambda))/(3*math::pow(J,2./3.)*lambda)+2*lambda/math::pow(J,2./3.))+0.25*K*(2*math::pow(J,2)/lambda-2./lambda))/J;
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==1 && !Compressibility)
        {
            San = mu * (lambda*lambda - 1/lambda);
            Lan = math::pow(1./lambda,0.5);
        }
        else if (material==3 && Compressibility)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = UAT_solveJ(lambda,K,std::vector<real_t>{c1*mu,-c2*mu},std::vector<real_t>{2.0,-2.0});
            UAT_pinLegacyJ(J,lambda,1.099905842);   // was hardcoded, "specific for lambda==2!!"
            San = lambda*(0.5*c1*mu*(-(2*(math::pow(lambda,2)+2*J/lambda))/(3*math::pow(J,2./3.)*lambda)+2*lambda/math::pow(J,2./3.))+0.5*c2*mu*(-(4*(2*lambda*J+math::pow(J,2)/math::pow(lambda,2)))/(3*math::pow(J,4./3.)*lambda)+4/math::pow(J,1./3.))+0.25*K*(2*math::pow(J,2)/lambda-2/lambda))/J;
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==3 && !Compressibility)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            San =-mu*(c2*lambda*lambda+c2/lambda+c1)/lambda+lambda*(c1*lambda*mu+2*c2*mu);
            Lan = math::pow(1./lambda,0.5);
        }
        else if (material==4 && Compressibility)
        {
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = UAT_solveJ(lambda,K,std::vector<real_t>{mu1,mu2,mu3},
                                    std::vector<real_t>{alpha1,alpha2,alpha3});
            UAT_pinLegacyJ(J,lambda,1.088778638);   // was hardcoded, "specific for lambda==2!!"
            San = 1./J* (lambda *( mu1*(2*math::pow(lambda/math::pow(J,1./3.),alpha1)*alpha1/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha1)*alpha1/(3*lambda))/alpha1+mu2*(2*math::pow(lambda/math::pow(J,1./3.),alpha2)*alpha2/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha2)*alpha2/(3*lambda))/alpha2+mu3*(2*math::pow(lambda/math::pow(J,1./3.),alpha3)*alpha3/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha3)*alpha3/(3*lambda))/alpha3+0.25*K*(2*math::pow(J,2)/lambda-2/lambda) ) );
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==4 && !Compressibility)
        {
            San =-mu1*math::pow((1./lambda),0.5*alpha1)-mu2*math::pow((1./lambda),0.5*alpha2)-mu3*math::pow((1./lambda),0.5*alpha3)+mu1*math::pow(lambda,alpha1)+mu2*math::pow(lambda,alpha2)+mu3*math::pow(lambda,alpha3);
            Lan = math::pow(1./lambda,0.5);
        }
        else
            GISMO_ERROR("Material not treated");

        std::pair<real_t,real_t> result;
        result.first = Lan;
        result.second = San;
        return result;
    }

    void UAT_CHECK(const index_t material, const index_t impl, const bool Compressibility)
    {
        // See balloon_CHECK above: the skip must RETURN.
        if (material==4 && impl!=3)
        {
          CHECK(true);
          return;
        }

        real_t Lnum, Snum, Lana, Sana;
        std::tie(Lnum,Snum) = UAT_numerical(material,impl,Compressibility);
        std::tie(Lana,Sana) = UAT_analytical(material,impl,Compressibility);
        CHECK_CLOSE(std::abs(Lnum-Lana)/Lana,0,1e-9);
        /*  (task 65) THE STRESS IS NOW ASSERTED TOO -- task 61's named gap (HN-2)
            closed, and it is not optional here.

            For every INCOMPRESSIBLE branch the analytic lateral stretch is
                Lan = pow(1/lambda,0.5)
            with NO material dependence at all (see UAT_analytical: the material==1,
            3 and 4 incompressible branches all set exactly that). So the L check
            alone CANNOT distinguish Neo-Hookean from Mooney-Rivlin from Ogden -- it
            gates incompressibility, not the material. Seven of these fourteen tests
            would therefore have stayed vacuous with respect to their own name even
            after the arguments were corrected. San is the only material-carrying
            quantity in the incompressible half.

            TOLERANCE, chosen on measured evidence and NOT fitted to the paths this
            task newly exercises: on the two honest pre-existing tests (NH /
            Analytical, the ONLY UAT combination that was ever really run) the
            measured relative stress residual is 1.8e-16 incompressible and 1.6e-12
            compressible. 1e-9 is the same gate the stretch already carries, and
            leaves ~600x headroom on the worst honest path. It was fixed BEFORE any
            argument was re-armed. The deformation here is homogeneous and exactly
            representable in the spline space, which is why these residuals are at
            solver level rather than at discretisation level.
        */
        CHECK_CLOSE(std::abs(Snum-Sana)/std::abs(Sana),0,1e-9);
        gsInfo << "[UAT_CHECK] material "<<material<<" impl "<<impl
               << (Compressibility ? " compressible" : " incompressible")
               << " : L "<<Lnum<<" vs "<<Lana<<" (rel "<<std::abs(Lnum-Lana)/Lana
               << ") ; S "<<Snum<<" vs "<<Sana<<" (rel "<<std::abs(Snum-Sana)/std::abs(Sana)
               << ")\n";
    }

    gsVector<real_t> Modal_numerical(bool composite)
    {
        // Input options
        int numElevate  = 2;
        int numHref     = 4;

        real_t thickness = 0.01;
        real_t E_modulus = 1e5;
        real_t Density = 1e0;
        real_t PoissonRatio = 0.3;

        gsMultiPatch<> mp;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);

        for(index_t i = 0; i< numElevate; ++i)
            mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
            mp.patch(0).uniformRefine();

        gsMultiBasis<> dbasis(mp);

        // Boundary conditions
        gsBoundaryConditions<> BCs;

        // Plate
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.setGeoMap(mp);

        // Initialise solution object
        gsMultiPatch<> mp_def = mp;

        // Linear isotropic material model
        gsVector<> tmp(3);
        tmp << 0, 0, 0;
        gsConstantFunction<> force(tmp,3);
        gsFunctionExpr<> t(std::to_string(thickness), 3);
        gsFunctionExpr<> E(std::to_string(E_modulus),3);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
        gsConstantFunction<> rho(Density,3);

        // Linear anisotropic material model
        real_t pi = math::atan(1)*4;
        index_t kmax = 5;

        std::vector<gsFunctionSet<> * > Gs(kmax);
        std::vector<gsFunctionSet<> * > Ts(kmax);
        std::vector<gsFunctionSet<> * > Rs(kmax);
        std::vector<gsFunctionSet<> * > Phis(kmax);

        Rs[0] = Rs[1] = Rs[2] = Rs[3] = Rs[4] = &rho;


        gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
        Gmat.resize(Gmat.rows()*Gmat.cols(),1);
        gsConstantFunction<> Gfun(Gmat,3);
        Gs[0] = Gs[1] = Gs[2] = Gs[3] = Gs[4] = &Gfun;

        gsConstantFunction<> phi1, phi2, phi3, phi4, phi5;
        phi1.setValue(0/kmax * pi / 2.0,3);
        phi2.setValue(1/kmax * pi / 2.0,3);
        phi3.setValue(2/kmax * pi / 2.0,3);
        phi4.setValue(3/kmax * pi / 2.0,3);
        phi5.setValue(4/kmax * pi / 2.0,3);

        Phis[0] = &phi1;
        Phis[1] = &phi2;
        Phis[2] = &phi3;
        Phis[3] = &phi4;
        Phis[4] = &phi5;

        gsConstantFunction<> thicks(thickness/kmax,3);
        Ts[0] = Ts[1] = Ts[2] = Ts[3] = Ts[4] = &thicks;

        std::vector<gsFunctionSet<>*> parameters;
        gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

        gsOptionList options;

        if (composite)
        {
            materialMatrix = memory::make_unique(new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis,Rs));
        }
        else
        {
            parameters.resize(2);
            parameters[0] = &E;
            parameters[1] = &nu;
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }

        // options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        // options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",static_cast<int>(!composite));
        // materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

        assembler->assemble();
        gsSparseMatrix<> K =  assembler->matrix();
        assembler->assembleMass();
        gsSparseMatrix<> M =  assembler->matrix();

        gsEigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;
        eigSolver.compute(K,M);
        gsMatrix<> values  = eigSolver.eigenvalues();
        gsMatrix<> vectors = eigSolver.eigenvectors();

        values = values.cwiseSqrt();
        values = values.col(0).head(10);

        delete assembler;

        return values;
    }

    gsVector<real_t> Modal_analytical()
    {
        real_t thickness = 0.01;
        real_t E_modulus = 1e5;
        real_t Density = 1e0;
        real_t PoissonRatio = 0.3;

        real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));

        std::vector<real_t> omegas;
        for (index_t m=1; m!=10; m++)
          for (index_t n=1; n!=10; n++)
            omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        std::sort(omegas.begin(),omegas.end());
        omegas.resize(10);
        gsAsVector<> analytical(omegas);

        return analytical;
    }

    void Modal_CHECK(const bool composite)
    {
     gsVector<real_t> num = Modal_numerical(composite);

     gsVector<real_t> ana = Modal_analytical();

     gsVector<real_t> relError = (num - ana).array()/ana.array();
     gsVector<real_t> zeros = gsVector<real_t>::Zero(relError.rows());
     CHECK_MATRIX_CLOSE(relError,zeros,1e-3);
    }
}
