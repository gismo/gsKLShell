/** @file example_pointload.cpp

    @brief Asymmetric loading applied with an offset angle (pi/50)
    from the crown of a semicircular arch.

    Example 7 from Yang et al 2006

    Yang, Y.B., Lin, S.P., & Leu,L.J. (2006) Solution strategy and rigid element for nonlinear analysis of elastically structures based on updated Lagrangian formulation. Engineering Structures, 29, 1189-1200. https://doi.org/10.1016/j.engstruct.2006.08.015

    Author(s): J. Li

    TODO: Not able to find the density information
 **/
#include <gismo.h>

#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsThinShellAssembler.h>


#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
using namespace gismo;


template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

template <class T>
gsMultiPatch<T> Strip(T Lb, T Hw, T x = 0, T y = 0, T z = 0);


int main(int argc, char* argv[]) {
    // Input options
    bool plot = true; // If set to true, paraview file is generated and launched on exit.
    bool write = true;
    bool xml = false;
    bool stress = false;
    bool membrane = false;
    bool nonlinear = false;
    index_t numRefine = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
//    bool Compressibility = false;
    bool verbose = false;
    int step = 10;
//    int method = 2；// (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)


    // Material properties
    index_t material = 0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 200;
    real_t PoissonRatio = 0.4;
    real_t Density = 1.0; //Couldn't find

    real_t ifcDirichlet = 1.0;
    real_t ifcClamped = 1.0;

    int result = 0;
    index_t maxit = 50;


    // Arc length method options
    real_t dL = 2e-1; // General arc length
    real_t tol = 1e-6;
    real_t tolU = 1e-6;
    real_t tolF = 1e-3;
    real_t relax = 1.0;
    bool adaptive = false;
    int method = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)



    // Read solver option file
    std::string assemberOptionsFile("options/solver_options.xml");
    std::string wn("data.csv");


    // Arc length method options
    gsCmdLine cmd(
            "Arc-length analysis of a semicircular arch with pinned-pinned boundary subjected to an asymmetrical force.");
    cmd.addInt("e", "degreeElevation",
               "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)",
               numElevate);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("m", "Material", "Material law", material);
//    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt("I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral", impl);
    cmd.addSwitch("composite", "Composite material", composite);
    cmd.addReal("A", "relaxation", "Relaxation factor for arc length method", relax);
    cmd.addInt("q", "QuasiNewtonInt", "Use the Quasi Newton method every INT iterations", quasiNewtonInt);


//
//    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
//    cmd.addSwitch("weak", "Impose boundary conditions weakly", weak);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("xml", "Write geometry into XML files", xml);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }


    //! [Read input file]
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;


    //! [Set material parameters]
    real_t I = 1;
    real_t A = 1;
    real_t PI = 3.1415926535;
    real_t aDim = 50;
    real_t bDim = pow(8 * I / PI, 1 / 4.0);
    real_t thickness = A / bDim;
    gsInfo << thickness << bDim;


    gsReadFile<>("surface/half_cylinder.xml", mp);
    mp.embed(3);
    mp.patch(0).coefs().row(0) << 0.00, 0.00, -aDim;
    mp.patch(0).coefs().row(1) << bDim, 0.00, -aDim;
    mp.patch(0).coefs().row(2) << 2 * bDim, 0.00, -aDim;
    mp.patch(0).coefs().row(3) << 0.00, aDim, -aDim;
    mp.patch(0).coefs().row(4) << bDim, aDim, -aDim;
    mp.patch(0).coefs().row(5) << 2 * bDim, aDim, -aDim;
    mp.patch(0).coefs().row(6) << 0.00, aDim, aDim;
    mp.patch(0).coefs().row(7) << bDim, aDim, aDim;
    mp.patch(0).coefs().row(8) << 2 * bDim, aDim, aDim;
    mp.patch(0).coefs().row(9) << 0.00, 0.00, aDim;
    mp.patch(0).coefs().row(10) << bDim, 0.00, aDim;
    mp.patch(0).coefs().row(11) << 2 * bDim, 0.00, aDim;
    gsWrite<>(mp, "Strip_example");

    //![Make geometry and refine/elevate]
    //p-refine
    if (numElevate != 0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r = 0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    mp_def.patch(0).coefs() -= mp.patch(0).coefs();
    gsWriteParaview(mp_def, "mp", 1000, true);

    //! [Refine and elevate]
    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: " << mp.nPatches() << ", degree: " << dbasis.minCwiseDegree() << "\n";
    gsInfo << dbasis.basis(0) << "\n";

    //![Set boundary conditions]
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    gsPiecewiseFunction<> force(mp.nPatches());
    gsPiecewiseFunction<> t(mp.nPatches());

    gsVector<> tmp(3);
    tmp << 0, 0, 0;


    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();


    real_t load = -2.0;
    gsVector<> neu(3);
    neu << 0.0, 0.0, load;
    gsFunctionExpr<> neuDataFun1;
    gsConstantFunction<> neuData(neu, 3);
    real_t pressure = 0.0;

    gsVector<> refPoint(2);
    real_t refPatch = 0;

    // Pinched-pinched
//    BCs.addCondition(boundary::north, condition_type::neumann, &neuData);
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1);
    BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2);

    for (index_t d = 0; d != 3; d++) {
        BCs.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);
        BCs.addCondition(0, boundary::west, condition_type::dirichlet, 0, 0, false, d);

    }

    gsMatrix<> writePoints(2, 1);
    writePoints.col(0) << 0.5, 0.54;

    std::string dirname = "ArcLengthResults";
    dirname = dirname + "/" + "Point_load";
    std::string output = "solution";
    wn = output + "data.txt";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    int systemRet = system(command);
    GISMO_ASSERT(systemRet != -1, "Something went wrong with calling the system argument");

    //thickness
    gsFunctionExpr<> t0(std::to_string(thickness), 3);
    t.addPiece(t0);
    // Point loads
    gsVector<> point(2);
    point << 0.5, 0.54;
    refPoint = point;
        gsVector<> newload(3);
        newload << 0.0,0.0,-2.0;
    pLoads.addLoad(point, newload, 0);
    // Surface forces
    tmp << 0, 0, 0;
    gsConstantFunction<> force0(tmp, 3);
    force.addPiece(force0);

    // plot geometry
    if (plot)
        gsWriteParaview(mp, dirname + "/" + "mp", 1000, true);
    if (write)
        initStepOutput(dirname + "/" + wn, writePoints);


//
    //! [Make material functions]
    // Linear isotropic material model
    gsConstantFunction<> pressFun(pressure, 3);
//    gsConstantFunction<> t(thickness,3);
    gsFunctionExpr<> E(std::to_string(E_modulus), 3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio), 3);
    gsFunctionExpr<> rho(std::to_string(Density), 3);

    // Linear anisotropic material model (only one layer for example purposes)
    index_t kmax = 1; // number of layers
    std::vector<gsFunctionSet<> *> Gs(kmax); // Material matrices
    std::vector<gsFunctionSet<> *> Ts(kmax); // Thickness per layer
    std::vector<gsFunctionSet<> *> Phis(kmax); // Fiber angle per layer

    //Make material matrix
    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus, E_modulus, 0.5 * E_modulus / (1 + PoissonRatio), PoissonRatio,
                                        PoissonRatio);
    Gmat.resize(Gmat.rows() * Gmat.cols(), 1);
    gsConstantFunction<> Gfun(Gmat, 3);
    Gs[0] = &Gfun;

    // Define thickness
    gsConstantFunction<> thicks(thickness / kmax, 3);
    Ts[0] = &thicks;

    // Define fiber angle
    gsConstantFunction<> phi;
    phi.setValue(0, 3);
    Phis[0] = &phi;

    //! [Make assembler]
    std::vector<gsFunctionSet<> *> parameters;
    gsMaterialMatrixBase<real_t> *materialMatrix;
    gsOptionList options;

    if (composite) {
        materialMatrix = new gsMaterialMatrixComposite<3, real_t>(mp, Ts, Gs, Phis);
    } else {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
        options.addInt("Material", "Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden", 0);
        options.addInt("Implementation",
                       "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
        materialMatrix = getMaterialMatrix<3, real_t>(mp, t, parameters, rho, options);
    }

    gsMaterialMatrixContainer<real_t> materialMats(mp.nPatches());
    for (size_t p = 0; p != mp.nPatches(); p++)
        materialMats.add(materialMatrix);

    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t> *assembler;

    if (membrane) // no bending term
        assembler = new gsThinShellAssembler<3, real_t, false>(mp, dbasis, BCs, force, materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true>(mp, dbasis, BCs, force, materialMatrix);

    // Set the penalty parameter for the interface C1 continuty
    assembler->options().setReal("IfcDirichlet", ifcDirichlet);
    assembler->options().setReal("IfcClamped", ifcClamped);
    assembler->addWeakC0(mp.topology().interfaces());
    assembler->addWeakC1(mp.topology().interfaces());
    assembler->initInterfaces();

    assembler->setPointLoads(pLoads);

    if (pressure != 0.0)
        assembler->setPressure(pressFun);
    //! [Make assembler]

    //! [Define jacobian and residual]
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t>(gsVector<real_t> const &)> Jacobian_t;
    typedef std::function<gsVector<real_t>(gsVector<real_t> const &)> Residual_t;
    typedef std::function<gsVector<real_t>(gsVector<real_t> const &, real_t, gsVector<real_t> const &)> ALResidual_t;

    Jacobian_t Jacobian = [&assembler, &mp_def](gsVector<real_t> const &x) {
        ThinShellAssemblerStatus status;
        assembler->constructSolution(x, mp_def);
        status = assembler->assembleMatrix(mp_def);
        GISMO_ENSURE(status == ThinShellAssemblerStatus::Success, "Assembly failed");
        gsSparseMatrix<real_t> m = assembler->matrix();
        return m;
    };

    // Function for the Residual
    Residual_t Residual = [&assembler, &mp_def](gsVector<real_t> const &x) {
        ThinShellAssemblerStatus status;
        assembler->constructSolution(x, mp_def);
        status = assembler->assembleVector(mp_def);
        GISMO_ENSURE(status == ThinShellAssemblerStatus::Success, "Assembly failed");
        return assembler->rhs();
    };

    // Function for the Residual
    ALResidual_t ALResidual = [&assembler, &mp_def](gsVector<real_t> const &x, real_t lam,
                                                    gsVector<real_t> const &force) {
        assembler->constructSolution(x, mp_def);
        assembler->assembleVector(mp_def);
        gsVector<real_t> Fint = -(assembler->rhs() - force);
        gsVector<real_t> result = Fint - lam * force;
        return result; // - lam * force;
    };

    //! [Define jacobian and residual]
    ThinShellAssemblerStatus status;
    status = assembler->assemble();
    GISMO_ENSURE(status == ThinShellAssemblerStatus::Success, "Assembly failed");

    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> Force = assembler->rhs();

    //! [Assemble linear part]

    //! [Assemble ALM]
    gsALMBase<real_t> *arcLength;
    if (method == 0)
        arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method == 1)
        arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method == 2)
        arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
        GISMO_ERROR("Method " << method << " unknown");

    arcLength->options().setString("Solver", "SimplicialLDLT"); // LDLT solver
    arcLength->options().setInt("BifurcationMethod", 0); // 0: determinant, 1: eigenvalue
    arcLength->options().setReal("Length", dL);
    arcLength->options().setInt("AngleMethod", 0); // 0: step, 1: iteration
    arcLength->options().setSwitch("AdaptiveLength", adaptive);
    arcLength->options().setInt("AdaptiveIterations", 5);
    arcLength->options().setReal("Scaling", 0.0);
    arcLength->options().setReal("Tol", tol);
    arcLength->options().setReal("TolU", tolU);
    arcLength->options().setReal("TolF", tolF);
    arcLength->options().setInt("MaxIter", maxit);
    arcLength->options().setSwitch("Verbose", true);
    arcLength->options().setReal("Relaxation", relax);

    if (quasiNewtonInt > 0) {
        quasiNewton = true;
        arcLength->options().setInt("QuasiIterations", quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi", quasiNewton);


    gsDebug << arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();

    gsParaviewCollection collection(dirname + "/" + output);
//    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
//    gsParaviewCollection Sflexural(dirname + "/" + "flexural");
//    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR


    for (index_t k = 0; k < step; k++) {
        gsInfo << "Load step " << k << "\n";
        arcLength->step();

        if (!(arcLength->converged()))
            GISMO_ERROR("Loop terminated, arc length method did not converge.\n");

        solVector = arcLength->solutionU();
        Uold = solVector;

        assembler->constructSolution(solVector, mp_def);

        gsMatrix<> pts(2, 1);
        pts << 0.0, 1.0;

        deformation = mp_def;
        deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here


        if (plot) {
            gsField<> solField;
            solField = gsField<>(mp, deformation);

            std::string fileName = dirname + "/" + output + util::to_string(k);
            gsWriteParaview<>(solField, fileName, 1000, true);
            fileName = output + util::to_string(k) + "0";
            collection.addTimestep(fileName, k, ".vts");
            collection.addTimestep(fileName, k, "_mesh.vtp");
        }
        if (stress) {
            std::string fileName;

            gsField<> membraneStress, flexuralStress, membraneStress_p;

            gsPiecewiseFunction<> membraneStresses;
            assembler->constructStress(mp_def, membraneStresses, stress_type::membrane);
            membraneStress = gsField<>(mp, membraneStresses, true);

            fileName = dirname + "/" + "membrane" + util::to_string(k);
            gsWriteParaview(membraneStress, fileName, 1000);
            fileName = "membrane" + util::to_string(k) + "0";
//            Smembrane.addTimestep(fileName, k, ".vts");

            gsPiecewiseFunction<> flexuralStresses;
            assembler->constructStress(mp_def, flexuralStresses, stress_type::flexural);
            flexuralStress = gsField<>(mp, flexuralStresses, true);

            fileName = dirname + "/" + "flexural" + util::to_string(k);
            gsWriteParaview(flexuralStress, fileName, 1000);
            fileName = "flexural" + util::to_string(k) + "0";
//            Sflexural.addTimestep(fileName, k, ".vts");

            if (impl == 3) {
                gsPiecewiseFunction<> membraneStresses_p;
                assembler->constructStress(mp_def, membraneStresses_p, stress_type::principal_stress_membrane);
                membraneStress_p = gsField<>(mp, membraneStresses_p, true);

                fileName = dirname + "/" + "membrane_p" + util::to_string(k);
                gsWriteParaview(membraneStress_p, fileName, 1000);
                fileName = "membrane_p" + util::to_string(k) + "0";
//                Smembrane_p.addTimestep(fileName, k, ".vts");
            }

        }
    }

        //! [Assemble ALM]


        //! [Construct and evaluate solution]
        gsVector<> refVals = deformation.patch(refPatch).eval(refPoint);

        gsInfo << "Displacement at reference point: " << refVals << "\n";
        //! [Construct and evaluate solution]

        //! [Export visualization in ParaView]
        if (plot) {
            collection.save();
//            Smembrane.save();
//            Sflexural.save();
//            Smembrane_p.save();
        }

    if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints,1, 201);
        //! [Export visualization in ParaView]
        delete assembler;
        delete materialMatrix;
        delete arcLength;
        return EXIT_SUCCESS;
    }


template<class T>
gsMultiPatch<T> Strip(T Lb, T Hw, T x, T y, T z) {
    gsMultiPatch<T> result, tmp;

    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0) << 0, 0, 0;
    result.patch(0).coefs().row(1) << 0, Lb, 0;
    result.patch(0).coefs().row(2) << 0, 0, Hw;
    result.patch(0).coefs().row(3) << 0, Lb, Hw;

    for (size_t p = 0; p != result.nPatches(); p++) {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }
    result.addAutoBoundaries();
    return result;
}


template<class T>
void initStepOutput(const std::string name, const gsMatrix<T> &points) {
    std::ofstream file;
    file.open(name, std::ofstream::out);
    file << std::setprecision(20)
             << "Deformation norm" << ",";
    for (index_t k = 0; k != points.cols(); k++) {
        file << "point " << k << " - x" << ","
        << "point " << k << " - y" << ","
        << "point " << k << " - z" << ",";
    }

    file << "Lambda" << ","
    << "Indicator"
    << "\n";
    file.close();

    gsInfo << "Step results will be written in file: " << name << "\n";
}

template<class T>
void writeStepOutput(const gsALMBase<T> *arcLength, const gsMultiPatch<T> &deformation, const std::string name,
                         const gsMatrix<T> &points, const index_t extreme,
                         const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
    {
        gsMatrix<T> P(2, 1), Q(2, 1);
        gsMatrix<T> out(3, points.cols());
        gsMatrix<T> tmp;

        for (index_t p = 0; p != points.cols(); p++) {
            P << points.col(p);
            deformation.patch(0).eval_into(P, tmp);
            out.col(p) = tmp;
        }

        std::ofstream file;
        file.open(name, std::ofstream::out | std::ofstream::app);
        if (extreme == -1) {
            file << std::setprecision(6)
                 << arcLength->solutionU().norm() << ",";
            for (index_t p = 0; p != points.cols(); p++) {
                file << out(0, p) << ","
                     << out(1, p) << ","
                     << out(2, p) << ",";
            }

            file << arcLength->solutionL() << ","
                 << arcLength->indicator() << ","
                 << "\n";
        } else if (extreme == 0 || extreme == 1) {
            gsMatrix<T> out2(kmax, points.cols()); // evaluation points in the rows, output (per coordinate) in columns
            for (int p = 0; p != points.cols(); p++) {
                Q.at(1 - extreme) = points(1 - extreme, p);
                for (int k = 0; k != kmax; k++) {
                    Q.at(extreme) = 1.0 * k / (kmax - 1);
                    deformation.patch(0).eval_into(Q, tmp);
                    out2(k, p) = tmp.at(2); // z coordinate
                }
            }

            file << std::setprecision(6)
                 << arcLength->solutionU().norm() << ",";
            for (index_t p = 0; p != points.cols(); p++) {
                file << out(0, p) << ","
                     << out(1, p) << ","
                     << std::max(abs(out2.col(p).maxCoeff()), abs(out2.col(p).minCoeff())) << ",";
            }

            file << arcLength->solutionL() << ","
                 << arcLength->indicator() << ","
                 << "\n";
        } else
            GISMO_ERROR("Extremes setting unknown");

        file.close();
    }
