/** @file example_shell3D.cpp

    @brief Examples for the surface thin shells including the shell obstacle course

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

//#include <gsThinShell/gsNewtonIterator.h>

template<typename T>
index_t sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

using namespace gismo;

int main(int argc, char *argv[])
{
    // Number of adaptive refinement loops
    index_t RefineLoopMax;
    // Flag for refinemet criterion
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    index_t refCriterion;
    // Parameter for computing adaptive refinement threshold
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    real_t refParameter; // ...specified below with the examples

    //! [Parse command line]
    bool plot = false;
    index_t numRefine = 1;
    index_t numElevate = 1;
    index_t goal = 1;
    bool nonlinear = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.3;
    real_t Density = 1.0;
    real_t thickness = 0.01;

    index_t modeIdx = 0;

    refCriterion = GARU;
    refParameter = 0.85;
    RefineLoopMax = 1;

    int testCase = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt("i", "index", "index of mode", modeIdx);
    cmd.addInt("e", "degreeElevation",
               "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("R", "refine", "Maximum number of adaptive refinement steps to perform",
               RefineLoopMax);
    cmd.addInt("t", "testcase",
                "Test case: 0: Beam - pinned-pinned, 1: Beam - fixed-fixed, 2: beam - fixed-free, 3: plate - fully pinned, 4: plate - fully fixed, 5: circle - fully pinned, 6: 5: circle - fully fixed",
               testCase);
    cmd.addReal("T", "thickness", "thickness", thickness);
    cmd.addInt("g", "goal", "Goal function to use", goal);
    cmd.addString("f", "file", "Input XML file", fn);
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsVector<> pts(2);
    pts.setConstant(0.25);

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    if (testCase==0)
    {
        // Unit square
        mp.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1.0;
        Density = 1.0;
        PoissonRatio = 0.3;
        thickness = 0.01;
    }
    else if (testCase==1)
    {
        std::string fn = "planar/unitcircle.xml";
        gsReadFile<>(fn, mp);
        thickness = 0.01;
        PoissonRatio = 0.3;
        E_modulus     = 1e0;
        Density = 1e0;
    }


    // thickness = 1.0;

    // p-refine
    if (numElevate != 0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r = 0; r < numRefine; ++r)
        mp.uniformRefine();

    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: " << mp.nPatches() << ", degree: " << dbasis.minCwiseDegree() << "\n";
    gsInfo << dbasis.basis(0) << "\n";

    // // Cast all patches of the mp object to THB splines
    // gsTHBSpline<2,real_t> thb;
    // for (index_t k=0; k!=mp.nPatches(); ++k)
    // {
    //     gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
    //     thb = gsTHBSpline<2,real_t>(*geo);
    //     mp.patch(k) = thb;
    // }

    mp_def = mp;

    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH = basisL;
    basisH.degreeElevate(1);

    gsInfo << "Basis Primal: " << basisL.basis(0) << "\n";
    gsInfo << "Basis Dual:   " << basisH.basis(0) << "\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    real_t D = E_modulus * math::pow(thickness, 3) / (12 * (1 - math::pow(PoissonRatio, 2)));
    gsInfo<<"D = "<<D<<"\n";

    std::vector<real_t> omegas;
    if (testCase==0)
    {
        for (index_t i = 0; i != 3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i);
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i);
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i);
        }
        for (index_t m=1; m!=10; m++)
          for (index_t n=1; n!=10; n++)
            omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        std::sort(omegas.begin(),omegas.end());
    }
    else if (testCase==1)
    {
        // Circle
        // Pinned-Pinned-Pinned-Pinned
        for (index_t i = 0; i != 3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i);
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i);
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i);
        }
        bc.addCondition(boundary::north, condition_type::clamped,0,0,false,2);
        bc.addCondition(boundary::east, condition_type::clamped,0,0,false,2);
        bc.addCondition(boundary::south, condition_type::clamped,0,0,false,2);
        bc.addCondition(boundary::west, condition_type::clamped,0,0,false,2);
        gsVector<> gammas(10);
        gammas<<3.1962206158252, 4.6108999, 5.9056782, 6.3064370, 7.1435310, 7.7992738, 8.3466059, 9.1968826, 9.4394991, 9.5257014;
        for (index_t n=0; n!=gammas.size(); n++)
          omegas.push_back(math::pow(math::pow(gammas[n],4)*D/(Density*thickness),0.5));
    }
    else
        GISMO_ERROR("TESTCASE UNKNOWN!");

    //! [Refinement]

    gsConstantFunction<> force(tmp, 3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus), 3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio), 3);
    gsFunctionExpr<> rho(std::to_string(Density), 3);

    std::vector<gsFunction<> *> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrixBase<real_t> *materialMatrix;
    gsOptionList options;
    options.addInt("Material", "Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden", 0);
    options.addInt("Implementation", "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
    materialMatrix = getMaterialMatrix<3, real_t>(mp, t, parameters, rho, options);

    gsThinShellAssemblerDWR<3, real_t, true, GoalFunction::Modal> DWR(mp, basisL, basisH, bc, force, materialMatrix);

    gsSparseSolver<>::LU solver;
    gsVector<> solVector, solVectorDualL, solVectorDualH;
    gsMultiPatch<> primalL, dualL, dualH;

    gsMatrix<> points(2, 0);
    // points.col(0).setConstant(0.25);
    // points.col(1).setConstant(0.50);
    // points.col(2).setConstant(0.75);

    gsInfo << "Assembling primal... " << std::flush;
    DWR.assembleMatrixL();
    DWR.assembleMassL();
    gsInfo << "done\n";

    // Solve system
    gsInfo << "Solving primal, size =" << DWR.matrixL().rows() << "," << DWR.matrixL().cols() << "... " << std::flush;
    Eigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> eigSolver;
    eigSolver.compute(DWR.matrixL(), DWR.massL());
    gsDebugVar(math::sqrt(eigSolver.eigenvalues()[modeIdx]));

    solVector = solVectorDualL = eigSolver.eigenvectors().col(modeIdx);

    real_t eigvalL, dualvalL;
    eigvalL = dualvalL = eigSolver.eigenvalues()[modeIdx];

    real_t Mnorm;
    // Mass-normalize primal
    solVector = 1 / (solVector.transpose() * DWR.massL() * solVector) * solVector;
    DWR.constructMultiPatchL(solVector, primalL);
    DWR.constructSolutionL(solVector, mp_def);
    gsField<> primalLField(mp, primalL);
    gsWriteParaview(primalLField, "primalL", 1000);


    // mass-normalize w.r.t. primal
    DWR.constructMultiPatchL(solVectorDualL, dualL);
    Mnorm = DWR.matrixNorm(primalL, dualL);
    solVectorDualL *= 1. / Mnorm;
    DWR.constructMultiPatchL(solVectorDualL, dualL);
    gsField<> dualLField(mp, dualL);
    gsWriteParaview(dualLField, "dualL", 1000);

    // solVectorDualL = 1 / (solVectorDualL.transpose() * DWR.massL() * solVectorDualL) * solVectorDualL;
    // DWR.constructMultiPatchL(solVectorDualL, dualL);

    gsInfo << "done.\n";

    gsInfo << "Assembling dual matrix (H)... " << std::flush;
    DWR.assembleMatrixH();
    DWR.assembleMassH();
    gsInfo << "done.\n";

    gsInfo << "Solving dual (high), size = " << DWR.matrixH().rows() << "," << DWR.matrixH().cols() << "... " << std::flush;
    eigSolver.compute(DWR.matrixH(), DWR.massH());
    gsDebugVar(math::sqrt(eigSolver.eigenvalues()[modeIdx]));
    solVectorDualH = eigSolver.eigenvectors().col(modeIdx);
    real_t dualvalH = eigSolver.eigenvalues()[modeIdx];

    // mass-normalize w.r.t. primal
    DWR.constructMultiPatchH(solVectorDualH, dualH);
    gsField<> dualHField(mp, dualH);
    gsWriteParaview(dualHField, "dualH", 1000);
    Mnorm = DWR.matrixNorm(primalL, dualH);
    gsDebugVar(Mnorm);
    solVectorDualH *= 1. / Mnorm;
    DWR.constructMultiPatchH(solVectorDualH, dualH);

    // solVectorDualH = 1 / (solVectorDualH.transpose() * DWR.massH() * solVectorDualH) * solVectorDualH;
    // DWR.constructMultiPatchH(solVectorDualH, dualH);

    // // Swap multipatch
    // solVectorDualH *= sgn(DWR.matrixNorm(dualL, dualH));
    // DWR.constructMultiPatchH(solVectorDualH, dualH);

    gsInfo << "done.\n";

    gsField<> dualField(mp, dualH);
    gsWriteParaview(dualField, "dualH", 1000);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real_t approx = DWR.computeErrorEig(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL);

    real_t lambda_an = omegas[0];
    gsDebugVar(math::pow(lambda_an, 2));
    gsDebugVar(lambda_an);
    gsDebugVar(eigvalL);

    real_t exact = math::pow(lambda_an, 2) - eigvalL;

    gsInfo << "approx = " << approx << "\n";
    gsInfo << "Exact = " << exact << "\n";
    gsInfo << "Efficiency = " << approx / exact << "\n";

    return EXIT_SUCCESS;

} // end main
