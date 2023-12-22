/** @file example_shell3D.cpp

    @brief Examples for the surface thin shells including the shell obstacle course

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>
#include <typeinfo>

#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsThinShellAssemblerDWR.h>
#include <gsKLShell/src/gsThinShellUtils.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

template<typename T>
index_t sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool write = false;
    index_t numRefine = 1;
    index_t numRefineIni = 0;
    index_t numElevate = 1;
    bool last = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.3;
    real_t Density = 1.0;
    real_t thickness = 0.01;

    index_t modeIdx = 0;

    int testCase = 1;

    int adaptivity = 0;

    std::string mesherOptionsFile("options/shell_mesher_options.xml");
    std::string dirname  = "ModalResults";

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt("i", "index", "index of mode", modeIdx);
    cmd.addInt("e", "degreeElevation",
               "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate);
    cmd.addInt("R", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", numRefineIni);
    cmd.addInt("r", "refine", "Maximum number of adaptive refinement steps to perform",
               numRefine);
    cmd.addInt("t", "testcase",
                "Test case: 0: Beam - pinned-pinned, 1: Beam - fixed-fixed, 2: beam - fixed-free, 3: plate - fully pinned, 4: plate - fully fixed, 5: circle - fully pinned, 6: 5: circle - fully fixed",
               testCase);
    cmd.addInt("A", "adaptivity", "Adaptivity scheme: 0) uniform refinement, 1) adaptive refinement, 2) adaptive refinement and coarsening", adaptivity);

    cmd.addString("f", "file", "Input XML file", fn);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence to file", write);
    cmd.addSwitch("last", "Only last refinement", last);
    cmd.addString( "O", "mesherOpt", "Input XML file for mesher options", mesherOptionsFile );
    cmd.addString( "o", "output", "output directory", dirname );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

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
        numElevate -= 1;
    }
    else if (testCase==2)
    {
        mp.addPatch(gsNurbsCreator<>::BSplineTriangle(1,1));
        mp.addAutoBoundaries();
        mp.embed(3);
        mp.degreeElevate(1);

        thickness = 0.01;
        PoissonRatio = 0.3;
        E_modulus     = 1e0;
        Density = 1e0;
        numElevate -= 1;
    }
    else if (testCase==3)
    {
        gsMultiPatch<> mp_tmp;
        std::string fn = "planar/weirdShape.xml";
        gsReadFile<>(fn, mp_tmp);
        mp.addAutoBoundaries();

        gsKnotVector<> kv1(0,1,8,3);
        gsKnotVector<> kv2(0,1,0,3);
        gsTensorBSplineBasis<2,real_t> basis(kv2,kv1);
        gsMatrix<> coefs;
        gsQuasiInterpolate<real_t>::localIntpl(basis, mp_tmp.patch(0), coefs);
        gsTensorBSpline<2,real_t> bspline(kv2,kv1,coefs);

        gsDebugVar(mp_tmp.patch(0).coefs());
        gsDebugVar(coefs);
        gsDebugVar(bspline);

        mp.addPatch(bspline);
        mp.embed(3);
        mp.degreeElevate(1);

        thickness = 0.01;
        PoissonRatio = 0.3;
        E_modulus     = 1e0;
        Density = 1e0;

        gsTensorBSpline<2,real_t> *geo;
        gsDebug<<(geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(0)))<<"\n";

        gsDebug<<(geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&bspline))<<"\n";
    }

    // p-refine
    if (numElevate != 0)
        mp.degreeElevate(numElevate);

    // Cast all patches of the mp object to THB splines
    if (adaptivity!=0)
    {
        if (testCase!=1)
        {
            gsTHBSpline<2,real_t> thb;
            for (index_t k=0; k!=mp.nPatches(); ++k)
            {
                gsTensorBSpline<2,real_t> *geo;
                if (geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k)))
                {
                    thb = gsTHBSpline<2,real_t>(*geo);
                    gsMatrix<> bbox = geo->support();
                    for (index_t i = 0; i< numRefineIni; ++i)
                        thb.refineElements(thb.basis().asElements(bbox));
                    mp.patch(k) = thb;
                }
                else
                {
                    std::cout<<typeid(geo).name();
                    GISMO_ERROR("Cannot cast to THB spline");
                }
            }
        }
        else // Quasi-interpolate NURBS geometry
        {
            gsMultiPatch<> mp_thb;
            gsTHBSpline<2,real_t> thb;
            gsTHBSplineBasis<2,real_t> thbBasis;
            gsMatrix<> coefs;
            for (index_t k=0; k!=mp.nPatches(); ++k)
            {
                if(gsTensorNurbs<2,real_t> *geo = dynamic_cast< gsTensorNurbs<2,real_t> * > (&mp.patch(k)))
                {
                    thbBasis = gsTHBSplineBasis<2,real_t>(geo->basis().source());
                    gsQuasiInterpolate<real_t>::localIntpl(thbBasis, mp.patch(k), coefs);
                    thb = gsTHBSpline<2,real_t>(thbBasis,coefs);
                    gsMatrix<> bbox = geo->support();
                    for (index_t i = 0; i< numRefineIni; ++i)
                        thb.refineElements(thb.basis().asElements(bbox));
                    mp_thb.addPatch(thb);
                }
            }

            gsField<> field1(mp,mp);
            gsField<> field2(mp,mp_thb);

            gsInfo<<"THB Approximation error: "<<field1.distanceL2(field2)<<"\n";
            gsWriteParaview<>(mp_thb,"mp_thb",10000,true);

            mp = mp_thb;
            // GISMO_ERROR("Adaptivity not available for NURBS");
        }
    }
    else
    {
        for (index_t i = 0; i< numRefineIni; ++i)
            mp.uniformRefine();
    }

    // h-refine
    if (last)
    {
        for (index_t r =0; r < numRefine; ++r)
            mp.uniformRefine();
        numRefine = 0;
    }

    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: " << mp.nPatches() << ", degree: " << dbasis.minCwiseDegree() << "\n";
    gsInfo << dbasis.basis(0) << "\n";
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
        gsVector<> gammas(8);
        gammas<<3.1962206165825410939805274034037203415990811116202,
                4.6108998790490558272421332787278913722925337365292,
                4.6108998790490558272421332787278913722925337365292,
                5.9056782354205228796795933090304630477909550533454,
                5.9056782354205228796795933090304630477909550533454,
                6.3064370476884237158917750270645901577640872572062,
                7.1435310235048408654996353878889350576865598730402,
                7.1435310235048408654996353878889350576865598730402;
                // 7.7992738008112319024715415103761572702583726170654,
                // 8.3466059387507383472419351935644166660470955569012,
                // 9.1968825996353207360388779515890837406166652434270,
                // 9.4394991378764049051037582259363913219818055853008,
                // 9.5257013556717228655724356734908469832430657314660;
        for (index_t n=0; n!=gammas.size(); n++)
          omegas.push_back(math::pow(math::pow(gammas[n],4)*D/(Density*thickness),0.5));
    }
    else if (testCase==2)
    {
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

        omegas.resize(modeIdx+1);
        omegas[modeIdx] = 0;
    }
    else if (testCase==3)
    {
        for (index_t i = 0; i != 3; ++i)
        {
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i);
        }
        bc.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

        omegas.resize(modeIdx+1);
        omegas[modeIdx] = 0;
    }
    else
        GISMO_ERROR("TESTCASE UNKNOWN!");

    //   [Analytical solution]
    real_t lambda_an = omegas[modeIdx];
    // ! [Analytical solution]


    gsConstantFunction<> force(tmp, 3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus), 3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio), 3);
    gsFunctionExpr<> rho(std::to_string(Density), 3);

    std::vector<gsFunctionSet<> *> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrixBase<real_t> *materialMatrix;
    gsOptionList options;
    options.addInt("Material", "Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden", 0);
    options.addInt("Implementation", "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
    materialMatrix = getMaterialMatrix<3, real_t>(mp, t, parameters, rho, options);

    gsMatrix<> points(2, 0);
    // points.col(0).setConstant(0.25);
    // points.col(1).setConstant(0.50);
    // points.col(2).setConstant(0.75);

    // measures
    std::vector<real_t> exacts(numRefine+1);
    std::vector<real_t> approxs(numRefine+1);
    std::vector<real_t> efficiencies(numRefine+1);
    std::vector<real_t> numGoal(numRefine+1);
    std::vector<real_t> estGoal(numRefine+1);
    std::vector<real_t> exGoal(numRefine+1);
    std::vector<index_t> DoFs(numRefine+1);

    // solutions
    gsMultiPatch<> primalL, dualL, dualH;
    gsVector<> solVector, solVectorDualL, solVectorDualH;
    real_t eigvalL, dualvalL, dualvalH;

    // matrix norm
    real_t Mnorm;

    // matrices
    gsSparseMatrix<> K_L, K_NL;
    gsVector<> rhs;

    // DWR assembler
    gsThinShellAssemblerDWRBase<real_t> * DWR;

    std::string commands = "mkdir -p " + dirname;
    const char *command  = commands.c_str();
    system(command);

    gsParaviewCollection collection(dirname + "/" + "solution");
    gsParaviewCollection errors_elem(dirname + "/" + "error_elem_ref");
    gsParaviewCollection errors(dirname + "/" + "errors");

    std::vector<real_t> elErrors;
    std::vector<bool> refVec;

    gsAdaptiveMeshing<real_t> mesher;
    if (adaptivity!=0)
    {
        gsFileData<> fd_mesher(mesherOptionsFile);
        gsOptionList mesherOpts;
        fd_mesher.getFirst<gsOptionList>(mesherOpts);

        mesher = gsAdaptiveMeshing<real_t>(mp);
        mesher.options() = mesherOpts;
        mesher.getOptions();
    }

    for (index_t r=0; r!=numRefine+1; r++)
    {
        // -----------------------------------------------------------------------------------------
        // ----------------------------Prepare basis------------------------------------------------
        // -----------------------------------------------------------------------------------------
        // Set deformed multipatch
        mp_def = mp;

        // Set bases
        gsMultiBasis<> basisL(mp);
        gsMultiBasis<> basisH(mp);
        basisH.degreeElevate(1);

        gsInfo<<"Basis Primal: "<<mp.basis(0)<<"\n";
        gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n"; //// different than the one above!!
        gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";

        // -----------------------------------------------------------------------------------------
        // ----------------------------DWR method---------------------------------------------------
        // -----------------------------------------------------------------------------------------
        DWR = new gsThinShellAssemblerDWR<3, real_t, true>(mp, basisL, basisH, bc, force, materialMatrix);
        DWR->setGoal(GoalFunction::Modal);


        gsInfo << "Assembling primal... " << std::flush;
        DWR->assembleMatrixL();
        DWR->assembleMassL();
        gsInfo << "done\n";

        gsVector<> eigenvalues;
        gsMatrix<> eigenvectors;
        // Solve system
        gsInfo << "Solving primal, size =" << DWR->matrixL().rows() << "," << DWR->matrixL().cols() << "... " << std::flush;
#ifdef gsSpectra_ENABLED
        index_t numL = std::min(DWR->matrixL().cols()-1,10);
        gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solverL(DWR->matrixL(),DWR->massL(),numL,2*numL,0.0);
        solverL.init();
        solverL.compute(Spectra::SortRule::LargestMagn,1000,1e-30,Spectra::SortRule::SmallestMagn);
#else
        gsEigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  solverL;
        solverL.compute(DWR->matrixL(),DWR->massL());
#endif
        eigenvalues = solverL.eigenvalues();
        eigenvectors = solverL.eigenvectors();
        if (modeIdx > eigenvalues.size()-1)
        {
            gsWarn<<"No error computed because mode does not exist (system size)!\n";
            approxs[r] = 0;
            exacts[r] = 0;
            efficiencies[r] = 0;
            numGoal[r] = 0;
            estGoal[r] = 0;
            exGoal[r] = 0;
            DoFs[r] = 0;
            mp.uniformRefine();
            continue;
        }

        solVector = solVectorDualL = eigenvectors.col(modeIdx);

        eigvalL = dualvalL = eigenvalues(modeIdx,0);

        // Mass-normalize primal
        solVector = 1 / (solVector.transpose() * DWR->massL() * solVector) * solVector;
        DWR->constructMultiPatchL(solVector, primalL);
        DWR->constructSolutionL(solVector, mp_def);

        // mass-normalize w.r.t. primal
        DWR->constructMultiPatchL(solVectorDualL, dualL);
        Mnorm = DWR->matrixNorm(primalL, dualL);
        gsDebugVar(Mnorm);
        solVectorDualL *= 1. / Mnorm;
        DWR->constructMultiPatchL(solVectorDualL, dualL);

        // solVectorDualL = 1 / (solVectorDualL.transpose() * DWR->massL() * solVectorDualL) * solVectorDualL;
        // DWR->constructMultiPatchL(solVectorDualL, dualL);

        gsInfo << "done.\n";

        gsInfo << "Assembling dual matrix (H)... " << std::flush;
        DWR->assembleMatrixH();
        DWR->assembleMassH();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (H), size = " << DWR->matrixH().rows() << "," << DWR->matrixH().cols() << "... " << std::flush;
#ifdef gsSpectra_ENABLED
        index_t numH = std::min(DWR->matrixL().cols()-1,10);
        gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solverH(DWR->matrixH(),DWR->massH(),numH,2*numL,0.0);
        solverH.init();
        solverH.compute(Spectra::SortRule::LargestMagn,1000,1e-30,Spectra::SortRule::SmallestMagn);
#else
        gsEigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  solverH;
        solverH.compute(DWR->matrixH(),DWR->massH());
#endif
        eigenvalues = solverH.eigenvalues();
        eigenvectors = solverH.eigenvectors();

        solVectorDualH = eigenvectors.col(modeIdx);
        dualvalH = eigenvalues(modeIdx,0);

        // mass-normalize w.r.t. primal
        DWR->constructMultiPatchH(solVectorDualH, dualH);
        Mnorm = DWR->matrixNorm(primalL, dualH);
        gsDebugVar(Mnorm);
        solVectorDualH *= 1. / Mnorm;
        DWR->constructMultiPatchH(solVectorDualH, dualH);

        // solVectorDualH = 1 / (solVectorDualH.transpose() * DWR->massH() * solVectorDualH) * solVectorDualH;
        // DWR->constructMultiPatchH(solVectorDualH, dualH);

        // // Swap multipatch
        // solVectorDualH *= sgn(DWR->matrixNorm(dualL, dualH));
        // DWR->constructMultiPatchH(solVectorDualH, dualH);

        gsInfo << "done.\n";

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (plot)
        {
            gsField<> VMStress(mp,primalL, true);
            std::string fileName = dirname + "/" + "solution" + util::to_string(r);
            gsWriteParaview<>(VMStress, fileName, 5000, true);
            fileName = "solution" + util::to_string(r) + "0";
            collection.addTimestep(fileName,r,".vts");
            collection.addTimestep(fileName,r,"_mesh.vtp");
        }

        exacts[r] = 0;
        numGoal[r] = eigvalL;
        exGoal[r] = math::pow(lambda_an, 2);
        DoFs[r] = DWR->numDofsL();

        exacts[r] += exGoal[r];
        exacts[r] -= numGoal[r];
        approxs[r] = DWR->computeErrorEig(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL);

        estGoal[r] = numGoal[r]+approxs[r];

        efficiencies[r] = approxs[r]/exacts[r];

        elErrors = DWR->computeErrorEigElements(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL,dirname + "/" + "errors" + util::to_string(r),10000,true,true);
        errors.addTimestep("errors" + util::to_string(r) + "0",r,".vts");
        errors.addTimestep("errors" + util::to_string(r) + "0",r,"_mesh.vtp");
        // for (std::vector<real_t>::iterator it = elErrors.begin(); it != elErrors.end(); it++)
        // {
        //     *it = std::abs(*it);
        // }

        gsElementErrorPlotter<real_t> err_eh(mp.basis(0),elErrors);
        const gsField<> elemError_eh( mp.patch(0), err_eh, true );
        gsWriteParaview<>( elemError_eh, dirname + "/" + "error_elem_ref" + util::to_string(r), 1000, true);
        errors_elem.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,".vts");
        if (adaptivity==0)
        {
            mp.uniformRefine();
        }
        else if (adaptivity > 0)
        {
            gsFileData<> fd_mesher(mesherOptionsFile);
            gsOptionList mesherOpts;
            fd_mesher.getFirst<gsOptionList>(mesherOpts);

            // Make container of the boxes
            gsHBoxContainer<2,real_t> markRef, markCrs;

            mesher.markRef_into(elErrors,markRef);
            mesher.refine(markRef);

            if (adaptivity>1)
            {
                mesher.markCrs_into(elErrors,markRef,markCrs);
                mesher.refine(markRef);
                mesher.unrefine(markCrs);
                // gsDebugVar(markCrs);
            }
            mesher.rebuild();
        }
        mp_def = mp;
    }

    if (plot)
    {
        collection.save();
        errors.save();
        errors_elem.save();
    }

    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    gsInfo<<"Ref.\tApprox    \tExact     \tEfficiency\tNumGoal   \tEstGoal   \texGoal    \t#DoFs \n";
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    for(index_t r=0; r!=numRefine+1; r++)
    {
        gsInfo  <<std::setw(4 )<<std::left<<r<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<approxs[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<exacts[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<efficiencies[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<numGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<estGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<exGoal[r]<<"\t";
        gsInfo  <<std::setw(10)<<std::left<<DoFs[r]<<"\n";
    }
    gsInfo<<"-------------------------------------------------------------------------------------------------\n";

    if (write)
    {
        std::string filename;
        filename = "example_shell3D_DWR_modal_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "_I" + std::to_string(modeIdx) + "_t" + std::to_string(testCase);
        filename = filename + ".csv";
        std::ofstream file_out;
        file_out.open (filename);

        file_out<<"Ref,Approx,Exact,Efficiency,NumGoal,EstGoal,exGoal,DoFs\n";
        for(index_t r=0; r!=numRefine+1; r++)
        {
            file_out<<r<<","<<approxs[r]<<","<<exacts[r]<<","<<efficiencies[r]<<","<<numGoal[r]<<","<<estGoal[r]<<","<<exGoal[r]<<","<<DoFs[r]<<"\n";
        }

        file_out.close();
    }

    if (plot)
    {
        gsField<> fieldDL(mp, dualL);
        gsField<> fieldDH(mp, dualH);

        gsField<> fieldPL(mp, primalL);


        gsWriteParaview<>( fieldDL, "dualL", 1000);
        gsWriteParaview<>( fieldDH, "dualH", 1000);
        gsWriteParaview<>( fieldPL, "primalL", 1000);
    }

    delete DWR;
    return EXIT_SUCCESS;

} // end main
