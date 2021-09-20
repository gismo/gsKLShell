/** @file example_shell3D.cpp

    @brief Examples for the surface thin shells including the shell obstacle course

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#ifdef GISMO_WITH_SPECTRA
#include <gsSpectra/gsSpectra.h>
#endif

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsThinShellAssemblerDWR.h>
#include <gsKLShell/gsThinShellUtils.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template<typename T>
class gsElementErrorPlotter : public gsFunction<T>
{
public:
    gsElementErrorPlotter(const gsBasis<T>& mp, const std::vector<T>& errors ) : m_mp(mp),m_errors(errors)
    {

    }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
    {
        // Initialize domain element iterator -- using unknown 0
        res.setZero(1,u.cols());
        for(index_t i=0; i<u.cols();++i)
        {
            int iter =0;
            // Start iteration over elements

            typename gsBasis<T>::domainIter domIt = m_mp.makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                 bool flag = true;
                const gsVector<T>& low = domIt->lowerCorner();
                const gsVector<T>& upp = domIt->upperCorner();


                for(int d=0; d<domainDim();++d )
                {
                    if(low(d)> u(d,i) || u(d,i) > upp(d))
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                     res(0,i) = m_errors.at(iter);
                     break;
                }
                iter++;
            }
        }
    }

    short_t domainDim() const { return m_mp.dim();}

private:
    const gsBasis<T>& m_mp;
    const std::vector<T>& m_errors;
};


template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

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
    bool write = false;
    index_t numRefine = 1;
    index_t numElevate = 1;
    bool loop = false;
    bool adaptive = false;
    std::string fn;

    real_t E_modulus = 1e6;
    real_t PoissonRatio = 0.3;
    real_t thickness = 1e-2;

    real_t aDim = 1.0;
    real_t bDim = 1.0;

    index_t modeIdx = 0;

    refCriterion = GARU;
    refParameter = 0.85;
    RefineLoopMax = 1;

    int testCase = 0;

    int refExt = -1;
    int crsExt = -1;

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
    cmd.addInt("E", "refExt", "Refinement extension", refExt);
    cmd.addInt("C", "crsExt", "Coarsening extension", crsExt);
    cmd.addReal("T", "thickness", "thickness", thickness);

    cmd.addReal("a","adim", "dimension a", aDim);
    cmd.addReal("b","bdim", "dimension b", bDim);

    cmd.addString("f", "file", "Input XML file", fn);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence to file", write);
    cmd.addSwitch("loop", "Uniform Refinemenct loop", loop);

    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int rv)
    {
        return rv;
    }
    //! [Parse command line]

    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    gsVector<> pts(2);
    pts.setConstant(0.25);

    if (testCase==0)
    {
        // Unit square
        mp.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        mp.patch(0).coefs().col(0) *= aDim;
        mp.patch(0).coefs().col(1) *= bDim;
        E_modulus = 1e6;
        PoissonRatio = 0.3;
        thickness = 1e-2;
    }
    else if (testCase==1)
    {
        // Unit square
        mp.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        mp.patch(0).coefs().col(0) *= aDim;
        mp.patch(0).coefs().col(1) *= bDim;
        E_modulus = 1e6;
        PoissonRatio = 0.3;
        thickness = 1e-2;
    }
    else
        GISMO_ERROR("Test case" << testCase<<" unknown!");

    // Unit square
    mp = Rectangle(aDim,bDim);

    // p-refine
    if (numElevate != 0)
        mp.degreeElevate(numElevate);

    // h-refine
    if (!loop)
    {
        for (index_t r =0; r < numRefine; ++r)
            mp.uniformRefine();
        numRefine = 0;
    }

    // Cast all patches of the mp object to THB splines
    if (refExt!=-1 || crsExt!=-1)
    {
        gsTHBSpline<2,real_t> thb;
        for (index_t k=0; k!=mp.nPatches(); ++k)
        {
            gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
            thb = gsTHBSpline<2,real_t>(*geo);
            mp.patch(k) = thb;
        }
    }

    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: " << mp.nPatches() << ", degree: " << dbasis.minCwiseDegree() << "\n";
    gsInfo << dbasis.basis(0) << "\n";

    mp_def = mp;

    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH = basisL;
    basisH.degreeElevate(1);

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();


    real_t D = E_modulus * math::pow(thickness, 3) / (12 * (1 - math::pow(PoissonRatio, 2)));
    std::vector<real_t> gammas;
    real_t Load = 0;
    real_t lambda_an = 0;
    if (testCase==0)
    {
        Load = 1e-1;
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 1.0, 0.0 ; load << Load, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );

        // // Clamped-Clamped
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 0 - x

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        // bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z

        //   [Analytical solution]
        gammas.resize(4);
        gammas[0] = 0.7692307692;
        gammas[1] = 1.453488372;
        gammas[2] = 2.688172043;
        gammas[3] = 4.432515337;
        real_t sigma1 = 3.61523970735874E+02;
        lambda_an = gammas[modeIdx]*sigma1*thickness / (Load);
        // ! [Analytical solution]
    }
    else if (testCase==1)
    {
        real_t Load = 1e3;

        gsVector<> point(2);
        gsVector<> load (3);
        point<< 0.5, 1.0 ; load << Load, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );

        bc.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 0); // unknown 0 - x
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0); // unknown 0 - x
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z
        lambda_an = 0;
    }

    gsConstantFunction<> force(tmp, 3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus), 3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio), 3);

    std::vector<gsFunction<> *> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrixBase<real_t> *materialMatrix;
    gsOptionList options;
    options.addInt("Material", "Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden", 0);
    options.addInt("Implementation", "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
    materialMatrix = getMaterialMatrix<3, real_t>(mp, t, parameters, options);

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
    std::vector<real_t> DoFs(numRefine+1);

    // solvers
    gsSparseSolver<>::LU solver;
    Eigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> eigSolver;

    // solutions
    gsMultiPatch<> primalL, dualL, dualH;
    gsVector<> solVector, solVectorLinear, solVectorDualL, solVectorDualH;
    real_t eigvalL, dualvalL, dualvalH;

    // matrix norm
    real_t Mnorm;

    // matrices
    gsSparseMatrix<> K_L, K_NL, Kdiff;
    gsVector<> rhs;

    // DWR assembler
    gsThinShellAssemblerDWRBase<real_t> * DWR;

    gsParaviewCollection collection("solution");
    gsParaviewCollection errors("error_elem_ref");

    std::vector<real_t> elErrors;
    std::vector<bool> refVec;
    MarkingStrategy adaptRefCrit = BULK;
    real_t adaptRefParam;

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
        DWR->setPointLoads(pLoads);
        DWR->setGoal(GoalFunction::Buckling);

        gsInfo << "Computing load step... " << std::flush;
        DWR->assembleMatrixL();
        DWR->assemblePrimalL();
        K_L =  DWR->matrixL();
        rhs = DWR->primalL();
        solver.compute(K_L);
        solVectorLinear = solver.solve(rhs);

        DWR->constructSolutionL(solVectorLinear, mp_def);
        gsInfo << "done\n";

        gsMultiPatch<> deformed = mp_def;
        deformed.patch(0).coefs() -= mp.patch(0).coefs();
        gsField<> fielddef(mp, deformed);
        gsWriteParaview(fielddef,"deformed");

        gsInfo << "Assembling primal... " << std::flush;
        DWR->assembleMatrixL(mp_def);
        K_NL = DWR->matrixL();
        gsInfo << "done\n";

        // Solve system
        gsInfo << "Solving primal, size =" << DWR->matrixL().rows() << "," << DWR->matrixL().cols() << "... " << std::flush;

#ifdef GISMO_WITH_SPECTRA
        index_t num = std::min(K_L.cols()-1,10);
        Kdiff = K_NL - K_L;
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> eigsolverL(K_L,Kdiff,num,2*num,0.0);
        eigsolverL.init();
        eigsolverL.compute(Spectra::SortRule::LargestMagn,1000,1e-6,Spectra::SortRule::SmallestMagn);

        if (eigsolverL.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<eigsolverL.num_iterations()<<" iterations and with "<<eigsolverL.num_operations()<<"operations. \n"; }
        else if (eigsolverL.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
        else if (eigsolverL.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
        else if (eigsolverL.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
        else                                                      { GISMO_ERROR("No error code known"); }
#else
        Eigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> eigsolverL;
        eigsolverL.compute(K_L,K_NL-K_L);
#endif

        if (modeIdx > eigsolverL.eigenvalues().size()-1)
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

        gsDebugVar(eigsolverL.eigenvalues()[modeIdx]*Load);
        solVector = solVectorDualL = eigsolverL.eigenvectors().col(modeIdx);
        eigvalL = dualvalL = eigsolverL.eigenvalues()[modeIdx];

        // eigSolver.compute(K_L, K_NL-K_L);
        // gsDebugVar(eigSolver.eigenvalues()[modeIdx]*Load);
        // solVector = solVectorDualL = eigSolver.eigenvectors().col(modeIdx);
        // eigvalL = dualvalL = eigSolver.eigenvalues()[modeIdx];

        DWR->constructMultiPatchL(solVector, primalL);

        Mnorm = DWR->matrixNorm(primalL, primalL,mp_def);
        gsDebugVar(solVector.transpose() * (K_NL-K_L) * (solVector));
        gsDebugVar(Mnorm);

        // Mass-normalize primal
        solVector *= 1 / Mnorm;
        DWR->constructMultiPatchL(solVector, primalL);

        // mass-normalize w.r.t. primal
        DWR->constructMultiPatchL(solVectorDualL, dualL);
        Mnorm = DWR->matrixNorm(primalL, dualL,mp_def);
        gsDebugVar(solVector.transpose() * (K_NL-K_L) * (solVectorDualL));
        gsDebugVar(Mnorm);
        solVectorDualL *= 1. / Mnorm;
        DWR->constructMultiPatchL(solVectorDualL, dualL);

        gsInfo << "done.\n";

        gsDebugVar(solVector.transpose() * (K_NL-K_L) * (solVectorDualL));
        // gsDebugVar(solVector.transpose() * K_NL * (solVectorDualL - solVector));

        gsInfo << "Assembling dual matrix (H)... " << std::flush;
        DWR->assembleMatrixH();
        K_L = DWR->matrixH();
        DWR->assembleMatrixH(mp_def);
        K_NL = DWR->matrixH();
        gsInfo << "done.\n";

        gsInfo << "Solving dual (high), size = " << DWR->matrixH().rows() << "," << DWR->matrixH().cols() << "... " << std::flush;

#ifdef GISMO_WITH_SPECTRA
        Kdiff = K_NL - K_L;
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> eigsolverH(K_L,Kdiff,num,2*num,0.0);
        eigsolverH.init();
        eigsolverH.compute(Spectra::SortRule::LargestMagn,1000,1e-6,Spectra::SortRule::SmallestMagn);

        if (eigsolverH.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<eigsolverH.num_iterations()<<" iterations and with "<<eigsolverH.num_operations()<<"operations. \n"; }
        else if (eigsolverH.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
        else if (eigsolverH.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
        else if (eigsolverH.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
        else
#else
        Eigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> eigsolverH;
        eigsolverH.compute(K_L,K_NL-K_L);
#endif

        gsDebugVar(eigsolverH.eigenvalues()[modeIdx]*Load);
        solVectorDualH = eigsolverH.eigenvectors().col(modeIdx);
        dualvalH = eigsolverH.eigenvalues()[modeIdx];
        //
        // eigSolver.compute(K_L, K_NL-K_L);
        // gsDebugVar(eigSolver.eigenvalues()[modeIdx]*Load);
        // solVectorDualH = eigSolver.eigenvectors().col(modeIdx);
        // dualvalH = eigSolver.eigenvalues()[modeIdx];

        // mass-normalize w.r.t. primal
        DWR->constructMultiPatchH(solVectorDualH, dualH);
        Mnorm = DWR->matrixNorm(primalL, dualH,mp_def);
        gsDebugVar(Mnorm);
        solVectorDualH *= 1. / Mnorm;
        DWR->constructMultiPatchH(solVectorDualH, dualH);
        gsInfo << "done.\n";

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (plot)
        {
            gsField<> VMStress(mp,primalL, true);
            std::string fileName = "solution" + util::to_string(r);
            gsWriteParaview<>(VMStress, fileName, 5000, true);
            fileName = "solution" + util::to_string(r) + "0";
            collection.addTimestep(fileName,r,".vts");
            collection.addTimestep(fileName,r,"_mesh.vtp");
        }

        exacts[r] = 0;
        numGoal[r] = eigvalL;
        exGoal[r] = lambda_an;
        DoFs[r] = basisL.basis(0).numElements();

        exacts[r] += exGoal[r];
        exacts[r] -= numGoal[r];
        approxs[r] = DWR->computeErrorEig(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL,mp_def);

        estGoal[r] = numGoal[r]+approxs[r];

        efficiencies[r] = approxs[r]/exacts[r];

        if (refExt==-1 && crsExt==-1)
        {
            mp.uniformRefine();
        }
        else
        {
            elErrors = DWR->computeErrorEigElements(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL);
            for (std::vector<real_t>::iterator it = elErrors.begin(); it != elErrors.end(); it++)
            {
                *it = std::abs(*it);
                gsDebugVar(*it);
            }

            gsElementErrorPlotter<real_t> err_eh(mp.basis(0),elErrors);
            const gsField<> elemError_eh( mp.patch(0), err_eh, true );
            gsWriteParaview<>( elemError_eh, "error_elem_ref" + util::to_string(r), 1000, true);
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,".vts");
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,"_mesh.vtp");

            adaptRefParam = 0.25;
            std::vector<bool> elMarked( elErrors.size() );
            gsMarkElementsForRef( elErrors, adaptRefCrit, adaptRefParam, elMarked);

            if (refExt!=-1 && crsExt==-1)
            {
                gsRefineMarkedElements( mp, elMarked,refExt );

                gsInfo<<"Error\tRefined?\n";
                for (index_t k=0; k!=elMarked.size(); k++)
                    gsInfo<<elErrors[k]<<"\t"<<elMarked[k]<<"\n";
            }
            else if (refExt!=-1 && crsExt!=-1)
            {
                adaptRefParam = 0.1;
                // Invert errors for coarsening marking
                std::vector<real_t> elErrorsC = elErrors;
                for (index_t k=0; k!=elErrors.size(); k++)
                    elErrorsC[k] = -elErrors[k];

                std::vector<bool> elCMarked( elErrorsC.size() );
                gsMarkElementsForRef( elErrorsC, adaptRefCrit, adaptRefParam, elCMarked);
                gsProcessMarkedElements( mp, elMarked, elCMarked, 1, 0 );

                gsInfo<<"Error\tRefined?\tCoarsened?\n";
                for (index_t k=0; k!=elMarked.size(); k++)
                    gsInfo<<elErrors[k]<<"\t"<<elMarked[k]<<"\t"<<elCMarked[k]<<"\n";
            }
        }
        mp_def = mp;
    }

    if (plot)
    {
        collection.save();
        errors.save();
    }

    gsInfo<<"-------------------------------------------------------------------------------------------------\n";
    gsInfo<<"Ref.\tApprox    \tExact     \tEfficiency\tNumGoal   \tEstGoal   \texGoal    \t#elements \n";
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
        filename = "example_shell3D_DWR_buckling_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "_I" + std::to_string(modeIdx);
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

template <class T>
gsMultiPatch<T> Rectangle(T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,2,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,2,1);

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

