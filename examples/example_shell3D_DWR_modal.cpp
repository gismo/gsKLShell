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

template<typename T>
index_t sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

using namespace gismo;

int main(int argc, char *argv[])
{
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
    index_t numRefineIni = 0;
    index_t numElevate = 1;
    bool loop = false;
    bool adaptive = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.3;
    real_t Density = 1.0;
    real_t thickness = 0.01;

    index_t modeIdx = 0;

    refCriterion = GARU;
    refParameter = 0.85;

    int testCase = 0;

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
    cmd.addReal("T", "thickness", "thickness", thickness);
    cmd.addString("f", "file", "Input XML file", fn);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence to file", write);
    cmd.addSwitch("loop", "Uniform Refinemenct loop", loop);
    cmd.addSwitch("adaptive", "Adaptive Refinemenct loop", adaptive);

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

    // p-refine
    if (numElevate != 0)
        mp.degreeElevate(numElevate);

    for (index_t r =0; r < numRefineIni; ++r)
        mp.uniformRefine();

    // h-refine
    if (!loop)
    {
        for (index_t r =0; r < numRefine; ++r)
            mp.uniformRefine();
        numRefine = 0;
    }

    // Cast all patches of the mp object to THB splines
    if (adaptive)
    {
        if (testCase!=1)
        {
            gsTHBSpline<2,real_t> thb;
            for (index_t k=0; k!=mp.nPatches(); ++k)
            {
                gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
                thb = gsTHBSpline<2,real_t>(*geo);
                mp.patch(k) = thb;
            }
        }
        else
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
        gammas<<3.1962206158252, 4.6108999, 4.6108999, 5.9056782, 5.9056782, 6.3064370, 7.1435310, 7.1435310;//, 7.7992738, 8.3466059, 9.1968826, 9.4394991, 9.5257014;
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

        omegas.push_back(0);
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

    std::vector<gsFunction<> *> parameters(2);
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

    // solvers
    gsSparseSolver<>::LU solver;

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

    gsParaviewCollection collection("solution");
    gsParaviewCollection errors("error_elem_ref");

    std::vector<real_t> elErrors;
    std::vector<bool> refVec;
    MarkingStrategy adaptRefCrit = PUCA;
    const real_t adaptRefParam = 0.8;

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

        // Solve system
        gsInfo << "Solving primal, size =" << DWR->matrixL().rows() << "," << DWR->matrixL().cols() << "... " << std::flush;

#ifdef GISMO_WITH_SPECTRA
        index_t numL = std::min(DWR->matrixL().cols()-1,10);
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> eigsolverL(DWR->matrixL(),DWR->massL(),numL,2*numL,0);
        eigsolverL.init();
        eigsolverL.compute(Spectra::SortRule::LargestMagn,1000,1e-10,Spectra::SortRule::SmallestAlge);
        if (eigsolverL.info()==Spectra::CompInfo::Successful)           { gsDebug<<"Spectra converged in "<<eigsolverL.num_iterations()<<" iterations and with "<<eigsolverL.num_operations()<<"operations. \n"; }
        else if (eigsolverL.info()==Spectra::CompInfo::NumericalIssue)  { GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
        else if (eigsolverL.info()==Spectra::CompInfo::NotConverging)   { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
        else if (eigsolverL.info()==Spectra::CompInfo::NotComputed)     { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
        else                                                            { GISMO_ERROR("No error code known"); }
#else
        Eigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> eigsolverL;
        eigsolverL.compute(DWR->matrixL(), DWR->massL());
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

        solVector = solVectorDualL = eigsolverL.eigenvectors().col(modeIdx);

        eigvalL = dualvalL = eigsolverL.eigenvalues()[modeIdx];

        // Mass-normalize primal
        solVector = 1 / (solVector.transpose() * DWR->massL() * solVector) * solVector;
        DWR->constructMultiPatchL(solVector, primalL);
        DWR->constructSolutionL(solVector, mp_def);

        // mass-normalize w.r.t. primal
        DWR->constructMultiPatchL(solVectorDualL, dualL);
        Mnorm = DWR->matrixNorm(primalL, dualL);
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

#ifdef GISMO_WITH_SPECTRA
        index_t numH = std::min(DWR->matrixH().cols()-1,10);
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> eigsolverH(DWR->matrixH(),DWR->massH(),numH,2*numH,0);
        eigsolverH.init();
        eigsolverH.compute(Spectra::SortRule::LargestMagn,1000,1e-10,Spectra::SortRule::SmallestAlge);
        if (eigsolverH.info()==Spectra::CompInfo::Successful)           { gsDebug<<"Spectra converged in "<<eigsolverH.num_iterations()<<" iterations and with "<<eigsolverH.num_operations()<<"operations. \n"; }
        else if (eigsolverH.info()==Spectra::CompInfo::NumericalIssue)  { GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
        else if (eigsolverH.info()==Spectra::CompInfo::NotConverging)   { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
        else if (eigsolverH.info()==Spectra::CompInfo::NotComputed)     { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
        else                                                            { GISMO_ERROR("No error code known"); }
#else
        Eigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> eigsolverH;
        eigsolverH.compute(DWR->matrixH(), DWR->massH());
#endif

        solVectorDualH = eigsolverH.eigenvectors().col(modeIdx);
        dualvalH = eigsolverH.eigenvalues()[modeIdx];

        // mass-normalize w.r.t. primal
        DWR->constructMultiPatchH(solVectorDualH, dualH);
        Mnorm = DWR->matrixNorm(primalL, dualH);
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
            std::string fileName = "solution" + util::to_string(r);
            gsWriteParaview<>(VMStress, fileName, 5000, true);
            fileName = "solution" + util::to_string(r) + "0";
            collection.addTimestep(fileName,r,".vts");
            collection.addTimestep(fileName,r,"_mesh.vtp");
        }

        exacts[r] = 0;
        numGoal[r] = eigvalL;
        exGoal[r] = math::pow(lambda_an, 2);
        DoFs[r] = basisL.basis(0).numElements();

        exacts[r] += exGoal[r];
        exacts[r] -= numGoal[r];
        approxs[r] = DWR->computeErrorEig(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL);

        estGoal[r] = numGoal[r]+approxs[r];

        efficiencies[r] = approxs[r]/exacts[r];

        if (!adaptive)
            mp.uniformRefine();
        else
        {
            elErrors = DWR->computeErrorEigElements(eigvalL, dualvalL, dualvalH, dualL, dualH, primalL);
            for (std::vector<real_t>::const_iterator it = elErrors.begin(); it != elErrors.end(); it++)
                gsDebug<<*it<<"\n";

            // std::vector<real_t> errorsDofs = DWR->computeErrorDofs(dualL,dualH);
            // for (std::vector<real_t>::const_iterator it = errorsDofs.begin(); it != errorsDofs.end(); it++)
            //     gsDebug<<*it<<"\n";


            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( elErrors.size() );
            gsMarkElementsForRef( elErrors, adaptRefCrit, adaptRefParam, elMarked);

            // Invert errors for coarsening marking
            std::vector<real_t> elErrorsC = elErrors;
            for (index_t k=0; k!=elErrors.size(); k++)
                elErrorsC[k] = -elErrors[k];

            std::vector<bool> elCMarked( elErrorsC.size() );
            gsMarkElementsForRef( elErrorsC, adaptRefCrit, adaptRefParam, elCMarked);


            gsElementErrorPlotter<real_t> err_eh(mp.basis(0),elErrors);
            const gsField<> elemError_eh( mp.patch(0), err_eh, true );
            gsWriteParaview<>( elemError_eh, "error_elem_ref" + util::to_string(r), 1000, true);
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,".vts");
            errors.addTimestep("error_elem_ref" + util::to_string(r) + "0",r,"_mesh.vtp");

            gsProcessMarkedElements( mp, elMarked, elCMarked, 1, 0 );
            // gsRefineMarkedElements( mp, elMarked,1 );
            mp_def = mp;

            for (index_t k=0; k!=elMarked.size(); k++)
            {
                gsInfo<<elErrors[k]<<"\t"<<elMarked[k]<<"\n";
                // gsInfo<<elErrors[k]<<"\t"<<elMarked[k]<<"\t"<<elCMarked[k]<<"\n";
            }
        }
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
        filename = "example_shell3D_DWR_modal_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "_I" + std::to_string(modeIdx);
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
