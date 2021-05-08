/** @file unittest_UAT.cpp

    @brief Unit tests performs Uniaxial Tension Test for Neo-Hookean, Mooney-Rivlin and Ogden material models

    This file tests the following classes and functions:
    - gsMaterialMatrix          (dim=2, mat=1,3,4, impl=1,2,3)
    - gsMaterialMatrixIntegrate (dim=2)
    - gsMaterialMatrixBase
    - gsThinShellAssembler      (dim=2)
        - assemble(), assembleMatrix(), assembleVector()
        - boundaryFoceVector(), getArea()
        - constructSolution(), computePrincipalStretches()

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

using namespace gismo;

std::pair<real_t,real_t> numerical(index_t material, index_t impl, bool Compressibility)
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

    real_t lambda = 2.0;
    gsConstantFunction<> displx(lambda-1.0,2);

    GISMO_ASSERT(mp.targetDim()==2,"Geometry must be planar (targetDim=2)!");
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

    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;
    gsMaterialMatrixBase<real_t>* materialMatrix;

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

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2, real_t, false >(mp,dbasis,bc,force,materialMatrix);

    assembler->setPointLoads(pLoads);

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
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
    for (index_t it = 0; it != 100; ++it)
    {
        jacMat = Jacobian(solVector);
        solver.compute(jacMat);
        updateVector = solver.solve(resVec); // this is the UPDATE
        solVector += updateVector;

        resVec = Residual(solVector);

        if (updateVector.norm() < 1e-6)
            break;
        else if (it+1 == it)
            gsWarn<<"Maximum iterations reached!\n";
    }

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Check solutions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NOTE: all analytical solutions for compressible materials are fixed for displ=1; (lambda=2)

    // Compute stretches (should be the same everywhere)
    // Ordering: lambda(0) < lambda(1); lambda(2) is ALWAYS the through-thickness stretch
    gsVector<> pt(2);
    pt<<1,0;
    gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

    // Get the total force on the tension boundary
    patchSide ps(0,boundary::east);
    gsMatrix<> forceVector = assembler->boundaryForceVector(mp_def,ps,0);
    real_t sideForce = forceVector.sum();
    real_t S   = -sideForce / (thickness*lambdas(0)*lambdas(2));
    real_t L   = lambdas(0);

    std::pair<real_t,real_t> result;
    result.first = L;
    result.second = S;

    delete assembler;

    return result;
}

std::pair<real_t,real_t> analytical(index_t material, index_t impl, bool Compressibility)
{
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

    real_t lambda = 2.0;

    real_t San,J,K,Lan;
    if      (material==1 && Compressibility)
    {
        K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
        J = 1.105598565;// specific for lambda==2!!
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
        J = 1.099905842;// specific for lambda==2!!
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
        J = 1.088778638;// specific for lambda==2!!
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

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{

real_t S, San,L,Lan;

std::vector<index_t> materials{ 1,3,4 };
std::vector<index_t> implementations{ 1,2,3 };
std::vector<bool> compressibility{ true,false };

std::pair<real_t,real_t> num, an;
real_t tol = 1e-9;
for (std::vector<index_t>::iterator mat = materials.begin(); mat!=materials.end(); mat++)
{
    for (std::vector<index_t>::iterator impl = implementations.begin(); impl!=implementations.end(); impl++)
    {
        for (std::vector<bool>::iterator comp = compressibility.begin(); comp!=compressibility.end(); comp++)
        {
            if (*mat==4 && *impl!=3) continue;
            num = numerical(*mat,*impl,*comp);
            L = num.first;
            S = num.second;

            an = analytical(*mat,*impl,*comp);
            Lan = an.first;
            San = an.second;


            if ( (std::abs(L-Lan)/Lan < tol) && (std::abs(S-San)/San < tol) )
                gsInfo<<"Passed\n";
            else
                gsInfo<<"Failed; L error = "<<std::abs(L-Lan)/Lan<<"\t S error = "<<std::abs(S-San)/San<<"\n";
        }
    }
}

return EXIT_SUCCESS;

}// end main