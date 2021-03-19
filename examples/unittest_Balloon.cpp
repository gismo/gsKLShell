/** @file unittest_Balloon.cpp

    @brief inflated balloon unit test

    Test follower pressures

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

std::pair<real_t,real_t> numerical(index_t material, index_t impl)
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

    gsReadFile<>("surface/eighth_sphere.xml", mp);

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

    GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
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

    gsMatrix<> pt(2,2);
    pt.col(0)<<0.0,1.0;
    pt.col(1)<<1.0,1.0;


    gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

    real_t tol = 10e-10;
    GISMO_ASSERT((lambdas.col(0)-lambdas.col(1)).norm() < tol, "Stretches must be equal over the balloon");

    real_t r = (mp_def.patch(0).eval(pt).col(0)).norm();

    // Get the total force on the tension boundary
    real_t P = pressure * assembler->getArea(mp) / assembler->getArea(mp_def);

    delete assembler;

    std::pair<real_t,real_t> result;
    result.first = P;
    result.second = r;

    return result;
}

real_t analytical(index_t material, index_t impl, real_t r)
{
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


// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{

    real_t P, Pan, rnum;

    std::vector<index_t> materials{ 1,3,4 };
    std::vector<index_t> implementations{ 1,2,3 };

    std::pair<real_t,real_t> num;
    real_t tol = 1e-3;
    for (std::vector<index_t>::iterator mat = materials.begin(); mat!=materials.end(); mat++)
    {
        for (std::vector<index_t>::iterator impl = implementations.begin(); impl!=implementations.end(); impl++)
        {
            if (*mat==4 && *impl!=3) continue;

            num = numerical(*mat,*impl);
            P = num.first;
            rnum = num.second;

            Pan = analytical(*mat,*impl,rnum);


            if (std::abs(P-Pan)/Pan < tol)
                gsInfo<<"Passed\n";
            else
                gsInfo<<"Failed; L error = "<<std::abs(P-Pan)/Pan<<"\n";
        }
    }
return EXIT_SUCCESS;

}// end main