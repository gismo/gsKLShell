/** @file example_shell3D_CC.cpp

    @brief Static analysis of a pinned X-stiffened plate
           using embedded Euler-Bernoulli beams
           in a Kirchhoff-Love shell.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Authors: C. Chianese, H.M. Verhelst
*/

#include <gismo.h>
#include <gsKLShell/gsKLShell.h>
#include <gsAssembler/gsEmbeddingUtils.h>

using namespace gismo;

template <class T>
void assembleLinearEmbeddedBeam_into(gsThinShellAssemblerBase<T> *assembler,
                                     const gsMultiPatch<T> &mp_surf,
                                     const gsMultiPatch<T> &mp_curve,
                                     const gsMultiBasis<T> &mbasis_curve,
                                     const T &R,
                                     gsSparseMatrix<T> &matrix);

template <class T>
void assembleNonlinearEmbeddedBeam_into(gsThinShellAssemblerBase<T> *assembler,
                                        const gsMultiPatch<T> &mp_surf,
                                        const gsMultiPatch<T> &mp_surf_def,
                                        const gsMultiPatch<T> &mp_curve,
                                        const gsMultiBasis<T> &mbasis_curve,
                                        const T &R,
                                        gsSparseMatrix<T> &matrix,
                                        gsVector<T> &vector);

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    bool nonlinear = false;
    real_t R = 1.00; //beam-to-shell stiffness ratio
    index_t numRefine  = 0;
    std::string output;

    gsCmdLine cmd("Static analysis of a pinned X-stiffened shell.");
    cmd.addReal("R", "stiffnessRatio", "Stiffness ratio between shell and beam material", R);
    cmd.addInt( "r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addString("o", "output", "Name of the output file", output);
    cmd.addSwitch( "nl", "Enable nonlinear analysis", nonlinear );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Shell material properties]
    real_t E_modulus = 210e9; // [Pa]
    real_t PoissonRatio = 0.30;
    real_t Density = 7.85e3; // [kg/m^3]
    real_t thickness = 3.0e-3; // [m]
    //! [Shell material properties]

    //! [Shell geometry]
    gsMultiPatch<> mp_surf;
    real_t r = 10; // [m]
    mp_surf.addPatch(gsNurbsCreator<>::BSplineSquare(r));
    mp_surf.embed(3);
    mp_surf.addAutoBoundaries();
    mp_surf.degreeElevate(1,-1); //set degree to 2 in both directions

    for (int r =0; r < numRefine; ++r)
        mp_surf.uniformRefine(); //h-refine;
    
    gsGeometry<real_t> &surfgeo = mp_surf.patch(0);
    gsTensorBSpline<2, real_t>* surf = dynamic_cast< gsTensorBSpline<2, real_t>* >(&surfgeo);
    gsWriteParaview(mp_surf, "surf", 1000, true, false);

    gsMultiBasis<> mbasis_surf(mp_surf);
    gsInfo << "Patches: "<< mp_surf.nPatches() <<", degree: "<< mbasis_surf.minCwiseDegree() <<"\n";
    gsInfo << mbasis_surf.basis(0)<<"\n";

    gsMultiPatch<> mp_surf_def;
    mp_surf_def = mp_surf;
    //! [Shell geometry]

    //! [Embedded entities]
    gsKnotVector<real_t> kv_c = surf->knots(0);
    gsBSplineBasis<> basis_c(kv_c);

    gsEigen::ArrayXXd cpvec (surf-> knots(0).size()-surf-> degree(0)-1,1);
    cpvec = (surf->coefs().block(0,0,cpvec.rows(),1))/r;
    auto cpvec_flipped = cpvec.reverse();
    gsMatrix<real_t> coef_c1(basis_c.size(), surf->parDim()); //fill stiffener 1 control net
    coef_c1.col(0) = cpvec;
    coef_c1.col(1) = cpvec_flipped;
    gsMatrix<real_t> coef_c2(basis_c.size(), surf->parDim()); //fill stiffener 2 control net
    coef_c2.col(0) = cpvec;
    coef_c2.col(1) = cpvec;

    gsBSpline<> stiffener1(basis_c, coef_c1);
    gsBSpline<> stiffener2(basis_c, coef_c2);

    gsMultiPatch<> mp_curve;
    mp_curve.addPatch(stiffener1);
    mp_curve.addPatch(stiffener2);
    gsWriteParaview(mp_curve,  "curve",  1000, true,  true);
    gsMultiBasis<> mbasis_curve(mp_curve);
    //! [Embedded entities]

    //! [Set boundary conditions and loads]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp_surf);

    // Pressure under shell
    gsVector<> tmp(3);
    tmp << 0,0,-10; //[N/m^2]
    gsConstantFunction<> force(tmp,3);

    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, i);
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false,  i);
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false,  i);
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, i);
    }
    //! [Set boundary conditions and loads]

    //! [Make material functions]
    // Linear isotropic material model
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    //! [Make material functions]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    materialMatrix = getMaterialMatrix<3,real_t>(mp_surf,t,parameters,rho,options);

    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t,true>(mp_surf,mbasis_surf,bc,force,materialMatrix);
    //assembler->options().setInt("Continuity",0);
    //! [Make assembler]

    //! [Assemble linear part]
    ThinShellAssemblerStatus status;
    gsInfo<<"Setting up shell assembly\n";
    status = assembler->assemble();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    gsInfo<<"Shell assembly done\n";

    gsInfo <<"Setting up beam assembly\n";
    assembleLinearEmbeddedBeam_into<real_t>(assembler,mp_surf,mp_curve,mbasis_curve,R,matrix);
    gsInfo<<"Beam assembly done\n";
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler->numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    //! [Lambda functions for Jacobian and Residual definition]
    std::function<void(gsVector<real_t> const&, gsSparseMatrix<real_t> &, gsVector<real_t> &)> AssembleNonlinear = 
                 [&assembler,&mp_surf,&mp_curve,&mbasis_curve,&R](gsVector<real_t> const &x, gsSparseMatrix<real_t> &jacMat, gsVector<real_t> &rhsVec)
    {
        gsMultiPatch<> mp_surf_def;
        ThinShellAssemblerStatus status;

        assembler->constructSolution(x,mp_surf_def);

        status = assembler->assembleMatrix(mp_surf_def);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");
        status = assembler->assembleVector(mp_surf_def);
        GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");

        jacMat = assembler->matrix();
        rhsVec = assembler->rhs();
        assembleNonlinearEmbeddedBeam_into(assembler,mp_surf,mp_surf_def,mp_curve,mbasis_curve,R,jacMat,rhsVec);
    };
    //! [Lambda functions for Jacobian and Residual definition]

    //! [Solve nonlinear problem]
    if (nonlinear)
    {
        real_t residual = vector.norm();
        real_t residual0 = residual;
        real_t residualOld = residual;
        gsVector<real_t> updateVector, rhsVec;
        gsSparseMatrix<real_t> jacMat;
        for (index_t it = 0; it != 100; ++it)
        {
            AssembleNonlinear(solVector,jacMat,rhsVec);
            solver.compute(jacMat);
            updateVector = solver.solve(rhsVec);
            solVector += updateVector;
            residual = rhsVec.norm();

            gsInfo<<"Iteration: "<< it
               <<", residue: "<< residual
               <<", rel. residue: "<< residual/residual0
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

    }
    //! [Solve nonlinear problem]

    //! [Construct and evaluate solution]
    mp_surf_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_surf_def;
    for (size_t k = 0; k != mp_surf_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp_surf.patch(k).coefs();

    gsMatrix<> pts(2,2);
    pts.col(0)<<5,5;
    pts.col(1)<<7.5,5;

    gsMatrix<> defpts = deformation.patch(0).eval(pts);
    gsMatrix<> physpts = mp_surf.patch(0).eval(pts);
    gsInfo<<"Deformation at point "<<pts.col(0).transpose()<<": "<<defpts.col(0).transpose()<<"\n";
    gsInfo<<"Deformation at point "<<pts.col(1).transpose()<<": "<<defpts.col(1).transpose()<<"\n";

    if (!output.empty())
    {
        std::ofstream file;
        file.open(output);
        for (index_t k=0; k!=physpts.cols(); k++)
            file<<"point "<<k<<" - x,"<<"point "<<k<<" - y,"<<"point "<<k<<" - z,"
                <<"point "<<k<<" - dx,"<<"point "<<k<<" - dy,"<<"point "<<k<<" - dz,";
        file<<"\n";
        for (index_t k=0; k!=physpts.cols(); k++)
            file<<physpts(0,k)<<","<<physpts(1,k)<<","<<physpts(2,k)<<","
                <<defpts(0,k)<<","<<defpts(1,k)<<","<<defpts(2,k)<<",";
        file<<"\n";
        file.close();
    }

    gsInfo <<"Maximum deformation coef: "
            << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
            << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    //! [Construct and evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_surf_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, false);
        gsWriteParaview<>( mp_surf_def, "mp_surf_def", 1000, true,true);
    }

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;
}

template <class T>
void assembleLinearEmbeddedBeam_into(gsThinShellAssemblerBase<T> *assembler,
                                     const gsMultiPatch<T> &mp_surf,
                                     const gsMultiPatch<T> &mp_curve,
                                     const gsMultiBasis<T> &mbasis_curve,
                                     const T &R,
                                     gsSparseMatrix<T> &matrix)
{
    GISMO_ASSERT(assembler,"Assembler not defined");
    GISMO_ASSERT(matrix.rows()==assembler->numDofs() && matrix.cols()==assembler->numDofs(),"Matrix size mismatch");

    // Embedded entity mechanical properties
    real_t E_modulus_b = R * 210e9; // [Pa]
    real_t PoissonRatio_b = 0.30;
    real_t thickness_b = 3e-3; // [m]
    real_t height_b = 30e-3; // [m]
    real_t G_modulus_b = 0.5 * E_modulus_b / (1+PoissonRatio_b);
    real_t area_b = height_b * thickness_b;                   //cross-sectional area
    real_t inertiamin_b = (height_b * pow(thickness_b,3))/12; //minimum principal inertia moment of beam cross-section
    real_t inertiamax_b = (thickness_b * pow(height_b,3))/12; //maximum principal inertia moment of beam cross-section
    real_t inertiap_b = inertiamin_b + inertiamax_b;          //polar inertia moment of beam cross-section

    // Declare the expression assembler
    gsExprAssembler<> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();

    // Register expressions inside assembler
    auto G_surf  = exprAssembler.getMap(mp_surf);
    auto u_surf = exprAssembler.trialSpace(0);
    index_t N = exprAssembler.numDofs();

    for (size_t p = 0; p != mp_curve.nPatches(); ++p)
    {
        auto G_curve = exprAssembler.getMap(mp_curve.patch(p));
        gsSparseMatrix<T> matrix_curve(N,N);
 
        auto eps_der = ctv_var1(u_surf,G_surf,G_curve) * ctv(G_surf, G_curve);
        auto k21_der = cnv_vara_var1_normalized(u_surf,G_surf,G_curve)*ctv(G_surf, G_curve) + ctv_var1(u_surf,G_surf,G_curve)*cnv_vara_normalized(G_surf, G_curve);
        auto k31_der = cbv_vara_var1_normalized(u_surf,G_surf,G_curve)*ctv(G_surf, G_curve) + ctv_var1(u_surf,G_surf,G_curve)*cbv_vara_normalized(G_surf, G_curve);
        auto k23_der = cnv_vara_var1_normalized(u_surf,G_surf,G_curve)*cbv(G_surf, G_curve) + cbv_var1(u_surf,G_surf,G_curve)*cnv_vara_normalized(G_surf, G_curve);
        auto k32_der = cbv_vara_var1_normalized(u_surf,G_surf,G_curve)*cnv(G_surf, G_curve) + cnv_var1(u_surf,G_surf,G_curve)*cbv_vara_normalized(G_surf, G_curve);

        //Beam stiffness matrix (BendingAxial + Torsional)
        auto stiff = E_modulus_b/pow(ctv(G_surf, G_curve).norm(),3) *
                    (area_b * eps_der * eps_der.tr() + inertiamax_b * k21_der * k21_der.tr() + inertiamin_b * k31_der * k31_der.tr())
                                                                   +
                    G_modulus_b*inertiap_b/(4*pow(ctv(G_surf, G_curve).norm(),3)) *
                    (-k32_der * k32_der.tr() + k23_der * k23_der.tr());

        // Setup domain iterator
        const gsBasis<T> &basis = mbasis_curve.basis(p);
        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

        // Define quadrature rule
        typename gsQuadRule<T>::uPtr QuRule = gsQuadrature::getPtr(basis,assembler->options().getGroup("ExprAssembler"));

        // Loop over curve elements
        gsMatrix<T> quPointsCurve, quPointsSurface, quPointsPhysical;
        gsVector<T> quWeights;

        // Make expression data
        gsExprEvaluator<T> exprEvaluator(exprAssembler);

        index_t pt=0;
        for (; domIt->good(); domIt->next(), pt++ )
        {
            // Map the quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(), quPointsCurve, quWeights);

            // Map the quadrature points to the surface
            mp_curve.patch(p).eval_into(quPointsCurve, quPointsSurface);

            // Loop over quadrature points
            for (index_t k = 0; k != quWeights.rows(); ++k)
            {
                gsMatrix<T> evalMat = exprEvaluator.eval(stiff,quPointsCurve.col(k));
                gsMatrix<T> localMat = quWeights[k] * evalMat;

                // Push the local element matrix inside the big system
                const expr::gsFeSpace<real_t> & v = stiff.rowVar();
                const expr::gsFeSpace<real_t> & u = stiff.colVar();
                const index_t rd                  = v.dim();
                const index_t cd                  = u.dim();
                const gsDofMapper  & rowMap       = v.mapper();
                const gsDofMapper  & colMap       = u.mapper();
                gsMatrix<index_t>   colInd0       = assembler->getSpaceBasis().basis(0).active(quPointsSurface.col(k));
                gsMatrix<index_t>   rowInd0       = assembler->getSpaceBasis().basis(0).active(quPointsSurface.col(k));

                for (index_t r = 0; r != rd; ++r)
                {
                    const index_t rls = r * rowInd0.rows();     
                    for (index_t i = 0; i != rowInd0.rows(); ++i)
                    {
                        const index_t ii = rowMap.index(rowInd0.at(i),0,r); 
                        if ( rowMap.is_free_index(ii) )
                        {
                            for (index_t c = 0; c != cd; ++c)
                            {
                                const index_t cls = c * colInd0.rows();     
                                for (index_t j = 0; j != colInd0.rows(); ++j)
                                {
                                    if ( 0 == localMat(rls+i,cls+j) ) continue;
                                    const index_t jj = colMap.index(colInd0.at(j),0,c);
                                    if ( colMap.is_free_index(jj) )
                                    {
                                        matrix_curve.coeffRef(ii, jj) += localMat(rls+i,cls+j);
                                    }
                                }
                            }
                        }
                    }
                }

            }
        }
        matrix += matrix_curve; // TODO: Write everything to sparse entries and build the sparse matrix later.
    }
    return;
};

template <class T>
void assembleNonlinearEmbeddedBeam_into(gsThinShellAssemblerBase<T> *assembler,
                                        const gsMultiPatch<T> &mp_surf,
                                        const gsMultiPatch<T> &mp_surf_def,
                                        const gsMultiPatch<T> &mp_curve,
                                        const gsMultiBasis<T> &mbasis_curve,
                                        const T &R,
                                        gsSparseMatrix<T> &matrix,
                                        gsVector<T> &vector)
{
    GISMO_ASSERT(assembler,"Assembler not defined");
    GISMO_ASSERT(matrix.rows()==assembler->numDofs() && matrix.cols()==assembler->numDofs(),"Matrix size mismatch");
    GISMO_ASSERT(vector.size()==assembler->numDofs(),"Vector size mismatch");

    // Embedded entity mechanical properties
    real_t E_modulus_b = R * 210e9; // [Pa]
    real_t PoissonRatio_b = 0.30;
    real_t thickness_b = 3e-3; // [m]
    real_t height_b = 30e-3; // [m]
    real_t G_modulus_b = 0.5 * E_modulus_b / (1+PoissonRatio_b);
    real_t area_b = height_b * thickness_b;                   //cross-sectional area
    real_t inertiamin_b = (height_b * pow(thickness_b,3))/12; //minimum principal inertia moment of beam cross-section
    real_t inertiamax_b = (thickness_b * pow(height_b,3))/12; //maximum principal inertia moment of beam cross-section
    real_t inertiap_b = inertiamin_b + inertiamax_b;          //polar inertia moment of beam cross-section

    // Declare the expression assembler
    gsExprAssembler<> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();

    // Register expressions inside assembler
    auto G_surf  = exprAssembler.getMap(mp_surf);
    auto defG_surf = exprAssembler.getMap(mp_surf_def);
    auto u_surf = exprAssembler.trialSpace(0);
    index_t N = exprAssembler.numDofs();

    for (size_t p = 0; p != mp_curve.nPatches(); ++p)
    {
        auto G_curve = exprAssembler.getMap(mp_curve.patch(p));
        gsSparseMatrix<T> matrix_curve(N,N);
        gsVector<> vector_curve(N);
        vector_curve.setZero();

        auto eps = 0.5 * (ctv(defG_surf, G_curve).tr()*ctv(defG_surf, G_curve) - ctv(G_surf, G_curve).tr()*ctv(G_surf, G_curve));
        auto eps_der = ctv_var1(u_surf,defG_surf,G_curve) * ctv(defG_surf, G_curve);
        auto eps_der2 = ctv_var1(u_surf,defG_surf,G_curve) * ctv_var1(u_surf,defG_surf,G_curve).tr();
        
        auto k21 = cnv_vara_normalized(defG_surf, G_curve).tr()*ctv(defG_surf, G_curve) - cnv_vara_normalized(G_surf, G_curve).tr()*ctv(G_surf, G_curve);
        auto k21_der = cnv_vara_var1_normalized(u_surf,defG_surf,G_curve)*ctv(defG_surf, G_curve) + ctv_var1(u_surf,defG_surf,G_curve)*cnv_vara_normalized(defG_surf, G_curve);
        auto k21_der2 = cnv_vara_var2dot(u_surf,u_surf,defG_surf,G_curve,ctv(defG_surf, G_curve)) 
                        + var1_dot_var1(u_surf,u_surf,G_curve,cnv_vara_var1_normalized(u_surf,defG_surf,G_curve),ctv_var1(u_surf,defG_surf,G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,ctv_var1(u_surf,defG_surf,G_curve),cnv_vara_var1_normalized(u_surf,defG_surf,G_curve));

        auto k31 = cbv_vara_normalized(defG_surf, G_curve).tr()*ctv(defG_surf, G_curve) - cbv_vara_normalized(G_surf, G_curve).tr()*ctv(G_surf, G_curve);
        auto k31_der = cbv_vara_var1_normalized(u_surf,defG_surf,G_curve)*ctv(defG_surf, G_curve) + ctv_var1(u_surf,defG_surf,G_curve)*cbv_vara_normalized(defG_surf, G_curve);
        auto k31_der2 = cbv_vara_var2dot(u_surf,u_surf,defG_surf,G_curve,ctv(defG_surf, G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,cbv_vara_var1_normalized(u_surf,defG_surf,G_curve),ctv_var1(u_surf,defG_surf,G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,ctv_var1(u_surf,defG_surf,G_curve),cbv_vara_var1_normalized(u_surf,defG_surf,G_curve));
        
        auto k23 = cnv_vara_normalized(defG_surf, G_curve).tr()*cbv(defG_surf, G_curve) - cnv_vara_normalized(G_surf, G_curve).tr()*cbv(G_surf, G_curve);
        auto k23_der = cnv_vara_var1_normalized(u_surf,defG_surf,G_curve)*cbv(defG_surf, G_curve) + cbv_var1(u_surf,defG_surf,G_curve)*cnv_vara_normalized(defG_surf, G_curve);
        auto k23_der2 = cnv_vara_var2dot(u_surf,u_surf,defG_surf,G_curve,cbv(defG_surf, G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,cnv_vara_var1_normalized(u_surf,defG_surf,G_curve),cbv_var1(u_surf,defG_surf,G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,cbv_var1(u_surf,defG_surf,G_curve),cnv_vara_var1_normalized(u_surf,defG_surf,G_curve))
                        + cbv_var2dot(u_surf,u_surf,defG_surf,G_curve,cnv_vara_normalized(defG_surf,G_curve));

        auto k32 = cbv_vara_normalized(defG_surf, G_curve).tr()*cnv(defG_surf, G_curve) - cbv_vara_normalized(G_surf, G_curve).tr()*cnv(G_surf, G_curve);
        auto k32_der = cbv_vara_var1_normalized(u_surf,defG_surf,G_curve)*cnv(defG_surf, G_curve) + cnv_var1(u_surf,defG_surf,G_curve)*cbv_vara_normalized(defG_surf, G_curve);
        auto k32_der2 = cbv_vara_var2dot(u_surf,u_surf,defG_surf,G_curve,cnv(defG_surf,G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,cbv_vara_var1_normalized(u_surf,defG_surf,G_curve),cnv_var1(u_surf,defG_surf,G_curve))
                        + var1_dot_var1(u_surf,u_surf,G_curve,cnv_var1(u_surf,defG_surf,G_curve),cbv_vara_var1_normalized(u_surf,defG_surf,G_curve))
                        + cnv_var2dot(u_surf,u_surf,G_curve,cbv_vara_normalized(defG_surf,G_curve));

        //Beam stiffness matrix (BendingAxial + Torsional)
        auto stiff = E_modulus_b/pow(ctv(G_surf, G_curve).norm(),3) * 
                     (area_b * eps_der * eps_der.tr() + inertiamax_b * k21_der * k21_der.tr() + inertiamin_b * k31_der * k31_der.tr() +
                      area_b * eps.val() * eps_der2 + inertiamax_b * k21.val() * k21_der2 + inertiamin_b * k31.val() * k31_der2)
                                                                   +
                     G_modulus_b*inertiap_b/(4*pow(ctv(G_surf, G_curve).norm(),3)) *
                     (-k32_der * k32_der.tr() + k23_der * k23_der.tr() + k23.val() * k23_der2 + (-k32.val()) * k32_der2);

        //Beam residual vector (BendingAxial + Torsional) (rhs = -F_int) <<<<<<<< IS MINUS SIGN REALLY NEEDED? >>>>>>>
        auto rhs =  -((E_modulus_b/pow(ctv(G_surf, G_curve).norm(),3)) * 
                      (area_b * eps.val() * eps_der + inertiamax_b * k21.val() * k21_der + inertiamin_b * k31.val() * k31_der)
                                                                   +
                      (G_modulus_b*inertiap_b/(4*pow(ctv(G_surf, G_curve).norm(),3))) * (-k32.val() * k32_der + k23.val() * k23_der));

        // Setup domain iterator
        const gsBasis<T> &basis = mbasis_curve.basis(p);
        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

        // Define quadrature rule
        typename gsQuadRule<T>::uPtr QuRule = gsQuadrature::getPtr(basis,assembler->options().getGroup("ExprAssembler"));

        // Loop over curve elements
        gsMatrix<T> quPointsCurve, quPointsSurface, quPointsPhysical;
        gsVector<T> quWeights;

        // Make expression data
        gsExprEvaluator<T> exprEvaluator(exprAssembler);

        index_t pt=0;
        for (; domIt->good(); domIt->next(), pt++ )
        {
            // Map the quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(), quPointsCurve, quWeights);

            // Map the quadrature points to the surface
            mp_curve.patch(p).eval_into(quPointsCurve, quPointsSurface);

            // Loop over quadrature points
            for (index_t k = 0; k != quWeights.rows(); ++k)
            {
                gsMatrix<T> evalMat = exprEvaluator.eval(stiff,quPointsCurve.col(k));
                gsMatrix<T> localMat = quWeights[k] * evalMat;
                gsMatrix<> evalRhs = exprEvaluator.eval(rhs,quPointsCurve.col(k));
                gsMatrix<> localRhs = quWeights[k] * evalRhs;

                // Push the local element matrix inside the big system
                const expr::gsFeSpace<real_t> & v = stiff.rowVar();
                const expr::gsFeSpace<real_t> & u = stiff.colVar();
                const index_t rd                  = v.dim();
                const index_t cd                  = u.dim();
                const gsDofMapper  & rowMap       = v.mapper();
                const gsDofMapper  & colMap       = u.mapper();
                gsMatrix<index_t>   colInd0       = assembler->getSpaceBasis().basis(0).active(quPointsSurface.col(k));
                gsMatrix<index_t>   rowInd0       = assembler->getSpaceBasis().basis(0).active(quPointsSurface.col(k));

                for (index_t r = 0; r != rd; ++r)
                {
                    const index_t rls = r * rowInd0.rows();     //local stride
                    for (index_t i = 0; i != rowInd0.rows(); ++i)
                    {
                        // const index_t ii = rowMap.index(rowInd0.at(i),v.data().patchId,r); //N_i
                        const index_t ii = rowMap.index(rowInd0.at(i),0,r); //N_i
                        if ( rowMap.is_free_index(ii) )
                        {
                            for (index_t c = 0; c != cd; ++c)
                            {
                                const index_t cls = c * colInd0.rows();     //local stride
                                for (index_t j = 0; j != colInd0.rows(); ++j)
                                {
                                    if ( 0 == localMat(rls+i,cls+j) ) continue;

                                    // const index_t jj = colMap.index(colInd0.at(j),u.data().patchId,c); // N_j
                                    const index_t jj = colMap.index(colInd0.at(j),0,c); // N_j
                                    if ( colMap.is_free_index(jj) )
                                    {
                                        matrix_curve.coeffRef(ii, jj) += localMat(rls+i,cls+j);
                                    }
                                }
                            }
                            vector_curve[ii] += localRhs(rls+i,0);
                        }
                    }
                }
            }
        }
        matrix += matrix_curve; // TODO: Write everything to sparse entries and build the sparse matrix later
        vector += vector_curve;
    }
    return;
};