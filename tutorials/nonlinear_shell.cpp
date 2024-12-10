/** @file gsKLShell/tutorials/nonlinear_shell.cpp

    @brief Tutorial for assembling the Non-Linear Kirchhoff-Love shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

//! [Includes]
#include <gsKLShell/gsKLShell.h>
//! [Includes]

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    gsCmdLine cmd("Nonlinear shell tutorial.");
    cmd.addInt( "e", "degreeElevation","Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    // Initialize [ori]ginal and [def]ormed geometry
    gsMultiPatch<> ori;
    gsReadFile<>("surfaces/paraboloid.xml", ori);
    ori.computeTopology();
    //! [Read geometry]

    //! [Initialize geometry]
    // p-refine
    if (numElevate!=0)
        ori.degreeElevate(numElevate);
    // h-refine
    for (int r =0; r < numRefine; ++r)
        ori.uniformRefine();
    //! [Initialize geometry]

    //! [Construct basis]
    gsMultiBasis<> bases(ori);
    //! [Construct basis]

    //! [Set boundary conditions]
    // Define the boundary conditions object
    gsBoundaryConditions<> bc;
    // Set the geometry map for computation of the Dirichlet BCs
    bc.setGeoMap(ori);

    // Set the boundary conditions
    bc.addCornerValue(boundary::southwest, 0.0, 0, -1); // (corner,value, patch, unknown)
    bc.addCornerValue(boundary::southeast, 0.0, 0, -1); // (corner,value, patch, unknown)
    bc.addCornerValue(boundary::northwest, 0.0, 0, -1); // (corner,value, patch, unknown)
    bc.addCornerValue(boundary::northeast, 0.0, 0, -1); // (corner,value, patch, unknown)
    //! [Set boundary conditions]

    //! [Loads]
    // Point loads are defined using a point in the parametric domain, and a 3D load.
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsVector<> point(2);
    gsVector<> load (3);
    point<< 0.5, 0.5 ; load << 0, 0.0, -1e4 ;
    pLoads.addLoad(point, load, 0 );

    // The surface force is zero, and defined in the physical space, i.e. 3D
    gsFunctionExpr<> force("0","0","0",3);
    //! [Loads]


    //! [Define the material matrix class]
    // Define material parameters
    // The material parameters are defined in the physical domain as well (!)
    gsConstantFunction<> E (1.E9,3);
    gsConstantFunction<> nu(0.45,3);
    gsConstantFunction<> t (1E-2,3);

    // Define a linear material, see \ref gsMaterialMatrixLinear.h
    // The first parameter is the physical domain dimension

    // gsMaterialMatrixLinear<3,real_t> materialMatrix(ori,t,E,nu);

    std::vector<gsFunctionSet<>*> parameters{&E,&nu};
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix = getMaterialMatrix<3,real_t>(ori,t,parameters,options);
    //! [Define the material matrix class]

    //! [Define the assembler]
    // Define the assembler.
    // The first parameter is the physical domain dimension
    // The second parameter is the real type
    // The third parameter controls the bending stiffness
    gsThinShellAssembler<3, real_t, true > assembler(ori,bases,bc,force,materialMatrix);
    assembler.setPointLoads(pLoads);
    //! [Define the assembler]


    //! [Define nonlinear residual functions]
    typedef std::function<bool (gsVector<real_t> const &, gsSparseMatrix<real_t>&)>    Jacobian_t;
    typedef std::function<bool (gsVector<real_t> const &, gsVector<real_t>      &) >   Residual_t;

    Jacobian_t Jacobian = [&assembler](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        gsMultiPatch<> def;
        assembler.constructSolution(x,def);
        ThinShellAssemblerStatus status = assembler.assembleMatrix(def);
        m = assembler.matrix();
        return status==ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    Residual_t Residual = [&assembler](gsVector<real_t> const &x, gsVector<real_t> & v)
    {
        gsMultiPatch<> def;
        assembler.constructSolution(x,def);
        ThinShellAssemblerStatus status = assembler.assembleVector(def);
        v = assembler.rhs();
        return status==ThinShellAssemblerStatus::Success;
    };
    //! [Define nonlinear residual functions]


    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> F = assembler.rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(F);
    //! [Solve linear problem]

    //! [Solve nonlinear problem]
    gsVector<> updateVector = solVector;
    gsVector<> R;
    GISMO_ENSURE(Residual(solVector,R),"Assembly of the residual failed");
    for (index_t it = 0; it != 100; ++it)
    {
        GISMO_ENSURE(Jacobian(solVector,matrix),"Assembly of the Jacobian failed");
        solver.compute(matrix);
        updateVector = solver.solve(R); // this is the UPDATE
        solVector += updateVector;

        GISMO_ENSURE(Residual(solVector,R),"Assembly of the residual failed");

        gsInfo<<"Iteration: "<< it
           <<", residue: "<< R.norm()/F.norm()
           <<", update norm: "<<updateVector.norm()/solVector.norm()
           <<"\n";

        if (updateVector.norm() < 1e-6)
            break;
        else if (it+1 == it)
            gsWarn<<"Maximum iterations reached!\n";
    }
    //! [Solve nonlinear problem]

    //! [Construct solution and deformed geometry]
    gsMultiPatch<> displ = assembler.constructDisplacement(solVector);
    gsMultiPatch<> def = assembler.constructSolution(solVector);
    //! [Construct solution and deformed geometry]


    //! [Evaluate solution]
    // Evaluate the solution on a reference point (parametric coordinates)
    // The reference points are stored column-wise
    gsMatrix<> refPoint(2,1);
    refPoint<<0.5,0.5;
    // Compute the values
    gsVector<> physpoint = def.patch(0).eval(refPoint);
    gsVector<> refVals = displ.patch(0).eval(refPoint);
    // gsInfo << "Displacement at reference point: "<<numVal<<"\n";
    gsInfo  << "Displacement at reference point ("
            <<ori.patch(0).eval(refPoint).transpose()<<"): "
            <<": "<<refVals.transpose()<<"\n";
    //! [Evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        // Plot the displacements on the deformed geometry
        gsField<> solField(def, displ);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);

        // Plot the membrane Von-Mises stresses on the geometry
        gsPiecewiseFunction<> VMm;
        assembler.constructStress(def,VMm,stress_type::von_mises_membrane);
        gsField<> membraneStress(def,VMm, true);
        gsWriteParaview(membraneStress,"MembraneStress",1000);
    }
    // ! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main