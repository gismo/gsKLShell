/** @file gsKLShell/tutorials/linear_shell.cpp

    @brief Tutorial for assembling the Kirchhoff-Love shell

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

    gsCmdLine cmd("Linear shell tutorial.");
    cmd.addInt( "e", "degreeElevation","Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    // Initialize [ori]ginal and [def]ormed geometry
    gsMultiPatch<> ori;
    std::string fileName = "surfaces/scordelis_lo_roof.xml";
    gsReadFile<>(fileName, ori);
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
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 1 );
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 2 );

    bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 1 );
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 2 );
    //! [Set boundary conditions]

    //! [Set surface force]
    // The surface force is defined in the physical space, i.e. 3D
    gsFunctionExpr<> force("0","0","-90",3);
    //! [Set surface force]


    //! [Define the material matrix class]
    // Define material parameters
    // The material parameters are defined in the physical domain as well (!)
    gsConstantFunction<> E (4.32E8,3);
    gsConstantFunction<> nu(0.0   ,3);
    gsConstantFunction<> t (0.25  ,3);

    // Define a linear material, see \ref gsMaterialMatrixLinear.h
    // The first parameter is the physical domain dimension

    // gsMaterialMatrixLinear<3,real_t> materialMatrix(ori,t,E,nu);

    std::vector<gsFunctionSet<>*> parameters{&E,&nu};
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix = getMaterialMatrix<3,real_t>(ori,t,parameters,options);
    //! [Define the material matrix class]

    //! [Define the assembler]
    // Define the assembler.
    // The first parameter is the physical domain dimension
    // The second parameter is the real type
    // The third parameter controls the bending stiffness
    gsThinShellAssembler<3, real_t, true > assembler(ori,bases,bc,force,materialMatrix);
    //! [Define the assembler]

    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    //! [Construct solution and deformed geometry]
    gsMultiPatch<> displ = assembler.constructDisplacement(solVector);
    gsMultiPatch<> def = assembler.constructSolution(solVector);
    //! [Construct solution and deformed geometry]


    //! [Evaluate solution]
    // Evaluate the solution on a reference point (parametric coordinates)
    // The reference points are stored column-wise
    gsMatrix<> refPoint(2,1);
    refPoint<<0.5,1;
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