# To Do

# Changelog

## v21.12:
  - Compatibility with the new `gsExprAssembler` in `gismo/gsAssembler`
## v21.6:
  - Initial version!

# gsKLShell

Module for the isogeometric Kirchhoff-Love shell element. The module is based on `gismo`'s Expression Assembler `gsExprAssembler`.

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsKLShell"```|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Build status| [CDash](link) |
|Repository|[gismo/gismo](https://github.com/gismo/gismo)|
|Status|completed|
|Developer|[Hugo Verhelst](https://github.com/hverhelst)|
|Maintainer|[h.m.verhelst@tudelft.nl](mailto:h.m.verhelst@tudelft.nl)|

#### Dependencies
No dependencies

#### Installation
```
cd path/to/build/dir
cmake . -DGISMO_OPTIONAL="<other submodules>;gsKLShell"
make
```

***

#### Overview of the `gsKLShell` module
`gsThinShellAssembler`
* Linear and Non-Linear kinematics
* Follower pressures and elastic foundation stiffness
* Supports B-spline, NURBS, H-Spline and THB-Spline bases
* Membrane or shell elements via template parameters

`gsMaterialMatrixNonlinear`
* Linear materials via Saint-Venant Kirchhoff model
* (In)Compressible non-linear materials: Neo-Hookean, Mooney-Rivlin and Ogden materials.
* Direct implementation for other material models possible
* Generalized formulations given the derivatives of the Strain Energy Density Function w.r.t components of the deformation tensor possible
* Stretch-based implementations given the derivatives of the Strain Energy Density Function w.r.t. the stretches
* Material and compressibility flags via template parameters


#### Use of the `gsKLShell` module
The `gsKLShell` module consists of the classes
* `gsThinShellAssembler`: class that resolves the kinematics of the Kirchhoff-Love shells
* `gsThinShellAssemblerDWR`: same as the above, but contains extra functions for the Dual Weighted Residual (DWR) method
* `gsMaterialMatrixNonlinear`: class that handles the (non-)linear constitutive relations of the shells.
* `gsMaterialMatrixLaminate`: class that handles constitutive relations for linear composites.

See the doxygen manuals for more information about the classes. (**to do: add link**)

To use the `gsKLShell` module, one should always define a `gsMaterialMatrixNonlinear` and a `gsThinShellAssembler`. The `gsMaterialMatrixNonlinear` object is used to compute the constitutive relations in the `gsThinShellAssembler` and should therefore be defined upon initialisation of this class.

Additionally, the geometry and the deformed geometry (both `gsMultiPatch`) together with a basis (`gsMultiBasis`) and the boundary conditions (`gsBoundaryConditions`) and a surface force (`gsFunctionExpression`) should be provided in the definition of the class.

The template parameters of the class are the dimension of the geometry (`dim`) which is 2D (planar) or 3D (surface) and a flag for the computation of bending stiffness term (`bending`) which is only relevant if `dim==3`. Other options that can be set are follower pressures (`setPressure`), elastic foundation stiffness (`setFoundation`) and point loads (`setPointLoads`).

#### Not (yet) supported:
* Multipatch coupling
* Error estimation via DWR
