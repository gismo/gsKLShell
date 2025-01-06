![GitHub commits since latest release](https://img.shields.io/github/commits-since/gismo/gsKLShell/latest?color=008A00)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/gismo/gsKLShell?color=008A00)

# gsKLShell

Module for the isogeometric Kirchhoff-Love shell element. The module is based on `gismo`'s Expression Assembler `gsExprAssembler`.

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsKLShell"```|
|--:|---|
|License|![GitHub License](https://img.shields.io/github/license/gismo/gismo?color=008A00)|
|OS support|Linux, Windows, macOS|
|Build status|[![ci](https://github.com/gismo/gsKLShell/actions/workflows/ci.yml/badge.svg)](https://github.com/gismo/gsKLShell/actions/workflows/ci.yml)|
|Repository|[gismo/gismo](https://github.com/gismo/gismo)|
|Developers/maintainers| [![Static Badge](https://img.shields.io/badge/@hverhelst-008A00)](https://github.com/hverhelst) [![Static Badge](https://img.shields.io/badge/@Crazy--Rich--Meghan-008A00)](https://github.com/Crazy-Rich-Meghan)|

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

`gsThinShellAssemblerDWR`
* Error estimation via the Dual-Weighted Residual method

`gsMaterialMatrixBase` and derivatives
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
* `gsMaterialMatrixLinear`: class that handles the linear constitutive relations of the shells.
* `gsMaterialMatrixNonlinear`: class that handles the non-linear constitutive relations of the shells.
* `gsMaterialMatrixComposite`: class that handles constitutive relations for linear composites.
* `gsMaterialMatrixTFT`: class that handles tension-field-theory-based constitutive models for membranes.

See the doxygen manuals for more information about the classes. (**to do: add link**)

To use the `gsKLShell` module, one should always define a `gsMaterialMatrixBase` and a `gsThinShellAssembler`. The `gsMaterialMatrixBase` object is used to compute the constitutive relations in the `gsThinShellAssembler` and should therefore be defined upon initialisation of this class.

Additionally, the geometry and the deformed geometry (both `gsMultiPatch`) together with a basis (`gsMultiBasis`) and the boundary conditions (`gsBoundaryConditions`) and a surface force (`gsFunctionExpression`) should be provided in the definition of the class.

The template parameters of the class are the dimension of the geometry (`dim`) which is 2D (planar) or 3D (surface) and a flag for the computation of bending stiffness term (`bending`) which is only relevant if `dim==3`. Other options that can be set are follower pressures (`setPressure`), elastic foundation stiffness (`setFoundation`) and point loads (`setPointLoads`).
