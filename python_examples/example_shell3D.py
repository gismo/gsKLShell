#!/usr/bin/python

""""
    @file BSpline curve example

    @brief Play with a B-spline curve in Python

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np
import scipy.sparse.linalg as la

## See gismo/filedata/surfaces/simple.xml for the geometry
c1 = np.array([0.,0.,1.,1.])
c2 = np.array([0.,0.,1.,1.])
ku1 = gs.nurbs.gsKnotVector(c1,1)
ku2 = gs.nurbs.gsKnotVector(c2,1)

coefs = np.array([
                    [0     ,0    ,0   ],
                    [1     ,0    ,0   ],
                    [0     ,1    ,0   ],
                    [1     ,1    ,0   ],
                        ])


# Construct basis using knot vectors
tbasis1 = gs.nurbs.gsTensorBSplineBasis2(ku1,ku2)
tspline1 = gs.nurbs.gsTensorBSpline2(tbasis1,coefs)

print("Coefficients:\n", tspline1.coefs())

mp = gs.core.gsMultiPatch()
mp.addPatch(tspline1)

mp.degreeElevate()
mp.uniformRefine()

mb = gs.core.gsMultiBasis(mp)

t = gs.core.gsFunctionExpr("1",3)
E = gs.core.gsFunctionExpr("1",3)
nu = gs.core.gsFunctionExpr("0.0",3)
f = gs.core.gsFunctionExpr("0","0","-1",3)

null = gs.core.gsFunctionExpr("0",3)
side = gs.core.side.west
bcs = gs.pde.gsBoundaryConditions();

for d in range(0,3):
    for side in [gs.core.side.west, gs.core.side.east, gs.core.side.south, gs.core.side.north]:
        bcs.addCondition(0,side,gs.pde.bctype.dirichlet,null,0,False,d)

bcs.setGeoMap(mp);

mm = gs.klshell.gsMaterialMatrixLinear3(mp,t)
mm.setYoungsModulus(E)
mm.setPoissonsRatio(nu)

assembler = gs.klshell.gsThinShellAssembler3(mp,mb,bcs,f,mm)
assembler.assemble()
matrix = assembler.matrix()
vector = assembler.rhs()
print(assembler.matrix());
print(assembler.rhs());
solution = la.spsolve(matrix,vector[:,0])
sol = assembler.constructSolution(solution)
assembler.assembleVector(solution)
assembler.assembleMatrix(solution)
print(assembler.matrix());
print(assembler.rhs());

print("Solution coefficients:\n",sol.patch(0).coefs())

