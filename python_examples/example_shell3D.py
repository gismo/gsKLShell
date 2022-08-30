#!/usr/bin/python

""""
    @file example_shell3D.py

    @brief Replicates example_shell3D from the gsKLShell module

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../../../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np
import scipy.sparse.linalg as la
from matplotlib import cm
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D



mp = gs.core.gsMultiPatch()
fd = gs.io.gsFileData(os.path.join(os.getcwd() , "../../../filedata/3dm/BB2_clean.3dm"))
# fd.getId(0,mp)
fd.dump()


exit(0);

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
# mp.uniformRefine()

mb = gs.core.gsMultiBasis(mp)

t = gs.core.gsFunctionExpr("1",3)
E = gs.core.gsFunctionExpr("1",3)
nu = gs.core.gsFunctionExpr("0.3",3)
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
print(solution);
sol = assembler.constructSolution(solution)
# assembler.assembleVector(solution)
# assembler.assembleMatrix(solution)
# print(assembler.matrix());
# print(assembler.rhs());

print("Solution coefficients:\n",sol.patch(0).coefs())
u = np.array([0.5,0.5])
print(f"Linear solution on {u[0],u[1]}:\n", sol.patch(0).eval(u))


def Residual(resvec):
    sol = assembler.constructSolution(resvec)
    assembler.assembleVector(sol)
    return assembler.rhs()

def Jacobian(resvec):
    sol = assembler.constructSolution(resvec)
    assembler.assembleMatrix(sol)
    return assembler.matrix()



print("Nonlinear solve...")
residual = np.linalg.norm(vector)
residual0 = residual
residualOld = residual
update = solution
resvec = Residual(solution)

itmax = 100
tol = 1e-6
for it in range(0,itmax):
    jacmat = Jacobian(solution)
    update = la.spsolve(jacmat,resvec[:,0])
    solution += update

    resvec = Residual(solution)
    residual = np.linalg.norm(resvec)

    print("Iteration ",it,end="")
    print(", residue %0.5e" %residual,end="")
    print(", update norm %0.5e" %np.linalg.norm(update),end="")
    print(", log(Ri/R0) %0.5e" %np.log(residualOld/residual0),end="")
    print(", log(Ri+1/R0) %0.5e" %np.log(residual/residual0),end="")
    print("")

    residualOld = residual

    if (np.linalg.norm(update) < tol):
        break
    elif (it+1==itmax):
        print("Maximum iterations reached")

N = 20
X = Y = np.linspace(0,1,N)
XX,YY = np.meshgrid(X,Y)

XX = XX.reshape(N*N)
YY = YY.reshape(N*N)

u = np.stack([XX,YY])

sol = assembler.constructSolution(solution)

solution = sol.patch(0).eval(u)

fig = plt.figure()
ax = fig.add_subplot() # projection='3d'

XX = solution[0,:].reshape((N,N))
YY = solution[1,:].reshape((N,N))
ZZ = solution[2,:].reshape((N,N))
p = ax.contourf(XX,YY,ZZ,cmap=cm.coolwarm)
fig.colorbar(p)
plt.show()