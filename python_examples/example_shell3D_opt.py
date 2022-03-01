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
from mpl_toolkits.mplot3d import Axes3D

import copy
import scipy.optimize as opt

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


mp = gs.core.gsMultiPatch()
mp.addPatch(tspline1)

mp.degreeElevate()
# mp.degreeElevate()
mp.uniformRefine()
mp.uniformRefine()

coefs = mp.patch(0).coefs()

# np.where(mp.patch(0).coefs()[:,0].any() < 1e-12 and mp.patch(0).coefs()[:,1].any() < 1e-12)

mask_00     = (coefs[:,0] < 1e-12  ) & (coefs[:,1] < 1e-12  )
mask_01     = (coefs[:,0] < 1e-12  ) & (coefs[:,1] > 1-1e-12)
mask_10     = (coefs[:,0] > 1-1e-12) & (coefs[:,1] < 1e-12  )
mask_11     = (coefs[:,0] > 1-1e-12) & (coefs[:,1] > 1-1e-12)
cornermask  = mask_00 | mask_01 | mask_10 | mask_11
interiormask= ~cornermask

# print("Coefficients:\n", mp.patch(0).coefs())

mb = gs.core.gsMultiBasis(mp)

t = gs.core.gsFunctionExpr("1",3)
E = gs.core.gsFunctionExpr("1",3)
nu = gs.core.gsFunctionExpr("0.3",3)
f = gs.core.gsFunctionExpr("0","0","0",3)
pload = gs.pde.gsPointLoads()
pload.addLoad(np.array([0.5,0.5]),np.array([0,0,-10]),0,True)
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


low = coefs[interiormask]
upp = coefs[interiormask]
L = np.max(low[:,0]) - np.min(low[:,0])
W = np.max(upp[:,1]) - np.min(upp[:,1])

low[:,0] = low[:,0] - 0.5*L
low[:,1] = low[:,1] - 0.5*W
low[:,2] = low[:,2] - 0.5*L

upp[:,0] = upp[:,0] + 0.5*L
upp[:,1] = upp[:,1] + 0.5*W
upp[:,2] = upp[:,2] + 0.5*L

low = low.flatten()
upp = upp.flatten()

u = coefs[interiormask]
u[:,2] = 0.01*L*np.sin(u[:,0] * (np.pi))*np.sin(u[:,1] * (np.pi))
shape = np.shape(u)
u = u.flatten()

assembler = gs.klshell.gsThinShellAssembler3(mp,mb,bcs,f,mm)
assembler.setPointLoads(pload)

def makeGeometry(design):
    design = np.resize(design, shape)
    mp_tmp = gs.core.gsMultiPatch(mp)
    tmp_coefs = mp_tmp.patch(0).coefs()
    tmp_coefs[interiormask] = design
    mp_tmp.patch(0).setCoefs(tmp_coefs)
    return mp_tmp

def constructDisplacement(solution):
    return assembler.constructDisplacement(solution)

def constructSolution(solution):
    return assembler.constructSolution(solution)

def computeDeformation(mp_tmp):
    assembler.setGeometry(mp_tmp)
    assembler.assemble()
    matrix = assembler.matrix()
    vector = assembler.rhs()
    solution = la.spsolve(matrix,vector[:,0])
    # def Residual(resvec):
    #     sol = assembler.constructSolution(resvec)
    #     assembler.assembleVector(sol)
    #     return assembler.rhs()

    # def Jacobian(resvec):
    #     sol = assembler.constructSolution(resvec)
    #     assembler.assembleMatrix(sol)
    #     return assembler.matrix()
    #
    # residual = np.linalg.norm(vector)
    # residual0 = residual
    # residualOld = residual
    # update = solution
    # resvec = Residual(solution)

    # itmax = 100
    # tol = 1e-6
    # for it in range(0,itmax):
    #     jacmat = Jacobian(solution)
    #     update = la.spsolve(jacmat,resvec[:,0])
    #     solution += update

    #     resvec = Residual(solution)
    #     residual = np.linalg.norm(resvec)

    #     print("Iteration ",it,end="")
    #     print(", residue %0.5e" %residual,end="")
    #     print(", update norm %0.5e" %np.linalg.norm(update),end="")
    #     print(", log(Ri/R0) %0.5e" %np.log(residualOld/residual0),end="")
    #     print(", log(Ri+1/R0) %0.5e" %np.log(residual/residual0),end="")
    #     print("")

    #     residualOld = residual

    #     if (np.linalg.norm(update) < tol):
    #         break
    #     elif (it+1==itmax):
    #         print("Maximum iterations reached")

    return solution

def computeObjective(design):
    mp_tmp = makeGeometry(design)
    solution = computeDeformation(mp_tmp)
    sol = constructDisplacement(solution)

    nx = ny = 100
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    xv, yv = np.meshgrid(x,y,indexing='xy')
    pts = np.stack((xv.flatten(),yv.flatten()))

    deformation = -sol.patch(0).eval(pts)
    return np.max(deformation[2,:])

def computeConstraint(design):
    design = np.resize(design, shape)
    mp_tmp = gs.core.gsMultiPatch(mp)
    tmp_coefs = mp_tmp.patch(0).coefs()
    tmp_coefs[interiormask] = design
    mp_tmp.patch(0).setCoefs(tmp_coefs)
    return assembler.getArea(mp_tmp) - assembler.getArea(mp)

def plotGeometry(design,ax):
    mp_tmp = makeGeometry(design)
    nx = ny = 100
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x,y,indexing='xy')
    pts = np.stack((xv.flatten(),yv.flatten()))
    geom = mp_tmp.patch(0).eval(pts)
    x = geom[0,:].reshape(nx,ny)
    y = geom[1,:].reshape(nx,ny)
    z = geom[2,:].reshape(nx,ny)
    ax.plot_surface(x,y,z)
    return

def plotDeformation(design,ax):
    mp_tmp = makeGeometry(design)
    nx = ny = 100
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x,y,indexing='xy')
    pts = np.stack((xv.flatten(),yv.flatten()))

    solution = computeDeformation(mp_tmp)
    sol = constructSolution(solution)

    deformed = sol.patch(0).eval(pts)
    XX = deformed[0,:].reshape((nx,ny))
    YY = deformed[1,:].reshape((nx,ny))
    ZZ = deformed[2,:].reshape((nx,ny))
    ax.plot_surface(XX,YY,ZZ,cmap=cm.coolwarm)
    return

nlc = opt.NonlinearConstraint(computeConstraint, 0, 0)
bnd = opt.Bounds(low,upp,keep_feasible=False)

sol = opt.minimize(computeObjective, u,
    method = 'trust-constr',
    bounds=bnd,
    constraints=nlc,
    options={   'verbose':3,
                'maxiter':1000,
                # 'xtol':1e-7,
                # 'gtol':1e-5,
                # 'barrier_tol':1e-5,
                })


# PLOTTING -------------------------------------------------
fig = plt.figure(figsize =(14, 9))
ax11 = plt.subplot(221,projection ='3d')
ax12 = plt.subplot(222,projection ='3d')
ax21 = plt.subplot(223,projection ='3d')
ax22 = plt.subplot(224,projection ='3d')

plotGeometry(u,ax11)
plotDeformation(u,ax12)

plotGeometry(sol.x,ax21)
plotDeformation(sol.x,ax22)

ax11.set_title('Initial geometry')
ax12.set_title('Initial deformation')
ax21.set_title('Final geometry')
ax22.set_title('Final deformation')

plt.show()