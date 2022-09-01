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

# Required for finding pygismo
import os, sys
# Obtain pygismo
gismo_path=os.path.join(os.getcwd() , "../../../")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path+"build/lib")

# Import pygismo
import pygismo as gs ## If this line gives an error, check your path or check if pygismo is compiled

# Import other modules
import numpy as np
import scipy.sparse.linalg as la
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# for creating a responsive plot
# %matplotlib widget

import scipy.optimize as opt

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

# elevate the degree (starts at 1)
mp.degreeElevate()     #degree = 2
# mp.degreeElevate()     #degree = 3
mp.uniformRefine()     # 4 elements
# mp.uniformRefine()     # 16 elements

# multibasis
mb = gs.core.gsMultiBasis(mp)

# store the coefficients
coefs = mp.patch(0).coefs()

mask_00     = (coefs[:,0] < 1e-12  ) & (coefs[:,1] < 1e-12  )
mask_01     = (coefs[:,0] < 1e-12  ) & (coefs[:,1] > 1-1e-12)
mask_10     = (coefs[:,0] > 1-1e-12) & (coefs[:,1] < 1e-12  )
mask_11     = (coefs[:,0] > 1-1e-12) & (coefs[:,1] > 1-1e-12)
cornermask  = mask_00 | mask_01 | mask_10 | mask_11
interiormask= ~cornermask

shape = np.shape(coefs[interiormask])

# thickness
t = gs.core.gsFunctionExpr("1",3)
# Young's modulus
E = gs.core.gsFunctionExpr("1",3)
# Poisson's ratio
nu = gs.core.gsFunctionExpr("0.3",3)
# Surface force
f = gs.core.gsFunctionExpr("0","0","0",3)
# Point loads
pload = gs.pde.gsPointLoads()
pload.addLoad(np.array([0.5,0.5]),np.array([0,0,-1]),0,True)


bcs = gs.pde.gsBoundaryConditions();
for corner in [gs.core.corner.northwest,
             gs.core.corner.northeast,
             gs.core.corner.southwest,
             gs.core.corner.southeast]:
    bcs.addCornerValue(corner,0.0,0,0,-1) # corner, value, patch, unknown, component

# null = gs.core.gsFunctionExpr("0",3)
# side = gs.core.side.west
# for side in [gs.core.side.west, gs.core.side.east, gs.core.side.south, gs.core.side.north]:
#     bcs.addCondition(0,side,gs.pde.bctype.dirichlet,null,0,False,-1)

# assign the geometry to the boundary conditions
bcs.setGeoMap(mp)

mm = gs.klshell.gsMaterialMatrixLinear3(mp,t)
mm.setYoungsModulus(E)
mm.setPoissonsRatio(nu)

assembler = gs.klshell.gsThinShellAssembler3(mp,mb,bcs,f,mm)
assembler.setPointLoads(pload)

# Makes a deformed geometry given a vector of design variables
def makeGeometry(design):
    design = np.resize(design, shape)
    mp_tmp = gs.core.gsMultiPatch(mp)
    tmp_coefs = mp_tmp.patch(0).coefs()
    tmp_coefs[interiormask] = design
    mp_tmp.patch(0).setCoefs(tmp_coefs)
    return mp_tmp


# Computes the deformation given a deformed multipatch
# (NB: the commented part is a nonlinear solver)
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

# Constructs the displacement
def constructDisplacement(solution):
    return assembler.constructDisplacement(solution)

# Constructs the displaced shell
def constructSolution(solution):
    return assembler.constructSolution(solution)

# Computes the optimization objective
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

# Computes the area constraint
def computeConstraint(design):
    design = np.resize(design, shape)
    mp_tmp = gs.core.gsMultiPatch(mp)
    tmp_coefs = mp_tmp.patch(0).coefs()
    tmp_coefs[interiormask] = design
    mp_tmp.patch(0).setCoefs(tmp_coefs)
    return assembler.getArea(mp_tmp) - assembler.getArea(mp)

# Plots the geometry
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

# Plots the deformed geometry
def plotDeformation(design,ax):
    mp_tmp = makeGeometry(design)
    nx = ny = 100
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x,y,indexing='xy')
    pts = np.stack((xv.flatten(),yv.flatten()))

    sol = constructSolution(design)

    deformed = sol.patch(0).eval(pts)
    XX = deformed[0,:].reshape((nx,ny))
    YY = deformed[1,:].reshape((nx,ny))
    ZZ = deformed[2,:].reshape((nx,ny))
    ax.plot_surface(XX,YY,ZZ,cmap=cm.coolwarm)
    return

## Static solve
# Make vector of design variables
u = coefs[interiormask]
u = u.flatten()

mp_tmp = makeGeometry(u)
soln = computeDeformation(mp_tmp)
fig = plt.figure(figsize =(14, 9))
ax = fig.add_subplot(projection ='3d')
plotDeformation(soln,ax)
plt.show()

## Shape optimization
# Define nonlinear constraint
nlc = opt.NonlinearConstraint(computeConstraint, 0, 0)
# Define bounds
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

bnd = opt.Bounds(low,upp,keep_feasible=False)

# Set initial design
u = coefs[interiormask]
u[:,2] = 0.01*L*np.sin(u[:,0] * (np.pi))*np.sin(u[:,1] * (np.pi))
u = u.flatten()

# Shape optimization
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

# Plotting
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