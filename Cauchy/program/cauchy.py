#!/usr/bin/python

from dolfin import *
import numpy as np
from numpy import linalg as la
from scipy.optimize import minimize
# from scipy import optimize


""" serial do all interpolation operations, slow but no need to store files"""
""" have to test parallel vs. serial """
parameters["reorder_dofs_serial"] = False

# Create mesh and define function space
DIMX, DIMY = 320,320
mesh = UnitSquareMesh(DIMX, DIMY)
V = FunctionSpace(mesh, 'Lagrange', 1)
# TRUE SOLUTION
u0 = Expression('x[0] + x[1]')
# your guess have to be extension of the boundary condition, otherwise the 
# program won't give a accurate solution.
u1 = Expression('x[0]')

# g = psi
g = Expression('1.0')
f = Expression('0.0')

""" what is the best parameter here ?"""
reg = 1e-5

exact = interpolate(u0,V).vector().array()
# settings for coefficient
# psi = Expression('x[0]*x[0] + 1')
psi = Expression('1.0')
phi = Expression('1.0')
# phi = Expression('x[1]*x[1] + 1')
"""psi, phi arrays for computing"""
psi_array = interpolate(psi,V).vector().array()
phi_array = interpolate(phi,V).vector().array()



def set2bd(init_array, p):
  init_array[DIMX+1:(DIMY+1)*(DIMX+1):DIMX+1] = p[0:DIMY]
def set4bd(init_array,p):
  # init_array[DIMX:(DIMY+1)*(DIMX+1)-1:DIMX+1] = p[DIMY+DIMX:DIMY+DIMX+DIMY]
  init_array[2*DIMX+1:(DIMY+1)*(DIMX+1)-1:DIMX+1] = p[DIMY+DIMX:DIMY+DIMX+DIMY-1]
  """release control to the Dirichlet bd"""
def set3bd(init_array,p):
  init_array[DIMY*(DIMX+1)+1:(DIMY+1)*(DIMX+1)] = p[DIMY:DIMY+DIMX]

def ex2bd(u_array,p):
  p[0:DIMY] = u_array[DIMX+1:(DIMY+1)*(DIMX+1):DIMX+1] 
def ex4bd(u_array,p):
  # p[DIMY+DIMX:DIMY+DIMX+DIMY] = u_array[DIMX:(DIMY+1)*(DIMX+1)-1:DIMX+1] 
  p[DIMY+DIMX:DIMY+DIMX+DIMY-1] = u_array[2*DIMX+1:(DIMY+1)*(DIMX+1)-1:DIMX+1] 
def ex3bd(u_array,p):
  p[DIMY:DIMY+DIMX] = u_array[DIMY*(DIMX+1)+1:(DIMY+1)*(DIMX+1)] 



# function ----> bd array
def extractbd(u):
  # p=np.zeros(2*DIMY+DIMX)
  p=np.zeros(2*DIMY+DIMX-1)
  u_array = u.vector().array()
  ex2bd(u_array,p)
  ex3bd(u_array,p)
  ex4bd(u_array,p) 
  return p
# function ----> bd array
def extractinnerbd(u):
  # p=np.zeros(2*DIMY+DIMX)
  p=np.zeros(2*DIMY+DIMX-1)
  u_array = u.vector().array()
  p[0:DIMY] = u_array[DIMX+2:(DIMY+1)*(DIMX+1):DIMX+1] 
  p[DIMY:DIMY+DIMX] = u_array[(DIMY-1)*(DIMX+1)+1:DIMY*(DIMX+1)] 
  # p[DIMY+DIMX:DIMY+DIMX+DIMY] = u_array[DIMX-1:(DIMY+1)*(DIMX+1)-2:DIMX+1] 
  p[DIMY+DIMX:DIMY+DIMX+DIMY-1] = u_array[2*DIMX:(DIMY+1)*(DIMX+1)-2:DIMX+1] 
  return p
# bd array ----> function(interpolate with 0)
def bd2func(p):
  init = u1 # not u0
  # init = u0
  init_array = interpolate(init,V).vector().array()
  # assemble
  set2bd(init_array,p)
  set3bd(init_array,p)
  set4bd(init_array,p)

  data0 = array2func(V,init_array)
  return data0

def goal(p):
  # p as an array of 2DIMY + DIMX
  data0 = bd2func(p)
  sol = solver(mesh,V,data0, g,f,psi,phi)
  sol_array = sol.vector().array()
  # NORM = la.norm(sol_array[0:DIMX] - exact[0:DIMX],2)
  NORM = la.norm(sol_array[0:DIMX+1] - exact[0:DIMX+1],2)
  # regularization
  p2 = np.array(p[0:DIMY])
  p3 = np.array(p[DIMY:DIMY+DIMX])
  p4 = np.array(p[DIMY+DIMX:DIMY+DIMX+DIMY-1])
  NORMP2 = la.norm(p2[0:-1]-p2[1::],2)
  NORMP3 = la.norm(p3[0:-1]-p3[1::],2)
  NORMP4 = la.norm(p4[0:-1]-p4[1::],2)
  return 0.5*NORM*NORM/DIMX + 0.5*reg*(NORMP2*NORMP2 + NORMP3*NORMP3 + NORMP4*NORMP4)/DIMX


def diff(p):
  p1 = np.concatenate([p,[0]])
  p2 = np.concatenate([[0],p])
  return np.array(p1) - np.array(p2)

def derivative(p):
  # p as an array
  data0 = bd2func(p)
  sol = solver(mesh,V,data0, g,f,psi,phi)     
  sol_array = sol.vector().array()
  """ setup short names """
  gg = array2func(V,sol_array - exact) # Neumann
  # for later use
  # multiplier = np.zeros(2*DIMY+DIMX)
  multiplier = np.zeros(2*DIMY+DIMX-1)

  ex2bd(psi_array,multiplier)
  ex3bd(phi_array,multiplier)
  ex4bd(psi_array,multiplier)
  
  u2 = Constant((0.0)) # Dirichlet
  ff = Constant((0.0)) # force
  w = solver(mesh,V,u2, gg, ff, psi, phi)
  w_barray = extractbd(w)
  w_iarray = extractinnerbd(w)
  diffw = (w_barray - w_iarray)
  grad = np.array(diffw)*np.array(multiplier)
  """ regularization term use int |grad p|^2"""
  p2 = np.array(p[0:DIMY])
  p3 = np.array(p[DIMY:DIMY+DIMX])
  p4 = np.array(p[DIMY+DIMX:DIMY+DIMX+DIMY-1])
  pp2 = p2[0:-1] - p2[1::]
  pp3 = p3[0:-1] - p3[1::]
  pp4 = p4[0:-1] - p4[1::]
  ppp2 = diff(pp2)
  ppp3 = diff(pp3)
  ppp4 = diff(pp4)
  diffp = np.concatenate([ppp2,ppp3,ppp4])
  # TODO HERE ########################
  reggrad = reg*np.array(diffp)/DIMX
  return grad + reggrad

def bfgs_solver():
  # u1 = Expression('x[0] + x[1]')
  x0 = extractbd(interpolate(u1,V))

  res = minimize(goal,x0, method = 'L-BFGS-B', jac=derivative, \
     options = {'disp': True,
            'iprint': 0,
            'maxcor': 10,
            'ftol': (1e-30)* np.finfo(float).eps,
            'gtol': 1e-30,
            'eps': 1e-30,
            'maxfun': 150000,
            'maxiter': 150000
          })
  
  return res

def solver(mesh, V, udb, g, f, psi, phi):
  class DirichletBoundary(SubDomain):
      def inside(self, x, on_boundary):
          tol = 1E-15   # tolerance for coordinate comparisons
          return on_boundary and \
                 (abs(x[0]) < tol or abs(x[0] - 1) < tol or abs(x[1] - 1) < tol)

  u0_boundary = DirichletBoundary()
  bc = DirichletBC(V, udb, u0_boundary)


  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  # f = Expression('2.0')
  # set up coefficients

  # variational form
  mgrad = as_vector((psi*u.dx(0),phi*u.dx(1)))
  a = inner(mgrad, grad(v))*dx
  L = -f*v*dx - phi*g*v*ds # here negative on Neumann term, since outward normal vector means negative gradient

  # Compute solution
  A = assemble(a)
  b = assemble(L)
  u = Function(V)
  U = u.vector()
  bc.apply(A,b)
  # solve(a == L, u, bc)
  solve(A,U,b)
  return u

#array ----> func
def array2func(V,phi_array):
  psi = Function(V)
  psi.vector()[:] = phi_array
  return psi



def main():
  # plot(interpolate(u0,V))
  res = bfgs_solver()
  p = res.x
  data0 = bd2func(p)
  sol = solver(mesh,V,data0, g,f,psi,phi)
  sol_array = sol.vector().array()
  NORM = la.norm(sol_array -exact,2)

  print 'L2 error:', NORM/la.norm(exact,2)
  plot(sol)
  interactive()



if __name__ == '__main__':
  main()
