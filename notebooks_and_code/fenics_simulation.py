import os, sys

import numpy as np
import matplotlib.pyplot as plt

from dolfin import *
import ufl_legacy as ufl
from ufl_legacy import Index

from mpl_toolkits.mplot3d import Axes3D

mesh_idx = int(sys.argv[1])
mesh=Mesh()
with XDMFFile("meshes/mesh_2d_"+str(mesh_idx).zfill(4)+".xdmf") as infile:
#with XDMFFile("/home/csaba/Documents/only_spiral.xdmf") as infile:
    infile.read(mesh)
    
mvc = MeshValueCollection("size_t", mesh, 1) 
with XDMFFile("meshes/mesh_1d_"+str(mesh_idx).zfill(4)+".xdmf") as infile:
    infile.read(mvc)
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

parameters["form_compiler"]["quadrature_degree"] = 4

output_dir = "output/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
rho_g = (1200-1000)*9.81 / ((1e2)**2)   #1200 * 9.81 kg /m^2 s = 1200*9.81 e-4 kg / cm^2, density SU-8
#L = 3.048
E, nu = 2.0E9/1e2, 0.22 #CONVERSION TO CM #Young's Modulus SU-8
mu = E/(2.0*(1.0 + nu))
lmbda = 2.0*mu*nu/(1.0 - 2.0*nu)
t = Constant(1e-4) #1e-4 cm = 1e-6 = 1 micron

initial_shape = Expression(('x[0]','x[1]','0'), degree = 4)
V_phi =  FunctionSpace(mesh, VectorElement("P", triangle, degree = 2, dim = 3))
phi0 = project(initial_shape, V_phi)

# Given the midplane, we define the corresponding unit normal as below and
# project on a suitable function space (here :math:`P_1` but other choices
# are possible)::

def normal(y):
    n = cross(y.dx(0), y.dx(1))
    return n/sqrt(inner(n,n))

V_normal = FunctionSpace(mesh, VectorElement("P", triangle, degree = 1, dim = 3))
n0 = project(normal(phi0), V_normal)

# The kinematics of the Nadghi shell model is defined by the following
# vector fields :
# 
# - :math:`\phi`: the position of the midplane, or the displacement from the reference configuration :math:`u = \phi - \phi_0`:
# - :math:`d`: the director, a unit vector giving the orientation of the microstructure
# 
# We parametrize the director field by two angles, which correspond to spherical coordinates,
# so as to explicitly resolve the unit norm constraint (see [3])::

def director(beta):
    return as_vector([sin(beta[1])*cos(beta[0]), -sin(beta[0]), cos(beta[1])*cos(beta[0])])

# We assume that in the initial configuration the director coincides with
# the normal. Hence, we can define the angles :math:`\beta`: for the initial
# configuration as follows: ::

beta0_expression = Expression(["atan2(-n[1], sqrt(pow(n[0],2) + pow(n[2],2)))",
                               "atan2(n[0],n[2])"], n = n0, degree=4)

V_beta = FunctionSpace(mesh, VectorElement("P", triangle, degree = 2, dim = 2))
beta0 = project(beta0_expression, V_beta)

# The director in the initial configuration is then written as ::

d0 = director(beta0)

P2 = FiniteElement("Lagrange", triangle, degree = 2)
bubble = FiniteElement("B", triangle, degree = 3)
enriched = P2 + bubble

element = MixedElement([VectorElement(enriched, dim=3), VectorElement(P2, dim=2)])

Q = FunctionSpace(mesh, element)

# Then, we define :py:class:`Function`, :py:class:`TrialFunction` and :py:class:`TestFunction` objects to express the variational forms and we split them into each individual component function::
    
q_, q, q_t = Function(Q), TrialFunction(Q), TestFunction(Q)
u_, beta_ = split(q_)

# The gradient of the transformation and the director in the current
# configuration are given by::

F = grad(u_) + grad(phi0)
d = director(beta_ + beta0)

# With the following definition of the natural metric and curvature ::

a0 = grad(phi0).T*grad(phi0)
b0 = -0.5*(grad(phi0).T*grad(d0) + grad(d0).T*grad(phi0))

# The membrane, bending, and shear strain measures of the Naghdi model are defined by::

e = lambda F: 0.5*(F.T*F - a0)
k = lambda F, d: -0.5*(F.T*grad(d) + grad(d).T*F) - b0
gamma = lambda F, d: F.T*d - grad(phi0).T*d0

# Using curvilinear coordinates,  and denoting by ``a0_contra`` the
# contravariant components of the metric tensor :math:`a_0^{\alpha\beta}` (in the initial curved configuration)
# the constitutive equation is written in terms of the matrix :math:`A` below,
# representing the contravariant components of the constitutive tensor
# for isotropic elasticity in plane stress (see *e.g.* [4]).
# We use the index notation offered by UFL to express
# operations between tensors::

a0_contra = inv(a0)
j0 = det(a0)

i, j, l, m = Index(), Index(), Index(), Index()
A_ = as_tensor((((2.0*lmbda*mu)/(lmbda + 2.0*mu))*a0_contra[i,j]*a0_contra[l,m]
                + 1.0*mu*(a0_contra[i,l]*a0_contra[j,m] + a0_contra[i,m]*a0_contra[j,l]))
                ,[i,j,l,m])

# The normal stress :math:`N`, bending moment :math:`M`, and shear stress :math:`T` tensors are (they are purely Lagrangian stress measures,
# similar to the so called 2nd Piola stress tensor in 3D elasticity)::

N = as_tensor(t*A_[i,j,l,m]*e(F)[l,m], [i,j])
M = as_tensor((t**3/12.0)*A_[i,j,l,m]*k(F,d)[l,m],[i,j])
T = as_tensor(t*mu*a0_contra[i,j]*gamma(F,d)[j], [i])

# Hence, the contributions to the elastic energy density due to membrane, :math:`\psi_m`,
# bending, :math:`\psi_b`, and shear, :math:`\psi_s` are
# (they are per unit surface in the initial configuration)::

psi_m = 0.5*inner(N, e(F))
psi_b = 0.5*inner(M, k(F,d))
psi_s = 0.5*inner(T, gamma(F,d))

# Shear and membrane locking is treated using the partial reduced
# selective integration proposed in Arnold and Brezzi [2]. In this approach
# shear and membrane energy are splitted as a sum of two contributions
# weighted by a factor :math:`\alpha`. One of the two contributions is
# integrated with a reduced integration. While [1] suggests a 1-point
# reduced integration, we observed that this leads to spurious modes in
# the present case. We use then :math:`2\times 2`-points Gauss integration
# for a portion :math:`1-\alpha` of the energy, whilst the rest is
# integrated with a :math:`4\times 4` scheme. We further refine the
# approach of [1] by adopting an optimized weighting factor
# :math:`\alpha=(t/h)^2`, where :math:`h` is the mesh size. ::

dx_h = dx(metadata={'quadrature_degree': 2})
h = CellDiameter(mesh)
alpha = project(t**2/h**2, FunctionSpace(mesh,'DG',0))

Pi_PSRI = psi_b*sqrt(j0)*dx + alpha*psi_m*sqrt(j0)*dx + alpha*psi_s*sqrt(j0)*dx + (1.0 - alpha)*psi_s*sqrt(j0)*dx_h + (1.0 - alpha)*psi_m*sqrt(j0)*dx_h

# Hence the total elastic energy and its first and second derivatives are ::

Pi = Pi_PSRI
dPi = derivative(Pi, q_, q_t)
J = derivative(dPi, q_, q)

bc_clamped = DirichletBC(Q, project(q_, Q), mf, 1) #from facet mesh, from meshvaluefunction (mvc->mf)
bcs=[bc_clamped]

class NonlinearProblemPointSource(NonlinearProblem):

    def __init__(self, L, a, bcs):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.bcs = bcs
        self.P = 0.0

    def F(self, b, x):
        weight=(1-alpha)*inner(Constant((0.0,0.0,-P*rho_g*1e-4)),u_)*dx_h
        w=MeshCoordinates(mesh)
        #weight=weight+(1-alpha)*Constant((-P/5*mg_org_per_column))*exp(-sqrt(w[0]**2+w[1]**2)**4/0.06**4)*u_[2]*dx_h #gravitational work
        assemble(derivative(self.L-weight, q_, q_t), tensor=b)
        for bc in self.bcs:
            bc.apply(b, x)

    def J(self, A, x):
        assemble(self.a, tensor=A)
        for bc in self.bcs:
            bc.apply(A, x)

problem = NonlinearProblemPointSource(Pi, J, bcs)
# We use a standard Newton solver and setup the files for the writing the
# results to disk::

solver = NewtonSolver()
output_dir = "output/"
file_phi = File(output_dir + str(mesh_idx).zfill(4)+"_configuration.pvd")
file_strain = File(output_dir + str(mesh_idx).zfill(4)+"_strain.pvd")
# Finally, we can solve the quasi-static problem, incrementally increasing
# the loading from :math:`0` to :math:`2000` N::

#P_values = np.linspace(0.0, 0.02, 40)
P_values = np.linspace(0.0, 1, 100)    
print(P_values[-1])
displacement = 0.*P_values
q_.assign(project(Constant((0,0,0,0,0)), Q))

for (i, P) in enumerate(P_values):
    problem.P = P
    (niter,cond) = solver.solve(problem, q_.vector())

    phi = project(u_ + phi0, V_phi)
    #displacement[i] = phi(0.0, 0.0)[2] - phi0(0.0, 0.0)[2]
    print(np.min(project(u_[2]).vector()[:]))
    phi.rename("phi", "phi")
    file_phi << (phi, P)
    print("Increment %d of %s. Converged in %2d iterations. P:  %.2f, Displ: %.2f" %(i, P_values.size,niter,P, displacement[i]))

    en_function = project(psi_m + psi_b + psi_s, FunctionSpace(mesh, 'Lagrange', 1))
    strain_function = project(k(F,d))
    
    strain_function.rename("Strain", "Strain")
    en_function.rename("Elastic Energy", "Elastic Energy")
    file_strain << (strain_function,P)

