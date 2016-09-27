
import numpy
from matplotlib import pyplot

def norm(A, B=numpy.array([0.0, 0.0, 0.0]), C=numpy.array([0.0, 0.0, 0.0])):
    nA = numpy.linalg.norm(A)
    nB = numpy.linalg.norm(B)
    nC = numpy.linalg.norm(C)
    if nA >= nB:
        if nA >= nC:
            return nA, A/nA, B, C
        else:
            return nC, C/nC, A, B
    else:
        if nB >= nC:
            return nB, B/nB, A, C
        else:
            return nC, C/nC, A, B

def orthogonalize(U, A):
    return A - numpy.dot(U, A)*U

class RK4:
    def __init__(self, tmin, tmax, n, renorm, A):
        self.h = (tmax - tmin)/float(n)
        self.tmin = tmin
        self.tmax = tmax
        self.nmax = n
        self.renorm = renorm
        self.A = A
        
    def compute(self, y0):
        n1, n2, n3 = 0.0, 0.0, 0.0
        y = numpy.empty((self.nmax, len(y0)))
        y[0] = y0
        for i in range(0, self.nmax - 1):
            y[i + 1] = self._iteration(y[i])
            if (i + 1) % self.renorm == 0:
                A = y[i+1, 3:6]
                B = y[i+1, 6:9]
                C = y[i+1, 9:12]
                N1, V1, A, B = norm(A, B, C)
                A = orthogonalize(V1, A)
                B = orthogonalize(V1, B)
                N2, V2, A, B = norm(A, B)
                A = orthogonalize(V2, A)
                N3, V3 = numpy.linalg.norm(A), A/numpy.linalg.norm(A)
                y[i+1, 3:6] = V1
                y[i+1, 6:9] = V2
                y[i+1, 9:12] = V3
                n1 = n1 + numpy.log(N1)
                n2 = n2 + numpy.log(N2)
                n3 = n3 + numpy.log(N3)
        print(n1/(i*self.h), n2/(i*self.h), n3/(i*self.h))
        return y
    
    def tscale(self):
        return numpy.linspace(self.tmin, self.tmin + self.h*self.nmax, self.nmax)
        
    def _iteration(self, yn):
        k1 = self.A(yn)
        k2 = self.A(yn + self.h*k1/2.0)
        k3 = self.A(yn + self.h*k2/2.0)
        k4 = self.A(yn + self.h*k3)
        return (yn + self.h*(k1 + k4 + 2*(k2 + k3))/6.0)

def set_A(R):
    def A(yn):
        sigma = 16 
        rho = 45.92
        beta = 4
        x=yn[0]
        y=yn[1]
        z=yn[2]
        x_pert1=yn[3]
        y_pert1=yn[4]
        z_pert1=yn[5]
        x_pert2=yn[6]
        y_pert2=yn[7]
        z_pert2=yn[8]
        x_pert3=yn[9]
        y_pert3=yn[10]
        z_pert3=yn[11]
#        return numpy.array([y-x, -x*z, x*y-R,
#                            y_pert1 - x_pert1, -x*z_pert1 - z*x_pert1, x*y_pert1 + y*x_pert1,
#                            y_pert2 - x_pert2, -x*z_pert2 - z*x_pert2, x*y_pert2 + y*x_pert2,
#                            y_pert3 - x_pert3, -x*z_pert3 - z*x_pert3, x*y_pert3 + y*x_pert3])
        return numpy.array([sigma*(y-x), x*(rho - z) - y, x*y-beta*z,
                            sigma*(y_pert1 - x_pert1), -x*z_pert1 + (rho - z)*x_pert1 - y_pert1, x*y_pert1 + y*x_pert1 - beta*z_pert1,
                            sigma*(y_pert2 - x_pert2), -x*z_pert2 + (rho - z)*x_pert2 - y_pert2, x*y_pert2 + y*x_pert2 - beta*z_pert2,
                            sigma*(y_pert3 - x_pert3), -x*z_pert3 + (rho - z)*x_pert3 - y_pert3, x*y_pert3 + y*x_pert3 - beta*z_pert3])
    return A

class lorenz:
    def __init__(self, resolution=20000):
        self.resolution = resolution
        
    def run(self, R, speed, inital_conditions, renorm=10):
        solver = RK4(0.0, 1.0*speed, self.resolution, renorm, set_A(R))
        path = solver.compute(numpy.array(inital_conditions + [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]))
#        path = solver.compute(numpy.array(inital_conditions))
        return path

from matplotlib import pyplot
n = 100
solver = lorenz(10000*n)
path = solver.run(1.0, 1*n, [-1, 2, 3], 10) #8415.96,8412.00,8399.43
pyplot.plot(path[:, 0], path[:, 1])
pyplot.show()
A = path[:, 3:6]
B = path[:, 6:9]
C = path[:, 9:12]
pyplot.plot(numpy.sqrt(A[:, 0]**2 + A[:, 1]**2 + A[:, 2]**2))
pyplot.plot(numpy.sqrt(B[:, 0]**2 + B[:, 1]**2 + B[:, 2]**2))
pyplot.plot(numpy.sqrt(C[:, 0]**2 + C[:, 1]**2 + C[:, 2]**2))
pyplot.show()
