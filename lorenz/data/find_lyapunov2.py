
import numpy
from matplotlib import pyplot
import sys

def progress(x, max, div=20):
    i = div*x//max + 1
    sys.stdout.write("\r")
    sys.stdout.write(("[%-" + str(div) + "s] %d%%") % ("="*i, (100/div)*i))
    sys.stdout.flush()

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
        dff = 0.000001
        n1, n2, n3 = [], [], []
        y = numpy.empty((self.nmax, len(y0)))
        y[0] = y0

        #initalize pertubation
        S = y[0, :3]
        A = y[0, 3:6]
        B = y[0, 6:9]
        C = y[0, 9:12]
        N1, V1, A, B = norm(A, B, C)
        A = orthogonalize(V1, A)
        B = orthogonalize(V1, B)
        N2, V2, A, _ = norm(A, B)
        A = orthogonalize(V2, A)
        N3, V3, _, _ = norm(A)

        y[0, 3:6] = S + V1*dff
        y[0, 6:9] = S + V2*dff
        y[0, 9:12] = S + V3*dff
        

        for i in range(0, self.nmax - 1):
            y[i + 1] = self._iteration(y[i])
            if (i + 1) % self.renorm == 0:
                progress(i, self.nmax - 1, 50)
                S = y[i+1, :3]
                A = y[i+1, 3:6] - S
                B = y[i+1, 6:9] - S
                C = y[i+1, 9:12] - S
                N1, V1, A, B = norm(A, B, C)
                A = orthogonalize(V1, A)
                B = orthogonalize(V1, B)
                N2, V2, A, _ = norm(A, B)
                A = orthogonalize(V2, A)
                N3, V3, _, _ = norm(A)

                y[i+1, 3:6] = S + V1*dff
                y[i+1, 6:9] = S + V2*dff
                y[i+1, 9:12] = S + V3*dff
                n1 = n1 + [numpy.log(N1/dff)/self.h]
                n2 = n2 + [numpy.log(N2/dff)/self.h]
                n3 = n3 + [numpy.log(N3/dff)/self.h]
        n1 = numpy.array(n1)/self.renorm
        n2 = numpy.array(n2)/self.renorm
        n3 = numpy.array(n3)/self.renorm
        print("\n")
        print("N1:", numpy.mean(n1), numpy.std(n1))
        print("N2:", numpy.mean(n2), numpy.std(n2))
        print("N3:", numpy.mean(n3), numpy.std(n3))
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
        sigma = 10
        rho = 28
        beta = 8/3.0
        x1=yn[0]
        y1=yn[1]
        z1=yn[2]
        x2=yn[3]
        y2=yn[4]
        z2=yn[5]
        x3=yn[6]
        y3=yn[7]
        z3=yn[8]
        x4=yn[9]
        y4=yn[10]
        z4=yn[11]
#        return numpy.array([y-x, -x*z, x*y-R,
#                            y_pert1 - x_pert1, -x*z_pert1 - z*x_pert1, x*y_pert1 + y*x_pert1,
#                            y_pert2 - x_pert2, -x*z_pert2 - z*x_pert2, x*y_pert2 + y*x_pert2,
#                            y_pert3 - x_pert3, -x*z_pert3 - z*x_pert3, x*y_pert3 + y*x_pert3])
        return numpy.array([sigma*(y1-x1), x1*(rho - z1) - y1, x1*y1-beta*z1,
                            sigma*(y2-x2), x2*(rho - z2) - y2, x2*y2-beta*z2,
                            sigma*(y3-x3), x3*(rho - z3) - y3, x3*y3-beta*z3,
                            sigma*(y4-x4), x4*(rho - z4) - y4, x4*y4-beta*z4])
    return A

class lorenz:
    def __init__(self, resolution=20000):
        self.resolution = resolution
        
    def run(self, R, speed, inital_conditions, renorm=10):
        solver = RK4(0.0, 1.0*speed, self.resolution, renorm, set_A(R))
        path = solver.compute(numpy.array(inital_conditions + [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]))
        return path

from matplotlib import pyplot
n = 10
solver = lorenz(10000*n)
path = solver.run(1.0, 1*n, [-1, 2, 1], 10) #8415.96,8412.00,8399.43

print("True values:", 0.9056, 0, -14.5723)
pyplot.plot(path[:, 0], path[:, 1])
pyplot.show()
#S = path[:, :3]
#A = path[:, 3:6]
#B = path[:, 6:9]
#C = path[:, 9:12]

#pyplot.plot(A[:, 0] - S[:, 0])
#pyplot.plot(B[:, 0] - S[:, 0])
#pyplot.plot(C[:, 0] - S[:, 0])
#pyplot.show()
