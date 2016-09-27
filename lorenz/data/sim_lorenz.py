import numpy
from matplotlib import pyplot

class RK4:
    def __init__(self, tmin, tmax, n, A):
        self.h = (tmax - tmin)/float(n)
        self.tmin = tmin
        self.tmax = tmax
        self.nmax = n
        self.A = A
        
    def compute(self, y0):
        y = numpy.empty((self.nmax, len(y0)))
        y[0] = y0
        for i in range(0, self.nmax - 1):
            y[i + 1] = self._iteration(y[i])
        return y
    
    def tscale(self):
        return numpy.linspace(self.tmin, self.tmin + self.h*self.nmax, self.nmax)
        
    def _iteration(self, yn):
        k1 = self.A(yn)
        k2 = self.A(yn + self.h*k1/2.0)
        k3 = self.A(yn + self.h*k2/2.0)
        k4 = self.A(yn + self.h*k3)
        return (yn + self.h*(k1 + k4 + 2*(k2 + k3))/6.0)

def set_A(R, speed):
    def A(yn):
        x=yn[0]
        y=yn[1]
        z=yn[2]
#        return speed*numpy.array([sig*(y-x), x*(r - z), x*y-b*z])
        return speed*numpy.array([y-x, -x*z, x*y-R])
    return A

class lorenz:
    def __init__(self, resolution=20000):
        self.resolution = resolution
        
    def run(self, R, speed, inital_conditions, skip=0.0):
        solver = RK4(0.0, 1.0, self.resolution, set_A(R, speed))
        path = solver.compute(numpy.array(inital_conditions))
        iskip = int(self.resolution*skip)
        return path[iskip:]
    
    def field(self, R, speed, range, res, skip=0.0):
        paths = []
        for x in numpy.linspace(range[0][0], range[0][1], res):
            for y in numpy.linspace(range[1][0], range[1][1], res):
                for z in numpy.linspace(range[2][0], range[2][1], res):
                    paths.append(self.run(R, speed, [x, y, z]))
        return numpy.array(paths)
