import numpy
from matplotlib import pyplot

data = numpy.loadtxt("solution.txt")
pyplot.plot(data[:, 0])
pyplot.plot(data[:, 1])
pyplot.plot(data[:, 2])
pyplot.show()
pyplot.plot(data[:, 0], data[:, 1])
pyplot.show()
