import sim_lorenz
import numpy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

lorenz = sim_lorenz.lorenz(20000)
data = lorenz.field(1.5, 400.0, ((-1, 1), (-1, 1), (-1, 1)), 4, 0.1)

#fig = pyplot.figure()
#ax = fig.add_subplot(111, projection="3d")
for run in data:
    x = run[:,0]
    y = run[:,1]
    z = run[:,2]
    x_dot = numpy.abs(numpy.gradient(x))
    y_dot = numpy.abs(numpy.gradient(y))
    z_dot = numpy.abs(numpy.gradient(z))
    v = numpy.sqrt(x_dot**2 + y_dot**2 + z_dot**2)
    v = v/numpy.max(v)
    c = numpy.zeros(shape=(len(x_dot), 3))
    c[:, 0] = v
    c[:, 1] = v
    c[:, 2] = v
#    pyplot.plot(x)
#    pyplot.plot(y)
#    pyplot.plot(z)
#    pyplot.show()
#    pyplot.plot(x_dot)
#    pyplot.plot(y_dot)
#    pyplot.plot(z_dot)
#    pyplot.show()

    f, ((a1, a2), (a3, a4)) = pyplot.subplots(2, 2, sharex="col", sharey="row")
    a1.set_title("The Lorenz system")
    a1.plot(x,y, color="black")
    a1.set_xlabel("x")
    a1.set_xlabel("y")
    a2.plot(x,z, color="black")
    a2.set_xlabel("x")
    a2.set_xlabel("z")
    a3.plot(y,z, color="black")
    a3.set_xlabel("y")
    a3.set_xlabel("z")
    a4 = f.add_subplot(224, projection="3d")
    a4.plot(x, y, z, color="black")
    a4.set_xlabel("x")
    a4.set_ylabel("y")
    a4.set_zlabel("z")
    pyplot.show()

#    fig = pyplot.figure()
#    a4 = fig.add_subplot(2,1,2, projection="3d")
#    ax.scatter(x, y, z, c=c, s=5**2, edgecolors="none")
#    pyplot.show()
