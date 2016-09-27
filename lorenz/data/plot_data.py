import loggerpro_csv
import numpy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

labpro = loggerpro_csv.loggerpro(4)
data = labpro.load("lorenz_deep.csv")

#fig = pyplot.figure()
#ax = fig.add_subplot(111, projection="3d")
for run in data:
    t = run[0]
    x = run[1]
    y = run[2]
    z = run[3]
    dt = numpy.gradient(t)
    x_dot = numpy.abs(numpy.gradient(x)/dt)
    y_dot = numpy.abs(numpy.gradient(y)/dt)
    z_dot = numpy.abs(numpy.gradient(z)/dt)
    v = numpy.sqrt(x_dot**2 + y_dot**2 + z_dot**2)
    v = v/numpy.max(v)
    c = numpy.zeros(shape=(len(x_dot), 3))
    c[:, 0] = v
    c[:, 1] = v
    c[:, 2] = v
    R, Rerr = numpy.mean(run[4]), numpy.std(run[4])
    print(R, Rerr)
#    pyplot.plot(t, x)
#    pyplot.show()
#    pyplot.plot(t, y)
#    pyplot.show()
#    pyplot.plot(t, z)
#    pyplot.show()
#    pyplot.plot(t, x_dot, t, y_dot, t, z_dot)
#    pyplot.show()

    f, ((a1, a2), (a3, a4)) = pyplot.subplots(2, 2, sharex="col", sharey = 'row')
    a1.set_title("Circuit Voltage")
    a1.plot(x,y, color="black")
    a1.set_xlabel("x in $v$")
    a1.set_ylabel("y in $v$")
    a2.plot(x,z, c="black")
    a2.set_xlabel("x in $v$")
    a2.set_ylabel("z in $v$")
    a3.plot(y,z, c="black")
    a3.set_xlabel("y in $v$")
    a3.set_ylabel("z in $v$")
    a4 = f.add_subplot(224, projection="3d")
    a4.plot(x,y,z, c="black")
    a4.set_xlabel("x")
    a4.set_ylabel("y")
    a4.set_zlabel("y")
    pyplot.show()

#    fig = pyplot.figure()
#    ax = fig.add_subplot(111, projection="3d")
#    ax.scatter(x, y, z, c=c, s=5**2, edgecolors="none")
#pyplot.show()
