import sim_lorenz
import numpy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

lorenz = sim_lorenz.lorenz(2000)
data = lorenz.field(3, 40.0, ((-1, 1), (-1, 1), (-1, 1)), 2, 0.1)

shape = data.shape
data = data.reshape(shape[0]*shape[1], shape[2])
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]
hist, bins = numpy.histogramdd((x, y, z), bins=[100, 100, 100])
hist = hist/numpy.max(hist)

time = 10000#ms
frames = []
fig = pyplot.figure()
for i in range(len(hist)):
    frame = pyplot.imshow(hist[i], cmap=pyplot.get_cmap("gray"), interpolation="nearest", animated=True, vmin=0.0, vmax=1.0)
    frames.append([frame])
ani = animation.ArtistAnimation(fig, frames, interval=time/len(hist))
pyplot.show()
