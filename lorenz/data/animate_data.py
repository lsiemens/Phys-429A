import loggerpro_csv
import numpy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

labpro = loggerpro_csv.loggerpro(4)
data = labpro.load("lorenz_deep.csv")
data_continued = labpro.load("lorenz_deep_continued.csv")

data = numpy.concatenate((data, data_continued))
if True:
    data = numpy.swapaxes(data, 1, 2)
    shape = data.shape
    data = data.reshape(shape[0]*shape[1], shape[2])
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]
    if True:
 #   for run in [data[0]]:
 #       t = run[0]
 #       x = run[1]
 #       y = run[2]
 #       z = run[3]
        frames = 360
    
        fig = pyplot.figure()
        ax = Axes3D(fig)
        
        def init():
#            ax.scatter(x, y, z, c='k', s=2**2, marker=".", alpha=0.3)
            ax.scatter(x, y, z, c='k', s=1**2, marker=".", alpha=0.015)
            ax.view_init(elev=10., azim=0.0)
            ax.grid(False)
            ax.autoscale_view("tight")
            ax.margins()
            ax.set_xlabel("$X$ in volts")
            ax.set_ylabel("$Y$ in volts")
            ax.set_zlabel("$Z$ in volts")
            ax.set_title("Simplified Lorenz attractor")

        def animate(i):
            print(frames, i)
            ax.view_init(elev=10., azim=i*360.0/frames)
        
        width = 9
        fig.set_size_inches(width, width)
        dpi = int(1000/width)
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(frames))
        writer = animation.writers["ffmpeg"](fps=30, bitrate=2000)
        anim.save("animation.mp4", writer=writer, dpi=dpi)
#        pyplot.show()
    
if False:
    data = numpy.swapaxes(data, 1, 2)
    shape = data.shape
    data = data.reshape(shape[0]*shape[1], shape[2])
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]
    
    hist, bins = numpy.histogramdd((x, y, z), bins=[100,100,100])
    hist = hist/numpy.max(hist)
    
    time = 10000#ms
    frames = []
    fig = pyplot.figure()
    for i in range(len(hist)):
        frame = pyplot.imshow(hist[i], cmap=pyplot.get_cmap("gray"), interpolation="nearest", vmin=0.0, vmax=1.0, animated=True)
        frames.append([frame])
    ani = animation.ArtistAnimation(fig, frames, interval=time/len(hist))
    ani.save("animation.mp4", fps=30)
    #pyplot.show()
    
