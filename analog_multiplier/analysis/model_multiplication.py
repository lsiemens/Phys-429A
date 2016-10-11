from matplotlib import pyplot
import tektronix
import numpy
import scipy.optimize
import os

def theory(vin, bias, I1, I2):
    return bias - 2*R2*I1*numpy.sqrt(1 + (vin/R1*I2)**2)

def general(vin, a, b, c):
    return a + b*numpy.sqrt(1 + (vin/c)**2)

def mult(vx, vy, a):
    return a*(numpy.sqrt(1 + (vx + vy)**2) - numpy.sqrt(1 + (vx - vy)**2))

tek = tektronix.tektronix()

dir = "../data/"
files = os.listdir(dir)
fnames = [file for file in files if file.startswith("s1")]

for fname in fnames:
    data = tek.load_file(dir + fname)
    
    time = data[0, 0]
    
    Vx = data[0, 1]
    Vy = data[1, 1]
    Vo = data[2, 1]
    
    pyplot.plot(time, Vx, label="Vx")
    pyplot.plot(time, Vy, label="Vy")
    pyplot.plot(time, Vo, label="Vout")
    pyplot.plot(time, mult(Vx, Vy, -0.7), label="Best fit")
    pyplot.legend(loc=0)
    pyplot.show()
    
