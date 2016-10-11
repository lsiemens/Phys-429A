from matplotlib import pyplot
import tektronix
import numpy
import scipy.optimize
import os

def get_resistor_units(value):
    if value[-1] == "k":
        value = int(value[:-1]) * 1000
    elif value[-1] == "m":
        value = int(value[:-1]) * 1000000
    else:
        value = int(value)
    return value

def theory(vin, bias, I1, I2):
    return bias - 2*R2*I1*numpy.sqrt(1 + (vin/R1*I2)**2)

def general(vin, a, b, c):
    return a + b*numpy.sqrt(1 + (vin/c)**2)

tek = tektronix.tektronix()

dir = "../data/"
files = os.listdir(dir)
fnames = [file for file in files if file.startswith("R2")]

for fname in fnames:
    data = tek.load_file(dir + fname)
    
    print(fname)
    fname = fname[2:-6].split("R1")
    fname[1], amplitude = fname[1].split("V")
    R2 = get_resistor_units(fname[0])
    R1 = get_resistor_units(fname[1])
    amplitude = int(amplitude)
    
    time = data[0, 0]
    
    Vin = data[0, 1]
    Vo = data[1, 1]
    V1 = data[2, 1]
    Vout = data[3, 1]
    
    popt, pcov = scipy.optimize.curve_fit(general, Vin, Vout, [0.7, 0.05, 0.05])
    print(popt)
    print(numpy.sqrt(numpy.diag(pcov)))
    
    pyplot.title("R1:$" + str(R1) + "\Omega$, R2:$" + str(R2) + "\Omega$, V:$" + str(amplitude) + "mv$")
    pyplot.plot(time, Vin, label="Vin")
    pyplot.plot(time, Vo, label="Vo")
    pyplot.plot(time, V1, label="V1")
    pyplot.plot(time, Vout, label="Vout")
    pyplot.plot(time, general(Vin, *popt), label="Best fit")
    pyplot.legend(loc=0)
    pyplot.show()
    
