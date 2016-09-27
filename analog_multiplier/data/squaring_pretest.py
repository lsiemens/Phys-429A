from matplotlib import pyplot
import tektronix
from scipy import optimize
import numpy
import os

def get_resistance(str):
    if "k" in str:
        return float(str[:-1])*1000.0
    elif "m" in str:
        return float(str[:-1])*1000000.0
    else:
        return float(str)

tek = tektronix.tektronix()

a, b, c = [], [], []
r1 = []
r2 = []
fnames = os.listdir()
fnames = [fname for fname in fnames if fname[:2] == "R2"]
for fname in fnames[2:]:
    data = tek.load_file(fname)
    fname = fname[2:]
    fname = fname.split("R1",1)
    fname = [fname[0]] + fname[1].split("V",1)
    fname[2] = fname[2][:-4]
    R2 = get_resistance(fname[0])
    R1 = get_resistance(fname[1])
    
    print("R2",R2,"R1",R1)
    
    def theory(vin, bias, I1, I2):
        return bias - 2*R2*I1*numpy.sqrt(1 + (vin/R1*I2)**2)

    def general(vin, a, b, c):
        return a + b*numpy.sqrt(1 + (vin/c)**2)
    
#    popt, pcov = optimize.curve_fit(general, data[0,1], data[3,1], [1.0, 0.01, 0.01])
    popt, pcov = optimize.curve_fit(theory, data[0,1], data[3,1], [10.0, 0.01, 0.0001])
    if popt[1] > 0.0:
        a.append(popt[0])
        b.append(popt[1])
        c.append(popt[2])
        r1.append(R1)
        r2.append(R2)
        print(pcov)

    perr = numpy.sqrt(numpy.diag(pcov))
    print(popt, perr)
    
    pyplot.plot(data[0,0], theory(data[0,1], *popt), c="k")
    pyplot.plot(data[3,0], data[3,1], c="r")
    pyplot.show()

r1 = numpy.array(r1)
r2 = numpy.array(r2)

pyplot.scatter(a,b, c="g")
pyplot.scatter(a,c, c="b")
pyplot.show()

pyplot.scatter(b,a, c="r")
pyplot.scatter(b,c, c="b")
pyplot.show()

pyplot.scatter(c,b, c="g")
pyplot.scatter(c,a, c="r")
pyplot.show()

a = sorted(a)
b = sorted(b)
#c = sorted(c)

for fname in fnames[2:]:
    data = tek.load_file(fname)

    def theory(vin, bias, Is):
        return bias - 2*R2*Is*numpy.sqrt(1 + (vin/R1*Is)**2)

    def general(vin, a, b, c):
        return a + b*numpy.sqrt(1 + (vin/c)**2)
    
    pyplot.plot(data[0,0], general(data[0,1], numpy.mean(a), numpy.mean(b), numpy.mean(c)), c="k")
    pyplot.plot(data[3,0], data[3,1], c="r")
    pyplot.show()
