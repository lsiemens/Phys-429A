from matplotlib import pyplot
import tektronix
from scipy import optimize
import numpy
import os

tek = tektronix.tektronix()

fnames = os.listdir()
fnames = [fname for fname in fnames if fname[:2] == "s2"]
for fname in fnames:
    data = tek.load_file(fname)
    
    T = data[0,0]
    X = data[0,1]
#    Y = data[1,1]
    Z = data[1,1]
    
    def F1(X, Y):
        def func(x, a):
            return a*X*Y
        return func

    def S1(X):
        def func(x, a, w, off, yof):
            return a*numpy.sin((x+off)/w)+yof
        return func

    def S2(X):
        def func(x, a, w, off, yof):
            return a*numpy.abs(numpy.sin((x+off)/w))+yof
        return func

    def F2(X, Y):
        def func(x, a, b):
            return a*X*Y*(1 + b*(X**2 + Y**2))
        return func

    def F3(X, Y):
        def func(x, a, b, c):
            return a*X*Y*(1 + b*(X**2 + Y**2) + c*(3*X**4 + 5*X**2*Y**2 + 3*Y**4))
        return func

    function = S1

    popt, pcov = optimize.curve_fit(function(T), T, X, [1, 0.005, 0, 0])
    print(popt, numpy.sqrt(numpy.diag(pcov)))

    pyplot.plot(T,X)
    pyplot.plot(T,function(T)(T, *popt))

    function = S2

    popt, pcov = optimize.curve_fit(function(T), T, Z, [0.1, 0.005, 0, 0])
    print(popt, numpy.sqrt(numpy.diag(pcov)))

    pyplot.plot(T,function(T)(T, *popt))
    pyplot.plot(T,Z)
#    pyplot.plot(function(X, X)(X, *popt), Z)

    pyplot.show()
