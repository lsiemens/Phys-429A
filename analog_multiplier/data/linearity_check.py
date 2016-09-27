from matplotlib import pyplot
import tektronix
from scipy import optimize
import numpy
import os

tek = tektronix.tektronix()

a = []
aerr = []
fnames = os.listdir()
fnames = [fname for fname in fnames if fname[:2] == "s1"]
for fname in fnames[20:]:
    data = tek.load_file(fname)
    
    X = data[0,1]
    Y = data[1,1]
    Z = data[2,1]
    
    def F1(X, Y):
        def func(x, a):
            return a*X*Y
        return func

    def F2(X, Y):
        def func(x, a, b):
            return a*X*Y*(1 + b*(X**2 + Y**2))
        return func

    def F3(X, Y):
        def func(x, a, b, c):
            return a*X*Y*(1 + b*(X**2 + Y**2) + c*(3*X**4 + 5*X**2*Y**2 + 3*Y**4))
        return func

    function = F1

    try:
        popt, pcov = optimize.curve_fit(function(X, Y), X, Z)
        perr = numpy.sqrt(numpy.diag(pcov))
        print(popt[0], perr[0])
        a.append(popt[0])
        aerr.append(perr[0])
    
#    pyplot.plot(data[0,0], function(X, Y)(X, *popt))
#    pyplot.plot(data[0,0], Z)
#    pyplot.show()

        pyplot.scatter(Z/0.1, function(X, Y)(X, *popt)/0.1, color="k", s=2, alpha=0.05)
#        pyplot.plot(Y*X, Z)
    except:
        print("fit faild")

pyplot.title("model ($Z=Axy$) vs circuit output")
pyplot.xlabel("Circuit output in $v$")
pyplot.ylabel("model in $v$")
pyplot.show()
#print("average: ", numpy.mean(a), numpy.mean(aerr))
#pyplot.show()
