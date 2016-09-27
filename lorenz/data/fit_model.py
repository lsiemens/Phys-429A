import loggerpro_csv
import numpy
from matplotlib import pyplot
import scipy.optimize as optimize

labpro = loggerpro_csv.loggerpro(4)

def fstr(value):
    return str(value) + "_rp"

def fortran_polynomial(popt):
    data = fstr(popt[0]) + "+"
    for i in range(3):
        data = data + fstr(popt[i+1]) + "*xn("+str(i+1) + ")+"
#    for i in range(3):
#        data = data + fstr(popt[i+4]) + "*xn("+str(i+1)+")**2+"
    data = data + fstr(popt[4]) + "xn(1)*xn(2)+"
    data = data + fstr(popt[5]) + "xn(1)*xn(3)+"
    data = data + fstr(popt[6]) + "xn(2)*xn(3)"
    return data

def fortran_model(fname, poptx, popty, poptz):
    data = ""
    data = data + "general_polynomial = [" + fortran_polynomial(poptx) + "," + fortran_polynomial(popty) + "," + fortran_polynomial(poptz) + "]"
    print(str(data))
    with open(fname, "w") as fout:
        fout.write(data)

def rsquared(data, model, weight=1.0):
    return 1 - numpy.sum(((data - model)*weight)**2)/numpy.sum(((data - numpy.mean(data))*weight)**2)

sets = []

fnames = [#"lorenz_s1.csv",
#          "lorenz_s2.csv",
#          "lorenz_s3.csv",
#          "lorenz_s4.csv",
#          "lorenz_s5.csv",
#          "lorenz_s6.csv",
#          "lorenz_s7.csv",
#          "lorenz_s8.csv",
#          "lorenz_s9.csv",
          "lorenz_deep.csv",
          "lorenz_deep_continued.csv"]

for fname in fnames:
    data = labpro.load(fname)
    for set in data:
        sets.append(set)
        
sets = numpy.array(sets)
poptxs = []
poptys = []
poptzs = []
varx = []
vary = []
varz = []

print(len(sets))

t = []
x = []
y = []
z = []
R = []
x_dot = []
y_dot = []
z_dot = []
for run in sets:
    t.extend(run[0])
    x.extend(run[1])
    y.extend(run[2])
    z.extend(run[3])
    R.extend(run[4])
    x_dot.extend(numpy.gradient(run[1])/numpy.gradient(run[0]))
    y_dot.extend(numpy.gradient(run[2])/numpy.gradient(run[0]))
    z_dot.extend(numpy.gradient(run[3])/numpy.gradient(run[0]))
t = numpy.array(t)
x = numpy.array(x)
y = numpy.array(y)
z = numpy.array(z)
R = numpy.array(R)
x_dot = numpy.array(x_dot)
y_dot = numpy.array(y_dot)
z_dot = numpy.array(z_dot)
print(numpy.mean(R), numpy.std(R)/numpy.sqrt(len(R)))

if True:
#for run in sets:
#    t = run[0]
#    x = run[1]
#    y = run[2]
#    z = run[3]
#    R = run[4]
    dt = numpy.gradient(t)
#    x_dot = numpy.gradient(x)/dt
#    y_dot = numpy.gradient(y)/dt
#    z_dot = numpy.gradient(z)/dt
#    pyplot.plot(x, y)
#    pyplot.show()
#    pyplot.plot(t, x_dot, t, y_dot, t, z_dot)
#    pyplot.show()

    def linear(x, m, b):
        return m*x + b

    def poly_first(t, a, b, c, d, e):
        return a + b*x + c*y + d*z + e*R
        
    def poly_first_cross(t, a, b, c, d, e, f, g, h):
        return a + b*x + c*y + d*z + e*x*y + f*x*z + g*y*z + h*R

    def poly_second(t, a, b, c, d, e, f, g, h, i, j, k):
        return a + b*x + c*y + d*z + e*x*x + f*y*y + g*z*z + h*x*y + i*x*z + j*y*z + k*R

    def iweight(radius):
        return numpy.exp((x**2 + y**2 + z**2)*radius**2)

    def weight(radius):
        return numpy.exp(-(x**2 + y**2 + z**2)*radius**2)

    iradius = 1/5.0
    model = poly_first_cross
    poptx, pcov = optimize.curve_fit(model, t, x_dot, sigma=iweight(iradius))
    lin, lincov = optimize.curve_fit(linear, x_dot, model(t, *poptx), sigma=iweight(iradius))
#    print(lin, numpy.sqrt(numpy.diag(lincov)), rsquared(x_dot, model(t, *poptx), weight(iradius)))
#    poptxs.append(poptx)
#    print(pcov)
    varx.append(numpy.diag(pcov))
    perr = numpy.sqrt(numpy.diag(pcov))
#    print("x", poptx, "\n", perr)
#    pyplot.plot(model(t, *poptx))
#    pyplot.plot(model(t, *poptx) + model(t, *perr))
#    pyplot.plot(x_dot)
    pyplot.scatter(x_dot, model(t, *poptx))
    pyplot.show()

    popty, pcov = optimize.curve_fit(model, t, y_dot, sigma=iweight(iradius))
    lin, lincov = optimize.curve_fit(linear, y_dot, model(t, *popty), sigma=iweight(iradius))
#    print(lin, numpy.sqrt(numpy.diag(lincov)), rsquared(y_dot, model(t, *popty), weight(iradius)))
#    poptys.append(popty)
#    vary.append(numpy.diag(pcov))
    perr = numpy.sqrt(numpy.diag(pcov))
    print("y", popty/11.1175, "\n", perr/11.1175)
#    print("chi:", chi_squared(y_dot, model(t, *popty), numpy.sqrt(len(t))*model(t, *perr)), "N", len(y_dot))
#    pyplot.plot(model(t, *popty))
#    pyplot.plot(y_dot)
    pyplot.scatter(y_dot, model(t, *popty))
    pyplot.show()

    poptz, pcov = optimize.curve_fit(model, t, z_dot, sigma=iweight(iradius))
    lin, lincov = optimize.curve_fit(linear, z_dot, model(t, *poptz), sigma=iweight(iradius))
#    print(lin, numpy.sqrt(numpy.diag(lincov)), rsquared(z_dot, model(t, *poptz), weight(iradius)))
#    poptzs.append(poptz)
#    varz.append(numpy.diag(pcov))
    perr = numpy.sqrt(numpy.diag(pcov))
    print("z", poptz/11.1175, "\n", perr/11.1175)
#    print("chi:", chi_squared(z_dot, model(t, *poptz), numpy.sqrt(len(t))*model(t, *perr)), "N", len(z_dot))
#    pyplot.plot(model(t, *poptz))
#    pyplot.plot(z_dot)
    pyplot.scatter(z_dot, model(t, *poptz))
    pyplot.show()

    fortran_model("./../fortran/ode_polynomial.f90", poptx/11.1175, popty/11.1175, poptz/11.1175)
#poptxs = numpy.array(poptxs)
#poptys = numpy.array(poptys)
#poptzs = numpy.array(poptzs)
#varx = numpy.sqrt(numpy.array(varx))
#vary = numpy.sqrt(numpy.array(vary))
#varz = numpy.sqrt(numpy.array(varz))

#poptx = []
#popty = []
#poptz = []
#for i in range(len(poptxs)):
#    if not ((float("Nan") in varx) or (float("Nan") in vary) or (float("Nan") in varz)):
#        print("good")
#        poptx.append(poptxs[i])
#        popty.append(poptys[i])
#        poptz.append(poptzs[i])
#    else:
#        print(varx[i], vary[i], varz[i])
#        print("bad")
#print(numpy.mean(poptxs, axis=0), numpy.std(poptxs, axis=0))
#print(numpy.mean(poptys, axis=0), numpy.std(poptys, axis=0))
#print(numpy.mean(poptzs, axis=0), numpy.std(poptzs, axis=0))
