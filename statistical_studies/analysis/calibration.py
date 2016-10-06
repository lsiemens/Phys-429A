from matplotlib import pyplot
import load_data
import numpy
import scipy.stats
import distributions
import statistics

def slope(x, params):
    m, b = params
    return m*x + b

resolution = 100

guess = (0, 0)
m_bounds = (-1, 1.01)
b_bounds = (-20, 20)

data = numpy.array([0.511, 1.17, 1.274, 1.33])
yerr = numpy.array([0.0005, 0.005, 0.0005, 0.005])
x = numpy.array([337.159, 755.298, 822.709, 856.509])
xerr = numpy.array([0.017, 0.041, 0.071, 0.041])
model = slope

m, b, r, p, std = scipy.stats.linregress(x, data)
print(m, b)
mx = x.mean()
sx2 = ((x-mx)**2).sum()
std_intercept = std*numpy.sqrt(1.0/len(x) + mx**2/sx2)
std_slope = std*numpy.sqrt(1.0/sx2)
print(std_intercept, std_slope)
pyplot.errorbar(x, data, xerr=xerr, yerr=yerr)
pyplot.plot(x, slope(x, (m, b)))
pyplot.show()
#statistics.two2_parameter_fit(x, xerr, data, yerr, model, guess, [m_bounds, b_bounds])
