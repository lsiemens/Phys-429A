from matplotlib import pyplot
import load_data
import numpy
import distributions
import statistics

def poisson_gaussian(x, mu):
    return distributions.gaussian(x, mu, numpy.sqrt(mu))

def gaussian_slope(x, A, mu, sigma, m, b):
    return A*numpy.exp(-(x - mu)**2/sigma) + m*(x - mu) + b

resolution = 100
fnames = ["../data/sample_31.mca", "../data/sample_51.mca", "../data/sample_8.mca"]
for i in [0]:
    fname = fnames[i]
    raw_data = load_data.load_spectra(fname)
    pyplot.plot(raw_data)
    pyplot.show()
    x = 1.0*numpy.arange(len(raw_data))
    model = gaussian_slope
    
    x, data = x[250:400], raw_data[250:400]
    statistics.five_parameter_fit(x, data, model, (21000, 340, 200, 0, 1000), [(0, 40000), (250, 400), (0, 1000), (-10, 10), (0, 10000)])
    
