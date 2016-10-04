from matplotlib import pyplot
import load_data
import numpy
import distributions
import statistics

def poisson_gaussian(x, mu):
    return distributions.gaussian(x, mu, numpy.sqrt(mu))

resolution = 100
fnames = ["../data/sample_31.mca", "../data/sample_51.mca", "../data/sample_8.mca"]
for i in [0, 1, 2]:
    fname = fnames[i]
#    min, max = ranges[i]
#    mu_min, mu_max = mu_range[i]
    raw_data = load_data.load_spectra(fname)
    pyplot.plot(raw_data)
    pyplot.show()
#    data, bins = statistics.binning(raw_data, min, max)
#    x = numpy.linspace(min - 0.5, max + 0.5, resolution)
#    model = poisson_gaussian
#    statistics.one_parameter_fit(bins, data, model, mu_min, mu_max, 2000)
    
