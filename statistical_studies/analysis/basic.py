from matplotlib import pyplot
import load_data
import numpy
import distributions
import statistics

max = 10
resolution = 100
fnames = ["../data/Run_1.dat", "../data/Run_2.dat", "../data/Run_3.dat", "../data/Run_4.dat"]
ranges = [(90, 200), (0, 15), (0, 15), (50, 150)]
for i in [0, 1, 2, 3]:
    fname = fnames[i]
    min, max = ranges[i]
    raw_data = load_data.load_counts(fname)
    data, bins = statistics.binning(raw_data, min, max)
    samples = numpy.sum(data)
    x = numpy.linspace(min - 0.5, max + 0.5, resolution)
    
    mu = numpy.linspace(min + 1, max, 10000)
    chi = []
    for m in mu:
        #theory = samples*distributions.poisson(bins, m)
        theory = samples*distributions.gaussian(bins, m, numpy.sqrt(m))
        chi.append(statistics.chi_squared_statistic(data, theory, numpy.sqrt(theory)))
    mu_best = mu[chi.index(numpy.array(chi).min())]
    pyplot.plot(mu, chi)
    pyplot.show()
    pyplot.plot(mu, numpy.log(chi))
    pyplot.show()
    
    pyplot.plot(bins, data)
    pyplot.plot(x, samples*distributions.gaussian(x, mu_best, numpy.sqrt(mu_best)))
    #pyplot.plot(bins, samples*distributions.poisson(bins, mu_best))
    pyplot.show()
