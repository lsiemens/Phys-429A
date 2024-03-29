from matplotlib import pyplot
import load_data
import numpy
import distributions
import statistics

def poisson_gaussian(x, mu):
    return distributions.gaussian(x, mu, numpy.sqrt(mu))

max = 10
resolution = 100
fnames = ["../data/Run_1.dat", "../data/Run_2.dat", "../data/Run_3.dat", "../data/Run_4.dat"]
ranges = [(100, 187), (0, 13), (0, 13), (68, 129)]
mu_range = [(100, 187), (0.5, 13), (0.5, 13), (68, 129)]
sigma_range = [(10, 14), (0.5, 4), (0.5, 4), (8, 11)]
print("Gaussian Fit")
for i in [0, 1, 2, 3]:
    fname = fnames[i]
    min, max = ranges[i]
    mu_min, mu_max = mu_range[i]
    sigma_min, sigma_max = sigma_range[i]
    raw_data = load_data.load_counts(fname)
    data, bins = statistics.binning(raw_data, min, max)
    x = numpy.linspace(min - 0.5, max + 0.5, resolution)
    
    model = distributions.gaussian

    statistics.two_parameter_fit(bins, data, model, mu_min, mu_max, sigma_min, sigma_max, 200)
