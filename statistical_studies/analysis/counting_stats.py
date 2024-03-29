from matplotlib import pyplot
import load_data
import numpy
import distributions
import statistics

def poisson_gaussian(x, mu):
    return distributions.gaussian(x, mu, numpy.sqrt(mu))

resolution = 100
fnames = ["../data/Run_1.dat", "../data/Run_2.dat", "../data/Run_3.dat", "../data/Run_4.dat"]
ranges = [(102, 185), (0, 12), (0, 12), (70, 127)]
mu_range = [(140, 150), (0.5, 12), (0.5, 12), (90, 100)]
sigma_range = [(10, 16), (0.5, 4), (0.5, 4), (10.5, 11)]
print("Gaussian constrained Fit")
for i in []:
    fname = fnames[i]
    min, max = ranges[i]
    mu_min, mu_max = mu_range[i]
    raw_data = load_data.load_counts(fname)
    data, bins = statistics.binning(raw_data, min, max)
    x = numpy.linspace(min - 0.5, max + 0.5, resolution)
    
    model = poisson_gaussian

    statistics.one_parameter_fit(bins, data, model, mu_min, mu_max, 2000)
    
print("Poisson Fit")
for i in []:
    fname = fnames[i]
    min, max = ranges[i]
    mu_min, mu_max = mu_range[i]
    raw_data = load_data.load_counts(fname)
    data, bins = statistics.binning(raw_data, min, max)
    x = numpy.linspace(min - 0.5, max + 0.5, resolution)
    
    model = distributions.poisson

    statistics.one_parameter_fit(bins, data, model, mu_min, mu_max, 2000)

print("Gaussian Fit")
for i in [3]:
    fname = fnames[i]
    min, max = ranges[i]
    mu_min, mu_max = mu_range[i]
    sigma_min, sigma_max = sigma_range[i]
    raw_data = load_data.load_counts(fname)
    data, bins = statistics.binning(raw_data, min, max)
    x = numpy.linspace(min - 0.5, max + 0.5, resolution)
    
    model = distributions.gaussian

    statistics.two_parameter_fit(bins, data, model, mu_min, mu_max, sigma_min, sigma_max, 200)

ranges = [(112, 175), (0, 12), (0, 12), (70, 122)]
mu_range = [(140, 150), (0.5, 12), (0.5, 12), (70, 127)]
sigma_range = [(10, 14), (0.5, 4), (0.5, 4), (8, 11)]
