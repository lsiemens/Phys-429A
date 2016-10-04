from matplotlib import pyplot
import load_data
import numpy
import distributions
import statistics

def poisson_gaussian(x, mu):
    return distributions.gaussian(x, mu, numpy.sqrt(mu))

def gaussian_slope(x, params):
    A, mu, sigma, m, b = params
    return A*numpy.exp(-(x - mu)**2/sigma) + m*(x - mu) + b

resolution = 100
sample_31, sample_51, sample_8 = "../data/sample_31.mca", "../data/sample_51.mca", "../data/sample_8.mca"

channel_range = [(250, 400), (750, 900), (650, 805), (800, 950), (350, 500), (6, 42)]
fnames = [sample_31, sample_31, sample_51, sample_51, sample_8, sample_8]
guess = [(21000, 340, 200, 0, 1000), (2200, 820, 300, 0, 200),
         (7000, 750, 600, -10, 1000), (6000, 850, 600, -5, 1000),
         (20000, 430, 200, 0, 300), (20000, 20, 20, 0, 3500)]
A_bounds = [(0, 40000), (0, 10000), (0, 15000), (0, 12000), (0, 40000), (0, 40000)]
mu_bounds = [(250, 400), (750, 900), (650, 805), (800, 950), (350, 500), (6, 42)]
sigma_bounds = [(0, 1000), (0, 1000), (0, 2000), (0, 2000), (0, 1000), (0, 200)]
m_bounds = [(-10, 10), (-10, 10), (-15, 15), (-15, 15), (-10, 10), (-100, 100)]
b_bounds = [(0, 10000), (0, 1000), (0, 4000), (0, 3000), (0, 2000), (0, 10000)]

for i in [5]:
    fname = fnames[i]
    raw_data = load_data.load_spectra(fname)
    pyplot.plot(raw_data)
    pyplot.show()
    xaxis = 1.0*numpy.arange(len(raw_data))
    model = gaussian_slope
    
    x, data = xaxis[channel_range[i][0]:channel_range[i][1]], raw_data[channel_range[i][0]:channel_range[i][1]]
    statistics.five_parameter_fit(x, data, model, guess[i], [A_bounds[i], mu_bounds[i], sigma_bounds[i], m_bounds[i], b_bounds[i]])
