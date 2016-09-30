import numpy
import scipy.misc

def poisson(x, mu):
    return numpy.exp(-mu)*numpy.power(mu, x)/scipy.misc.factorial(x)

def gaussian(x, mu, sigma):
    return numpy.exp(-(x - mu)**2/(2*sigma**2))/numpy.sqrt(2*numpy.pi*sigma**2)
