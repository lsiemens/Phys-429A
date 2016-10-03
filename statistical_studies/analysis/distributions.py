import numpy
import scipy.misc
import scipy.special

#def poisson(x, mu):
#    return numpy.exp(-mu)*numpy.power(mu, x)/scipy.misc.factorial(x)

def poisson(x, mu):
    # this expression for the poisson distribution is equivelent to e^(-mu)*mu^(k)/(k)!
    # but it is numericaly even when k is large. Note scipy.special.gammaln(x + 1) = ln(Gamma(x + 1)) = ln(x!)
    return numpy.exp(x*numpy.log(mu) - mu - scipy.special.gammaln(x + 1.0))

def poisson_CPF(x, mu):
    # scipy.special.gammainc(x, mu) = GammaInc(x, mu)/Gamma(x) where GammaInc(x, mu) is the incomplete gamma function
    # note for integer x where x > 0 and mu >= 0 then GammaInc(x, mu) = (x - 1)!*SUM_{k=0}^{x - 1}e^(-mu)*mu^k/k!
        
    value = 1.0 - scipy.special.gammainc(x, mu)
    try:
        value[x <= 0] = 0.0
    except TypeError:
        if x <= 0:
            value = 0.0
    return value

def gaussian(x, mu, sigma):
    return numpy.exp(-(x - mu)**2/(2*sigma**2))/numpy.sqrt(2*numpy.pi*sigma**2)
