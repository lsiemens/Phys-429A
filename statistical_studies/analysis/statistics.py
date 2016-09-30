import numpy
import scipy.stats

def binning(x, min, max, width=1.0, density=False):
    bins = numpy.arange(min - width, max + 1, width) + 0.5*width
    data, _ = numpy.histogram(x, bins=bins, density=density)
    return data, numpy.arange(min, max + 1, width)

def chi_squared_statistic(observed, expected, std):
    return numpy.sum(((observed - expected)/std)**2, axis=0)

def one_parameter_fit(xaxis, data, model, pmin, pmax, steps):
    from matplotlib import pyplot
    import distributions
    samples = numpy.sum(data)
    parameter = numpy.linspace(pmin, pmax, steps)
    
    theory = samples*model(xaxis[:, numpy.newaxis], parameter)
    chi2 = chi_squared_statistic(data[:, numpy.newaxis], theory, numpy.sqrt(theory))

    min_chi2_i = numpy.argmin(chi2)
    p_min, chi2_min = parameter[min_chi2_i], chi2[min_chi2_i]
    dof = len(xaxis) - 2
    
    print("dof: " + str(dof))
    print("parameter: " + str(p_min) + ", Chi sqrd: " + str(chi2_min))
    print("probability of finding larger chi2: " + str(1 - scipy.stats.chi2.cdf(chi2_min, dof)))

    pyplot.plot(xaxis, data)
    pyplot.errorbar(xaxis, samples*model(xaxis, p_min), yerr=numpy.sqrt(samples*model(xaxis, p_min)))
    pyplot.show()

def two_parameter_fit(xaxis, data, model, amin, amax, bmin, bmax, steps):
    from matplotlib import pyplot
    import distributions
    samples = numpy.sum(data)
    a = numpy.linspace(amin, amax, steps)
    b = numpy.linspace(bmin, bmax, steps)
    
    A, B = numpy.meshgrid(a, b)
    
    theory = samples*model(xaxis[:, numpy.newaxis, numpy.newaxis], A[numpy.newaxis, :, :], B[numpy.newaxis, :, :])
    chi2 = chi_squared_statistic(data[:, numpy.newaxis, numpy.newaxis], theory, numpy.sqrt(theory))

    min_chi2_i = numpy.argmin(chi2.flat)
    a_min, b_min, chi2_min = A.flat[min_chi2_i], B.flat[min_chi2_i], chi2.flat[min_chi2_i]
    dof = len(xaxis) - 3
    
    print("dof: " + str(dof))
    print("log10(chi squared): " + str(numpy.log10(chi2_min)))
    print("a: " + str(a_min) + ", b: " + str(b_min) + ", Chi sqrd: " + str(chi2_min))
    print("probability of finding larger chi2: " + str(1 - scipy.stats.chi2.cdf(chi2_min, dof)))

#    A, B = A[::-1, :], B[::-1, :]
#    levels = numpy.power(10, numpy.linspace(0, numpy.log10(numpy.log10(chi2.max())), 50))
#    pyplot.imshow(numpy.log10(chi2), aspect="auto", extent=(a.min(), a.max(), b.min(), b.max()))
#    pyplot.contour(A, B, numpy.log10(chi2), levels, colors="k", linewidths=0.5)
#    pyplot.colorbar()
#    pyplot.show()

    pyplot.plot(xaxis, data)
    pyplot.errorbar(xaxis, samples*model(xaxis, a_min, b_min), yerr=numpy.sqrt(samples*model(xaxis, a_min, b_min)))
    pyplot.show()
