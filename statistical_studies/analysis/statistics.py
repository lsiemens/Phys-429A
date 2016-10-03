import numpy
import scipy.stats

def binning(x, min, max, width=1.0, density=False):
    bins = numpy.arange(min - width, max + 1, width) + 0.5*width
    data, _ = numpy.histogram(x, bins=bins, density=density)
    return data, numpy.arange(min, max + 1, width)

def CPFinv_one_sigma(CPF, guess=(1,1), step=1, res=0.01):
    # This method will break down for highly asymmetric distributions
    # where CPF(x_mean) < 0.15 or CPF(x_mean) > 0.85
    
    # for the poisson distribution the errors are within 10% of sqrt(x_mean) for mu >= 30
    # and if x_mean < 0.085033 then CPF(x_mean) < 0.15
    p_target = (1.0 - scipy.special.erf(1/numpy.sqrt(2.0)))/2.0 # =~ 0.15
    
    errors = numpy.zeros(shape=(2))
    errors[:] = numpy.float("nan")

    x = guess[0]
    p = CPF(x)
    # get x to within one "step" of p_target
    while(True):
        if p > p_target:
            x = x - step
            p = CPF(x)
            if p < p_target:
                break
        elif p < p_target:
            x = x + step
            p = CPF(x)
            if p > p_target:
                break
        else:
            break

    dx = 2*step
    # get x to within one "res" of p_target
    while(numpy.abs(p - p_target) > res):
        dx = 0.75*dx
        if p > p_target:
            x = x - dx
            p = CPF(x)
        elif p < p_target:
            x = x + dx
            p = CPF(x)
        else:
            break

    errors[0] = x
    p_target = 1.0 - p_target # =~ 0.85

    x = guess[1]
    p = CPF(x)
    # get x to within one "step" of p_target
    while(True):
        if p > p_target:
            x = x - step
            p = CPF(x)
            if p < p_target:
                break
        elif p < p_target:
            x = x + step
            p = CPF(x)
            if p > p_target:
                break
        else:
            break

    dx = 2*step
    # get x to within one "res" of p_target
    while(numpy.abs(p - p_target) > res):
        dx = 0.75*dx
        if p > p_target:
            x = x - dx
            p = CPF(x)
        elif p < p_target:
            x = x + dx
            p = CPF(x)
        else:
            break

    errors[1] = x
    return errors
    
def chi_squared_statistic(observed, expected, std):
    return numpy.sum(((observed - expected)/std)**2, axis=0)

def chi_squared(observed, expected, std):
    return ((observed - expected)/std)**2

def one_parameter_fit(xaxis, data, model, pmin, pmax, steps):
    from matplotlib import pyplot
    import distributions
    samples = numpy.sum(data)
    parameter = numpy.linspace(pmin, pmax, steps)
    chi2 = numpy.zeros(shape=parameter.shape)
    
    for i, param in enumerate(parameter):
        print(param, pmax)
        theory = samples*model(xaxis, param)
        theory_err = numpy.sqrt(samples*model(xaxis, param))
        for i in range(len(theory)):
            # if less then 30 events are predicted use full poisson errors
            if theory[i] < 30:
                errors = (theory[i] - theory_err[i], theory[i] + theory_err[i])
                errors = CPFinv_one_sigma(lambda x: distributions.poisson_CPF(x, param), guess=errors, step=theory_err[i]/3, res=0.001)
                if data[i] > theory[i]:
                    theory_err[i] = errors[1] - theory[i]
                else:
                    theory_err[i] = theory[i] - errors[0]
        chi2[i] = chi_squared_statistic(data, theory, theory_err)

    min_chi2_i = numpy.argmin(chi2)
    p_min, chi2_min = parameter[min_chi2_i], chi2[min_chi2_i]
    dof = len(xaxis) - 2
    
    print("dof: " + str(dof))
    print("parameter: " + str(p_min) + ", Chi sqrd: " + str(chi2_min))
    print("probability of finding larger chi2: " + str(1 - scipy.stats.chi2.cdf(chi2_min, dof)))
    
    pyplot.plot(parameter, chi2)
    pyplot.show()

    pyplot.plot(xaxis, chi_squared(data, samples*model(xaxis, p_min), numpy.sqrt(samples*model(xaxis, p_min))))
    pyplot.show()

    pyplot.plot(xaxis, data)
    pyplot.title("best fit by minimizing $\chi^2$")
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

    A, B = A[::-1, :], B[::-1, :]
    levels = numpy.power(10, numpy.linspace(0, numpy.log10(numpy.log10(chi2.max())), 50))
    pyplot.imshow(numpy.log10(chi2), aspect="auto", extent=(a.min(), a.max(), b.min(), b.max()))
    pyplot.title("$\log_{10}(\chi^2)$ surface")
    pyplot.contour(A, B, numpy.log10(chi2), levels, colors="k", linewidths=0.5)
    pyplot.colorbar()
    pyplot.show()

    pyplot.plot(xaxis, data)
    pyplot.title("best fit by minimizing $\chi^2$")
    pyplot.errorbar(xaxis, samples*model(xaxis, a_min, b_min), yerr=numpy.sqrt(samples*model(xaxis, a_min, b_min)))
    pyplot.show()
