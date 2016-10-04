import numpy
import scipy.stats
import scipy.interpolate
import cPickle as pickle

def binning(x, min, max, width=1.0, density=False):
    bins = numpy.arange(min - width, max + 1, width) + 0.5*width
    data, _ = numpy.histogram(x, bins=bins, density=density)
    return data, numpy.arange(min, max + 1, width)

def precompute_CPFinv(mean_values):
# takes mean_values to precompute
# ---------------------------
#    for mean in mean_values:
#        errors = (mean - numpy.sqrt(mean), x + numpy.sqrt(mean))
#        errors = CPFinv_one_sigma(lambda x: distributions.poisson_CPF(x, mean), guess=errors, step=numpy.sqrt(mean)/3, res=1.0E-9)

    import distributions

    errors = numpy.empty(shape=mean_values.shape + (2,))
    Merr_min = mean_values - numpy.sqrt(mean_values)
    Merr_max = mean_values + numpy.sqrt(mean_values)
    step=numpy.sqrt(mean_values)/3.0
    for i in range(len(mean_values)):
        if i % int(len(mean_values) / 100) == 0:
            print(i, len(mean_values)) 
        err_guess = (Merr_min[i], Merr_max[i])
        errors[i, :] = CPFinv_one_sigma(lambda x: distributions.poisson_CPF(x, mean_values[i]), mean_values[i], guess=err_guess, step=step[i], res=1.0E-9)
    
    data = (mean_values, errors)
    
    with open("CPFinv_cache.dat", "wb") as fout:
        pickle.dump(data, fout)

def load_CPFinv_dump(zero_zero=False):
# returns interpolation functions for lower and upper bounds
    data = None
    with open("CPFinv_cache.dat", "rb") as fin:
        data = pickle.load(fin)
    if zero_zero:
        lower = scipy.interpolate.interp1d([0] + data[0].tolist(), [0] + data[1][:, 0].tolist())
        upper = scipy.interpolate.interp1d([0] + data[0].tolist(), [0] + data[1][:, 1].tolist())
    else:
        lower = scipy.interpolate.interp1d(data[0], data[1][:, 0])
        upper = scipy.interpolate.interp1d(data[0], data[1][:, 1])
    return lower, upper

def CPFinv_one_sigma(CPF, mean, guess, step=1, res=0.01):
    # This method will break down for highly asymmetric distributions
    # where CPF(x_mean) < 0.15 or CPF(x_mean) > 0.85
    
    # for the poisson distribution the errors are within 10% of sqrt(x_mean) for mu >= 30
    # and if x_mean < 0.085033 then CPF(x_mean) < 0.15
    p_int = scipy.special.erf(1/numpy.sqrt(2.0)) # =~ 0.68
    p_mid = CPF(mean)
    p_min = p_mid - p_int/2
    p_max = p_mid + p_int/2

    if p_min < 0.0:
        p_min = 0.0
        p_max = p_int
        guess = (0.0, guess[1])
    elif p_max > 1.0:
        p_max = 1.0
        p_min = 1.0 - p_int
    
    errors = numpy.zeros(shape=(2))
    errors[:] = numpy.float("nan")

    p_target = p_min
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
    p_target = p_max

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
    
def find_zero(func, minima, step=1, res=0.01):
    zeros = numpy.zeros(shape=(2))
    zeros[:] = numpy.float("nan")
    
    x = minima
    value = func(x)
    # get x to within one "step" of p_target
    while(True):
        if value < 0.0:
            x = x + step
            value = func(x)
            if value > 0:
                break
        else:
            break

    dx = 2*step
    # get x to within one "res" of p_target
    while(numpy.abs(value) > res):
        dx = 0.75*dx
        if value > 0:
            x = x - dx
            value = func(x)
        elif value < 0:
            x = x + dx
            value = func(x)
        else:
            break

    zeros[1] = x

    x = minima
    value = func(x)
    # get x to within one "step" of p_target
    while(True):
        if value < 0.0:
            x = x - step
            value = func(x)
            if value > 0:
                break
        else:
            break

    dx = 2*step
    # get x to within one "res" of p_target
    while(numpy.abs(value) > res):
        dx = 0.75*dx
        if value > 0:
            x = x + dx
            value = func(x)
        elif value < 0:
            x = x - dx
            value = func(x)
        else:
            break

    zeros[0] = x
    return zeros
    
def chi_squared_statistic(observed, expected, std):
    return numpy.sum(((observed - expected)/std)**2, axis=0)

def chi_squared(observed, expected, std):
    return ((observed - expected)/std)**2

def one_parameter_fit(xaxis, data, model, pmin, pmax, steps):
    lower, upper = load_CPFinv_dump(zero_zero=True)

    from matplotlib import pyplot
    samples = numpy.sum(data)
    parameter = numpy.linspace(pmin, pmax, steps)
    chi2 = numpy.zeros(shape=parameter.shape)
    
    theory = samples*model(xaxis[:, numpy.newaxis], parameter)
    theory_err = numpy.empty(shape=theory.shape)
    try:
        theory_err[data[:, numpy.newaxis] < theory] = lower(theory[data[:, numpy.newaxis] < theory])
        theory_err[data[:, numpy.newaxis] >= theory] = upper(theory[data[:, numpy.newaxis] > theory])

    except ValueError:
        print(theory.max())
        raise
    chi2 = chi_squared_statistic(data[:, numpy.newaxis], theory, theory_err - theory)

    min_chi2_i = numpy.argmin(chi2)
    p_min, chi2_min = parameter[min_chi2_i], chi2[min_chi2_i]
    dof = len(xaxis) - 2

    chi2 = scipy.interpolate.interp1d(parameter, chi2)
    confidence_interval = find_zero(lambda x:chi2(x) - chi2_min - 1.0, p_min, step=(pmax - pmin)/100.0, res=1.0E-5)
    print(confidence_interval)
    
    print("dof: " + str(dof))
    print("parameter: " + str(p_min) + ", Chi sqrd: " + str(chi2_min))
    print("probability of finding larger chi2: " + str(1 - scipy.stats.chi2.cdf(chi2_min, dof)))

    pyplot.plot(parameter, chi2(parameter))
    pyplot.show()

    theory = samples*model(xaxis, p_min)
    theory_err = numpy.empty(shape=theory.shape)
    theory_err[data < theory] = lower(theory[data < theory])
    theory_err[data >= theory] = upper(theory[data > theory])

    pyplot.plot(xaxis, chi_squared(data, theory, theory_err - theory))
#    pyplot.show()

    pyplot.plot(xaxis, data)
    pyplot.title("best fit by minimizing $\chi^2$")
    pyplot.errorbar(xaxis, theory, yerr=[theory - lower(theory), upper(theory) - theory])
    pyplot.show()

def two_parameter_fit(xaxis, data, model, amin, amax, bmin, bmax, steps):
    lower, upper = load_CPFinv_dump(zero_zero=True)

    from matplotlib import pyplot
    samples = numpy.sum(data)
    a = numpy.linspace(amin, amax, steps)
    b = numpy.linspace(bmin, bmax, steps)
    
    A, B = numpy.meshgrid(a, b)
    
    theory = samples*model(xaxis[:, numpy.newaxis, numpy.newaxis], A[numpy.newaxis, :, :], B[numpy.newaxis, :, :])
    theory_err = numpy.empty(shape=theory.shape)
    try:
        theory_err[data[:, numpy.newaxis, numpy.newaxis] < theory] = lower(theory[data[:, numpy.newaxis, numpy.newaxis] < theory])
        theory_err[data[:, numpy.newaxis, numpy.newaxis] >= theory] = upper(theory[data[:, numpy.newaxis, numpy.newaxis] >= theory])
    except ValueError:
        print(theory.max())
        raise
    
    chi2 = chi_squared_statistic(data[:, numpy.newaxis, numpy.newaxis], theory, theory_err - theory)

    min_chi2_i = numpy.argmin(chi2.flat)
    a_min, b_min, chi2_min = A.flat[min_chi2_i], B.flat[min_chi2_i], chi2.flat[min_chi2_i]
    dof = len(xaxis) - 3

    a_i, b_i = numpy.unravel_index(min_chi2_i, chi2.shape)
    print(A[a_i, b_i], B[a_i, b_i])

    chi2_b = scipy.interpolate.interp1d(b, chi2[:, b_i])
    chi2_a = scipy.interpolate.interp1d(a, chi2[a_i, :])
    confidence_interval_b = find_zero(lambda x:chi2_b(x) - chi2_min - 1.0, b_min, step=(bmax - bmin)/100.0, res=1.0E-5)
    confidence_interval_a = find_zero(lambda x:chi2_a(x) - chi2_min - 1.0, a_min, step=(amax - amin)/100.0, res=1.0E-5)
    print(confidence_interval_a)
    print(confidence_interval_b)

    pyplot.plot(b, chi2_b(b) - chi2_min - 1)
    pyplot.show()
    pyplot.plot(a, chi2_a(a) - chi2_min - 1)
    pyplot.show()

    
    print("dof: " + str(dof))
    print("log10(chi squared): " + str(numpy.log10(chi2_min)))
    print("a: " + str(a_min) + ", b: " + str(b_min) + ", Chi sqrd: " + str(chi2_min))
    print("probability of finding larger chi2: " + str(1 - scipy.stats.chi2.cdf(chi2_min, dof)))

    levels = numpy.power(10, numpy.linspace(-2, numpy.log10(numpy.log10(chi2.max())), 50))
    pyplot.imshow(numpy.log10(chi2), aspect="auto", extent=(a.min(), a.max(), b.max(), b.min()))
    pyplot.title("$\log_{10}(\chi^2)$ surface")
    pyplot.contour(A, B, numpy.log10(chi2), levels, colors="k", linewidths=0.5)
    pyplot.colorbar()
    pyplot.show()

    theory = samples*model(xaxis, a_min, b_min)
    theory_err = numpy.empty(shape=theory.shape)
    theory_err[data < theory] = lower(theory[data < theory])
    theory_err[data >= theory] = upper(theory[data >= theory])

    pyplot.plot(xaxis, chi_squared(data, theory, theory_err - theory))
#    pyplot.show()

    pyplot.plot(xaxis, data)
    pyplot.title("best fit by minimizing $\chi^2$")
    pyplot.errorbar(xaxis, theory, yerr=[theory - lower(theory), upper(theory) - theory])
    pyplot.show()
