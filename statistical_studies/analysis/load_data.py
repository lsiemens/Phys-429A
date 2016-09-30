import numpy

def load_counts(fname):
    return numpy.genfromtxt(fname)

def load_spectra(fname):
    return numpy.genfromtxt(fname, skip_header=12, skip_footer=1)
    
