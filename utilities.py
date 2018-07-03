import numpy as np

def normPDF(x, mean, sigma):
    return exp(-1.*(x-mean)**2/(2.*sigma**2))
    # return (1/np.sqrt(2.*np.pi*sigma**2))*exp(-1.*(x-mean)**2/(2.*sigma**2))

def normSigma(x, mean, probability):
    return np.sqrt(-0.5*(x-mean)**2/np.log(probability))

def normInverse(p, mu, sigma):
    """
    return the offset to the Gaussian center given the Gaussian informaton

    Keyword arguments:
    p -- the Gaussian probabilty
    mu -- the center of the Gaussian
    sigma -- the deviation  of the Gaussian

    return:
    offset -- the offset to the Gaussian center given the parameters

    """

    x = np.sqrt(-2.*sigma**2.*np.log(p))+mu
    return x

