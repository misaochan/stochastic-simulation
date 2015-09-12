from math import *
from scipy.special import *
from numpy import *

def variance():
    """Calculates the variance of N lengths. Exactly the same as
    the function in the main code."""
    differences = []  
    for i in n_lengths:
        diff2 = (i - mean_length) **2
        differences.append(diff2)
    
    print "Sum of differences is ", sum(differences)
    s2 = 1.0/(n_count - 1) * sum(differences)  
    print "Unbiased estimate of variance, s2 = %f." % s2
    return s2

def delta_N(conf_level, s2):
    """Calculates delta(N), the half-width of the confidence interval.
    Exactly the same as the function in the main code."""
    alpha = 1 - conf_level
    p = 1 - alpha/2
    df = n_count - 1
    t = stdtrit(df,p)
    s = sqrt(s2)
    delta_N = t * s / sqrt(n_count)
    print "Half width of confidence interval, delta_N, is %f " % delta_N
    return delta_N

"""Calculates the mean length of the LS after 10 repetitions with 
different sequences of uniformly distributed numbers produced by the 
generator. Then calculates the variance, half-width of confidence 
interval, and confidence interval at confidence level 0.95."""
conf_level = 0.95
n_lengths = [14459, 14689, 15205, 14444, 13555, 14855, 14137, 15269, 14285, 14439]
n_count = 10
mean_length = float(sum(n_lengths))/float(n_count)
print "Mean length of the longest simulation run is %f" % (mean_length)
var = variance()
delta_val = delta_N(conf_level, var)
print "Confidence interval at %0.2f confidence level is (%d - %f, %d + %f)." \
      % (conf_level, mean_length, delta_val, mean_length, delta_val)
print "Confidence interval at %0.2f confidence level is (%f, %f)." \
      % (conf_level, float(mean_length)-delta_val, \
         float(mean_length)+delta_val)