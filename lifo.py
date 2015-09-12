"""Sequential simulation model for LIFO at confidence level 0.9 and max error 0.05"""
"""#
# LastInFirstOut simulation
#
# r[i] = (107374182 * r[i-1] + 104480 * r[i-5]) mod (2**31 - 1)
#
"""
from SimPy.Simulation import *
from SimPy.SimPlot import *
from math import *
from scipy.special import *
from numpy import *

##Random number generator. 

class RandomGenerator(object):
    """Generates uniformly-distributed random numbers that are transformed to the appropriate CDF that we require in exponential_transform() and pareto_transform(), which are called from the main simulation code to produce interarrival times and printing times respectively. """
    
    A = 107374182
    B = 104480
    M = 2**31 - 1
    state = [4,3,2,1,0]
    
    def random(self):
        """Returns the next number in the sequence as an integer in the range 
        0 to 2**31 - 2. This method is called each time a realization of 
        interarrival times or printing times is required in the 
        simulation."""
        r = (self.A * self.state[-1] + self.B * self.state[-5]) % self.M
        del self.state[0]
        self.state.append(r)
        return r
    
    def seed(self, values):
        """Initialises the state of the prng from a sequence of 5 integers. 
        Set for each simulation run with last five numbers from prev run."""
        if len(values) <> 5:
            raise ValueError("Seed must have length 5")
        for i in values:
            if not isinstance(i, int):
                raise TypeError("Seed values are not all integers")
        self.state[:] = values
    
    def exponential_transform(self):
        """Calls random() to generate a random number from the uniform 
        distribution in the range (0,1). Then transforms the number to obtain 
        a realization of an exponential random variable."""
        number = float(float(self.random())/(2**31-2))
        expo_number = -1.32 * log(1 - number) 
        return expo_number
    
    def pareto_transform(self):
        """Calls random() to generate a custom number in the uniform 
        distribution range (0,1). Then transforms the number to obtain a 
        realization of a Pareto random variable."""
        number = float(float(self.random())/(2**31-2)) 
        number2 = pow((1 - number), 1/2.025)
        pareto_number = 1/number2
        return pareto_number
    

## Model components

class Source(Process):
    """Inherits from the SimPy.Simulation Process class. 
    Creates PrintJob objects to be used in the simulation.
    Each PrintJob generated operates independently of the others."""

    def generate(self, resource, monitor):
        i = 0
        # Creates a sequence of numbered PrintJobs until 
        # stopped (when max simulation time is reached)
        while(1):
            i = i + 1
            c = PrintJob(name="PrintJob%02d" % (i))
            # PrintJob is activated at the current simulation time 
            # and submits a printing task to 'resource'. 
            # This is monitored to record the data, and priority 
            # increases with each PrintJob generated because in LIFO,
            # the most recent job submitted is given priority.
            activate(c, c.visit(res=resource, moni=monitor, priority=i))
            ta = random_gen.exponential_transform()
            # Pauses for simulation time 'ta', 
            # which is a realization of an exponential RV, before
            # generating the next Print Job
            yield hold, self, ta


class PrintJob(Process):
    """Inherits from the SimPy.Simulation Process Class. 
    Models a PrintJob that is submitted to the buffer,
    is placed in queue (if printers are not free), and 
    exits the system when it is completed."""

    def visit(self, res, moni, priority):
        """Models the activity of a PrintJob that is generated 
        by the Source()."""
        # Records the time of submission of printing task
        arrive = now()

        # Requests for an available Resource (printer). If the printer 
        # is free then the PrintJob can start service immediately and 
        # the program execution moves on to the next line. 
        # If the printer is busy, the PrintJob is automatically queued 
        # by the Resource according to his priority.
        # Print jobs with higher priority are serviced first.
        tp = random_gen.pareto_transform()
        yield request, self, res, priority

        #Operation of the printer is modelled by the yield hold statement. 
        #During the simulated time period of tp, the print job is being served.
        yield hold, self, tp
        wait = now() - arrive
        
        #The Monitor observes and records the waiting time for each print 
        #job for use in statistical analyses.
        moni.observe(wait)  

        #The current print job completes service and the printer becomes 
        #available for any remaining jobs in the queue.
        yield release, self, res


## Main simulation components

maxTime = 600.0    #In terms of minutes. Simulate for ten hours
random_gen = RandomGenerator()
# PRNG is initialized outside of replication loop so the state of the
# generator is saved over replications and the last 5 numbers from one
# replication is used as the seed for the next.
random_gen.seed((371949050, 1907301114, 1662645853, 2070256321, 1977108697)) 
d_max = 0.01
n_means = []
n_count = 0
conf_level = 0.95

def model():
    """Creates the objects needed for one simulation run and initializes 
    its state. Generates one replication of the simulation each time it 
    is called."""
    # Creates a Resource object that models two printers with a common 
    # buffer. The PriorityQ qType allows us to set the priorities for each 
    # print job, which is done during the generation of
    # print jobs in Source. Higher priority jobs are served first. When 
    # both printers are available, Resource() inherently chooses one printer 
    # at random so there is no need to add extra code for that.
    k = Resource(capacity=2, name="Printer", qType=PriorityQ)
    
    # Creates a Monitor object to record the data obtained
    wM = Monitor()                                     

    # Sets up the simulation system ready to receive activate calls
    initialize()
    
    # Creates a Source object and activates it, specifying as parameters the
    # object to be activated, the call of the action routine(s.generate())
    # and that it is to be activated at simulation time 0.0
    s = Source('Source')
    activate(s, s.generate(resource=k, monitor=wM), at=0.0)
    
    # Starts the simulation, which runs until the simulation time is maxTime
    simulate(until=maxTime)
    
    # Returns the number of print jobs served in one simulation run 
    # and the mean of the waiting time for jobs.
    return (wM.count(), wM.mean())                    

def variance():
    """Calculates the variance of mean waiting times obtained during the
    past N replications."""
    differences = []  
    for i in n_means:
        diff2 = (i - total_mean)**2
        differences.append(diff2)
    
    print "Sum of differences is ", sum(differences)
    s2 = 1.0/(n_count - 1) * sum(differences)
    
    print "Variance of N means is %f." % s2
    return s2

def delta_N(conf_level, s2):
    """Calculates the half-width of the confidence interval at
    a given confidence level, after knowing the variance."""
    alpha = 1 - conf_level
    p = 1 - alpha/2
    df = n_count - 1
    t = stdtrit(df,p)
    
    s = sqrt(s2)
    delta_N = t * s / sqrt(n_count)
    print "Half width of confidence interval, delta_N, is %f " % delta_N

    return delta_N

##Main simulation loop

print "Running nmin = 3 replications."

# Generates replications of the model with non-overlapping sequences of
# numbers generated by the PRNG. Runs(almost) indefinitely until stopping
# condition is fulfilled
for i in range (1, 9999999):
    result = model()
    (count, mean) = result
    print "Replication #%d. MeanWaitingTime(i) for %d printing \
    tasks was %f min." % (i, count, mean)
    # Saves a list of mean waiting times obtained by current and
    # previous replications 
    n_means.append(mean)
    n_count = n_count + 1
    
    # Executes at least nmin replications, then calculates the 
    # current value of relative statistical error d(n)
    if i >= 3:
        total_mean = sum(n_means)/n_count
        print "MeanWaitingTime(n) for %d replications is %f" % \
              (n_count, total_mean)
        var = variance()
        delta_val = delta_N(conf_level, var)
        d_N = delta_val / total_mean
        print "Statistical error for %d replications is %f . \n" % \
              (n_count, d_N)
    
        # If the relative error d(n) is not greater than the max relative 
        # error allowed, stop the simulation. Otherwise go to (n+1)st 
        # replication
        if d_N <= d_max:
            print "Simulation completed with d_N = %f . Number of \
            replications required: %d. Final MeanWaitingTime(%d) \
            is %f" % (d_N, i, i, total_mean)
            break

# Print out last 5 numbers generated by PRNG to use as seed in next run
print "Last five pseudorandom numbers generated from previous run, \
in order of generation", random_gen.state
