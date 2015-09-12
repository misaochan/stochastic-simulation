# stochastic-simulation


A LIFO (last in first out) stochastic simulation model using the SimPy library, to estimate the mean waiting time for printing tasks during 10 hours of simulated time.

**Hypothetical task:** Evaluate options that exist for  upgrading  printing facilities in a computer laboratory at a university. Administration is interested how to install two identical laser printers for a class  of 1200 students, in a laboratory equipped with 600 workstations linked by a local area network.  The plan is that jobs submitted for printing will be spooled in one common buffer 1  of these two  printers, but the order in which jobs stored in that buffer will be selected for printing has not been decided yet. 

With this framework in place, alternative models to LIFO can be simulated by altering the `def visit(self, res, moni, priority)` method. 
