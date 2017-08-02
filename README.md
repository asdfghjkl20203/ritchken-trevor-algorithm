This repo contains an implementation of the Ritchken-Trevor algorithm to price American put options.
# Running the Code
[NumPy](http://www.numpy.org/) package is required. Then, just type `python hw4.py` to execute the code.
# Inputs
The following inputs are hard-coded at the beginning of `hw4.py`:
* E (days before expiration)
* r (%) (annual interest rate)
* S (stock price at time 0)
* h0, b0, b1, b2, c
* X (strike price)
* n1 (number of partitions per day)
* n2 (number of variances per node)
# Warning
If you are looking for a *reference* for your homework (probably of the course Principles of Financial Computing in NTU), be careful! The probability of escaping from the suspection is slim, if possible at all.
