##QNC HW 2
#import libraries
import numpy as np 
import random as rnd
import collections
import matplotlib.pyplot as plt
import time
import scipy.stats as st
from scipy.stats import bernoulli, binom, poisson, chi2
#from IPython.display import clear_output 
from operator import itemgetter
from statsmodels.stats import proportion
from numpy import matlib

#EXCERCISE 1
n = 10   # total number of available quanta
p = 0.2  # probability a single quantum is released

for x in range(n + 1):  # We loop over all possible outcomes: 0 quanta, 1 quantum, 2 quanta, … up to all 10 quanta, x in our case ranges from 0 to 10
    prob = st.binom.pmf(x, n, p)
    print(f"P({x} quanta released) = {prob:.5f}")


#EXCERCISE 2
n = 14      # number of quanta available
x = 8       # observed quanta released

# (1) & (2): check p=0.1 and p=0.7
p01 = st.binom.pmf(x, n, 0.1)
p07 = st.binom.pmf(x, n, 0.7)
print(f"P(X=8 | p=0.1) = {p01:.8f}")
print(f"P(X=8 | p=0.7) = {p07:.8f}")

# (3): the 
probs = [] #emplty list to store (p, likelihood) 
for p in [i/10 for i in range(1, 11)]:  # Loops through all the probabilities from 0.1 to 1.0 in steps of 0.1
    pr = st.binom.pmf(x, n, p) #Calls our binomial probability function
    probs.append((p, pr)) #putting back into the list
    print(f"P(X=8 | p={p:.1f}) = {pr:.8f}")

# (4): find which decile is most likely
best_p, best_val = max(probs, key=lambda t: t[1])
print(f"Most likely decile p given X=8: p={best_p:.1f} (likelihood={best_val:.8f})")


#EXERCISE 3
#You now have two independent experiments:
#Trial 1: observed 8/14 released
#Trial 2: observed 5/14 released
#You want to know: If the true release probability was p=0.1 how likely are these results together?
n = 14
p = 0.1   # assumed release probability
x1, x2 = 8, 5  # your two observations

# individual likelihoods
L1 = st.binom.pmf(x1, n, p)
L2 = st.binom.pmf(x2, n, p)

# total likelihood = product
L_total = L1 * L2

# total log-likelihood = sum of logs
logL_total = np.log(L1) + np.log(L2)

print(f"Likelihood for X=8: {L1:.8e}")
print(f"Likelihood for X=5: {L2:.8e}")
print(f"Total likelihood:   {L_total:.8e}")
print(f"Total log-likelihood: {logL_total:.8f}")

#EXERCISE 4
counts = {0:0,1:0,2:3,3:7,4:10,5:19,6:26,7:16,8:16,9:5,10:5,11:0,12:0,13:0,14:0}
n = 14

total_experiments = sum(counts.values())
total_releases = sum(k * c for k, c in counts.items())

p_hat = total_releases / (n * total_experiments)
print(f"Total experiments: {total_experiments}")
print(f"Total releases: {total_releases}")
print(f"p-hat = {p_hat:.2f}")


