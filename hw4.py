##################################################################
# Inputs
##################################################################

# E (days before expiration)
days    = int(30)
# r (annual interest rate)
rate    = float(0.01)
# S (stock price at time 0)
stock   = float(100.0)
# h0, b0, b1, b2, c
h0      = float(0.010469)
b0      = float(0.000006575)
b1      = float(0.9)
b2      = float(0.04)
ccc     = float(0.0)
# X (strike price)
strike  = float(100.0)
# n1 (number of partitions per day)
n1      = int(2)
# n2 (number of variances per node)
n2      = int(2)




import numpy as np


def get_eta(h2, rate, gamma, n1):
    eta = int(np.ceil((h2 ** .5) / gamma))
    while True:
        low_bound = abs(rate - (h2 / 2.)) / (2. * eta * gamma * (n1 ** .5))
        up_bound = min(1 - low_bound, .5)
        mid = h2 / (2. * eta * eta * gamma * gamma)

        if low_bound <= mid and mid <= up_bound:
            return eta
        elif mid < low_bound:
            print 'Valid eta does not exsist!'
            ##################################################################
            # This is important!
            ##################################################################
            print 'You motherf*cker should never get "here" :('
            print 'Sorry it was not "there" hahaha...'
            return 0

        eta += 1

def get_next_h2(b0, b1, b2, ccc, h2, rate, l, eta, gamma_n):
    eps = (l * eta * gamma_n - rate + h2 / 2.) / (h2 ** .5)
    return b0 + b1 * h2 + b2 * h2 * ((eps - ccc) ** 2.)

def get_pu_pm_pd(h2, eta, gamma, n1, rate):
    tmp1 = h2 / (eta * eta * gamma * gamma)
    tmp2 = (rate - h2 / 2.) / (2. * eta * gamma * (n1 ** .5))
    return tmp1 * .5 + tmp2, 1. - tmp1, tmp1 * .5 - tmp2


# Daily interest rate
rate /= 365
# Well, you know :)
gamma = h0
gamma_n = h0 / (n1 ** .5)
# Eta (jump parameter) tree
eta_tree = [{} for i in range(days)]
# Variance tree
h2_tree = [{} for i in range(days + 1)]
h2_tree[0][0] = [h0 ** 2 for k in range(n2)]
# Probability tree
p_tree = [{} for i in range(days)]


# Forward trees buildup using RT algorithm
for i in range(days):
    # Compute eta_tree[i] and p_tree[i] using h2_tree[i].
    for j in sorted(h2_tree[i].keys()):
        eta_tree[i][j] = [0 for k in range(n2)]
        p_tree[i][j] = [[0 for l in range(2 * n1 + 1)] for k in range(n2)]
        for k in range(n2):
            eta = get_eta(h2_tree[i][j][k], rate, gamma, n1)
            if eta == 0:
                continue

            # 
            eta_tree[i][j][k] = eta

            # (2 * n1 + 1) coefficients of (pu * x ** 2 + pm * x + pd) ** n1
            pu, pm, pd = get_pu_pm_pd(h2_tree[i][j][k], eta, gamma, n1, rate)
            poly = np.poly1d([pu, pm, pd]) ** n1
            coefs = poly.c
            for l in range(-n1, n1 + 1):
                p_tree[i][j][k][l + n1] = coefs[n1 - l]

    # Compute max./min. variances in h2_tree[i + 1].
    for j in sorted(h2_tree[i].keys()):
        for k in range(n2):
            eta = eta_tree[i][j][k]
            if eta == 0:
                continue
            h2 = h2_tree[i][j][k]

            # 
            for l in range(-n1, n1 + 1):
                next_j = j + eta * l
                next_h2 = get_next_h2(b0, b1, b2, ccc, h2, rate, l, eta, gamma_n)
                if not h2_tree[i + 1].has_key(next_j):
                    h2_tree[i + 1][next_j] = [next_h2 for k1 in range(n2)]
                else:
                    min_ = h2_tree[i + 1][next_j][0]
                    h2_tree[i + 1][next_j][0] = min(next_h2, min_)
                    max_ = h2_tree[i + 1][next_j][-1]
                    h2_tree[i + 1][next_j][-1] = max(next_h2, max_)

    # Interpolation of variances (for n2 > 2)
    for next_j in sorted(h2_tree[i + 1].keys()):
        min_ = h2_tree[i + 1][next_j][0]
        max_ = h2_tree[i + 1][next_j][-1]
        for k in range(n2):
            h2_tree[i + 1][next_j][k] = min_ + k * (max_ - min_) / (n2 - 1.)


# Pricing at the last day
put_tree = [{} for i in range(days + 1)]
for j in sorted(h2_tree[-1].keys()):
    # put = max(stock * np.exp(h0 * j) - strike, 0.)
    put = max(strike - stock * np.exp(gamma_n * j), 0.)
    put_tree[-1][j] = [put for k in range(n2)]


# Backward induction
for i in range(days - 1, -1, -1):
    for j in sorted(h2_tree[i].keys()):
        put_tree[i][j] = [put for k in range(n2)]
        for k in range(n2):
            eta = eta_tree[i][j][k]
            if eta == 0:
                continue
            h2 = h2_tree[i][j][k]

            # 
            put = 0.
            for l in range(-n1, n1 + 1):
                next_j = j + eta * l
                next_h2 = get_next_h2(b0, b1, b2, ccc, h2, rate, l, eta, gamma_n)

                # Find the next (k, k+1) interval bounding next_h2.
                for next_k in range(n2 - 1):
                    low = h2_tree[i + 1][next_j][next_k]
                    up = h2_tree[i + 1][next_j][next_k + 1]
                    if low <= next_h2 and next_h2 <= up:
                        break

                x = (next_h2 - up) / (low - up) if low - up != 0 else 0
                put_ = x * put_tree[i + 1][next_j][next_k] +\
                        (1. - x) * put_tree[i + 1][next_j][next_k + 1]

                put += p_tree[i][j][k][l + n1] * put_

            # American put pricing
            exercise = strike - stock * np.exp(gamma_n * j)
            put_tree[i][j][k] = max(put / np.exp(rate), exercise)

print 'Price: %f' % (put_tree[0][0][0])
