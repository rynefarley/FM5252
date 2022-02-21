import numpy as np
import sys
from scipy.stats import norm
sys.setrecursionlimit(10**7)

# Initialise parameters
S0 = 50      # initial stock price
K = 50       # strike price
T = 1         # time to maturity in years
r = 0.06      # annual risk-free rate
V = .5      # Volatility
N = 10000       # number of time steps
q = .01        # dividend



def binomial_tree(S0,K,r,V,T,N,q):
    # Allow more rescursion than standard
    sys.setrecursionlimit(10**7)
    # Set constants
    t = T/N
    u = np.exp(V*np.sqrt(t))
    d = 1/u
    p = (np.exp((r-q)*t) - d) / (u-d)
    disc = np.exp(-r*t)
    #adjust for dividends
    r= r-q
    #setting up the tree
    tree = S0 * d ** (np.arange(N,-1,-1)) * u ** (np.arange(0,N+1,1))

    # European options
    def European(S0,K,r,V,T,N,q,tree):    
        EC = np.maximum( tree - K , np.zeros(N+1) )
        EP = np.maximum(K - tree, np.zeros(N+1))
        i = N
       # iterate backwards through tree recursively
        def NotALoop(S0,K,r,V,T,N,q,i,tree,EC,EP):
            if i > 0:
                EC = disc * ( p * EC[1:i+1] + (1-p) * EC[0:i] )
                EP = disc * ( p * EP[1:i+1] + (1-p) * EP[0:i] )
                NotALoop(S0,K,r,V,T,N,q,i-1,tree,EC,EP)
            else:
               print('Europoean Call:',EC,'European Put:',EP)
        NotALoop(S0,K,r,V,T,N,q,i,tree,EC,EP)
    European(S0,K,r,V,T,N,q,tree)

    # American options
    def American(S0,K,r,V,T,N,q,tree):
        AP = np.maximum(0, K - tree)
        AC = np.maximum(0, tree - K)
        j = N - 1
        # iterating through tree recursively
        def NotALoop2(j,AC,AP,p,disc,u,d,S0,K):
            if j  > -1:
                S = S0 * d**(np.arange(j,-1,-1)) * u**(np.arange(0,j+1,1))
                AC [:j+1] = disc * ( p*AC[1:j+2] + (1-p)*AC[0:j+1] )
                AC = AC[:-1]
                AC = np.maximum(AC, S - K)
                AP[:j+1] = disc * ( p*AP[1:j+2] + (1-p)*AP[0:j+1] )
                AP = AP[:-1]
                AP = np.maximum(AP,K - S)
                NotALoop2(j-1,AC,AP,p,disc,u,d,S0,K)
            else:
                # printing results
                print('American Call:',AC,'American Put:',AP)
        NotALoop2(j,AC,AP,p,disc,u,d,S0,K)                        
    American(S0,K,r,V,T,N,q,tree)
    
    # d1 equation
    d1=(np.log(S0/K)+(r +(V ** 2)/2)*T)/(V *np.sqrt(T))
    # d2 equation
    d2 = d1 - V*np.sqrt(T)
    # call delta
    calldelta = norm.cdf(d1, 0, 1)
    # put delta
    putdelta = -norm.cdf(-d1, 0, 1)
    # gamma
    gamma = norm.pdf(d1, 0, 1)/(S0*V*np.sqrt(T))
    # vega
    vega  = S0*norm.pdf(d1, 0, 1)*np.sqrt(T)
    # call theta
    calltheta = -S0*norm.pdf(d1, 0, 1)*V/(2*np.sqrt(T)) - r*K*np.exp(-r*T)*norm.cdf(d2, 0, 1)
    # put theta
    puttheta = -S0*norm.pdf(d1, 0, 1)*V/(2*np.sqrt(T)) + r*K*np.exp(-r*T)*norm.cdf(-d2, 0, 1)
    # call rho
    callrho = K*T*np.exp(-r*T)*norm.cdf(d2, 0, 1)
    # put rho
    putrho = -K*T*np.exp(-r*T)*norm.cdf(-d2, 0, 1)
    
    print('Call Delta:',calldelta,'Put Delta:',putdelta,'Gamma:',gamma,'Vega:',vega/100)
    print('Call Theta:',calltheta/100,'Put Theta:',puttheta/100,'Call Rho',callrho/100,'Put Rho:',putrho/100)

# ASX Portfolio was a great resource for this assignment
binomial_tree(S0,K,r,V,T,N,q)

