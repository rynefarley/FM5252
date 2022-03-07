import numpy as np
from scipy.stats import norm
N = norm.cdf
import pandas as pd
!pip install yfinance
import yfinance as yf
import datetime
import scipy.optimize
import matplotlib.pyplot as plt


def BS_CALL(S, K, T, r, sigma):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return S * N(d1) - K * np.exp(-r*T)* N(d2)

def BS_PUT(S, K, T, r, sigma):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma* np.sqrt(T)
    return K*np.exp(-r*T)*N(-d2) - S*N(-d1)

def bisection(S,K,r,t,price,otype):
    sigma = .5
    lowerbound = 0 
    upperbound = 10
    imp  = -9.999
    a = False
    
    while abs(price - imp) > .00001:
        if a == True:
            if imp > price:
                upperbound = sigma
                sigma = (sigma + lowerbound) / 2
            else:
                lowerbound = sigma
                sigma = (sigma + upperbound) / 2
        
        a = True
        
        if otype == 'c':
            imp = BS_CALL(S,K, t, r, sigma)
        elif otype == 'p':
            imp = BS_PUT(S, K, t, r, sigma)
        else:
            return print("Please enter 'c' for call or 'p' for put")
            
    print('Implied vol: ',sigma)


def vega(r, S, K, T, sigma):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    vega = S*norm.pdf(d1, 0, 1)*np.sqrt(T)
    return vega

def Newton(S,K,r,t,price,otype):
    sigma = .5
    imp = -9.999
    a = False
    
    while abs(price - imp) > .001:
        if a == True:
            if vega(r,S,K,t,sigma) != 0:
                sigma = sigma + ((price - imp) / vega(r,S,K,t,sigma))
            else:
                return sigma
        a = True
        
        if otype == 'c':
            imp = BS_CALL(S,K, t, r, sigma)
        elif otype == 'p':
            imp = BS_PUT(S, K, t, r, sigma)
        else:
            return print("Please enter 'c' for call or 'p' for put")
        
    return sigma

tk = yf.Ticker('AAPL')
exps = tk.options

options = pd.DataFrame()
for e in exps:
    opt = tk.option_chain(e)
    opt = pd.DataFrame().append(opt.calls).append(opt.puts)
    opt['expirationDate'] = e
    options = options.append(opt, ignore_index=True)

options['expirationDate'] = pd.to_datetime(options['expirationDate']) + datetime.timedelta(days = 1)
options['dte'] = (options['expirationDate'] - datetime.datetime.today()).dt.days / 365    
options['CALL'] = options['contractSymbol'].str[4:].apply(lambda x: "C" in x) 
options[['bid', 'ask', 'strike']] = options[['bid', 'ask', 'strike']].apply(pd.to_numeric)
options['mark'] = (options['bid'] + options['ask']) / 2 # Calculate the midpoint of the bid-ask
options = options.drop(columns = ['contractSize', 'currency', 'change', 'percentChange', 'lastTradeDate', 'lastPrice'])
options = options.loc[(options['expirationDate']) == options['expirationDate'].median()]
options = options.loc[(options['CALL']) == True] 
    
ticker = 'AAPL'
ticker_yahoo = yf.Ticker(ticker)
data = ticker_yahoo.history()
price = (data.tail(1)['Close'].iloc[0])
impliedVol = np.arange(len(options))
options['impliedVol'] = 0
for row in options.iterrows():
    impliedVol = Newton(price,options['strike'].astype(float)[row[0]],.03,1,
                                   ((options['ask'].astype(float)[row[0]]) + options['bid'].astype(float)[row[0]]) / 2 ,'c')
    plt.scatter(options['strike'][row[0]],impliedVol)
    
options['impliedVol']
plt.ylim(-.5, 1.5)
plt.show

S,K,r,t = 50,50,.05,1
print(bisection(S,K,r,t,9.98,'c'))
print(Newton(S,K,r,t,9.98,'c'))

