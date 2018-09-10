# Momentum-Trading-Optimization
### Final Project for Optimization Modeling for Finance

In this project we combined simple momentum trading strategy and Markovitz portfolio optimization. On each rebalance day after we determine the list of stock to long/short, and we put this list of stock in the markovitz ooptimzation algorithm. In order to make our strategy more robust, we tested the parameters for simple momentum trading as well as the appropriate expected return and covariance matrix to be used in portfolio optimization.

To view source pdf please follow the link:
https://v2.overleaf.com/read/cqzyvnwntcyb

The code folder includes all the dependencies of codes that we used in our project

### Main files:
#### longshort.m  
It is the main matlab file to conduct all backtests, includes the zero financing portfolio rebalancing algorithm. Choose parameters and time period to conduct backtest, return portfolio performance and a plot of cumulative return.
#### muV.m  
It is the dependency file that include the algorithm for calculating different expected return and covariance matrix to be plugged in the markvovitz algorithm.
#### cvx_markowitz.m  
It is the markovitz algorithm for our portfolio optimization
### Data files:
#### ffdata_m.mat & ffdata_d.mat  
fama french factor data, monthly and daily
#### p2data.mat  
contains the stock data, 302 stocks, monthly data and daily data
#### stocklist.txt  
list of companies that meet our requirement.
