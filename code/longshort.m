function  longshort
load('p2data.mat');%from 1/1/2000 to 12/31/2017 daily and monthly stock data of 302
                   %stocks, variables:
                   %date_d,date_mo,Price_d,Price_mo,tickers_d,tickers_mo
spx = hist_stock_data('01012000','12312017','^GSPC','frequency','mo');
%%%%%%%%%%%%%%%
%%%% PARAMETERS
%%%% (the following choices of parameters can easily be changed)
lambda = 1 ;%parameter for penalization
w = 0.15; %restriction for each abs(weight)
meanmode = 1; %1:simple mean, 2:101mean(1 for hi,-1 for lo,0 for others) 3: sum of means
varmode = 2;%1:simple cov, 2:fama french 3 factor, 3:single factor, 4:shrinkage model

start_date='01012007';
end_date='12312009';
decile = 10; %long top x% stocks, short top x%, enter 50 if include all stocks 
horizon = 1; % rebalance period (months)

number_of_samples = 6; %length of lookback period (month)

sample_frequency = 1; %sample every x month
						
r_w_f_o_y_e = .4; % "relative weight for one year earlier" 
				   % -- a value .4 means that for the (exponential) weights 
				    % used in computing return averages and covariances, 
				     % the weight assigned to the time period one year ago
				      % should be .4 times the weight assigned 
				       % to the most recent period.    	 
allowable_risk = 1;
	% This is the level of risk relative to the benchmark portfolio,
	% where risk is measured as standard deviation of portfolio returns.
	% Choosing this value to equal 1 means exactly the same amount of risk is allowed,
	% whereas choosing 2 means twice as much risk is allowed as the benchmark, and so on.
							
	
cumureturn = 0; % initial point

Price = Price_mo; %set data to be monthly price history
date = date_mo; %get date
[~,ind1] = min(abs(datenum(date)-datenum(start_date,'mmddyyyy'))); %get index for month of startdate
start_date = date(ind1);
[~,ind2] = min(abs(datenum(date)-datenum(end_date,'mmddyyyy'))); %index for month of enddate
end_date = date(ind2);
included_period = ind2-ind1; %months in the period

number_rebalances = round(included_period/horizon)-1;% the number of times the portfolio will be rebalanced

start = ind1;

rf = hist_stock_data(datenum(start_date),datenum(end_date),'^FVX','frequency','mo');
risk_free_rate = mean(rf.AdjClose);

tickers = tickers_mo;


d = decile;
n = length(tickers);% n = number of risky assets
q = round(n*d/100); %numeber of top/bottom stocks in the portfolio

e = ones(2*q,1);	


%get index for selected stocks               
[val, index] = sort((Price(start-1,:)-Price(start-1-number_of_samples,:))./Price(start-1-number_of_samples,:),'descend');               
lookbackreturn = (Price(start-1,:)-Price(start-1-number_of_samples,:))./Price(start-1-number_of_samples,:);
lookbackreturn(index(1:q)) = abs(lookbackreturn(index(1:q)));
lookbackreturn(index(end-q+1:end)) = -1*abs(lookbackreturn(index(end-q+1:end)));
index = index([1:q end-q+1:end]);
signal = sign(lookbackreturn(index));

bmindex = index;	                
               
               

x = 1/q*signal'; % proportions in risky assets initially



benchmark_x = x;  % proportions for benchmark portfolio
									
%%%% END OF PARAMETERS	
%%%%%%%%%%%%%%%%%%%%%%

rate_of_decay = 1 - r_w_f_o_y_e^(sample_frequency/12);
initial_wealth = cumureturn;
benchmark_wealth = cumureturn;	


rebalance_dates = start + horizon*(0:number_rebalances-1);


cumu_hist = [cumureturn];
port_return = [0];
bmwealth_hist = [cumureturn];
bm_return = [0];
for i = 1:length(rebalance_dates)

	trade_date = rebalance_dates(i);
		
    %%%%% REBALANCE PORTFOLIO AND PAY TRANSACTION COSTS %%%%%%
    [mu,V] = muV(Price(:,index),trade_date, ...
    horizon,sample_frequency,number_of_samples,rate_of_decay,meanmode,varmode,decile,index);
    
    benchmark_risk = sqrt(e'*V*e)/length(index); 		
    sigma = allowable_risk*benchmark_risk; 
 
    [x] =  cvx_markowitz(mu,V,sigma,lambda,w);
    	
    %%%%%% PROCEED TO END OF TIME PERIOD AND ACCOUNT FOR GAINS, LOSSES %%%%%	
    
    %%%find everyday progress of the portfolio%%%
    x_index = index;
    x_returns = zeros(n,1); % x = all assets(index)
    x_returns(x_index) = x;
    for j = 1:horizon  
        dailyreturns = Price(trade_date+j,:)./Price(trade_date+j-1,:);
        x_returns(x_index) = dailyreturns(x_index)'.*x_returns(x_index);
        port_return = [port_return sum(x_returns(x_index))];
        cumu_hist = [cumu_hist (1+cumureturn)*(1+port_return(end))-1];        
    end    
    cumureturn = cumu_hist(end)
      
    %%%benchmark progression, equally weighted portfolio%%% 	  
    bmreturns = zeros(n,1);
    bmreturns(bmindex) = benchmark_x;
        for j = 1:horizon  
            dailyreturns = Price(trade_date+j,:)./Price(trade_date+j-1,:);
            bmreturns(bmindex) = dailyreturns(bmindex)'.*bmreturns(bmindex);
            bm_return = [bm_return sum((bmreturns(bmindex)))];
            bmwealth_hist = [bmwealth_hist (1+benchmark_wealth)*(1+bm_return(end))-1];
        end  
    benchmark_wealth = bmwealth_hist(end)
    %bm_delta = abs(benchmark_x0/(1+benchmark_x0+sum(benchmark_x))-1/(q+1))+sum(abs(benchmark_x/(1+benchmark_x0+sum(benchmark_x))-((1-benchmark_x0)/q)*ones(q,1)));
    benchmark_x = 1/q*signal';;
    	% At the beginning of the time period, the benchmark portfolio was
    	% rebalanced 
    

    
if i < length(rebalance_dates)   
%reissue the index before the last rebalance date
[~, index] = sort((Price(rebalance_dates(i+1)-1,:)-Price(rebalance_dates(i+1)-number_of_samples-1,:))./Price(rebalance_dates(i+1)-number_of_samples-1,:),'descend');               
lookbackreturn = (Price(rebalance_dates(i+1)-1,:)-Price(rebalance_dates(i+1)-number_of_samples-1,:))./Price(rebalance_dates(i+1)-number_of_samples-1,:);
lookbackreturn(index(1:q)) = abs(lookbackreturn(index(1:q)));
lookbackreturn(index(end-q+1:end)) = -1*abs(lookbackreturn(index(end-q+1:end)));
index = index([1:q end-q+1:end]);
signal = sign(lookbackreturn(index));
               

x = 1/q*signal'; % proportions in risky assets initially
    
bmindex = index;

end




    
       
    
%%% sp500 as another benchmark         
spx_ret =  [spx.AdjClose(start+1:trade_date+horizon)./spx.AdjClose(start:trade_date+horizon-1)-1];      
spx_hist = [0];
for i = 1:length(spx_ret)
    spx_hist = [spx_hist (1+spx_ret(i))*(1+spx_hist(end))-1];
end
      
end

fprintf('portfolio final wealth %f\n',cumureturn);
fprintf('Simple Momentum final wealth %f\n',benchmark_wealth);
fprintf('sp500 final wealth %f\n',spx.AdjClose(rebalance_dates(end)+horizon)/spx.AdjClose(start)-1);

fprintf('portfolio annualized rate of return %f\n', (cumureturn)^(12/(horizon*number_rebalances))-1);
fprintf('Simple Momentum annualized rate of return %f\n' , (benchmark_wealth)^(12/(horizon*number_rebalances))-1);
fprintf('sp500 annualized rate of return %f\n' , (spx.AdjClose(rebalance_dates(end)+horizon)/spx.AdjClose(start))^(12/(horizon*number_rebalances))-1);


figure()
plot(datetime(date_mo(ind1:ind2-1)),cumu_hist)
hold on
plot(datetime(date_mo(ind1:ind2-1)),bmwealth_hist)
plot(datetime(date_mo(ind1:ind2-1)),spx_hist)
legend('Portfolio return','Simple Momentum return','S&P500 return')
hold off

exprtn = mean(port_return)-risk_free_rate/100/12
risk = std(port_return)
portfolio_sharpe = exprtn/risk*sqrt(12) %sharpe ratio
fama3 = timeseries(cumu_hist,date_mo(ind1:ind2-1),'Name','FamaFrenchCov lambda = 3');
bm = timeseries(bmwealth_hist,date_mo(ind1:ind2-1),'Name','Equally weighted BM');

