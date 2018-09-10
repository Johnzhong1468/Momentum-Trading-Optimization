function [Price_all, tickers_all,date] = stockdownload(evalstart,evalend,freq)
start_date=evalstart;
end_date=evalend;
tickers_all_name = "stocklist.txt";

stocks_price_all = hist_stock_data(start_date,end_date,tickers_all_name,'frequency',freq);
%rf = hist_stock_data(start_date,end_date,'^FVX');

Price_all = [];
tickers_all = [];
date = [stocks_price_all(1).Date];

%get the close price and name of each stock
for i = 1:length(stocks_price_all)
    if length(stocks_price_all(i).AdjClose) == length(stocks_price_all(1).AdjClose)
        Price_all = [Price_all stocks_price_all(i).AdjClose];
        tickers_all = [tickers_all {stocks_price_all(i).Ticker}];
    end
end

end