start_date='01012000';
end_date='01012018';

[Price_mo, tickers_mo,date_mo] = stockdownload(start_date,end_date,'mo');
[Price_d, tickers_d,date_d] = stockdownload(start_date,end_date,'d');

save('p2data.mat','Price_mo', 'tickers_mo','date_mo','Price_d', 'tickers_d','date_d')