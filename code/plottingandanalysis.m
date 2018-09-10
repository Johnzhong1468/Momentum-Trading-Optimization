figure()
plot(timeseries(mean(Price_d./mean(Price_mo),2),date_d,'Name','Equally Weighted Mean Stock Price'));
hold on
[val, index] = sort(mean(diff(Price_mo)));
plot(timeseries(mean(Price_d(:,index(1:100))./mean(Price_d(:,index(1:100))),2),date_d,'Name','Top 100 Highest Return Mean Stock Price'));
plot(timeseries(mean(Price_d(:,index(end-100:end))./mean(Price_d(:,index(end-100:end))),2),date_d,'Name','Low 100 Return Mean Stock Price'));
legend('Arithmetic Mean Stock Price','Top 100 Highest Return Mean Stock Price','Low 100 Return Mean Stock Price')
hold off

figure()
plot(fama1);
hold on
plot(fama2);
plot(fama3);
plot(bm);
legend('FamaFrenchCov lambda = 0.5','FamaFrenchCov lambda = 1','FamaFrenchCov lambda = 3','Equally Weighted Benchmark')
hold off

tsdata = horzcat(decile_6_01,decile_9_01,vig_6_01,vig_9_01);
dlmwrite('tsdata.txt', tsdata);
writetable(tsdata,'data.csv','Delimiter',' ')