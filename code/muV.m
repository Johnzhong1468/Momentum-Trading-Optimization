function [mu,V] = muV(Prices,trade_date, ...
    horizon,sample_frequency,number_of_samples,rate_of_decay,meanmode,varmode,decile,index)

%mu: mean vector of historical returns

%V: covariance matrix of historical returns

%Prices: Matrix of asset prices, with rows corresponding to time periods
% (ordered chronologically), and with columns corresponding to assets.

%trade_date: The number of the period prior to which data can be
% accessed for calculating mu and V.  Think of a new portfolio as being
% constructed at the beginning of the period whose number is trade_date.
% The only prices then known would be ones from earlier periods.)

%horizon: The number of periods which the new portfolio will be held.
% (For example, if the portfolio is to be held one month and the data
% is weekly, then horizon=4, whereas if data is daily, horizon=20.)

%sample_frequency: The number of periods between sampled returns in 
% calculating mu and V.  For example, if the historical returns are
% to be sampled weekly, then sample_frequency=1 if data is weekly,
% whereas sample_frequency=5 if data is daily.

%number_of_samples: The number of historical returns to compute. Must
% satisfy the condition 
%          sample_frequency*number_of_samples < trade_date  
% since otherwise there will insufficient historical data.

%rate_of_decay: Each of the historical returns will be assigned a weight.
% The most recent return will be given the greatest weight.  The
% weight of the next most recent return will be (1-rate_of_decay) times
% the weight of the most recent return.  Likewise, the weight of the
% second most recent return will be (1-rate_of_decay)^2 times the
% weight of the most recent return.  And so on.
% The weights will sum to 1, thus giving a probability measure on the 
% samples.  This is the probability measure used in computing mu and V.
% Be careful to satisfy the condition 0 <= rate_of_decay < 1.
% (Clearly, if all samples are to be given the same weight, choose
% rate_of_decay=0.
% BE CAREFUL to change rate_of_decay when changing sample_frequency.
% In particular, if r1 is the rate chosen for frequency f1, then for
% frequency f2 the equivalent rate is 
%                  r2 = 1 - (1-r1)^(f2/f1)
% These two rates would result in giving the same relative weight 
% to a sample return from one month in the past, one year in the past, 
% or however long ago. In other words, if r1 causes samples from
% one year ago to be given half as much weight as the most recent sample,
% so will r2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = horizon; P = Prices; t_d = trade_date; s_f = sample_frequency;
n_s = number_of_samples; r_d = rate_of_decay;

s_d = t_d - 1 - s_f*[0:n_s];    % -2 to account for reversal effect
s_d = fliplr(s_d);              % chronological
% sample dates, a row vector

S_P = P(s_d,:);
% Sample Prices, a matrix

S_C_R = log(S_P(2:end,:)./S_P(1:(end-1),:));
% Sample Compound Returns, a matrix

% now it ’s time to construct the weight
w=(1-r_d).^(1:n_s);
w=fliplr(w)'/sum(w);
% weights, a non-negative vector that sums to 1

mean_c_r = S_C_R'*w;
% mean vector of compound returns

Cov_C_R = (S_C_R'*diag(w)*S_C_R) - mean_c_r*mean_c_r';
% covariance matrix of compound returns

adapted_mean_c_r = (h/s_f)*mean_c_r;
% adapting mean vector to reflect length of holding period

Adapted_Cov_C_R = (h/s_f)*Cov_C_R;
% adapting covariance matrix to length of holding period

muu = exp(adapted_mean_c_r + .5*diag(Adapted_Cov_C_R)); 
mu = muu - 1;
% resulting mean vector of ARITHMETIC returns

V = (muu*muu').*(exp(Adapted_Cov_C_R) - 1);
% resulting covariance matrix of ARITHMETIC returns

len = length(mu);
thre_l=prctile(mu,decile,2);
thre_h=prctile(mu,100-decile,2);


switch meanmode
    case  1
        %simple return
        mu=mu;
    case 2
        %101 momentum
        mu=zeros(len,1);
        mu(mu>thre_h)=1;
        mu(mu<thre_l)=-1;
    case 3
        mu = sum(exp(S_C_R)-1)';
end


switch varmode
    
    case  1
        %simple covariance (may be out of rank)
        V = V;
    case {2,3,4}
        load('p2data.mat');
        Price = Price_d;
        [~,ind1] = min(abs(datenum(date_d)-datenum(date_mo(t_d))));
        [~,ind2] = min(abs(datenum(date_d)-datenum(date_mo(t_d-6))));
        t = ind1-ind2;
        return_d = Price_d(ind2:ind1-1,index)./Price_d(ind2-1:ind1-2,index)-1;
        n = length(index);
        load('ffdata_d.mat')
        b_m = zeros(n,1);
        b_s = zeros(n,1);
        b_v = zeros(n,1);
        res = zeros(n,1);
        delta = zeros(n,n);

        x1 = mkt_rf(ind2:ind1-1);
        x2 = smb(ind2:ind1-1);
        x3 = hml(ind2:ind1-1);
        for j = 1:n
            y = return_d(1:end,j);
            X = [ones(size(x1)) x1 x2 x3];
            [b,bint,r] = regress(y,X);
            b_m(j) = b(2);
            b_s(j) = b(3);
            b_v(j) = b(4);
            res(j) = var(r);
            delta(j,j) = res(j);
        end
            bs = [b_m,b_s,b_v];
        if varmode == 2
        %fama french covariance
            factor = [mkt_rf,smb,hml];
            factor_cov = cov(factor);
            V_ff = bs*factor_cov*bs'+delta;
            V = V_ff;
        elseif varmode == 3
        %singe factor
            factor_cov = cov(mkt_rf);
            V_sim = b_m*factor_cov*b_m'+delta;
            V = V_sim;
        else
        %shrinkage
        factor = [mkt_rf,smb,hml];
        factor_cov = cov(factor);
        n_bs = [];
        for i = 1:n
            n_bs =[n_bs (bs'*bs/std(res)+cov(bs)^(-1))^(-1)*(bs'*bs*bs(i,:)'/std(res)+cov(bs)^(-1)*mean(bs)')];
        end
        n_bs = n_bs';
        V_sim = n_bs*factor_cov*n_bs'+delta;
        V = V_sim;        
    
        end
end
end




