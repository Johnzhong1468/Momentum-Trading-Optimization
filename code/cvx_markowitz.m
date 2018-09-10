function [x] =  cvx_markowitz(mu,V,sigma,lambda,w)

n = length(mu);
U = chol(V);

    cvx_begin quiet
    
        variables x(n) z1(n) z2(n)
              
        maximize(mu'*x-lambda*x'*V*x)
        
        subject to
                    sum(x)==0;
                    x == z1+z2;
                    sum(z1) == 1;
                    abs(x)<=w;
                    z1 >= 0;
                    z2 <= 0;
    cvx_end
end