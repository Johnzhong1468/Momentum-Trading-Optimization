function W = semidef(U,V)
 n = length(U);
 x = ones(n,1)/n;
 
 cvx_begin quiet
    variable W(n,n) symmetric
    minimize (norm(V-W)) 
    subject to
        W<In>semidefinite(n)
        %lambda_max(W) <= 15^2*lambda_min(W)
        lambda_min(W) > 0
        diag(W) == diag(U)
        x'*W*x == x'*U*x
        
cvx_end
end