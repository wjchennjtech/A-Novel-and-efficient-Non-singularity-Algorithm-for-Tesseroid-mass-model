function [x, w] = GaussLegendre(a, b, n)

    beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2)); 
    T = diag(beta, 1) + diag(beta, -1);    
    [V, D] = eig(T);                    
    x0 = diag(D);                        
    w0 = 2*V(1,:).^2;                      

    [x0, idx] = sort(x0);
    w0 = w0(idx);

    x = 0.5*(b - a)*x0 + 0.5*(a + b);
    w = 0.5*(b - a)*w0;
    
end


