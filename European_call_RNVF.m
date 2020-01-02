K = 100;
T = 1;
sigma = 0.2;
S = 90;

%p = 3;
pricevector = zeros(1,100);
for p = 1:100
    
    u = exp(sigma*sqrt(T/p));
    d = 1/u;
    R = exp(0.05/p);
    qd = (u-R)/(u-d);
    qu = (R-d)/(u-d);
    p_triag = pt(p+1);

    last_row = zeros(p+1,1);
    probabilities = zeros(p+1,1);

    for i = 1:p+1
        last_row(i) = S*u^(p+1-i)*d^(i-1);
        probabilities(i) = p_triag(end,i)*qu^(p+1-i)*qd^(i-1);
    end

    sum = 0;

    for ii = 1:length(last_row)
        if last_row(ii) > K
            sum = sum + probabilities(ii)*(last_row(ii)-K); 
        else

        end
    end

    price_ec_RNVF = R^(-p)*sum;
    pricevector(p) = price_ec_RNVF;
end

d1 = 1/sigma * (log(S/K) + (0.05 + sigma^2/2));
d2 = d1 - sigma;
price_BS = S * normcdf(d1, 0, 1) - K * exp(-0.05) * normcdf(d2, 0, 1);

figure
hold on
plot(pricevector)
plot(price_BS*ones(1,p))
hold off

%Notice that the call price converges to the price calculated by the
%Black-Scholes model (5.0912).
