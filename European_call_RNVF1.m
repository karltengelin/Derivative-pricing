K = 100;
T = 1;
sigma = 0.2;
S = 90;

p = 3;
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
            sum = sum + probabilities(ii)*(last_row(ii)-K); %ska det vara -K eller bara inget?
        else

        end
    end

    price_ec_RNVF = R^(-p)*sum;
    pricevector(p) = price_ec_RNVF;
end

plot(pricevector)