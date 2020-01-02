
K = 100;                    %strike price
T = 1;                      %end time
sigma = 0.2;                %volatility
p = 3;                      %number of periods
u = exp(sigma*sqrt(T/p));   %up-coefficient
d = 1/u;                    %down-coefficient
S = 90;                     %stock price at time 0
R = exp(0.05/3);            %risk free rate coefficient
qd = (u-R)/(u-d);           %down probability
qu = (R-d)/(u-d);           %up probability

%% European call
values = zeros(p+1,p);

for i = 1:p
    for ii = 0:i
        values(ii+1,i) = S*u^(i-ii)*d^(ii);
    end
end

values_stock = values;
values(:,1:end-1) = 0;
values_contract = max(values-K,0);

for col = 1:p-1
    for row = 1:p+1-col
        values_contract(row,p-col) = R^(-1)*(qu*values_contract(row,p-col+1)+qd*values_contract(row+1,p-col+1));  
    end
end

price_ec = R^(-1)*(qu*values_contract(1,1)+qd*values_contract(2,1))

%% European call RNVF

p_triag = pt(p+1); %pt is a function that constructs a pascal triangle, this we do in order to get the binomial coefficients (we didn't know about the "binomcdf"-command...)

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

price_ec_RNVF = R^(-p)*sum

%same answer as for the european call using tree method

%% American put

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

