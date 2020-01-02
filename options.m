clear all
close all

K = 100;
T = 1;
sigma = 0.2;
p = 3;
u = exp(sigma*sqrt(T/p));
d = 1/u;
S = 90;
R = exp(0.05/p);
qd = (u-R)/(u-d);
qu = (R-d)/(u-d);

values = zeros(p+1);
for i = 1:p+1
    for ii = 1:i
        values(ii,i) = S*u^((i-1)-(ii-1))*d^(ii-1);      
    end
end
values_stock = values;

%% European Call Option (price is identical to American Call Option)

values_call = zeros(p+1);
values_call(:,end) = max(values(:,end)-K,0);

for col = 1:p
    col = p + 1 - col
    for row = 1:col
        values_call(row,col) = R^(-1)*(qu*values_call(row,col+1)+qd*values_call(row+1,col+1));  
    end
end

price_ec = values_call(1,1)

%Price of an European call using three periods: 4.5603

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

%%Price of an European call with RNVF using three periods: 4.5603, same answer as for the european call using the tree method

%% American Put Option

values_put = zeros(p+1);
values_put(:,end) = max(K-values(:,end),0);

for col = 1:p
    col = p + 1 - col;
    for row = 1:col
        values_put(row,col) = max(R^(-1)*(qu*values_put(row,col+1)+qd*values_put(row+1,col+1)),(K-values_stock(row,col)));  
    end
end

price_ap = values_put(1,1)

%price using 1 period: 10.5757
%price using 3 period: 11.2900
%price using 10 period: 11.4437

%As period increases, there are more opportunities for a non-zero payoff,
%thus the option price increases.

%% Forward

%Here we use the same assumptions as before, except for forward contracts
%there is no option to not execute, thus payoff = Stock price @ t - strike

values_fwd = zeros(p+1);
values_fwd(:,end) = values(:,end)-K;

for col = 1:p
    col = p + 1 - col;
    for row = 1:col
        values_fwd(row,col) =R^(-1)*(qu*values_fwd(row,col+1)+qd*values_fwd(row+1,col+1));  
    end
end

%% Put-call parity inequality test

%Here we check if the put-call parity inequality is satisfied for each
%node, a 1 signifies a pass and a zero a fail (in matrices test1 and test2, test1 tests lower bound and test2 tests upper bound)

test1 = zeros(p);
test2 = zeros(p);
results = zeros(p);

for i = 1:p
    for j = 1:p
        if values_stock(i, j) - K <= values_call(i, j) - values_put(i, j)
            test1(i,j) = 1;
        end
        if values_fwd(i, j) >= values_call(i, j) - values_put(i, j)
            test2(i, j) = 1;
        end
        if test1(i,j) == test2(i,j) == 1
            results(i,j) = 1;
        end
    end
end

test1
test2

%Observe that the inequality holds true for all nodes

%% Convergence of European call RNVF


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





%% Function to get binomial coefficients
function temp = pt(n) 

temp(1, 1) = 1;
temp(2,1:2) = [1 1]; 

if n < 3
    return
end 

for r = 3 : n
    temp(r, 1) = 1;   

    for c = 2 : r-1
        temp(r, c) = temp(r-1, c-1) + temp(r-1, c);
    end   
    
    temp(r, r) = 1;
end
end