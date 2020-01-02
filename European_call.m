K = 100;
T = 1;
sigma = 0.2;
p = 3;
u = exp(sigma*sqrt(T/p));
d = 1/u;
S = 90;
R = exp(0.05/3);
qd = (u-R)/(u-d);
qu = (R-d)/(u-d);
%p_triag = pt(p+1);

values = zeros(p+1,p);
%probabilities = zeros(p+1,p);

for i = 1:p
    for ii = 0:i
        values(ii+1,i) = S*u^(i-ii)*d^(ii);
        %probabilities(ii+1,i) = p_triag(i+1,ii+1)*qu^(i-ii)*qd^(ii);
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


%sum = 0;

%for ii = 1:length(values)
%    if values(ii) > K
%        sum = sum + probabilities(ii)*values(ii);
%    else
        
%    end
%end

%price = rf^(-p)*sum



