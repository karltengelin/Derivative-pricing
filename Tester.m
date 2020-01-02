K = 100;
u = 1.2;
d = 0.8;
rf = 1.05; %borde det vara 1.016396... eftersom att det är tre perioder men det är på ett år? och riskfria räntan räknas på år?
S = 90;
p = 3;
qu = 0.375;
qd = 0.625;
p_triag = pt(p+1);

values = zeros(p,p);
probabilities = zeros(p,p);
%%
for i = 2:p
    for ii = 1:i
        values(ii,i) = S*u^(
%        probabilities(ii,i) = p_triag(i,ii+1)*qu^(i-1-ii)*qd^(ii);
    end
end
%%
values_stock = values;
values(:,1:end-1) = 0;
values_contract = max(values-100,0);


for col = 1:p-1
    for row = 1:p+1-col
        values_contract(row,p-col) = rf^(-1)*(qu*values_contract(row,p-col+1)+qd*values_contract(row+1,p-col+1));  
    end
end

%%
sum = 0;

for ii = 1:length(values)
    if values(ii) > K
        sum = sum + probabilities(ii)*values(ii);
    else
        
    end
end

price = rf^(-p)*sum


