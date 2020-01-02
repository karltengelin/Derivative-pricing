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