warning('off','all');

p = 0.2219;
len = 16569;
result = zeros(11,1);

for d = 0:0.1:1
    EL = 0;
    for i = 1:len
        EP = 0;
        EP2 = 0;
        for j = 0:i
            EP = EP + (1-exp(-i*d))*2^(i)*nchoosek(i,j)*p^j*(0.5-p)^(i-j)*(1-p^j*(0.5-p)^(i-j))^len;
        end
        for j = 0:i-1
            EP2 = EP2 + (1-exp(-(i-1)*d))*2^(i-1)*nchoosek(i-1,j)*p^j*(0.5-p)^(i-1-j)*(1-p^j*(0.5-p)^(i-1-j))^len;
        end
%         disp(EP);
        EL = EL + i * (EP-EP2);
    end
    disp(EL);
    result(d*10+1) = EL;
end