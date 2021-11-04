function n = zetadist(z)

p = rand();
n = 0;

while(1)
    n = n + 1;
    if p < n^(-z)/zeta(z)
        break;
    else
        p = p - n^(-z)/zeta(z);
    end
end

end