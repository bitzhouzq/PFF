function F =  gainv(M)
Mmin = min(M(:));
if Mmin < 0 
    M = M-Mmin;
end

c=0.8; 
F = 1./(1*(M+1).^c)+0.00;

end