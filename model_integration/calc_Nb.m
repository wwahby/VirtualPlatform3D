function Nb = calc_Nb(Nc)

Nb = zeros(1,length(Nc));
Nb(1) = 0;
Nb(2) = 0; % Nb summation starts with l' = 1 and ends with l' = l-1 so first nonzero Nb is Nb(l=2) (lind=3)
for lind = 3:length(Nc)
    Nb(lind) = Nb(lind-1) + Nc(lind-1);
end