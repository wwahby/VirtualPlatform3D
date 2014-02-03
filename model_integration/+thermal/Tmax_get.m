function T = Tmax_get(x, die, chip)   
    T = zeros(die.N,1);
    for k=1:1:die.N
        left = chip.Nx*chip.Ny*(2+(k-1)*die.model)+1;
        right = chip.Nx*chip.Ny*(2+(k-1)*die.model+1);
        
        T(k) = max(x(left:right))-273;
    end
end

