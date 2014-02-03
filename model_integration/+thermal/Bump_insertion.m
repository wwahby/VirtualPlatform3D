function [Cl, Cv] = Bump_insertion(Cdt_L, Cdt_V, chip, bump, die, K_Material)
% change the conductivity of meshes which contain bumps
    bump.N = bump.Nx*bump.Ny;
    ratio = bump.N/(chip.Nx*chip.Ny);
    
    area = chip.Xgrid*chip.Ygrid;
    
    bump.r = ratio*bump.d^2/area;
    
    id = 1+(die.N-1)*die.model+3;
    left = 1+chip.Nx*chip.Ny*(id-1);
    right = chip.Nx*chip.Ny*id;

    temp = Cdt_L(left:right);
    Cdt_L(left:right) = 1./ ((1-bump.r)./temp +bump.r/K_Material(6));

    temp = Cdt_V(left:right);
    Cdt_V(left:right) = (1-bump.r).*temp +bump.r*K_Material(6);
    
    Cl = Cdt_L;
    Cv = Cdt_V;
end

