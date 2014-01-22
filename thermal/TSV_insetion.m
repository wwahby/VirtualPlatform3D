function [Cl, Cv] = TSV_insetion(Cdt_L, Cdt_V, chip, tsv, die, K_Material)
% change the conductivity of meshes which contain TSVs
    tsv.N = tsv.Nx*tsv.Ny;
    ratio = tsv.N/(chip.Nx*chip.Ny);
    
    area = chip.Xgrid*chip.Ygrid;
    
    tsv.r1 = ratio*(tsv.d-tsv.liner)^2/area;
    tsv.r2 = ratio*(tsv.d^2-(tsv.d-tsv.liner)^2)/area;
    
    for i=1:1:die.N
        for k=2:1:3
            id = 1+(i-1)*die.model+k;
            left = 1+chip.Nx*chip.Ny*(id-1);
            right = chip.Nx*chip.Ny*id;
            temp = Cdt_L(left:right);
            Cdt_L(left:right) = 1./ ((1-tsv.r1-tsv.r2)./temp ...
                                +tsv.r1/K_Material(4)+tsv.r2/K_Material(5));

            temp = Cdt_V(left:right);
            Cdt_V(left:right) = ((1-tsv.r1-tsv.r2).*temp ...
                                +tsv.r1*K_Material(4)+tsv.r2*K_Material(5));
        end
    end
    
    Cl = Cdt_L;
    Cv = Cdt_V;
end

