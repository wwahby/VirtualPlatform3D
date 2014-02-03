function [ cXmesh cYmesh, pXmesh, pYmesh, ...
           cNx, cNy, pNx, pNy] = mesh(chip, pack)
    %chip uniform meshing
    cXmesh = [0 : chip.Xgrid : chip.Xsize, chip.Xsize];
    cYmesh = [0 : chip.Ygrid : chip.Ysize, chip.Xsize];
    
    cNx = length(cXmesh);
    if abs(cXmesh(cNx-1) - chip.Xsize)<1e-9
        cXmesh = cXmesh(1:cNx-1);
        cNx = cNx - 1;
    end
    
    cNy = length(cYmesh);
    if abs(cYmesh(cNy-1) - chip.Ysize)<1e-9
        cYmesh = cYmesh(1:cNy-1);
        cNy = cNy - 1;
    end  
    
    %interposer uniform meshing
    left = (pack.Xsize - chip.Xsize)/2;
    right = (pack.Xsize + chip.Xsize)/2;
    
    left_part = 0 : pack.Xgrid : left-pack.Xgrid/100;
    middle_part = left : chip.Xgrid : right;
    right_part = sort(pack.Xsize-left_part);
    
    pXmesh = [left_part, middle_part, right_part];

    left = (pack.Ysize - chip.Ysize)/2;
    right = (pack.Ysize + chip.Ysize)/2;
    
    left_part = 0 : pack.Ygrid : left-pack.Ygrid/100;
    middle_part = left : chip.Ygrid : right;
    right_part = sort(pack.Ysize-left_part);
    
    pYmesh = [left_part, middle_part, right_part];
    
    %caculate the mesh number
    
    pNx = length(pXmesh);
    pNy = length(pYmesh);   
end

