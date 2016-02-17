function [chip_, TSV_] = mesh(chip, TSV)
%{
%%%%%%%%%%%%%%%%%%for TSV domain%%%%%%%%%%%%%%%%%%%%%
TSV.d = 10e-6; %TSV diameter
TSV.px = 400e-6; %TSV X pitch
TSV.py = 400e-6; %TSV Y pitch

TSV.map = [1.4e-3 3.0e-3 1.4e-3 5.6e-3 200e-6 200e-6
           7.3e-3 3.0e-3 1.4e-3 5.6e-3 200e-6 200e-6];
%this is the block which uses different TSV pitch
%the non-uniform meshing is based on the map

TSV.model_level = 1; %TSV model level, 2^TSV_model meshes between P and G pad

chip.Xsize = 1e-2;
chip.Ysize = 1e-2;

%}

%this function is used for meshing the chip
%TSV here is equivalent to PAD;
%we may add a function to differentiate the TSV and PAD, because P/G TSV
%may not be aligned to the PAD in the bottom die
    chip.Xgrid = (TSV.px/2)/2^TSV.model_level;
    %TSV.model_level is a parameter to control the mesh number between a
    %power PAD and a ground PAD
    chip.Ygrid = (TSV.py/2)/2^TSV.model_level;
    
    chip.Xmesh = 0:chip.Xgrid:chip.Xsize;
    chip.Ymesh = 0:chip.Ygrid:chip.Ysize;
    
    [N, ~] = size(TSV.map);
    
    for i=1:N
        xl = TSV.map(i,1);
        xr = xl+TSV.map(i,3);
        yb = TSV.map(i,2);
        yt = yb+TSV.map(i,4);
        block_Xgrid = (TSV.map(i,5)/2)/2^TSV.model_level;
        block_Ygrid = (TSV.map(i,6)/2)/2^TSV.model_level;
        
        block_Xmesh = xl:block_Xgrid:xr;
        block_Ymesh = yb:block_Ygrid:yt;
        
        chip.Xmesh = [chip.Xmesh, block_Xmesh];
        chip.Ymesh = [chip.Ymesh, block_Ymesh];
    end
    
    [chip.Xmesh, chip.Nx]= pdnt.my_unique(chip.Xmesh, chip.Xgrid/1000);
    [chip.Ymesh, chip.Ny]= pdnt.my_unique(chip.Ymesh, chip.Ygrid/1000);
    
    chip.type = zeros(chip.Ny, chip.Nx);
    
    figure(chip.N*10);
    for j=1:1:chip.Ny
        for i=1:1:chip.Nx
            x = chip.Xmesh(i);
            y = chip.Ymesh(j);
            if N > 0
                [id, sign] = pdnt.block_check(x, y, TSV.map(:,1:4));            
            else
                sign = 0;
            end
            if sign == 1
                %this step is important, the location determines the 
                %block TSV positions
                P_x = 0;
                P_y = 0;
                TSV_px = TSV.px;
                TSV_py = TSV.py;
                G_x = TSV_px/2;
                G_y = TSV_py/2;
                %the global TSV position information
                xl = TSV.map(id,1);
                yb = TSV.map(id,2);
                
                [N1, res1] = pdnt.mymod(xl-P_x, TSV_px);
                [N2, res2] = pdnt.mymod(yb-P_y, TSV_py);
                [N3, res3] = pdnt.mymod(xl-G_x, TSV_px);
                [N4, res4] = pdnt.mymod(yb-G_y, TSV_py);
                
                P_x = N1*TSV_px+res1;
                P_y = N2*TSV_py+res2;
                G_x = N3*TSV_px+res3;
                G_y = N4*TSV_py+res4;              
                TSV_px = TSV.map(id,5);
                TSV_py = TSV.map(id,6);                
                if P_x < G_x                                    
                    G_x = P_x + TSV_px/2;
                    G_y = P_y + TSV_py/2;
                else
                    P_x = G_x + TSV_px/2;
                    P_y = G_y + TSV_py/2;
                end
            else
                P_x = 0;
                P_y = 0;    
                G_x = TSV.px/2;
                G_y = TSV.py/2;
                TSV_px = TSV.px;
                TSV_py = TSV.py;                
            end
            
            [~, res1] = pdnt.mymod(x-P_x, TSV_px);
            [~, res2] = pdnt.mymod(y-P_y, TSV_py);
            [~, res3] = pdnt.mymod(x-G_x, TSV_px);
            [~, res4] = pdnt.mymod(y-G_y, TSV_py);            
            
            if  abs(res1) < 1e-12 && abs(res2) < 1e-12
                chip.type(j,i) = 1; % mark the power pad
                TSV.P = TSV.P+1;
                plot(chip.Xmesh(i)*100, chip.Ymesh(j)*100, 'o');
                hold on;                
            elseif abs(res3) < 1e-12 && abs(res4) < 1e-12
                chip.type(j,i) = 2; % mark the ground pad
                TSV.G = TSV.G+1;
                plot(chip.Xmesh(i)*100, chip.Ymesh(j)*100, 'x');
                hold on;                
            else
                chip.type(j,i) = 0; %common nodes
            end
        end
    end
%     axis equal;
    set(gca,'FontSize',16);
    xlabel('x(cm)');
    ylabel('y(cm)');set(gca,'FontSize',16);    
    
    chip_ = chip;
    TSV_ = TSV;
end

