function current  = assign_current(system, chip, var, TSV, type, sign)
%assign the current density for use;
    die_num = chip.N;    
    
    current = zeros(var, 1);
    %for each current excitation, the total # of unknown is var
    
    draw_P = system.drawP;
    if draw_P == 1
        drawP_die = zeros(chip.Ny, chip.Nx, chip.N);
        %this variable is for ploting
    end
    
    if type == 1
        tsvN = TSV.P;
    else
        tsvN = TSV.G;
    end
    
    const = chip.Ny*chip.Nx;
    current(const*chip.N+1:const*chip.N+tsvN, 1) = -system.Vdd/2;
    %keep in mind, the variable above chip domain (1:const*chip.N) is the
    %unknown of the volt of each node in chip
    %for the variable above chip domain, is the current flowing into each
    %pad, the excitation for them is the voltage
    
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;
    
    chip_xmesh = chip.Xmesh;
    chip_ymesh = chip.Ymesh;    
    
    for DIE = 1:1:die_num
        %this is for assigning the background power excitation
        background = 0;
        if DIE == 1
            start = 1;
        else
            start = sum(chip.blk_num(1:DIE-1))+1;
        end
        if DIE == die_num
            End = sum(chip.blk_num(1:DIE));
        else
            End = start + chip.blk_num(DIE)-1;
        end

        %calculate whether need to assign background power        
        if End-start < 0
            background = chip.power(DIE)/(chip.Xsize*chip.Ysize);
        else
            if abs(chip.power(DIE) - sum(chip.map(start:End, 5))) > 1e-8
                power_back = chip.power(DIE) - sum(chip.map(start:End, 5));
                tmp = chip.Xsize*chip.Ysize - chip.map(start:End, 3)'*chip.map(start:End, 4);
                if tmp <= 10e-12
                    area_back = chip.Xsize*chip.Ysize;
                else
                    area_back = tmp;
                end
                background = power_back / (area_back*system.Vdd);
            end
        end
        %disp(background);
        if background > 1e-8 
            for i = 1:1:gridNx_chip
                for j=1:1:gridNy_chip
                    id = i+(j-1)*gridNx_chip+const*(DIE-1);
                    if i == 1
                        boundary(1) = chip_xmesh(i);
                    else
                        boundary(1) = (chip_xmesh(i-1)+chip_xmesh(i))/2;
                    end
                    if i == gridNx_chip
                        boundary(2) = chip_xmesh(i);
                    else
                        boundary(2) = (chip_xmesh(i)+chip_xmesh(i+1))/2;
                    end

                    if j == 1
                        boundary(3) = chip_ymesh(j);
                    else
                        boundary(3) = (chip_ymesh(j-1)+chip_ymesh(j))/2;
                    end
                    if j == gridNy_chip
                        boundary(4) = chip_ymesh(j);
                    else
                        boundary(4) = (chip_ymesh(j)+chip_ymesh(j+1))/2;
                    end
                    gridx = boundary(2) - boundary(1);
                    gridy = boundary(4) - boundary(3);
                    area = gridx*gridy;
                    current(id) = current(id) + background*area;
                    if draw_P == 1
                        drawP_die(j,i,DIE) = background;                 
                    end
                end
            end  
        end

        %assign the block power maps
        if End < start
            continue
        else
            density = zeros(1:End-start);
        end
        for k=start:1:End
            density(k) = chip.map(k,5)/(chip.map(k,3)*chip.map(k,4)*system.Vdd) - background;
            xl = sum(sum(chip_xmesh<chip.map(k,1)));
            if xl <= 0
                xl = 1;
            end
            xr = sum(sum(chip_xmesh<chip.map(k,3)+chip.map(k,1)))+1;
            if xr >= gridNx_chip
                xr = gridNx_chip;
            end 
            yb = sum(sum(chip_ymesh<chip.map(k,2)));
             if yb <= 0
                yb = 1;
            end           
            yt = sum(sum(chip_ymesh<chip.map(k,4)+chip.map(k,2)))+1;
            if yt >= gridNy_chip
                yt = gridNy_chip;
            end            
            boundary_blk = [chip.map(k,1) chip.map(k,3)+chip.map(k,1) ...
                            chip.map(k,2) chip.map(k,4)+chip.map(k,2)];
            for i=xl:1:xr
                for j=yb:1:yt                    
                    id = i+(j-1)*gridNx_chip+const*(DIE-1);
                    if i == 1
                        boundary(1) = chip_xmesh(i);
                    else
                        boundary(1) = (chip_xmesh(i-1)+chip_xmesh(i))/2;
                    end
                    if i == gridNx_chip
                        boundary(2) = chip_xmesh(i);
                    else                        
                        boundary(2) = (chip_xmesh(i)+chip_xmesh(i+1))/2;
                    end

                    if j == 1
                        boundary(3) = chip_ymesh(j);
                    else
                        boundary(3) = (chip_ymesh(j-1)+chip_ymesh(j))/2;
                    end
                    if j == gridNy_chip
                        boundary(4) = chip_ymesh(j);
                    else
                        boundary(4) = (chip_ymesh(j)+chip_ymesh(j+1))/2;
                    end                                        
                    grid_area = cal_overlap (boundary, boundary_blk);
                    
                    gridx = boundary(2) - boundary(1);
                    gridy = boundary(4) - boundary(3);
                    area = gridx*gridy;
                    
                    current(id) = current(id) + density(k)*grid_area;
                    
                    if draw_P == 1
                        drawP_die(j,i,DIE) = current(id)/area;
                    end                                        
                end
            end
        end
    end

    if draw_P == 1 && sign == 1
        for k=1 : 1 : chip.N
            if chip.blk_num(k) > 0
                figure(k);
                contourf(chip.Xmesh*100, chip.Ymesh*100,drawP_die(:,:,k)*1e-6, max(chip.blk_num(k))*2,'Linestyle','none');
                h=colorbar;
                set(get(h,'Title'),'string','A/mm2','FontSize',16)
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
            end
        end
    end

end

