function C = assign_cap(system, chip, var)
%assign the current density for use;
    die_num = chip.N;    
    
    C = zeros(var, 1);
    %for each current excitation, the total # of unknown is var
    
    draw_P = system.drawC;
    if draw_P == 1
        drawP_die = zeros(chip.Ny, chip.Nx, chip.N);
        %this variable is for ploting
    end
        
    const = chip.Ny*chip.Nx;
    
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;
    
    chip_xmesh = chip.Xmesh;
    chip_ymesh = chip.Ymesh;    
    
    for DIE = 1:1:die_num
        %this is for assigning the background power excitation
        background = chip.cap_per(DIE)*chip.c_gate;
        drawP_die(:,:,DIE) = ones(chip.Ny, chip.Nx)*chip.cap_per(DIE);
        st = 1+const*(DIE-1);
        ed = const*DIE;
        C(st:ed) = background;
        
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

        %assign the block decap maps
        if End < start
            continue
        else
            cap = zeros(1:End-start);
        end
        for k=start:1:End
            cap(k) = chip.map(k,6) - background;
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
                    
                    C(id) = C(id) + cap(k)*grid_area/area;
                    
                    if draw_P == 1
                        drawP_die(j,i,DIE) = C(id)/chip.c_gate;
                    end                                        
                end
            end
        end
    end

    if draw_P == 1
        for k=1 : 1 : chip.N
            if chip.blk_num(k) > 0
                figure(k+chip.N);
                contourf(chip.Xmesh*100, chip.Ymesh*100,drawP_die(:,:,k), max(chip.blk_num(k))*2,'Linestyle','none');
                h=colorbar;
                set(get(h,'Title'),'string','A/mm2','FontSize',16)
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
            end
        end
    end

end

