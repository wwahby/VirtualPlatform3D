function Pa = add_heat(power, map, blk_num, P, chip, die, var, draw_P)
    heat = zeros(var,1);
    %caculate power
    die_num = die.N;
    
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;
    const = chip.Nx*chip.Ny;
    
    chip_xmesh = chip.Xmesh;
    chip_ymesh = chip.Ymesh;
    
    if draw_P == 1
        drawP_die = zeros(chip.Ny, chip.Nx, die.N);
    end
        
    for DIE = 1:1:die_num
        background = 0;
        if DIE == 1
            start = 1;
        else
            start = sum(blk_num(1:DIE-1))+1;
        end
        if DIE == die_num
            End = sum(blk_num(1:DIE));
        else
            End = start + blk_num(DIE)-1;
        end
        
        %calculate whether need to assign background power        
        if End-start < 0
            background = power(DIE)/(chip.Xsize*chip.Ysize);
        else
            if abs(power(DIE) - sum(map(start:End, 5))) > 1e-8
                power_back = power(DIE) - sum(map(start:End, 5));
                tmp = chip.Xsize*chip.Ysize - map(start:End, 3)'*map(start:End, 4);
                if tmp <= 10e-12
                    area_back = chip.Xsize*chip.Ysize;
                else
                    area_back = tmp;
                end
                background = power_back / area_back;
            end
        end
        %disp(background);
        if background > 1e-8 
            for i = 1:1:gridNx_chip
                for j=1:1:gridNy_chip
                    id = i+(j-1)*gridNx_chip+const*(2+(DIE-1)*die.model);
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
                    grid_area = gridx*gridy;
                    heat(id) = heat(id) + background*grid_area;
                    if draw_P == 1
                        drawP_die(j,i,DIE) = heat(id)/(gridx*gridy*10^4);                 
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
            density(k) = map(k,5)/(map(k,3)*map(k,4))-background;
            xl = sum(sum(chip_xmesh<map(k,1)));
            if xl <= 0
                xl = 1;
            end
            xr = sum(sum(chip_xmesh<map(k,3)+map(k,1)))+1;
            if xr >= gridNx_chip
                xr = gridNx_chip;
            end 
            yb = sum(sum(chip_ymesh<map(k,2)));
             if yb <= 0
                yb = 1;
            end           
            yt = sum(sum(chip_ymesh<map(k,4)+map(k,2)))+1;
            if yt >= gridNy_chip
                yt = gridNy_chip;
            end            
            boundary_blk = [map(k,1) map(k,3)+map(k,1) ...
                            map(k,2) map(k,4)+map(k,2)];
            for i=xl:1:xr
                for j=yb:1:yt                    
                    id = i+(j-1)*gridNx_chip+const*(2+(DIE-1)*die.model);               
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
                    grid_area = thermal.cal_overlap (boundary, boundary_blk);
                    heat(id) = heat(id) + density(k)*grid_area;    
                    if draw_P == 1
                        gridx = boundary(2) - boundary(1);
                        gridy = boundary(4) - boundary(3);                   
                        drawP_die(j,i,DIE) = heat(id)/(gridx*gridy*10^4);                   
                    end
                end
            end
        end
    end
    Pa = P+heat;
    if draw_P == 1
        for k=1 : 1 : die.N
            if blk_num(k) > 0
                figure(die.N+1+k);
                contourf(chip.Xmesh*100, chip.Ymesh*100,drawP_die(:,:,k), max(blk_num)*2,'Linestyle','none');
                h=colorbar;
                set(get(h,'Title'),'string','Tjunc-Tamb','FontSize',16)
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
            end
        end
    end    
end

