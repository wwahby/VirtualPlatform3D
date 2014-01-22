function Pa = add_heat(power, P, chip, die, var)
    heat = zeros(var,1);
    %caculate power

    power = power/(chip.Xsize*chip.Ysize*1e4);
    die_num = die.N;
    
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;
    const = chip.Nx*chip.Ny;
    
    chip_xmesh = chip.Xmesh;
    chip_ymesh = chip.Ymesh;
    
    bgd_die = power;
    
    for k = 1:1:die_num
        for i = 1:1:gridNx_chip
            for j=1:1:gridNy_chip
                id = i+(j-1)*gridNx_chip+const*(2+(k-1)*die.model);
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
                heat(id) = heat(id) + bgd_die(k)*10^4*grid_area;
            end
        end
    end

    Pa = P+heat;
end

