function T = ThermSim( die, thick, chip, pack, ...
                   tsv, bump, portion, power, ...
                   granularity, draw, h, displayT)
    disp(['start time: ', num2str(fix(clock))]);
    start = tic;
%%%%%%%%%%%%%%%%%%%%material information for the chip%%%%%%%%%%%%%%%%%%%%%%
%thermal conductivity
    K_Material = [3  %TIM 1
                  149 %CHIP 2
                  0.9 %UNDER 3
                  400 %COPPER 4
                  1.38 %SIO2 5
                  60 %UBUMP 6
                  149 %SILICON INTERPOSER 7
                  0.024 %Air 8
                  400*portion+1.38*(1-portion) %ILD vertical 9
                  1/(portion/400+(1-portion)/1.38) %ILD lateral 10
                  ]; 
    
    Layer = stack_build(die, thick);
    %Layer.material, Layer.thick, Layer.N
%%%%%%%%%%%%%%%%%%%%finish material information for%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%start meshing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('starting meshing the package');
    tic;
    [ chip.Xmesh chip.Ymesh, pack.Xmesh, pack.Ymesh, ...
      chip.Nx, chip.Ny, pack.Nx, pack.Ny ] = mesh(chip, pack);
    disp('finish meshing the package');
    toc;   
%%%%%%%%%%%%%%%%%%%%%%%%%finish meshing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%array initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Cdt_L, Cdt_V, row, column, value, var] ...
    = initial(chip, pack, Layer, K_Material);
%%%%%%%%%%%%%%%%%%%%%%%%%fnish initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%TSV & bump insertion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Cdt_L, Cdt_V] = TSV_insetion(Cdt_L, Cdt_V, chip, ...
                                  tsv, die, K_Material);
    [Cdt_L, Cdt_V] = Bump_insertion(Cdt_L, Cdt_V, chip, ... 
                                    bump, die, K_Material);
%%%%%%%%%%%%%%%%%%%%%%%%finish insertion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%build thermal conductance materix A%%%%%%%%%%%%%%%%%%%%
    %layers beyond board
    [A P] = Conductance(chip, pack, h, Cdt_L, Cdt_V, Layer, ...
                    row, column, value, var);
%%%%%%%%%%%%%%%%%%%%finish thermal conductance materix A%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%add heat matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = add_heat(power, P, chip, die, var);
%%%%%%%%%%%%%%%%%%%%finish adding heat%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%solve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('calculate temperature');
    tic;
    x = A\P;
    toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%finish solving%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%draw map%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%in 3D system tool, this function should be disabled%%%%%%%%%%%%
    if draw == 1
        draw_map( chip, die, pack, Layer, h, x, granularity);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%finish drawing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%return the maximum temperature of each die%%%%%%%%%%%%%%%%%
    T = Tmax_get(x, die, chip);
    if displayT == 1
        disp(T);
    end
%%%%%%%%%%%%%%%%%temperature returned%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp(['total run time:', num2str(toc(start))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%end here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




















