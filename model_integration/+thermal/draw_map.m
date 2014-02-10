function draw_map( chip, die, pack, Layer, h, x, granularity)
    %%%%%%%%%%variable conversion%%%%%%%%%%%%%%%%%%
    Ta = h.Ta;
    drawT_die = zeros(chip.Ny, chip.Nx, die.N);
    drawT_pack= zeros(pack.Ny, pack.Nx);
    
    const = chip.Nx*chip.Ny;
    for k=1:1:die.N
        for i = 1:1:chip.Nx
            for j=1:1:chip.Nx
                id = i+(j-1)*chip.Nx+const*(2+die.model*(k-1));
                drawT_die(j,i,k) = x(id);  
            end
        end
    end
    const = chip.Nx*chip.Ny*(Layer.N-1);
    for i = 1:1:pack.Nx
        for j=1:1:pack.Ny
            id = i+(j-1)*pack.Nx+const;
            drawT_pack(j,i) = x(id);      
        end
    end
    
    for k=1 : 1 : die.N
        figure(k);
        contourf(chip.Xmesh*100, chip.Ymesh*100,drawT_die(:,:,k)-Ta, granularity,'Linestyle','none');
        h=colorbar;
        set(get(h,'Title'),'string','Tjunc-Tamb','FontSize',16)
        set(gca,'FontSize',16);
        xlabel('x(cm)');
        ylabel('y(cm)');set(gca,'FontSize',16);
    end
    
    figure(die.N+1);
    contourf(pack.Xmesh*100, pack.Ymesh*100,drawT_pack(:,:)-Ta, granularity,'Linestyle','none');
    h=colorbar;
    set(get(h,'Title'),'string','Tjunc-Tamb','FontSize',16)
    set(gca,'FontSize',16);
    xlabel('x(cm)');
    ylabel('y(cm)');set(gca,'FontSize',16);    
end

