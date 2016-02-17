function write_map(system, chip, x, t)
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;    
    const = gridNx_chip * gridNy_chip;
    drawP_die = zeros(gridNy_chip, gridNx_chip, chip.N);
    file_name = 'Tmesh';
    if t == 0
        fid=fopen(file_name,'w+');
    else
        fid=fopen(file_name,'a+');
    end
    fprintf(fid,'%g\r\n',t);
    fclose(fid);      
    for k=1 : 1 : chip.N
        for i = 1:1:gridNx_chip
            for j=1:1:gridNy_chip
                id = i + (j-1)*gridNx_chip + const*(k-1);
                drawP_die(j,i,k) = x(id);
            end
        end
        st = const*(k-1)+1;
        ed = st + chip.Nx * chip.Ny - 1;
        file_name = ['chip', num2str(k)];

        if t == 0
            fid=fopen(file_name,'w+');
            fprintf(fid,'%g\r\n',x(st : ed)-system.Vdd);
            fclose(fid);
        else
            fid=fopen(file_name,'a+');
            fprintf(fid,'%g\r\n',x(st : ed)-system.Vdd);
            fclose(fid);
        end              
        if system.tran ~= 1            
            if system.map == 1
                figure(k+chip.N);
                contourf(chip.Xmesh*100, chip.Ymesh*100,abs(drawP_die(:,:,k)-system.Vdd),30, 'Linestyle','none');
                h=colorbar;
                set(get(h,'Title'),'string','Noise','FontSize',16);
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
                disp(['Die ', num2str(k), '; Max: ', num2str(max(max(drawP_die(:,:,k))))]);
                disp(['Die ', num2str(k), '; Min: ', num2str(min(min(drawP_die(:,:,k))))]);
            end
        else
            if system.map == 1
                figure(k+chip.N);
                range = system.range;
                contourf(chip.Xmesh*100, chip.Ymesh*100,drawP_die(:,:,k)-system.Vdd,30,'Linestyle','none');
                caxis([range(1),range(2)]);
                h=colorbar;                    
                set(h,'Ylim', [range(1),range(2)]);
                axis equal;
                %axis off;
                set(get(h,'Title'),'string','Vdd','FontSize',16);
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
                temp = ['Time: ', num2str(t), 's.'];
                title(temp, 'FontSize', 16);
                disp(['Die ', num2str(k), '; Max: ', num2str(max(max(drawP_die(:,:,k))))]);
                disp(['Die ', num2str(k), '; Min: ', num2str(min(min(drawP_die(:,:,k))))]);  

                filename = ['die', num2str(k), '.gif'];
                frame = getframe(k+chip.N);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);                
                if t == 0;
                    imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.2);
                end
            end
        end
    end
end

