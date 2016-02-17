function [max_noise, max_noise_time, time_mat, voltage_mat] = draw_map(system, chip, Tmesh)
max_noise = zeros(1,chip.N);
max_noise_time = zeros(1,chip.N);

time_mat = zeros(chip.N, length(Tmesh));
voltage_mat = zeros(chip.N, length(Tmesh));
    if system.draw == 0 
        return;
    else
        gridNx_chip = chip.Nx;
        gridNy_chip = chip.Ny;
        const = chip.Nx * chip.Ny;
        drawT_die = zeros(gridNy_chip, gridNx_chip, chip.N);        
        for k=1 : 1 : chip.N
            file_name = ['chip', num2str(k)];
            data = load(file_name);
            [~, id] = max(abs(data));
            tid = floor((id-1)/const);
            St = 1+const*(tid);
            Ed = St + const -1;                        
            drawT_die(:,:,k) = reshape(data(St:Ed), gridNx_chip, gridNy_chip)';
            disp(['maximum noise occurs in ', num2str(Tmesh(tid+1)*1e9), 'ns'])
            max_noise_time(k) = Tmesh(tid+1);
            value = max(max(abs(drawT_die(:,:,k))*1000));
            disp(['maximum nose: ', num2str(value), 'mV']);
            max_noise(k) = value*1e-3;
%             figure(k+chip.N*2);
%             contourf(chip.Xmesh*100, chip.Ymesh*100,abs(drawT_die(:,:,k))*1000,30, 'Linestyle','none');
%             h=colorbar;
%             set(get(h,'Title'),'string','Noise/(mV)','FontSize',16);
%             set(gca,'FontSize',16);
%             xlabel('x(cm)');
%             ylabel('y(cm)');set(gca,'FontSize',16);
%             
            time_mat(k,:) = Tmesh(1,:);
            
            St = id - tid*const;
            Vmax = data(St:const:length(data));
            voltage_mat(k, :) = Vmax(1,:);
            
            if system.tran == 1
                figure(k+chip.N*3);
                St = id - tid*const;
                Vmax = data(St:const:length(data));
%                 plot(Tmesh*1e9, Vmax*1e3, 'linewidth', 3);
%                 set(gca,'FontSize',16);
%                 xlabel('Time(ns)');
%                 ylabel('Power Delivery Noise(mV)');set(gca,'FontSize',16);
%                 file_name = ['max_chip', num2str(k)];
% 
%                 fid=fopen(file_name,'w+');
%                 fprintf(fid,'%g\r\n',Vmax);
%                 fclose(fid);
            end
        end
    end
end

