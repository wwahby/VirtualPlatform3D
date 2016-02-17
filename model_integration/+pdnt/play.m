        Tmesh = load('Tmesh');
        gridNx_chip = 101;
        gridNy_chip = 101;
        const = 101 * 101;
        drawT_die = zeros(gridNy_chip, gridNx_chip, 2);    
        chip.Xmesh = [0:100:10000]*1e-6;
        chip.Ymesh = [0:100:10000]*1e-6;
        for k=1 : 1 : 2
            file_name = ['chip', num2str(k)];
            data = load(file_name);
            [~, id] = max(abs(data));
            tid = floor((id-1)/const);
            St = 1+const*(tid);
            Ed = St + const -1;                        

            figure(k+2*3);
            St = id - tid*const;
            Vmax = data(St:const:length(data));
            h1 = plot(Tmesh*1e9, Vmax*1e3, 'r', 'linewidth', 2);
            l1 = 'PDN';
            set(gca,'FontSize',16);
            xlabel('Time(ns)');
            ylabel('Power Delivery Noise(mV)');set(gca,'FontSize',16);
            hold on;
            if k == 1
                h2 = plot(T_P2*1e9, T_P_P1*1e3, 'b', 'linewidth', 2);
                l2 = 'PDN\_POWER';
                hold on;
                h3 = plot(T_T2*1e9, T_P_P_T1*1e3, 'k', 'linewidth', 2);
                l3 = 'PDN\_POWER\_THERMAL';                
            else
                h2 = plot(T_P2*1e9, T_P_P2*1e3, 'b', 'linewidth', 2);
                l2 = 'PDN\_POWER';                
                h3 = plot(T_T2*1e9, T_P_P_T2*1e3, 'k', 'linewidth', 2);
                l3 = 'PDN\_POWER\_THERMAL';                
            end        
            legend([h3 h2 h1 ], l3, l2, l1);
        end