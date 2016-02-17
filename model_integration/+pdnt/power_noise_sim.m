function [max_noise, max_noise_time, time_mat, voltage_mat] = power_noise_sim(chip, system, TSV)
%this is the main PDN simulator
    P = 1; G = 2; %two signs
    disp(['start time: ', num2str(fix(clock))]);
    start = tic;
    
    disp('meshing');tic;
    [chip, TSV] = pdnt.mesh(chip, TSV); %mesh the chip    
    [len, Tmesh] = pdnt.mesh_T(system);
    toc;
    
    disp('initialization');tic;
    [row_P, column_P, value_P, var_P, var_chip] = pdnt.initial(chip, TSV, P);
    [row_G, column_G, value_G, var_G, ~] = pdnt.initial(chip, TSV, G);
    toc;
    %initial the solution matrix parameters. we solve Ptsv and Gtsv, need
    %two matrix later
    
    disp('assign cap and power'); tic;
    chip.c = pdnt.assign_cap(system, chip, var_chip);
    %assign the decap distribution to the corresponding array
    
    current_full_P  = pdnt.assign_current(system, chip, var_P, TSV, P, 1);
    current_full_G  = pdnt.assign_current(system, chip, var_G, TSV, G, 0);
    %assign the current requirement to each mesh
    %due to some specific reason, need to assign according to Ptsv and Gtsv
    toc;
    power_full = chip.power;
    power_map = chip.map(:,5);
    
    if system.tran == 0                
        %this is for IR drop simulations
        [~, Y_P] = pdnt.Matrix_build(chip, system, TSV, row_P, column_P, value_P, var_P, P);
        %Y matrix is the conductance matrix
        [~, Y_G] = pdnt.Matrix_build(chip, system, TSV, row_G, column_G, value_G, var_G, G);

        current_P  = pdnt.assign_current(system, chip, var_P, TSV, P, 0);
        current_G  = pdnt.assign_current(system, chip, var_G, TSV, G, 0);
        disp(sum(system.Vdd*current_P(1:var_chip)));
        disp('solve the IR drop');
        tic;
        x_P = pdnt.Noise_solver_ss(Y_P, current_P);
        %matrix solver, call the systematical function \
        x_G = pdnt.Noise_solver_ss(Y_G, current_G);
        x = x_P(1:var_chip)+x_G(1:var_chip);
        %becasue from var_chip+1 to var, is for the current rows
        toc;      
        disp(sum(x.*current_P(1:var_chip)));
        pdnt.write_map(system, chip, x, 0);
        if system.display == 1
            pdnt.display_T(chip, x)
        end
        [max_noise, max_noise_time, time_mat, voltage_mat] = pdnt.draw_map(system, chip, Tmesh);
    else      
        [C_P, Y_P] = pdnt.Matrix_build(chip, system, TSV, row_P, column_P, value_P, var_P, P);
        [C_G, Y_G] = pdnt.Matrix_build(chip, system, TSV, row_G, column_G, value_G, var_G, G);
        %%%%%%%%%%%%%%%%%%solve the initial noise%%%%%%%%%%%%%%%%%%%%%%%%%%
        chip.map(:, 5) = power_map*system.start_ratio;
        chip.power = power_full*system.start_ratio;
        
        disp('assign power'); tic;
        current_P  = pdnt.assign_current(system, chip, var_P, TSV, P, 0);
        current_G  = pdnt.assign_current(system, chip, var_G, TSV, G, 0);
        toc;
        
%         disp('solve the IR drop');
%         tic;
%         x_P = Noise_solver_ss(Y_P, current_P);
%         x_G = Noise_solver_ss(Y_G, current_G);
%         x = x_P(1:var_chip)+x_G(1:var_chip);
%         toc;
        x_P = zeros(var_P, 1); x_P(1:var_chip) = ones(var_chip, 1)*system.Vdd/2;
        x_G = zeros(var_G, 1); x_G(1:var_chip) = ones(var_chip, 1)*system.Vdd/2;
        x = x_P(1:var_chip)+x_G(1:var_chip);
        pdnt.write_map(system, chip, x, 0);
        
        if system.display == 1
            pdnt.display_T(chip, x)
        end
        %%%%%%%%%%%%%%%%%end the initial noise calculation%%%%%%%%%%%%%%%%%    
        i = 2;
        %transient solving steps
        while i <= len
            current_prev_P = current_P;
            current_prev_G = current_G;            
            if Tmesh(i) < system.Tr + 1e-15
                disp('assign power'); tic;
                chip.map(:,5) = power_map * (system.start_ratio + system.step);
                chip.power = power_full*(system.start_ratio + system.step);
                system.start_ratio = system.start_ratio + system.step;
                current_P  = pdnt.assign_current(system, chip, var_P, TSV, P, 0);
                current_G  = pdnt.assign_current(system, chip, var_G, TSV, G, 0);
                toc;
            else                
                current_prev_P = current_P;
                current_prev_G = current_G;
                current_P = current_full_P;
                current_G = current_full_G;
            end

            disp(['Time: ', num2str(Tmesh(i)*1e9), 'ns, ' , 'solving']);            
            x_prev_P = x_P;
            x_prev_G = x_G;
            x_P = pdnt.Noise_solver(C_P, Y_P, current_P, current_prev_P, x_prev_P, Tmesh(i)-Tmesh(i-1));
            %need to use previous step of current excitation, current step
            %of current excitation; previously got solution, x_prev
            x_G = pdnt.Noise_solver(C_G, Y_G, current_G, current_prev_G, x_prev_G, Tmesh(i)-Tmesh(i-1));
            x = x_P(1:var_chip)+x_G(1:var_chip);
            
            pdnt.write_map(system, chip, x, Tmesh(i));
            
            if system.display == 1
                pdnt.display_T(chip, x)
            end
            
            %temperature returned
            i = i+1;
        end
        [max_noise, max_noise_time, time_mat, voltage_mat] = pdnt.draw_map(system, chip, Tmesh);
    end
    toc(start);
end
