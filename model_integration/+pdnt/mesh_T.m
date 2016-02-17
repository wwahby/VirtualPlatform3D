function [Len, mesh] = mesh_T(system)
%this function is to generate time domain non-uniform meshes
    len = floor(system.T/system.dt)+1;
    temp = zeros(len,1);

    mesh_Tr = 0:system.dt:system.Tr; 
    len_Tr = length(mesh_Tr);
    
    temp(1:len_Tr,1) = mesh_Tr;
    %rise time domain is key, use the smallest time step
    
    step = system.dt*2;
    threshold = system.dt * system.ratio;
    for i = len_Tr+1:1:len      
        if step < threshold
            temp(i) = temp(i-1) + step;
            step = system.dt + step;
        else
            temp(i) = temp(i-1) + threshold;
        end
        
        if temp(i) - system.T < -1e-12
            continue;
        else
            temp(i) = system.T;
            break;
        end
    end
    
    mesh = temp(1:i,1);
    Len = length(mesh);
end

