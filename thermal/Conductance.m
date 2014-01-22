function [A, Pa] = Conductance(chip, pack, hc, Cdt_L, Cdt_V, Layer, ...
                             row, column, value, var)
    tic
    chip_xmesh = chip.Xmesh;
    chip_ymesh = chip.Ymesh;
    pack_xmesh = pack.Xmesh;
    pack_ymesh = pack.Ymesh;
    
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;
    gridNx_pack = pack.Nx;
    gridNy_pack = pack.Ny;
    
    Thick_layer = Layer.thick;
    
    h_s = hc.up;
    h_mf = hc.down;
    h_d = hc.d;
    h = hc.side;
    Ta = hc.Ta;
    
    layer_id = 1;
    disp(['analyze: layer' num2str(layer_id)]);    
    
    P = zeros(var,1); 
    
    pointer = 1;
    
    xl = (pack.Nx-chip.Nx)/2+1;
    xr = xl+chip.Nx-1;
    yb = (pack.Ny-chip.Ny)/2+1;   
    yt = yb+chip.Ny-1;
    for i = 1:1:gridNx_chip
        for j = 1:1:gridNy_chip
            % z information
            lz_u = 0;        
            lz_d = Thick_layer(layer_id);
            lz = (lz_d+lz_u)/2;
            
            % x direction
            if(i>1) 
                x1 = chip_xmesh(i)-chip_xmesh(i-1);
            else
                x1 = 0;
            end
            
            if(i<gridNx_chip)
                x2 = chip_xmesh(i+1) - chip_xmesh(i);
            else
                x2 = 0;
            end
            
            % y direction
            if(j>1) 
                y1 = chip_ymesh(j)-chip_ymesh(j-1);
            else
                y1 = 0;
            end
            
            if(j<gridNy_chip)
                y2 = chip_ymesh(j+1) - chip_ymesh(j);
            else
                y2 = 0;
            end
            
            lx = (x1+x2)/2;
            ly = (y1+y2)/2;
            
            const = gridNx_chip*gridNy_chip*(layer_id-1);            
            id = i+(j-1)*gridNx_chip+const;
            w = id-1;
            e = id+1;
            s = id - gridNx_chip;
            n = id + gridNx_chip;
            b = id+gridNx_chip*gridNy_chip;

            %A is diagonal dominant
            %variable temp is for diagonal element calculation
            if( ( i==1 || i==gridNx_chip ) && ( j==1 || j==gridNy_chip ))
                if i == 1 && j == 1 
                        ke_d = Cdt_L(id);kn_d = Cdt_L(id);kz_d = Cdt_V(id);

                        P(id)=Ta*h*ly*lz+Ta*h*lx*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d*ly*lz/x2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d*lx*lz/y2;
                        temp = temp-value(pointer);                   
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d*lx*ly/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;

                elseif i == gridNx_chip && j == 1
                        kw_d = Cdt_L(id-1);kn_d = Cdt_L(id-1);kz_d = Cdt_V(id-1);                

                        P(id)=Ta*h*ly*lz+Ta*h*lx*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d*ly*lz/x1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d*lx*lz/y2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d*lx*ly/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;

                elseif i == 1 && j == gridNy_chip
                        ke_d = Cdt_L(s);ks_d = Cdt_L(s);kz_d = Cdt_V(s);                 

                        P(id)=Ta*h*ly*lz+Ta*h*lx*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d*ly*lz/x2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d*lx*lz/y1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d*lx*ly/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                         

                else
                        kw_d = Cdt_L(s-1);ks_d = Cdt_L(s-1);kz_d = Cdt_V(s-1);                 
                        P(id)=Ta*h*ly*lz+Ta*h*lx*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d*ly*lz/x1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d*lx*lz/y1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;


                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d*lx*ly/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1; 

                end
            elseif i ~= 1 && i~= gridNx_chip && j ~= 1 && j~=gridNy_chip
                kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz/2;
                ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz/2;
                ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz/2;
                kn_d = (Cdt_L(id)*x1+Cdt_L(id-1)*x2)*lz/2;            
                kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;     

                P(id)=Ta*h_s*lx*ly;
                temp = P(id)/Ta;

                row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/x2;
                temp = temp-value(pointer);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/y2;
                temp = temp-value(pointer);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/x1;
                temp = temp-value(pointer);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/y1;
                temp = temp-value(pointer);
                pointer=pointer+1;           

                row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                temp = temp-value(pointer);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=id; value(pointer)=temp; 
                pointer=pointer+1;

            else
                if i == 1
                        ke_d = (Cdt_L(id)*y2+Cdt_L(s)*y1)*lz/2;
                        kn_d = Cdt_L(id)*lx*lz;ks_d = Cdt_L(s)*lx*lz;
                        kz_d = (Cdt_V(id)*y2*x2+Cdt_V(s)*y1*x2)/4;

                        P(id)=Ta*h*ly*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/y1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/y2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/x2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                     

                elseif i == gridNx_chip
                        kw_d = (Cdt_L(id-1)*y2+Cdt_L(s-1)*y1)*lz/2;
                        kn_d = Cdt_L(id-1)*lx*lz;ks_d = Cdt_L(s-1)*lx*lz;
                        kz_d = (Cdt_V(id-1)*y2+Cdt_V(s-1)*y1)*x1/4;                

                        P(id)=Ta*h*ly*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/y1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/y2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/x1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                    

                elseif j == gridNy_chip
                        ks_d = (Cdt_L(s)*x2+Cdt_L(s-1)*x1)*lz/2;
                        kw_d = Cdt_L(s-1)*ly*lz;ke_d = Cdt_L(s)*ly*lz;
                        kz_d = (Cdt_V(s)*x2+Cdt_V(s-1)*x1)*y1/4;                 

                        P(id)=Ta*h*lx*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/x2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/x1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/y1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                     

                else
                        kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz/2;
                        kw_d = Cdt_L(id-1)*ly*lz;ke_d = Cdt_L(id)*ly*lz;
                        kz_d = (Cdt_V(id)*x2+Cdt_V(id-1)*x1)*y2/4;                

                        P(id)=Ta*h*lx*lz+Ta*h_s*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/x2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/x1;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/y2;
                        temp = temp-value(pointer);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                        temp = temp-value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                       

                end
            end        
        end
    end
    disp(['finish: layer' num2str(layer_id)]);
    toc;
    for layer_id = 2:1:Layer.N-1
        disp(['analyze: layer' num2str(layer_id)]);
        tic;
        for i = 1:1:gridNx_chip
            for j = 1:1:gridNy_chip
                lz_u = Thick_layer(layer_id-1);       
                lz_d = Thick_layer(layer_id);
                lz = (lz_u+lz_d)/2;
                
                % x direction
                if(i>1) 
                    x1 = chip_xmesh(i)-chip_xmesh(i-1);
                else
                    x1 = 0;
                end

                if(i<gridNx_chip)
                    x2 = chip_xmesh(i+1) - chip_xmesh(i);
                else
                    x2 = 0;
                end

                % y direction
                if(j>1) 
                    y1 = chip_ymesh(j)-chip_ymesh(j-1);
                else
                    y1 = 0;
                end

                if(j<gridNy_chip)
                    y2 = chip_ymesh(j+1) - chip_ymesh(j);
                else
                    y2 = 0;
                end
                
                lx = (x1+x2)/2;
                ly = (y1+y2)/2;
                
                const = gridNx_chip*gridNy_chip*(layer_id-1);            
                id = i+(j-1)*gridNx_chip+const;
                w = id-1;
                e = id+1;
                s = id - gridNx_chip;
                n = id + gridNx_chip;
                if layer_id == Layer.N-1
                    b = (i-1+xl)+(j-2+yb)*gridNx_pack+const+gridNx_chip*gridNy_chip;
                else b = id+gridNx_chip*gridNy_chip;
                end
                t = id-gridNx_chip*gridNy_chip;
                t_s = t - gridNx_chip;

                if( ( i==1 || i==gridNx_chip ) && ( j==1 || j==gridNy_chip ))
                    if i == 1 && j == 1
                            ke_d = Cdt_L(id)*y2*lz_d;kn_d = Cdt_L(id)*x2*lz_d;kz_d = Cdt_V(id)*x2*y2;  
                            ke_u = Cdt_L(t)*y2*lz_u;kn_u = Cdt_L(t)*x2*lz_u;kz_u = Cdt_V(t)*x2*y2;

                            P(id)=1/2*(Ta*h*y2*lz+Ta*h*x2*lz);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1;                    

                    elseif i == gridNx_chip && j == 1
                            kw_d = Cdt_L(id-1)*y2*lz_d;kn_d = Cdt_L(id-1)*x1*lz_d;kz_d = Cdt_V(id-1)*x1*y2;
                            kw_u = Cdt_L(t-1)*y2*lz_u;kn_u = Cdt_L(t-1)*x1*lz_u;kz_u = Cdt_V(t-1)*x1*y2;

                            P(id)=1/2*(Ta*h*y2*lz+Ta*h*x1*lz);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1;                      

                    elseif i == 1 && j == gridNy_chip
                            ke_d = Cdt_L(s)*y1*lz_d;ks_d = Cdt_L(s)*x2*lz_d;kz_d = Cdt_V(s)*y1*x2;
                            ke_u = Cdt_L(t_s)*y1*lz_u;ks_u = Cdt_L(t_s)*x2*lz_u;kz_u = Cdt_V(t_s)*y1*x2;

                            P(id)=1/2*(Ta*h*y1*lz+Ta*h*x2*lz);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1;  

                    else
                            kw_d = Cdt_L(s-1)*y1*lz_d;ks_d = Cdt_L(s-1)*x1*lz_d;kz_d = Cdt_V(s-1)*x1*y1; 
                            kw_u = Cdt_L(t_s-1)*y1*lz_u;ks_u = Cdt_L(t_s-1)*x1*lz_u;kz_u = Cdt_V(t_s-1)*x1*y1;

                            P(id)=1/2*(Ta*h*y1*lz+Ta*h*x1*lz);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u/lz_u;               
                            temp = temp - value(pointer);
                            pointer=pointer+1;                    

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;        
                            pointer=pointer+1;                      
                    end

                elseif i ~= 1 && i~= gridNx_chip && j ~= 1 && j~=gridNy_chip
                    kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                    ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                    ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                    kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                    kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                    kw_u = (Cdt_L(t_s-1)*y1+Cdt_L(t-1)*y2)*lz_u/2;
                    ks_u = (Cdt_L(t_s-1)*x1+Cdt_L(t_s)*x2)*lz_u/2;
                    ke_u = (Cdt_L(t_s)*y1+Cdt_L(t)*y2)*lz_u/2;
                    kn_u = (Cdt_L(t)*x2+Cdt_L(t-1)*x1)*lz_u/2;            
                    kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                    row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_u/x2-1/2*ke_d/x2;
                    temp = - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_u/y2-1/2*kn_d/y2;
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_u/x1-1/2*kw_d/x1;
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_u/y1-1/2*ks_d/y1;
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/lz_u;
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                    pointer=pointer+1;                      

                else
                    if i == 1
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = Cdt_L(id)*x2*lz_d;ks_d = Cdt_L(s)*x2*lz_d;
                            kz_d = (Cdt_V(id)*y2+Cdt_V(s)*y1)*x2/2;

                            ke_u = (Cdt_L(t_s)*y1+Cdt_L(t)*y2)*lz_u/2;
                            kn_u = Cdt_L(t)*x2*lz_u;ks_u = Cdt_L(t_s)*x2*lz_u;
                            kz_u = (Cdt_V(t)*y2+Cdt_V(t_s)*y1)*x2/2;                        
                            
                            P(id)=Ta*h*ly*lz;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_u/x2-1/2*ke_d/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1;                    

                    elseif i == gridNx_chip
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            kn_d = Cdt_L(id-1)*x1*lz_d;ks_d = Cdt_L(s-1)*x1*lz_d;
                            kz_d = (Cdt_V(id-1)*y2+Cdt_V(s-1)*y1)*x1/2;  

                            kw_u = (Cdt_L(t_s-1)*y1+Cdt_L(t-1)*y2)*lz_u/2;
                            kn_u = Cdt_L(t-1)*x1*lz_u;ks_u = Cdt_L(t_s-1)*x1*lz_u;
                            kz_u = (Cdt_V(t-1)*y2+Cdt_V(t_s-1)*y1)*x1/2;                         

                            P(id)=Ta*h*ly*lz;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_u/x1-1/2*kw_d/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;                              
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1;                      

                    elseif j == gridNy_chip
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            kw_d = Cdt_L(s-1)*y1*lz_d;ke_d = Cdt_L(s)*y1*lz_d;
                            kz_d = (Cdt_V(s)*x2+Cdt_V(s-1)*x1)*y1/2;  

                            ks_u = (Cdt_L(t_s-1)*x1+Cdt_L(t_s)*x2)*lz_u/2;
                            kw_u = Cdt_L(t_s-1)*y1*lz_u;ke_u = Cdt_L(t_s)*y1*lz_u;
                            kz_u = (Cdt_V(t_s)*x2+Cdt_V(t_s-1)*x1)*y1/2;                         

                            P(id)=Ta*h*lx*lz;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_u/y1-1/2*ks_d/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1; 

                    else
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;
                            kw_d = Cdt_L(id-1)*y2*lz_d;ke_d = Cdt_L(id)*y2*lz_d;
                            kz_d = (Cdt_V(id)*x2+Cdt_V(id-1)*x1)*y2/2;  

                            kn_u = (Cdt_L(t)*x2+Cdt_L(t-1)*x1)*lz_u/2;
                            kw_u = Cdt_L(t-1)*y2*lz_u;ke_u = Cdt_L(t)*y2*lz_u;
                            kz_u = (Cdt_V(t)*x2+Cdt_V(t-1)*x1)*y2/2;                        

                            P(id)=Ta*h*lx*lz;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_u/y2-1/2*kn_d/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;                            
                            temp = temp - value(pointer);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                            pointer=pointer+1;                     

                    end
                end
            end
        end
        disp(['finish: layer' num2str(layer_id)]);
        toc;
    end
    %top of board layer
    tic;
    layer_id = Layer.N;
    disp(['analyze: layer' num2str(layer_id)]);   
    for i = 1:1:gridNx_pack
        for j=1:1:gridNy_pack
            lz_u = Thick_layer(layer_id-1);
            lz_d = Thick_layer(layer_id);
            
            if(i>1) 
                x1 = pack_xmesh(i)-pack_xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<gridNx_pack)
                x2 = pack_xmesh(i+1) - pack_xmesh(i);
            else
                x2 = 0;
            end

            % y direction
            if(j>1) 
                y1 = pack_ymesh(j)-pack_ymesh(j-1);
            else
                y1 = 0;
            end

            if(j<gridNy_pack)
                y2 = pack_ymesh(j+1) - pack_ymesh(j);
            else
                y2 = 0;
            end

            lx = (x1+x2)/2;
            ly = (y1+y2)/2;
            const = gridNx_chip*gridNy_chip*(layer_id-1);            
            id = i+(j-1)*gridNx_pack+const;
            w = id-1;
            e = id+1;
            s = id - gridNx_pack;
            n = id + gridNx_pack;
            b = id+gridNx_pack*gridNy_pack;
            t = ((i-xl+1)+(j-yb)*gridNx_chip)+(layer_id-2)*gridNx_chip*gridNy_chip; 
            t_s = t - gridNx_chip;

            if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ))
                if i == 1 && j == 1
                        ke_d = Cdt_L(id);kn_d = Cdt_L(id);kz_d = Cdt_V(id);

                        P(id)=1/4*(Ta*h*y2*lz_d+Ta*h*x2*lz_d+Ta*h*x2*y2);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_d*y2*lz_d/x2;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_d*x2*lz_d/y2;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d*x2*y2/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;

                elseif i == gridNx_pack && j == 1
                        kw_d = Cdt_L(id-1);kn_d = Cdt_L(id-1);kz_d = Cdt_V(id-1);

                        P(id)=1/4*(Ta*h*y2*lz_d+Ta*h*x1*lz_d+Ta*h*x1*y2);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_d*y2*lz_d/x1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_d*x1*lz_d/y2;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d*x1*y2/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;                 

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                     

                elseif i == 1 && j == gridNy_pack
                        ke_d = Cdt_L(s);ks_d = Cdt_L(s);kz_d = Cdt_V(s);

                        P(id)=1/4*(Ta*h*y1*lz_d+Ta*h*x2*lz_d+Ta*h*x2*y1);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_d*y1*lz_d/x2;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_d*x2*lz_d/y1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d*x2*y1/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                       

                else
                        kw_d = Cdt_L(s-1);ks_d = Cdt_L(s-1);kz_d = Cdt_V(s-1);

                        P(id)=1/4*(Ta*h*y1*lz_d+Ta*h*x1*lz_d+Ta*h*x1*y1);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_d*y1*lz_d/x1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_d*x1*lz_d/y1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/4*kz_d*x1*y1/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                                           

                end
            elseif i ~= 1 && i ~= gridNx_pack && j ~= 1 && j ~= gridNy_pack
                if(  i >= xl && i<= xr && j >= yb && j <= yt )
                    if (i==xl || i==xr) && (j==yb || j==yt)
                        if i==xl && j==yb
                            ke_u = Cdt_L(t)*y2*lz_u;kn_u = Cdt_L(t)*x2*lz_u;kz_u = Cdt_V(t)*x2*y2;
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;
 
                            P(id)=1/4*(Ta*h*y2*lz_u+Ta*h*x2*lz_u+Ta*h*x1*y2+Ta*h*x1*y1+Ta*h*x2*y1);
                            temp = P(id)/Ta;                            

                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_u/(4*x2)-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_u/(4*y2)-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/(lz_d);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(4*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)= temp;
                            pointer=pointer+1;
                        elseif i==xl && j==yt
                            ke_u = Cdt_L(t_s)*y1*lz_u;ks_u = Cdt_L(t_s)*x2*lz_u;kz_u = Cdt_V(t_s)*x2*y1;                        
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;                        

                            P(id)=1/4*(Ta*h*y1*lz_u+Ta*h*x2*lz_u+Ta*h*x1*y2+Ta*h*x1*y1+Ta*h*x2*y2);
                            temp = P(id)/Ta;                            
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_u/(4*x2)-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_u/(4*y1)-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(4*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;             
                            pointer=pointer+1;

                        elseif i==xr && j==yb
                            kw_u = Cdt_L(t-1)*y2*lz_u;kn_u = Cdt_L(t-1)*x1*lz_u;kz_u = Cdt_V(t-1)*x1*y2;
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;                        

                            P(id)=1/4*(Ta*h*y2*lz_u+Ta*h*x1*lz_u+Ta*h*x2*y1+Ta*h*x1*y1+Ta*h*x2*y2);
                            temp = P(id)/Ta;                            
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_u/(4*y2)-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_u/(4*x1)-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(4*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;          
                            pointer=pointer+1;
                        else
                            kw_u = Cdt_L(t_s-1)*y1*lz_u;ks_u = Cdt_L(t_s-1)*x1*lz_u;kz_u = Cdt_V(t_s-1)*x1*y1;
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            P(id)=1/4*(Ta*h*y1*lz_u+Ta*h*x1*lz_u+Ta*h*x2*y1+Ta*h*x1*y2+Ta*h*x2*y2);
                            temp = P(id)/Ta;                            
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_u/(4*x1)-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_u/(4*y1)-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(4*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)= temp;            
                            pointer=pointer+1;
                        end

                    elseif i~=xl && i~=xr && j~=yb && j~=yt
                        kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                        ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                        ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                        kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                        kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                        kw_u = (Cdt_L(t_s-1)*y1+Cdt_L(t-1)*y2)*lz_u/2;
                        ks_u = (Cdt_L(t_s-1)*x1+Cdt_L(t_s)*x2)*lz_u/2;
                        ke_u = (Cdt_L(t_s)*y1+Cdt_L(t)*y2)*lz_u/2;
                        kn_u = (Cdt_L(t)*x2+Cdt_L(t-1)*x1)*lz_u/2;            
                        kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_u/x2-1/2*ke_d/x2;
                        temp = - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_u/y2-1/2*kn_d/y2;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_u/x1-1/2*kw_d/x1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_u/y1-1/2*ks_d/y1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/lz_u;
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;
                        pointer=pointer+1;                          
                    else
                        if j == yb
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            kn_u = (Cdt_L(t)*x2+Cdt_L(t-1)*x1)*lz_u/2;
                            kw_u = Cdt_L(t-1)*y2*lz_u;ke_u = Cdt_L(t)*y2*lz_u;
                            kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t-1)*x1*y2)/2;                       

                            P(id)=1/2*(Ta*h*lx*y1+Ta*lx*lz_u);
                            temp = P(id)/Ta;                              
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_u/(4*x2)-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_u/(2*y2)-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_u/(4*x1)-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(2*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;          
                            pointer=pointer+1;                  

                        elseif j == yt
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            ks_u = (Cdt_L(t_s)*x2+Cdt_L(t_s-1)*x1)*lz_u/2;
                            kw_u = Cdt_L(t_s-1)*y1*lz_u;ke_u = Cdt_L(t_s)*y1*lz_u;
                            kz_u = (Cdt_V(t_s)*x2*y1+Cdt_V(t_s-1)*x1*y1)/2;                        

                            P(id)=1/2*(Ta*h*lx*y2+Ta*lx*lz_u);
                            temp = P(id)/Ta;                            
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_u/(4*x2)-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_u/(4*x1)-kw_d/(2*x1);
                            temp = temp - value(pointer); 
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_u/(2*y1)-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(2*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;        
                            pointer=pointer+1;                  

                        elseif i == xl
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            ke_u = (Cdt_L(t)*y2+Cdt_L(t_s)*y1)*lz_u/2;
                            kn_u = Cdt_L(t)*x2*lz_u;ks_u = Cdt_L(t_s)*x2*lz_u;
                            kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t_s)*x2*y1)/2;                        

                            P(id)=1/2*(Ta*h*ly*x1+Ta*ly*lz_u);
                            temp = P(id)/Ta;                             
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/(2*x2)-ke_u/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_u/(4*y2)-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_u/(4*y1)-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(2*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;          
                            pointer=pointer+1;                                           

                        else
                            kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            kw_u = (Cdt_L(t-1)*y2+Cdt_L(t_s-1)*y1)*lz_u/2;
                            kn_u = Cdt_L(t-1)*x1*lz_u;ks_u = Cdt_L(t_s-1)*x1*lz_u;
                            kz_u = (Cdt_V(t-1)*x1*y2+Cdt_V(t_s-1)*x1*y1)/2;                        

                            P(id)=1/2*(Ta*h*ly*x2+Ta*ly*lz_u);
                            temp = P(id)/Ta;                            
                            
                            row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/(2*x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-kn_u/(4*y2)-kn_d/(2*y2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-kw_u/(2*x1)-kw_d/(2*x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-ks_u/(4*y1)-ks_d/(2*y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/(2*lz_u);            
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;             
                            pointer=pointer+1;  
                        end
                    end
                else
                        kw_d = (Cdt_L(s-1)*y1+Cdt_L(id-1)*y2)*lz_d/2;
                        ks_d = (Cdt_L(s-1)*x1+Cdt_L(s)*x2)*lz_d/2;
                        ke_d = (Cdt_L(s)*y1+Cdt_L(id)*y2)*lz_d/2;
                        kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;            
                        kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;                    

                        P(id)=Ta*h*lx*ly;
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-ke_d/(2*x2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-kn_d/(2*y2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-kw_d/(2*x1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-ks_d/(2*y1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=b; value(pointer)=-kz_d/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;              
                        pointer=pointer+1;                
                end        
            else
                if i == 1
                        ke_d = (Cdt_L(id)*y2+Cdt_L(s)*y1)*lz_d/2;
                        kn_d = Cdt_L(id)*x2*lz_d;ks_d = Cdt_L(s)*x2*lz_d;
                        kz_d = (Cdt_V(id)*x2*y2+Cdt_V(s)*x2*y1)/2;

                        P(id)=1/2*(Ta*h*ly*lz_d+Ta*h*x2*ly);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_d/(2*y1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_d/(2*y2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_d/(x2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;              
                        pointer=pointer+1;                    

                elseif i == gridNx_pack
                        kw_d = (Cdt_L(id-1)*y2+Cdt_L(s-1)*y1)*lz_d/2;
                        kn_d = Cdt_L(id-1)*x1*lz_d;ks_d = Cdt_L(s-1)*x1*lz_d;
                        kz_d = (Cdt_V(id-1)*x1*y2+Cdt_V(s-1)*x1*y1)/2;                

                        P(id)=1/2*(Ta*h*ly*lz_d+Ta*h*x1*ly);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_d/(2*y1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_d/(x1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_d/(2*y2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;              
                        pointer=pointer+1;                          

                elseif j == gridNy_pack
                        ks_d = (Cdt_L(s)*x2+Cdt_L(s-1)*x1)*lz_d/2;
                        kw_d = Cdt_L(s-1)*y1*lz_d;ke_d = Cdt_L(s)*y1*lz_d;
                        kz_d = (Cdt_V(s)*x2*y1+Cdt_V(s-1)*x1*y1)/2;                

                        P(id)=1/2*(Ta*h*lx*lz_d+Ta*h*lx*y1);
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_d/(2*x2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_d/(2*x1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_d/y1;
                        temp = temp - value(pointer);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;                   

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;              
                        pointer=pointer+1;                                              

                else
                        kn_d = (Cdt_L(id)*x2+Cdt_L(id-1)*x1)*lz_d/2;
                        kw_d = Cdt_L(id-1)*y2*lz_d;ke_d = Cdt_L(id)*y2*lz_d;
                        kz_d = (Cdt_V(id)*x2*y2+Cdt_V(id-1)*x1*y2)/2;                

                        P(id)=1/2*(Ta*h*lx*lz_d+Ta*h*lx*y2);  
                        temp = P(id)/Ta;

                        row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_d/(2*x2);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_d/(2*x1);
                        temp = temp - value(pointer);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_d/y2;
                        temp = temp - value(pointer);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer)=-1/2*kz_d/lz_d;
                        temp = temp - value(pointer);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer)=temp;              
                        pointer=pointer+1;                                          
                end
            end       
        end
    end
    disp(['finish: layer' num2str(layer_id)]);
    toc;
    %bottom of the board
    tic;
    for layer_id = Layer.N+1
        disp(['analyze: layer' num2str(layer_id)]);
        for i = 1:1:gridNx_pack
            for j = 1:1:gridNy_pack
                if(layer_id < Layer.N)
                    lz_d = Thick_layer(layer_id);
                else
                    lz_d = 0;
                end
                
                lz_u = Thick_layer(layer_id-1);
                lz = (lz_u+lz_d)/2;
                
                if(i>1) 
                    x1 = pack_xmesh(i)-pack_xmesh(i-1);
                else
                    x1 = 0;
                end

                if(i<gridNx_pack)
                    x2 = pack_xmesh(i+1) - pack_xmesh(i);
                else
                    x2 = 0;
                end

                % y direction
                if(j>1) 
                    y1 = pack_ymesh(j)-pack_ymesh(j-1);
                else
                    y1 = 0;
                end

                if(j<gridNy_pack)
                    y2 = pack_ymesh(j+1) - pack_ymesh(j);
                else
                    y2 = 0;
                end

                lx = (x1+x2)/2;
                ly = (y1+y2)/2;
                
                const = gridNx_chip*gridNy_chip*(Layer.N-1)+gridNx_pack*gridNy_pack;
                id = i+(j-1)*gridNx_pack+const;
                w = id-1;
                e = id+1;
                s = id - gridNx_pack;
                n = id + gridNx_pack;

                t = id-gridNx_pack*gridNy_pack;
                t_s = t - gridNx_pack; 

                if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ) )
                    if i == 1 && j == 1
                            ke_u = Cdt_L(t);kn_u = Cdt_L(t);kz_u = Cdt_V(t);
                            P(id)=1/4*(Ta*h*y2*lz_u+Ta*h*x2*lz_u+Ta*h_d*x2*y2);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u*y2*lz_u/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u*x2*lz_u/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u*x2*y2/lz_u;         
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;          
                            pointer=pointer+1;                    

                    elseif i == gridNx_pack && j == 1
                            kw_u = Cdt_L(t-1);kn_u = Cdt_L(t-1);kz_u = Cdt_V(t-1);
                            P(id)=1/4*(Ta*h*y2*lz_u+Ta*h*x1*lz_u+Ta*h_d*x1*y2);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u*y2*lz_u/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u*x1*lz_u/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u*x1*y2/lz_u;         
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;          
                            pointer=pointer+1;                      

                    elseif i == 1 && j == gridNy_pack
                            ke_u = Cdt_L(t_s);ks_u = Cdt_L(t_s);kz_u = Cdt_V(t_s);
                            P(id)=1/4*(Ta*h*y1*lz_u+Ta*h*x2*lz_u+Ta*h_d*x2*y1);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u*y1*lz_u/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u*x2*lz_u/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*kz_u*x2*y1/lz_u;     
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;           
                            pointer=pointer+1;                     

                    else
                            kw_u = Cdt_L(t_s-1);ks_u = Cdt_L(t_s-1);kz_u = Cdt_V(t_s-1);
                            P(id)=1/4*(Ta*h*y1*lz_u+Ta*h*x1*lz_u+Ta*h_d*x1*y1);
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*(kw_u*y1*lz_u/x1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*(ks_u*x1*lz_u/y1);
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/4*(kz_u*x1*y1/lz_u);
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;          
                            pointer=pointer+1;                        

                    end
                elseif i ~= 1 && i~= gridNx_pack && j ~= 1 && j~=gridNy_pack
                    kw_u = (Cdt_L(t_s-1)*y1+Cdt_L(t-1)*y2)*lz_u/2;
                    ks_u = (Cdt_L(t_s-1)*x1+Cdt_L(t_s)*x2)*lz_u/2;
                    ke_u = (Cdt_L(t_s)*y1+Cdt_L(t)*y2)*lz_u/2;
                    kn_u = (Cdt_L(t)*x2+Cdt_L(t-1)*x1)*lz_u/2;            
                    kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;                

                    if i >= xl && i<= xr && j>=yb && j<=yt                   
                        P(id)=Ta*h_mf*lx*ly;
                        temp = P(id)/Ta;
                    else
                        P(id)=Ta*h_d*lx*ly;
                        temp = P(id)/Ta; 
                    end                    

                    row(pointer)=id; column(pointer)=e; value(pointer)=-ke_u/(2*x2);
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=n; value(pointer)=-kn_u/(2*y2);
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=w; value(pointer)=-kw_u/(2*x1);
                    temp = temp - value(pointer);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=s; value(pointer)=-ks_u/(2*y1);
                    temp = temp - value(pointer);
                    pointer=pointer+1;           

                    row(pointer)=id; column(pointer)=t; value(pointer)=-kz_u/lz_u;
                    temp = temp - value(pointer);
                    pointer=pointer+1;           

                    row(pointer)=id; column(pointer)=id; value(pointer)=temp;              
                    pointer=pointer+1;  

                else
                    if i == 1
                            ke_u = (Cdt_L(t_s)*y1+Cdt_L(t)*y2)*lz_u/2;
                            kn_u = Cdt_L(t)*x2*lz_u;ks_u = Cdt_L(t_s)*x2*lz_u;
                            kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t_s)*x2*y1)/2;                    

                            P(id)=Ta*h*ly*lz+Ta*h_d*lx*ly;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/2*ke_u/(x2);
                            temp = temp - value(pointer);
                            pointer=pointer+1;         

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;             
                            pointer=pointer+1;                      

                    elseif i == gridNx_pack
                            kw_u = (Cdt_L(t-1)*y2+Cdt_L(t_s-1)*y1)*lz_u/2;
                            kn_u = Cdt_L(t-1)*x1*lz_u;ks_u = Cdt_L(t_s-1)*x1*lz_u;
                            kz_u = (Cdt_V(t-1)*x1*y2+Cdt_V(t_s-1)*x1*y1)/2;                     

                            P(id)=Ta*h*ly*lz+Ta*h_d*lx*ly;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/4*ks_u/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/2*kw_u/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/4*kn_u/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;         

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;             
                            pointer=pointer+1;                       

                    elseif j == gridNy_pack
                            ks_u = (Cdt_L(t_s)*x2+Cdt_L(t_s-1)*x1)*lz_u/2;
                            kw_u = Cdt_L(t_s-1)*y1*lz_u;ke_u = Cdt_L(t_s-1)*y1*lz_u;
                            kz_u = (Cdt_V(t_s)*x2*y1+Cdt_V(t_s-1)*x1*y1)/2;                    

                            P(id)=Ta*h*lx*lz+Ta*h_d*lx*ly;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer)=-1/2*ks_u/y1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;         

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;             
                            pointer=pointer+1;                        

                    else
                            kn_u = (Cdt_L(t)*x2+Cdt_L(t-1)*x1)*lz_u/2;
                            kw_u = Cdt_L(t-1)*y2*lz_u;ke_u = Cdt_L(t)*y2*lz_u;
                            kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t-1)*x1*y2)/2;                    

                            P(id)=Ta*h*lx*lz+Ta*h_d*lx*ly;
                            temp = P(id)/Ta;

                            row(pointer)=id; column(pointer)=e; value(pointer)=-1/4*ke_u/x2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer)=-1/4*kw_u/x1;
                            temp = temp - value(pointer);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer)=-1/2*kn_u/y2;
                            temp = temp - value(pointer);
                            pointer=pointer+1;         

                            row(pointer)=id; column(pointer)=t; value(pointer)=-1/2*kz_u/lz_u;
                            temp = temp - value(pointer);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer)=temp;             
                            pointer=pointer+1;                      

                    end
                end
            end
        end
        disp(['finish: layer' num2str(layer_id)]);   
    end
    toc;
    row = row(1:pointer-1,1);
    column = column(1:pointer-1,1);
    value = value(1:pointer-1,1);
    A = sparse(row, column, value, var, var);
    Pa = P;
end

