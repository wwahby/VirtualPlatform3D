function [C, Y] = Matrix_build(chip, system, TSV, row, column, value, var, type)
    %if type is 1, simulate the power nodes
    %if type is 2, simulate the gnd nodes.
    %clear; 
    %clc;
    %this is the core function of the solver
    %it depends on the schemes we are using
    %for a general node which is connected with a pad, the equation is:
    % [x(i,j)-x(i-1,j)]/Rx + [x(i,j)-x(i+1,j)]/Rx +[x(i,j)-x(i,j+1)]/Ry +
    % [x(i,j)-x(i,j-1)]/Ry + c*delta_V(i,j)/dt + current(i,j) = I_pad(i,j);
    % the scheme we are using is Trapezoid scheme (TR)
    % what we do is, use (x(i,j,n)+x(i,j,n-1))/2 replace x(i,j) (same to I_pad)
    % delta_V(i,j) = V(i,j,n)-V(i,j,n-1);
    % for the above scheme, we will see the final equation is:
    % [x(i,j,n)-x(i-1,j,n)+(n-1 term)]/(2*Rx) + [x(i,j,n)-x(i+1,j,n)+(n-1 term)]/(2*Rx) +[x(i,j,n)-x(i,j+1,n)+(n-1 term)]/(2*Ry) +
    % [x(i,j,n)-x(i,j-1, n)+(n-1 term)]/(2*Ry) + c*(V(i,j,n)-V(i,j,n-1)/dt + (current(i,j,n)+current(i,j,n-1))/2  = (I_pad(i,j,n)+I_pad(i,j,n-1))/2;
    % move the n-1 term to the right side, leave the n term to the left
    % side, using matrix to rewrite the equation it is
    % ( -1/2Rx -1/2Rx 1/Rx+1/Ry -1/2Ry -1/2Ry)*X' + (0 0 C 0 0)/dt*X' = 
    % -(-1/2Rx -1/2Rx 1/Rx+1/Ry -1/2Ry -1/2Ry)'* X(n-1)'+ (0 0 C 0 0)/dt*X(n-1)'-(current(i,j,n)+current(i,j,n-1))/2
    % that is Y*X' + C*X' = (C-Y)*X(n-1) + Current
    % for each node, we should build such correlation Matrix, this is what
    % this funtion does
    
    %for each PAD node, the equation is:
    % Vdd/2 - Lp*di/dt - i*Rp = V
    % discretize the equation, we get:
    %Vdd/2 = Lp*(I(i,j,n)-I(i,j,n-1))/dt - (I(i,j,n)+I(i,j,n-1))/2*Rp =
    %(V(i,j,n)+V(i,j,n-1))/2
    %build the same matrix as above
    
%%%%%%%%%%%%%%%%%%%%%%%%simulation mode%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
    gridNx = chip.Nx;
    gridNy = chip.Ny;
    const = gridNx*gridNy*chip.N;
    chip_xmesh = chip.Xmesh;
    chip_ymesh = chip.Ymesh;
    
    rv = TSV.R;
    lv = TSV.L;
    C_mat = chip.c;
    L = system.L;
    rp = system.R;
    if type == 1
        tsvN = TSV.P;
    else
        tsvN = TSV.G;
    end
    
%%%%assign the capacitance inductance
    P_tsv = 0;
    %pivoting the number of TSV, for referencing the current
    pointer = 1;
    for i=1:1:gridNx
        for j=1:1:gridNy
            if (chip.type(j,i) == type)
                P_tsv = P_tsv + 1;
                id = const + P_tsv;
                t = (j-1)*gridNx+i;
                
                row(pointer) = id;
                column(pointer) = id;
                value(pointer, 1) = 0.5*rp;
                value(pointer, 2) = L;
                pointer = pointer+1;              
                
                row(pointer) = id;
                column(pointer) = t;
                value(pointer, 1) = 0.5;
                value(pointer, 2) = 0;
                pointer = pointer+1;  
            end            
        end
    end
    
    for k = 2:1:chip.N
        P_tsv = 0;
        for i=1:1:gridNx
            for j=1:1:gridNy
                if (chip.type(j,i) == type)
                    P_tsv = P_tsv + 1;
                    id = const + P_tsv + (k-1)*tsvN;
                    t = (j-1)*gridNx+i + (k-1)*gridNx*gridNy;
                    b = (j-1)*gridNx+i + (k-2)*gridNx*gridNy;

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) = 0.5*rv;
                    value(pointer, 2) = lv;
                    pointer = pointer+1;              

                    row(pointer) = id;
                    column(pointer) = t;
                    value(pointer, 1) = 0.5;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    
                    row(pointer) = id;
                    column(pointer) = b;
                    value(pointer, 1) = -0.5;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                    
                end            
            end
        end        
    end
      
    P_tsv = 0;
    Layer_id = 1;
    for i=1:1:gridNx
        for j=1:1:gridNy
            if(i>1) 
                x1 = chip_xmesh(i)-chip_xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<gridNx)
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

            if(j<gridNy)
                y2 = chip_ymesh(j+1) - chip_ymesh(j);
            else
                y2 = 0;
            end

            gridx = (x1+x2)/2;
            gridy = (y1+y2)/2;
            
            rx1 = chip.Rs*(x1/chip.wire_width);
            rx2 = chip.Rs*(x2/chip.wire_width);
            ry1 = chip.Rs*(y1/chip.wire_width);
            ry2 = chip.Rs*(y2/chip.wire_width);
            
            id = (j-1)*gridNx+i+(Layer_id-1)*gridNx*gridNy;
            w = id-1;
            e = id+1;
            n = id+gridNx;
            s = id-gridNx;
            c = C_mat(id);
            
            if (chip.type(j,i) == type)
                P_tsv = P_tsv + 1;
                b = const + P_tsv;
                t = const + P_tsv + tsvN;
            end                                                
                        
            if (i==1 || i==gridNx) && (j==1 || j==gridNy)
                if i==1 && j==1
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = 0.5/rx2;
                    
                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1; 
                    temp = temp+0.5/ry2;
               
                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
               
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) = temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;
                    
                elseif i==1 && j==gridNy
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = 0.5/rx2;
           
                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = temp+0.5/ry1;

                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1; 
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) = temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;               
                elseif i==gridNx && j==1
                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = 0.5/ry2;

                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = temp + 0.5/rx1;

                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) = temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;
                else
                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = 0.5/ry1;
                    
                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    temp = temp + 0.5/rx1;

                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end;                       
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) = temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;
                end
            elseif i~=1 && i~=gridNx && j~=gridNy && j~=1
                row(pointer) = id;
                column(pointer) = e;
                value(pointer, 1) = -0.5/rx2;
                value(pointer, 2) = 0;
                pointer = pointer+1;

                row(pointer) = id;
                column(pointer) = n;
                value(pointer, 1) = -0.5/ry2;
                value(pointer, 2) = 0;
                pointer = pointer+1;                  

                row(pointer) = id;
                column(pointer) = s;
                value(pointer, 1) = -0.5/ry1;
                value(pointer, 2) = 0;
                pointer = pointer+1;

                row(pointer) = id;
                column(pointer) = w;
                value(pointer, 1) = -0.5/rx1;
                value(pointer, 2) = 0;
                pointer = pointer+1;
                
                temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry1 + 0.5/ry2;
                            
                if (chip.type(j,i) == type)
                    row(pointer) = id;
                    column(pointer) = b;
                    value(pointer, 1) = -1/2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    
                    if chip.N > 1
                        row(pointer) = id;
                        column(pointer) = t;
                        value(pointer, 1) = 1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;    
                    end
                end
                
                row(pointer) = id;
                column(pointer) = id;
                value(pointer, 1) =  temp;
                value(pointer, 2) = 2*c*gridx*gridy;
                pointer = pointer+1;
            else
                if i==1
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                  

                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    
                    temp = 0.5/rx2 + 0.5/ry1 + 0.5/ry2;
                              
                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) =  temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;   
                elseif j==gridNy
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                 

                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    
                    temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry1;
              
                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) =  temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;           
                elseif i==gridNx
                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                  

                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    
                    temp = 0.5/rx1 + 0.5/ry1 + 0.5/ry2;
               
                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                                            
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) = temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;
                else
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                  

                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;
                    
                    temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry2;
                             
                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -1/2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        
                        if chip.N > 1
                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;    
                        end
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) =  temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;    
                end
            end
        end
    end
    
    for Layer_id = 2:1:chip.N-1;
        P_tsv = 0;
        for i=1:1:gridNx
            for j=1:1:gridNy
                if(i>1) 
                    x1 = chip_xmesh(i)-chip_xmesh(i-1);
                else
                    x1 = 0;
                end

                if(i<gridNx)
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

                if(j<gridNy)
                    y2 = chip_ymesh(j+1) - chip_ymesh(j);
                else
                    y2 = 0;
                end

                gridx = (x1+x2)/2;
                gridy = (y1+y2)/2;

                rx1 = chip.Rs*(x1/chip.wire_width);
                rx2 = chip.Rs*(x2/chip.wire_width);
                ry1 = chip.Rs*(y1/chip.wire_width);
                ry2 = chip.Rs*(y2/chip.wire_width);

                id = (j-1)*gridNx+i+(Layer_id-1)*gridNx*gridNy;
                w = id-1;
                e = id+1;
                n = id+gridNx;
                s = id-gridNx;
                c = C_mat(id);
                
                if (chip.type(j,i) == type)
                    P_tsv = P_tsv + 1;
                    t = const + P_tsv + tsvN*(Layer_id);
                    b = const + P_tsv + tsvN*(Layer_id-1);
                end                                                

                if (i==1 || i==gridNx) && (j==1 || j==gridNy)
                    if i==1 && j==1
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/rx2;

                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1; 
                        temp = temp+0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;

                    elseif i==1 && j==gridNy
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/rx2;

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = temp+0.5/ry1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 1/2;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;                      
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;               
                    elseif i==gridNx && j==1
                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/ry2;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = temp + 0.5/rx1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;                      
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;
                    else
                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/ry1;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = temp + 0.5/rx1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;                      
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;
                    end
                elseif i~=1 && i~=gridNx && j~=gridNy && j~=1
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                  

                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry1 + 0.5/ry2;

                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -0.5;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = t;
                        value(pointer, 1) = 0.5;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) =  temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;
                else
                    if i==1
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                  

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx2 + 0.5/ry1 + 0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) =  temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;   
                    elseif j==gridNy
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                 

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) =  temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;           
                    elseif i==gridNx
                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                  

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx1 + 0.5/ry1 + 0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;
                    else
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                  

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;

                            row(pointer) = id;
                            column(pointer) = t;
                            value(pointer, 1) = 0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) =  temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;    
                    end
                end
            end
        end
    end
    
    P_tsv = 0;
    if chip.N > 1
        Layer_id = chip.N;
        for i=1:1:gridNx
            for j=1:1:gridNy
                if(i>1) 
                    x1 = chip_xmesh(i)-chip_xmesh(i-1);
                else
                    x1 = 0;
                end

                if(i<gridNx)
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

                if(j<gridNy)
                    y2 = chip_ymesh(j+1) - chip_ymesh(j);
                else
                    y2 = 0;
                end

                gridx = (x1+x2)/2;
                gridy = (y1+y2)/2;

                rx1 = chip.Rs*(x1/chip.wire_width);
                rx2 = chip.Rs*(x2/chip.wire_width);
                ry1 = chip.Rs*(y1/chip.wire_width);
                ry2 = chip.Rs*(y2/chip.wire_width);

                id = (j-1)*gridNx+i+(Layer_id-1)*gridNx*gridNy;
                w = id-1;
                e = id+1;
                n = id+gridNx;
                s = id-gridNx;
                c = C_mat(id);
                
                if (chip.type(j,i) == type)
                    P_tsv = P_tsv + 1;
                    b = const + P_tsv + tsvN*(Layer_id-1);
                end                                                

                if (i==1 || i==gridNx) && (j==1 || j==gridNy)
                    if i==1 && j==1
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/rx2;

                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1; 
                        temp = temp+0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;

                    elseif i==1 && j==gridNy
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/rx2;

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = temp+0.5/ry1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;                       
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;               
                    elseif i==gridNx && j==1
                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/ry2;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = temp + 0.5/rx1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;                     
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;
                    else
                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = 0.5/ry1;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                        temp = temp + 0.5/rx1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;
                    end
                elseif i~=1 && i~=gridNx && j~=gridNy && j~=1
                    row(pointer) = id;
                    column(pointer) = e;
                    value(pointer, 1) = -0.5/rx2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = n;
                    value(pointer, 1) = -0.5/ry2;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;                  

                    row(pointer) = id;
                    column(pointer) = s;
                    value(pointer, 1) = -0.5/ry1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    row(pointer) = id;
                    column(pointer) = w;
                    value(pointer, 1) = -0.5/rx1;
                    value(pointer, 2) = 0;
                    pointer = pointer+1;

                    temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry1 + 0.5/ry2;

                    if (chip.type(j,i) == type)
                        row(pointer) = id;
                        column(pointer) = b;
                        value(pointer, 1) = -0.5;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;
                    end

                    row(pointer) = id;
                    column(pointer) = id;
                    value(pointer, 1) =  temp;
                    value(pointer, 2) = 2*c*gridx*gridy;
                    pointer = pointer+1;
                else
                    if i==1
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                  

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx2 + 0.5/ry1 + 0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) =  temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;   
                    elseif j==gridNy
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                 

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry1;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) =  temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;           
                    elseif i==gridNx
                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                  

                        row(pointer) = id;
                        column(pointer) = s;
                        value(pointer, 1) = -0.5/ry1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx1 + 0.5/ry1 + 0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) = temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;
                    else
                        row(pointer) = id;
                        column(pointer) = e;
                        value(pointer, 1) = -0.5/rx2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        row(pointer) = id;
                        column(pointer) = n;
                        value(pointer, 1) = -0.5/ry2;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;                  

                        row(pointer) = id;
                        column(pointer) = w;
                        value(pointer, 1) = -0.5/rx1;
                        value(pointer, 2) = 0;
                        pointer = pointer+1;

                        temp = 0.5/rx1 + 0.5/rx2 + 0.5/ry2;

                        if (chip.type(j,i) == type)
                            row(pointer) = id;
                            column(pointer) = b;
                            value(pointer, 1) = -0.5;
                            value(pointer, 2) = 0;
                            pointer = pointer+1;
                        end

                        row(pointer) = id;
                        column(pointer) = id;
                        value(pointer, 1) =  temp;
                        value(pointer, 2) = 2*c*gridx*gridy;
                        pointer = pointer+1;    
                    end
                end
            end
        end
    end
    
    row = row(1:pointer-1,1);
    column = column(1:pointer-1,1);
    value = value(1:pointer-1,:);
    Y = sparse(row, column, value(:,1), var, var);
    C = sparse(row, column, value(:,2), var, var);
end

