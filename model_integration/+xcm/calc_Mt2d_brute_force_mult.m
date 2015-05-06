function [Mt2dc, Mt2d, term2, term3] = calc_Mt2d_brute_force_mult(Lx, Nuc_1d, w_tsv)


w_uc = Lx/Nuc_1d;
forbidden_min = round(w_uc/2 - w_tsv/2);
forbidden_max = forbidden_min + w_tsv-1;

T = Lx/Nuc_1d;
t = forbidden_min-1;

lmax = 2*Lx;
term2 = zeros(1,lmax+1);
term3 = zeros(1,lmax+1);
Mt2dc = zeros(1,lmax+1);
Mt2d = zeros(1,lmax+1);
Ly = Lx;
Nx = Lx;
Ny = Ly;

w = w_tsv;

for l_ind = 1:lmax+1
    l = l_ind - 1
    lxmax = min(l,Lx);
    lxmin = max(0,l-Nx+1);

    for lx_ind = (lxmin+1):(lxmax+1)
        lx = lx_ind - 1;
        ly = l - lx;
        ly_ind = ly + 1;

        for x1_ind = 1:Nx
            x1 = x1_ind - 1;
            x1_uc = mod(x1,T);
            gx = (x1_uc >= forbidden_min) && (x1_uc <= forbidden_max);
            x2 = x1 - lx;
            x2_uc = mod(x2,T);
            gx2 = (x2_uc >= forbidden_min) && (x2_uc <= forbidden_max) && (x2 >= 0);
            
            for y1_ind = 1:Ny
                y1 = y1_ind - 1;
                y1_uc = mod(y1,T);
                gy = (y1_uc >= forbidden_min) && (y1_uc <= forbidden_max);

                y2 = y1 - ly;
                y2_uc = mod(y2,T);
                gy2 = (y2_uc >= forbidden_min) && (y2_uc <= forbidden_max) && (y2 >=0);
                
                p2ok = (x2 >= 0 ) && (y2 >= 0);


                term2(l_ind) = term2(l_ind) + gx*gy;
                term3(l_ind) = term3(l_ind) + gx2*gy2;
                Mt2dc(l_ind) = Mt2dc(l_ind)+(1-gx*gy)*( (p2ok)*(1-gx2*gy2));
                Mt2d(l_ind) = Mt2d(l_ind) + p2ok;
            end
        end
    end
end



%% If we do this one we get the same result as Mt2dbfc, obviously
% for x1_ind = 1:Nx
%     x1 = x1_ind - 1
%     x1_uc = mod(x1,T);
%     gx = (x1_uc >= forbidden_min) && (x1_uc <= forbidden_max);
%     
% 
%     for y1_ind = 1:Ny
%         y1 = y1_ind - 1;
%         y1_uc = mod(y1,T);
%         gy = (y1_uc >= forbidden_min) && (y1_uc <= forbidden_max);
%         
%         for x2_ind = 1:Nx
%             x2 = x2_ind-1;
%             x2_uc = mod(x2,T);
%             gx2 = (x2_uc >= forbidden_min) && (x2_uc <= forbidden_max) && (x2 >= 0);
%             
%             for y2_ind = y1_ind:Ny
%                 y2 = y2_ind-1;
%                 y2_uc = mod(y2,T);
%                 gy2 = (y2_uc >= forbidden_min) && (y2_uc <= forbidden_max) && (y2 >=0);
% 
%                 p2ok = (x2 >= 0 ) && (y2 >= 0);
%                 
%                 l = abs(x1_ind - x2_ind) + abs(y1_ind - y2_ind);
%                 l_ind = l + 1;
% 
%                 same_row = (y2_ind <= y1_ind);
%                 before_x1 = (x2_ind < x1_ind);
%                 already_counted = (same_row && before_x1);
%                 
%                 if (~already_counted)
%                     term2(l_ind) = term2(l_ind) + gx*gy;
%                     term3(l_ind) = term3(l_ind) + gx2*gy2;
%                     Mt2dc(l_ind) = Mt2dc(l_ind)+(1-gx*gy)*( (p2ok)*(1-gx2*gy2));
%                     Mt2d(l_ind) = Mt2d(l_ind) + p2ok;
%                 end
%             end
%         end
%     end
% end


