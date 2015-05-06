function [term2, term3] = calc_mid_terms_brute_force(Lx, Nuc_1d, w_tsv)


w_uc = Lx/Nuc_1d;
forbidden_min = round(w_uc/2 - w_tsv/2);
forbidden_max = forbidden_min + w_tsv-1;

T = Lx/Nuc_1d;
t = forbidden_min-1;

r = @(x,start,width) (x >= start) && (x <= start + width);


fx = @(x,Lx) r(x,0,Lx-1);
q = @(x_uc,y_uc,t,w) r(x_uc,t,w) * r(y_uc,t,w);


lmax = 2*Lx;
term2 = zeros(1,lmax+1);
term3 = zeros(1,lmax+1);
Ly = Lx;
Nx = Lx;
Ny = Ly;

w = w_tsv;
% for l_ind = 1:lmax+1
%     l = l_ind - 1
%     lxmax = min(l,Lx);
%     for lx_ind = 1:lxmax+1
%         lx = lx_ind - 1;
%         for x1_ind = 1:Nx
%             x1 = x1_ind - 1;
%             x1_uc = mod(x1,T);
%             x2 = x1 - lx;
%             x2_uc = mod(x2,T);
%             x2_nogo = x2 < 0;
%             if (~x2_nogo)
%                 for y1_ind = 1:Ny
%                     y1 = y1_ind - 1;
%                     y1_uc = mod(y1,T);
%                     ly = l - lx;
%                     y2 = y1 - ly;
%                     y2_uc = mod(y2,T);
%                     y2_nogo = y2 < 0;
%                     
%                     p2_nogo = x2_nogo && y2_nogo;
% 
%                     if(~p2_nogo)
%                         %a1 = fx(x1,Lx)*fx(y1,Ly); % guaranteed to be fine
%                         %a2 = fx(x2,Lx)*fx(y2,Ly); % guaranteed to be fine
% 
%                         %b1 = q(x1_uc,y1_uc,t,w_tsv);
%                         b1x = (x1_uc < t) && (x1_uc > t + w);
%                         b1y = (y1_uc < t) && (y1_uc > t + w);
%                         b1 = b1x*b1y;
%                         
%                         b2x = (x2_uc < t) && (x2_uc > t + w);
%                         b2y = (y2_uc < t) && (y2_uc > t + w);
%                         b2 = b2x*b2y;
%                         %b2 = q(x2_uc,y2_uc,t,w_tsv);
% 
%                         term2(l_ind) = term2(l_ind) + b2;
%                         term3(l_ind) = term3(l_ind) + b1;
%                     end
%                 end
%             end
%         end
%     end
% end
 

for l_ind = 1:lmax+1
    l = l_ind - 1;
    lxmax = min(l,Lx);
    lxmin = max(0,l-Nx+1);

    for lx_ind = (lxmin+1):(lxmax+1)
        lx = lx_ind - 1;
        ly = l - lx;
        ly_ind = ly + 1;

        for x1_ind = lx_ind:Nx
            x1 = x1_ind - 1;
            x1_uc = mod(x1,T);
            gx = (x1_uc >= forbidden_min) && (x1_uc <= forbidden_max);
            x2 = x1 - lx;
            x2_uc = mod(x2,T);
            gx2 = (x2_uc >= forbidden_min) && (x2_uc <= forbidden_max) && (x2 >= 0);
            
            if ((gx ~= 0) || (gx2 ~=0))
                for y1_ind = ly_ind:Ny
                    y1 = y1_ind - 1;
                    y1_uc = mod(y1,T);
                    gy = (y1_uc >= forbidden_min) && (y1_uc <= forbidden_max);
                    
                    y2 = y1 - ly;
                    y2_uc = mod(y2,T);
                    gy2 = (y2_uc >= forbidden_min) && (y2_uc <= forbidden_max) && (y2 >=0);
                    
                    
                    term2(l_ind) = term2(l_ind) + gx*gy;
                    term3(l_ind) = term3(l_ind) + gx2*gy2;
                end
            end
        end
    end
end
