function [Mt2d, SFXC, DFXC] = Mt2d_brute_force_corrected(Lx, Nuc_1d, w_tsv)


w_uc = Lx/Nuc_1d;

forbidden_min = round(w_uc/2) - round(w_tsv/2);
forbidden_max = forbidden_min + (w_tsv - round(w_tsv/2));

Mt2d = zeros(1,2*Lx+1);
SFXC = zeros(1,2*Lx+1);
DFXC = zeros(1,2*Lx+1);
for x1_ind = 1:Lx
    x1_uc = mod(x1_ind,w_uc);
    if (x1_uc == 0)
        x1_uc = w_uc;
    end
    
    x1_nogo = (x1_uc >= forbidden_min) && (x1_uc <= forbidden_max);
    
    for y1_ind = 1:Lx
        y1_uc = mod(y1_ind,w_uc);
        if (y1_uc == 0)
            y1_uc = w_uc;
        end
        
        y1_nogo = (y1_uc >= forbidden_min) && (y1_uc <= forbidden_max);
        
        p1_nogo = x1_nogo && y1_nogo;
        
        %if (~p1_nogo)
            for x2_ind = 1:Lx
                x2_uc = mod(x2_ind,w_uc);
                if (x2_uc == 0)
                    x2_uc = w_uc;
                end

                x2_nogo = (x2_uc >= forbidden_min) && (x2_uc <= forbidden_max);

                for y2_ind = y1_ind:Lx
                    y2_uc = mod(y2_ind,w_uc);
                    if (y2_uc == 0)
                        y2_uc = w_uc;
                    end

                    y2_nogo = (y2_uc >= forbidden_min) && (y2_uc <= forbidden_max);
                    p2_nogo = x2_nogo && y2_nogo;

                    %if (~p2_nogo)
                        same_row = (y2_ind == y1_ind);
                        before_x1 = (x2_ind < x1_ind);
                        already_counted = (same_row && before_x1);
                        if(~already_counted)
                            l = abs(x1_ind - x2_ind) + abs(y1_ind - y2_ind);
                            if (p2_nogo && p1_nogo)
                                DFXC(l+1) = DFXC(l+1) + 1;
                            elseif ((p2_nogo) || (p1_nogo))
                                SFXC(l+1) = SFXC(l+1) + 1;
                            else
                                Mt2d(l+1) = Mt2d(l+1) + 1;
                            end
                        end
                    %end
                end
            end
        %end
    end
end
