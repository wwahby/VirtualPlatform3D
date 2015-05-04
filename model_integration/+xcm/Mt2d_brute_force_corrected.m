function [Mt2d, SFXC, DFXC] = Mt2d_brute_force_corrected(Lx, Nuc_1d, w_tsv)


w_uc = Lx/Nuc_1d;

forbidden_min = round(w_uc/2 - w_tsv/2);
forbidden_max = forbidden_min + w_tsv-1;

Mt2d = zeros(1,2*Lx+1);
SFXC = zeros(1,2*Lx+1);
DFXC = zeros(1,2*Lx+1);
for x1_ind = 1:Lx
    x1 = x1_ind - 1;
    x1_uc = mod(x1,w_uc);
    x1_nogo = (x1_uc >= forbidden_min) && (x1_uc <= forbidden_max);
    
    for y1_ind = 1:Lx
        y1 = y1_ind - 1;
        y1_uc = mod(y1,w_uc);
        y1_nogo = (y1_uc >= forbidden_min) && (y1_uc <= forbidden_max);
        
        p1_nogo = x1_nogo && y1_nogo;
        
        %if (~p1_nogo)
            for x2_ind = 1:Lx
                x2 = x2_ind - 1;
                x2_uc = mod(x2,w_uc);
                x2_nogo = (x2_uc >= forbidden_min) && (x2_uc <= forbidden_max);

                for y2_ind = y1_ind:Lx
                    y2 = y2_ind - 1;
                    y2_uc = mod(y2,w_uc);
                    y2_nogo = (y2_uc >= forbidden_min) && (y2_uc <= forbidden_max);
                    p2_nogo = x2_nogo && y2_nogo;

                    %if (~p2_nogo)
                        % all p1 points from (1,1) to (x1,y1) have already had
                        % all possible combinations of p1,p2 counted
                        % Only count p2s if they start at (x1,y1) or
                        % greater
                        same_row = (y2_ind <= y1_ind);
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
