function Mt2d = Mt_2d_brute_force(Lx)

% Mt2d = zeros(1,2*Lx+1);
% for x1_ind = 1:Lx
%     for y1_ind = 1:Lx
%         for x2_ind = x1_ind:Lx
%             for y2_ind = y1_ind:Lx
%                 l = abs(x1_ind - x2_ind) + abs(y1_ind - y2_ind);
%                 Mt2d(l+1) = Mt2d(l+1) + 1;
%             end
%         end
%     end
% end
%                 

Mt2d = zeros(1,2*Lx+1);
for x1_ind = 1:Lx
    for y1_ind = 1:Lx
        for x2_ind = 1:Lx
            for y2_ind = y1_ind:Lx
                same_row = (y2_ind == y1_ind);
                before_x1 = (x2_ind < x1_ind);
                already_counted = (same_row && before_x1);
                if(~already_counted)
                    l = abs(x1_ind - x2_ind) + abs(y1_ind - y2_ind);
                    Mt2d(l+1) = Mt2d(l+1) + 1;
                end
            end
        end
    end
end
