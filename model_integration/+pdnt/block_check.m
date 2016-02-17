function [id, sign] = block_check(x, y, map)
%check whether (x, y) falls into blocks specified by map
%map format: xl yb w h;
    [len, ~] = size(map);
    id = 0;
    sign = 0;
    for i=1:len
        xl = map(i,1);
        xr = xl+map(i,3);
        yb = map(i,2);
        yt = yb+map(i,4);
        if x >= xl-1e-12 && x <= xr+1e-12 && y >= yb-1e-12 && y <= yt+1e-12
            sign = 1;
            id = i;
            return;
        end
    end
end

