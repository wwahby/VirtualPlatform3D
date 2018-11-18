function twl = get_total_length(wld)
    ll = 0:(length(wld)-1);
    twl = sum(ll.*wld);
end