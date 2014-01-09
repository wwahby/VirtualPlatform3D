function [l_av] = get_average_wirelength(iidf)

lengths = 0:length(iidf)-1;

l_av = sum(lengths.*iidf)/sum(iidf);