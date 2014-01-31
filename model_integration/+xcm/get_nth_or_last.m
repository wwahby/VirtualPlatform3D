function val = get_nth_or_last(data,n)
% If data is a vector, grab the nth element
% If n is larger than the last element, grab the last element

if (n <= length(data))
    val = data(n);
else
    val = data(end);
end