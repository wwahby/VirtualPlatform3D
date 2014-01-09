function ind = binary_search_le(LHSF,RHSF,lbnd,rbnd)
% binary search to find the last point at which LHSF <= RHSF
% LHSF and RHSF are function handles which accept one input (an index
% between lbnd and rbnd)
% LHSF and RHSF are assumed to be monotonic functions, so there is at most
% only one intersection point
% lbnd and rbnd are the smallest and largest indices to consider in this
% binary search

% Have we found an answer yet?
found = 0;
ind = -1;

% First, check to see if lbnd satisfies the condition
% if lbnd doesn't satisfy the condition, we'll never find an answer
if(LHSF(lbnd) > RHSF(lbnd)) % Even lbnd doesn't satisfy this condition
    ind = -1; % Return an invalid index
    found = 1;
end

while (found == 0)
    % First, do we have more than two valid points?
    % If lbnd and rbnd are adjacent, test rbnd
    % At this point lbnd will be known good, so...
    %   If rbnd fails, just return lbnd
    %   If rbnd passes, return rbnd
    if (rbnd - lbnd == 1)
        if (LHSF(rbnd) > RHSF(rbnd))
            ind = lbnd;
        else
            ind = rbnd;
        end
        
        found = 1;
    end
    
    % We'll always round the center point up
    c = ceil(1/2*(lbnd + rbnd));
    
    % check LHSF and RHSF at the center point, and move the search interval
    % up or down accordingly
    if (LHSF(c) <= RHSF(c))
        lbnd = c;
    else
        rbnd = c;
    end
end
    