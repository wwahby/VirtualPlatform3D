function fixfigs(figures,linewidth,text_size,axis_size)
% ==================================================
% fixfigs(figures,linewidth,fontsize,axis_size)
% ==================================================
%  Update linewidths and fontsizes in a list of figures
%   figures is a list of figure indices
%   text_size controls the size of all text except the axis numbering
%   axis size controls the axis numbering
%   Good defaults
%       linewidth = 3
%       fontsize = 14
%       axis_size = 12

% if no arguments passed in, fix all figures with default sizes
if (nargin < 1)
    figure_handles = findobj('Type','figure');
end

for i=1:length(figures)
    if nargin < 1
        figh = figure_handles(i); % get handle
        figure(figh); % set current figure so gca works correctly
    else
        figure(figures(i)); % set current figure
        figh = gcf; % get handle
    end
    
    % Check to see if we have a legend.
    % If we have a legend, the first four line objects will be parts of the
    % legend, and we'll have to ignore them when we set curve properties
    legend_check = size(findobj(figh,'Type','axes','Tag','legend'));
    if legend_check(1) == 0
        legend_exists = 0;
    else
        legend_exists = 1;
    end

    hline = findobj(figh,'type','line');
    htext = findall(figh,'type','text');
    set(hline,'Linewidth',linewidth)
    set(htext,'fontSize',text_size)
    set(gca,'fontsize',axis_size)
end
