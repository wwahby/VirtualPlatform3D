% plot power and power density

power_vec = zeros(1,num_stacks);
power_density_vec = zeros(1,num_stacks);
for nind = 1:num_stacks
    power_vec(nind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    power_density_vec(nind) = power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                 
end

f1 = figure(1)
clf
hold on
[ax h1 h2] = plotyy(1:num_stacks, power_vec, 1:num_stacks, power_density_vec/1e4);
set(h1,'color','b')
set(ax(1),'ycolor','k')
set(h2,'color','r')
set(ax(2),'ycolor','k')
ax(1).YLabel.String = 'Power (W)';
ax(2).YLabel.String = 'Power Density (W/cm^2)';
% ax(1).Position(3) = ax(1).Position(3) - 0.01;
ax(2).Position(3) = ax(2).Position(3) - 0.05;
ax(1).Position(2) = ax(1).Position(2) + 0.03;
ax(2).Position(2) = ax(2).Position(2) + 0.03;
xlabel('Number of Tiers')

figh = ax(2);
linewidth = 2;
text_size = 12;
axis_size = 14;
ax(2).FontSize = text_size;
ax(1).YLim = [8 22];
ax(1).YTick = [8:2:22];
% ax(2).YLim = [8 22];
ax(2).YTick = [0:50:200];
hline = findobj(figh,'type','line');
htext = findall(figh,'type','text');
set(hline,'Linewidth',linewidth)
set(htext,'FontSize',text_size)
% set(htext,'Fontsize',axis_size)
fixfigs(1,2,14,12)


figure(2)
clf
% hold on
% dummy axes for the box and background color
ax0 = axes;
set (ax0, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
ax0.Position(3) = ax0.Position(3) - 0.05;
ax0.Position(2) = ax0.Position(2) + 0.03;
ax0.Position(2) = ax0.Position(2) + 0.03;

% first axes for left y-axis
ax1 = axes ('Position', get (ax0, 'Position'));
h1 = plot(ax1, 1:num_stacks, power_vec);
set (ax1, 'Box', 'off', 'Color', 'none', 'YAxisLocation', 'left');
ax1.XLabel.String = 'Number of Tiers';
ax1.YLabel.String = 'Power (W)';
h1.Color = 'b';

% second axes for right y-axis assuming common x-axis controlled by ax1
ax2 = axes ('Position', get (ax0, 'Position'));
h2 = plot(ax2, 1:num_stacks, power_density_vec/1e4);
h2.Color = 'r';
ax2.YLabel.String = 'Power Density (W/cm^2)';
ax2.FontSize = text_size;
figh = ax2;
hline = findobj(figh,'type','line');
htext = findall(figh,'type','text');
set(hline,'Linewidth',linewidth)
set(htext,'FontSize',text_size)
% ax(1).Position(3) = ax(1).Position(3) - 0.01;
set (ax2, 'Box', 'off', 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right');
axes(ax1)
fixfigs(2,2,14,12)