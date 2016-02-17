    t = 0:0.1:20;
    x1 = (0.1:0.1:1)*75;
    x2 = ones(191,1)*75;
    x = [x1'; x2]';
    plot(t, x, 'Linewidth',5);
    set(gca,'FontSize',16);
    xlabel('T(ns)');
    ylabel('Die #1 Power(w)');set(gca,'FontSize',20);