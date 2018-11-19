clear all
close all

colors = {'k', 'b', 'r', 'g'};

designs = {'des_perf', 'cf_rca_16', 'cf_fft_256_8', 'mul_256_256'};

footprint_mm2 = [0.0655, 0.445, 1.690, 5.198];

Ng_vec = [33024, 146542, 288145, 1639050];

miv_num = [ 0   1800    2738    3823    ;
            0   1747    2925    3936    ;
            0   1050   1921    2475    ;
            0   48513   79382   102994  ];
 
miv_tf_num = [ 0    124     482     1098    ;
               0    324     609     850     ;
               0    105     210     518     ;
               0    13162   20955   24260   ];
twl = [ 563293      478166      432728      415356      ;
        1578160     1499774     1466258     1451201     ;
        4927746     4754600     4745069     4759862     ;
        30224686    26169716    23746536    22470562    ]/1e6;
    
twl_tf = [ 566977   581311      563714      448844      ;
           1578160  1534249     1491711     1473710     ;
           4927746  4848700     4829858     4801278     ;
           30224686 28482288    27610390    27478991    ]/1e6;

       
rent_exps =     [ 0.7  0.6  0.55   0.6    ];
rent_consts =   [ 1.8    1    1    1.25      ];
alpha = 0.5;
k = 1;
tiers = [1,2,3,4];
r=10;

gate_pitch_m = 465e-9*2;
gate_pitch_mm = gate_pitch_m*1e3;
area_per_gate = footprint_mm2./Ng_vec;
gp_per_side = round(sqrt(area_per_gate)/gate_pitch_m);

num_designs = length(designs);

twl_vec = zeros(num_designs, length(tiers));
tsv_vec = zeros(num_designs, length(tiers));

for dind = 1:num_designs
    p = rent_exps(dind);
    k = rent_consts(dind);
    Ng = Ng_vec(dind);
    
    for nind = 1:length(tiers)
        S = tiers(nind);
        footprint_now_mm2 = footprint_mm2(dind)/S;
        gp_per_side = round( sqrt(footprint_now_mm2)/gate_pitch_mm );
        Iidf = xcm.calc_Iidf(alpha,k,p,gp_per_side,S,r);
        twl_vec(dind, nind) = get_total_length(Iidf);
        [nt_max, nt_tot, nt_to, nt_through, Tacmat] = xcm.estimate_tsvs_required(Ng,S,k,p,alpha);
        tsv_vec(dind,nind) = sum(nt_tot);
    end
end

twl_vec_m = twl_vec*gate_pitch_m;


for dind = 1:num_designs
    figure(dind)
    clf
    hold on
    grid on
    plot(twl(dind,:), colors{dind}, 'linestyle', '-', 'linewidth', 1.5)
    plot(twl_tf(dind,:), colors{dind}, 'linestyle', '--', 'linewidth', 1.5)
    plot(twl_vec_m(dind,:), colors{dind}, 'linestyle', ':', 'linewidth', 1.5)
end

for dind = 1:num_designs
    figure(10 + dind)
    clf
    hold on
    grid on
    plot(miv_num(dind,:), colors{dind}, 'linestyle', '-', 'linewidth', 1.5)
    plot(miv_tf_num(dind,:), colors{dind}, 'linestyle', '--', 'linewidth', 1.5)
    plot(tsv_vec(dind,:), colors{dind}, 'linestyle', ':', 'linewidth', 1.5)
end

%WL Normalized
for dind = 1:num_designs
    figure(20 + dind)
    clf
    hold on
    grid on
    plot(twl(dind,:)/twl(dind,1), colors{dind}, 'linestyle', '-', 'linewidth', 1.5)
    plot(twl_tf(dind,:)/twl_tf(dind,1), colors{dind}, 'linestyle', '--', 'linewidth', 1.5)
    plot(twl_vec_m(dind,:)/twl_vec_m(dind,1), colors{dind}, 'linestyle', ':', 'linewidth', 1.5)
end

%WL Normalized bar
for dind = 1:num_designs
    figure(30 + dind)
    clf
    hold on
    grid on
    twl_norm = twl(dind,:)/twl(dind,1);
    twl_calc_norm = twl_vec_m(dind,:)/twl_vec_m(dind,1);
    bar(tiers, [twl_norm; twl_calc_norm]', 1)
    set(gca,'xtick', tiers)
    xlabel('Number of tiers')
    ylabel('Total WL (Normalized)')
end


%WL NON-Normalized bar
for dind = 1:num_designs
    figure(40 + dind)
    clf
    hold on
    grid on
    bar(tiers, [twl(dind,:); twl_vec_m(dind,:)]', 1)
    set(gca,'xtick', tiers)
    xlabel('Number of tiers')
    ylabel('Total WL (m)')
end
