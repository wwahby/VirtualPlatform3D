clear all
close all

Nt = 86e6/4;
S = 4;
Ntsv = 40e3;
r = 10;
pitch_gate_m = 930e-9;

w_tsv = 11;


Ns_nom = Nt/S;

Nuc_1d = round(sqrt(Ntsv));
Ncells = Nuc_1d^2;
w_cell = round(sqrt(Ns_nom/Ncells + w_tsv^2));
Luc_1d = w_cell;

Ns_eff = Ncells*w_cell^2;
Ns_act = Ncells*(w_cell^2-w_tsv^2);

%%

area_ratio = (w_tsv/Luc_1d)^2;

Lx = Luc_1d*Nuc_1d;
Ns = Lx^2;

Lx_orig = round(sqrt(Ns_act));
Ns_eff_orig = Lx_orig^2;



 Mt2dc = xcm.Mt2d_corrected(Lx, Nuc_1d, w_tsv);
 Mt2d = xcm.Mt_2d_joyner(Lx);
 Mt2d_orig = xcm.Mt_2d_joyner(Lx_orig);
 
 err_rel = abs(Mt2d-Mt2dc)./Mt2d;
 
 figure(1)
 clf
 hold on
 plot(err_rel,'k', 'linewidth',1.5)
% set(gca, 'yscale','log')

figure(11)
clf
hold on
plot(Mt2d_orig, 'k', 'linewidth',1.5)
plot(Mt2d, 'b', 'linewidth',1.5)
plot(Mt2dc, 'r', 'linewidth',1.5)
 
 
 %%
 
 Mt3dc = xcm.Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);
 Mt3d = xcm.Mt_3d_joyner(Lx,S,r);
 Mt3d_orig = xcm.Mt_3d_joyner(Lx_orig, S, r);
 
 
 err_rel = abs(Mt3d-Mt3dc)./Mt3d;
 
 figure(2)
 clf
 hold on
 plot(err_rel,'k', 'linewidth',1.5)
% set(gca, 'yscale','log')

figure(21)
clf
hold on
plot(Mt3d_orig, 'k', 'linewidth',1.5)
plot(Mt3d, 'b', 'linewidth',1.5)
plot(Mt3dc, 'r', 'linewidth',1.5)


%%

alpha = 0.5;
k = 3;
p = 0.6;

Iidf = xcm.calc_Iidf(alpha,k,p,Lx,S,r);
Iidfc = xcm.calc_Iidf_corrected(alpha,k,p,Lx,S,r,Nuc_1d,w_tsv);
Iidf_orig = xcm.calc_Iidf(alpha, k, p, Lx_orig, S, r);

err_rel = abs(Iidf-Iidfc)./Iidf;
 
figure(3)
clf
hold on
grid on
plot(err_rel(Iidf>1),'k', 'linewidth',1.5)
%set(gca, 'yscale','log')
xlabel('Interconnect Length (GP)')
ylabel('Relative Error (-)')

figure(31)
clf
hold on
plot(Iidf_orig, 'k', 'linewidth',1.5)
plot(Iidf, 'b', 'linewidth',1.5)
plot(Iidfc, 'r', 'linewidth',1.5)
set(gca, 'yscale', 'log')
set(gca, 'xscale', 'log')


%% figure(4)

figure(4)
clf
hold on
grid on
plot(Iidf(Iidf>0), 'k', 'linewidth', 1.5)
%plot(Iidfc(Iidf>1), 'r', 'linewidth', 1.5)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Interconnect Length (GP)')
ylabel('Number of Interconnects')

tsv_width_vec = 1:2:20;

twl_vec = zeros(1,length(tsv_width_vec) );
twlc_vec = zeros(1, length(tsv_width_vec) );
twl_orig_vec = zeros(1, length(tsv_width_vec) );

for ind = 1:length(tsv_width_vec)
    w_tsv = tsv_width_vec(ind);
    
    
    Ns_nom = Nt/S;
    
    Nuc_1d = round(sqrt(Ntsv));
    Ncells = Nuc_1d^2;
    w_cell = round(sqrt(Ns_nom/Ncells + w_tsv^2));
    Luc_1d = w_cell;
    
    Ns_eff = Ncells*w_cell^2;
    Ns_act = Ncells*(w_cell^2-w_tsv^2);
    
    
    area_ratio = (w_tsv/Luc_1d)^2;
    
    Lx = Luc_1d*Nuc_1d;
    Ns = Lx^2;
    
    Lx_orig = round(sqrt(Ns_act));
    Ns_eff_orig = Lx_orig^2;
    
    
    
    Mt2dc = xcm.Mt2d_corrected(Lx, Nuc_1d, w_tsv);
    Mt2d = xcm.Mt_2d_joyner(Lx);
    Mt2d_orig = xcm.Mt_2d_joyner(Lx_orig);
    
    err_rel = abs(Mt2d-Mt2dc)./Mt2d;
    
    
    Mt3dc = xcm.Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);
    Mt3d = xcm.Mt_3d_joyner(Lx,S,r);
    Mt3d_orig = xcm.Mt_3d_joyner(Lx_orig, S, r);
    
    
    %err_rel = abs(Mt3d-Mt3dc)./Mt3d;
    
    
    alpha = 0.5;
    k = 3;
    p = 0.6;
    
    Iidf = xcm.calc_Iidf(alpha,k,p,Lx,S,r);
    Iidfc = xcm.calc_Iidf_corrected(alpha,k,p,Lx,S,r,Nuc_1d,w_tsv);
    Iidf_orig = xcm.calc_Iidf(alpha, k, p, Lx_orig, S, r);
    
    twl = get_total_length(Iidf);
    twlc = get_total_length(Iidfc);
    twl_orig = get_total_length(Iidf_orig);
    
    twl_vec(ind) = twl;
    twlc_vec(ind) = twlc;
    twl_orig_vec(ind) = twl_orig;
    
end

%%
figure(5)
clf
hold on
grid on
plot(tsv_width_vec*pitch_gate_m*1e6, twl_orig_vec, 'k', 'linewidth', 1.5)
plot(tsv_width_vec*pitch_gate_m*1e6, twl_vec, 'b', 'linewidth', 1.5)
plot(tsv_width_vec*pitch_gate_m*1e6, twlc_vec, 'r', 'linewidth', 1.5)

tsv_width_um = tsv_width_vec*pitch_gate_m*1e6;

err_rel = abs(twl_orig_vec - twlc_vec)./twl_orig_vec;
err_rel_b = abs(twl_vec - twlc_vec)./twl_vec;
figure(6)
clf
hold on
grid on
plot(tsv_width_um, err_rel*100, 'k', 'linewidth', 1.5)
plot(tsv_width_um, err_rel_b*100, 'b', 'linewidth', 1.5)
xlabel('TSV Width (microns)')
ylabel('Relative Error (%)')
ylim([0 100])

