function Ngpyr = calc_Ngpyr(rst,Nst,r)

% derived version. not quite right
Ngpyr = 1/3*r^2*(Nst-1).*(2*Nst-1).*Nst - r*Nst.*(Nst-1).*(2*rst+1) + 2*Nst.*rst.*(rst+1);

% correct version?
%Ng_pyr = 2*r_step.*(1+r_step).*N_step + 1/2*r.^2.*(N_step-1).*(2*N_step-1).*N_step - r*(1+2*r_step).*N_step.*(N_step-1);
%Ngpyr = 1/2*r.^2.*(Nst-1).*(2*Nst-1).*Nst - r*Nst.*(Nst-1).*(2*rst+1) + 2*Nst.*rst.*(rst+1);