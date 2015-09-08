function [C_top, C_side, Ctot_lumped, Cc_lumped, Ctot_lumped_no_coupling] = gnr_get_capacitance_ns2015(epsr_dielectric, h1, h2, num_layers, Ef, gnr_width, temp)
    % Get quantum capacitance for a graphene nanoribbon
    % Uses the method described in Nishad and Sharma 2015
    %.
    % == Inputs ==
    % epsr_dielectric (-) : Relative permittivity of surrounding dielectric
    % h1 (m) : Distance from top of MLGNR to bottom of next conducting level above
    % h2 (m) : Distance from bottom of MLGNR to top of next lowest conducting level
    % Ef (eV) : Fermi level of doped MLGNRs
    % Width (m) : Width of MLGNR
    % Temp (K) : Ambient temperature
    %.
    % == Outputs ==
    % C_top (F/m) : PUL Capacitance of top-contacted MLGNR interconnect
    % C_side (F/m) : PUL Capacitance of side-contacted MLGNR interconnect


    %%  Constants
    q = 1.602e-19; %% electron charge in coulomb
    h=6.62e-34; %%% Plancks constant in J-s
    vf=8e5; %% Fermi velocity in m/s
    eps0 = 8.854e-12; % (F/m) Vacuum permittivity

    epsd = epsr_dielectric;% * eps0;
    w = gnr_width;
    delta = 0.35e-9; % (m) Interlayer spacing in MLGNRs

    %% Quantum capacitance element
    Nch = xcm.gnr_get_num_channels( Ef, gnr_width, temp );
    Cq_el = 4*q^2*Nch / (h*vf);

    %% Capacitance elements
    % Upper and lower capacitances
    C_10 = epsd*eps0 * calc_M( tanh(pi*w/4/h1) );
    C_N0 = epsd*eps0 * calc_M( tanh(pi*w/4/h2) );

    % Internal capacitances
    k_diel = 1; % (-) relative permittivity of medium between the graphene sheets (assuming vacuum per NS2015)
    C_int = k_diel*eps0 * (w/delta + 4*log(2)/pi);

    %% Capacitance matrices

    cq_diag = ones(1,num_layers) * Cq_el;

    ce_diag = ones(1,num_layers) * 2*C_int;
    ce_diag(1) = C_10 + C_int;
    ce_diag(end) = C_N0 + C_int;

    ce_diag_p1 = ones(1,num_layers - 1) * (-C_int);

    Cq = diag(cq_diag, 0);
    Ce = diag(ce_diag, 0) + diag(ce_diag_p1, -1) + diag(ce_diag_p1, +1);

    %% Overall capacitance
    C = (Cq^-1 + Ce^-1)^-1;

    C_top = C(1,1);
    C_side = sum(sum(C));
    
    [Ctot_lumped, Cc_lumped, Ctot_lumped_no_coupling] = calc_lumped_capacitance(epsr_dielectric, Cq_el, C_N0, gnr_width, gnr_width);


end



function M = calc_M(gamma)
    A = (1+ (1-gamma^2)^(1/4) )/ (1 - (1-gamma^2)^(1/4) );

    B = (1 + sqrt(gamma)) / (1-sqrt(gamma));

    if ( (gamma < 1/sqrt(2) ) && (gamma > 0) )
        M = 2*pi/log(2*A);
    elseif ( (gamma >= 1/sqrt(2) ) && (gamma <= 1) )
        M = 2/pi * log(2*B);
    else
        M = 0;
    end
end


function [Ctot, Cc, Ctot_no_coupling] = calc_lumped_capacitance(epsr_dielectric, Cq, Ce, gnr_width, gnr_space)
     % not based on NS2015 paper. Simple lumped element formulae to include
     % impact of coupling capacitances to the side.
    eps0 = 8.854e-12; % (F/m) Vacuum Permittivity
    
    Cc = epsr_dielectric*eps0/4 * calc_M( sqrt(1 - (1+2*gnr_width/gnr_space)^-2) );
    CA = (1/(Cq + Ce) + 1/Cc )^(-1);
    
    Ctot = (1/Cq + 1/(2*CA + Ce))^-1;
    
    CA2 = (1/(Cq + Ce) )^(-1);
    Ctot_no_coupling = (1/Cq + 1/(2*CA2 + Ce))^-1;
end

