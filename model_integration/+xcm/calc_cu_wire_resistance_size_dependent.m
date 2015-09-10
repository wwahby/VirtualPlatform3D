function [R, rho_cu, R_cu, R_barrier] = calc_cu_wire_resistance_size_dependent(resistivity_bulk,width,height,barrier_thickness,rho_barrier,wire_length,electron_mfp,specularity_coeff,reflection_coeff)

% Dimensions of entire wire
width_tot = width;
height_tot = height;

% Get dimensions of actual Cu wire
width_barrier = 2*barrier_thickness; % barrier on both sides of wire
height_barrier = barrier_thickness; % barrier only needed on bottom edge of wire
width_cu = width_tot - width_barrier;
height_cu = height_tot - height_barrier;

% Parameters for Cu size effect models
film_thickness = min(width_cu,height_cu);
grain_size = film_thickness;

if ((width_cu > 0) && (height_cu > 0)) % Only find Cu resistivity if there's actually space to put cu in the wire!
    [rho_cu] = xcm.cu_resistivity_fsms(resistivity_bulk,film_thickness,grain_size,electron_mfp,specularity_coeff,reflection_coeff);
    R_cu = rho_cu*wire_length/width_cu/height_cu;
else % otherwise the wire is made entirely of barrier material
    width_barrier = width_tot;
    height_barrier = height_tot;
    R_cu = inf;
    rho_cu = inf;
end

if (barrier_thickness > 0)
    R_barrier = rho_barrier * wire_length/width_barrier/height_barrier;
else
    R_barrier = inf;
end

R = 1/(1/R_cu + 1/R_barrier);