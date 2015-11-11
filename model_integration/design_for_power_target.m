
function core = design_for_power_target(power_target_W, temperature_target_C, core, simulation)

temperature_C = temperature_target_C;
leakage_reference_temperature_K = core.transistor.leakage_reference_temperature;
Vdd = core.chip.Vdd;
Ioff_per_um = core.transistor.leakage_current_per_micron;
w_trans_min_m = core.transistor.gate_length;
transistor_widths_um = w_trans_min_m *1e6 * 1.5;
transistor_numbers = 4 * core.chip.num_gates;
Vth = core.transistor.Vt;

[Plk_tot, Plk_per_transistor_per_um] = calc_transistor_leakage_at_temp(temperature_C, leakage_reference_temperature_K, Vdd, Vth, Ioff_per_um, w_trans_min_m, transistor_widths_um, transistor_numbers);

P_avail_non_leakage = power_target_W - Plk_tot;


% Need to set the power target
simulation.power_binsearch_target = P_avail_non_leakage;
core = find_power_limited_max_frequency(core, simulation);
