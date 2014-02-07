3D Virtual Integration Platform (3DVIP)
==========================================
This package is a system prediction tool for 2D and
3DIC systems. It is composed of several different modules,
which work together to give an estimate of the wiring
properties, system power, system temperature, and power
supply noise for a system implemented in 2D or 3D

The main function is codesign_system.m, which can be found in the
model_integration folder. Several test scripts are available
to illustrate the use of the codesign_system wrapper function.

This function is implemented in MATLAB and makes significant use of
vectorization for functional speedup. Matlab 2012a or better is
recommended.

=========
Modules
=========

XCM - Interconnect determination
Determines the wire length distribution (WLD) in the system,
and uses that to predict the number and pitch of the metal
layers required for signal routing, as well as the number
of repeaters required to meet timing requirements.
The wire length distribution used is based on the Joyner/Davis
distribution for 3DICs, but is modified to account for the finite
lateral size of TSVs. Dynamic, leakage, repeater, and wiring power
are all determined in this module.

Thermal - Thermal evaluation
Coarse finite-difference thermal solver which figures out the
system temperature based on the output power and the user-defined
heatsinking parameters

Power Noise - PSN resistive/inductive droop
Power noise simulator for 3DICs. Uses an analytic periodic model
for the power delivery network (PDN), accounting for both
temperature-dependent leakage current as well as the use of
TSVs for power delivery.

