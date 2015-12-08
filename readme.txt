3D Virtual Integration Platform (3DVIP)
==========================================
This package is a system prediction tool for 2D and
3DIC systems. It is composed of several different modules,
which work together to give an estimate of the wiring
properties, system power, system temperature, and power
supply noise for a system implemented in 2D or 3D

The main function is codesign_block.m, which can be found in the
model_integration folder. Several test scripts are available
to illustrate the use of the codesign_system wrapper function.
The most useful of these is sweep_sandy_bridge_all_the_things.m,
which includes setup functions to simulate a 32nm Sandy Bridge
CPU core with a range of different material and technological parameters.

This program is implemented in MATLAB and makes significant use of
vectorization for functional speedup. Matlab 2012a or better is
recommended. The program is verified to work with Matlab 2015a or newer.

This program is a constant work in progress, and the version available on
the I3DS website may not always be the most recent. You can find the
latest version at github.com/wwahby/VirtualPlatform3D.

If you have any comments, concerns, questions, or suggestions, please
contact William Wahby at w.wahby@gmail.com


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

