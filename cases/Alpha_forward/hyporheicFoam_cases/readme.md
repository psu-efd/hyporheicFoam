# Readme

This folder contains the 12 cases for the fully coupled simulations using hyporheicFoam.
 - Before run the simulation, the user needs to run blockmesh to generate the mesh files.
 - The user then needs to pre-simulate the flow to steady state and then use the pre-simulated flow field to simulate the solute transport.
 - Each case solves scalar transport in both the surface and subsurface domains. And the dynamic coupling at SWI is fully considered.
