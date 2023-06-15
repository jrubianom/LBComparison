# LBComparison
# Comparison of two Lattice Boltzmann models for Electrodynamics

This repository contains code and resources for comparing the MM (Mendoza-Muñoz) model and HV (Hauser-Verhey) model for electromagnetic simulation. The goal of this project is to evaluate the performance and accuracy of these two models in simulating electromagnetic phenomena using a Lattice Boltzmann method.

## Overview

Electromagnetic simulation plays a crucial role in various fields such as physics, engineering, and telecommunications. The MM model and HV model are two approaches proposed to simulate electromagnetic wave propagation and interactions in linear media using a Lattice Boltzmann method subjected to a BGK scheme.

The MM model published by M. Mendoza and D. Muñoz in [2008](https://arxiv.org/abs/0806.2678) incorporates the relative permittivity and permeability of the medium being simulated, furthermore sources such as currents can be implemented. It provides a reliable framework for accurately representing electromagnetic phenomena in different materials.

On the other hand, the HV model is an alternative approach proposed by A. Hauser and L. Verhey in [2017](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.96.063306). The authors claim that the HV model address certain limitations of the MM model, particularly regarding stability and smooth transitions at interfaces with different permeability and permittivity. Here, we modified the HV model and we added a source term as part of the collision term.

## Repository Structure

This repository contains code implementations of the Modified Hauser-Verhey (HV2017model) and Mendoza-Muñoz (MMmodel) models. The directory structure is organized as follows:

- HV2017model: Contains the code for the Modified Hauser-Verhey model.
- MMmodel: Contains the code for the Mendoza-Muñoz model.

The directories HV2017model and MMmodel have four subdirectories:

- DielectricPulse2Interfaces: Simulates a dielectric Gaussian pulse propagating in one dimension through one interface with different permittivities.
- DielectricPulseSeveralInterfaces: Simulates a dielectric Gaussian pulse propagating in one dimension through several interfaces with different permittivities.
- HertzDipoleRad: Analyzes the radiation pattern and dynamics of a punctual Hertz dipole in three dimensions.
- SkinEffect: Simulates an infinite wave in one dimension penetrating a conducting wall with characteristic skin depth.

Each subdirectory consists of the following components:
- GnuplotFiles or plotfiles folder: Contains Gnuplot scripts for generating plots.
- AuxiliarLibraries folder: Contains header files and auxiliary libraries.
- Makefile: Automates the compilation process.
- main.cpp: Main source code file.

To run the code in each of the four subdirectories, simply execute the make command. This will create folders containing graphs and text files as simulation results.

Additionally, there are two extra directories: MMmodelTimes and HVmodelTimes. These directories include similar cases but also measure the execution times.

## Getting Started

To get started with the repository, follow these steps:

- Clone the repository to your local machine.
- Navigate to the desired model directory (HV2017model, MMmodel, MMmodelTimes or HVmodelTimes.).
- Choose one of the subdirectories based on the simulation you want to run (e.g., DielectricPulse2Interfaces).
- Open the chosen subdirectory and review the available files, including Gnuplot scripts, auxiliary libraries, Makefile, and main.cpp.
- Compile and execute the code by running the make command in the terminal.
- The simulation will produce output files and graphs, which can be found in the created folders.
- Feel free to explore different subdirectories and modify the parameters in the main.cpp file to customize the simulations according to your requirements.
- Feel free to delete the output files by running make clean.

## Contributing

Contributions to this project are welcome! If you have any improvements, bug fixes, or new features to propose, please submit a pull request. Be sure to follow the project's guidelines and coding conventions.

## Acknowledgements

We would like to acknowledge the contributions of the original authors of the MM model and HV model, whose research and insights form the foundation of this project.

## Contact

For any inquiries or further information about this project, please contact [jrubianom@unal.edu.co](mailto:jrubianom@unal.edu.co).

We look forward to collaborating and advancing the understanding and implementation of electromagnetic simulation using the MM model and HV model, and other alterantives Lattice Boltzmann methods.
