# 1. Introduction
Simulation of a Cortico-Basal Ganglia-Thalamocortical network model.
The architecture and dynamics of the network can be seen in the reference.
The objective is to analyze the effect of the intensities of the synaptic connections on the dynamics of the network.
In particular, associating pathological and physiological states in Parkinson's disease.
Furthermore, it is desired to study the effect of deep brain stimulation on the pathological dynamics of the network.

# 2. References
https://doi.org/10.1371/journal.pone.0182884

# 3. About code
How to compile: gcc -o main BGThNetwork_1DModel.c
How to execute: 
./main
A_DBS Period_DBS Width_Pulse Slope

where A_DBS: amplitude of stimulation (e.g 10)
	  Period_DBS: Period of stimulation (ms) (e.g 7) 
	  Width_Pulse: Width of pulse train (ms) (e.g 0.5)
	  Slope: Slope (e.g 0)

# 4. About files
- 4.1: Ganglios.dat
		Numeric matrix (number of populations x 4). 
		For each population, threshold, constant external input, noise level and number of neurons are stored.

- 4.2: Interaction.dat
		Numeric matrix (number of "conductores" x 4)
		For each "conductor", pre-synaptic, post-synaptic, average of number of connections, level of divergence of connections are stored. 

- 4.3: Conductores.dat
		Numeric matrix (number of "conductores" x 3)
		For each "conductor", synaptic efficac, time constant and delay are stored.

- 4.4: ElectrodosLFP.dat
		Numeric matrix (number of populations x 4)
		For each population, The parameters of the spatial modulation of the measurement (electrode position, threshold, volume of tissue involved) are saved.

- 4.5: ElectrodosEstim.dat
		Numeric array (4)
		The parameters of spatial modulation of stimulation are saved (target population, position, threshold, volumen of activated tissue)