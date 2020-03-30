This folder contains the simulation code used for the results in the manuscript 'Pathological responses to single pulse electrical stimuli in epilepsy: the role of feedforward inhibition'.

**Content**

_read_me.txt		This file
create_NM.m         	Matlab function used to create structs containing all parameters of the neural masses
create_NM_network.m	Matlab function used to crate a struct containing all information of a set of coupled neural masses
init_NM_network.m	Matlab script used to initialize the computation
sim_NM_network.m	Matlab function numerically solves the ODEs of a set of coupled neural masses


**Design note**

The code provides a framework to simulate arbitrary networks of neural masses and not only two feedforward coupled neural masses as used in the manuscript. The idea is that each neural mass is represented by a struct. This struct contains the parameters of that neural mass. This allows to chose the parameters of each neural mass differently. The complete network is represented by another struct. This struct contains a matrix describing how the neural masses are coupled and a vector of neural mass structs.


(c) 2019 Jurgen Hebbink (University of Twente, University Medical Centre Utrecht), Geertjan Huiskamp (University Medical Centre Utrecht), Frans Leijten (University Medical Centre Utrecht), Stephan van Gils (University of Twente) and Hil Meijer (University of Twente)