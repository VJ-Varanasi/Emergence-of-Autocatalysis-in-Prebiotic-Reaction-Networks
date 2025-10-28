# Emergence-of-Autocatalysis-in-Prebiotic-Reaction-Networks
Accompanying code repository for 'Emergence of Autocatalysis in Prebiotic Reaction Networks' by V. Varanasi and J. Korenaga


This repository includes all code used to generate figures and simulations in the aforementioned manuscript. A breakdown of the files included can be found below: 

**figures.ipynb:** The code used to produce the data and accompanying figures for the manuscript. All functions can be found within the script or in the corresponding helper files. 
**RAF_Emergence_Sim.py:** Implementation of the RAF algorithm as introduced in 'W. Hordijk and M. Steel, Detecting autocatalytic, self-sustaining sets in chemical reaction systems'. This script is used to simulate the probability of observing an RAF set across multiple network sizes and average catalytic number, f, as shown in Fig 1 of the manuscript.  
**catalytic_schemes.py:** Defines functions associated with the four catalytic schemes studied in the work. Code includes classes for theoretical models and functions used to generate catalytic sets in simulation. 
**kauffman_network_helper.py:** Defines functions used to create Kauffman Networks and run the RAF detection algorithm
**theoretical_estimate_helper.py:** Defines functions used to implement the theory introduced in the paper. Code is labeled with the associated equation numbers of the manuscript. 


