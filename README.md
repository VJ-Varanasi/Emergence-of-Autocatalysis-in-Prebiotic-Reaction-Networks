# Emergence-of-Autocatalysis-in-Prebiotic-Reaction-Networks
Accompanying code repository for 'Emergence of Autocatalysis in Prebiotic Reaction Networks' by V. Varanasi and J. Korenaga

## Overview

This repository contains all code used to generate the figures and simulations in the manuscript.  
It includes implementations of the RAF detection algorithm, catalytic scheme generators, and theoretical estimates introduced in the paper.


## Repository Structure

| File | Description |
|------|--------------|
| **`figures.ipynb`** | Main Jupyter notebook used to produce all figures. Each cell is labeled with the corresponding figure number from the paper. |
| **`RAF_Emergence_Sim.py`** | Implements the RAF detection algorithm following Hordijk & Steel (2004). Used to estimate the probability of observing an RAF set across network sizes `n` and mean catalysis `f`. |
| **`catalytic_schemes.py`** | Defines the four catalytic schemes analyzed in the paper. Contains both analytical models and random generators for simulation. |
| **`kauffman_network_helper.py`** | Functions for generating Kauffman networks and mapping molecule–reaction relations. |
| **`theoretical_estimate_helper.py`** | Implements the closed-form theory described in Section III, with each function annotated by the corresponding equation number. |
| **`single_solutions.json`** | JSON containing the enumrated single catalyst irrRAFs. Reactions are indexed by number as implemented in kauffman_network_helper.py. |



