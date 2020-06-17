# Linking Minimal and Detailed Hippocampal CA1 Networks

Melisa Gumus - May 2018 (uploaded June 2020)

Languages: MATLAB

<img src="https://github.com/FKSkinnerLab/CA1_SimpleDetailed/blob/master/SimpleDetailedMovie.gif" height="300px" width="500px" >

## Summary ##
This github repository includes complete scripts for comparing the minimal (Ferguson et al., 2017) and the detailed (Beziare et al., 2016) microcircuit models of hippocampus CA1 in terms of theta rhythm generation. For more details, please see Chatzikalymniou and colleagues (2020).

This repository includes several folders listed below; each contains scripts for a specific purpose and calculations. 

1. **Excitatory_Inhibitory_Ratios**

     Excitatory and inhibitory currents, that each cell receives (specifically, PYR, BC, BiC and AAC) in the detailed model, are obtained using SimTracker and Network Clamp tool (Bezaire et al., 2016). Each simulation was run for 1000 msec. Each script helps you calculate the following (Please note that this example is based on PYR cells onlys , but the same concept applies to other cells/scripts):

     1. Find Mean EPSCs or IPSCs and standard deviation from each cell (e.g. BC, BiC, AAC separately) or combination of cells (e.g. BC and BiC together) onto PYR cells.
     2. Identify the peaks in EPSCs and IPSCs onto PYR, and find the mean and stardard deviation.
     3. Calculate the EPSCs/IPSCs ratios on PYR cells.
     4. Voltage recordings

_**Note:** IPSCs and EPSCs in detailed model obtained from Network Clamp are saved as .dat files and can be found in /Excitatory_Inhibitory_Ratios/Network_Clamp_Results/mytrace_'cell_ID'_syns.dat. The voltage recordings for each cell can be found in /Excitatory_Inhibitory_Ratios/Network_Clamp_Results/mytrace_'cell_ID'_soma.dat_

2. **Anton's stuff (another repository)**



## Status ##
This project is completed, and the relevant publication (Chatzikalymniou et al., 2020) is listed below.


## Download the Repository ##
Download a copy of the repo somewhere locally by running:

`git clone https://github.com/FKSkinnerLab/CA1_SimpleDetailed.git`


## Contact ##
If there are any difficulties or errors in running any of the scripts or something is unclear, please feel to contact me at melisa.gumus@mail.utoronto.ca.

Project link: https://github.com/FKSkinnerLab/CA1_SimpleDetailed. To give credit, please mention this GitHub repository by the name and the link.


## References ##
Bezaire, M.J., Raikov, I., Burk, K., Vyas, D., Soltesz, I., 2016. Interneuronal mechanisms of hippocampal theta oscillations in a full-scale model of the rodent CA1 circuit. eLife 5, e18566. https://doi.org/10.7554/eLife.18566

Chatzikalymniou, A. P.; Gumus, M.; Lunyov, A; Rich, S.; Lefebvre, J.; Skinner, F.K (submitted, 2020). 
Translating mechanisms from minimal to detailed models of CA1 hippocampal microcircuits reveals how theta rhythms emerge and how their frequencies are controlled - bioarxiv 

Ferguson, K.A., Chatzikalymniou, A.P., Skinner, F.K., 2017. Combining Theory, Model, and Experiment to Explain How Intrinsic Theta Rhythms Are Generated in an In Vitro Whole Hippocampus Preparation without Oscillatory Inputs. eNeuro 4. https://doi.org/10.1523/ENEURO.0131-17.2017
