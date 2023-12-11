REINFORCEMENT LEARNING FOR SELF-ASSEMBLY OF PATCHY PARTICLES

This file is to show how to use the codes for Reinforcement Learning for self-assemblies of patchy particles, by employing Q-table.
According to [I], the Q-table Q(sigma,T,a), where sigma, T are states and a is action, is trained with N_e epochs, each epoch contains N_step of applying action a (changing T), in each step of N_step, the Brownian Dynamics simulation is performed in a certain time and the last snapshot is used for the next step.  

If use the code, please refer to articles
   (I)  Uyen Tu Lieu, and Natsuhiko Yoshinaga, "Dynamic control of self-assembly of quasicrystalline structures through reinforcement learning", (arXiv:2309.06869)


### 1-STRUCTURE OF THE FILES ### 
   (1) 	RF_patchy.m   		main file (MATLAB)
   (2)  CallSimulation.m   	function to call the simulation, used in (1) (MATLAB)
   (3)  CalculateState.m	function to calculate the sigma ratio of the structure formed by the simulation, used in (1) (MATLAB)
   (4)  PatchyRL_step1.f90	for simulation when the initial configuration is random in position and orientation (1st step of an epoch), used in (2), (Fortran90) 
   (5)  PatchyRL.f90		for simulation when the initial configuration is taken from the previous step, used in (2) (Fortran90)  	


### 2-TESTED ENVIRONMENT ###
We tested the codes in the following environments:
   - Windows 10 ver22H2, Ubuntu 22.04.2 LTS  
   - MATLAB R2018b
   - oneAPI toolkit		 
 	

### 3-HOW TO USE ###
*INPUT:
   - For RF_patchy.m (1), the input for RL is in line #9-19, input for the self-assembly is in line #22-24
   - For CallSimulation.m (2): 
	+ user must choose the lines compatible to either linux (line #17-21, 61-71) or windows (line #23-27, #74-80).  
	+ Line #62,66 and #75,78 can be modified so that the .f90 files in (4) and (5) can be compile by a Fortran compiler. 	
*OUTPUT: 
After each epoch i, the output data is
   - Q-table after each epoch (train_Q_epoch*.dat)       
   - Simulation data of each step in the epoch
   - file saving the input of RL (00input.dat)	
After each epoch i, the output figures are
   - the policy after each epoch, determined by argmax_a Q(sigma, T, a). The blue, grey, red elements on the sigma, T plane corresponds to action decrease T, keep T constant, increase T. (fig_epoch*_policy.fig)
   - the trajectories of sigma and T in each epoch (fig_epoch*_sT.fig)

  	