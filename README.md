# ReciPSIICOS

Data: https://drive.google.com/drive/folders/1Behww7qRw2fMDG5jxIShGwfZRBUTSg2H?usp=sharing Here you can find forward model sparse and reduced, channels, generated noise, C_re, Upwr and Apwr matrices for whitened ReciPSIICOS reconstruction.

Simulations:

  1. Distances:
	1.1 save_simulations_dist.m — Generate simulations for different methods and different distances for synchronous and asynchronous sources, the resulting activation maps saved at Z_totalDist.mat, and source locations in pickedSrcDist.mat
	1.2 plot_metrics_simulations_dist.mat — calculate point spreading and bias for precomputed Z_totalDist and picked_src and draw them on a graph

  2. Correlations:
	2.1 save_simulations_corr.m — Generate simulations for different methods and different correlations, the resulting activation maps saved at Z_totalCorr.mat, and source locations in pickedSrcCorr.mat
	2.2 plot_metrics_simulations_corr.mat — calculate point spreading and bias for precomputed Z_totalCorr and picked_src and draw them on a graph
