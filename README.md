# ReciPSIICOS

Data source: https://drive.google.com/drive/folders/1Behww7qRw2fMDG5jxIShGwfZRBUTSg2H?usp=sharing Here you can find forward model sparse and reduced, channels, generated noise, "C_re50" (which was renamed from "C_re") and "C_re70" matrices for whitened ReciPSIICOS reconstruction.

Simulations:

	1. Distances:
		1.1 save_simulations_dist.m — Generate simulations for different methods and different distances for synchronous and asynchronous sources, the resulting activation maps saved at Z_totalDist.mat, and source locations at pickedSrcDist.mat. Uses sparse sensor matrix.
		1.2 plot_metrics_simulations_dist.m — Calculate point spreading and bias for precomputed Z_totalDist and picked_src and draw them on a graph
			1.2.1 DistPlot.m — sub function that plots calculated point spreading, bias and detection. Can't be used alone, it is useful to seperate it from plot_metrics_simulations_dist.m for convenient debug.

	2. Correlations:
		2.1 save_simulations_corr.m — Generate simulations for different methods and different correlations, the resulting activation maps saved at Z_totalCorr.mat, and source locations at pickedSrcCorr.mat Uses sparse sensor matrix.
		2.2 plot_metrics_simulations_corr.m — Calculate point spreading and bias for precomputed Z_totalCorr and picked_src and draw them on a graph
			2.2.1 CorrPlot.m — sub function that plots calculated point spreading, bias and detection. Can't be used alone, it is useful to seperate it from plot_metrics_simulations_corr.m for convenient debug.

	3. Grouped Sources correlations:
		3.1 SaveSimulationsGroupedSources_Corr.m — Generate simulations for N=2 sources we will call main and Ngr=2 additional for each main source, which we will call paired. Simulations are generated for different methods and different correlation between main and paired sources. Resulting activation maps are saved at "ZtotalGrpsCorr DD-mm-YY HH-MM-SS.mat", and source locations at "pickedSrcGrpsCorr DD-mm-YY HH-MM-SS.mat". Uses dense sensor matrix.
		3.2 PlotMetricsGroupedSimulations_Corr.m — Calculate point spreading and bias for precomputed ZtotalGrpsCorr and pickedSrcGrpsCorr and draw them on a graph.
			3.2.1 GroupedSourcesCorrPlot.m — sub function that plots calculated point spreading, bias and detection. Can't be used alone, it is useful to seperate it from PlotMetricsGroupedSimulations_Corr.m for convenient debug.