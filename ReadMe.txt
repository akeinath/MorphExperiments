All custom code for analysis of the data from:

	Keinath, Mosser, Brandon. (2022). The hippocampal representation of context is 	preserved despite neural drift. 

Also requires the circular statistics toolbox (2012a) to generate some of the figures. 

To run in MATLAB, set this folder to your root directory, add the Functions folder and subfolders to the path, replace 'FOLDER_CONTAINING_DATA' in main.m with the folder containing the data (downloaded as per the specs in the Methods), and run main.m.

Originally ran in MATLAB version R2021a. Typically takes ~5-10min to output all the figures to run, depending on the particular hardware. All figures will be written to the 'Plots' folder. Some statistics output to the Command Window, others saved as text files accompanying the figures in Plots.