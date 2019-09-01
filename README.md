# gutmicrobiota
Repository for the manuscript "Quantifying the impact of treatment history on plasmid-mediated resistance evolution in human gut microbiota" by Burcu Tepekule, Pia Abel zur Wiesch, Roger Kouyos, Sebastian Bonhoeffer

**STEP 0 : REQUIRED LIBRARIES**

Two libraries are required for the code to run properly. Please refer to the websites of the libraries to install them to your local environment.

1) GSL - GNU Scientific Library : https://www.gnu.org/software/gsl/
2) Boost : https://www.boost.org/

**STEP 1 : RUNNING C++ CODE TO GENERATE DATA**

**FILES**

GUT_BIOTA_SIM.cpp : This is the main script which has the hybrid deterministic-stochastic simulation for the model proposed in the manuscript. Please refer to the comments in the file for more information on how the script works.


GUT_BIOTA_SIM.sh : This is the bash script to generate the executable and run it.

**HOW IT MAKE IT WORK**

1) Open the Terminal
2) Go to the directory you have downloaded the GUT_BIOTA_SIM.cpp and GUT_BIOTA_SIM.sh (cd ./Downloads/gutbiota/)
3) Give permission to the bash script by typing the following : chmod +x ./GUT_BIOTA_SIM.sh
4) Run the bash script by typing the following : ./GUT_BIOTA_SIM.sh

**WHAT IT DOES**

1) The bash script will generate a subfolder called "SIM_RESULTS", where all the simulation results will be saved.

2) Simulation results for different maximum time for treatment period, for each number of treatment courses, and for each number of random recolonizations will be saved in a different folder. Folders are named as follows : "TOTAL_(number of max days for treatment period) _ COUNT _ (number of treatment courses) _ INFS _ (number of recolonizations)"

3) Since no recolonization is considered in the results in the manuscript for now, the maximum time for treatment period is 1000 days, and the number of treatment courses vary from 1 to 20, the folders will be named as "TOTAL_1000_COUNT_1_INFS_0", "TOTAL_1000_COUNT_1_INFS_0", ... , "TOTAL_1000_COUNT_20_INFS_0".

4) In each folder, you will see four different text files indexed with a "_ 0.txt" at the end. 

("_ 0" ending represents the loop index. If the user wants to run the simulations in parallel, the user can submit multiple loops with same number of simulations in it. In this case, the code will generate files with "_ 1", "_ 2", until "_ (number of loops -1)". This can be done via having another for loop in the bash script. For now, the script only submits one loop, and therefore only files with "_ 0" ending is created)  

Each simulation result is kept in a different row (If you run 10 simulations, then you will have 10 rows for each of the files below, with different number of columns depending on the file).T he files are named as follows,

	- c0counter_0.txt : Two column text file keeping the information about the time step when C_0+ and C_1+ goes extinct. 

	- samplePops_0.txt : Sampled population abundances  [t,C_0,C_0^(+),C_1,C_1^(+)] at sampling times (t) for the drug-free time after treatment. As an example, if the last treatment was at day 120, and the last treatment lasts for 5 days, and the time step of the simulation is 0.01, you will see a text file like,

	day (T_df=0) | C_0 | C_0^(+) | C_1 | C_1^(+) | day (T_df=15) | C_0 | C_0^(+) | C_1 | C_1^(+) | ... | day (T_df=360) | C_0 | C_0^(+) | C_1 | C_1^(+)
	125.01       | ... | ...     | ... | ...     | 140.01        | ... | ...     | ... | ...     | ... | 485.01         | ... | ...     | ... | ...     
	
where each row will represent another simulation with different values for day (T_df=0), day (T_df=15), ... , day (T_df=360).

	- schedule_trtInit_0.txt : Initiation day of treatments. (number of columns will depend on the number of treatments)

	- schedule_trtLen_0.txt : Duration of treatments. (number of columns will depend on the number of treatments)


**STEP 2 : RUNNING MATLAB CODE TO READ AND SAVE DATA**

When all simulations are done, you will have the following folders and files,

./SIM_RESULTS/TOTAL_1000_COUNT_1_INFS_0/c0counter_0.txt
                                       /samplePops_0.txt
				       /schedule_trtInit_0.txt
                                       /schedule_trtLen_0.txt

             /TOTAL_1000_COUNT_2_INFS_0/c0counter_0.txt
                                       /samplePops_0.txt
				       /schedule_trtInit_0.txt
                                       /schedule_trtLen_0.txt
			.....

             /TOTAL_1000_COUNT_20_INFS_0/c0counter_0.txt
                                        /samplePops_0.txt
				        /schedule_trtInit_0.txt
                                        /schedule_trtLen_0.txt


These files need to be read and analyzed. This is done via MATLAB, using GUT_BIOTA_SAVEDATA.m This script goes through all the files generated, and creates two .mat files, 

- allData.mat : where all the raw data including the binary sequence transformation of each simulation (explained in the Supplementary information) is also saved.

- predictorMat.mat : where each column is a predictor for the predictor importance analysis, in the following order

	Number of Treatments (N) | Total days of treatment (sum_di) | Duration of last treatment (d_N) | Days to first treatment (t_1-T_I) | Drug-free time after last treatment (T_df) | Coefficient of variation (c_v) | Prevalence of resistance 

predictorMat.mat will can used to generate datasets for predictor importance analysis

**STEP 3 : RUNNING MATLAB CODE FOR PREPARING DATA FOR THE PREDICTOR IMPORTANCE ANALYSIS**

For predictor importance analysis, first the user should create datasets suitable for building classification and regression forests. 

This is done via CREATE_DATASETS_PI.m, which generates 3 different folder with .csv files with different number of samples (this is customizable in the matlab code. These folders and files are named as follows,


./R/PI_DATA_5E2/DATA_REG_TRAIN_CLS_BLNCD.csv (balanced classification dataset with sample size of 5E2 data points)
		     /DATA_REG_TRAIN_POSONLY.csv   (regression dataset with sample size of 5E2 data points)

./R/PI_DATA_5E3/DATA_REG_TRAIN_CLS_BLNCD.csv (balanced classification dataset with sample size of 5E3 data points)
		     /DATA_REG_TRAIN_POSONLY.csv   (regression dataset with sample size of 5E3 data points)

./R/PI_DATA_5E4/DATA_REG_TRAIN_CLS_BLNCD.csv (balanced classification dataset with sample size of 5E4 data points)
		     /DATA_REG_TRAIN_POSONLY.csv   (regression dataset with sample size of 5E4 data points)


These .csv files will be used for the predictor importance analysis in the next step

**STEP 4 : RUNNING R CODE FOR THE PREDICTOR IMPORTANCE ANALYSIS**

Required R code is run from a bash file called PI.sh as follows, 

1) Open the Terminal
2) Go to the directory you have downloaded the PI.sh (cd ./Downloads/gutbiota/)
3) Give permission to the bash script by typing the following : chmod +x ./PI.sh
4) Run the bash script by typing the following : ./PI.sh

The bash script will generate the following text files in a loop for different sample sizes (5E2, 5E3, 5E4) and different number of trees (50 100 150 200) named as follows, 


./R/PI_RESULTS/tableNormal_T<number of trees>_CLS_S_5E<power of sample size>.txt ("classical" classification predictor importance results)
              /tableCond_T<number of trees>_CLS_S_5E<power of sample size>.txt   ("conditional" classification predictor importance results - used in the manuscript to avoid bias due to correlations between variables)
              /tableNormal_T<number of trees>_REG_S_5E<power of sample size>.txt ("classical" regression predictor importance results)
              /tableCond_T<number of trees>_REG_S_5E<power of sample size>.txt   ("conditional" regression predictor importance results - used in the manuscript to avoid bias due to correlations between variables)

*IMPORTANT* -> Results in the manuscript can be reproduced ONLY by using the conditional predictor importance results, since they account for the biases for the correlations among the variables.

**STEP 5 : RUNNING MATLAB CODE TO NORMALIZE AND PLOT PREDICTOR IMPORTANCE RESULTS**

run PLOT_PREDICTOR_IMP.m to generate the predictor importance results, which is saved in the same directory as a figure.


**END**

Please feel free to contact me in case of any questions by sending an email to burcu.tepekule@env.ethz.ch



