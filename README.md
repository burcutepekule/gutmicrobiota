# Gut Microbiota Model
Repository for the manuscript _"Quantifying the impact of treatment history on plasmid-mediated resistance evolution in human gut microbiota"_  by Burcu Tepekule, Pia Abel zur Wiesch, Roger Kouyos, Sebastian Bonhoeffer.

# Step 0 : Getting started
### Required Libraries

Two libraries are required for the code to run properly. Please refer to the websites of the libraries to install them to your local environment.

1) GSL - GNU Scientific Library : https://www.gnu.org/software/gsl/
2) Boost : https://www.boost.org/

### Cloning the repository
1) Open the Terminal
2) Go to the directory you wish to clone the repository. As an example, if you want to clone the repository to your ``Desktop`` folder

```sh 
$ cd ./Desktop/
```

3) Clone the repository by typing the following

```sh 
$ git clone https://github.com/burcutepekule/gutmicrobiota.git
```
You will have a folder named ``gutmicrobiota'' on your desktop including all files necessary to run the simulations.

# Step 1 : Running Simulations & Data Generation
### Related Scripts

``GUT_BIOTA_SIM.cpp`` : Main script for the hybrid deterministic-stochastic simulations of the model proposed in the manuscript. Please refer to the comments in the file for more information on how the script works.

``GUT_BIOTA_SIM.sh`` : Bash script to generate the executable and run it.

### Workflow

1) Open the Terminal
2) Go to the directory you have cloned the repository. As an example, if you cloned the repository to your ``Desktop`` folder

```sh 
$ cd ./Desktop/gutmicrobiota/
```
3) Give permission to the bash script by typing the following

```sh
$ chmod +x ./GUT_BIOTA_SIM.sh
```

4) Run the bash script by typing the following 

```sh
$ ./GUT_BIOTA_SIM.sh
```
The bash script will generate a subfolder called ``SIM_RESULTS``, where all the simulation results will be saved.

Simulation results for,
  - Different maximum time for treatment period (denoted by T_L (Figure 7) in the manuscript), 
  - Each number of treatment courses (denoted by N in the manuscript),  

will be saved in a different folder with the name indexed as ``TL_<number of max days for treatment period>_N_<number of treatment courses>``. Since the maximum time for treatment period is T_L=1000 days, and the number of treatment courses N vary from 1 to 20, the folders will be named as ``TL_1000_N_1``, ``TL_1000_N_2``, ..., ``TL_1000_N_3``.

In each folder, you will see four different text files, where each simulation result is kept in a different row (If you run 10 simulations, then you will have 10 rows for each of the files below, with different number of columns depending on the file),

- ``extCounter.txt`` : Two column text file keeping the information about the time step when C_0^{+} and C_1^{+} goes extinct.

- ``samplePops.txt`` : Sampled population abundances  [t,C_0,C_0^{+},C_1,C_1^{+}] at sampling times t for the drug-free period after the treatment period (refer to Figure 7 in the manuscript). Each row represents another simulation with As an example, for the first simulation, if the last treatment is at day 120, and the last treatment lasts for 5 days, and the time step of the simulation is 0.01, you will see a text file with the following columns for the first row,

| day (T_{df}=0)  | C_0 | C_0^{+} | C_1 | C_1^{+} | day (T_{df}=0)  | C_0 | C_0^{+} | C_1 | C_1^{+} | ... | day (T_{df}=360)  | C_0 | C_0^{+} | C_1 | C_1^{+} | 
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
125.01       | ... | ...     | ... | ...     | 140.01        | ... | ...     | ... | ...     | ... | 485.01         | ... | ...     | ... | ...     

- ``schedule_trtInit.txt`` : Initiation day of treatments. Number of columns will depend on the number of treatments, where i^{th} column keeps the initiation day of the i^{th} treatment, denoted by t_{i} in the manuscript (Figure 7).

- ``schedule_trtLen.txt`` : Duration of treatments. Number of columns will depend on the number of treatments, where i^{th} column keeps the duration of the i^{th} treatment, denoted by d_{i} in the manuscript (Figure 7).

# Step 2 : Reading and Saving Data
When all simulations are done, you will have the following folders and text files,

``./SIM_RESULTS/TL_1000_N_<number of treatment courses>/extCounter.txt``

``./SIM_RESULTS/TL_1000_N_<number of treatment courses>/samplePops.txt``

``./SIM_RESULTS/TL_1000_N_<number of treatment courses>/schedule_trtInit.txt``

``./SIM_RESULTS/TL_1000_N_<number of treatment courses>/schedule_trtLen.txt``

These text files need to be read and converted to ``.mat`` files for further processing via MATLAB. 
### Related Scripts

``GUT_BIOTA_SAVEDATA.m`` : \text{MATLAB} script that goes through all the text files generated, and creates two ``.mat`` files, 

- ``allData.mat`` : Raw data including the binary sequence transformation of each simulation (explained in the Supplementary information) is saved.

- ``predictorMat.mat`` : Predictor matrix generated for the predictor importance analysis. Each column represents a different predictor, in the following order

| Number of Treatments (N) | Total days of treatment (\sum d_i) | Duration of last treatment (d_N) | Days to first treatment (t_1-T_I) | Drug-free time after last treatment (T_{df}) | Coefficient of variation (c_v)  | Prevalence of resistance (C_0^{+}/(C_0+C_0^{+})) | 
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |

``CREATE_DATASETS_PI.m`` : Creates dataset tables (.csv files) suitable for building classification and regression forests for the predictor importance analysis. Currently, this script generates subfolders ``PI_DATA_5E<p>``  for each sample size 5x10^{p}. Currently, only sample size 5 x10^{3} is used, but the choice of sample size(s) is customizable. Dataset tables will be named as follows, 

``./R/PI_DATA_5E3/DATA_REG_TRAIN_CLS_BLNCD.csv`` : Classification dataset with sample size of 5x10^3 data points.

``./R/PI_DATA_5E3/DATA_REG_TRAIN_POSONLY.csv`` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; : Regression dataset with sample size of 5x10^3 data points.

These .csv files will be used for the predictor importance analysis in the next step.

# STEP 3 : Predictor Importance Analysis

### Related Scripts

``PI_CLS.r`` : R script running predictor importance analysis using classification forests.

``PI_REG.r`` : R script running predictor importance analysis using regression forests.

``PI.sh `` : Bash script to run the R scripts in loops for different sample and tree sizes.

``PLOT_PREDICTOR_IMP.m`` : MATLAB script to plot the predictor importance results, which is saved in the same directory as an .eps file.

### Workflow

1) Open the Terminal
2) Go to the directory you have cloned the repository. Most likely, this is going to be your ``Downloads`` folder

```sh 
$ cd ./Desktop/gutbiota/
```
3) Give permission to the bash script by typing the following

```sh
$ chmod +x ./PI.sh
```

4) Run the bash script by typing the following 

```sh
$ ./PI.sh
```
The bash script will generate a subfolder called ``PI_RESULTS``, where all the predictor importance analysis results will be saved. Text files will be generated for two different types of predictor importance analysis (classical and conditional), four different number of tree sizes (50,\,100,\,150,\,200), and one sample size (currently only 5x10^3) used in growing the random forests. These files will be named as

``./R/PI_RESULTS/tableNormal_T<number of trees>_CLS_S_5E3.txt `` : Classical classification predictor importance results with sample size of 5x10^3 data points. 

``./R/PI_RESULTS/tableCond_T<number of trees>_CLS_S_5E3.txt`` &nbsp;&nbsp;&nbsp;&nbsp; : Conditional classification predictor importance results with sample size of 5x10^3 data points. 

Results in the manuscript can be reproduced **ONLY** by using the conditional predictor importance results, since they account for the biases for the correlations among the variables. Classical predictor importance results are provided for comparison. **Note that this step can take time due to the computational cost of calculating conditional predictor importances.**

To plot the predictor importance analysis results, open MATLAB and run ``PLOT_PREDICTOR_IMP.m``.

_Note : This whole repository is modified to run locally and serially. Unfortunately, this causes simulations and predictor importance analysis to take considerable amount of computational time. Many parts of this pipeline can be parallelized if the user has access to a cluster (you can contact me for more details on parallelizing the code, easiest way to start with is to run the simulations in parallel for different values of N). In case the user has to run everything locally but has a time constraint, we provided the results used in the manuscript in the_ ``./gutmicrobiota/R/`` _folder of the repository. You can delete this folder if you wish to produce your own results._

**Please feel free to contact me in case of any questions by sending an email to burcu.tepekule@env.ethz.ch**

## SUMMARY

- To clone the repository to your desktop, open the terminal and type,

```sh 
$ cd ./Desktop/
$ git clone https://github.com/burcutepekule/gutmicrobiota.git
```
- To run the simulations and generate data, stay in the terminal and type,

```sh 
$ cd ./gutmicrobiota/
$ chmod +x ./GUT_BIOTA_SIM.sh
$ ./GUT_BIOTA_SIM.sh
```
- After all the simulations are complete, open MATLAB and run ``GUT_BIOTA_SAVEDATA.m``.
- After the .mat files are generated, stay in MATLAB and run ``CREATE_DATASETS_PI.m``.
- After the .csv files are generated, go back to the terminal and type,

```sh
$ chmod +x ./PI.sh
$ ./PI.sh
```
 - To plot the predictor importance analysis results, open MATLAB and run ``PLOT_PREDICTOR_IMP.m``.

## APPENDIX

#### Customizing simulations for recolonization events

Currently, the model used in the manuscript only has **one** event of colonization, and it is prior to all the treatment courses. User can add multiple colonization events at random times during the treatment period, by adjusting ``totalInfs`` and ``infCountVector``, and set the population abundances for colonization events by adjusting ``infAbundVec`` in ``GUT_BIOTA_SIM.cpp``. Please refer to the comments in ``GUT_BIOTA_SIM.cpp`` for more information.
