Source code and data for the paper: "Estimating causal effects in observational studies for survival data with a cure fraction using propensity score adjustment"

Authors: Ziwen Wang, wzwstat@163.com

GUIDELINE FOR REPRODUCING THE RESULTS IN THE PAPER

1.- Table 1 and Table 2 in the paper (simulation_study).
    (a) Obtain "MSE_ACE_SPCEcure_weak_alpxx_PExx_xxxx.RData" in any of the following ways:
        - Running the code in "XXXX_mse.R" in the folder "simulation_study/Methods_MSE" to fit various comparison methods
        - Running load("simulation_study/MSE_intermediate_results/weak/MSE_ACE_SPCEcure_weak_alpxx_PExx_xxxx.RData")
    (b) Represent the results by running the master script "MSE_table1&2.R"
                                             
2.- Table 3 and Table 4 in the paper (simulation_study).
    (a) Obtain "boots_weak_alpxx_PExx_PSBeran.RData" in any of the following ways:
        - Running the code in "PSBeran_boots.R" in the folder "simulation_study/PSBeran_bootstrap"
        - Running load("simulation_study/bootstrap_intermediate_results/boots_weak_alpxx_PExx_PSBeran.RData") 
    (b) Represent the results by running the master scripts "table3.R" and "table4.R"

3.- Table 1 in Supplementary Material (simulation_study).
    Running the master script "Table1_Supplementary.R"

4.- Table 5-6 and Figure 1-2 in the paper, and Table 2 in the Supplementary Material (case_data):
    Running the master script "table5&6_figure1&2_table2supplementary.R"

5.- Figure 3 (case_data): 
    (a) Obtain "bandwidth.RData" and "bandwidth_Beran.RData" in any of the following ways:
        - Running "PSBeran_hn.R" and "Beran_hn.R" in the folder "case_data/methods"
        - Running load("case_data/intermediate_results/bandwidth.RData")
                  load("case_data/intermediate_results/bandwidth_Beran.RData")
    (b) Running the master script "plot_methods_figure3.R"
    
6.- Figure 1 in the Supplementary Material (case_data): 
    (a) Obtain "bandwidth.RData" in any of the following ways:
        - Running "PSBeran_hn.R" in the folder "case_data/methods"
        - Running load("case_data/intermediate_results/bandwidth.RData")
    (b) Running the master script "plot_figure1_supplementary.R"  

7.- Table 7 (case_data)
    (a) Obtain "boots_PHMC_stomach.RData", "boots_PSBeran_stomach.RData",
               "boots_Beran_stomach.RData", "boots_PSKM_stomach.RData" in any of the following ways:
        - Running "XXXX_stomach.R" in the folder "case_data/methods"
        - Running load("case_data/intermediate_results/boots_XXXX_stomach.Rdata")
    (b) Running the master script "table7.R"  


------------------------------------------------------------------------------------------------------

CONTENT

The Code_and_Data folder contains R scripts and data that can be used to reproduce all analyses in the paper and in the Supplementary Material.
An outline of the subfolder structure, which details the name and location of the R scripts and data, is:

/case_data
      table5&6_figure1&2_table2supplementary.R
      table7.R
      plot_methods_figure3.R
      plot_figure1_supplementary.R
      ./data
          stomach.csv	
      ./function
          function_stomach.R        
      ./methods
          Beran_stomach.R
          PHMC_stomach.R
          PSBeran_stomach.R
          PSKM_stomach.R  
          PSBeran_hn.R
          Beran_hn.R
      ./intermediate_results
          boots_Beran_stomach.RData
          boots_PHMC_stomach.RData
          boots_PSBeran_stomach.RData
          boots_PSKM_stomach.RData
          bandwidth.RData
          bandwidth_Beran.RData
       ./results
          Figure1.eps
          Figure2.eps
          Figure3.eps
          Figure1_supplementary.eps
          table5.RData
          table6.RData
          table7.RData
          table2_supplementary.RData
     
/functions
      function_definitions.R

/true
      True_value.R
      True_value.RData
 
/simulation_study
      MSE_table1&2.R
      table3.R
      table4.R
      Table1_Supplementary.R
      ./results
          table1.Rdata
          table2.Rdata
          table3.Rdata
          table4.Rdata
          Table1_supplementary.RData
      ./Methods_MSE
          1Native_mse.R
          2PHMC_mse.R 
          3Beran_mse.R
          4PSKM_mse.R
          5MSPBeran_mse.R
          6PSBeran_mse.R
      ./MSE_intermediate_results
           ./weak
               Native_mse_weak.RData
               PHMC_mse_weak.RData
               PSKM_mse_weak.RData
               MSE_ACE_SPCEcure_weak_alp0_PE20_Beran.RData 
               MSE_ACE_SPCEcure_weak_alp0_PE35_Beran.RData
               MSE_ACE_SPCEcure_weak_alp0_PE50_Beran.RData
               MSE_ACE_SPCEcure_weak_alp05_PE20_Beran.RData 
               MSE_ACE_SPCEcure_weak_alp05_PE35_Beran.RData
               MSE_ACE_SPCEcure_weak_alp05_PE50_Beran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE20_Beran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE35_Beran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE50_Beran.RData
               MSE_ACE_SPCEcure_weak_alp0_PE20_PSBeran.RData 
               MSE_ACE_SPCEcure_weak_alp0_PE35_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp0_PE50_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp05_PE20_PSBeran.RData 
               MSE_ACE_SPCEcure_weak_alp05_PE35_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp05_PE50_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE20_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE35_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE50_PSBeran.RData
               MSE_ACE_SPCEcure_weak_alp0_PE20_MPSBeran.RData 
               MSE_ACE_SPCEcure_weak_alp0_PE35_MPSBeran.RData
               MSE_ACE_SPCEcure_weak_alp0_PE50_MPSBeran.RData
               MSE_ACE_SPCEcure_weak_alp05_PE20_MPSBeran.RData 
               MSE_ACE_SPCEcure_weak_alp05_PE35_MPSBeran.RData
               MSE_ACE_SPCEcure_weak_alp05_PE50_MPSBeran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE20_MPSBeran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE35_MPSBeran.RData
               MSE_ACE_SPCEcure_weak_alp1_PE50_MPSBeran.RData
           ./strong
               Native_mse_strong.RData
               PHMC_mse_strong.RData
               PSKM_mse_strong.RData
               MSE_ACE_SPCEcure_strong_alp0_PE20_Beran.RData 
               MSE_ACE_SPCEcure_strong_alp0_PE35_Beran.RData
               MSE_ACE_SPCEcure_strong_alp0_PE50_Beran.RData
               MSE_ACE_SPCEcure_strong_alp05_PE20_Beran.RData 
               MSE_ACE_SPCEcure_strong_alp05_PE35_Beran.RData
               MSE_ACE_SPCEcure_strong_alp05_PE50_Beran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE20_Beran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE35_Beran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE50_Beran.RData
               MSE_ACE_SPCEcure_strong_alp0_PE20_PSBeran.RData 
               MSE_ACE_SPCEcure_strong_alp0_PE35_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp0_PE50_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp05_PE20_PSBeran.RData 
               MSE_ACE_SPCEcure_strong_alp05_PE35_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp05_PE50_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE20_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE35_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE50_PSBeran.RData
               MSE_ACE_SPCEcure_strong_alp0_PE20_MPSBeran.RData 
               MSE_ACE_SPCEcure_strong_alp0_PE35_MPSBeran.RData
               MSE_ACE_SPCEcure_strong_alp0_PE50_MPSBeran.RData
               MSE_ACE_SPCEcure_strong_alp05_PE20_MPSBeran.RData 
               MSE_ACE_SPCEcure_strong_alp05_PE35_MPSBeran.RData
               MSE_ACE_SPCEcure_strong_alp05_PE50_MPSBeran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE20_MPSBeran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE35_MPSBeran.RData
               MSE_ACE_SPCEcure_strong_alp1_PE50_MPSBeran.RData
     ./PSBeran_bootstrap
          PSBeran_boots.R
     ./bootstrap_intermediate_results
          boots_weak_alp0_PE20_PSBeran_n150.RData
          boots_weak_alp0_PE20_PSBeran_n300.RData
	  boots_weak_alp0_PE20_PSBeran_n500.RData
          boots_weak_alp0_PE35_PSBeran_n150.RData
          boots_weak_alp0_PE35_PSBeran_n300.RData
	  boots_weak_alp0_PE35_PSBeran_n500.RData
          boots_weak_alp0_PE50_PSBeran_n150.RData
          boots_weak_alp0_PE50_PSBeran_n300.RData
	  boots_weak_alp0_PE50_PSBeran_n500.RData
          boots_weak_alp05_PE20_PSBeran_n150.RData
          boots_weak_alp05_PE20_PSBeran_n300.RData
	  boots_weak_alp05_PE20_PSBeran_n500.RData
          boots_weak_alp05_PE35_PSBeran_n150.RData
          boots_weak_alp05_PE35_PSBeran_n300.RData
	  boots_weak_alp05_PE35_PSBeran_n500.RData
          boots_weak_alp05_PE50_PSBeran_n150.RData
          boots_weak_alp05_PE50_PSBeran_n300.RData
	  boots_weak_alp05_PE50_PSBeran_n500.RData
          boots_weak_alp1_PE20_PSBeran_n150.RData
          boots_weak_alp1_PE20_PSBeran_n300.RData
	  boots_weak_alp1_PE20_PSBeran_n500.RData
          boots_weak_alp1_PE35_PSBeran_n150.RData
          boots_weak_alp1_PE35_PSBeran_n300.RData
	  boots_weak_alp1_PE35_PSBeran_n500.RData
          boots_weak_alp1_PE50_PSBeran_n150.RData
          boots_weak_alp1_PE50_PSBeran_n300.RData
	  boots_weak_alp1_PE50_PSBeran_n500.RData
          boots_strong_alp0_PE20_PSBeran_n150.RData
          boots_strong_alp0_PE20_PSBeran_n300.RData
	  boots_strong_alp0_PE20_PSBeran_n500.RData
          boots_strong_alp0_PE35_PSBeran_n150.RData
          boots_strong_alp0_PE35_PSBeran_n300.RData
	  boots_strong_alp0_PE35_PSBeran_n500.RData
          boots_strong_alp0_PE50_PSBeran_n150.RData
          boots_strong_alp0_PE50_PSBeran_n300.RData
	  boots_strong_alp0_PE50_PSBeran_n500.RData
          boots_strong_alp05_PE20_PSBeran_n150.RData
          boots_strong_alp05_PE20_PSBeran_n300.RData
	  boots_strong_alp05_PE20_PSBeran_n500.RData
          boots_strong_alp05_PE35_PSBeran_n150.RData
          boots_strong_alp05_PE35_PSBeran_n300.RData
	  boots_strong_alp05_PE35_PSBeran_n500.RData
          boots_strong_alp05_PE50_PSBeran_n150.RData
          boots_strong_alp05_PE50_PSBeran_n300.RData
	  boots_strong_alp05_PE50_PSBeran_n500.RData
          boots_strong_alp1_PE20_PSBeran_n150.RData
          boots_strong_alp1_PE20_PSBeran_n300.RData
	  boots_strong_alp1_PE20_PSBeran_n500.RData
          boots_strong_alp1_PE35_PSBeran_n150.RData
          boots_strong_alp1_PE35_PSBeran_n300.RData
	  boots_strong_alp1_PE35_PSBeran_n500.RData
          boots_strong_alp1_PE50_PSBeran_n150.RData
          boots_strong_alp1_PE50_PSBeran_n300.RData
	  boots_strong_alp1_PE50_PSBeran_n500.RData

    
Next, an explanation of the R scripts follows.

1.- function_definitions.R

The code of this script is sourced by all .R scripts in the folder "simulation_study"
It contains a collection of functions needed in the computations performed by those scripts. 

2.- XXXX_mse.R

This script contains the code for carrying out the computations needed for obtaining the results shown in Table 1 and Table 2 of the paper 
As indicated in the script, some parameters must be manually set to specify:
(a) PE (gamma0_true = 1.5, 1, 0)
(b) the probability of known cure (alpha1_true = 0, 0.5, 1) and
(c) weak association and strong association will be simulated.
The results obtained under each set of conditions must be saved as a separate workspace in the folder "./simulation_study/MSE_intermediate_results" to be further processed by the master script "MSE_table1&2.R". 
Since the simulations are quite time-consuming, the workspaces are provided already in the folder "./simulation_study/MSE_intermediate_results".

3.- PSBeran_boots.R

This script contains the code for carrying out the computations needed for obtaining the results shown in Table 3 and Table 4 of the paper. 
As indicated in the script, some parameters must be manually set to specify:
(a) the sample size (nn=150, 300, 500)
(b) PE (gamma0_true = 1.5, 1, 0)
(c) the probability of known cure (alpha1_true = 0, 0.5, 1) and
(d) weak association and strong association will be simulated.

As in the previous script, the simulations are quite time-consuming.
To save time, the workspaces are provided already in the folder "./simulation_study/bootstrap_intermediate_results"
Using these intermediate results, tables 3 - 4 in the paper
can be represented using the master scripts in "table3.R" and "table4.R".

4.- function_stomach.R

The code of this script is sourced by all .R scripts in the folder "case_data"
It contains a collection of functions needed in the computations performed by those scripts. 

5.- XXXX_stomach.R

This script contains the code for carrying out the computations needed for obtaining the results shown in Table 7  of the paper 
The results obtained must be saved as a separate workspace in the folder "./case_data/intermediate_results" to be further processed by the master script "table7.R". 

6.- XXXX_hn.R

This is an intermediate process to calculate and save the bandwidth of Beran. 
The results obtained have been saved in the folder "./case_data/intermediate_results" 
to be further processed by the master script "plot_methods_figure3.R" and "plot_figure1_supplementary.R". 


------------------------------------------------------------------------------------------------------

Output of sessionInfo()

> sessionInfo()
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ResourceSelection_0.3-5 forcats_0.5.2           stringr_1.4.1           dplyr_1.0.9            
 [5] purrr_0.3.4             readr_2.1.3             tidyr_1.2.0             tibble_3.1.8           
 [9] ggplot2_3.3.6           tidyverse_1.3.2         npcure_0.1-5            condSURV_2.0.2         
[13] doParallel_1.0.17       iterators_1.0.14        foreach_1.5.2           foreign_0.8-82         
[17] smcure_2.1              survival_3.4-0         

loaded via a namespace (and not attached):
 [1] httr_1.4.4          jsonlite_1.8.2      splines_4.2.1       modelr_0.1.9       
 [5] assertthat_0.2.1    np_0.60-14          doRNG_1.8.2         googlesheets4_1.0.1
 [9] cellranger_1.1.0    pillar_1.8.1        backports_1.4.1     lattice_0.20-45    
[13] quantreg_5.94       glue_1.6.2          quadprog_1.5-8      digest_0.6.29      
[17] rvest_1.0.3         colorspace_2.0-3    Matrix_1.5-1        pkgconfig_2.0.3    
[21] broom_1.0.1         SparseM_1.81        haven_2.5.1         scales_1.2.1       
[25] tzdb_0.3.0          cubature_2.0.4.5    MatrixModels_0.5-1  googledrive_2.0.0  
[29] generics_0.1.3      ellipsis_0.3.2      withr_2.5.0         cli_3.4.1          
[33] crayon_1.5.2        magrittr_2.0.3      readxl_1.4.1        fs_1.5.2           
[37] fansi_1.0.3         MASS_7.3-57         xml2_1.3.3          tools_4.2.1        
[41] hms_1.1.2           gargle_1.2.1        lifecycle_1.0.3     munsell_0.5.0      
[45] reprex_2.0.2        rngtools_1.5.2      compiler_4.2.1      rlang_1.0.6        
[49] grid_4.2.1          rstudioapi_0.14     boot_1.3-28         gtable_0.3.1       
[53] codetools_0.2-18    DBI_1.1.3           R6_2.5.1            zoo_1.8-11         
[57] lubridate_1.8.0     utf8_1.2.2          KernSmooth_2.23-20  permute_0.9-7      
[61] stringi_1.7.8       Rcpp_1.0.9          vctrs_0.4.1         dbplyr_2.2.1       
[65] tidyselect_1.2.0 
