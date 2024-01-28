###############################################################################
#
#  Title:         Estimating causal effects in observational studies for survival 
#                 data with a cure fraction using propensity score adjustment
#
#  Script:        table7.R
#
#  Authors:       Ziwen Wang | Chenguang Wang | Xiaoguang Wang
#
#  Purpose:       Create results of Table 7 in the paper 
# 
###############################################################################

##====================================================================================
#  Remark:  To save time, intermediate results have been saved in the folder
#           "case_data/intermediate_results/...":
#           including the bootstrap bandwidth for covariates 
#           and the results of estimated different methods.
#
#           The saved results can be loaded with load().
#  Example: load("case_data/intermediate_results/boots_PSBeran_stomach")
##====================================================================================

load("case_data/intermediate_results/boots_PHMC_stomach.RData")
Summary_table_PHMC <- Get_table()

load("case_data/intermediate_results/boots_Beran_stomach.RData")
Summary_table_Beran <- Get_table()

load("case_data/intermediate_results/boots_PSKM_stomach.RData")
Summary_table_PSKM <- Get_table()

load("case_data/intermediate_results/boots_PSBeran_stomach.RData")
Summary_table_PSBeran <- Get_table()

race <-
  rbind(
    Summary_table_PHMC$table_result[1, ],
    Summary_table_Beran$table_result[1, ],
    Summary_table_PSKM$table_result[1, ],
    Summary_table_PSBeran$table_result[1, ]
  )
row_name <- c("PHMC", "Beran", "PS_KM", "PS_Beran")
col_name <- c("Est", "SE", "Lower", "Upper", "p_value")
dimnames(race) <- list(row_name, col_name)

spce <-
  rbind(
    Summary_table_PHMC$table_result[2, ],
    Summary_table_Beran$table_result[2, ],
    Summary_table_PSKM$table_result[2, ],
    Summary_table_PSBeran$table_result[2, ]
  )
row_name <- c("PHMC", "Beran", "PS_KM", "PS_Beran")
col_name <- c("Est", "SE", "Lower", "Upper", "p_value")
dimnames(spce) <- list(row_name, col_name)

## ========================================================
## Summary of the results for table 7
## ========================================================
table7 <- list(race = race, spce = spce)
print(table7)


#save(table7, file = "case_data/results/table7.RData")

