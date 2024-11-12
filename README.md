Two-Step Kernel Method contains the following folder and contents.
1. Functions folder: This contains the functions used in power analysis simulations. The functions include
   I) CODAK_functions.R - Has dcor.comp function used in computing KDC Statistics and P-values for the Zero-Imputation method.
   ii) final_addeffects_fn_version6.R - Has add.new.effects function used to add effects sizes while maintaining the sum to 1 constraint.
   iii)  match.indices.R - This is used by final_addeffects_fn_version6 to identify proportions where effect sizes will be added to.
   iv) new.effects.R - This is used by final_addeffects_fn_version6 to generate new effects for any given vector of proportions.
   v) simulation_functions.R - The function contains the simulation code for the two methods being compared among other functions used in the power analysis simulation.
   vi) two_step_KDC_v2.R - It contains dcor.comp.2step used in computing KDC Statistics and P-values for Two-Step Kernel method.
3. Power Analysis folder: This contains three subfolders;
   i) Codes - Has power simulation code and report generation code. In the report file, you will need to change file names when generating a report for different data.
   ii) Reports - Contains pdf reports from the power analysis report. These are the results to be reported for discussions.
   iii) Simulation Data - This contains data from the power analysis simulation. The data can be used to generate reports in the Reports folder using the report code in the Codes folder.
4. Similarity Plots folder: Contains reports with similarity plots used for preliminary check before running the power analysis simulation.
