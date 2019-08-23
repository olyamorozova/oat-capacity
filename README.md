# Replication materials

This repository contains replication materials for "Cost-effectiveness of expanding the capacity of opioid agonist treatment in Ukraine: Dynamic modeling analysis" by 
[Olga Morozova](http://campuspress.yale.edu/omorozova/), 
[Forrest W. Crawford](http://www.crawfordlab.io),
[Ted Cohen](https://medicine.yale.edu/lab/tcohen/),
[A. David Paltiel](https://publichealth.yale.edu/people/david_paltiel.profile),
and
[Frederick L. Altice](https://medicine.yale.edu/intmed/people/frederick_altice.profile). 


# Citation 

Morozova O., Crawford F.W., Cohen T., Paltiel A.D., Altice F.L. Cost-effectiveness of expanding the capacity of opioid agonist treatment in Ukraine: Dynamic modeling analysis. _Addiction_ 2019; in press. 

# How to run the code

The code is written in [R statistical computing environment](https://www.R-project.org/) (R Foundation for Statistical Computing; Vienna, Austria 2017) and uses functions from the packages [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) and [deSolve](https://cran.r-project.org/web/packages/deSolve/index.html). 

The functions used to conduct the data analysis presented in the paper are included in the file `oat_capacity_functions.R` To use these functions, download, open and run this file. Comments in the code provide details about the inputs and outputs of each function.

The R dataframe files `smpl_joint_same50K.Rdata`, `smpl_joint_kyiv50K.Rdata`, `smpl_joint_lviv50K.Rdata`, and `smpl_joint_mykolaiv50K.Rdata` contain samples from the joint distribution of model parameters and can be used to replicate the results presented in the paper.

# View the replication document

Download and open the file `oat_capacity_replication_code.pdf` for replication guidelines and examples.
