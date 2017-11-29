#Starting the program... 

options(warn = -1)
cat("Ejecutando... Puede tardar un poco... \n\n")

#starttime <- Sys.time()

##########  Biostring instalation   @#########################
# source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("Biostrings")


#Installation of Packages that we need. 

list.of.packages <- c("Maeswrap", "data.table", "ggplot2", "plotrix", 
                      "reshape", "MASS", "gridExtra", "grid", "plotrix", 
                      "rJava", "xlsx", "sqldf", "tcltk2", "dplyr", "oce", 
                      "qpcR", "gdata", "tcR", "stats", "XLConnect", "readr", "readxl", "tidyr", "rJava", "progress")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

print("Loading Packages...")
#Load Packages
library(Maeswrap);
library(data.table);
library(ggplot2);
library(plotrix);
library(reshape);
library(MASS);
library(gridExtra);
library(grid);
library(plotrix);
library(rJava);
library(xlsx);
library(sqldf);
library(tcltk2);
library(dplyr);
library(oce);
library(qpcR);
library(gdata);
library(tcR);
library(readxl);
library(stats);
library(XLConnect);
library(Biostrings);
library(readr);
library(readxl);
library(tidyr);
library(progress);



#Clean console
shell(cat("\014")  )

