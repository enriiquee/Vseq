#
#This program allows to run all the programs, avoiding execution by parts.
#

#Starting the program... 

print("Ejecutando")

#starttime <- Sys.time()

##########  Biostring instalation   @#########################
# source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("Biostrings")


#Installation of Packages that we need. 

list.of.packages <- c("Maeswrap", "data.table", "ggplot2", "plotrix", 
                      "reshape", "MASS", "gridExtra", "grid", "plotrix", 
                      "rJava", "xlsx", "sqldf", "tcltk2", "dplyr", "oce", 
                      "qpcR", "gdata", "tcR", "stats", "XLConnect", "readr", "readxl", "tidyr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

print("Loading Packages...")
#Load Packages
library(Maeswrap);library(data.table); library(ggplot2); library(plotrix); library(reshape); library(MASS); library(gridExtra); library(grid); library(plotrix); library(rJava);
library(xlsx); library(sqldf); library(tcltk2); library(dplyr); library(oce); library(qpcR); library(gdata); library(tcR); library(readxl);
library(stats); library(XLConnect); library(Biostrings);library(readr);library(readxl); library(tidyr)





repeat{
  
  #Creamos la carpeta Datos, donde almacenamos los archivos que son necesarios. 
  output_dir <- file.path(getwd(), "Datos2")
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("La carpeta Datos ya existe")
  }
  
  
  repeat{
    cat("¿Que tipo de datos tienes, DiS o DdS? Por defecto se suele seleccionar DiS.\n ")
    cat("1. DiS | 2. DdS\n")
    input_answer <- readLines(file("stdin"),1)#enter "yes"
    
    if (input_answer=="DiS" | input_answer==1 )  {
      readline(prompt="Mete en la carpeta Datos el archivo xlsx exportado de cualquier motor de busqueda y el mgf. 
               \n Press [enter] to continue")
      
      data_type <- "DiS" #CAMBIAR SI ES UN DiS o DdS
      experimento <- "SDR" 
      
      #Load files 
      infile <- list.files(file.path(getwd(), "Datos2"), pattern = "\\.mgf$") 
      infile2 <- list.files(file.path(getwd(), "Datos2"), pattern = "\\.xlsx$") 
      
      #This test if the file is in the folder. If not BREAK. 
      if (length(infile)<1){
        cat("El archivo mgf no está en la carpeta Datos")
        break
      }
      if (length(infile2)<1){
        cat("El archivo Excel no existe en la carpeta Datos")
        break
      }
      
      
      #Load data frame. 
      sql <-  read_excel(file.path(getwd(), "Datos2", infile2[1]))
      
      fr_ns <- read_delim(file.path(getwd(), "Datos2", infile[1]), 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)
      
      #################        Get the query file   ##################################
      #################################################################################
      ###################################################################################
      #Extract SCANS 
      squery <- filter(fr_ns, grepl("SCANS",  X1))
      #Remove the word SCAN=
      squery2 <- data.frame(SCAN=gsub('SCANS=', '', squery$X1)) 
      
      
      #Take PEPMASS
      mquery <- filter(fr_ns, grepl("PEPMASS", X1))
      #Remove name PEPMASS=
      mquery <- data.frame(MASS=gsub('PEPMASS=', '', mquery$X1)) #Eliminamos PEPMASS
      #Splited in two, having MZ and Intensity. 
      mquery <- separate(data = mquery, col = MASS, into = c("MASS", "INT"), sep = "\\s")
      
      

      
      #### Testing
      
      
      #hsquery <- as.data.frame(substr(squery[,'X1'], 7, 20))
      cquery <- filter(fr_ns, grepl("CHARGE", X1))
      
      
      
      
      
      tquery <- cbind(hsquery, mquery, cquery)
      
      names(tquery)[1] <- "SCAN"
      names(tquery)[2] <- "MASS"
      
      
      
      
      
      
      
      
      
      
      break
    }
    
    else{
      data_type <- "DdS" #CAMBIAR SI ES UN DiS o DdS
      experimento <- "SDR" 
      break
    }
      
    
    
    
    
  }
}