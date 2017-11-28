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
                      "qpcR", "gdata", "tcR", "stats", "XLConnect", "readr", "readxl", "tidyr", "rJava")

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
library(tidyr)


#Clean console
shell(cat("\014")  )



  
  #Creamos la carpeta Datos, donde almacenamos los archivos que son necesarios. 
  output_dir <- file.path(getwd(), "Datos3")
  output_dir2 <- file.path(getwd(), "Vseq_Graphs")
  
  if (!dir.exists(output_dir) & !dir.exists(output_dir2)){
    dir.create(output_dir)
    dir.create(output_dir2)
    
  } else {
    print("La carpeta Datos y R_graphs ya existe")
  }
  
  
  repeat{
    cat("¿Que tipo de datos tienes, DiS o DdS? Por defecto se suele seleccionar DiS.\n ")
    cat("1. DiS | 2. DdS\n")
    input_answer <- readLines(file("stdin"),1)#enter "yes"
    
    if (input_answer=="DiS" | input_answer==1 )  {
      
      #readline(prompt="\n Mete en la carpeta Datos el archivo xlsx exportado de cualquier motor de busqueda y el mgf. 
              # \n Press [enter] to continue")
      
      repeat{
      cat("\n Mete en la carpeta Datos el archivo .xlsx exportado de cualquier motor de busqueda y el mgf. 
               \n Press [enter] to continue")
      word <- readLines(file("stdin"),1) 
      print(word)
      if (word==""){
        break
      }
      
      }
      data_type <- "DiS" #CAMBIAR SI ES UN DiS o DdS
      experimento <- "SDR" 
      
      #Load files 
      infile <- list.files(file.path(getwd(), "Datos3"), pattern = "\\.mgf$") 
      infile2 <- list.files(file.path(getwd(), "Datos3"), pattern = "\\.xlsx$") 
      
      #This test if the file is in the folder. If not BREAK. 
      repeat{ 
      if (length(infile)<1){
        cat("El archivo mgf no está en la carpeta Datos \n")
        
      }
      if (length(infile2)<1){
        cat("El archivo Excel no existe en la carpeta Datos \n")
        
      }
      
      if (length(infile)>1 | length(infile2)>1){
        cat("Hay más de un archivo mgf o xlsx en la carpeta datos \n")
        
      }
      if (length(infile)==1 | length(infile2==1)){
        break
      } 
      }
      
      #Load data frame. 
      sql <-  read_excel(file.path(getwd(), "Datos3", infile2[1]))
      
      fr_ns <- read_delim(file.path(getwd(), "Datos3", infile[1]), 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)
      
      #################        Get the query file   ##################################
      #################################################################################
      ###################################################################################
      #Extract SCANS 
      squery <- filter(fr_ns, grepl("SCANS",  X1))
      #Remove the word SCAN=
      squery <- data.frame(SCAN=gsub('SCANS=', '', squery$X1)) 
      
      
      #Take PEPMASS
      mquery <- filter(fr_ns, grepl("PEPMASS", X1))
      #Remove name PEPMASS=
      mquery <- data.frame(MASS=gsub('PEPMASS=', '', mquery$X1)) #Eliminamos PEPMASS
      #Splited in two, having MZ and Intensity. 
      mquery <- separate(data = mquery, col = MASS, into = c("MASS", "INT"), sep = "\\s")
      
            #Take Change: 
      cquery <- filter(fr_ns, grepl("CHARGE", X1))
      
      #Finally we merge last columns
      tquery <- cbind(squery, mquery, cquery)
      colnames(tquery) <- c("SCAN", "MZ", "INTENSITY", "CHARGE")
      
      
      tquery <- transform(tquery, MZ = as.numeric(MZ), 
                     INTENSITY = as.numeric(INTENSITY), SCAN= as.numeric(SCAN))

      #Clean JAVA menory:
      xlcFreeMemory()
      #Export
      write.xlsx2(tquery, file=file.path(getwd(), "Datos3/tquery.xlsx"), sheetName="sheet1", row.names=FALSE)
      #We clean again: 
      xlcFreeMemory()
      
      
      #Run the rest of the programs. 
      setwd(file.path(getwd(),"Programa_R/"))
      source('Vseq_pre.R')
      
      #Open the tquery table and filter by MS
      break
      
      
      readline(prompt="El archivo tquery se ha creado. Abrelo y selecciona la masa que quieres buscar. \n 
               PULSE ENTER TO CONTINUE [ENTER]")  
      
      
      
      cat("Introduzca el número de MZ seleccionado: ")
      input_MZ <- readLines(file("stdin"),1)#enter "yes"
      
      cat("Introduzca la tolerancia de +- MZ: ")
      input_tolerance <- readLines(file("stdin"),1)#enter "yes"
      
      tquery2 <- data.frame(subset(tquery, tquery$MZ>(input_MZ-input_input_tolerance) & tquery$MZ<(input_MZ+input_tolerance)))
      
      
      break
    }
    
    else{
      data_type <- "DdS" #CAMBIAR SI ES UN DiS o DdS
      experimento <- "SDR" 
      break
    }
      
    
    
    
    
  }
  
  
  cat("WELL DONE!!!")
