#
#This program allows to run all the programs, avoiding execution by parts.
#


#Suppress warnings globally. Avoid error printed
options(warn = -1)
cat("Ejecutando... Puede tardar un poco... \n\n")

#Load packages. 
#Run the rest of the programs. 
setwd(file.path(getwd(),"Programa_R/"))
source('Vseq_Packages.R' )
setwd('..')


#Creamos la carpeta Datos, donde almacenamos los archivos que son necesarios. 
output_dir <- file.path(getwd(), "Datos3")
output_dir2 <- file.path(getwd(), "Vseq_Graphs")
  
if (!dir.exists(output_dir) | !dir.exists(output_dir2)){
  dir.create(output_dir)
  dir.create(output_dir2)
  
} else {
  cat("La carpeta Datos y R_graphs ya existe")
}
  
  
repeat{
  cat("?Que tipo de datos tienes, DiS o DdS? Por defecto se suele seleccionar DiS.\n ")
  cat("1. DiS | 2. DdS\n")
  #input_answer <- readline()#enter "yes"
  input_answer <- readline()#enter "yes"
  
  if (input_answer=="DiS" | input_answer==1 )  {

    #readline(prompt="\n Mete en la carpeta Datos el archivo xlsx exportado de cualquier motor de busqueda y el mgf. 
    # \n Press [enter] to continue")

    repeat{
      cat("\n Mete en la carpeta Datos el archivo .xlsx exportado de cualquier motor de busqueda y el mgf. 
      \n Press [enter] to continue")
      word <- readline() 
      if (word==""){
      break
      }
    }
    #Establecemos las variables.
    data_type <- "DiS" #CAMBIAR SI ES UN DiS o DdS
    experimento <- "SDR" 

    #This test if the file is in the folder. 

    repeat{
      #Load files 
      infile <- list.files(file.path(getwd(), "Datos3"), pattern = "\\.mgf$") 
      infile2 <- list.files(file.path(getwd(), "Datos3"), pattern = "\\.xlsx$") 
      if (length(infile)<1){
        cat("El archivo mgf no est? en la carpeta Datos. Mete el archivo y pulsa [ENTER] \n")
        repeat{
          word <- readline() 
          if (word==""){
            break
          }
        }
      } else if (length(infile2)<1){
      cat("El archivo Excel no existe en la carpeta Datos. Mete el archivo y pulsa [ENTER] \n")
        repeat{
          word <- readline() 
          if (word==""){
            break
          }
        }
      }

      else if (length(infile)>1 | length(infile2)>1){
        cat("Hay m?s de un archivo mgf o xlsx en la carpeta datos. Elimina los archivos extra y pulsa [ENTER] \n")
        repeat{
          word <- readline() 
          if (word==""){
            break
          }
        }
      }

      else if (length(infile)==1 | length(infile2==1)){
        break
      }

      else{
        cat("Hay alg?n problema con los archivos. Revisa que todo est? bien")
      }
    }

    cat("Ejecutando... Puede tardar un poco")

    #Load data frame. 
    sql <-  read_excel(file.path(getwd(), "Datos3", infile2[1]))

    fr_ns <- read_delim(file.path(getwd(), "Datos3", infile[1]), 
            "\t", escape_double = FALSE, col_names = FALSE, 
            trim_ws = TRUE, col_types = cols())

    
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

    #Transformamos el tipo de dato para que Excel lo coja adecuadamente. 
    tquery <- transform(tquery, MZ = as.numeric(MZ), INTENSITY = as.numeric(INTENSITY), SCAN= as.numeric(SCAN))

    #Clean JAVA menory:
    xlcFreeMemory()
    #Export
    #write.xlsx2(tquery, file=file.path(getwd(), "Datos3/tquery.xlsx"), sheetName="sheet1", row.names=FALSE)
    write.xlsx(tquery, file=file.path(getwd(), "Datos3/tquery.xlsx"), sheetName="sheet1", row.names=FALSE)
    #We clean again: 
    xlcFreeMemory()


    #Run the rest of the programs. 
    setwd(file.path(getwd(),"Programa_R/")) #Seleccionamos el Path donde se encuentra el archivo
    source('Vseq_pre.R') #Ejecutamos Vseq_pre
    setwd('..') #Volvemos al anterior directorio


    #Open the tquery table and filter by MS

    cat("\n El archivo Tquery se ha creado en la carpeta Datos. \n")
    
    #input_MZ <- as.numeric(1530)

    #CStdin mete todo como texto. Por eso 
    repeat{
      cat("\n Selecciona la masa Target que quieres buscar. Escribela y pulsa [ENTER]: \n")
      input_MZ <- readline()

      isnumber <- grepl("^[+-]?(\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+)?$",input_MZ)

      if (!isnumber){
        cat("Esto no es un n?mero")
      }else if(input_MZ==""){
        cat("Esto no es un n?mero")
      }else{
        input_MZ<- as.numeric(input_MZ)
        break
      }
    }

    #input_tolerance <- as.numeric(1)
    
    repeat{
      cat("Selecciona la tolerancia (Ventana de MS) que quieres buscar: \n")
      input_tolerance <- as.numeric(readline()) 

      isnumber2 <- grepl("^[+-]?(\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+)?$",input_tolerance)

      if (!isnumber2){
        cat("Esto no es un n?mero")
      }else if(input_MZ==""){
        cat("Esto no es un n?mero")
      }else {
        input_tolerance<- as.numeric(input_tolerance)
        break
      }
    }


  tquery_test <- data.frame(subset(tquery, tquery$MZ>(input_MZ-input_tolerance) & tquery$MZ<(input_MZ+input_tolerance)))
  
  repeat{
    cat("Selecciona la carga que tiene ese peptido: \n")
    charge <- as.numeric(readline()) 
    
    isnumber2 <- grepl("^[+-]?(\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+)?$",charge)
    
    if (!isnumber2){
      cat("Esto no es un n?mero")
    }else if(input_MZ==""){
      cat("Esto no es un n?mero")
    }else {
      input_tolerance<- as.numeric(input_tolerance)
      break
    }
  }
  
  if(charge==1){
    tquery_test <- tquery_test[tquery_test$CHARGE == 'CHARGE=1+',]
  }else if(charge==2){
    tquery_test <- tquery_test[tquery_test$CHARGE == 'CHARGE=2+',]
  }else if(charge==3){
    tquery_test <- tquery_test[tquery_test$CHARGE == 'CHARGE=3+',]
  }else if(charge==4){
    tquery_test <- tquery_test[tquery_test$CHARGE == 'CHARGE=4+',]
  }else if(charge==5){
    tquery_test <- tquery_test[tquery_test$CHARGE == 'CHARGE=5+',]
  }else if(charge==6){
    tquery_test <- tquery_test[tquery_test$CHARGE == 'CHARGE=6+',]
  }
  


  xlcFreeMemory()
  #Export

  write.xlsx(tquery_test, file=file.path(getwd(), "Datos3/tquery_explorer.xlsx"), sheetName="sheet1", row.names=FALSE)
  #We clean again: 
  xlcFreeMemory()
  
  ## Iniciamos el siguiente programa Vseq_explorer.R 
  
  x <- data.frame(X1=tquery_test$SCAN)
  
  #Run the rest of the programs. 
  
  #Path where the result will be saved. 
  varNamePath <- file.path(getwd(), "Vseq_Graphs/")
  
  
  
  setwd(file.path(getwd(),"Programa_R/")) #Seleccionamos el Path donde se encuentra el archivo
  source('Vseq_explorer.R') #Ejecutamos Vseq_pre
  setwd('..') #Volvemos al anterior directorio

  break

}
  else{
    data_type <- "DdS" #CAMBIAR SI ES UN DiS o DdS
    experimento <- "SDR" 
    break
  }





}


cat("WELL DONE!!!")
