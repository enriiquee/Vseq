
#This program allows to run all the programs, avoiding execution by parts.



repeat{

  output_dir <- file.path(getwd(), "Datos2")
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("La carpeta datos ya existe")
  }
  
  
  repeat{
    cat("¿Que tipo de datos tienes, DiS o DdS? Por defecto se suele seleccionar DiS.\n ")
    cat("1. DiS | 2. DdS\n")
    input_answer <- readLines(file("stdin"),1)#enter "yes"
    
    if (input_answer=="DiS" | input_answer==1 )  {
      data_type <- "DiS" #CAMBIAR SI ES UN DiS o DdS
      experimento <- "SDR" 
      
      datapath <- list.files(file.path(getwd(), "Datos"), pattern = "\\.mgf$") 
      
      
      datapath <- list.files(file.path(getwd(), "Datos"), pattern = "\\.xlsx$") 
      
      break
    }
    
    else{
      data_type <- "DdS" #CAMBIAR SI ES UN DiS o DdS
      experimento <- "SDR" 
      break
    }
      
    
    
    
    
  }
}