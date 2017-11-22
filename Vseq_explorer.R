

##      ###   ##   ###        ##
##     #    # #  #  #       ##
##    ###  ###  #  #      ##
##     #  #    #  #     ##         ###################################
##  ###  ###   ###   ##           ###################################
##              #  ##            #######  Vseq_EXPLORER  ###########
##            # ##             ###################################
##           ##               ###################################
##       ##
##  ##
##




#write.csv(tquery, "C:\\Users\\proteomica\\Desktop\\R_projects\\PESA\\tquery.csv")

####################################################################################3
starttime <- Sys.time()
####################################################################################

##      ###   ##   ###        ##
##     #    # #  #  #       ##
##    ###  ###  #  #      ##
##     #  #    #  #     ##
##  ###  ###   ###   ##
##              #  ##
##            # ##
##           ##
##       ##
##  ##
##

seq <- "LSETVAICR"
# seq <- "GLPDQMLYR"
# seq <- "LTSSVTAYDYSGK"
# seq <- "ADGVPIHLK"
# seq <- "EDPNLVPSISNK"
# seq <- "SATYVNTEGR"

#####   calculo de la masa de la seq    ###############
len<-nchar(seq)

yy <- substr(seq, 1, len)
print(yy)
yn <- strsplit(yy, "")
z <- as.data.frame(yn)
z<-as.matrix(z)

pas<- paste0("c(", as.character(z), ")")

tvar<-unlist(lapply(pas, 
                    function(x) sum(eval(parse(text=x)))))


parental <-sum(tvar)+19.0178
parental2c <- (parental + 1.0078)/2
parental3c <- (parental + 2.0156)/3
parental4c <- (parental + 3.0234)/4
parental5c <- (parental + 4.0312)/5
parental6c <- (parental + 5.0390)/6
parental7c <- (parental + 6.0468)/7

print(parental)
print(parental2c)
print(parental3c)
print(parental4c)
print(parental5c)
print(parental6c)
print(parental7c)



###########################################################################
############################################################################
##########################################################################
##################    ####################################

x = t(read.table("E:\\Vseq\\SDR_eca_explorer3.csv"))


###################################################################################
snvar = x

varconc3<- paste0(varNamePath,"SCAN=",snvar[1],"_",
                  "SCAN=",snvar[length(snvar)],"_",seq,experimento,".xls")


pb <- winProgressBar(title=paste("Vseq_EXPLORER", "-", seq), 
                     label="0% done", min=0, max=100, initial=0)

count = 1

for (ww in x)
{
  #################################################################################
  ##############################################################################
  ################################################################################
  ##################################################################################
  sn <- paste0("SCANS=", ww)
  source("C:\\Users\\Administrador\\Desktop\\R_projects\\Vseq_explorer_project2c.R")
  #params <- data.frame(nrow=length(x))
  params[[ww]] <-  data.frame(seq,mimMZ,ww,Escore,DeltaMass2c, DeltaMass3c, DeltaMass4c,
                              DeltaMass5c, DeltaMass6c, matched_ions, out_of)
  
  
  #Sys.sleep(0.1) # slow down the code for illustration purposes
  info <- sprintf("%d%% done", round((100/length(x)*count)))
  
  setWinProgressBar(pb, round((100/length(x)*count)), label=info)
  
  count=count +1
  
}  



uu = 1
dtapar <- data.frame()

while(uu <=length(x)){
  dtapar <- rbind(dtapar, params[[x[uu]]])
  uu = uu+1
}

write.xlsx(dtapar, varconc3, sheetName = "params")

close(pb)

endtime <- Sys.time()

print(paste(length(x),"-scans-from-",starttime,"-","to", "-", endtime))
