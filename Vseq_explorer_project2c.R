

#xxx =  4
#ww = xxx

seq2 <- strsplit(seq, NULL)[[1]]
seq2 <- paste(rev(seq2), collapse='')
#############################################            

fr <- fread(infile,skip=sn); as.data.frame(fr)
snc <- paste0(ww ,",")

fr<-as.data.frame(fr)

len<-nchar(seq)
out_of = len*2

yy <- substr(seq, 1, len)
print(yy)
yn <- strsplit(yy, "")
z <- as.data.frame(yn)
z<-as.matrix(z)

pas<- paste0("c(", as.character(z), ")")

tpar<-unlist(lapply(pas, 
                 function(x) sum(eval(parse(text=x)))))




mim <- subset(tquery, SCAN==snc)
mimMZ <- mim[,2]

DeltaMass2c = mimMZ - parental2c
DeltaMass3c = mimMZ - parental3c
DeltaMass4c = mimMZ - parental4c
DeltaMass5c = mimMZ - parental5c
DeltaMass6c = mimMZ - parental6c
DeltaMass7c = mimMZ - parental7c

# print(DeltaMass2c)
# print(DeltaMass3c)
# print(DeltaMass4c)
# print(DeltaMass5c)
# print(DeltaMass6c)
# print(DeltaMass7c)

###    preparando matriz de fragmentos teóricos    ###


######      calculate y series     #############


ipar=c(2:len-1)
outy <- matrix(NA, nrow=len, ncol= nrow(fr)-1)

for (ipar in 1:len)
{
  yy <- substr(seq, ipar, len)
  print(yy)
  yn <- strsplit(yy, "")
  z <- as.data.frame(yn)
  z<-as.matrix(z)
  tpar<-unlist(lapply(paste0("c(", as.character(z), ")"), 
                   function(x) sum(eval(parse(text=x)))))
  fragy <-sum(tpar)+19.0178
  
  #print (parent)
  #yions<-matrix(parent, ncol = len)
  outy[ipar,] <- fragy
}
###################################################################
###################################################################

###################    calculate b series             ######### ###


j=c(2:len-1)

outb <- matrix(NA, nrow=len, ncol= nrow(fr)-1)

for (j in 1:len)
{
  
  bb <- substr(seq2, j, len)
  print(bb)
  bn <- strsplit(bb, "")
  z <- as.data.frame(bn)
  z<-as.matrix(z)
  tpar<-unlist(lapply(paste0("c(", as.character(z), ")"), 
                   function(x) sum(eval(parse(text=x)))))
  fragb <-sum(tpar)+1.0078
  
  #print (parent)
  #yions<-matrix(parent, ncol = len)
  outb[j,] <- fragb
}

################################################################

#####            generando la matriz de fragmentos teoricos     ##########
yions<-t(outy)
outb <- apply(outb, 2, rev)
bions<-t(outb)

matrix_frags <- cbind(bions,yions)

#####################################################################

###  preparando la matriz de errores ppm   TRAGE   #############

frv1 <- matrix(fr[,1])
rr <- nrow(frv1)
exp <- matrix(frv1[-rr,]) 
exp <- matrix(rep(exp,ncol(matrix_frags)), ncol = ncol(matrix_frags))
exp <- apply(exp, 2, as.numeric)

###  preparando la matriz de masas experimentales con carga 2   TARGET    ######


double <- matrix(2,nrow(exp),1)
proton <- matrix(1.0078,nrow(exp),1)
frv1dc <- matrix(frv1[-rr,])
frv1dc <- apply(frv1dc, 2,as.numeric)
expdc = frv1dc * double - proton
expdc <- matrix(rep(expdc,ncol(matrix_frags)), ncol = ncol(matrix_frags))

###  preparando la matriz de masas experimentales con carga 3  TARGET   ######


triple <- matrix(3,nrow(exp),1)
dproton <- matrix(2.0156,nrow(exp),1)
frv1tc <- matrix(frv1[-rr,])
frv1tc <- apply(frv1tc, 2,as.numeric)
exptc = frv1tc * triple - dproton
exptc <- matrix(rep(exptc,ncol(matrix_frags)), ncol = ncol(matrix_frags))

####  Preparando la matriz de ppm para las cargas 1 2 3 TARGET   ######


ppm <-((exp - matrix_frags)/matrix_frags)*1000000
ppm <- abs(ppm)
ppmdc <- ((expdc - matrix_frags)/matrix_frags)*1000000
ppmdc <- abs(ppmdc)
ppmtc <- ((exptc - matrix_frags)/matrix_frags)*1000000
ppmtc <- abs(ppmtc)

ppmfinal <- pmin(ppm,ppmdc,ppmtc)

min <- transform(ppmfinal, minv=apply(ppmfinal[,],1, min, na.rm = TRUE))
min <- min[,ncol(min)]
min <- as.matrix(min)

zoom <- ifelse(min>50,sample(50:90,replace=T),min)


fppm <- ppmfinal[rowSums(ppmfinal < 1300) >= 0.01*ncol(ppmfinal),] 

fppmFALSE <- matrix(50, nrow = len*2, ncol = len*2)
# write.matrix(ppm, file = "C:\\Users\\proteomica\\Desktop\\R_projects\\ppm.txt", sep = ",")
# write.matrix(ppmdc, file = "C:\\Users\\proteomica\\Desktop\\R_projects\\ppmdc.txt", sep = ",")
# write.matrix(ppmtc, file = "C:\\Users\\proteomica\\Desktop\\R_projects\\ppmtc.txt", sep = ",")
#write.matrix(ppmfinal, file = "C:\\Users\\proteomica\\Desktop\\R_projects\\ppmfinal.txt", sep = ",")
#write.matrix(min, file = "C:\\Users\\proteomica\\Desktop\\R_projects\\min.txt", sep = ",")


fppm <-t(fppm)

if(ncol(fppm)==0)
{
  fppm = fppmFALSE
}

fppm <- fppm[,ncol(fppm):1]

####    selecciona el máximo error en ppm que se va a tener en cuenta     #####

fppm <- ifelse(fppm>50,50,fppm)  ### convierte todos los valores mayores de 50 en 50



if(is.null(dim(fppm)))
{
  fppm = fppmFALSE
}

##  preparando la matriz de intensidades experimentales  TARGET  ###

frv2 <- matrix(fr[,2])
int <- matrix(frv2[-rr,])
int2 <- matrix(frv2[-rr,])

#############      calculando Qscore       ######################
err = 9

qss <- transform(ppmfinal, minv=apply(ppmfinal[,],1, min, na.rm = TRUE))
qss <- qss[,ncol(qss)]
qss <- as.matrix(qss)

qss[qss>err] <- 0

qss <- as.numeric(qss)
int2 <- as.numeric(int2)

qscore <- cbind(qss, int2) 

qscoreFALSE <- matrix(21, nrow = 2, ncol = 2)

row_sub = apply(qscore, 1, function(row) all(row !=0 ))
qscore <- qscore[row_sub,]

if(length(qscore) == 2) 
{
  qscore <- qscoreFALSE
}

Escore <- sum(qscore[,2])
matched_ions = length(qscore)/2


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

