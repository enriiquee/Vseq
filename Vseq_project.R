
#x = 76933
#sn <- paste0("SCANS=", x)
##########  extracting Charge and XCorr    ############

sub <- subset(sql, FirstScan == x)
if(nrow(sub)==0){
  sub <- data.frame(1, Sequence=1, ScoreValue= 1.11,RetentionTime=1,
                    Charge=1, ProteinDescriptions=1, ScoreID=9, MonoisotopicMass=1000)
}
sub_xc <- subset(sub, ScoreID == 9 )
sub_xc_max <- sub_xc[,"ScoreValue"]
Xcorr <- max(sub_xc_max)

sub_ch_max <- subset(sub_xc, ScoreValue == Xcorr)

if (nrow(sub_ch_max)>1){
  sub_ch_max = sub_ch_max[1,]
}

Charge <- sub_ch_max[,"Charge"]

###########################################

###   Encontrar la secuencia del SCAN en SQL   y el RT    ####  
seql <- sub_ch_max[,"Sequence"]
seql <- as.character(seql)
head(seql)
rt <- sub_ch_max[,"RetentionTime"]
prot <- sub_ch_max[,"ProteinDescriptions"]
#########################################################

###    seleccionar UNA secuencia     ####
seq <- seql
if(length(seql)==2){
  seq = seql[1]
}
seq2 <- strsplit(seq, NULL)[[1]]
seq2 <- paste(rev(seq2), collapse='')

#############################################            


fr <- fread(infile,skip=sn); as.data.frame(fr)

fr<-as.data.frame(fr)
frnames <- c("mass", "intens")
colnames(fr) <- frnames

len<-nchar(seq)
out_of = len*2

yy <- substr(seq, 1, len)
print(yy)
yn <- strsplit(yy, "")
z <- as.data.frame(yn)
z<-as.matrix(z)

pas<- paste0("c(", as.character(z), ")")

tre<-unlist(lapply(pas, 
                   function(x) sum(eval(parse(text=x)))))

parental <-sum(tre)+19.0178

if (z[1]=="k") {
  
  parental <-sum(tre)+19.0178+isobLab
  
}

mim <- sub_ch_max[,"MonoisotopicMass"]  

DeltaMass = mim - parental
DMSS = DeltaMass


parentaldm = parental + DeltaMass
DeltaMassdm = mim - parentaldm

### aqui estaba el calculo de secuencia DECOY

###### preparando el spectrum ###################



frn1 <- matrix(fr[,1])
rr <- nrow(frn1)
exp <- matrix(frn1[-rr,]) 

exp <- apply(exp, 2, as.numeric)

frn2 <- matrix(fr[,2])
rr <- nrow(frn2)
exp2 <- matrix(frn2[-rr,]) 

exp2 <- apply(exp2, 2, as.numeric)
zero <- matrix(0, nrow=nrow(exp2), ncol=ncol(exp2))
ccu <- matrix(0.01, nrow=nrow(exp), ncol= ncol(exp))
idea <- exp - ccu
m_z <- c(matrix(c(idea, exp), 2, byrow = T)) 
RelInt <- c(matrix(c(zero, exp2), 2, byrow= T))

bind <- cbind(m_z, RelInt)
ccu2 <- matrix(0.01, nrow = nrow(bind))

bind2 <- bind[,1] + ccu2
morezero <- matrix(0, nrow=nrow(bind2))
bind3 <- c(matrix(c(bind[,1], bind2), 2, byrow = T)) 
bind4 <- c(matrix(c(RelInt, morezero), 2, byrow= T))
spec <- cbind(bind3, bind4)


rownamesspec <- c(1:nrow(spec))      #####################         NEW, for spec intensity correction
rownames(spec) <- rownamesspec


maxRelInt <- max(RelInt)
meanRelInt <- median(as.numeric(frn2[1:nrow(frn2)-1]))
stdRelInt <- sd(as.numeric(frn2[1:nrow(frn2)-1]))
#spec_correction = maxRelInt/meanRelInt


normRelInt <- (as.numeric(frn2[1:nrow(frn2)-1])-meanRelInt)/stdRelInt

pRelInt <- pnorm(normRelInt, mean = 0, sd = 1)

normspec <- cbind(intensity=as.numeric(frn2[1:nrow(frn2)-1]), pRelInt)
normspec <- subset(normspec, pRelInt > 0.81)


maxRelIntpos <- subset(spec, bind4 == maxRelInt)
maxposition <- rownames(maxRelIntpos)
spec_correction = maxRelInt/mean(normspec[,1])

spec2 <- spec_correction*spec[,2]
spec2 <- cbind(spec[,1], spec2)
spec2[as.numeric(maxposition),2] = bind4[as.numeric(maxposition)]

spec = spec2
spec <- ifelse(spec > maxRelInt,maxRelInt-13, spec)


###############################################################################
#####    preparando matriz de fragmentos tesricos    ###


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

###  preparando la matriz de errores ppm   TARGET   #############

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

#############################################################################
##############################################################################
#############################################################################
############################################################################
##################  DATOS  PARA  LOS  CALCULOS  DELTAMASS   #################

######      calculate y series with DeltaMass #############


ipar=c(2:len-1)
outydm <- matrix(NA, nrow=len, ncol= nrow(fr)-1)

for (ipar in 1:len)
{
  yydm <- substr(seq, ipar, len)
  #print(yy)
  yndm <- strsplit(yydm, "")
  zdm <- as.data.frame(yndm)
  zdm<-as.matrix(zdm)
  tdm<-unlist(lapply(paste0("c(", as.character(zdm), ")"), 
                     function(x) sum(eval(parse(text=x)))))
  fragydm <-sum(tdm)+19.0178+DeltaMass
  
  #print (parent)
  #yions<-matrix(parent, ncol = len)
  outydm[ipar,] <- fragydm
}



#################################################################

###################    calculate b series con deltamass######### 


j=c(2:len-1)

outbdm <- matrix(NA, nrow=len, ncol= nrow(fr)-1)

for (j in 1:len)
{
  
  bbdm <- substr(seq2, j, len)
  #print(bb)
  bndm <- strsplit(bbdm, "")
  zdm <- as.data.frame(bndm)
  zdm<-as.matrix(zdm)
  tdm<-unlist(lapply(paste0("c(", as.character(zdm), ")"), 
                     function(x) sum(eval(parse(text=x)))))
  fragbdm <-sum(tdm)+1.0078+DeltaMass
  
  #print (parent)
  #yions<-matrix(parent, ncol = len)
  outbdm[j,] <- fragbdm
}

################################################################

#####            generando la matriz de fragmentos teoricos   
###################################con deltamass   ##########

yionsdm<-t(outydm)
outbdm <- apply(outbdm, 2, rev)
bionsdm<-t(outbdm)

matrix_fragsdm <- cbind(bionsdm,yionsdm)

#####################################################################

###  preparando la matriz de errores ppm   TARGET  with Deltamass 

frv1 <- matrix(fr[,1])
rr <- nrow(frv1)
expdm <- matrix(frv1[-rr,]) 
expdm <- matrix(rep(expdm,ncol(matrix_fragsdm)), ncol = ncol(matrix_fragsdm))
expdm <- apply(expdm, 2, as.numeric)

###  preparando la matriz de masas experimentales con carga 2   
###########################################TARGET  Deltamass  ######

double <- matrix(2,nrow(exp),1)
proton <- matrix(1.0078,nrow(exp),1)
frv1dcdm <- matrix(frv1[-rr,])
frv1dcdm <- apply(frv1dcdm, 2,as.numeric)
expdcdm = frv1dcdm * double - proton
expdcdm <- matrix(rep(expdcdm,ncol(matrix_fragsdm)), ncol = ncol(matrix_fragsdm))

###  preparando la matriz de masas experimentales con carga 3  TARGET  
####################################################Deltamass ######

triple <- matrix(3,nrow(exp),1)
dproton <- matrix(2.0156,nrow(exp),1)
frv1tcdm <- matrix(frv1[-rr,])
frv1tcdm <- apply(frv1tcdm, 2,as.numeric)
exptcdm = frv1tcdm * triple - dproton
exptcdm <- matrix(rep(exptcdm,ncol(matrix_fragsdm)), ncol = ncol(matrix_fragsdm))

####  Preparando la matriz de ppm para las cargas 1 2 3 TARGET  Deltamass ######

ppmdm <-((expdm - matrix_fragsdm)/matrix_fragsdm)*1000000
ppmdm <- abs(ppmdm)
ppmdcdm <- ((expdcdm - matrix_fragsdm)/matrix_fragsdm)*1000000
ppmdcdm <- abs(ppmdcdm)
ppmtcdm <- ((exptcdm - matrix_fragsdm)/matrix_fragsdm)*1000000
ppmtcdm <- abs(ppmtcdm)




####preparando los ablines......con Deltamass
# masasdm <- cbind(expdm[,1], mindm)
# matchesdm <- masasdm[rowSums(masasdm < 30) >= 0.01*ncol(masasdm),]

ppmfinaldm <- pmin(ppmdm,ppmdcdm,ppmtcdm)
parcialdm = ppmfinaldm

fppmdm <- ppmfinaldm[rowSums(ppmfinaldm < 300) >= 0.01*ncol(ppmfinaldm),] 
fppmFALSEdm <- matrix(50, nrow = len*2, ncol = len*2)

### selecciona el maximo error en ppm que se va a tener en cuenta  Dmass #####

fppmdm <- ifelse(fppmdm>50,50,fppmdm)  ### convierte todos los valores mayores de 50 en 50




if(is.null(dim(fppmdm)))
{
  fppmdm = fppmFALSEdm
}


################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
################################################################


frags1 = list=c(paste0("b", 1:len))
frags2 = list=c(paste0("y", len:1))
frags1dc = list=c(paste0("b", 1:len, "++"))
frags2dc = list=c(paste0("y", len:1, "++"))
frags1dm = list=c(paste0("b", 1:len, "*"))
frags2dm = list=c(paste0("y", len:1, "*"))
frags1dcdm = list=c(paste0("b", 1:len,"*", "++"))
frags2dcdm = list=c(paste0("y", len:1,"*", "++"))
frags <- append(frags1,frags2)
fragsdc <- append(frags1dc, frags2dc)
fragsdm <- append(frags1dm,frags2dm)
fragsdcdm <- append(frags1dcdm, frags2dcdm)

frags <- as.vector(frags)
colnames(ppmdm) <- frags
colnames(ppmdcdm) <- frags
colnames(ppmtcdm) <- frags


# 

############# asignacion de los iones en el espectro



##############################
###############################
#################################

asign <- cbind(matrix_frags[1,], frags)
asigndc <- cbind(((matrix_frags[1,]+1.0078)/2), fragsdc)
asigndm <- cbind(matrix_fragsdm[1,], fragsdm)
asigndcdm <- cbind(((matrix_fragsdm[1,]+1.0078)/2), fragsdcdm)

c_asign <- rbind(asign, asigndc, asigndm, asigndcdm)

if (DeltaMass<3){
  c_asign <- rbind(asign, asigndc)
}
c_asignc1 <- as.numeric(c_asign[,1])
c_asignc2 <- as.character(c_asign[,2])
c_asign <- data.frame(c_asignc1, c_asignc2)
####################################sigue despues cuando se declare el parametro matches


ppmfinal <- pmin(ppm,ppmdc,ppmtc)
parcial = ppmfinal

if (DeltaMass>=3){
  DeltaMass=DeltaMassdm
  ppmfinal <- pmin(ppm,ppmdc,ppmtc,ppmdm,ppmdcdm,ppmtcdm)
}
#fppmy <-exp[,1]     ################   new
#rownames(ppmfinal) <- fppmy   ###############   new

min <- transform(ppmfinal, minv=apply(ppmfinal[,],1, min, na.rm = TRUE))
min <- min[,ncol(min)]
min <- as.matrix(min)

zoom <- ifelse(min>50,sample(50:90,replace=T),min)


fppm <- ppmfinal[rowSums(ppmfinal < 300) >= 0.01*ncol(ppmfinal),] 

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



####preparando los ablines......
masas <- cbind(exp[,1], min)

matches <- as.data.frame(masas[rowSums(masas < 30) >= 0.01*ncol(masas),])
if (nrow(matches)==2){
  matches <- matrix(1:4, ncol = 2, nrow = 2)
}

##################sigue la asignacisn###############
matches_ions <- data.frame()
for (mi in 1:length(matches[,1])){
  for (ci in 1:length(c_asign[,1])){
    if (abs((matches[mi,1]-c_asign[ci,1])) <= 0.3){
      asigned <- data.frame(matches[mi,1], c_asign[ci,2])
      matches_ions<-rbind(matches_ions,asigned)
    }
  }
}
if (length(matches_ions)==0){
  matches_ions <- data.frame( c1=c(1.2,1.3), c2=c("bNI", "yNI"))
}


proofnames1 <- c("dta1","dta2")
colnames(matches_ions) <- proofnames1
fr[,1] <- as.numeric(fr[,1])

proof <- merge(matches_ions, fr, by.x = 'dta1', by.y='mass')
if (nrow(proof)==0){
  proof <- data.frame(matches_ions, dta3=c(fr[1,1],fr[2,1]))
}

proof[,3] <- as.numeric(proof[,3])*spec_correction

proof$intens[proof$intens > maxRelInt] <- maxRelInt-3

proofb <- filter(proof, grepl("b", dta2))
proofy <- filter(proof, grepl("y", dta2))

####    selecciona el maximo error en ppm que se va a tener en cuenta     #####

fppm <- ifelse(fppm>50,50,fppm)  ### convierte todos los valores mayores de 50 en 50
########################
##########################
##########################
deltascolname <- c(1:ncol(ppmfinal))
deltasrowname <- c(1:nrow(ppmfinal))
colnames(parcial) <- deltascolname
colnames(parcialdm) <- deltascolname
colnames(ppmfinal) <- deltascolname
rownames(parcial) <- deltasrowname
rownames(parcialdm) <- deltasrowname
rownames(ppmfinal) <- deltasrowname



parcial <- ifelse(parcial<50,2,parcial)
parcialdm <- ifelse(parcialdm<50,3,parcialdm)
pppmfinal <- ifelse(ppmfinal<=300,1,ppmfinal)
pppmfinal <- ifelse(pppmfinal>300,0,pppmfinal)
parcial <- ifelse(parcial>50,0,parcial)
parcialdm <- ifelse(parcialdm>50,0,parcialdm)



deltamplot <- pmax(parcialdm, parcial, pppmfinal)

deltamplot <- deltamplot[rowSums(deltamplot > 0) >= 0.01*ncol(deltamplot),] 

if(is.null(nrow(deltamplot))){
  deltamplot = parcial
}
if(nrow(deltamplot)==0){
  deltamplot = parcial
}

filterrowname <- c(1:nrow(deltamplot))
rownames(deltamplot) <- filterrowname



rplot <- data.frame(NA)
cplot <- data.frame(NA)
for (ki in 1:nrow(deltamplot)){
  for (kj in 1:ncol(deltamplot)){
    if (deltamplot[ki,kj] == 3){
      rplot <- data.frame(cbind(rplot,rownames(deltamplot)[ki]))
      cplot <- data.frame(cbind(cplot,colnames(deltamplot)[kj]))
    }}}
rplot <- t(as.matrix(rplot[,-1]))
cplot <- t(as.matrix(cplot[,-1]))

deltaplot <- cbind(rplot,cplot)
deltav1 <- as.numeric(deltaplot[,1])
uves1 <- rep(nrow(deltamplot), length(deltav1))
deltav1 = uves1 - deltav1
deltav2 <- as.numeric(deltaplot[,2])

if(length(deltav1)==0){
  deltav1=0
  deltav2=0
}

if(length(deltav1)>0){
  if(deltav1==0){
    deltav1=1
  }
}


if(is.null(dim(fppm)))
{
  fppm = fppmFALSE
}




z <-max(fppm)


#Aqui estaba el calculo de matriz DECOY

##  preparando la matriz de intensidades experimentales  TARGET  ###

frv2 <- matrix(fr[,2])
int <- matrix(frv2[-rr,])
int2 <- matrix(frv2[-rr,])

#fmint <- mint[rowSums(mint > 3000000) >= 0.01*ncol(mint),] 


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

Escore <- sum(qscore[,2]/1000000)
matched_ions = length(qscore)/2

tqueryscan <- paste0( x, ",") 
pepmass <- subset(tquery, SCAN == tqueryscan)

massspec <- pepmass[,"MZ"]
chargespec <- pepmass[,"CHARGE"]
specpar <- paste( "MZ=", massspec,",", "Charge=", chargespec)

#############################################
#############################################

#####  Ascore calculation   ####################

# rAsc <- data.frame(NA)
# cAsc <- data.frame(NA)
# for (ki in 1:nrow(deltamplot)){
#   for (kj in 1:ncol(deltamplot)){
#     if (deltamplot[ki,kj] == 2){
#       rAsc <- data.frame(cbind(rAsc,rownames(deltamplot)[ki]))
#       cAsc <- data.frame(cbind(cAsc,colnames(deltamplot)[kj]))
#     }}}
# rAsc <- t(as.matrix(rAsc[,-1]))
# cAsc <- t(as.matrix(cAsc[,-1]))
# 
# Ascoor <- cbind(rAsc,cAsc)      ###############    posiciones de matriz V sin delta mass
# uves2 <- rep(ncol(fppm), length(coor1))
# coor1 <- uves2 - as.numeric(Ascoor[,1])
# coor2 <-  as.numeric(Ascoor[,2])

AsB <- data.frame()
AsY <- data.frame()

for (dvar0 in 1:length(deltav2)){
  if (deltav2[dvar0] <= len){
    AsBpar <- data.frame(cbind(row= deltav2[dvar0], col= deltav1[dvar0]))
    AsB <- rbind(AsBpar, AsB)
  }
  if (deltav2[dvar0] > len){
    AsYpar <- data.frame(cbind(row= deltav2[dvar0], col= deltav1[dvar0]))
    AsY <- rbind(AsYpar, AsY)
  }
}

if (nrow(AsB)==0){
  AsB <- data.frame(row=0, col=0)
}
if (nrow(AsY)==0){
  AsY <- data.frame(row=0, col=0)
}

AsB <- AsB[order(AsB$row),c(1,2)]
AsY <- AsY[order(AsY$row),c(1,2)]
BDAG <- data.frame(cbind(AsB), dist= 0, value= 0,consec= 0, int= 0, error= 0, count=c(1:nrow(AsB)))

for (dvar0 in nrow(BDAG):1){
  BDAG[nrow(BDAG),3] = BDAG[nrow(BDAG),2]
  dvar1 = dvar0-1
  BDAG[dvar1,3] = abs(BDAG[dvar0,2] - BDAG[(dvar1),2])
  for (dvar1 in 1:dvar0){
    if (BDAG[dvar0,3] <= 7){
      BDAG[dvar1,4] = BDAG[dvar0,4]+1
    }
    else if (BDAG[dvar0,3] > 7){
      BDAG[dvar1,4] = 0
    }
    BDAG[dvar0,5] = BDAG[dvar0,4]+1
  }}
if (nrow(AsB)==1){
  BDAG <- data.frame(cbind(AsB), dist= 0, value= 0,consec= 0, int= 0, error= 0, count=c(1:nrow(AsB)))
}



YDAG <- data.frame(cbind(AsY), dist= 0, value= 0,consec= 0, int= 0, error= 0, count=c(nrow(AsY):1))
YDAG <- YDAG[order(YDAG$count),c(1:8)]

for (dvar0 in nrow(YDAG):1){
  YDAG[nrow(YDAG),3] = YDAG[nrow(YDAG),2]
  dvar1 = dvar0-1
  YDAG[dvar1,3] = abs(YDAG[dvar0,2] - YDAG[(dvar1),2])
  for (dvar1 in 1:dvar0){
    if (YDAG[dvar0,3] <= 7){
      YDAG[dvar1,4] = YDAG[dvar0,4]+1
    }
    else if (YDAG[dvar0,3] > 7){
      YDAG[dvar1,4] = 0
    }
    YDAG[dvar0,5] = YDAG[dvar0,4]+1
  }}
if (nrow(AsY)==1){
  YDAG <- data.frame(cbind(AsY), dist= 0, value= 0,consec= 0, int= 0, error= 0, count=c(1:nrow(AsY)))
}



PTMprob <- strsplit(seq, NULL)[[1]]


BDAGmax <- subset(BDAG, BDAG[,4]==max(BDAG[,4]))
if (nrow(BDAG)==1){
  BDAGmax[1,1] = 0
}
YDAGmax <- subset(YDAG, YDAG[,4]==max(YDAG[,4]))
YDAGmax <- YDAGmax[,1] - len

##############################################################
############################################################

####################   Getting Survey scan information   ########################


#dtapath <- "E:\\NT_scaf1\\4730\\"
# dtafiles <- data.frame(as.character(list.files(dtapath)))

if (data_type == "DdS"){
  
  pauta <- as.data.frame(gregexpr("\\.", dtafiles[1,1]))  
  remove <- as.data.frame(substr(dtafiles[,1], as.numeric(pauta[1,1]+1), 50))
  remove <- data.frame(as.numeric(sub('\\..*', '', remove[,1])))
  dtamatrix <- cbind(dtafiles, remove)
  dtacolnames <- c("1","2")
  colnames(dtamatrix) <- dtacolnames
  dtamatrix <- dtamatrix[order(dtamatrix$`2`),]
  
  dtamatrix <- subset(dtamatrix, dtamatrix[,2]< x)
  dtamatrix <- dtamatrix[nrow(dtamatrix),]
  
  dta_prec <- read.csv(paste0(dtapath, dtamatrix[1,1]), sep = " ")
  
  prec_bar <- subset(dta_prec, as.numeric(dta_prec[,1]) < massspec+3)
  prec_bar <- subset(prec_bar, as.numeric(prec_bar[,1]) > massspec-3)
  colnames(prec_bar) <- c("mass", "intensity")
  
  if (nrow(prec_bar) == 0){
    prec_bar <- data.frame("mass"= 1, "intensity"=1)
  }
  
  zero_dta <- matrix(0, nrow=nrow(prec_bar), ncol=1)
  ccu_dta <- matrix(0.01, nrow=nrow(prec_bar), ncol= 1)
  
  
  idea_dta <- prec_bar[,1] - ccu_dta
  m_z_zoom <- c(matrix(c(idea_dta, prec_bar[,1]), 2, byrow = T)) 
  RelInt_dta <- c(matrix(c(zero_dta, prec_bar[,2]), 2, byrow= T))
  
  bind_dta <- cbind(m_z_zoom, RelInt_dta)
  ccu2_dta <- matrix(0.01, nrow = nrow(bind_dta))
  
  bind2_dta <- bind_dta[,1] + ccu2_dta
  morezero_dta <- matrix(0, nrow=nrow(bind2_dta))
  bind3_dta <- c(matrix(c(bind_dta[,1], bind2_dta), 2, byrow = T)) 
  bind4_dta <- c(matrix(c(RelInt_dta, morezero_dta), 2, byrow= T))
  spec_dta <- cbind(bind3_dta, bind4_dta)
  
  specpar_dta <- paste( "MZ=", massspec,",", "Charge=", chargespec, ",Scan=", dtamatrix[,2])
}

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


###  pintando las matrices  ppm de TARGET y spectrum    #####

rownames(fppm) <- frags
#rownames(fppmdec) <- frags


mainT <- paste(seq, DMSS, sep = "+")
#mainD <- paste("DECOY", decoy, sep = "---")


z <-max(fppm)


#plot.new()

varNamePath <- "E:\\R_graphs\\"
varconc<- paste0(varNamePath,seq,"_",x,".png")
#tiff(filename = varconc, width = 54 ,height = 36, units = c("cm"), res = 90)
png(file = varconc, width = 54 ,height = 36, units = c("cm"), res = 90)

#par(mfrow=c(2,3), oma=c(1,1,1,3))

m_paroma <- matrix(c(1,4,1,4,2,4,2,5,3,5,3,5), nrow = 2, ncol = 6)
##set up the plot
layout(m_paroma)


plot(zoom,int, log = "y",xlim = c(-3,90),asp = 1,
     xlab = "error in ppm______________________  >50", ylab = "intensity(log)",
     pch=21, bg="lightblue", cex=3, col = "purple", cex.lab=2)
abline(v=err, col="blue", lty=2)
title(prot, cex.main=1.5, col.main= "darkblue")


my_palette <- colorRampPalette(c("red","green","blue","orange","grey"))(z)
levels<-c(0,1,3,5,7,9,12,15,18,21,30:z)
#ndm <-  1:nrow(namesdm)
#ndm <- matrix(ndm)
#ndm2 = c(3,5)

p1 <- levelplot(fppm, col.regions=my_palette, levels=levels, main= list(label=mainT, 
                                                                        cex= 2),
                xlab=list(label="b series --------- y series", cex=2),
                ylab=list(label="large--Exp.masses--small", cex= 2), aspect = "fill",
                panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = 1:len, col = "white", lty = 2)
                  panel.abline(v = len+1:len, col = "white", lty = 2)
                  if(DMSS > 1){
                    grid.points(deltav2,deltav1, pch=24)
                  }
                })

if (data_type == "DdS"){
  plot(m_z_zoom ,RelInt_dta, cex=0, lines(spec_dta, col="darkblue", lty=1, lwd=2),
       xlab="m/z", ylab = "Intensity", cex.lab=2) 
  title(specpar_dta, cex.main=2, col.main= "darkblue")
  abline(v=massspec, col="green", lty=3, lwd=5)
}

if (data_type == "DiS"){
  plot(c(100, 200), c(300, 450), xlab = "", ylab = "")
  title("DiS experiment", cex.main=2, col.main= "darkblue")
}

# scales=list(x=list(cex=1), y=list(cex=1),
#                 xlab=list(cex=2), ylab=list(cex=2)))


print(p1, split=c(2,2,2,2), more=TRUE)

plot.new()

textbox(c(0,0.6), 1, c("FirstScan=",x,"Xcorr=",Xcorr,
                       "Charge=",Charge,"SQL-Sequence=",seql,"RT=",rt,
                       "DeltaM=", DeltaMass, "Escore=", Escore, "m.Mass=", mim)
        ,cex= 2.1,col="blue", border = FALSE)
if (DMSS > 1){
  textbox(c(0.6,0.9), 1, c("PTM pinpointing:","Bseries=",PTMprob[BDAGmax[,1]],BDAGmax[,1],"Yseries=",PTMprob[YDAGmax],YDAGmax)
          ,cex= 2.1,col="red", border = FALSE)}
if (DMSS < 1) {
  textbox(c(0.6,0.9), 1, c("unmodified peptide")
          ,cex= 2.1,col="red", border = FALSE)}


plot(m_z ,RelInt, cex=0, lines(spec, col="darkblue", lty=1, lwd=2),
     xlab="m/z", ylab = "Relative Intensity", cex.lab=2) 
title(specpar, cex.main=2, col.main= "darkblue")
abline(v=matches[,1], col="orange", lty=2)
for (pri in 1:nrow(proof)){
  
  text(proofb[pri,1], y= proofb[pri,3], labels = proofb[pri,2], cex = 2, col = "red")
  text(proofy[pri,1], y= proofy[pri,3], labels = proofy[pri,2], cex = 2, col = "blue")
  
}


dev.off()

params <- data.frame(x, seq,mim, Charge, Xcorr, Escore, rt, DeltaMass, 
                     prot, matched_ions, out_of, massspec)



