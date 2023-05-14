library(lubridate)
library(magrittr)
library(MSnbase)
library(xcms)
library(plotly)
library(dplyr)
library(openxlsx)

extract_pre <- function(x, sdmz, ppm) {
  if(!is.na(x)){
    if (any (
      sapply (sdmz, function(mz) {
        any(mz-(ppm * mz/1e+06) < x 
            & x < mz +(ppm * mz/1e+06))}))){
      TRUE
    } else { FALSE }
  } else { FALSE }
}
classfinder<-function(mz,ClassPeak,intensity,err.ppm=5){
  matches<-massmatch(mz,ClassPeak,err.ppm=err.ppm)
  if(any(!is.na(matches))){
    if(sum(!is.na(matches))==1){
      return(matches[!is.na(matches)])
    }
    return(PCfilter(mz,intensity,ClassPeak))
  }
  return("NOTHING")
}
massmatch<-function(observed,expected,err.ppm=5){
  sapply(observed,function(x){
    hit=x>expected*(1-err.ppm*10^-6)&x<expected*(1+err.ppm*10^-6)
    if(sum(hit)==1){
      return(names(expected)[hit])
    }
    return(NA)
  })
}
PCfilter<-function(mz,intensity,ClassPeak){
  mzintensitydf<-data.frame(mz, intensity)
  mzintensitydf$matches<-massmatch(mz,ClassPeak)
  mzintensitydf<-mzintensitydf[!is.na(mzintensitydf$matches),]
  notPCdf<-mzintensitydf[mzintensitydf$matches!="PC",]
  PCdf<-mzintensitydf[mzintensitydf$matches=="PC",]
  
  if(nrow(PCdf)==0){
    return(mzintensitydf$matches[which.max(mzintensitydf$intensity)])
  }
  if(any(notPCdf$intensity/PCdf$intensity>=1)){
    return(mzintensitydf$matches[mzintensitydf$matches!="PC"][which.max(notPCdf$intensity)])
  }
  if(any(notPCdf$intensity/PCdf$intensity < 1 & notPCdf$intensity/PCdf$intensity < 0.01)){
    return("PCPE")
  }
  return("PC")
}
PLoverTMT<-function(intensity,Chain, class, Channel){
  PLintensity<-suppressWarnings(max(intensity[which(is.na(Chain) & is.na(Channel) & class != "PC")]))
  PCintensity<-suppressWarnings(max(intensity[which(is.na(Chain) & is.na(Channel) & class == "PC")]))
  TMTintensity<- suppressWarnings(max(intensity[which(Channel %in% names(TMT))]))
  
  if(TMTintensity<0){
    return("N")
  }
  
  if (!any(PCintensity > 0)) {
    if( PLintensity > TMTintensity/100) {
      return("Y")
    } else {
      return("N")
    }
  }
  if(any(PCintensity > 0) & any(PLintensity/PCintensity > 0.1) & any(PLintensity > TMTintensity/50)) {
    return("Y")
  }
  if(PCintensity >= max(intensity)*0.7) {
    return("Y")
  }
  return("N")
  
}
PrecursorMatch <- function(precursorMz,PLdatabaseElite,class) {
  
  if (class %in% c("PC", "PE", "PG", "PS", "PA", "PI")) {
    return(group_pre(precursorMz, 
                     sd_vector = PLdatabaseElite,
                     class = class))
  }
  return("NOTHING")
  
}
group_pre <- function(pre_vector, sd_vector, class, err.ppm = 5) {
  
  sapply(pre_vector,function(x){
    res<-sd_vector$Exactmz[x < sd_vector$Exactmz*(1+err.ppm*10^-6) &
                             x > sd_vector$Exactmz*(1-err.ppm*10^-6) &
                             sd_vector$Class == class]
    if(length(res)!=1){
      if(length(res) > 1) {
        res<-res[which.min(abs(res-x))]
        Mol <- sd_vector$Molename[sd_vector$Class == class & sd_vector$Exactmz == res]
        
        unique(Mol)
      }else{
        NA
      }
    }else{
      Mol <- sd_vector$Molename[sd_vector$Class == class & sd_vector$Exactmz == res]
      unique(Mol)
    }
  }
  )
}

TMT <- c(TMT126 = 126.1276, TMT127 = 127.1245, TMT128 = 128.1344,
         TMT129 = 129.1312, TMT130 = 130.1412, TMT131 = 131.1380)
Reporterion = TMT

ClassPeak <- c(PC = 278.1767, PE = 513.2994, PS = 557.2893, 
               PA = 841.5298, PG = 544.2939, PI = 632.3098)


PLdatabaseElite <- read.xlsx("~/Documents/R script/PLdatabase.xlsx", colNames = F)
SidedatabasePL <- read.xlsx("~/Documents/R scipt/SidedatabasePL.xlsx", colNames = T)
FAfragment <- read.xlsx("~/Documents/R script/FAfragment_complement.xlsx", colNames = T)
subFAfragment<- FAfragment$Fmass
names(subFAfragment) <- FAfragment$Chain

#----------start here-----------------------
data <- readMSData("~/Documents/MS data/20230427_Liver_6plex_111111.mzXML", pdata = NULL, verbose = F,
                   centroided. = F, smoothed. = NA, mode = "onDisk")
dataTF <- sapply(precursorMz(data), extract_pre, 
                 sdmz = PLdatabaseElite[,2], ppm = 5)
dataMS1match <- as(data[dataTF], "MSnExp")
dataMS1match = clean(dataMS1match)

dflist<-lapply(dataMS1match@assayData,function(x)
  data.frame(precursorMz=x@precursorMz,mz=x@mz,intensity=x@intensity))

dfID<-Reduce(rbind,dflist)
dfID$IDnumber<-rep(names(dflist),sapply(dflist,nrow))
dfID$Charge <- precursorCharge(dataMS1match)[match(dfID$IDnumber, names(precursorCharge(dataMS1match)))]
dfID$RT <- round(rtime(dataMS1match)[match(dfID$IDnumber, 
                                           names(rtime(dataMS1match)))]/60, digits = 1)
dfID$PreInt <- round(precursorIntensity(dataMS1match)[match(dfID$IDnumber, 
                                                            names(precursorIntensity(dataMS1match)))])

dfIDclass<-dfID%>%group_by(IDnumber)%>%
  mutate(class=classfinder(mz,ClassPeak,intensity))%>%ungroup()
dfIDclass$Chain = massmatch(dfIDclass$mz, subFAfragment)
dfIDclass$Channel = massmatch(dfIDclass$mz, Reporterion)

dfIDclass = dfIDclass[which(dfIDclass$class != "NOTHING"),] 

dfIDforPC = dfIDclass[which(dfIDclass$class == "PC"),]
dfIDforPC = dfIDforPC%>%group_by(IDnumber)%>%
  mutate(MaxIT = mz[which(intensity == max(intensity))])
dfIDforPC = dfIDforPC[round(dfIDforPC$MaxIT) %in% c(278, 184),]
dfIDforPC = dfIDforPC[,-11]

dfIDclass = rbind(dfIDforPC, dfIDclass[-which(dfIDclass$class == "PC"),])

dfIDclean = dfIDclass[which(massmatch(dfIDclass$mz,ClassPeak)==dfIDclass$class | 
                              !is.na(dfIDclass$Chain) | !is.na(dfIDclass$Channel)),]
dfIDclean <- dfIDclean[order(dfIDclean$intensity, decreasing = T),]
dfIDclean <- dfIDclean[!duplicated(dfIDclean[, c("IDnumber", "Chain", "Channel")]), ]
dfIDclean <- dfIDclean[dfIDclean$class %in% c("PC", "PE", "PA") & dfIDclean$Charge == 2 |
                         dfIDclean$class %in% c("PG", "PI", "PS") & dfIDclean$Charge == 1,]

dfIDfilter <- dfIDclean%>%group_by(IDnumber)%>%
  mutate(TF = PLoverTMT(intensity,Chain, class, Channel))

dfIDfilter <- dfIDfilter[-which(dfIDfilter$TF == "N"),]
dfIDfilter <- dfIDfilter%>%group_by(IDnumber)%>%
  mutate(Mol = PrecursorMatch(precursorMz,PLdatabaseElite,unique(class)))
dfIDfilter$precursorMz <- round(dfIDfilter$precursorMz, digits = 3)
dfIDfilter <- dfIDfilter[!is.na(dfIDfilter$Mol),]
dfIDfilter = dfIDfilter[, -11]
dfIDfilter <- dfIDfilter %>% group_by(class, Channel, Mol, Chain, RT, PreInt, IDnumber)%>%
  summarise(intensity = sum(intensity))
sort(unique(dfIDfilter$Mol))

PGdatabase = PLdatabaseElite[PLdatabaseElite$Class == "PG", ]
PGtest <- sapply(dfID$precursorMz, extract_pre, 
                 sdmz = PGdatabase[, 2], ppm = 5)

PGtestdf <- dfID[PGtest,]                                                                                                                   
PGtestdf$Chain = massmatch(PGtestdf$mz, subFAfragment)
PGtestdf$Channel = massmatch(PGtestdf$mz, Reporterion)
PGtestdf$class <- "PG"

PGtestdf = PGtestdf[which(!is.na(PGtestdf$Chain) | !is.na(PGtestdf$Channel)),] 
PGMol <- PGtestdf%>%group_by(IDnumber)%>%
  mutate(Mol = PrecursorMatch(precursorMz,PLdatabaseElite,
                              unique(class)))

PGMol <- PGMol %>% group_by(PreInt, class, Channel, Mol, Chain, RT, IDnumber)%>%
  summarise(intensity = sum(intensity))
sort(unique(PGMol$Mol))

PLtotal = rbind(dfIDfilter[-which(dfIDfilter$class == "PG"),], PGMol)

sideTMTre2 = PLtotal %>% group_by(PreInt, class, Mol, RT, Chain, IDnumber)%>%
  summarise(TMT126 = sum(intensity[which(Channel == "TMT126")]),
            TMT127 = sum(intensity[which(Channel == "TMT127")]),
            TMT128 = sum(intensity[which(Channel == "TMT128")]),
            TMT129 = sum(intensity[which(Channel == "TMT129")]),
            TMT130 = sum(intensity[which(Channel == "TMT130")]),
            TMT131 = sum(intensity[which(Channel == "TMT131")]))

collapTMT = sideTMTre2 %>% group_by(IDnumber) %>%
  mutate(TMT126 = max(TMT126),
         TMT127 = max(TMT127),
         TMT128 = max(TMT128),
         TMT129 = max(TMT129),
         TMT130 = max(TMT130),
         TMT131 = max(TMT131))

collapChain = merge(collapTMT,SidedatabasePL,by.x=c("Mol","Chain"),by.y=c("Mol", "FAchain"))
HightPeak = collapChain %>% arrange(desc(PreInt)) %>% 
  group_by(class, Mol, Chain)%>%
  slice(1:10)
collapchain = HightPeak %>% group_by(IDnumber, PreInt, class, Mol, RT, 
                                     TMT126, TMT127, TMT128, 
                                     TMT129, TMT130, TMT131)%>%
  summarise(Chain = paste(unique(Chain), collapse = ', '))

collapRT = collapchain %>% group_by(class, Mol, round(RT), Chain)%>%
  summarise(PreInt = sum(PreInt),
            TMT126 = sum(TMT126),
            TMT127 = sum(TMT127),
            TMT128 = sum(TMT128),
            TMT129 = sum(TMT129),
            TMT130 = sum(TMT130),
            TMT131 = sum(TMT131))

colnames(collapRT)[3] = "RT"

collapRT = collapRT[order(collapRT$Mol, collapRT$RT),]

write.xlsx(collapRT, "~/Dropbox/PL R scipt/DennyWD/liver128_collap.xlsx")


