#Supporting Information
#Open source feature detection for non-target LC-MS analytics
#Christian Dietrich, Arne Wick, Thomas A. Ternes

#licensed under the GNU GPL version 3 or later

#This algorithm is designed to process centroided high-resolution MS data provided as mzML or mzXML

#before running this code:

#please make sure that you have installed the following packages:
#xcms: https://bioconductor.org/packages/release/bioc/html/xcms.html
#Rcpp: https://cran.r-project.org/web/packages/Rcpp/index.html or type: install.packages("Rcpp")
#RcppArmadillo: https://cran.r-project.org/web/packages/RcppArmadillo/index.html or type: install.packages("RcppArmadillo")
#make sure that this file and the "peakpicking2.cpp" file are in the same folder, otherwise adjust the file path in line 22

#please adjust the peak picking parameters in lines 106-117

library(xcms)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("peakpicking2.cpp")



#open .mzML or .mzXML using xcms: 
dataFile <- file.choose()
rawData <- xcmsRaw(dataFile,includeMSn=TRUE)



peakpicking_SingleXIC <- function(mz, rawData, mz_step,rt_min_scan,rt_max_scan,sn,int_threshold,peakwidth_min,peakwidth_max,maxPeaksPerSignal,precursormzTol){
  
  maxima <- NULL

  XIC_intensity <- rawEIC(rawData, mzrange = c(mz,mz+mz_step))
  XIC_intensity <- XIC_intensity$intensity
  
  maxima <- peakPicking_cpp(intensity = XIC_intensity, scantime = rawData@scantime, min_intensity = int_threshold, sn = sn, peakwidth_min = peakwidth_min, peakwidth_max = peakwidth_max, rt_min_scan = rt_min_scan, rt_max_scan = rt_max_scan, maxPeaksPerSignal = maxPeaksPerSignal)
  
  if (nrow(maxima) > 0) {
    for (j in 1:nrow(maxima)) {
      exactmass <- 0
      mass_spectrum <- NULL
      mass_spectrum_j <- NULL
      
      #extract mass spectra for each scan of the peak and calculate the intensity-weighted average of the exact m/z
      for (jj in (maxima[j,8]):(maxima[j,9])) {
        mass_spectrum_j <- getScan(rawData, jj, mzrange = c(mz,mz+mz_step))
        if (jj == maxima[j,5]) {
          maxima[j,3] <- max(mass_spectrum_j[,2])
          if (nrow(mass_spectrum_j) > 0) maxima[j,3] <- maxima[j,3]-maxima[j,15] 
        }
        if (nrow(mass_spectrum_j) > 0) {
          mass_spectrum_j <- mass_spectrum_j[which.max(mass_spectrum_j[,2]),,drop=FALSE]
          mass_spectrum <- rbind(mass_spectrum,mass_spectrum_j)
          exactmass <- exactmass+mass_spectrum_j[1,1]*mass_spectrum_j[1,2]
        }
      }
      
      if (nrow(mass_spectrum) > 0){
        exactmass <- exactmass/sum(mass_spectrum[,2])
      }
      #only keep those peaks with m/z in the "center" of the analyzed m/z range
      if ((exactmass == 0) | (exactmass < mz+mz_step/4) | (exactmass > mz+mz_step/4*3)) exactmass = 0 
      maxima[j,1] <- exactmass
      ms2 <- which((rawData@msnRt > maxima[j,6]) & (rawData@msnRt < maxima[j,7]) & (abs(rawData@msnPrecursorMz-exactmass) <= exactmass/1000000*precursormzTol))
      if (length(ms2) > 0) maxima[j,17] <- ms2[which.min(abs(rawData@msnRt[ms2]-maxima[j,2]))]
    }
  }
  
  maxima <- maxima[maxima[,1] > 0,,drop = FALSE]
  maxima[,16] <- mz
  if (nrow(maxima) > 0) return(maxima)
}


getmzbins <- function(mzmin,mzmax) {
  i <- 1
  mzbin <- mzmin
  while (mzbin[i] < mzmax) {
    i <- i+1
    #example for Orbitrap systems (mass resolution ~ sqrt(1/mz))
    mzbin[i] <- mzbin[i-1]+mzbin[i-1]*sqrt(1E-6/mzbin[i-1])
  }
  return(mzbin)
}


deriveoptimummzstep <- function(rawData, probs) {
  mzdifferences <- NULL
  for (i in 1:length(rawData@scanindex)) {
    msspektrum <- getScan(rawData, i)
    mzdifferences <- c(mzdifferences,-quantile(-diff(msspektrum[,1],lag=1),probs=probs))
  }
  return(min(mzdifferences))
}




peakPicking <- function(rawData) {
  
  peaklist <- NULL
  
  #peak picking parameters (to be adjusted individually):
  mz_min <- 100 #minimum m/z
  mz_max <- 1200 #maximum m/z
  mz_step <<- 0.02 #m/z bin size
  rt_min <- 120 #minimum retention time in seconds
  rt_max <- 1200 #maximum retention time in seconds
  sn <- 3 #signal/noise threshold
  int_threshold <- 1 #intensity threshold
  peakwidth_min <- 5 #minimum peak width in seconds
  peakwidth_max <- 60 #maximum peak width in seconds
  maxPeaksPerSignal <- 10 #signals containing more consecutive maxima of similar size than this
  #threshold will be regarded as noise
  precursormzTol <- 20 #mass tolerance (ppm) to find ms2
  
  
  
  #convert min and max retention times into scans
  rt_min_scan <- min(which(rawData@scantime > rt_min)) 
  rt_max_scan <- max(which(rawData@scantime < rt_max))
  if (is.infinite(rt_min_scan)) rt_min_scan <- 1
  if (is.infinite(rt_max_scan)) rt_max_scan <- length(daten@scantime)
  
  linear_binning <- TRUE #TRUE: fixed mz_step (line 101) is used. Check "deriveoptimummzstep(rawData,quantiles)", e.g. with quantiles=0.9
                         #FALSE: non-linear mz_step is used, e.g. as a function of m/z. Adjust "getmzbins" function accordingly.
  
  if (linear_binning == TRUE) {
  #loop over all m/z bins and pick peaks
  
   peaklist_list <- lapply(X = seq(mz_min,mz_max, by = mz_step*0.5), FUN = function(X)
     {peakpicking_SingleXIC(mz = X, rawData =  rawData, mz_step = mz_step,rt_min_scan = rt_min_scan,rt_max_scan = rt_max_scan, sn=sn,int_threshold=int_threshold,peakwidth_min=peakwidth_min,peakwidth_max=peakwidth_max,maxPeaksPerSignal=maxPeaksPerSignal,precursormzTol=precursormzTol)
   })
   peaklist <- as.data.frame(rbind(peaklist, do.call(rbind, peaklist_list)))
   
  } else {
    mzbins <- getmzbins(mz_min,mz_max)
    mz_step <- diff(mzbins)
    for (i in 1:length(mzbins)-1) {
      peaklist <- rbind(peaklist,peakpicking_SingleXIC(mz = mzbins[i], rawData =  rawData, mz_step = mz_step[i],rt_min_scan = rt_min_scan,rt_max_scan = rt_max_scan, sn=sn,int_threshold=int_threshold,peakwidth_min=peakwidth_min,peakwidth_max=peakwidth_max,maxPeaksPerSignal=maxPeaksPerSignal,precursormzTol=precursormzTol))
    }
  }
  
  
 
  
  colnames(peaklist) <- c("mz","RT","Intensity","XIC_Intensity","Scan","LeftendRT","RightendRT","Leftendscan","Rightendscan","NoiseScans", "NoiseDeviation","Area","FWHM_left","FWHM_right","Baseline","XICStartMass","MS2scan")
  
  peaklist <- peaklist[order(peaklist[,1]),]
  
  peaklist <- peaklist[peaklist[,1] > 0,]
  return(peaklist)
}


#run peak picking:
system.time({
  peaklist <- peakPicking(rawData)
})

#Plot a random feature from the peaklist:
peaknumber <- sample(nrow(peaklist),1)

noiseRTleft <- rawData@scantime[peaklist$Leftendscan[peaknumber]-(peaklist$Rightendscan[peaknumber]-peaklist$Leftendscan[peaknumber])*10]/60
noiseRTright <- rawData@scantime[peaklist$Rightendscan[peaknumber]+(peaklist$Rightendscan[peaknumber]-peaklist$Leftendscan[peaknumber])*10]/60
if (is.na(noiseRTleft)) noiseRTleft <- 0
if (length(noiseRTleft) > 1) noiseRTleft <- 0
if (is.na(noiseRTright)) noiseRTright <- max(rawData@scantime)/60
if (length(noiseRTright) > 1) noiseRTright <- 0
XIC <- rawEIC(rawData,mzrange=c(peaklist$mz[peaknumber]-mz_step/2,peaklist$mz[peaknumber]+mz_step/2))

plot(x=rawData@scantime/60,y=XIC$intensity, xlim = c(noiseRTleft,noiseRTright), ylim = c(0,(peaklist$XIC_Intensity[peaknumber]+peaklist$Baseline[peaknumber])*1.1),type = "l", main = paste0("XIC of m/z = ",round(peaklist$mz[peaknumber],4)," (",round(peaklist$mz[peaknumber]-mz_step/2,4),"-",round(peaklist$mz[peaknumber]+mz_step/2,4),")"),xlab="RT / min",ylab="Intensity",xaxs = "i",yaxs="i")
  
lines(x = c(peaklist$LeftendRT[peaknumber]/60,peaklist$LeftendRT[peaknumber]/60), y = c(peaklist$Baseline[peaknumber],(peaklist$XIC_Intensity[peaknumber])+peaklist$Baseline[peaknumber]), type = "l", col = "red")
lines(x = c(peaklist$RightendRT[peaknumber]/60,peaklist$RightendRT[peaknumber]/60), y = c(peaklist$Baseline[peaknumber],(peaklist$XIC_Intensity[peaknumber])+peaklist$Baseline[peaknumber]), type = "l", col = "red")
text(x = peaklist$RightendRT[peaknumber]/60, y = peaklist$XIC_Intensity[peaknumber]*0.9, labels = "Peak boundaries", col = "red", pos = 4)
lines(x = c(noiseRTleft,noiseRTright), y=c(peaklist$Baseline[peaknumber],peaklist$Baseline[peaknumber]), type = "l", col = "red")
text(x = noiseRTleft, y = peaklist$Baseline[peaknumber], labels = "Baseline", col = "red", adj = c(0,0))
lines(x = c(noiseRTleft,noiseRTright), y=c(peaklist$Baseline[peaknumber]+peaklist$NoiseDeviation[peaknumber]/2,peaklist$Baseline[peaknumber]+peaklist$NoiseDeviation[peaknumber]/2), type = "l", col = "blue")
lines(x = c(noiseRTleft,noiseRTright), y=c(peaklist$Baseline[peaknumber]-peaklist$NoiseDeviation[peaknumber]/2,peaklist$Baseline[peaknumber]-peaklist$NoiseDeviation[peaknumber]/2), type = "l", col = "blue")
text(x = noiseRTright, y = peaklist$Baseline[peaknumber]+peaklist$NoiseDeviation[peaknumber]/2, labels = "Noise", col = "blue", adj = c(1,0))
lines(x = c(peaklist$FWHM_left[peaknumber]/60,peaklist$FWHM_left[peaknumber]/60), y = c(peaklist$Baseline[peaknumber],(peaklist$XIC_Intensity[peaknumber]/2)+peaklist$Baseline[peaknumber]), type = "l", col = "green")
lines(x = c(peaklist$FWHM_right[peaknumber]/60,peaklist$FWHM_right[peaknumber]/60), y = c(peaklist$Baseline[peaknumber],(peaklist$XIC_Intensity[peaknumber]/2)+peaklist$Baseline[peaknumber]), type = "l", col = "green")
text(x = peaklist$FWHM_right[peaknumber]/60, y = peaklist$XIC_Intensity[peaknumber]/2, labels = "FWHM", col = "green", pos = 4)
  
