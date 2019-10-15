## Idea:
# Make new directory, simulations_v4, put in all the Rscripts, and then read in all the data from _v3

# Folder structure:
# /output
# /gen.data/signals
# /gen.data/wave/ for the cleaned simulation data

## Functions required:
# Change effect size to proportion script, and cater for some constant movement to 1
# Integrate the simulation script inside (make it a function)
# Integrate WaveQTL pre-processing (noting that a QT IS performed here.)

## Trim out any fat that doesn't need to be in the script.

## Remember to create directories as required, or check that they exist


# Functions ---------------------------------------------------------------

# Copied from '/src/R/my_utils.R' from HJ's repo, multiscale_analysis.
estBetaParams <- function(mu.orig, var) {
  
  del.ix = ((mu.orig <= 0) | (mu.orig >= 1))
  if(sum(del.ix) > 0){
    mu = mu.orig[!del.ix]
  }else{
    mu = mu.orig
  }
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  
  if(sum(del.ix) > 0){
    alpha.orig = beta.orig = rep(NA, length(mu.orig))
    alpha.orig[!del.ix] = alpha
    beta.orig[!del.ix] = beta
  }else{
    alpha.orig = alpha
    beta.orig = beta
  }
  
  ## handle non-positive alpha or beta 
  invalid.para = which((alpha.orig <= 0) | (beta.orig <= 0))
  if(length(invalid.para) > 0){
    alpha.orig[invalid.para] = rep(NA, length(invalid.para))
    beta.orig[invalid.para] = rep(NA, length(invalid.para))
  }
  
  return(params = list(alpha = alpha.orig, beta = beta.orig))
}

sample.from.Binomial.with.Overdispersion <- function(num.sam, total.count, mu.sig, over.dispersion=NULL){
  
  invalid.entry = ((mu.sig < 0) | (mu.sig > 1))
  if(sum(invalid.entry) > 0){ stop("ERROR, mu.sig have some values outside of valid range [0, 1]")}
  
  if(is.null(over.dispersion)){
    return(matrix(data=rbinom(length(mu.sig)*num.sam, total.count, mu.sig), nr = num.sam, byrow = TRUE))
  }else{
    
    
    final.dat = matrix(data=NA, nr = num.sam, nc = length(mu.sig))
    
    # get alpha and beta
    resBeta = estBetaParams(mu.sig, over.dispersion)
    alpha = resBeta$alpha
    beta = resBeta$beta
    
    # for valid alpha and beta, sample data 
    del.ix = is.na(alpha)
    p.sig = rbeta(sum(!del.ix)*num.sam, alpha[!del.ix], beta[!del.ix])
    dat.V = rbinom(sum(!del.ix)*num.sam, total.count[!del.ix], p.sig) 
    final.dat[,!del.ix] = matrix(data=dat.V, nr = num.sam, byrow = TRUE)
    
    # for invalid alpha and beta, sample without over dispersion
    if(sum(del.ix) > 0){
      dat.IV = matrix(data=rbinom(sum(del.ix)*num.sam, total.count[del.ix], mu.sig[del.ix]), nr = num.sam, byrow = TRUE)
      final.dat[,del.ix] = matrix(data=dat.IV, nr = num.sam, byrow = TRUE)
    }
    
    return(final.dat)
    
  }
  
}

# Copied from WaveQTL's 'WaveQTL/R/WaveQTL_preprocess_funcs.R'
require("wavethresh")
fiter.WCs <- function(Data, meanR.thresh){
  
  if(is.vector(Data)){
    dim(Data) <- c(1, length(Data))
  }        
  numWCs = dim(Data)[2]
  J = log2(numWCs)
  
  Mean_R = rep(NA, numWCs)
  Mean_R[1] = mean(apply(Data, 1, sum))
  Mean_R[2] = Mean_R[1]
  
  posi = 3
  for(ss in 1:(J-1)){
    num_loc = 2^ss
    size_int = numWCs/num_loc
    st = (0:(num_loc-1))*size_int + 1
    en = st + size_int -1
    
    for(ll in 1:num_loc){
      Mean_R[posi] = mean(apply(Data[,st[ll]:en[ll]], 1, sum))
      posi = posi + 1
    }
  }
  
  filtered.WCs = rep(0, numWCs)
  wh = which(Mean_R > meanR.thresh)
  
  if(length(wh) > 0){
    filtered.WCs[wh] = rep(1, length(wh))
  }
  
  return(filtered.WCs)
}

FWT <- function(Data, filter.number=1, family="DaubExPhase"){
  
  if(is.vector(Data)){
    dim(Data) <- c(1, length(Data))
  }    
  T = dim(Data)[2]
  J = log2(T)
  N = dim(Data)[1]
  
  dat_D = matrix(data=NA, nr = N, nc = (T - 1))
  dat_C = rep(NA, N)
  
  dat_W = matrix(data=NA, nr = N, nc = T)
  
  for(j in 1:N){
    each_WT	= wd(Data[j,], filter.number=filter.number ,family=family) 
    dat_D[j,] = each_WT$D
    dat_C[j] = each_WT$C[length(each_WT$C)]
  }
  
  dat_W[,1] = dat_C
  dat_W[,2] = -dat_D[,(T -1)]
  
  st_input = 3
  en_posi = T - 2
  for(k in 1:(J-1)){
    st_posi = en_posi - 2^k + 1
    en_input = st_input + 2^k - 1
    dat_W[,st_input:en_input] = -dat_D[,st_posi:en_posi]
    en_posi = st_posi - 1
    st_input = en_input + 1
  }
  
  return(list(WCs = dat_W))
  
}

QT_randomTie <- function(x) {
  
  x.rank = rank(x, ties.method="random")
  return(qqnorm(x.rank,plot.it = F)$x)
}

corrected_forCovariates <- function(x, Covariates){
  return(lm(x~Covariates)$residuals)
}

Normalize.WCs <- function(WCs, Covariates=NULL){
  
  # QT to a standard normal distribution
  QT_dat = apply(WCs, 2, QT_randomTie)
  
  # correct for covariates and QT to a standard normal distribution. 
  if(!is.null(Covariates)){
    corrected_QT.dat = apply(QT_dat, 2, corrected_forCovariates, Covariates)
    QT_dat = apply(corrected_QT.dat, 2, QT_randomTie)
  }
  
  return(list(QT_WCs = QT_dat))
  
}

WaveQTL_preprocess <- function(Data, library.read.depth = NULL, Covariates = NULL, meanR.thresh = 2, no.QT = FALSE, filter.number=1, family="DaubExPhase"){
  
  
  if(is.vector(Data)){dim(Data)<- c(1,length(Data))} #change Data to matrix
  if(nrow(Data)==1){Covariates = NULL} #if only one observation, don't correct for covariates
  
  if(!is.null(Covariates)){
    if(is.vector(Covariates)){dim(Covariates)<- c(1,length(Covariates))} #change C to matrix
  }
  
  
  
  T = dim(Data)[2]
  J = log2(T)
  if(!isTRUE(all.equal(J,trunc(J)))){stop("Error: number of columns of Data must be power of 2")}
  N = dim(Data)[1]
  
  
  ### generate filtered.WCs
  if(!is.null(meanR.thresh)){
    filtered.WCs = fiter.WCs(Data, meanR.thresh)				
  }else{
    filtered.WCs = NULL
  }	
  
  
  ## corrected for read depth
  if(!is.null(library.read.depth)){
    DataC = Data/library.read.depth
  }else{
    DataC = Data
  }
  
  ## Wavelet Transform
  WCs = FWT(DataC, filter.number=filter.number, family=family)$WCs
  
  if(!no.QT){ # Normalize phenotype for testing
    ## Normalize WCs    
    if(N > 1){
      WCs = Normalize.WCs(WCs, Covariates)
    }
    WCs = WCs$QT_WCs
  }else{  # Normalize for effect size estimation 
    ## correct for Covariates 
    if(!is.null(Covariates)){
      WCs = apply(WCs, 2, corrected_forCovariates, Covariates)
    }
  }
  
  return(list(WCs = WCs, filtered.WCs = filtered.WCs))
} 

# Effect to prop script ---------------------------------------------------

effect_to_prop <- function(shrink_factor
                           ,effect.size.path
                           ,dir.path
                           ,datasetNum){
  numPBs = 1024
  for(ss in datasetNum){
    ## read effect size
    effect.size.path.this = paste0(effect.size.path, ss)
    if(file.exists(effect.size.path.this)==TRUE){
      effect.size = scan(effect.size.path.this, what = double())
      
      # Amend effect size by shrink factor
      effect.size <- (effect.size - 1) * shrink_factor + 1
      
      ## alt
      mu0 = 2/70/(effect.size + 1)
      mu1 = 2/70/(effect.size + 1)*effect.size
      this.path = paste0(dir.path, "alt.sig0.", ss)
      cat(mu0, file=this.path)
      this.path = paste0(dir.path, "alt.sig1.", ss)
      cat(mu1, file=this.path)
      ## null
      effect.size = rep(1, numPBs)
      mu0 = 2/70/(effect.size + 1)
      mu1 = 2/70/(effect.size + 1)*effect.size
      this.path = paste0(dir.path, "null.sig0.", ss)
      cat(mu0, file=this.path)
      this.path = paste0(dir.path, "null.sig1.", ss)
      cat(mu1, file=this.path)
      
    }
  }
}


# Simulation script -------------------------------------------------------

simulation_waveqtl <- function(
  seed # should be set to be equal to the dataset number
  ,geno.path
  ,raw.dat.path
  ,sig0.path
  ,sig1.path
  ,read.depth.ratio
  ,over.dispersion
  ,multipleSig
  ,wavelet.preprocess
  ,DESeq.preprocess
  ,wd.path
  ,output.dir.name
){
  
  if(multipleSig == 1){
    sig0.path = paste0(sig0.path, ".", seed)
    sig1.path = paste0(sig1.path, ".", seed)
    if(!is.null(raw.dat.path)){
      raw.dat.path = paste0(raw.dat.path, ".", seed)
    }
  }
  
  setwd(wd.path)
  
  #####################
  # sample data
  #####################
  
  # read genotype data
  genoD = round(as.numeric(scan(geno.path, what=double())))
  numSam = length(genoD)
  
  # read signal
  mu0 = scan(file = sig0.path, what = double())
  mu1 = scan(file = sig1.path, what = double())
  numBPs = length(mu0)
  
  # phenotype data
  phenoD = matrix(data=NA, nr= length(genoD), nc = numBPs)
  
  # let's sample!!!
  set.seed(seed)
  
  # Amended to only sample beta-binomial from our raw data
  ## read raw data from which we will sample
  raw.data = read.table(raw.dat.path, as.is = TRUE)
  raw.data.T = ceiling(as.numeric(apply(raw.data, 2, sum)))
  
  ## change read depth
  if(!is.null(read.depth.ratio)){
    raw.data.T = floor(raw.data.T*read.depth.ratio)
  }
  
  mu0.sig = mu0
  mu1.sig = mu1
  
  ## ok!!
  ## upper and lower bound!
  trunc.fun = function(x){
    x = max(0, x)
    return(min(1,x))
  }
  mu0.sig = sapply(mu0.sig, trunc.fun)
  mu1.sig = sapply(mu1.sig, trunc.fun)
  
  ## geno = 0
  wh0 = which(genoD == 0)
  if(length(wh0) > 0){
    phenoD[wh0,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh0), total.count = raw.data.T, mu.sig = mu0.sig, over.dispersion = over.dispersion)
  }  
  ## geno = 1
  wh1 = which(genoD == 1)
  if(length(wh1) > 0){
    phenoD[wh1,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh1), total.count = raw.data.T, mu.sig = mu1.sig, over.dispersion = over.dispersion)
  }
  
  #########################################
  # data preprocessing for wavelet analysis
  #########################################
  if(wavelet.preprocess){
    meanR.thresh = 2
    res = WaveQTL_preprocess(Data = phenoD, library.read.depth =NULL , Covariates = NULL, meanR.thresh = meanR.thresh)
    
    filteredWCs = res$filtered.WCs
    norm.DNase = res$WCs
    
    ##save normaized data and useWCs information in output.path
    # out.dir.path = paste0(wd.path, "wave/", output.dir.name, ".data/") 
    out.dir.path = paste0(wd.path, output.dir.name, ".data/")
    # if(!file.exists(out.dir.path)){
    if(!dir.exists(out.dir.path)){
      dir.create(out.dir.path)
    }
    
    this.path = paste0(out.dir.path, "DNase.", seed, ".txt")
    write.table(norm.DNase, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(out.dir.path, "use.", seed, ".txt")
    cat(filteredWCs, file = this.path)
  }
}


# Altogether now ----------------------------------------------------------

effect_to_prop_and_sim <- function(shrink_factor = 1
                                   ,datasetNum = 1
                                   ,effect.size.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/smooth.ratio.3."
                                   ,geno.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.dat"
                                   ,raw.dat.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/pheno.dat"
                                   ,output.data.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v4/data/"
                                   ,read.depth.ratio = 1
                                   ,over.dispersion = 1/70/70/10
                                   ,multipleSig = 1
                                   ,wavelet.preprocess = TRUE
                                   ,DESeq.preprocess = TRUE
                                   ,wd.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v4/"
                                   ,output.dir.name = "full"){
  effect_to_prop(shrink_factor
                 ,effect.size.path = effect.size.path
                 ,dir.path = output.data.path
                 ,datasetNum = datasetNum)
  
  # Alt
  simulation_waveqtl(seed = datasetNum
                     ,geno.path = geno.path
                     ,raw.dat.path = raw.dat.path
                     ,sig0.path = paste0(output.data.path,"alt.sig0")
                     ,sig1.path = paste0(output.data.path,"alt.sig1")
                     ,read.depth.ratio = read.depth.ratio
                     ,over.dispersion = over.dispersion
                     ,multipleSig = multipleSig
                     ,wavelet.preprocess = wavelet.preprocess
                     ,DESeq.preprocess = DESeq.preprocess
                     ,wd.path = paste0(wd.path,"alt/")
                     ,output.dir.name = output.dir.name
  )

  # Null
  simulation_waveqtl(seed = datasetNum
                     ,geno.path = geno.path
                     ,raw.dat.path = raw.dat.path
                     ,sig0.path = paste0(output.data.path,"null.sig0")
                     ,sig1.path = paste0(output.data.path,"null.sig1")
                     ,read.depth.ratio = read.depth.ratio
                     ,over.dispersion = over.dispersion
                     ,multipleSig = multipleSig
                     ,wavelet.preprocess = wavelet.preprocess
                     ,DESeq.preprocess = DESeq.preprocess
                     ,wd.path = paste0(wd.path,"null/")
                     ,output.dir.name = output.dir.name
  )
}
