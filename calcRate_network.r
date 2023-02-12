calcRate_network <- function(t, Cpool, parameters) {
  
  #Function to calculate the rate of C change
  #DV 1 May 2013, modified from CTS 2 Feb 2013
  #Arguments are:
  #  Cpool - CO2 and DOC masses in the lake at the current time step (g)
  #  parameters - a vector giving the values of the fixed parameters
  
  #convert  state variables and parameters as lists
  with(as.list(c(Cpool,parameters)),{
    
    #splitting pools
    Doc <- Cpool[1:nNodes] 
    Co2 <- Cpool[(nNodes+1):(2*nNodes)]
    
    #DOC and CO2 concentrations
    DocConc <- Doc/vol #g C m-3
    Co2Conc <- Co2/vol #g C m-3
    Co2m <- Co2Conc*1000/12 #mmol m-3
    
    
######################   
##Hydrological inputs#
######################
    
  
    #output from outflow, g C d-1
    DocOut <- Q*DocConc
    Co2Out <- Q*Co2Conc
    
    #input from upstream, g C d-1
    DocUp=apply(DocOut*W,2,sum)
    Co2Up=apply(Co2Out*W,2,sum)
 
   
#####################################
#Net atmospheric exchange of CO2
#positive sign indicates gain by lake
#####################################
    
    #Henry constant as function of temp
    pKH=(-60.2409+93.4517*(100/(273.15+Tw))+23.3585*log((273.15+Tw)/100))/(log(10))
    Kh <- 10^(pKH)
 
    #Piston velocity for CO2 (m d-1)
    #k600
    #mixed model using Ulseth for streams and Vachon and Prairie (2013) for lakes
    
    k600=rep(0,nNodes)
    for (k in 1:nNodes){
      
      if (eD[k]>0.02) {
        k600[k]=exp(1.18*log(eD[k])+6.43)
      } else {
        k600[k]=exp(0.35*log(eD[k])+3.10) 
      }
    }
    
    #k600[lakes]*((2.51+1.48*U10+0.39*U10*log10(SA[lakes]/1000000))*24/100) #Vachon and Prairie 2013
    
    k600[lakes] <- ((1.266+(0.328*log10(SA[lakes]/1000000)+1.581)*U10-(0.066*6.906755))*24/100) #Klaus and Vachon 2020
    
    ScCO2=1923.6-125.06*Tw+4.3773*Tw^2-0.085681*Tw^3+0.00070284*Tw^4 #Wanninkhof 2014
    kCO2=k600*(ScCO2/600)^(-0.5)
    
    #Calculate flux of CO2, mmol m-2 d-1
    CO2Sat <- 400*Kh  #mmol m-3
    pCO2 <- Co2m/Kh    #uatm
    
    flux <- kCO2*(Co2m-CO2Sat) #mmol m-2 d-1
    atm <- (flux*12/1000)*SA  #g d-1
    
    
##################
##DOC degradation#
##################
    
    #DOC age
    tauDoc <- Doc/(DocUp + DocL)
    docAge <- rep(0,nNodes)
    #for (n in 1:nNodes){
    #  tauDocCumul[n] <- sum(tauDoc*Z[,n])
    #}
    
    docAge <- apply(Doc*Z,2,sum)/apply(DocL*Z,2,sum)
    
    ##
    #DOC loss by bacterial consumption 
    #function of age
  
    k_doc <- kv/(ka+docAge)
 
     
    #DOC removal rates
    #Volumetric
    DOC_R <- Doc*k_doc
    
    #temperature dependence 
    #Arrhenius equation
    r_doc <- DOC_R*exp((-Ae/8.314)*((1/(Tw+273.15))-(1/293.15))) #g d-1
    
   
    
    ########################
    #Differencial equations#
    ########################
    
    #Differential equation A: DOC pool
    dDoc.dt <- DocUp + DocL - r_doc - DocOut
    
    #Differential equation A: CO2 pool
    dCo2.dt <- Co2Up + Co2L  + r_doc - atm - Co2Out
    
    #return rates for each pool
    list(c(dDoc.dt, dCo2.dt),pCO2=pCO2,k_doc=k_doc,DocUp=DocUp,Co2Up=Co2Up,r_doc=r_doc,DocOut=DocOut,atm=atm,Co2Out=Co2Out)
  })
}