##
#Network model for DOC and CO2 in an Optimized Channel Network with lakes addition
#version 3 May 2021
#domvachon@gmail.com
##


#loading the required libraries and other scripts
library(deSolve)
library(rootSolve)
library(OCNet)
library(ggplot2)
library('igraph')
source("calcRate_network.r")


##################################
#Network configuration/generation#
##################################

  #create landscape with one outlet, positioned in the center, each pixel are 500 m wide
  set.seed(5) #change seed to explore other network configurations
  OCN_lake <- create_OCN(50, 50, outletPos = 25, cellsize = 250, typeInitialState = "T")
  draw_simple_OCN(OCN_lake)

  #slope affect stream gas exchange, zMin could potentially affect the temperature if(altitude ~ temperature function)
  OCN_lake <- landscape_OCN(OCN_lake, slope0 = 0.005, zMin = 0)

  #aggregate network
  OCN_lake <- aggregate_OCN(OCN_lake, thrA = 500000)
  
  draw_subcatchments_OCN(OCN_lake)
  points(OCN_lake$AG$X,OCN_lake$AG$Y, pch = 21, col = "blue", bg = "blue")
  
  #Catchment areas
  A=OCN_lake$AG$AReach #total CA at node i
  a=OCN_lake$SC$ALocal #direct CA at node i
  Atot=OCN_lake$CM$A #total catchment area in m2
  Atot/1000000 #total catchment area in km2
  pA=a/Atot #proportion of each sub-catchment to total catchment area
  
  slope <- OCN_lake$AG$slope #approximation of each stream slope
  elev <- OCN_lake$AG$ZReach #system elevation (m)
  streamOrder <- OCN_lake$AG$streamOrder #Strahler order

  draw_thematic_OCN(streamOrder, OCN_lake, chooseAggregation = NULL, colPalette=hcl.colors, colLevels=c(1,max(streamOrder),max(streamOrder)))
  title("Strahler order")
  
  ##
  #network attributes
  ##

  nNodes <- OCN_lake$AG$nNodes #number of nodes
  outlet <- OCN_lake$AG$outlet #ID of the outlet node

  ##connectivity matrix
  #determines which is the upstream system
  W <- as.matrix(OCN_lake$AG$W)

  ##cumulative connectivity matrix
  #determines which are all the upstream systems
  Z <- matrix(0, nNodes, nNodes)
  for (n in 1:nNodes){Z[,n][unlist(OCN_lake$AG$upstream[n])]=1}
  

##############
#Adding lakes#
##############
  
  lakeProp <- 0.15 #proportion of nodes that are lakes
    
  #assigning lakes per stream order
  SO1 <- which(streamOrder==1); SO2 <- which(streamOrder==2); SO3 <- which(streamOrder==3); SO4 <- which(streamOrder==4)
  #randomly assigning lakes in nodes of SO 1 to 3
  lakes <- c(sample(SO1,round(lakeProp*length(SO1))),sample(SO2,round(lakeProp*length(SO2))),sample(SO3,round(lakeProp*length(SO3))))
  
  #identifying system type
  syst <- rep("stream",nNodes); syst[lakes]="lake"
  
  #lake positions
  lakesX=OCN_lake$AG$X[lakes]
  lakesY=OCN_lake$AG$Y[lakes]
 
  #Show where lakes are
  draw_simple_OCN(OCN_lake,thrA = 500000,easyDraw=TRUE)
  points(OCN_lake$AG$XReach[lakes],OCN_lake$AG$YReach[lakes], pch = 21, col = "blue", bg = "blue")
  

###############
### Climate ### 
###############

  Tw=20 #temperature
  U10=2 #wind speed ms-1

  #discharge as a function of specific runoff and catchment area
  Kq=0.0025 #specific runoff in m d-1
  
  Q=A*Kq #discharge
  qL=a*Kq #lateral discharge
  qUp=Q-qL #upstream discharge
  qProp=qL/Q #proportion between laterla and upstream
   
#################################
#### morphology and hydrology ###
#################################

  #volume and CA of each node
  length=OCN_lake$AG$leng
  #Relationships from Raymond et al. 2012
  width=exp(2.56+0.423*log(Q/86400)) 
  depth=exp(-0.895+0.294*log(Q/86400))
  
  SA=length*width #surface area 
  vol=length*width*depth #volume
  
  draw_thematic_OCN(length, OCN_lake, chooseAggregation = NULL, colPalette=hcl.colors,colLevels=c(1,max(length),max(length)))
  title("Stream length")

  #enlarging and deepening for lakes
  SA[lakes]=10^((log10(A[lakes]/10000)-1.23)/1.05)*10000 #Scaling lake area with their catchment size (Walter et al. 2020, Aquatic Sciences)
  vol[lakes]=exp(1.39+1.12*log(SA[lakes]/1000000))*1000000 #Sobek et al. 2011 Inland Water
  depth[lakes]=vol[lakes]/SA[lakes]
  width[lakes]=SA[lakes]/length[lakes]
  
  #Water residence time (days)
  tau=vol/Q
  tauCumul=sum(vol)/Q[OCN_lake$AG$outlet]   #average residence time of the whole network
  alfa=qL/Q #relative contribution of lateral to upstream water loading
  WRTcumul<-(apply(vol*Z,2,sum))/Q #cumulative residence time, i.e., the average residence time of all the upstream systems
  

  #unzeroing the outlet SA and vol to avoid errors.
  SA[OCN_lake$AG$outlet] <- 0.001
  vol[OCN_lake$AG$outlet] <- 0.001
  

#################################
# gas exchange parameterization #
#################################
  
  #hydrological parameters
  vel=length/tau  #m d-1
  V=vel/(24*60*60)
  eD = 8.31*slope*V
  eD[is.nan(eD)] <- 0
  
  #Piston velocity for CO2 (m d-1)
  #k600
  #mixed model using Ulseth for streams and Klaus and Vachon (2020) for lakes
  
  k600=rep(0,nNodes)
  for (k in 1:nNodes){
    
    if (eD[k]>0.02) {
      k600[k]=exp(1.18*log(eD[k])+6.43)
    } else {
      k600[k]=exp(0.35*log(eD[k])+3.10) 
    }
  }
  
  k600[lakes] <- ((1.266+(0.328*log10(SA[lakes]/1000000)+1.581)*U10-(0.066*6.906755))*24/100) #Klaus and Vachon 2020
  
  ScCO2=1923.6-125.06*Tw+4.3773*Tw^2-0.085681*Tw^3+0.00070284*Tw^4 #Wanninkhof 2014
  kCO2=k600*(ScCO2/600)^(-0.5)
  
 
#######################################
#C input concentrations and reactivity#
#######################################

 
  #terrestrial C input
  DOCin <- 25 #input C concentration, g C m-3 
  CO2in <- 10 #input C concentration, g C m-3 
  
  #Lateral inputs, g C d-1
  #fixed
  #DocL <- propDOC*tCexport*a #by catchment export
  DocL <- qL*DOCin #by stream C concentration
 
  #baseline CO2 lateral C inputs, g C d-1
  #Henry constant as function of temp
  pKH=(-60.2409+93.4517*(100/(273.15+Tw))+23.3585*log((273.15+Tw)/100))/(log(10))
  Kh <- 10^(pKH)
  CO2Sat <- 400*Kh  #mmol m-3
  
  #Co2L <- (1-propDOC)*tCexport*a
  Co2L <- qL*CO2in
  exCo2L <- Co2L-qL*(CO2Sat*12.01/1000)
  
 
  #Doc reactivity 1: Slow 
  kv=0.5
  ka=0.5
  
  #Doc reactivity 2: Fast
  #kv=0.25
  #ka=0.0025
  
  #activation energy (J mol-1)
  Ae <- 50 #J/mol
  
  
###################
#running the model#
###################


  #defining fixed parameters
  parameters <- list(A=A,a=a,Z=Z,W=W,Kq=Kq,kv=kv,ka=ka,slope=slope,vol=vol,length=length,width=width,depth=depth,
                     SA=SA,U10=U10,Tw=Tw,Ae=Ae,nNodes=nNodes,lakes=lakes)

  Cpool <- c(5*vol,5*vol)

  #Run the simulation using deSolve

  out <- steady.1D(y = Cpool, method = "runsteady", jactype = 'sparse', func = calcRate_network, parms = parameters, nspec=2, names=c("Doc","Co2"))
  
#########
#RESULTS#
#########


  #results output
  Doc <- out$y[,1]/vol
  Co2 <- out$y[,2]
  pCO2 <- out$pCO2
  CO2flux <- out$atm/SA #g C m-2 d-1
  DOCresp <- out$r_doc/SA #g C m-2 d-1
  
  #plot results
  draw_thematic_OCN(Doc, OCN_lake, chooseAggregation = NULL, colPalette=hcl.colors)
  title("DOC concentration (mg L-1)")
  
  draw_thematic_OCN(pCO2, OCN_lake, chooseAggregation = NULL, colPalette=hcl.colors)
  title("pCO2 (uatm)")
  
  draw_thematic_OCN(CO2flux, OCN_lake, chooseAggregation = NULL, colPalette=hcl.colors)
  title("CO2 flux (g C m-2 d-1)")
  
  draw_thematic_OCN(DOCresp, OCN_lake, chooseAggregation = NULL, colPalette=hcl.colors)
  title("DOC mineralization (g C m-2 d-1)")
  
  
