#CALCULATING GROWING DEGREE DAYS
#FUNCTION: degreedays() v 0.3 has been tested for accuracy but has not been optimized for performance.
#Elise Larsen, 2013-11-05, SESYNC
#Contact: ealarsen.umd [at] gmail.com

##CHANGES: Default LDT set to 10 (assumed temperatures given in degrees Celsius).

#Single sine wave approximation for growign degree days with lower & upper development thresholds.
#Based on Degree-Days: The Calculation and Use of Heat Units in Pest Management FROM Zalom et al. Uc Coop. Extension 1983 (2m-12/83-WC/FB), Baskerville & Emin 1969.

##Define Terms

#tmin is minimum recorded temperature - required parameter
#tmax is maximum recorded temperature - required parameter
#ldt is [species-specific] lower developmental threshold - optional parameter; if no ldt is given then the lower development threshold is set to 10 degrees (10C is a commonly used threshold based on corn.)
#udt is [species-specific]upper developmental threshold - optional parameter; if no udt is given then the upper development threshold is set to whichever is higher, the maximum recorded temperature or the lower development threshold.

#Define subfunctions used in degreedays()
myalpha = function(tmin,tmax) {
  return((tmax - tmin) / 2)
}
paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )

mybeta = function(tmin,tmax) {
  return((tmax + tmin) / 2)
}
paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )

t1function = function(tmin,tmax,ldt) {
  alpha1<-myalpha(tmin,tmax)
  beta1<-mybeta(tmin,tmax)
  return(asin((ldt - beta1) / alpha1))
}
paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )

t2function = function(tmin,tmax,udt) {
  alpha1<-myalpha(tmin,tmax)
  beta1<-mybeta(tmin,tmax)
  return(asin((udt - beta1) / alpha1))
}
paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )


##BEGIN degreedays() FUNCTION


#Calculated degree days using single sine wave approximation
degreedays=function(tmin,tmax,ldt,udt){
  
  ## CHECK FOR APPROPRIATE PARAMETERS
  
  if(missing(tmin) | missing(tmax)) {
    warning("No calculation: Missing Temperature Parameter(s)", immediate. = TRUE)
    return(NA)
    stop
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  ## Check that temps are numeric
  if(!is.numeric(tmin) | !is.numeric(tmax)){
    warning("No calculation: Temperature(s) not numeric", immediate. = TRUE)
    return(NA)
    stop
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  ## Check that tmax > tmin
  if(tmax < tmin){
    warning("No calculation: Maximum Temperature Less Than Minimum Temperature", immediate. = TRUE)
    return(NA)
    stop
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  ### SET LDT AND UDT IF MISSING
  
  # 'missing()' is used to test whether a value was specified as an argument.  
  
  if(missing(ldt)) {  #Use tmin for lower development threshold is given
    ldt<-10 #assumes Celsius
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  if(missing(udt)) { # 'missing()' is used to test whether a value was specified as an argument. Code for Degree days when no upper development threshold is given
    udt<-max(tmax,ldt)
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  ##CHECK LDT, UDT PARAMETERS
  
  if(udt < ldt){
    warning("Upper Temperature Threshold Less Than Lower Temperature Threshold", immediate. = TRUE)
    return(NA)
    stop
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  
  ##BEGIN CALCULATIONS
  
  dd<-NA
  
  ## BEGIN WITH LDT <= TMIN
  if(ldt<=tmin) {
    ##Calculation: Case 1:  LDT < UDT < Tmin < Tmax: d = UDT-LDT
    if(udt<=tmin) {
      dd<-udt-ldt
    } 
    paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
    
    #Calculation: case 5: LDT < Tmin < UDT < Tmax: 
    if(udt>tmin & udt<tmax) {
      theta2 = t2function(tmin,tmax,udt)
      dd<-(1/pi) * ((mybeta(tmin,tmax) - ldt) * (pi/2 + theta2) + (udt - ldt) * (pi/2 - theta2) - myalpha(tmin,tmax) * (cos(theta2)))
    } 
    paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
    
    ##Calculation: case 3: LDT < Tmin < Tmax < UDT: d = beta - LDT  
    if(udt>=tmax) {
      dd<-mybeta(tmin,tmax) - ldt
    }
    paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
    
  } else { ##  LDT > TMIN
    
    ##Calculation: Case 2:  Tmin < Tmax < LDT < UDT: d = 0
    if(ldt>tmax) {
      dd<-0
    } 
    paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
    
    #Calculation: case 4: Tmin < LDT < Tmax < UDT: 
    if(ldt<=tmax & udt>tmax) {
      theta1<-t1function(tmin,tmax,ldt)
      dd<-(1/pi) *( (mybeta(tmin,tmax) - ldt) * (pi/2 - theta1) + myalpha(tmin,tmax) * (cos(theta1)) )
    } 
    paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
    
    #Calculation: case 6: Tmin < LDT < UDT < Tmax: 
    if(ldt<=tmax & udt<=tmax) {
      theta1 = t1function(tmin,tmax,ldt)
      theta2 = t2function(tmin,tmax,udt)
      dd<-(1/pi) * ((mybeta(tmin,tmax) - ldt) * (theta2 - theta1) + myalpha(tmin,tmax) * (cos(theta1) - cos(theta2)) + (udt - ldt) * (pi/2 - theta2) )
    }
    paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
    
  }
  paste(unclass(as.POSIXlt(Sys.time()))$min,unclass(as.POSIXlt(Sys.time()))$sec )
  
  return(dd)
  stop 
  
}



#END FUNCTION



######Pseudocode used in development
#Check that LDT < UDT (if given); Check that Tmax < Tmin

##alpha = (Tmax - Tmin) / 2
##beta = (Tmax + Tmin) / 2


#case 1: LDT < UDT < Tmin < Tmax: UDT is less than Tmin: d = UDT-LDT
#case 2: Tmin < Tmax < LDT < UDT: LDT is greater than Tmax: d = 0
#case 3: LDT < Tmin < Tmax < UDT: d = beta - LDT
#case 4: Tmin < LDT < Tmax < UDT: 
##theta1 = sin-1 ( (LDT - beta) / alpha)
## d = (1/pi) ( (beta - LDT) (pi/2 - theta1) + alpha (cos(theta1)) )
#case 5: LDT < Tmin < UDT < Tmax: 
##theta2 = sin-1 ( (UDT - beta) / alpha)
## d = (1/pi) ( (beta - LDT) (pi/2 + theta2) + (UDT - LDT) (pi/2 - theta2) - alpha (cos(theta2)) )
#case 6: Tmin < LDT < UDT < Tmax: 
##theta1 = sin-1 ( (LDT - beta) / alpha)
##theta2 = sin-1 ( (UDT - beta) / alpha)
## d = (1/pi) ( (beta - LDT) (theta2 - theta1) + alpha (cos(theta1) - cos(theta2)) + (UDT - LDT)(pi/2 - theta2) )


