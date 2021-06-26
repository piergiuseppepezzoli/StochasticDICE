
# per installare gdxrrw seguire queste istruzioni
# https://support.gams.com/gdxrrw:interfacing_gams_and_r
library(gdxrrw)

library(fda)

library(viridisLite)

library(ggplot2)

library(ggpubr)

library(gganimate)

library(animation)

library(purrr)

#-----------------FILENAMES-----------------------
setwd("/Users/piergiuseppepezzoli/Climate Change/Project/Our Own Code") 
#cartella del file gdx

igdx("/Library/Frameworks/GAMS.framework/Versions/34/Resources") 
#per trovare questa directory lanciare un file gams con   display "%gams.sysdir%";  

filename = "Stochastic7.gdx"
#nome file .gdx


#--------------------FUNCTIONS-----------------------
#varie helping functions

#trasforma la lista caricata con rgdx in una matrice più ordinata 

isolate = function(variable,as){
  
  if(dim(variable$val)[2] == 4){
    res = matrix(nrow = dim(as$val)[1], ncol = 4 )
    for(x in 1:dim(variable$val)[1]){
      for(y in 1:dim(as$val)[1]){
        if(all(variable$val[x,1:3] == as$val[y,])){
          
          res[y,1] = as.numeric(as$uels[[1]][(variable$val[x,1])])
          res[y,2] = as.numeric(as$uels[[2]][(variable$val[x,2])])
          res[y,3] = as.numeric(as$uels[[3]][(variable$val[x,3])])
          res[y,4] = as.numeric((variable$val[x,4]))
        }
      }
    }
  }
  if(dim(variable$val)[2] == 3){
    res = matrix(nrow = dim(variable$val)[1], ncol = 3 )
    for(x in 1:dim(variable$val)[1]){
      res[x,1] = as.numeric(as$uels[[2]][(variable$val[x,1])])
      res[x,2] = as.numeric(as$uels[[3]][(variable$val[x,2])])
      res[x,3] = as.numeric((variable$val[x,3]))
    }
  }
  
  return(res)
}

#a partire da una variabile "pulita" isola lo scenario specificato con st1 st2

isolate_scenario = function(cleanvar, st1, st2){
  res = NULL
  for(i in 1:dim(cleanvar)[1]){
    if(cleanvar[i,1] < st1 && cleanvar[i,1] < st2){
      if(all(cleanvar[i,2:3] == c(0,0))){
        res = c(res,cleanvar[i,4])
      }
    }
    if(cleanvar[i,1] >= st1 && cleanvar[i,1] < st2){
      if(all(cleanvar[i,2:3] == c(st1,0))){
        res = c(res,cleanvar[i,4])
      }
    }
    if(cleanvar[i,1] < st1 && cleanvar[i,1] >= st2){
      if(all(cleanvar[i,2:3] == c(0,st2))){
        res = c(res,cleanvar[i,4])
      }
    }
    if(cleanvar[i,1] >= st1 && cleanvar[i,1] >= st2){
      if(all(cleanvar[i,2:3] == c(st1,st2))){
        res = c(res,cleanvar[i,4])
      }
    }
  }
  return(res)
}

# restituisce valori max e min e scenario (0,0) in una matrice più ordinata

banda = function(clean_var){
  res = matrix(nrow = 50,ncol=2)
  for(t in 1:50){
    index = which(clean_var[,1] == t)
    res[t,1] = max(clean_var[index,4])
    res[t,2] = min(clean_var[index,4])
  }
  return(res)
}

#fa smoothing con basi spline

smooth_data = function(clean_var,nbasis){
  res = matrix(nrow = 50,ncol = 25)
  iter = 0
  for(st1 in seq(0,40,by=10)){
    for(st2 in seq(0,40,by=10)){
      iter = iter +1 
      res[,iter] = isolate_scenario(clean_var,st1,st2)
    }
  }
  basis.1 <- create.bspline.basis(rangeval=c(1,50),nbasis=10)
  data_W.fd.1 <- Data2fd(y = res,argvals = 1:50,basisobj = basis.1)
  return(data_W.fd.1)
}

# pulisce il report dopo montecarlo

clean_report = function(report){
  res = matrix(nrow = dim(report$val)[1], ncol = dim(report$val)[2])
  for(i in 1:dim(report$val)[1]){
    for(ii in 1:(dim(report$val)[2]-1)){
      res[i,ii] = (report$uels[[ii]][report$val[i,ii]])
    }
    res[i,6] = as.numeric(report$val[i,6])
  }
  
  return(res)
}

# estrae una variabile dal report pulito di montecarlo

variable_extract = function(rep_clean,variablename){
  res = rep_clean[which(rep_clean[,1] == variablename),2:6]
  return(res)
}

# estrae una specifica iterazione dalla variabile pulita

iteration_extract = function(variable_clean,iteration,total_iterations){
  iteration_string = paste('i', iteration, sep = '')
  res = variable_clean[which(variable_clean[,4] == iteration_string),c(1,2,3,5)]
  
  
  if(dim(variable_clean)[1] == 574*total_iterations){
    true_res = matrix(nrow = dim(res)[1],ncol = dim(res)[2])
    for(i in 1:dim(res)[1]-24){
      true_res[i,1] = as.numeric(res[i,1])
      true_res[i,2] = as.numeric(res[i,2])
      true_res[i,3] = as.numeric(res[i,3])
      true_res[i,4] = as.numeric(res[i,4])
    }
    return(true_res)
  }else{
    true_res = matrix(nrow = dim(res)[1]-24,ncol = dim(res)[2])
    for(i in 1:dim(res)[1]-24){
      true_res[i,1] = as.numeric(res[24+i,1])
      true_res[i,2] = as.numeric(res[24+i,2])
      true_res[i,3] = as.numeric(res[24+i,3])
      true_res[i,4] = as.numeric(res[24+i,4])
    }
    return(true_res)
  }
  
}
