rm(list = ls())

#-------------SUPPORT CODE --------------------

source('Stochastic_functions.R')


#--------------------ACTUALCODE-------------------

#per leggere una variabile usare rgdx come qui sotto
TATM = rgdx(filename,list(name='TATM'))
CPRICE = rgdx(filename,list(name='CPRICE'))
Y = rgdx(filename,list(name='Y'))
t = rgdx(filename,list(name='t'))

#importante leggere sempre as
as = rgdx(filename,list(name='as'))


#la funzione isolate lascia solo i nodi "giusti" e pulisce tutto
temperature = isolate(TATM,as)
temperature

#la funzione isolate_scenario richiede una variabile "pulita" e gli indici dello scenario
#restituisce il vettore corrispondente al singolo scenario
temperature1020 = isolate_scenario(temperature,10,20)
temperature1020





PROB = rgdx(filename,list(name='PROB'))
H1 = rgdx(filename,list(name='H1'))
H2 = rgdx(filename,list(name='H2'))

probabilities = isolate(PROB,as)
h1 = isolate(H1,as)
h2 = isolate(H2,as)



E_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'E'))
E = rgdx(filename,list(name='E'))
emissions = isolate(E,as)
iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(emissions,st1,st2),type = 'l',col='red',lwd = 5)
    }else{
      lines(isolate_scenario(emissions,st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(E_vanilla$val[,2],col = 'blue')


TATM_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'TATM'))
TATM = rgdx(filename,list(name='TATM'))
temperature = isolate(TATM,as)

iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(temperature,st1,st2),type = 'l',col='red',lwd = 5, ylim = c(0,5))
      
    }else{
      lines(isolate_scenario(temperature,st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(TATM_vanilla$val[,2],col = 'blue')



MIU_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'MIU'))
MIU = rgdx(filename,list(name = 'MIU'))
miu = isolate(MIU,as)

par(mfrow=c(1,1))
iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(miu[1:549,],st1,st2),type = 'l',col='red',lwd = 5)
    }else{
      lines(isolate_scenario(miu[1:549,],st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(MIU_vanilla$val[,2],col = 'blue')


Y_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'Y'))
y = isolate(Y,as)
iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(y,st1,st2),type = 'l',col='red',lwd = 5)
      
    }else{
      lines(isolate_scenario(y,st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(Y_vanilla$val[,2],col = 'blue')

ABAT_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'ABATECOST'))
ABATECOST = rgdx(filename,list(name = 'ABATECOST'))
abate = isolate(ABATECOST,as)
iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(abate[1:549,],st1,st2),type = 'l',col='red',lwd = 5)
      
    }else{
      lines(isolate_scenario(abate[1:549,],st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(ABAT_vanilla$val[,2],col = 'blue')

DAMAGE_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'DAMAGES'))
DAMAGES = rgdx(filename,list(name = 'DAMAGES'))
damages = isolate(DAMAGES,as)
iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(damages,st1,st2),type = 'l',col='red',lwd = 5, ylim =c(0,70))
      
    }else{
      lines(isolate_scenario(damages,st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(DAMAGE_vanilla$val[,2],col = 'blue')



DAMFRAC_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'DAMFRAC'))
DAMFRAC = rgdx(filename,list(name = 'DAMFRAC'))
damfrac = isolate(DAMFRAC,as)
iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(damfrac,st1,st2),type = 'l',col='red',lwd = 5, ylim = c(0,0.1))
      
    }else{
      lines(isolate_scenario(damfrac,st1,st2),col=viridis(16)[iter])
    }
  }
}
lines(DAMFRAC_vanilla$val[,2],col = 'blue')




#------------------BAND PLOT----------------

TATM = rgdx(filename, list(name = 'TATM'))
as = rgdx(filename, list(name ='as'))
temperature = isolate(TATM,as)


plot(banda(temperature)[,1],type = 'l', col= 'blue')
lines(banda(temperature)[,2], col = 'red')


lines(smooth.spline(banda(temperature)[,1]))




A = banda(temperature)
A = as.data.frame(A)
colnames(A) = c('max','min')
A['base'] = isolate_scenario(temperature,0,0)


basis = create.bspline.basis(rangeval=c(1,50), nbasis=20)
abscissa = 1:50

B = matrix(nrow = length(abscissa), ncol = dim(A)[2])
for(i in 1:dim(A)[2]){
  temp <- smooth.basis(argvals=abscissa, y=A[,i], fdParobj=basis)
  B[,i] <- eval.fd(abscissa, temp$fd) 
}

B = as.data.frame(B)
colnames(B) = c('max','min','base')



fig1 = ggplot(data=A, aes(x=abscissa, y=base, ymin=min, ymax=max)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5, fill = 'blue')

fig2 = ggplot(data=B, aes(x=abscissa, y=base, ymin=min, ymax=max)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5, fill = 'blue')

figure = ggarrange(fig1,fig2, labels = c('normal','smoothed'), ncol = 2,nrow = 1)
figure


#-----------------SMOOTHING and FPCA------------------------

temperature_smooth = smooth_data(temperature)


par(mfrow=c(1,2))
plot.fd(temperature_smooth) 


iter = 0
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    iter = iter+1
    if(st1 == 0 && st2 == 0){
      plot(isolate_scenario(temperature,st1,st2),type = 'l',col='red',lwd = 5, ylim = c(0.5,4))
      
    }else{
      lines(isolate_scenario(temperature,st1,st2),col=viridis(16)[iter])
    }
  }
}


# lo smoothing ha senso
temperature00 = isolate_scenario(temperature,20,0)

basis.1 <- create.bspline.basis(rangeval=c(1,50),nbasis=10)
smooth00 <- Data2fd(y = temperature00,argvals = 1:50,basisobj = basis.1)


par(mfrow = c(1,2))
plot.fd(smooth00, ylim = c(0.7,3.5))
lines(temperature00, col ='red')
plot.fd(temperature_smooth, ylim = c(0.7,3.5))

#la pca bisognerebbe usarla dopo montecarlo

pca_W.1 <- pca.fd(temperature_smooth,nharm=5,centerfns=TRUE)

par(mfrow=c(1,2))
plot(pca_W.1$values[1:10],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:10]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

par(mfrow=c(1,2))
plot.pca.fd(pca_W.1, nx=100, pointplot=FALSE, harm=c(1,2), expand=0, cycle=FALSE)





#-------------MONTECARLO---------------

#lanciare GAMS con montecarlo = 1

report = rgdx('rep.gdx',list(name = 'report'))

rep_clean = clean_report(report)


temperature_betavar = matrix(nrow= 50, ncol = 100)
for(i in 1:100){
  temperature_betavar[,i] = isolate_scenario(iteration_extract(variable_extract(rep_clean,'TATM'),i,100),0,0)
}





par(mfrow=c(1,1))

basis.1 <- create.bspline.basis(rangeval=c(1,50),nbasis=10)
smoothed <- Data2fd(y = temperature_betavar,argvals = 1:50,basisobj = basis.1)
plot.fd(smoothed) 



emission_betavar = matrix(nrow= 50, ncol = 100)
for(i in 1:100){
  emission_betavar[,i] = isolate_scenario(iteration_extract(variable_extract(rep_clean,'E'),i,100),0,0)
}

par(mfrow=c(1,1))
smoothed2 <- Data2fd(y = emission_betavar,argvals = 1:50,basisobj = basis.1)
plot.fd(smoothed2) 


pca_W.1 <- pca.fd(smoothed,nharm=5,centerfns=TRUE)
par(mfrow=c(1,2))
plot.pca.fd(pca_W.1, nx=100, pointplot=FALSE, harm=c(1,2), expand=0, cycle=FALSE)

pca_W.1 <- pca.fd(smoothed2,nharm=5,centerfns=TRUE)
par(mfrow=c(1,2))
plot.pca.fd(pca_W.1, nx=100, pointplot=FALSE, harm=c(1,2), expand=0, cycle=FALSE)



#------------BAND PLOT WITH MONTECARLO----------------


report = rgdx('rep.gdx',list(name = 'report'))

rep_clean = clean_report(report)

time = seq(2010,by=5,len = 50)


total_it = 100
basis = create.bspline.basis(rangeval=c(1,50), nbasis=20)
abscissa = 1:50
B = matrix(nrow = length(abscissa)*total_it, ncol = 4)
for (i in 1:total_it){
  temperature_mc = iteration_extract(variable_extract(rep_clean,'TATM'),i,total_it)

  A = banda(temperature_mc)
  A = as.data.frame(A)
  colnames(A) = c('max','min')
  A['base'] = isolate_scenario(temperature_mc,0,0)
  
  for(ii in 1:dim(A)[2]){
    temp <- smooth.basis(argvals=abscissa, y=A[,ii], fdParobj=basis)
    B[((i-1)*50+1):(i*50),ii] <- eval.fd(abscissa, temp$fd) 
  }
  B[((i-1)*50+1):(i*50),dim(A)[2]+1] = rep(i,50)
}
B = as.data.frame(B)
colnames(B) = c('max','min','base','iteration')


for(i in 1:total_it){
  assign(paste('p',i,sep = ''),
         ggplot(data=B[which(B$iteration == i),], aes(x=time, y=base, ymin=min, ymax=max)) + 
           geom_line() + 
           geom_ribbon(alpha=0.5, fill = heat.colors(total_it)[total_it-i+1]) +
           ylim(0.7,4) +
           xlab('Year') +
           ylab('TATM') +
           ggtitle(paste('eind =', 1.1+0.003*i))
  )
}


ggplots_in_env =NULL
for(i in 1:total_it){
  ggplots_in_env[i] = paste('p',i,sep = '')
}

animation::saveGIF(
  expr = {
    purrr::walk(
      ggplots_in_env,
      ~ plot(get(.))
    )
  },
  movie.name = "temp_eindvar.gif"
)


#-------------PRETTY PLOTS--------------------

as = rgdx(filename,list(name = 'as'))
TATM_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'TATM'))
TATM = rgdx(filename,list(name='TATM'))
temperature = isolate(TATM,as)



iter = 1
data = matrix(nrow = 26*50,ncol = 6)
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    data[(((iter-1)*50+1):(iter*50)),] = cbind(rep(iter,50),rep(st1,50),rep(st2,50),rep(0,50),1:50,isolate_scenario(temperature,st1,st2))
    iter = iter+1
  }
}
data[(((iter-1)*50+1):(iter*50)),] = cbind(rep(iter,50),rep(0,50),rep(0,50),rep(1,50),1:50,TATM_vanilla$val[,2])

data = as.data.frame(data)

colnames(data) = c('iter','st1','st2','vanilla','t','val')

time = seq(2010,by=5,len = 50)

highdata = data
highdata$high = ifelse((data$st1 == 0 & data$st2 == 0 & data$vanilla == 0), "3", "1")
highdata$high[which(data$st1 == 0 & data$st2 == 10)] = '2'
highdata$high[which(data$st1 == 10 & data$st2 == 0)] = '4'
highdata$high[which(data$st1 == 10 & data$st2 == 10)] = '5'
highdata$high[which(data$vanilla == 1)] = '6'

#highdata = highdata[1:1250,]

highdata %>%
  ggplot(aes(x=rep(time,26),y=val,group = iter,color=high, size=high)) +
  geom_line() +
  scale_color_manual(labels = c('Others','Tipping 1','No Tipping', 'Tipping 2', 'Both','Vanilla'),
                     values = c( "grey",viridis(4)[1],viridis(4)[2],viridis(4)[3],viridis(4)[4],"red")) +
  scale_size_manual(labels = c('Others','Tipping 1','No Tipping', 'Tipping 2', 'Both','Vanilla'),
                    values = c(0.2,1.5,1.5,1.5,1.5,1.5)) +
  ggtitle('Tipping points at step 10') +
  xlab('Year') +
  ylab('TATM') +
  labs(color = 'Legend', size = 'Legend')
  
  
  



as = rgdx(filename,list(name = 'as'))
UTIL_vanilla = rgdx('DICE_vanilla.gdx',list(name = 'Y'))
UTIL = rgdx(filename,list(name='Y'))
utility = isolate(UTIL,as)

time = seq(2010,by=5,len = 50)

iter = 1
data = matrix(nrow = 26*50,ncol = 6)
for(st1 in seq(0,40,by=10)){
  for(st2 in seq(0,40,by=10)){
    data[(((iter-1)*50+1):(iter*50)),] = cbind(rep(iter,50),rep(st1,50),rep(st2,50),rep(0,50),1:50,isolate_scenario(utility,st1,st2))
    iter = iter+1
  }
}
data[(((iter-1)*50+1):(iter*50)),] = cbind(rep(iter,50),rep(0,50),rep(0,50),rep(1,50),1:50,UTIL_vanilla$val[,2])

data = as.data.frame(data)

colnames(data) = c('iter','st1','st2','vanilla','t','val')

highdata = data
highdata$high = '0'
highdata$high[which(data$st1 == 0 & data$st2 == 0 & data$vanilla == 0)] = '1'
highdata$high[which(data$vanilla == 1)] = '2'

highdata %>%
  ggplot(aes(x=rep(time,26),y=val,group = iter,color=high, size=high)) +
  geom_line() +
  scale_color_manual(values = c( "red",'green','blue')) +
  scale_size_manual(values=c(1.5,1.5,1.5)) +
  theme(legend.position="none") +
  ggtitle('GDP') +
  ylab('Y') +
  xlab('Year') +
  geom_label(x = 2200,y= 1000,label = 'Stochastic', color = 'red', size = 4) +
  geom_label(x = 2100,y= 1500,label = 'Vanilla', color = 'blue', size = 4)

