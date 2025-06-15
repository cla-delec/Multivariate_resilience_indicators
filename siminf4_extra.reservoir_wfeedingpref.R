library(SimInf)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(hrbrthemes)
setwd("C:/Users/delec001/OneDrive - Wageningen University & Research/PhD_EWSinVBD/Multivariate_analysis/code_Yi")



################################
## WNV MODEL with import rate ##
################################

#model with increasing k over time, extra reservoir species



transitions <- c("M_s -> (delta_M*p_M*(k_s+(k_e-k_s)*t/n_steps)*(pref_Ba*Ba_i+pref_Bb*Bb_i)*M_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> M_e", "M_e -> gamma_M*M_e -> M_i",
                 "@ -> b_M*(M_s+M_e+M_i) -> M_s", "M_s -> M_s*b_M -> @",  "M_e->b_M*M_e -> @","M_i->M_i*b_M -> @",
                 "Ba_s -> (delta_M*p_M*(k_s+(k_e-k_s)*t/n_steps)*pref_Ba*M_i*Ba_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> Ba_e", "Ba_e -> gamma_B*Ba_e -> Ba_i", "Ba_i -> (1-nu_Ba)*alpha_B*Ba_i -> Ba_r", 
                 "@ -> rate_import -> Ba_i","@ -> b_B*(Ba_s+Ba_e+Ba_i+Ba_r) -> Ba_s", "Ba_s -> b_B*Ba_s -> @", "Ba_e -> b_B*Ba_e-> @", "Ba_i -> b_B*Ba_i + nu_Ba*alpha_B*Ba_i -> Ba_d", "Ba_r -> b_B*Ba_r-> @",
                 "Bb_s -> (delta_M*p_M*(k_s+(k_e-k_s)*t/n_steps)*pref_Bb*M_i*Bb_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> Bb_e", "Bb_e -> gamma_B*Bb_e -> Bb_i", "Bb_i -> (1-nu_Bb)*alpha_B*Bb_i -> Bb_r", 
                 "@ -> rate_import -> Bb_i","@ -> b_B*(Bb_s+Bb_e+Bb_i+Bb_r) -> Bb_s", "Bb_s -> b_B*Bb_s -> @", "Bb_e -> b_B*Bb_e-> @", "Bb_i -> b_B*Bb_i + nu_Bb*alpha_B*Bb_i -> @", "Bb_r -> b_B*Bb_r-> @",
                 "H_s -> (delta_M*p_M*(k_s+(k_e-k_s)*t/n_steps)*pref_H*M_i*H_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> H_e", "H_e -> gamma_H*H_e -> H_i", "H_i -> alpha_H*H_i -> H_r",
                 "E_s -> (delta_M*p_M*(k_s+(k_e-k_s)*t/n_steps)*pref_E*M_i*E_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> E_e", "E_e -> gamma_E*E_e -> E_i", "E_i -> (1-nu_E)*alpha_E*E_i -> E_r",
                 "@ -> b_E*(E_s+E_e+E_i+E_r) -> E_s", "E_s -> b_E*E_s -> @", "E_e -> b_E*E_e-> @", "E_i -> b_E*E_i + nu_E*alpha_E*E_i -> @", "E_r -> b_E*E_r-> @")




compartments <- c("M_s","M_e","M_i",
                  "Ba_s","Ba_e","Ba_i","Ba_r", "Ba_d",
                  "Bb_s","Bb_e","Bb_i","Bb_r", 
                  "H_s","H_e","H_i","H_r",
                  "E_s","E_e","E_i","E_r")

#number replicates
n <- 100

#initial conditions
N_m=300000
N_ba=200000
N_bb=150000
N_H=2000000
N_E=40000

u0 <- data.frame(M_s = rep(N_m, n), M_e = rep(0, n), M_i = rep(0, n),
                 Ba_s = rep(N_ba, n), Ba_e = rep(0, n), Ba_i = rep(1, n), Ba_r = rep(0,n), Ba_d = rep(0,n),
                 Bb_s = rep(N_bb, n), Bb_e = rep(0, n), Bb_i = rep(1, n), Bb_r = rep(0,n), 
                 H_s = rep(N_H, n), H_e = rep(0, n), H_i = rep(0, n), H_r = rep(0,n),
                 E_s = rep(N_E, n), E_e = rep(0, n), E_i = rep(0, n), E_r = rep(0,n))


#other parameters
#Temp=22 #for temp dependent parameters
b_M_val= 0.22 #(0.7998/10*(1+1.231*exp(-0.184*(Temp-20))))
gamma_B_val=0.55
gamma_H_val=0.25
gamma_E_val=0.05
gamma_M_val=0.11
nu_Ba_val=0.3
nu_Bb_val=0.05
nu_E_val=0.04
alpha_B_val=0.31
alpha_H_val=0.0714
alpha_E_val=0.2
b_B_val=0.002
b_E_val=0.00016
rate_import_val=0.15
delta_M_val=1
pref_Ba_val=5
pref_Bb_val=10
pref_H_val=1
pref_E_val=1
p_M_val=0.9
k_s_val=0.66
k_e_val=0.94
#k_s_val=0.75
#k_e_val=0.75

k_s<-seq(0.1,2,length.out=100)

# calculate R0

calc_R0<-function(k_val) {
  big_D=N_ba*pref_Ba_val+N_bb*pref_Bb_val+N_H*pref_H_val+N_E*pref_E_val #denominator: sum of hosts * corresponding feeding pref
  
  Tr<-matrix(0, 10, 10) #transmission matrix
  Tr[1,4]=N_m*delta_M_val*p_M_val*k_val*pref_Ba_val/big_D #transmission to mos by bird A
  Tr[1,6]=N_m*delta_M_val*p_M_val*k_val*pref_Bb_val/big_D #transmission to mos by bird B
  Tr[3,2]=N_ba*delta_M_val*p_M_val*k_val*pref_Ba_val/big_D #transmission to bird A by mos
  Tr[5,2]=N_bb*delta_M_val*p_M_val*k_val*pref_Bb_val/big_D #transmission to bird B by mos
  Tr[7,2]=N_H*delta_M_val*p_M_val*k_val*pref_H_val/big_D #transmission to human by mos
  Tr[9,2]=N_E*delta_M_val*p_M_val*k_val*pref_E_val/big_D #transmission to horse by mos
  
  Sigma<-matrix(0, 10, 10) #transition matrix
  Sigma[1,1]=-(gamma_M_val+b_M_val)
  Sigma[2,1]=gamma_M_val
  Sigma[2,2]=-b_M_val
  Sigma[3,3]=-(gamma_B_val+b_B_val) 
  Sigma[4,3]=gamma_B_val
  Sigma[4,4]=-(b_B_val+alpha_B_val)
  Sigma[5,5]=-(gamma_B_val+b_B_val)
  Sigma[6,5]=gamma_B_val
  Sigma[6,6]=-(b_B_val+alpha_B_val)
  Sigma[7,7]=-gamma_H_val
  Sigma[8,7]=gamma_H_val
  Sigma[8,8]=-alpha_H_val
  Sigma[9,9]=-(gamma_E_val+b_E_val)
  Sigma[10,9]=gamma_E_val
  Sigma[10,10]=-(alpha_E_val+b_E_val)
  
  K_L<- -Tr %*% solve(Sigma) #NGM w large domain
  eig<-eigen(K_L)
  R_0<-max(abs(eig$values))
  return(R_0)
  
  # E<-matrix(0, 10, 5)
  # E[1,1]<-1
  # E[3,2]<-1
  # E[5,3]<-1
  # E[7,4]<-1
  # E[9,5]<-1
  # 
  # K<-t(E) %*% K_L %*% E #(real) NGM 
  # eig2<-eigen(K)
}

R0_vals<-lapply(k_s, calc_R0)
plot(k_s,R0_vals)
abline(h=0.7,col='red')
abline(h=1,col='red')

id_07<-min(which(R0_vals>=0.7))
k_s[id_07]

id_1<-min(which(R0_vals>1))
k_s[id_1]

id_08<-min(which(R0_vals>=0.8))
k_s[id_08]

num_steps=3000;

model <- mparse(transitions = transitions, compartments = compartments,
                gdata = c(gamma_M=gamma_M_val, gamma_B=gamma_B_val, gamma_H=gamma_H_val, gamma_E=gamma_E_val,
                          nu_Ba=nu_Ba_val, nu_Bb=nu_Bb_val, nu_E=nu_E_val,
                          alpha_B=alpha_B_val, alpha_H=alpha_H_val, alpha_E=alpha_E_val,
                          b_M=b_M_val, b_B=b_B_val, b_E=b_E_val, 
                          rate_import=rate_import_val,
                          delta_M=delta_M_val, k_s=k_s_val, k_e=k_e_val,
                          pref_Ba=pref_Ba_val, pref_Bb=pref_Bb_val, pref_H=pref_H_val, pref_E=pref_E_val,
                          p_M=p_M_val, n_steps=num_steps),
                u0 = u0, tspan = 1:num_steps)



set_num_threads(123)
result <- run(model = model)


# PLOTS -------------------------------------------------------------------

plot(result@U[3,],type="l", col='blue', xlab="time", ylab="Incidence")
lines(result@U[6,], col='red')
new_deads<-diff(result@U[8,])
lines(c(0,new_deads), col='pink')
lines(result@U[11,], col='purple')
lines(result@U[15,], col='green')
lines(result@U[19,], col='black')
legend("topleft",legend=c("Infected mosquitoes","Infected birds 1","Dead birds 1","Infected birds 2","Infected humans","Infected horses"),col=c('blue','red','pink','purple','green','black'),lty=rep(1,4))

plot(result@U[23,],type="l", col='blue', xlab="time", ylab="Incidence")
lines(result@U[26,], col='red')
new_deads<-diff(result@U[8,])
lines(c(0,new_deads), col='pink')
lines(result@U[31,], col='purple')
lines(result@U[35,], col='green')
lines(result@U[39,], col='black')
legend("topleft",legend=c("Infected mosquitoes","Infected birds 1","Dead birds 1","Infected birds 2","Infected humans","Infected horses"),col=c('blue','red','pink','purple','green','black'),lty=rep(1,4))

plot(result@U[43,],type="l", col='blue', xlab="time", ylab="Incidence")
lines(result@U[46,], col='red')
new_deads<-diff(result@U[8,])
lines(c(0,new_deads), col='pink')
lines(result@U[51,], col='purple')
lines(result@U[55,], col='green')
lines(result@U[59,], col='black')
legend("topleft",legend=c("Infected mosquitoes","Infected birds 1","Dead birds 1","Infected birds 2","Infected humans","Infected horses"),col=c('blue','red','pink','purple','green','black'),lty=rep(1,4))


plot(result@U[63,],type="l", col='blue', xlab="time", ylab="Incidence")
lines(result@U[66,], col='red')
new_deads<-diff(result@U[8,])
lines(c(0,new_deads), col='pink')
lines(result@U[71,], col='purple')
lines(result@U[75,], col='green')
lines(result@U[79,], col='black')
legend("topleft",legend=c("Infected mosquitoes","Infected birds 1","Dead birds 1","Infected birds 2","Infected humans","Infected horses"),col=c('blue','red','pink','purple','green','black'),lty=rep(1,4))


# SAVE DATA -------------------------------------------------------------------


data_final=tail(result@U,10)
data_final=result@U[c(18, 21, 22, 25, 26, 29, 30),];

n_ts=dim(result@U)[1]
indexes=c(seq(3,n_ts,20), seq(6,n_ts,20), seq(7,n_ts,20), seq(8,n_ts,20), seq(11,n_ts,20), seq(12,n_ts,20), seq(15,n_ts,20), seq(16,n_ts,20), seq(19,n_ts,20), seq(20,n_ts,20));
indexes=sort(indexes)

data_final=result@U[indexes,]

#data_final[seq(4,n*10,10),]=rbind(rep(0,n),apply(data_final[seq(4,n*10,10),], 1, diff))
#matplot(t(data_final[seq(4,n*10,10),]), type='l')

write.csv(data_final, "data/siminf_extrares_feedingprefB10_propBA0.75_ts3000_100rep_k0.66-0.94_R00.7-1_rateimport0.15.csv", row.names=FALSE)
#write.csv(data_final, "data/siminf_extrares_feedingprefB10_propBA0.75_ts3000_100rep_k0.75_R00.8_rateimport0.15.csv", row.names=FALSE)



