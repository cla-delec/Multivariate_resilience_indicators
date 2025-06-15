library(SimInf)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(hrbrthemes)
#setwd("C:/Users/delec001/OneDrive - Wageningen University & Research/PhD_EWSinVBD/Multivariate_analysis/code_Yi")

rm(list=ls())

# Initialize model and simu params ----------------------------------------

transitions <- c("M_s -> (delta_M*p_M*k_fixed*(pref_Ba*Ba_i+pref_Bb*Bb_i)*M_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> M_e", "M_e -> gamma_M*M_e -> M_i",
                 "@ -> b_M*(M_s+M_e+M_i) -> M_s", "M_s -> M_s*b_M -> @",  "M_e->b_M*M_e -> @","M_i->M_i*b_M -> @",
                 "Ba_s -> (delta_M*p_M*k_fixed*pref_Ba*M_i*Ba_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> Ba_e", "Ba_e -> gamma_B*Ba_e -> Ba_i", "Ba_i -> (1-nu_Ba)*alpha_B*Ba_i -> Ba_r", 
                 "@ -> rate_import -> Ba_i","@ -> b_B*(Ba_s+Ba_e+Ba_i+Ba_r) -> Ba_s", "Ba_s -> b_B*Ba_s -> @", "Ba_e -> b_B*Ba_e-> @", "Ba_i -> b_B*Ba_i + nu_Ba*alpha_B*Ba_i -> Ba_d", "Ba_r -> b_B*Ba_r-> @",
                 "Bb_s -> (delta_M*p_M*k_fixed*pref_Bb*M_i*Bb_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> Bb_e", "Bb_e -> gamma_B*Bb_e -> Bb_i", "Bb_i -> (1-nu_Bb)*alpha_B*Bb_i -> Bb_r", 
                 "@ -> rate_import -> Bb_i","@ -> b_B*(Bb_s+Bb_e+Bb_i+Bb_r) -> Bb_s", "Bb_s -> b_B*Bb_s -> @", "Bb_e -> b_B*Bb_e-> @", "Bb_i -> b_B*Bb_i + nu_Bb*alpha_B*Bb_i -> @", "Bb_r -> b_B*Bb_r-> @",
                 "H_s -> (delta_M*p_M*k_fixed*pref_H*M_i*H_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> H_e", "H_e -> gamma_H*H_e -> H_i", "H_i -> alpha_H*H_i -> H_r",
                 "E_s -> (delta_M*p_M*k_fixed*pref_E*M_i*E_s)/(pref_Ba*(Ba_s+Ba_e+Ba_i+Ba_r)+pref_Bb*(Bb_s+Bb_e+Bb_i+Bb_r)+pref_H*(H_s+H_e+H_i+H_r)+pref_E*(E_s+E_e+E_i+E_r)) -> E_e", "E_e -> gamma_E*E_e -> E_i", "E_i -> (1-nu_E)*alpha_E*E_i -> E_r",
                 "@ -> b_E*(E_s+E_e+E_i+E_r) -> E_s", "E_s -> b_E*E_s -> @", "E_e -> b_E*E_e-> @", "E_i -> b_E*E_i + nu_E*alpha_E*E_i -> @", "E_r -> b_E*E_r-> @")




compartments <- c("M_s","M_e","M_i",
                  "Ba_s","Ba_e","Ba_i","Ba_r", "Ba_d",
                  "Bb_s","Bb_e","Bb_i","Bb_r", 
                  "H_s","H_e","H_i","H_r",
                  "E_s","E_e","E_i","E_r")

#number replicates
n <- 10

#initial conditions
N_m=300000
N_ba=200000
N_bb=150000
N_H=2000000
N_E=40000

u0 <- data.frame(M_s = rep(N_m, n), M_e = rep(0, n), M_i = rep(0, n),
                 Ba_s = rep(N_ba, n), Ba_e = rep(0, n), Ba_i = rep(0, n), Ba_r = rep(0,n), Ba_d = rep(0,n),
                 Bb_s = rep(N_bb, n), Bb_e = rep(0, n), Bb_i = rep(0, n), Bb_r = rep(0,n), 
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
rate_import_val=0
delta_M_val=1
pref_Ba_val=5
pref_Bb_val=10
pref_H_val=1
pref_E_val=1
p_M_val=0.9
k_s_val=0.66
k_e_val=0.94

num_steps=20;


# R0 values ---------------------------------------------------------------

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
  
}

R0_vals<-lapply(k_s, calc_R0)
plot(k_s,R0_vals)
abline(h=0.7,col='red')
abline(h=1,col='red')

id_03<-min(which(R0_vals>=0.3))
k_s[id_03]

id_08<-min(which(R0_vals>=0.8))
k_s[id_08]

id_02<-min(which(R0_vals>=0.2))
k_s[id_02]

id_1<-min(which(R0_vals>=1))
k_s[id_1]

# Example figure ----------------------------------------------------------

#add event at t=5
E<-matrix(0,20,1, dimnames = list(c("M_s","M_e","M_i",
                                    "Ba_s","Ba_e","Ba_i","Ba_r", "Ba_d",
                                    "Bb_s","Bb_e","Bb_i","Bb_r", 
                                    "H_s","H_e","H_i","H_r",
                                    "E_s","E_e","E_i","E_r"), c("1")))
E[11,1]<-1

size_perturb=5
add <- data.frame(event = "enter", time = 5,
                  node = 1, dest = 0, n = size_perturb, proportion = 0, select = 1, shift = 0)

#low R0
model_lowr0 <- mparse(transitions = transitions, compartments = compartments,
                      events = add, E=E,
                               gdata = c(gamma_M=gamma_M_val, gamma_B=gamma_B_val, gamma_H=gamma_H_val, gamma_E=gamma_E_val,
                                         nu_Ba=nu_Ba_val, nu_Bb=nu_Bb_val, nu_E=nu_E_val,
                                         alpha_B=alpha_B_val, alpha_H=alpha_H_val, alpha_E=alpha_E_val,
                                         b_M=b_M_val, b_B=b_B_val, b_E=b_E_val, 
                                         rate_import=rate_import_val,
                                         delta_M=delta_M_val, k_fixed=k_s[id_03],
                                         pref_Ba=pref_Ba_val, pref_Bb=pref_Bb_val, pref_H=pref_H_val, pref_E=pref_E_val,
                                         p_M=p_M_val, n_steps=num_steps),
                               u0 = u0, tspan = 1:num_steps)
set_num_threads(123)
result_lowr0 <- run(model = model_lowr0)
plot(result_lowr0@U[11,1:20], type="l", col='red', xlab="time", ylab="Prevalence")
#return rate
df<-data.frame(time=0:7,inc=log(c(size_perturb,result_lowr0@U[11,6:12]+0.01)))
# lm<-lm(inc~time,data=df)
# ret_rate_lowr0<-abs(lm$coefficients[2])
# intercept_lowr0<-abs(lm$coefficients[1])

lm2<-lm(inc~0+time, offset=rep(log(size_perturb),length(df$inc)), data=df)
plot(df$time, df$inc)
ret_rate_lowr0<-abs(lm2$coefficients[1])
intercept_lowr0<-log(size_perturb)

#high R0
model_highr0 <- mparse(transitions = transitions, compartments = compartments,
                       events = add, E=E,
                       gdata = c(gamma_M=gamma_M_val, gamma_B=gamma_B_val, gamma_H=gamma_H_val, gamma_E=gamma_E_val,
                                 nu_Ba=nu_Ba_val, nu_Bb=nu_Bb_val, nu_E=nu_E_val,
                                 alpha_B=alpha_B_val, alpha_H=alpha_H_val, alpha_E=alpha_E_val,
                                 b_M=b_M_val, b_B=b_B_val, b_E=b_E_val, 
                                 rate_import=rate_import_val,
                                 delta_M=delta_M_val, k_fixed=k_s[id_08],
                                 pref_Ba=pref_Ba_val, pref_Bb=pref_Bb_val, pref_H=pref_H_val, pref_E=pref_E_val,
                                 p_M=p_M_val, n_steps=num_steps),
                       u0 = u0, tspan = 1:num_steps)

result_highr0 <- run(model = model_highr0)
plot(result_highr0@U[11,1:20], type="l", col='red', xlab="time", ylab="Prevalence")
#return rate
df<-data.frame(time=0:15,inc=log(c(size_perturb,result_highr0@U[11,6:20])+0.01))
# lm<-lm(inc~time,data=df)
# ret_rate_highr0<-abs(lm$coefficients[2])
# intercept_highr0<-abs(lm$coefficients[1])
lm2<-lm(inc~0+time, offset=rep(log(size_perturb),length(df$inc)), data=df)
plot(df$time, df$inc)
abline(lm2,col="red")
ret_rate_highr0<-abs(lm2$coefficients[1])
intercept_highr0<-log(size_perturb)

#plot all
pdf("Figures/perturb_recov_examples4_10Bi.pdf") 
plot(0:20,c(0,result_lowr0@U[11,1:20]), xlab="Time (days)", ylab="Number of infected birds",ylim=c(0,10),pch=15,bty="n", cex=1.5,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col='#D04848') #LOW
lines(0:20,c(rep(NA,5),exp(-(0:15)*ret_rate_lowr0+intercept_lowr0)),lty=2, col='#D04848',lwd = 1) #FIT LOw
points((0:20)+0.1,c(0,result_highr0@U[11,1:20]),pch=17, col='#6895D2', cex=1.5) #HIGH
lines(0:20,c(rep(NA,5),exp(-(0:15)*ret_rate_highr0+intercept_highr0)),lty=2, col='#6895D2') #FIT HIGH
abline(v=5,lw=2)
legend_text <-c(as.expression(bquote('R'['0']*'=0.3, high resilience')), paste('Fitted line R0=0.3, return rate=',as.character(round(ret_rate_lowr0,2))), as.expression(bquote('R'['0']*'=0.8, low resilience')), paste('Fitted line R0=0.8, return rate=',as.character(round(ret_rate_highr0,2))))
legend("topright", legend_text,pch=c(15,NA,17,NA),lty=c(NA,2,NA,2), col=c('#D04848','#D04848','#6895D2','#6895D2'),cex=1.5)
dev.off() 


# Perturbation - recovery - with log method -birds -------------------------------

#number replicates
n <- 100

#num steps simulation
num_steps=100

#size perturbation
size_perturb<-10

#initial conditions
u0 <- data.frame(M_s = rep(N_m, n), M_e = rep(0, n), M_i = rep(0, n),
                 Ba_s = rep(N_ba, n), Ba_e = rep(0, n), Ba_i = rep(0, n), Ba_r = rep(0,n), Ba_d = rep(0,n),
                 Bb_s = rep(N_bb, n), Bb_e = rep(0, n), Bb_i = rep(size_perturb, n), Bb_r = rep(0,n), 
                 H_s = rep(N_H, n), H_e = rep(0, n), H_i = rep(0, n), H_r = rep(0,n),
                 E_s = rep(N_E, n), E_e = rep(0, n), E_i = rep(0, n), E_r = rep(0,n))

#Perturbation of the bird state
n_val_tested=25
k_vals<-seq(from = k_s[id_02], to = k_s[id_1], length.out = n_val_tested)
res_b_bperturb = matrix(, nrow = length(k_vals), ncol = n)
indexes_bi<-seq(11,n*20,20)

max_n_point<-100

#function to find 2 consecutive zeros
find_first_consecutive_zeros <- function(vec) {
  # Logical vector checking for two consecutive zeros
  consec_zeros <- (vec[-length(vec)] == 0) & (vec[-1] == 0)
  
  # Get the index of the first TRUE value (for the first zero in the pair)
  first_index <- which(consec_zeros)[1]
  
  # Return the index or NA if no consecutive zeros found
  if (is.na(first_index)) {
    return(NA)
  } else {
    return(first_index)
  }
}


#for loop for all results
for (i in 1:length(k_vals)) {
  cur_k=k_vals[i]
  
  model <- mparse(transitions = transitions, compartments = compartments,
                  gdata = c(gamma_M=gamma_M_val, gamma_B=gamma_B_val, gamma_H=gamma_H_val, gamma_E=gamma_E_val,
                            nu_Ba=nu_Ba_val, nu_Bb=nu_Bb_val, nu_E=nu_E_val,
                            alpha_B=alpha_B_val, alpha_H=alpha_H_val, alpha_E=alpha_E_val,
                            b_M=b_M_val, b_B=b_B_val, b_E=b_E_val, 
                            rate_import=rate_import_val,
                            delta_M=delta_M_val, k_fixed=cur_k,
                            pref_Ba=pref_Ba_val, pref_Bb=pref_Bb_val, pref_H=pref_H_val, pref_E=pref_E_val,
                            p_M=p_M_val, n_steps=num_steps),
                  u0 = u0, tspan = 1:num_steps)
  result <- run(model = model)
  
  for (j in 1:n) {
    cur_ts<-c(size_perturb,result@U[indexes_bi[j],])
    
    cur_max_point<-ifelse(is.na(find_first_consecutive_zeros(cur_ts)),max_n_point,find_first_consecutive_zeros(cur_ts))
    
    df<-data.frame(time=0:cur_max_point,inc=log(cur_ts[1:(cur_max_point+1)]+0.01))
    lm<-lm(inc~0+time, offset=rep(log(size_perturb),length(df$inc)), data=df)
    res_b_bperturb[i,j]<-abs(lm$coefficients[1])

    # lm<-lm(inc~time,data=df)
    # ret_rate_lowr0<-abs(lm$coefficients[2])
    # intercept_lowr0<-abs(lm$coefficients[1])
  }
  
}


#saveRDS(res_b_bperturb, file="data/res_b_bperturb_10Bi_100000reps_v2.Rda")

#res_b_bperturb<-readRDS("data/res_b_bperturb_10Bi_100000reps_v2.Rda")
res_b_bperturb<-readRDS("data/res_b_bperturb_10Bi_1000reps_25vals_max100_v3.Rda")

recov_per_k_b <- rowMeans(res_b_bperturb)
recov_per_k_b_lower <- apply(res_b_bperturb, 1, quantile, probs = 0.05)
recov_per_k_b_upper <- apply(res_b_bperturb, 1, quantile, probs = 0.95)


plot(k_vals,recov_per_k_b, type="l",col="blue")
plot(k_vals,recov_per_k_b_lower, lty=2)
plot(k_vals,recov_per_k_b_upper, lty=2)

plot(k_vals,1/recov_per_k_b, type="l",col="blue", ylim=c(0,8))
lines(k_vals,1/recov_per_k_b_lower, lty=2)
lines(k_vals,1/recov_per_k_b_upper, lty=2)

R0_vals<-seq(from = 0.2, to = 1, length.out = n_val_tested)

#pdf("Figures/plot_recov_time_birds_10Bi_1000reps_3.pdf") 
plot(R0_vals,1/recov_per_k_b, type="l", ylab="Recovery time", main="Recovery time after introducing infected birds", xlab="R0",bty="n", lwd=2,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,8))
lines(k_vals,1/recov_per_k_b_lower, lty=2)
lines(k_vals,1/recov_per_k_b_upper, lty=2)
dev.off() 

#plot 2
new_df_recov<-as.data.frame(cbind(R0_vals,res_b_bperturb))
colnames(new_df_recov)<-c("R0_vals",paste(rep('rep',n), as.character(1:n), sep = ""))
new_df_recov2<- new_df_recov %>% pivot_longer(cols= -R0_vals) %>% mutate(recov=1/value)

ggplot(new_df_recov2, aes(x=R0_vals, y=recov))+
  stat_summary(geom="ribbon", fun.data=mean_cl_normal,fun.args=list(conf.int=0.95), fill="lightblue")+
  stat_summary(geom="line", fun.y=mean, linetype="dashed")


# Perturbation - recovery - with log method -mosquitoes -------------------------------

#number replicates
n <- 1000

#size perturbation
size_perturb<-10

#initial conditions
u0 <- data.frame(M_s = rep(N_m, n), M_e = rep(0, n), M_i = rep(size_perturb, n),
                 Ba_s = rep(N_ba, n), Ba_e = rep(0, n), Ba_i = rep(0, n), Ba_r = rep(0,n), Ba_d = rep(0,n),
                 Bb_s = rep(N_bb, n), Bb_e = rep(0, n), Bb_i = rep(0, n), Bb_r = rep(0,n), 
                 H_s = rep(N_H, n), H_e = rep(0, n), H_i = rep(0, n), H_r = rep(0,n),
                 E_s = rep(N_E, n), E_e = rep(0, n), E_i = rep(0, n), E_r = rep(0,n))

#Perturbation of the mosquito state
res_m_mperturb = matrix(, nrow = length(k_vals), ncol = n)
indexes_mi<-seq(3,n*20,20)

max_n_point<-50

#for loop for all results
for (i in 1:length(k_vals)) {
  cur_k=k_vals[i]
  
  model <- mparse(transitions = transitions, compartments = compartments,
                  gdata = c(gamma_M=gamma_M_val, gamma_B=gamma_B_val, gamma_H=gamma_H_val, gamma_E=gamma_E_val,
                            nu_Ba=nu_Ba_val, nu_Bb=nu_Bb_val, nu_E=nu_E_val,
                            alpha_B=alpha_B_val, alpha_H=alpha_H_val, alpha_E=alpha_E_val,
                            b_M=b_M_val, b_B=b_B_val, b_E=b_E_val, 
                            rate_import=rate_import_val,
                            delta_M=delta_M_val, k_fixed=cur_k,
                            pref_Ba=pref_Ba_val, pref_Bb=pref_Bb_val, pref_H=pref_H_val, pref_E=pref_E_val,
                            p_M=p_M_val, n_steps=num_steps),
                  u0 = u0, tspan = 1:num_steps)
  result <- run(model = model)
  
  for (j in 1:n) {
    cur_ts<-c(size_perturb,result@U[indexes_mi[j],])
    
    cur_max_point<-ifelse(is.na(find_first_consecutive_zeros(cur_ts)),max_n_point,find_first_consecutive_zeros(cur_ts))
    
    df<-data.frame(time=0:cur_max_point,inc=log(cur_ts[1:(cur_max_point+1)]+0.01))
    lm<-lm(inc~0+time, offset=rep(log(size_perturb),length(df$inc)), data=df)
    res_m_mperturb[i,j]<-abs(lm$coefficients[1])
    
  }
  
}


#saveRDS(res_m_mperturb, file="data/res_m_mperturb_10Mi_1000reps_v3.Rda")

res_m_mperturb<-readRDS(file="data/res_m_mperturb_10Mi_1000reps_v2.Rda")

recov_per_k_m<-rowMeans(res_m_mperturb)
recov_per_k_m<-apply(res_m_mperturb, 1, median)

res_recov_m<-1/res_m_mperturb
recov_per_k_m_lower <- apply(res_m_mperturb, 1, quantile, probs = 0.05)
recov_per_k_m_upper <- apply(res_m_mperturb, 1, quantile, probs = 0.95)


plot(k_vals,recov_per_k_m, type="l",col="blue")

plot(k_vals,1/recov_per_k_m, type="l",col="blue", , ylim=c(1,15))
lines(k_vals,1/recov_per_k_m_lower, lty=2)
lines(k_vals,1/recov_per_k_m_upper, lty=2)

pdf("Figures/plot_recov_time_mosquitoes_10Mi_1000reps_3.pdf") 
#R0_vals<-c(0.19, 0.28, 0.37, 0.46, 0.56, 0.65, 0.74, 0.84,0.93, 1.02)
plot(R0_vals,1/recov_per_k_m, type="l", 
     ylab="Recovery time", main="Recovery time after introducing infected mosquitoes", 
     xlab="R0",bty="n", lwd=2,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
     cex.sub=1.5, ylim=c(1,15))
lines(k_vals,1/recov_per_k_m_lower, lty=2)
lines(k_vals,1/recov_per_k_m_upper, lty=2)
dev.off() 

