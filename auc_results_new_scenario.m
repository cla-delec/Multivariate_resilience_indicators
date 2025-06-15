%% load data

res_tot=load("data/results_roc_curve_siminf_extrares_feedingprefB10_length3000_100rep_tp_k0.65-0.92_R00.7-1_tn_k0.75_R00.8_rateimport0.15_agg.mat");

res_tp_multi=res_tot.final_results.all_taus_tp_inc_multi;
res_tn_multi=res_tot.final_results.all_taus_tn_stable_multi;
res_uni_tp=res_tot.final_results.all_taus_tp_inc_uni;
res_uni_tn=res_tot.final_results.all_taus_tn_stable_uni;


%% calculate systematically for each multi indicator and multi scenario

number_dim=10;
number_rep=size(res_tp_multi,3);

scenarios_multi=[1 2 7 9; 3 8 10 nan; 7 8 9 10; 1 2 3 4; 1 5 6 nan]; %scenarios for multivariate indicators
id_incidence=[1 2 4 5 7 9]; %also dead because are aggregated by summing
id_sero=[3 6 8 10];
tot_n_scenarios_multi=size(scenarios_multi,1);

label=[{'epidemics'},{'no epidemics'}];
label_tot=repelem(label,number_rep);

indicators_multi=["Deg Finger", "MAF AC", "MAF eig", "Mean AR","Max AR","MAF var","Max var","Mean var","PCA var","Expl var","Max abscorr","Max cov"];
n_indic_multi=size(indicators_multi,2);
id_autocor_based=[1,2,3,4,5]; n_indic_autocor=size(id_autocor_based,2);
id_var_based=[6,7,8,9,10,12]; n_indic_var=size(id_var_based,2);

name_scenarios=["Incidence scenario", "Seroprevalence scenario", "Anthro-equine scenario", "Wildlife scenario", "Hiding reservoir scenario"];

AUC_s=zeros(tot_n_scenarios_multi,n_indic_multi);
X_s=zeros(number_rep*2+1,tot_n_scenarios_multi,n_indic_multi);
Y_s=zeros(number_rep*2+1,tot_n_scenarios_multi,n_indic_multi);

for cur_scenario=1:tot_n_scenarios_multi
    subplot(3,2,cur_scenario)
    hold on
    title(name_scenarios(cur_scenario),'FontSize',14)
    xlabel('False positive rate','FontSize',12); ylabel('True positive rate','FontSize',12);
    for cur_indic=1:n_indic_multi
        res_tp_scen_indic=squeeze(res_tp_multi(cur_indic,cur_scenario,:));
        res_tn_scen_indic=squeeze(res_tn_multi(cur_indic,cur_scenario,:));
        [X_cur,Y_cur,T,AUC_s(cur_scenario,cur_indic)]=perfcurve(label_tot,[res_tp_scen_indic.',res_tn_scen_indic.'],'epidemics','ProcessNaN','addtofalse'); 
        X_s(1:size(X_cur,1),cur_scenario,cur_indic)=X_cur;
        Y_s(1:size(Y_cur,1),cur_scenario,cur_indic)=Y_cur;
        
        plot(X_cur,Y_cur)
    end
    legend(indicators_multi,'FontSize',12)
    hold off
end


%% calculate systematically for uni indicators and scenario

number_dim=10;
number_rep=size(res_tp_multi,3);

scenarios_uni=1:number_dim; %scenarios for multivariate indicators

tot_n_scenarios_uni=size(scenarios_uni,2);

indicators_uni=["AR", "std"];
n_indic_uni=size(indicators_uni,2);
name_scenarios_uni=["M_I", "Ba_I", "Ba_R", "Ba_D", "Bb_I", "Bb_R", "H_I", "H_R", "E_I", "E_R"];
name_scenarios_uni2=["Mosquito incidence", "Bird A incidence", "Bird A seroprevalence", "Bird A dead", "Bird B incidence", "Bird B seroprevalence", "Human incidence", "Human seroprevalence", "Horse incidence", "Horse seroprevalence"];

AUC_s_uni=zeros(tot_n_scenarios_uni,n_indic_uni);
X_s_uni=zeros(number_rep*2+1,tot_n_scenarios_uni,n_indic_uni);
Y_s_uni=zeros(number_rep*2+1,tot_n_scenarios_uni,n_indic_uni);

for cur_scenario=1:tot_n_scenarios_uni
    subplot(2,5,cur_scenario)
    hold on
    title(name_scenarios_uni(cur_scenario),'FontSize',14)
    xlabel('False positive rate','FontSize',12); ylabel('True positive rate','FontSize',12);
    for cur_indic=1:n_indic_uni
        res_tp_scen_indic_uni=squeeze(res_uni_tp(cur_indic,cur_scenario,:));
        res_tn_scen_indic_uni=squeeze(res_uni_tn(cur_indic,cur_scenario,:));
        [X_cur,Y_cur,T,AUC_s_uni(cur_scenario,cur_indic)]=perfcurve(label_tot,[res_tp_scen_indic_uni.',res_tn_scen_indic_uni.'],'epidemics','ProcessNaN','addtofalse'); 
        X_s_uni(1:size(X_cur,1),cur_scenario,cur_indic)=X_cur;
        Y_s_uni(1:size(Y_cur,1),cur_scenario,cur_indic)=Y_cur;
        
        plot(X_cur,Y_cur)
    end
    legend(indicators_uni,'FontSize',12)
    hold off
end

%% diff plot
% 
% AUC_s=abs(AUC_s-0.5);
% AUC_s_uni=abs(AUC_s_uni-0.5);

%% barplots for all scenarios

%sort by increasing order
avg_AUCs_allsc=mean(AUC_s,1);
[sorted_AUC,idx] = sort(avg_AUCs_allsc);

subplot(5,1,1)
bar(AUC_s(1,idx))
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('AUC'); title(strcat("Performance of the different indicators for the ", name_scenarios(1)))
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,1,2)
bar(AUC_s(2,idx))
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('AUC'); title(strcat("Performance of the different indicators for the ", name_scenarios(2)))
ax = gca; 
ax.FontSize = 10; 
box off


subplot(5,1,3)
bar(AUC_s(3,idx))
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('AUC'); title(strcat("Performance of the different indicators for the ", name_scenarios(3)))
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,1,4)
bar(AUC_s(4,idx))
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('AUC'); title(strcat("Performance of the different indicators for the ", name_scenarios(4)))
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,1,5)
bar(AUC_s(5,idx))
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('AUC'); title(strcat("Performance of the different indicators for the ", name_scenarios(5)))
ax = gca; 
ax.FontSize = 10; 
box off

%% barplots for all scenarios + uni

%sort by increasing order
avg_AUCs_allsc=mean(AUC_s,1);
[sorted_AUC,idx] = sort(avg_AUCs_allsc);

avg_AUCs_allsc_uni=mean(AUC_s_uni,1);
[sorted_AUC_uni,idx_uni] = sort(avg_AUCs_allsc_uni);

subplot(5,2,1)
bar(AUC_s(1,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(strcat("Performance of the different indicators for the ", name_scenarios(1)))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,2)
b=bar(AUC_s_uni(scenarios_multi(1,:),:))
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(1,:)))
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); %title("Performance of the univariate indicators")
%legend(["Autocorrelation","Variance"],'Location','eastoutside')
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,3)
bar(AUC_s(2,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(strcat("Performance of the different indicators for the ", name_scenarios(2)))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,4)
b=bar(AUC_s_uni(scenarios_multi(2,1:3),:))
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(2,1:3)))
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); %title("Performance of the univariate indicators")
%legend(["Autocorrelation","Variance"],'Location','eastoutside')
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,5)
bar(AUC_s(3,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(strcat("Performance of the different indicators for the ", name_scenarios(3)))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,6)
b=bar(AUC_s_uni(scenarios_multi(3,:),:))
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(3,:)))
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); %title("Performance of the univariate indicators")
%legend(["Autocorrelation","Variance"],'Location','eastoutside')
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,7)
bar(AUC_s(4,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(strcat("Performance of the different indicators for the ", name_scenarios(4)))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,8)
b=bar(AUC_s_uni(scenarios_multi(4,1:4),:))
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(4,1:4)))
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); %title("Performance of the univariate indicators")
%legend(["Autocorrelation","Variance"],'Location','eastoutside')
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,9)
bar(AUC_s(5,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(strcat("Performance of the different indicators for the ", name_scenarios(5)))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,10)
b=bar(AUC_s_uni(scenarios_multi(5,1:3),:))
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(5,1:3)))
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); %title("Performance of the univariate indicators")
%legend(["Autocorrelation","Variance"],'Location','eastoutside')
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

%% new figure

%sort by increasing order
avg_AUCs_allsc=mean(AUC_s,1);
[sorted_AUC,idx] = sort(avg_AUCs_allsc);

avg_AUCs_allsc_uni=mean(AUC_s_uni,1);
[sorted_AUC_uni,idx_uni] = sort(avg_AUCs_allsc_uni);

%subplot(1,4,1)
b=bar(AUC_s_uni(:,:))
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni)
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); title("Performance of the univariate indicators")
%legend(["Autocorrelation","Variance"],'Location','eastoutside')
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off
f = gcf;
%saveas(f,'Figures/barchart_uni.pdf')


%subplot(1,4,2)
bar(AUC_s(3,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(name_scenarios(3))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off
f = gcf;
%saveas(f,'Figures/barchart_multi1.pdf')

%subplot(1,4,3)
bar(AUC_s(4,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(name_scenarios(4))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off
f = gcf;
%saveas(f,'Figures/barchart_multi2.pdf')


%subplot(1,4,4)
bar(AUC_s(5,idx), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('|AUC-0.5|'); title(name_scenarios(5))
ylim([0.48 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off
f = gcf;
%saveas(f,'Figures/barchart_multi3.pdf')


%% summary barplot - mean per scenario

%sort by increasing order
avg_AUCs_allindic=mean(AUC_s,2);
[sorted_AUC2,idx2] = sort(avg_AUCs_allindic);

avg_AUCs_uni_allindic=mean(AUC_s_uni,2);
[sorted_AUC_uni2,idx_uni2] = sort(avg_AUCs_uni_allindic);

subplot(1,2,1)
bar(avg_AUCs_allindic(idx2), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',name_scenarios(idx2))
xlabel('Monitoring scenario'); ylabel('|AUC-0.5|'); title("Multivariate indicator")
ylim([0.49 1])
yticks(linspace(0,0.5,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 15; 
box off

subplot(1,2,2)
bar(avg_AUCs_uni_allindic(idx_uni2), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',name_scenarios_uni2(idx_uni2))
xlabel('Univariate time series'); ylabel('|AUC-0.5|'); title("Univariate indicator")
ylim([0.49 1])
yticks(linspace(0,0.5,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 15; 
box off

%% summary barplot - best indic (on avg)

id_autocor_based=[1,2,3,4,5]; n_indic_autocor=size(id_autocor_based,2);
id_var_based=[6,7,8,9,10,12]; n_indic_var=size(id_var_based,2);

%sort by increasing order
[sorted_AUC_bestindic,idx_bi] = sort(AUC_s(:,idx(end)));
[sorted_AUC_uni_bestindic,idx_uni_bi] = sort(AUC_s_uni(:,idx_uni(end)));

subplot(1,2,1)
bar(AUC_s(idx_bi,idx(end)), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',name_scenarios(idx_bi))
xlabel('Monitoring scenario'); ylabel('|AUC-0.5|'); title("Multivariate indicators")
ylim([0 0.52])
yticks(linspace(0,0.5,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 15; 
box off

subplot(1,2,2)
bar(AUC_s_uni(idx_uni_bi,idx_uni(end)), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',name_scenarios_uni2(idx_uni_bi))
xlabel('Data type'); ylabel('|AUC-0.5|'); title("Univariate indicators")
ylim([0 0.52])
yticks(linspace(0,0.5,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 15; 
box off

%% barplots for univariate

subplot(2,1,1)
bar(AUC_s_uni(:,1))
set(gca,'xticklabel',name_scenarios_uni)
xlabel('Indicator'); ylabel('AUC'); title("Performance of the autocorrelation for univariate time series")
ax = gca; 
ax.FontSize = 10; 
box off

subplot(2,1,2)
bar(AUC_s_uni(:,2))
set(gca,'xticklabel',name_scenarios_uni)
xlabel('Indicator'); ylabel('AUC'); title("Performance of the variance for univariate time series")
ax = gca; 
ax.FontSize = 10; 
box off

%% multibarplots for all scenarios

%sort by increasing order
avg_AUCs_allsc=mean(AUC_s,1);
[sorted_AUC,idx] = sort(avg_AUCs_allsc);

bar(AUC_s(:,idx).')
set(gca,'xticklabel',indicators_multi(idx))
xlabel('Indicator'); ylabel('AUC'); title("Performance of the different indicators");
ax = gca; 
ax.FontSize = 10; 
box off

%% calculations 

mean(AUC_s_uni,'all') %avg perf all uni
mean(AUC_s_uni,1) %avg uni for [ AR   var ]

mean(AUC_s(1:3,:),'all') %avg all multi
mean(AUC_s(1:3,id_var_based),'all') %avg multi var based indic
mean(AUC_s(1:3,id_autocor_based),'all') %avg multi autocor based indic