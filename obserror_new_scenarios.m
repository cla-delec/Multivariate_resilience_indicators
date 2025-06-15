  %% LOAD THE DATA

% data_tp=readtable('data/siminf_extrares_feedingprefB10_ts3000_100rep_k0.65-0.92_R00.7-1_rateimport0.15.csv');
% data_tp=table2array(data_tp);
% data_tp=data_tp.';
% 
% data_tn=readtable('data/siminf_extrares_feedingprefB10_ts3000_100rep_k0.75_R00.8_rateimport0.15.csv');
% data_tn=table2array(data_tn);
% data_tn=data_tn.';

data_tp=readtable('data/siminf_extrares_feedingprefB80_ts3000_100rep_k0.35-0.50_R00.7-1_rateimport0.15.csv');
data_tp=table2array(data_tp);
data_tp=data_tp.';

data_tn=readtable('data/siminf_extrares_feedingprefB80_ts3000_100rep_k0.40_R00.8_rateimport0.15.csv');
data_tn=table2array(data_tn);
data_tn=data_tn.';

%% SIMULATION PARAMETERS

number_dim=10;
number_rep=size(data_tp,2)/number_dim;

scenarios_multi=[1 2 7 9; 3 8 10 nan; 7 8 9 10; 1 2 3 4; 1 5 6 nan]; %scenarios for multivariate indicators
id_incidence=[1 2 4 5 7 9]; %also dead because are aggregated by taking max
id_sero=[3 6 8 10];  
n_scenarios_multi=size(scenarios_multi,1);

scenarios_uni=(1:1:number_dim).';
n_scenarios_uni=size(scenarios_uni,1);

indicators_uni=["AR","std"];
indicators_multi=["degfinger", "mafac", "mafeig", "mean_ar","max_ar","mafvar","max_var","mean_var","pcavar","explvar","max_abscorr","max_cov"];

id_autocor_based=[1,2,3,4,5]; n_indic_autocor=size(id_autocor_based,2);
%id_autocor_based=[4, 5, 2]; n_indic_autocor=size(id_autocor_based,2);
id_var_based=[6,7,8,9,10,12]; n_indic_var=size(id_var_based,2);
%id_var_based=[8, 9, 7]; n_indic_var=size(id_var_based,2);
n_indic_uni=size(indicators_uni,2); n_indic_multi=size(indicators_multi,2);

length_tot=length(data_tp);
length_tot_aggregated=ceil(length_tot/7);

%lengths_tested=linspace(100, length_tot_aggregated, n_lengths);
obs_rate_tested=[1, 0.1, 0.01, 0.001, 0.0001]; 
n_rates=size(obs_rate_tested,2);

name_scenarios=["Incidence scenario", "Seroprevalence scenario", "Anthro-equine scenario", "Wildlife scenario", "Hidden reservoir scenario"];
name_scenarios_uni=["M_I", "Ba_I", "Ba_R", "Ba_D", "Bb_I", "Bb_R", "H_I", "H_R", "E_I", "E_R"];
name_scenarios_uni2=["Mosquito incidence", "Bird V incidence", "Bird V seroprevalence", "Bird V dead", "Bird H incidence", "Bird H seroprevalence", "Human incidence", "Human seroprevalence", "Horse incidence", "Horse seroprevalence"];


%% calculate all taus

results_taus_agg_zsnans_tp_multi=zeros(n_indic_multi,n_scenarios_multi,number_rep,n_rates); 
results_taus_agg_zsnans_tn_multi=zeros(n_indic_multi,n_scenarios_multi,number_rep,n_rates); 
results_taus_agg_zsnans_tp_uni=zeros(n_indic_uni,n_scenarios_uni,number_rep,n_rates); 
results_taus_agg_zsnans_tn_uni=zeros(n_indic_uni,n_scenarios_uni,number_rep,n_rates); 

for ind_cur_rate=1:n_rates
    
    cur_rate=obs_rate_tested(ind_cur_rate);
                
    for cur_rep=1:number_rep
        
        for cur_scenario_multi=1:n_scenarios_multi
            
            %for tp data
            num_time_series=sum(~isnan(scenarios_multi(cur_scenario_multi,:)));
            index_cur_scenario=scenarios_multi(cur_scenario_multi,1:num_time_series);
            index_scenario_rep=index_cur_scenario+7*(cur_rep-1);
            cur_data_no_process_tp=data_tp(:,index_scenario_rep); %subset from the big dataset

            %find inc and sero ts in current scenario
            inc_ts=find(ismember(index_cur_scenario,id_incidence));

            cur_data_processed_tp=aggregate_week(cur_data_no_process_tp,inc_ts);
            cur_data_processed_tp_obs_process=poissrnd(cur_data_processed_tp*cur_rate);
            
            ews=generic_ews(cur_data_processed_tp_obs_process,'indicators',{'degfinger', 'mafac', 'mafeig', 'mean_ar','max_ar','mafvar','max_var','mean_var','pcavar','explvar','max_abscorr','max_cov'},'datacolumn',[],'silent',true,'nanflag','omitnan'); %,'ebisuzaki',100
            cur_taus=ews.taus;
            results_taus_agg_zsnans_tp_multi(:,cur_scenario_multi,cur_rep,ind_cur_rate)=cur_taus;

            %for tn data
            num_time_series=sum(~isnan(scenarios_multi(cur_scenario_multi,:)));
            index_cur_scenario=scenarios_multi(cur_scenario_multi,1:num_time_series);
            index_scenario_rep=index_cur_scenario+7*(cur_rep-1);
            cur_data_no_process_tn=data_tn(:,index_scenario_rep); %subset from the big dataset

            cur_data_processed_tn=aggregate_week(cur_data_no_process_tn,inc_ts);
            cur_data_processed_tn_obs_process=poissrnd(cur_data_processed_tn*cur_rate);


            ews=generic_ews(cur_data_processed_tn_obs_process,'indicators',{'degfinger', 'mafac', 'mafeig', 'mean_ar','max_ar','mafvar','max_var','mean_var','pcavar','explvar','max_abscorr','max_cov'},'datacolumn',[],'silent',true,'nanflag','omitnan'); %,'ebisuzaki',100
            cur_taus=ews.taus;
            results_taus_agg_zsnans_tn_multi(:,cur_scenario_multi,cur_rep,ind_cur_rate)=cur_taus;
        end
        
        for cur_scenario_uni=1:n_scenarios_uni
            
            %for tp data

            index_scenario_rep_uni=scenarios_uni(cur_scenario_uni)+7*(cur_rep-1);
            cur_data_tp_uni=data_tp(:,index_scenario_rep_uni); %subset from the big dataset
            
            %find inc and sero ts in current scenario and aggregate
            if ismember(scenarios_uni(cur_scenario_uni),id_incidence)
                cur_data_tp_uni_processed_tp=aggregate_week(cur_data_tp_uni,1);
            else
                cur_data_tp_uni_processed_tp=aggregate_week(cur_data_tp_uni,[]);
            end

            cur_data_tp_uni_processed_obs_process=poissrnd(cur_data_tp_uni_processed_tp*cur_rate);
            
            ews=generic_ews(cur_data_tp_uni_processed_obs_process,'indicators',{'var','ar'},'silent',true,'nanflag','omitnan'); %,'ebisuzaki',100
            cur_taus=ews.taus;
            results_taus_agg_zsnans_tp_uni(:,cur_scenario_uni,cur_rep,ind_cur_rate)=cur_taus;

            %for tn data
            cur_data_tn_uni=data_tn(:,index_scenario_rep_uni); %subset from the big dataset
            
            %find inc and sero ts in current scenario and aggregate
            if ismember(scenarios_uni(cur_scenario_uni),id_incidence)
                cur_data_tn_uni_processed_tn=aggregate_week(cur_data_tn_uni,1);
            else
                cur_data_tn_uni_processed_tn=aggregate_week(cur_data_tn_uni,[]);
            end

            cur_data_tn_uni_processed_obs_process=poissrnd(cur_data_tn_uni_processed_tn*cur_rate);
            
            ews=generic_ews(cur_data_tn_uni_processed_obs_process,'indicators',{'std','ar'},'silent',true,'nanflag','omitnan'); %,'ebisuzaki',100
            cur_taus=ews.taus;
            results_taus_agg_zsnans_tn_uni(:,cur_scenario_uni,cur_rep,ind_cur_rate)=cur_taus;
        end
        
    end

end


final_results.results_uni_tp=results_taus_agg_zsnans_tp_uni;
final_results.results_uni_tn=results_taus_agg_zsnans_tn_uni;
final_results.results_multi_tp=results_taus_agg_zsnans_tp_multi;
final_results.results_multi_tn=results_taus_agg_zsnans_tn_multi;

save("data/results_extrares_feedingpref80_obs_process_5rates1_0.0001_alltaus_siminf_length3000_100rep_tp_k0.35-0.50_R00.7-1_tn_k0.40_R00.8_rateimport0.15_agg.mat",'final_results');

%% load results

final_result=load("data/results_extrares_feedingpref80_obs_process_5rates1_0.0001_alltaus_siminf_length3000_100rep_tp_k0.35-0.50_R00.7-1_tn_k0.40_R00.8_rateimport0.15_agg.mat");


results_taus_agg_zsnans_tp_uni=final_result.final_results.results_uni_tp;
results_taus_agg_zsnans_tn_uni=final_result.final_results.results_uni_tn;
results_taus_agg_zsnans_tp_multi=final_result.final_results.results_multi_tp;
results_taus_agg_zsnans_tn_multi=final_result.final_results.results_multi_tn;



%% calculate AUCs

results_multi=zeros(n_indic_multi,n_rates,n_scenarios_multi); 
result_uni=zeros(n_indic_uni,n_rates,n_scenarios_uni);

label=[{'epidemics'},{'no epidemics'}];
label_tot=repelem(label,number_rep);

for ind_cur_rate=1:n_rates
    for cur_scenario_multi=1:n_scenarios_multi
        for cur_indic_multi=1:n_indic_multi
            res_tp_scen_indic_length=squeeze(results_taus_agg_zsnans_tp_multi(cur_indic_multi,cur_scenario_multi,:,ind_cur_rate));
            res_tn_scen_indic_length=squeeze(results_taus_agg_zsnans_tn_multi(cur_indic_multi,cur_scenario_multi,:,ind_cur_rate));

            [X_cur,Y_cur,T,results_multi(cur_indic_multi,ind_cur_rate,cur_scenario_multi)]=perfcurve(label_tot,[res_tp_scen_indic_length.',res_tn_scen_indic_length.'],'epidemics','ProcessNaN','addtofalse'); 
        end
    end
    
    for cur_scenario_uni=1:n_scenarios_uni    
        for cur_indic_uni=1:n_indic_uni
            res_tp_scen_indic_length=squeeze(results_taus_agg_zsnans_tp_uni(cur_indic_uni,cur_scenario_uni,:,ind_cur_rate));
            res_tn_scen_indic_length=squeeze(results_taus_agg_zsnans_tn_uni(cur_indic_uni,cur_scenario_uni,:,ind_cur_rate));

            [X_cur,Y_cur,T,result_uni(cur_indic_uni,ind_cur_rate,cur_scenario_uni)]=perfcurve(label_tot,[res_tp_scen_indic_length.',res_tn_scen_indic_length.'],'epidemics','ProcessNaN','addtofalse'); 
        end
    end
end

final_aucs.results_uni_aucs=result_uni;
final_aucs.results_multi_aucs=results_multi;

save("data/results2_aucs_extrares_obs_process_5rates1_0.0001_alltaus_siminf_feedingpref80_length3000_100rep_tp_k0.35-0.50_R00.7-1_tn_k0.40_R00.8_rateimport0.15_agg.mat",'final_aucs');

%% LOAD RESULTS

final_res=load("data/results2_aucs_extrares_obs_process_5rates1_0.0001_alltaus_siminf_feedingpref80_length3000_100rep_tp_k0.35-0.50_R00.7-1_tn_k0.40_R00.8_rateimport0.15_agg.mat");
 
result_uni=final_res.final_aucs.results_uni_aucs;
results_multi=final_res.final_aucs.results_multi_aucs;

%% PLOT RESULTS
 
indicators_multi=["Deg Finger", "MAF AC", "MAF eig", "Mean AR","Max AR","MAF var","Max var","Mean var","PCA var","Expl var","Max abscorr","Max cov"];

cols=["#0072BD","#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F","#FF00FF","#0000FF","#000000"];

for k=1:n_scenarios_multi
    %subplot(2,5,k)
    subplot(5,2,1+(k-1)*2)
    %title({name_scenarios(k),"AR-based indicators"})
    title(name_scenarios(k))
    xlabel('Observation rate'); ylabel('AUC');
    for i=1:n_indic_autocor
        id_indic_autocor=id_autocor_based(i);
        hold on
        plot(obs_rate_tested,results_multi(id_indic_autocor,:,k),'o-','LineWidth',2);
        
    end
    %xticklabels({'Weekly', 'Biweekly', ' ', 'Monthly'})
    
    id_uni_cur_multi=find(ismember(scenarios_uni,scenarios_multi(k,:)));
   
    for j=1:size(id_uni_cur_multi,1)
        id_cur_uni=id_uni_cur_multi(j);
        plot(obs_rate_tested,result_uni(2,:,id_cur_uni),'o--','LineWidth',2,'Color',cols(id_cur_uni));
    end
    %legend([indicators_multi(id_autocor_based) name_scenarios_uni2(id_uni_cur_multi)],'Location','southoutside')
    xlim([0.1 1.1])
    ylim([0.5 1])
    yticks(0.5:0.1:1);
    ax = gca; 
    ax.FontSize = 10; 
    grid on
    hold off
    
    subplot(5,2,2+(k-1)*2) 
    %subplot(5,2,k+5)
    %title({name_scenarios(k),"Var-based indicators"})
    title(name_scenarios(k))
    xlabel('Observation rate'); ylabel('AUC');
    for i=1:n_indic_var
        id_indic_var=id_var_based(i);
        hold on
        plot(obs_rate_tested,results_multi(id_indic_var,:,k),'o-','LineWidth',2);
    end
    %xticklabels({'Weekly', 'Biweekly', ' ', 'Monthly'})   
    for j=1:size(id_uni_cur_multi,1)
        id_cur_uni=id_uni_cur_multi(j);
        plot(obs_rate_tested,result_uni(1,:,j),'o--','LineWidth',2,'Color',cols(id_cur_uni));
    end
    %legend([indicators_multi(id_var_based) name_scenarios_uni2(id_uni_cur_multi)],'Location','southoutside')
    xlim([0.1 1.1])
    ylim([0.5 1])
    yticks(0.5:0.1:1);
    ax = gca; 
    ax.FontSize = 10; 
    grid on
    hold off
    
    
end


%% PLOT RESULTS CUR SCENARIO
 
indicators_multi=["Deg Finger", "MAF AC", "MAF eig", "Mean AR","Max AR","MAF var","Max var","Mean var","PCA var","Expl var","Max abscorr","Max cov"];

cols=["#0072BD","#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F","#FF00FF","#0000FF","#000000"];

k=5; %scenario

subplot(1,2,1)
title(name_scenarios(k))
xlabel('Observation rate'); ylabel('AUC');
for i=1:n_indic_autocor
    id_indic_autocor=id_autocor_based(i);
    hold on
    plot(obs_rate_tested,results_multi(id_indic_autocor,:,k),'^-','LineWidth',2);

end



id_uni_cur_multi=find(ismember(scenarios_uni,scenarios_multi(k,:)));

for j=1:size(id_uni_cur_multi,1)
    id_cur_uni=id_uni_cur_multi(j);
    plot(obs_rate_tested,result_uni(2,:,id_cur_uni),'^--','LineWidth',2,'Color',cols(id_cur_uni));
end
legend([indicators_multi(id_autocor_based) name_scenarios_uni2(id_uni_cur_multi)],'Location','southoutside')
xlim([-0.01 1.1])
ylim([0.5 1])
yticks(0.5:0.1:1);
set(gca, 'XScale', 'log')
ax = gca; 
ax.FontSize = 10; 
grid on
hold off

subplot(1,2,2) 
title(name_scenarios(k))
xlabel('Observation rate'); ylabel('AUC');
for i=1:n_indic_var
    id_indic_var=id_var_based(i);
    hold on
    plot(obs_rate_tested,results_multi(id_indic_var,:,k),'o-','LineWidth',2);
end

for j=1:size(id_uni_cur_multi,1)
    id_cur_uni=id_uni_cur_multi(j);
    plot(obs_rate_tested,result_uni(1,:,j),'o--','LineWidth',2,'Color',cols(id_cur_uni));
end
legend([indicators_multi(id_var_based) name_scenarios_uni2(id_uni_cur_multi)],'Location','southoutside')
xlim([0 1.1])
ylim([0.5 1])
yticks(0.5:0.1:1);
set(gca, 'XScale', 'log')
ax = gca; 
ax.FontSize = 10; 
grid on
hold off
    
    
%% NEW FIGURE - main
 
indicators_multi=["Deg Finger", "MAF AR", "MAF eig", "Mean AR","Max AR","MAF var","Max var","Mean var","PCA var","Expl var","Max abscorr","Max cov"];

k=5; %scenario

id_selected_indic=[8, 9, 7];
n_selected=size(id_selected_indic,2);

%add y-jitter to prevent overlap
results_multi(9,:,k)=results_multi(9,:,k)-0.005;

title(name_scenarios(k))
xlabel('Observation rate'); ylabel('AUC');
for i=1:n_selected
    id_indic_var=id_selected_indic(i);
    hold on
    plot(obs_rate_tested(1:4),results_multi(id_indic_var,1:4,k),'o-','LineWidth',2);
end


id_selected_ts=[1, 5, 6, 7];
n_selected_uni=size(id_selected_ts,2);

cols_uni=["#F26B6C","#FF9A00","#C78520","#6BCB77"];


for j=1:n_selected_uni
    id_cur_uni=id_selected_ts(j);
    plot(obs_rate_tested(1:4),result_uni(1,1:4,j),'o--','LineWidth',2,'Color',cols_uni(j));
end

legend([indicators_multi(id_selected_indic) name_scenarios_uni2(id_selected_ts)],'Location','southoutside')
set(gca, 'XScale', 'log')
xlim([0.0007 1.5])
ylim([0.37 1])
yticks(0.4:0.1:1);
ax = gca; 
ax.FontSize = 10; 
grid on
hold off
    

%% Figure supplementary - univariate time series
 
cols=["#FF6B6B","#FFD93D", "#FFD93D", "#FFD93D", "#F2BE22", "#F2BE22", "#6BCB77","#6BCB77","#4D96FF","#4D96FF"];
linetypes=["--^","--^","--square","--o","--^","--square","--^","--square","--^","--square"];

subplot(1,2,1)
title("Autocorrelation")
xlabel('Observation rate'); ylabel('AUC');
hold on
for j=1:n_scenarios_uni
    plot(obs_rate_tested(1:4),result_uni(2,1:4,j),linetypes(j),'LineWidth',2,'Color',cols(j));
end
hold off
legend(name_scenarios_uni2,'Location','southoutside')
xlim([0.0007 1.5])
xticks([0.001 0.01 0.1 1])
ylim([0.20 1])
yticks(0.2:0.1:1);
set(gca, 'XScale', 'log')
ax = gca; 
ax.FontSize = 10; 
grid on
hold off

subplot(1,2,2) 
title("Variance")
xlabel('Observation rate'); ylabel('AUC');
hold on
for j=1:n_scenarios_uni

    plot(obs_rate_tested(1:4),result_uni(1,1:4,j),linetypes(j),'LineWidth',2,'Color',cols(j));
end
hold off
legend(name_scenarios_uni2,'Location','southoutside')
xlim([0.0007 1.5])
xticks([0.001 0.01 0.1 1])
ylim([0.20 1])
yticks(0.2:0.1:1);
set(gca, 'XScale', 'log')
ax = gca; 
ax.FontSize = 10; 
grid on
hold off

%% Figure supplementary - multivariate indicators for all scenarios
 
indicators_multi=["Deg Finger", "MAF AC", "MAF eig", "Mean AR","Max AR","MAF var","Max var","Mean var","PCA var","Expl var","Max abscorr","Max cov"];

k=3; %scenario

subplot(1,2,1)
title(strcat(name_scenarios(k), " autocorrelation-based indicators"))
xlabel('Observation rate'); ylabel('AUC');
for i=1:n_indic_autocor
    id_indic_autocor=id_autocor_based(i);
    hold on
    plot(obs_rate_tested(1:4),results_multi(id_indic_autocor,1:4,k),'^-','LineWidth',2);

end

legend(indicators_multi(id_autocor_based),'Location','southoutside')
xlim([0.0007 1.5])
xticks([0.001 0.01 0.1 1])
ylim([0.10 1])
yticks(0.1:0.1:1);
set(gca, 'XScale', 'log')
ax = gca; 
ax.FontSize = 10; 
grid on
hold off

subplot(1,2,2) 
title(strcat(name_scenarios(k), " variance-based indicators"))
xlabel('Observation rate'); ylabel('AUC');
for i=1:n_indic_var
    id_indic_var=id_var_based(i);
    hold on
    plot(obs_rate_tested(1:4),results_multi(id_indic_var,1:4,k),'o-','LineWidth',2);
end

legend(indicators_multi(id_var_based),'Location','southoutside')
xlim([0.0007 1.5])
xticks([0.001 0.01 0.1 1])
ylim([0.10 1])
yticks(0.1:0.1:1);
set(gca, 'XScale', 'log')
ax = gca; 
ax.FontSize = 10; 
grid on
hold off


%% PLOT RESULTS summary
 
indicators_multi=["Deg Finger", "MAF AC", "MAF eig", "Mean AR","Max AR","MAF var","Max var","Mean var","PCA var","Expl var","Max abscorr","Max cov"];
id_var_based2=[6 7 8 9 12];

k=4;

res_AUCs_multi_arindic_st_ress=abs(mean(results_multi(id_autocor_based,:,k),1)-0.5);
res_AUCs_multi_varindic_st_ress=abs(mean(results_multi(id_var_based2,:,k),1)-0.5);
AUCs_uni_st_ress=abs(result_uni(:,:,:)-0.5);

 
xlabel('Sampling frequency', 'FontSize',24); ylabel('Prediction performance (AUC)', 'FontSize',24);
hold on
%plot both types of indic
plot(obs_rate_tested,mean(results_multi(id_autocor_based,:,k),1),'-o','LineWidth',2,'Color',"#0072BD",'MarkerSize',7);
plot(obs_rate_tested,mean(results_multi(id_var_based2,:,k),1),'-s','LineWidth',2,'Color',"#0072BD",'MarkerSize',10,'MarkerFaceColor',"#0072BD");
xticklabels({'Weekly','', 'Biweekly', ' ','','', 'Monthly'})
yticks(0.5:0.1:1)

cols=["#77AC30" "#D95319" "#EDB120"];
%plot ar uni
id_uni_cur_multi=find(ismember(scenarios_uni,scenarios_multi(k,:)));
for j=1:size(id_uni_cur_multi,1)
    plot(obs_rate_tested,result_uni(2,:,j),'--o','LineWidth',2,'Color', cols(j),'MarkerSize',7);
end

%plot var uni
for j=1:size(id_uni_cur_multi,1)
    plot(obs_rate_tested,result_uni(1,:,j),'--s','LineWidth',2,'Color', cols(j),'MarkerSize',10,'MarkerFaceColor',cols(j));
end

legend(["Multivariate indicators - AR-based" "Multivariate indicators - var-based" name_scenarios_uni2(id_uni_cur_multi)+repelem(" Autocorrelation", 3) name_scenarios_uni2(id_uni_cur_multi)+repelem(" variance", 3)],'Location','southoutside') %,'Orientation','horizontal'
legend(["Multivariate indicators - AR-based" "Multivariate indicators - var-based" name_scenarios_uni2(id_uni_cur_multi) name_scenarios_uni2(id_uni_cur_multi)],'Location','southoutside') %,'Orientation','horizontal'
ylim([0.5 1])
xlim([0.3 1.1])
ax = gca; 
ax.FontSize = 16; 
hold off

    

%% barplot tot - to check if consistent with results auc

%sort by increasing order
avg_AUCs_perindic=mean(squeeze(results_multi(:,1,:)),2);
[sorted_AUC2,idx_indic] = sort(avg_AUCs_perindic);


subplot(5,2,1)
bar(results_multi(idx_indic,1,1), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx_indic))
xlabel('Indicator'); ylabel('Prediction performance (AUC)'); title(name_scenarios(1))
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,2)
b=bar(squeeze(result_uni(:,1,scenarios_multi(1,:))).')
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(1,:)))
xlabel('Univariate time series'); ylabel('Prediction performance (AUC)'); %title("Performance of the univariate indicators")
%legend(["Variance","Autocorrelation"],'Location','eastoutside')
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,3)
bar(results_multi(idx_indic,1,2), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx_indic))
xlabel('Indicator'); ylabel('Prediction performance (AUC)'); title(name_scenarios(2))
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,4)
b=bar(squeeze(result_uni(:,1,scenarios_multi(2,1:3))).')
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(2,1:3)))
xlabel('Univariate time series'); ylabel('Prediction performance (AUC)'); %title("Performance of the univariate indicators")
%legend(["Variance","Autocorrelation"],'Location','eastoutside')
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,5)
bar(results_multi(idx_indic,1,3), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx_indic))
xlabel('Indicator'); ylabel('Prediction performance (AUC)'); title(name_scenarios(3))
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,6)
b=bar(squeeze(result_uni(:,1,scenarios_multi(3,:))).')
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(3,:)))
xlabel('Univariate time series'); ylabel('Prediction performance (AUC)'); %title("Performance of the univariate indicators")
%legend(["Variance","Autocorrelation"],'Location','eastoutside')
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,7)
bar(results_multi(idx_indic,1,4), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx_indic))
xlabel('Indicator'); ylabel('Prediction performance (AUC)'); title(name_scenarios(4))
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,8)
b=bar(squeeze(result_uni(:,1,scenarios_multi(4,1:4))).')
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(4,1:4)))
xlabel('Univariate time series'); ylabel('Prediction performance (AUC)'); %title("Performance of the univariate indicators")
%legend(["Variance","Autocorrelation"],'Location','eastoutside')
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off


subplot(5,2,9)
bar(results_multi(idx_indic,1,5), 'FaceColor','#6895D2')
grid on
set(gca,'xticklabel',indicators_multi(idx_indic))
xlabel('Indicator'); ylabel('Prediction performance (AUC)'); title(name_scenarios(5))
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off

subplot(5,2,10)
b=bar(squeeze(result_uni(:,1,scenarios_multi(5,1:3))).')
b(1).FaceColor='#D04848';
b(2).FaceColor='#F3B95F';
grid on
set(gca,'xticklabel',name_scenarios_uni(scenarios_multi(5,1:3)))
xlabel('Univariate time series'); ylabel('Prediction performance (AUC)'); %title("Performance of the univariate indicators")
%legend(["Variance","Autocorrelation"],'Location','eastoutside')
ylim([0.5 1])
yticks(linspace(0.5,1,6))
xtickangle(45)
ax = gca; 
ax.FontSize = 10; 
box off