function behav_taskswitch2()

dbstop if error;
datadir = 'Z:\Group\Duncan-lab\users\az01\task_switch2\Results';
filename  = dir(datadir);
addpath 'Z:\Group\Duncan-lab\users\az01\task_switch2\DataAnalysis'
addpath 'Z:\Group\Duncan-lab\users\dm01\MoreTools'
%store each subject's avg rt per switch condition in 2- or 4- domain blocks
prev_sub_num ='';
run=7;
sub=0;
%
allsubs={};
two_domain = nan(6,36);
four_domain = nan(7,36);

all_switch_names = {
    'rest'
    'task_stay'
    'within_domain'
    'within_chunk_between_domain'
    'between_domain'
    'between_chunk_between_domain'
    'restart'};

all_task_names = {'a1','a2','b1','b2','c1','c2','d1','d2','r'};
sub_mean_type_rt = nan(9,36);
sub_mean_type_acc= nan(9,36);
first_chunk_rt = nan(6,36);
second_chunk_rt = nan(6,36);
ci_switch=nan(1,36);
for num =1:length(filename)
    
    if strcmp(filename(num).name,'.') || strcmp(filename(num).name,'..')
        continue;
    end
    
    comp = split(filename(num).name,'_');
    subnum = comp{2};
    run_num = split(comp{4},'.');
    run_num = str2num(run_num{1});
    allsubs=[allsubs subnum];
    if strcmp(subnum,'230103')|| strcmp(subnum,'230231') 
        %fprintf('found and excluded %s',subnum);
        continue;
    end
    %fprintf('Retrieved data from subject %s',subnum);
    
    %initialize subject's rts per condition
    for i = 1:length(all_switch_names)
        sub_rts(run_num).(all_switch_names{i})=[];
        sub_acc(run_num).(all_switch_names{i})=[];
        
    end
    for i = 1:length(all_task_names)
        sub_type_rts(run_num).(all_task_names{i})=[];
        sub_type_acc(run_num).(all_task_names{i})=[];
    end
    
    if strcmp(subnum, prev_sub_num)
        run=run+1;
    else
        %makes sure each subject has 6 runs of data, outputs subject if not
        if run<6
            fprintf([subnum ' has only ' run_num ' runs']);
            
        end
        run = 1;
        prev_sub_num = subnum;
        sub=sub+1;
        
        if sub>1
            for i = 1:length(all_switch_names)
                %compute each run's mean rt and mean acc for each participant
                sub_mean_rt(sub-1).(all_switch_names{i}) = mean([sub_rts.(all_switch_names{i})]);
            end
            
            for i =  1:length(all_task_names)
                sub_mean_type_rt(i,sub) = mean([sub_type_rts.(all_task_names{i})]);
                sub_mean_type_acc(i,sub) = mean([sub_type_acc.(all_task_names{i})]);
            end
        end
        
    end
    
    %for each subject, average the rts and accuracy of each task type, and
    %each switch type, per run
    load(fullfile(datadir, filename(num).name),'result');
    %switch_names = unique({result.switch_type});
    task_names = unique({result.type});
    
    
    first_chunk_sub=[];
    second_chunk_sub=[];
    for trial=1:length(result)
        
        switch_type=regexprep(result(trial).switch_type,'-','_');
        
        if strcmp(switch_type,'dummy_trial')
            continue;
        end
        trial_rt = result(trial).rt - result(trial).stim_onset;
        trial_accuracy = result(trial).accuracy;
%         if any(contains(result(1).first_two_domains,type(1)))
%                 %add rt to first chunk sub avg
%                 first_chunk_sub = [first_chunk_sub trial_rt];
%             elseif any(contains(result(1).second_two_domains,type(1)))
%                 second_chunk_sub = [second_chunk_sub trial_rt];
% 
%             end
        sub_acc(run_num).(switch_type)=[sub_acc(run_num).(switch_type) trial_accuracy];
        sub_type_acc(run_num).(result(trial).type)=[sub_type_acc(run_num).(result(trial).type) trial_accuracy];
        if trial_accuracy
            sub_rts(run_num).(switch_type)=[sub_rts(run_num).(switch_type) trial_rt];
            sub_type_rts(run_num).(result(trial).type)=[sub_type_rts(run_num).(result(trial).type) trial_rt];
            
            ttype=(result(trial).type);
%             if any(contains(result(1).first_two_domains,ttype(1)))
%                 %add rt to first chunk sub avg
%                 first_chunk_sub = [first_chunk_sub trial_rt];
%             elseif any(contains(result(1).second_two_domains,ttype(1)))
%                 second_chunk_sub = [second_chunk_sub trial_rt];
% 
%             end
        end
    end
%     first_chunk_rt(run,sub)=nanmean(first_chunk_sub);
%     second_chunk_rt(run,sub)=nanmean(second_chunk_sub);
        
    
    %average the condition's rts per condition
    for i = 1:length(all_switch_names)
       
        %average the mean rt for two domain and four domain blocks
        if run==1
            two_domain_avg.(all_switch_names{i}) = sub_rts(run_num).(all_switch_names{i});
            two_domain_avg_acc.(all_switch_names{i}) = sub_acc(run_num).(all_switch_names{i});

        elseif run ==3
            four_domain_avg.(all_switch_names{i}) = sub_rts(run_num).(all_switch_names{i});
            four_domain_avg_acc.(all_switch_names{i}) = sub_acc(run_num).(all_switch_names{i});

        elseif run ==4
            four_domain_avg.(all_switch_names{i}) = mean([four_domain_avg.(all_switch_names{i}) sub_rts(run_num).(all_switch_names{i})]);
            four_domain_avg_acc.(all_switch_names{i}) = mean([four_domain_avg_acc.(all_switch_names{i}) sub_acc(run_num).(all_switch_names{i})]);

            %add average to 4-domain block
            four_domain(i,sub)=four_domain_avg.(all_switch_names{i});
            four_domain_acc(i,sub)=four_domain_avg_acc.(all_switch_names{i});
            %four_domain.(switch_names{i}) = [four_domain.(switch_names{i}) four_domain_avg];
        elseif run==6
            two_domain_avg.(all_switch_names{i}) =  mean([two_domain_avg.(all_switch_names{i}) sub_rts(run_num).(all_switch_names{i})]);
            two_domain_avg_acc.(all_switch_names{i}) =  mean([two_domain_avg_acc.(all_switch_names{i}) sub_acc(run_num).(all_switch_names{i})]);

            %add average to 2-domain block
            two_domain(i,sub)=two_domain_avg.(all_switch_names{i});
            two_domain_acc(i,sub)=two_domain_avg_acc.(all_switch_names{i});
            %two_domain.(switch_names{i}) = [two_domain.(switch_names{i}) two_domain_avg];
        else
            two_domain_avg.(all_switch_names{i}) =  [two_domain_avg.(all_switch_names{i}) sub_rts(run_num).(all_switch_names{i})];
            two_domain_avg_acc.(all_switch_names{i}) = [two_domain_avg_acc.(all_switch_names{i}) sub_acc(run_num).(all_switch_names{i})];

        end
    end
    
end
% save('two_domain_acc','two_domain_acc');
% save('four_domain_acc','four_domain_acc');

chunk_rt_bar = [nanmean(first_chunk_rt,2) nanmean(second_chunk_rt,2)];
%participants' average rt per block and per switch

 figure(100);clf(100);
bar(chunk_rt_bar);
[H,P,CI,stats]=ttest(first_chunk_rt')
error_first = abs(CI-repmat(mean(first_chunk_rt,2)',2,1));

[H,P,CI,stats]=ttest(second_chunk_rt')
error_two = abs(CI-repmat(mean(second_chunk_rt,2)',2,1));

hold on;

e1=errorbar([1:6]-0.15,mean(first_chunk_rt,2)',error_first(1,:),error_first(2,:));
e1.LineStyle='none';

e2=errorbar([1:6]+0.15,mean(second_chunk_rt,2)',error_two(1,:),error_two(2,:));
e2.LineStyle='none';

legend('First chunk','Second chunk');
ylabel('RT(s)');
xlabel('Run');

%comparing rt in 1st and 2nd learned tasks, collapsed across runs
[h,p,ci,tstat]=ttest(nanmean(first_chunk_rt,1),nanmean(second_chunk_rt,1))
fprintf('Ttest of 1st-learned tasks compared to 2nd learned tasks: t=%.3f, p=%.3f, BF=%.3e',tstat.tstat,p,t1smpbf(tstat.tstat,tstat.df+1)); 

%graph the mean rts per condition
figure(11); clf(11)
domain_means=[mean(two_domain(2:end,:),2),mean(four_domain(2:end,:),2)];
domain_means(3,:)=nanmean(domain_means(3:4,:));
domain_means=domain_means([1:3,5:6],:);
h=bar(domain_means);
hold on;
set(h(1),'FaceColor',[255/255 242/255 204/255]);
set(h(2),'FaceColor',[177/255 208/255 149/255]);
adj=four_domain'-repmat(nanmean(four_domain',2),1,7)+nanmean(four_domain(:));
adj2=two_domain'-repmat(nanmean(two_domain',2),1,7)+nanmean(two_domain(:));

scatter(reshape(repmat([1,2,3,5]-0.15,36,1),36*4,1),reshape(adj2(:,[2,3,5,7]),36*4,1),10,'MarkerFaceColor',[255/255 242/255 204/255],'MarkerFaceAlpha',0.7,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);
scatter(reshape(repmat([1:5]+0.15,36,1),36*5,1),reshape(adj(:,[2:4,6,7]),36*5,1),10,'MarkerFaceColor',[177/255 208/255 149/255],'MarkerFaceAlpha',0.7,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);

[H,P,CI,stats]=ttest(adj);
error_four = abs(CI-repmat(mean(four_domain,2)',2,1));

[H,P,CI,stats]=ttest(adj2);
error_two = abs(CI-repmat(mean(two_domain,2)',2,1));

box off;
e1=errorbar([1:5]-0.15,domain_means(:,1),error_two(1,[2,3,5:end]),error_two(2,[2,3,5:end]),'k');
e2=errorbar([1:5]+0.15,domain_means(:,2),error_four(1,[2:4,6:end]),error_four(2,[2:4,6:end]),'k');
e1.LineStyle='none';
e2.LineStyle='none';
ylim([1 2.5]);
legend('Two domain runs','Four domain runs');
legend boxoff;
switch_names={'task-repeat','within-\ndomain','between-\ndomain','between-\nchunk','restart'};
switch_names={'task-repeat','within-group\nwithin-domain','within-group\nbetween-domain','between-group\nbetween-domain','restart'};

set(gca,'xticklabel',switch_names);
%xtickangle(30);
ylabel('RT(s)');

figure(12);clf(12)
bar(nanmean(sub_mean_type_rt,2));
[~,p,CI,stats_rt]=ttest(sub_mean_type_rt');
error_rt = abs(CI-repmat(nanmean(sub_mean_type_rt,2)',2,1));
task_type_name = {'Living','Shoebox','Letter A','Letter I','Same shape','Same size','Sex','Age','rest'};
set(gca,'xticklabel',task_type_name);
xtickangle(30);
title('Avg rts for task types');
hold on;
e3=errorbar(1:9,nanmean(sub_mean_type_rt,2),error_rt(1,:),error_rt(2,:));
e3.LineStyle='none';
ylabel('RT(s)');
 xlabel('Task type');

figure(13);clf(13);
bar(nanmean(sub_mean_type_acc,2));
[~,p,CI,stats_rt]=ttest(sub_mean_type_acc');
error_rt = abs(CI-repmat(nanmean(sub_mean_type_acc,2)',2,1));
set(gca,'xticklabel',task_type_name);
xtickangle(30);
title('Avg accuracy for task types');
hold on;
e3=errorbar(1:9,nanmean(sub_mean_type_acc,2),error_rt(1,:),error_rt(2,:));
e3.LineStyle='none';
ylabel('Accuracy');
ylim([0.85 1.05]);
 xlabel('Task type');

 %analysis of behav data
dat = zeros(2,6,36);
dat(1,:,:)=first_chunk_rt;
dat(2,:,:)=second_chunk_rt;
dat=permute(dat,[3,1,2]); %sub x chunk x task type

dom_dat=zeros(2,2,36);  %domain x within/between x sub
%for anova, sub x domain x within/between task switch

dom_dat(1,1,:) = four_domain(3,:);
dom_dat(1,2,:) = mean(four_domain([4,6],:),1);

dom_dat(2,1,:)=two_domain(3,:);
dom_dat(2,2,:)=two_domain(5,:);

dom_dat = permute(dom_dat,[3,1,2]);
[tbl,rm]=simple_mixed_anova(dom_dat,[],{'Domains','Conditions'});
disp(tbl);
F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % roi
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nDomains BF=%.3e\n',BF10);
        F=tbl.F(5);x=tbl.DF(5);y=tbl.DF(5+1); % cond
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition BF=%.3e\n',BF10);
        F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction BF=%.3e\n',BF10);
%t-test of all task switches compared to task stay in four and two domains
[h,p,ci,tstat]=ttest(squeeze(mean(dom_dat(:,1,:),3))',four_domain(2,:));
fprintf('Ttest of four domain task switches compared to task stay: t=%.3f, p=%.3f, BF=%.3e',tstat.tstat,p,t1smpbf(tstat.tstat,tstat.df+1)); 

[h,p,ci,tstat]=ttest(squeeze(mean(dom_dat(:,2,:),3))',two_domain(2,:));
fprintf('Ttest of two domain task switches compared to task stay: t=%.3f, p=%.3f, BF=%.3e',tstat.tstat,p,t1smpbf(tstat.tstat,tstat.df+1)); 
return;