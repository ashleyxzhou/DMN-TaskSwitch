function figures_task_switch()
%Ashley Zhou, June 11
%Script for generating stats and figures for manuscript 
%'DMN activation when switching tasks reflects mental task-set structure'

dbstop if error;
datadir = 'Z:\Duncan-lab\users\az01\task_switch2\Results';
filename  = dir(datadir);

%Each subject's avg rt per switch condition in 2- or 4- domain blocks
load('rt_two_domain.mat');
load('rt_four_domain.mat');

all_switch_names = {
    'rest'
    'task_stay'
    'within_domain'
    'within_chunk_between_domain'
    'between_domain'
    'between_chunk_between_domain'
    'restart'};

load('rt_first_chunk.mat');
load('rt_second_chunk.mat');

chunk_rt_bar = [nanmean(first_chunk_rt,2) nanmean(second_chunk_rt,2)];

%comparing rt in 1st and 2nd learned tasks, collapsed across runs
[h,p,ci,tstat]=ttest(nanmean(first_chunk_rt,1),nanmean(second_chunk_rt,1));
fprintf('Ttest of 1st-learned tasks compared to 2nd learned tasks: t=%.3f, p=%.3f, BF=%.3e',tstat.tstat,p,t1smpbf(tstat.tstat,tstat.df+1)); 

%graph the mean rts per condition
figure(1); clf(1)
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

%% Core univariate figures
load('all_conditions.mat');
%data 40 contrasts x 1 roi x 36 subjects
%load('Z:\Duncan-lab\users\az01\task_switch2\DataAnalysis\aa5_analysis_030323_unsmoothed_GLMs\complexity_conditions\grouping_by_run_aamod_firstlevel_contrasts_djm_00003\all_conditions.mat');
core_data = squeeze(data(:,1,:));
con_names = dataspecs.contrasts;

all_comp.domain = {'two', 'four'};
all_comp.group_learnt = {'first', 'second'};
all_switch_conds={'task_stay','within_domain','between_domain','between_chunk','restart','rest'};
type= fieldnames(all_comp);
comb_anova=nan(36,2,2,2); %sub x cond x domain x chunk
bd_anova=nan(36,2,2,2);
figdata=nan(36,5);
for cond=1:6
    
    %find and average the two
    num_rows=find(startsWith(con_names,[all_switch_conds{cond} '_']));
    avg_cond = mean(core_data(num_rows,:));
    
    
    figdata(:,cond)=avg_cond;
    
end

%subtract between-subject variance 
%figdata=figdata-repmat(mean(figdata,2),1,6)+mean(figdata(:));
figdata=figdata-repmat(figdata(:,1),1,6);

[h,p,ci,stat]=ttest(figdata);
unvbar=mean(figdata);
error=abs(ci-repmat(mean(figdata),2,1));
figure(2);clf(2);b=bar(unvbar(2:end));hold on;
b.FaceColor = 'flat';
%color1=[216/255 217/255 122/255];
%color2=[0.58 0.76 0.43];
orange=[239,138,71]; color1=orange/255;
yellow=[247,170,88]; 
paleyellow=[255,208,111];color2=paleyellow/255;
paleblue=[70,220,224];
lightblue=[55,103,149];
blue=[114,188,213];darkblue=[30,70,110]; color3=[blue;lightblue;darkblue];
% color3=[0.45 0.78 0.76];
% Change the color of a specific bar (e.g., the second bar)
b.CData(1, :) = color1; %[255/255 242/255 204/255]; % RGB color value
b.CData(2, :) = color2; %[177/255 208/255 149/255]; % RGB color value
b.CData(3:5,:)=color3/255; %repmat(color3,3,1);%[224/255 234/255 246/255]
scatter(reshape(repmat(1:5,36,1),36*5,1),reshape(figdata(:,2:6),36*5,1),15,'k','MarkerEdgeAlpha',0.3);
ebar=errorbar(1:5,unvbar(2:end),error(1,2:end),error(2,2:end),'k');
ebar.LineStyle='none';box off;
labname=regexprep(all_switch_conds(2:end),'_','-');
labname{1}= ['within-group-\newline' labname{1}];
labname{2}= ['within-group-\newline' labname{2}];
labname{3}='between-group-\newlinebetween-domain';
set(gca,'xticklabel',labname);
ylabel('BOLD response relative to task-repeat trials');

[h,p,ci,stat]=ttest(mean(figdata(:,2:5),2)); %against task stay
fprintf('\nt-test for averaged task switch conditions compared to task repeats.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[h,p,ci,stat]=ttest(mean(figdata(:,2:5),2),figdata(:,6)); %against rest
fprintf('\nt-test for averaged task switch conditions compared to rest.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[h,p,ci,stat]=ttest(figdata(:,3),figdata(:,4)); %against rest
fprintf('\nt-test for between-domain-within-chunk compared to between-domain-between-chunk.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[h,p,ci,stat]=ttest(figdata(:,2),figdata(:,3)); %wd and bd
fprintf('\nt-test for within-domain compared to between-domain-within-chunk.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[h,p,ci,stat]=ttest(figdata(:,5),figdata(:,2)); %restart and wd
fprintf('\nt-test for restart compared to within-domain.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[h,p,ci,stat]=ttest(figdata(:,5),figdata(:,3)); %restart and bd
fprintf('\nt-test for restart compared to to between-domain-within-chunk.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[h,p,ci,stat]=ttest(figdata(:,5),figdata(:,4)); %restart and bgbd
fprintf('\nt-test for restart compared to between-domain-between-chunk.\nt=%.3f,p=%.3f,BF=%.3e\n',stat.tstat,p,t1smpbf(stat.tstat,36));

all_switch_conds_comb_bd_except_rest = {'within_domain','between_domain'};
bd_conds = {'within_domain','between_domain'};
bc_conds={'between_chunk_between_domain','between_domain'};
for cond=1:2
    for domain = 1:2
        for chunk =1:2
            
            cn=[all_switch_conds_comb_bd_except_rest{cond} '_' all_comp.domain{domain} '_' all_comp.group_learnt{chunk}];
            %find and average the two
            num_rows=find(contains(con_names,cn));
            avg_cond = mean(core_data(num_rows,:));
            
            comb_anova(:,cond,domain,chunk)=avg_cond;
            
            bdcn=[bd_conds{cond} '_' all_comp.domain{domain} '_' all_comp.group_learnt{chunk}];
            num_rows=find(startsWith(con_names,bdcn));
            avg_cond = mean(core_data(num_rows,:));
            bd_anova(:,cond,domain,chunk)=avg_cond;
            
            if domain==2
                bdcn=[bc_conds{cond} '_' all_comp.domain{domain} '_' all_comp.group_learnt{chunk}];
                num_rows=find(startsWith(con_names,bdcn));
                avg_cond = mean(core_data(num_rows,:));
                four_dom_anova(:,cond,chunk)=avg_cond;
                
            end
        end
    end
end

%between-subject factor
load('participant_groups.mat')
%[tbl_comb_anovaa,rm]=simple_mixed_anova(comb_anova,group_vector,{'Condition','Domains','Chunk'},{'Participant_subgroup'});
[tbl,rm]=simple_mixed_anova(bd_anova,group_vector,{'Condition','Domains','Chunk'},{'Participant_subgroup'});  %only within-domain and within-chunk-between-domain
disp(tbl);
% F=tbl.F(5);x=tbl.DF(5);y=tbl.DF(5+1); % domain
F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+2); % domain

BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nDomain F(%d,%d)=%.3e, p=%.3e, BF=%.3e\n',x,y,F,tbl.pValueGG(7),BF10);

%F=tbl.F(9);x=tbl.DF(9);y=tbl.DF(9+1); % interaction of cond and domain
F=tbl.F(13);x=tbl.DF(13);y=tbl.DF(13+2); % interaction of cond and domain

BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction of domain and Condition F(%d,%d)=%.3e, p=%.3e, BF=%.3e\n',x,y,F,tbl.pValueGG(13),BF10);

%F=tbl.F(11);x=tbl.DF(11);y=tbl.DF(11+1); % interaction of cond and chunk
F=tbl.F(16);x=tbl.DF(16);y=tbl.DF(16+2);
BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction of chunk and Condition  F(%d,%d)=%.3e, p=%.3e,BF=%.3e\n',x,y,F,tbl.pValueGG(16),BF10);

%F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % chunk
F=tbl.F(10);x=tbl.DF(10);y=tbl.DF(10+2); % chunk 
BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nChunk  F(%d,%d)=%.3e, p=%.3e,BF=%.3e\n',x,y,F,tbl.pValueGG(10),BF10);
% between-subject factor
F=tbl.F(2);x=tbl.DF(2);y=tbl.DF(3); % interaction of cond and chunk
BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nBetween-subject factor  F(%d,%d)=%.3e, p=%.3e,BF=%.3e\n',x,y,F,tbl.pValueGG(2),BF10);


[tbl,rm]=simple_mixed_anova(four_dom_anova,group_vector,{'Condition','Chunk'},{'Participant_subgroup'});
disp(tbl)
% F=tbl.F(1);x=tbl.DF(1);y=tbl.DF(1+1); % condition
F=tbl.F(4);x=tbl.DF(4);y=tbl.DF(4+2); % condition
BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition F(%d,%d)=%.3e, p=%.3e,BF=%.3e\n',x,y,F,tbl.pValueGG(4),BF10);

% F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % chunk
F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+2); % chunk
BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nChunk F(%d,%d)=%.3e, p=%.3e,BF=%.3e\n',x,y,F,tbl.pValueGG(7),BF10);

% F=tbl.F(5);x=tbl.DF(5);y=tbl.DF(5+1); % interaction of cond and chunk
F=tbl.F(10);x=tbl.DF(10);y=tbl.DF(10+2); % interaction of cond and chunk
BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction of cond and chunk: F(%d,%d)=%.3e, p=%.3e,BF=%.3e\n',x,y,F,tbl.pValueGG(10),BF10);

colorsets=nan(2,2,3);
colorsets(1,:,:)=[255/255 242/255 204/255;177/255 208/255 149/255];
%colorsets(2,:,:)=[0/255 32/255 96/255;237/255 125/255 49/255];
colorsets(2,:,:)=[255/255 242/255 204/255;224/255 234/255 246/255];

%color1=[255/255 242/255 204/255];
%color2=[177/255 208/255 149/255];
legend_names={{'two domain runs','four domain runs'},{'group learnt first','group learnt second'}};

for graph = 1:2
    graph_name = regexprep(type{graph},'_',' ');
    groups = all_comp.(type{graph});
    tobar=[];
    error=[];
    
    %averages the ones with same switch condition and grouping
    for g=1:2
        for cond =1:6
            num_rows=find(startsWith(con_names,strcat(all_switch_conds{cond},'_')).*contains(con_names,groups{g}));
            avg_cond = mean(core_data(num_rows,:));
            group_all_conditions(cond,:)=avg_cond;
            %plot four domain bd and bcbds separated by group
            if graph==2
                num_rows=find(startsWith(con_names,strcat(all_switch_conds{cond},'_')).*contains(con_names,['four_' groups{g}]));
                avg_cond = mean(core_data(num_rows,:));
                four_domain_conds(cond,:)=avg_cond;
            end
        end
        
        group_against_task_stay = group_all_conditions(2:end,:)-group_all_conditions(1,:);
        group_avg = mean(group_against_task_stay,2);
        [H,P,CI,stats]=ttest(group_against_task_stay');
        error1 =abs(CI-repmat(group_avg',2,1));
        
        tobar=[tobar group_avg];
        error =[error; error1];
        anovamatrix(:,:,g) = group_against_task_stay(1:2,:)';
        
        if graph==2 %chunk, no rest
            %average the task switches (bd,bcbd and restart)
            chunk_anova_matrix(:,:,g) = group_against_task_stay(2:4,:)';
            four_against_task_stay =four_domain_conds(2:end,:)-four_domain_conds(1,:);
            four_domain_matrix(:,:,g)=four_against_task_stay(2:3,:)';
        end
        
    end
    
    figure(3);subplot(1,2,graph);
    h=bar(tobar(1:2,:)'); hold on;
    set(h(1),'FaceColor',color1);
    set(h(2),'FaceColor',color2);
    
    anovamatrix(:,:,1)=anovamatrix(:,:,1)-repmat(mean(anovamatrix(:,:,1),2),1,2)+repmat(mean(reshape(anovamatrix(:,:,1),36*2,1)),36,2);
    anovamatrix(:,:,2)=anovamatrix(:,:,2)-repmat(mean(anovamatrix(:,:,2),2),1,2)+repmat(mean(reshape(anovamatrix(:,:,2),36*2,1)),36,2);
    %add individual points
    scatter(repmat(0.85,36,1),anovamatrix(:,1,1),10,'MarkerFaceColor',color1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);
    scatter(repmat(1.15,36,1),anovamatrix(:,2,1),10,'MarkerFaceColor',color2,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);
    scatter(repmat(1.85,36,1),anovamatrix(:,1,2),10,'MarkerFaceColor',color1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);
    scatter(repmat(2.15,36,1),anovamatrix(:,2,2),10,'MarkerFaceColor',color2,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);
    set(gca,'xticklabel',legend_names{graph});
    %     xlabel(graph_name);
    
    etest1=errorbar([0.85,1.15],tobar(1:2,1),error(1,1:2),error(2,1:2),'k');
    etest2=errorbar([1.85,2.15],tobar(1:2,2),error(3,1:2),error(4,1:2),'k');
    etest1.LineStyle='none';
    etest2.LineStyle='none';
    %ymin=min(reshape(anovamatrix,36*4,1));
    %ymax=max(reshape(anovamatrix,36*4,1));
    %ylim([ymin ymax]);
    ylim([-0.5 0.75]);
    
    box off;
   
end

%for first chunk and second chunk, ttest between wd and bd
[~,p,ci,stat]=ttest(anovamatrix(:,1,1),anovamatrix(:,2,1)); %first chunk
fprintf('\nFirst chunk: t-test for within-domain compared to between-domain.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));

[~,p,ci,stat]=ttest(anovamatrix(:,1,2),anovamatrix(:,2,2)); %s chunk
fprintf('\nSecond chunk: t-test for within-domain compared to between-domain.\nt=%.3f,p=%.3f,BF=%.3e',stat.tstat,p,t1smpbf(stat.tstat,36));


return

function bf10 = t1smpbf(t,n,r)
%
% bf10 = t1smpbf(t,n,[r=0.707])
%
% Calculates JZS Bayes Factor for a one-sample t-test given t and sample size n.
% The optional input r is the scale factor which defaults to 0.707.
% This quantifies the evidence in favour of the alternative hypothesis. 
% See Rouder et al, 2009, Psychon Bull Rev for details.
%

% Default scale factor
if nargin < 3
    r = 0.707;
end

% Function to be integrated
F = @(g,t,n,r) (1+n.*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+n.*g.*r.^2).*(n-1))).^(-n./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));

% Bayes factor calculation
bf01=nan(size(t));
for i=1:numel(t)
    bf01(i) = (1 + t(i).^2./(n-1)).^(-n/2) ./ integral(@(g) F(g,t(i),n,r),0,Inf);
end

% Invert Bayes Factor
bf10 = 1 ./ bf01;

function [tbl,rm] = simple_mixed_anova(datamat, varargin)

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names,
% between_factor_names)
%
% Repeated-measures or mixed ANOVA with any number of factors. 
%
% Function built on top of existing MATLAB functions in an attempt to
% simplify user inputs and manipulations while still being able to perform
% all common ANOVA designs.
%
% DATAMAT is a numerical matrix containing the values associated with each 
% level of each within-subject factor in each subject (responses). The 
% subjects should always be the first dimension of the matrix. Each 
% within-subject factor is another dimension.
%
% BETWEEN_FACTORS is a numerical matrix where each row represents a subject 
% and each column represents a between-subjects factor. A given value
% represents the level of the column's factor associated with the row's 
% subject. Optional.
%
% WITHIN_FACTOR_NAMES is a cell array of strings indicating the name of
% each within-subject factor (one for each dimension of the datamat
% variable except the first, in order). Optional.
%
% BETWEEN_FACTOR_NAMES is a cell array of strings indicating the name of
% each between-subjects factor (one for each column of the between_factors
% matrix, in order). These factors are assumed to be categorical (groups). 
% Optional.
%
% TABLE is a table indicating the F statistics, p-values and other
% statistics associated with each term of the model. The "intercept" terms
% can be ignored (e.g. "(Intercept):WS01" indicates the main effect of 
% WS01).
%
% RM is a structure with the repeated measures model parameters and
% statistics.
%
% Does not support covariates or partial models (without all interactions)
% for now.
%
%
% EXAMPLE
%
% A design with 24 subjects and 2 within-subject factors, the first one 
% having 3 levels (time: pre-test, 1st post-test, 2nd post-test) and the 
% second one 4 levels (experimental condition: A, B, C, D). The subjects 
% are grouped in 4 groups: 2 variables with 2 levels each (gender: male or
% female; age group: young or old).
%
% The input datamat should be a 24 x 3 x 4 matrix with each row
% corresponding to a subject. So the element datamat(1,1,1) will correspond
% to the response of the subject #1 in the pre-test in experimental
% condition A, the element datamat(2,3,2) will correspond to the response
% of the subject #2 in the 2nd post-test in experimental condition B, and
% so on.
%
% The input between_factors will be a 24 x 2 matrix with each row
% corresponding to a subject and each column to a between-subjects factor
% (gender and age). Each column will be filled with 1s and 2s, or
% 0s and 1s, or other numbers, indicating the gender/age group of the 
% respective subject.
%
% tbl = simple_mixed_anova(datamat, between_factors, {'Time', 'Exp_cond'},
% {'Gender', 'Age_group'})
%
% Copyright 2017, Laurent Caplette
% https://www.researchgate.net/profile/Laurent_Caplette

% Check if correct number of inputs
narginchk(1,4)

% Assign inputs to variables; if none, will be empty array
between_factors = [];
within_factor_names = [];
between_factor_names = [];
if nargin>1    
    between_factors = varargin{1};
    if nargin>2
        within_factor_names = varargin{2};
        if nargin>3
            between_factor_names = varargin{3};
        end
    end
end

% Determine numbers of variables and measures
nWithin = ndims(datamat)-1;
nBetween = size(between_factors,2);
nVars = size(datamat);
nVars = nVars(2:end); % don't use the nb of subjects on the first dim
nMeas = prod(nVars);

% Check if dimensions of matrices are ok
if size(datamat,1)<2
    error('There must be more than one subject.')
end
if ~isempty(between_factors)
    if size(between_factors,1)~=size(datamat,1)
        error('Both input matrices must have the same nb of subjects.')
    end
end

% Check if there is more than one unique value
if length(unique(datamat))<2
    error('The data matrix must contain more than one unique value.')
end
for ii = 1:size(between_factors,2)
    if length(unique(between_factors(:,ii)))<2
        error('Each between-subjects factor must contain more than one unique value.')
    end
end

% Error if more variable names than variables as input
if length(between_factor_names)>nBetween
    error('Too many between-subject factor names or not enough between-subject variables as input.')
end
if length(within_factor_names)>nWithin
    error('Too many within-subject factor names or not enough within-subject variables as input.')
end

% Check validity of variable names
for ii = 1:length(between_factor_names)
    if ~isvarname(between_factor_names{ii})
        error('Variable names must be continuous strings starting with a letter and without symbols.')
    end
end
for ii = 1:length(within_factor_names)
    if ~isvarname(within_factor_names{ii})
        error('Variable names must be continuous strings starting with a letter and without symbols.')
    end
end

% Assign variable names if not enough or empty
if length(between_factor_names)<nBetween
    nMissing = nBetween - length(between_factor_names);
    BS = repmat('BS', [nMissing 1]); % list of 'BS'
    missing_factor_names = cellstr([BS num2str([1:nMissing]', '%02.0f')]);
    between_factor_names = [between_factor_names missing_factor_names];
end
if length(within_factor_names)<nWithin
    nMissing = nWithin - length(within_factor_names);
    WS = repmat('WS', [nMissing 1]); % list of 'WS'
    missing_factor_names = cellstr([WS num2str([1:nMissing]', '%02.0f')]);
    within_factor_names = [within_factor_names missing_factor_names];
end

% Create table detailing within-subject design
withinVarLevels = fullfact(nVars); % all level combinations
within_table = array2table(withinVarLevels, 'VariableNames', within_factor_names);
for ii = 1:nWithin % ensure that each within-subject factor is categorical (levels==discrete)
    evalc(sprintf('within_table.%s = categorical(within_table.%s)', within_factor_names{ii}, within_factor_names{ii}));
end

% Vectorize all dimensions after first one of the data matrix
y = datamat(:,:);

% Create data table
yList = repmat('Y', [nMeas 1]); % list of 'Y'
numList = num2str([1:nMeas]', '%03.0f'); % support up to 999 measures
measureNames = cellstr([yList numList]); % create names for every measure
for ii = 1:nBetween % add between-subject factors
    measureNames{nMeas+ii} = between_factor_names{ii};
end
total_table = array2table([y between_factors],'VariableNames', measureNames);
for ii = 1:nBetween % ensure that each between-subject factor is categorical (levels/groups==discrete)
    evalc(sprintf('total_table.%s = categorical(total_table.%s)', between_factor_names{ii}, between_factor_names{ii}));
end

% Create between-subjects model using Wilkinson notation
betweenModel = '';
for ii = 1:nBetween
    betweenModel = [betweenModel,measureNames{nMeas+ii},'*'];
end
betweenModel = betweenModel(1:end-1); % remove last star
if isempty(betweenModel)
    betweenModel = '1'; % if no between-subjects factor, put constant term (usually implicit)
end

% Create within-subject model using Wilkinson notation
withinModel = '';
for ii = 1:nWithin
    withinModel = [withinModel,within_factor_names{ii},'*']; % stars for full model (all interactions)
end
withinModel = withinModel(1:end-1); % remove last star

% Fit repeated measures model
rm = fitrm(total_table, sprintf('%s-%s~%s', measureNames{1}, measureNames{nMeas}, betweenModel),...
    'WithinDesign', within_table);

% Run ANOVA
tbl = ranova(rm, 'WithinModel', withinModel);

