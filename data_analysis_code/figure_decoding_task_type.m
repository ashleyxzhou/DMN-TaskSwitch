function figure_decoding_task_type()

dbstop if error;

%datadir='Z:\Duncan-lab\users\az01\task_switch\DataAnalysis';
datadir = '/group/duncan-lab/users/az01/task_switch2/Results/';
%addpath Z:\Duncan-lab\users\dm01\MoreTools;
addpath /group/duncan-lab/users/dm01/MoreTools;

decoding_dir = '/group/duncan-lab/users/az01/task_switch2/DataAnalysis/aa5_analysis_030323_unsmoothed_GLMs/task_type/aamod_decoding5_az_00001';
%\CBU230093\DecodingAnalysis
measures = {'dprime','balanced_accuracy_minus_chance','signed_decision_values'};

names = dir(decoding_dir);
groups={'TwoDomains',  'FourDomains'};
tasks = {'a1','a2','b1','b2','c1','c2','d1','d2'};
load(fullfile(fullfile(decoding_dir,names(end).name),'DecodingAnalysis/AllResults.mat'),'cfg');
setnames =cfg.design.setnames;
load(fullfile(fullfile(decoding_dir,names(end).name),'DecodingAnalysis/Decoding_dprime.mat'),'results');
roi=results.roi_names;
avg_X.(roi{1})=struct(); avg_X.(roi{2}) =struct();avg_X.(roi{3})=struct();avg_X.(roi{4})=struct();
subjects=0;
for k = 1:length(names)
    subname = names(k).name;
    if strcmp(subname,'.')||strcmp(subname,'..')
        continue;
    end
    subjects = subjects +1;
    toload = fullfile(fullfile(decoding_dir,subname),'DecodingAnalysis');
    load(fullfile(toload,'AllResults.mat'),'cfg');
    
    subnum=split(subname,'CBU');
    load(fullfile(datadir,['participant_' subnum{2} '_run_1.mat']),'result');
    chunk = result.first_two_domains;
    chunk2=result.second_two_domains;
%     chunk_index = contains(tasks,chunk); %1s for all the rows in first chunk,
%     sec_chunk_index=contains(tasks,chunk2);
%     
%     %%To-DO
%     %use regexp to select the withinchunk indexes 
%     %regexp(setnames,[chunk{1} '[12]_vs_' chunk{2} '[12]']);
    indexes=~cellfun(@isempty,regexp(setnames,[chunk{1} '[12]_vs_' chunk{2} '[12]']));
    indexes=indexes+~cellfun(@isempty,regexp(setnames,[chunk{2} '[12]_vs_' chunk{1} '[12]']));
    second_chunk_indexes=~cellfun(@isempty,regexp(setnames,[chunk2{1} '[12]_vs_' chunk2{2} '[12]']));
    second_chunk_indexes=second_chunk_indexes+~cellfun(@isempty,regexp(setnames,[chunk2{2} '[12]_vs_' chunk2{1} '[12]']));
    within_domain_indexes=[1,2,27,28,45,46,55,56];
    between_chunk_indexes=~ismember(1:56,within_domain_indexes)-indexes-second_chunk_indexes;
    
    %     between_chunk_indexes1=~cellfun(@isempty,regexp(setnames,[chunk{1} '[12]_vs_[' chunk2{1} chunk2{2} '][12]']));
%     between_chunk_indexes2=~cellfun(@isempty,regexp(setnames,[chunk{2} '[12]_vs_[' chunk2{1} chunk2{2} '][12]']));
%     between_chunk_indexes3=~cellfun(@isempty,regexp(setnames,[chunk2{1} '[12]_vs_[' chunk{1} chunk{2} '][12]']));
%     between_chunk_indexes4=~cellfun(@isempty,regexp(setnames,[chunk2{2} '[12]_vs_[' chunk{1} chunk{2} '][12]']));

%     for o=1:length(measures)
%         load(fullfile(toload,['Decoding_' measures{o} '.mat']));
%         roiXset=[results.(measures{o}).set.output]; %56x4
%         roi_withinchunk=mean(roiXset(:,indexes+second_chunk_indexes),2);
%         %roi_between_chunk=mean(roiXset(:,~(indexes+second_chunk_indexes),2));
%         %between chunks exclude rest, restart,task-stay,within domain
%         all_sub_within=[all_sub_within roi_withinchunk];
%         all_sub_within
%         for r=1:results.n_decodings
%         roiXset(r,~cellfun(@isempty,regexp(setnames,[chunk{1} '[12]_vs_' chunk{2} '[12]'])))
%         end
%     end
    %%
    
    for o=1:length(measures)
        load(fullfile(toload,['Decoding_' measures{o} '.mat']));
        roiXset=[results.(measures{o}).set.output];
        two_domains(k).(measures{o})=mean(roiXset(:,[1:2:56]),2);
        four_domains(k).(measures{o}) = mean(roiXset(:,[2:2:56]),2);
        
        within_chunk(k).(measures{o})=mean(roiXset(:,find((indexes+second_chunk_indexes)==1)),2);
        between_chunk(k).(measures{o})=mean(roiXset(:,find(between_chunk_indexes)),2);
        
        for r=1:results.n_decodings
            % reshape by class and group
            testset=0;
            X=nan(results.n_cond,results.n_cond,2);
            for a=1:results.n_cond
                for b=(a+1):results.n_cond
                    for g=1:2
                        testset=testset+1;
                        X(a,b,g)=roiXset(r,testset);
                        
                        %if two domain %2=0
                    end
                end
            end
            X=nansum(cat(3,X(:,:,1),X(:,:,2)'),3);
            if isempty(fieldnames(avg_X.(roi{r})))
                
                avg_X.(roi{r}).(measures{o})=X;
            elseif ~contains(fieldnames(avg_X.(roi{r})),measures{o})
                avg_X.(roi{r}).(measures{o})=X;
            else
                avg_X.(roi{r}).(measures{o})=avg_X.(roi{r}).(measures{o})+X;
                
            end
            
            %here average subjects' X for each roi and each measure
        end
        
    end
    
end
figure(99);clf(99);
pos=fitplots2([results.n_decodings, numel(cfg.results.output)-2],'',4);

pos2=fitplots2([numel(cfg.results.output)-2,2],'',4);
for o=1:length(measures)
    
    for r=1:length(roi)
        figure(99);
        ax=subplot('position',pos{r,o});
        X=avg_X.(roi{r}).(measures{o})/subjects;
        imagesc(X); axis tight equal;
        cb=colorbar; cb.Label.String=measures{o}; cb.Label.Interpreter='none';
        colormap(colorcet('D1'));
%         caxis([-1 1]*max(abs(roiXset(:))));
        caxis([-1 1]*max(abs(X(:))));

        set(gca,'XTick',1:results.n_cond,'XTickLabel',tasks,'YTick',1:results.n_cond,'YTickLabel',tasks)
        
        
        
        if o==1, ylabel(regexprep(cfg.files.mask{r},'.*/',''),'Interpreter','none'); end
        
        if o==1 && r==1, title(sprintf('Test on %s at upper right\nTest on %s at lower left',groups{1},groups{2})); end
        drawnow

    end
    
    figure(1);
    ax=subplot('position',pos2{o,1});
    all_two=[two_domains.(measures{o})];
    all_four=[four_domains.(measures{o})];
        bar([mean(all_two,2),mean(all_four,2)]);
    [~,p,CI,stats]=ttest(all_two');
    error1 = abs(CI-repmat(mean(all_two,2)',2,1));
    [~,p,CI,stats]=ttest(all_four');
        error2 = abs(CI-repmat(mean(all_four,2)',2,1));
        
        hold on;
        e1=errorbar([1:length(roi)]-0.15,mean(all_two,2),error1(1,:),error1(2,:));
        e2=errorbar([1:length(roi)]+0.15,mean(all_four,2),error2(1,:),error2(2,:));

        e1.LineStyle='none';
        e2.LineStyle='none';
        
        set(gca,'XTick',1:results.n_cond,'XTickLabel',roi);
        legend('Two domains','Four domains')
        ylabel(measures{o})
        if o==1 title('Decoding accuracy of task pairs in rois by number of domains'); end;
        
        ax=subplot('position',pos2{o,2});
        

    withinchunk=[within_chunk.(measures{o})];
    betweenchunk=[between_chunk.(measures{o})];
        bar([mean(withinchunk,2),mean(betweenchunk,2)]);
    [~,p,CI,stats]=ttest(withinchunk');
    error1 = abs(CI-repmat(mean(withinchunk,2)',2,1));
    [~,p,CI,stats]=ttest(betweenchunk');
        error2 = abs(CI-repmat(mean(betweenchunk,2)',2,1));
        
        hold on;
        e1=errorbar([1:length(roi)]-0.15,mean(withinchunk,2),error1(1,:),error1(2,:));
        e2=errorbar([1:length(roi)]+0.15,mean(betweenchunk,2),error2(1,:),error2(2,:));

        e1.LineStyle='none';
        e2.LineStyle='none';
        
        set(gca,'XTick',1:results.n_cond,'XTickLabel',roi);
        legend('Withink Chunk','Between Chunk')
        ylabel(measures{o})
    if o==1 title('Decoding accuracy of task pairs in rois by chunk'); end;
end
set(findall(99,'-property','FontSize'),'FontSize',10);

% figure(1);


return;