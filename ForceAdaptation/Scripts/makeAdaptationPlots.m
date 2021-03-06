% Starting from scratch on the behavioral adaptation plotting code

monkeyColors = {'r','b','g'};
% Do two things:
%   1) For each session, plot a trace of all trials (aligned as percentage) on top of each other
%   2) Pool across sessions in blocks to match the neural data
epochs = {'BL','AD','WO'};
useMetrics = {'errors'};
numTrials = 200;
remOutliers = false;

for iTask = 1:length(useTasks)
    switch lower(useMetrics{iTask})
        case 'curvatures'
            valScale = 1;
            x_min = -0.2;
            x_max = 0.2;
        case 'errors'
            valScale = (180/pi);
            x_min = -35;
            x_max = 35;
    end
    
    h = figure('Position',[400 400 650 550]);
    hold all;
    for iMonkey = 1:length(useMonkeys)
        plot(-1,0,'o','Color',monkeyColors{iMonkey},'LineWidth',3);
    end
    
    for iMonkey = 1:length(useMonkeys)
        if strcmpi(useMonkeys{iMonkey},'MrT')
            useArray = 'PMd';
        else
            useArray = 'M1';
        end
        % start with just CO files
        useFiles = doFiles(strcmpi(doFiles(:,4),useTasks{iTask}) & strcmpi(doFiles(:,1),useMonkeys{iMonkey}),:);
        
        pertDir = zeros(1,size(useFiles,1));
        results = cell(size(useFiles,1),length(epochs));
        targInfo = cell(size(useFiles,1),length(epochs));
        for iFile = 1:size(useFiles,1)
            a = loadResults(root_dir,useFiles(iFile,:),'adaptation');
            
            % get direction of perturbation to flip the clockwise ones to align
            if flipClockwisePerts
                % gotta hack it
                dataPath = fullfile(root_dir,useFiles{iFile,1},'Processed',useFiles{iFile,2});
                expParamFile = fullfile(dataPath,[useFiles{iFile,2} '_experiment_parameters.dat']);
                a.BL.params.exp = parseExpParams(expParamFile);
                switch lower(a.BL.params.exp.angle_dir)
                    case 'ccw'
                        pertDir(iFile) = 1;
                    case 'cw'
                        pertDir(iFile) = -1;
                end
            else
                pertDir(iFile) = 1;
            end
            
            for iEpoch = 1:length(epochs)
                results{iFile,iEpoch} = pertDir(iFile).*a.(epochs{iEpoch}).(useMetrics{iTask}).*valScale;
                targInfo{iFile,iEpoch} = a.(epochs{iEpoch}).movement_table(:,1);
            end
        end
        
        % find minimum length of trials for each epoch
        %         trialMins = numTrials.*ones(1,length(epochs));
%         trialMins = min(cellfun(@(x) length(x),results),[],1);
        trialMins = [100, 128, 189];
        
        if remOutliers
            disp('Removing outliers at 3x std of baseline...');
        end
        
        % get the baseline error for each unique target direction
        blErr = zeros(size(useFiles,1),8);
        for iFile = 1:size(useFiles,1)
            t = targInfo{iFile,1};
            r = results{iFile,1};
            utheta = unique(t);
            for i = 1:length(utheta)
                idx = t == utheta(i);
                blErr(iFile,i) = mean(r(idx));
            end
        end
        
        % Plot first trials from each session as an average
        figure('Position',[400 400 650 550]);
        subplot1(1,length(epochs));
        for iEpoch = 1:length(epochs)
            subplot1(iEpoch);
            hold all;
            
            set(gca,'Box','off','TickDir','out','FontSize',14);
            if iEpoch ~= 1
                set(gca,'YTick',[]);
            else
                title(useMonkeys{iMonkey},'FontSize',16);
                ylabel('Angular Error (Deg)','FontSize',16);
            end
            
            sResults = cell(1,size(results,1));
            for iFile = 1:size(results,1)
                r = results{iFile,iEpoch};
                t = targInfo{iFile,iEpoch};
                
                % remove outliers
                if remOutliers
                    badInds = abs(r) > 3*std(results{iFile,1});
                    if ~isempty(badInds)
                        r(badInds) = [];
                        t(badInds) = [];
                    end
                end
                
                try
                    r = r(1:trialMins(iEpoch));
                    t = t(1:trialMins(iEpoch));
                catch
                    r = [r; NaN(trialMins(iEpoch)-length(r),1)];
                    t = [t; NaN(trialMins(iEpoch)-length(t),1)];
                end
                
                for j = 1:length(r)
                    idx = utheta == t(j);
                    if ~isnan(r(j))
                        r(j) = r(j) - blErr(iFile,idx);
                    end
                end
                sResults{iFile} = r;
            end
            
            % plot all of the session traces in gray
            %                 plot(cell2mat(sResults),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
            
            % plot mean across sessions for each trial
            %             plot(repmat(1:trialMins(iEpoch),size(sResults,1),1),cell2mat(sResults),'o','LineWidth',1,'Color',[0.7 0.7 0.7],'MarkerSize',4);
            
            fuck_matlab = nanmean(cell2mat(sResults),2);
            plot(1:trialMins(iEpoch),fuck_matlab,'b+','LineWidth',2);
            
            % fit trendline and plot it
            [fit,~,~,~,s] = regress(fuck_matlab,[ones(trialMins(iEpoch),1), (1:trialMins(iEpoch))']);
            hold all; 
            plot(1:trialMins(iEpoch),fit(1)+(1:trialMins(iEpoch))*fit(2),'k-','LineWidth',3);
            
            
            %             plot(1:trialMins(iEpoch),mean(cell2mat(sResults),2) - std(cell2mat(sResults),[],2)./sqrt(size(results,1)),'k--');
            %             plot(1:trialMins(iEpoch),mean(cell2mat(sResults),2) + std(cell2mat(sResults),[],2)./sqrt(size(results,1)),'k--');
            
            axis('tight');
            set(gca,'YLim',[x_min,x_max]);
            xlabel('Trial Number','FontSize',16);
        end
        
        % Take tuning blocks of trials and compute average metrics
        
        % load something to get the parameters used for the main tuning
        c = loadResults(root_dir,useFiles(iFile,:),'tuning',{'classes'},useArray,paramSetName,tuneMethod,tuneWindow);
        blocks = c.params.tuning.blocks; clear c;
        
%         blocks = {[0,0.5;0.5,1], ...
%             [0,0.1;0.1,0.2;0.2,0.3;0.3,0.4;0.4,0.5;0.5,0.6;0.6,0.7;0.7,0.8;0.8,0.9;0.9,1], ...
%             [0,0.1;0.1,0.2;0.2,0.3;0.3,0.4;0.4,0.5;0.5,0.6;0.6,0.7;0.7,0.8;0.8,0.9;0.9,1]};
        
        % deduce how many total blocks
        sResults = cell(size(results,1),sum(cellfun(@(x) length(x)-1,blocks)));
        
        for iFile = 1:size(results,1)
            count = 0;
            for iEpoch = 1:length(epochs)
                r = results{iFile,iEpoch};
                t = targInfo{iFile,iEpoch};
                for iBlock = 1:length(blocks{iEpoch})-1
                    count = count + 1;
                    % get relevant trials
                    idx1 = floor(blocks{iEpoch}(iBlock)*length(r))+1;
                    idx2 = floor(blocks{iEpoch}(iBlock+1)*length(r));
                    
                    temp1 = r(idx1:idx2);
                    temp2 = t(idx1:idx2);
                    for j = 1:length(temp1)
                        idx = utheta == temp2(j);
                        temp1(j) = temp1(j) - blErr(iFile,idx);
                    end
                    sResults{iFile,count} = temp1;
                end
            end
        end
        
%         figure('Position',[400 400 650 550]);
%         figure(h);
figure;
        hold all;
        
        m = zeros(1,size(sResults,2));
        s = zeros(2,size(sResults,2));
        for iBlock = 1:size(sResults,2)
            %temp = cell2mat(sResults(:,iBlock));
            temp = cellfun(@mean,sResults(:,iBlock));
            
            m(iBlock) = mean(temp);
            s(:,iBlock) = [mean(temp)-std(temp)./sqrt(length(temp)),mean(temp)+std(temp)./sqrt(length(temp))];
            
            %             plot(iBlock,m,'o','Color',monkeyColors{iMonkey},'LineWidth',2);
            %             plot([iBlock,iBlock],s,'-','Color',monkeyColors{iMonkey},'LineWidth',2);
            
        end
        
        all_x = 1:size(sResults,2);
        
        plot(all_x,m,'o','Color',monkeyColors{iMonkey},'LineWidth',2);
        plot([all_x; all_x],s,'-','Color',monkeyColors{iMonkey},'LineWidth',2);
        
        %         patch([all_x, fliplr(all_x)],[s(1,:), fliplr(s(2,:))],monkeyColors{iMonkey},'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',monkeyColors{iMonkey});
        %         h(i) = plot(all_x,m,'-','LineWidth',3,'Color',monkeyColors{iMonkey});
        
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[1 size(sResults,2)],'YLim',[-16,10]);
        %         title(useMonkeys{iMonkey},'FontSize',16);
        xlabel('Blocks','FontSize',16);
        ylabel(useMetrics{iTask},'FontSize',16);
    end
%     legend(useMonkeys,'FontSize',14);
end
