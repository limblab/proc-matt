%
clear; clc; close all;
% fname = 'F:\Chewie\CerebusData\2017-01-27\Chewie_M1_stimtest_01272017_001.nev';
% cds=commonDataStructure();
% cds.file2cds(fname,6,'arrayM1','monkeyChewie','taskRW','ranByMatt','mapFileZ:\lab_folder\Animal-Miscellany\Chewie 8I2\Chewie Left M1 SN 6250-001474.cmp','ignoreJumps')
% % fname = 'F:\Chewie\CerebusData\2016-10-21\Chewie_PMd_RT_CS_BL_10212016_002.nev';
% % cds.file2cds(fname,6,'arrayPMd','monkeyChewie','taskRW','ranByMatt','mapFileZ:\lab_folder\Animal-Miscellany\Chewie 8I2\Chewie Left PMd SN 6251-001469.cmp','ignoreJumps')
% save('F:\Chewie\CDS\2017-01-27\Chewie_M1_stimtest_01272017_001.mat','cds','-v7.3')

% clear; clc; close all;
% fname = 'F:\Mihili\CerebusData\2017-02-01\Mihili_M1_stim_02012017_001.nev';
% cds=commonDataStructure();
% cds.file2cds(fname,6,'arrayM1','monkeyMihili','taskRW','ranByMatt','mapFile\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left M1 SN 1025-001452.cmp','ignoreJumps')
% % fname = 'F:\Chewie\CerebusData\2016-10-21\Chewie_PMd_RT_CS_BL_10212016_002.nev';
% % cds.file2cds(fname,6,'arrayPMd','monkeyChewie','taskRW','ranByMatt','mapFileZ:\lab_folder\Animal-Miscellany\Chewie 8I2\Chewie Left PMd SN 6251-001469.cmp','ignoreJumps')
% save('F:\Mihili\CDS\2017-02-01\Mihili_M1_stim_02012017_001.mat','cds','-v7.3')

% fname = 'F:\Han\CerebusData\2017-01-26\Han_20170126_CObumpcurl_washout_area2EMG_003.nev';
% cds=commonDataStructure();
% cds.file2cds(fname,6,'arrayS1','monkeyHan','taskCObump','ranByRaeed','mapFile\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp','ignoreJumps')
% % fname = 'F:\Chewie\CerebusData\2016-10-21\Chewie_PMd_RT_CS_BL_10212016_002.nev';
% % cds.file2cds(fname,6,'arrayPMd','monkeyChewie','taskRW','ranByMatt','mapFileZ:\lab_folder\Animal-Miscellany\Chewie 8I2\Chewie Left PMd SN 6251-001469.cmp','ignoreJumps')
% save('F:\Han\CDS\2017-01-26\Han_20170126_CObumpcurl_washout_area2EMG_003.mat','cds','-v7.3')



%%
num_pulses = 1000; % per channel
win_size = [0,10.1]/1000; % in ms
bin_size = 1/1000;

plot_colors = distinguishable_colors(length(bin_size:bin_size:win_size(2)));

if 1 % get stim from spikes
    stim_ts = cds.units(end).spikes.ts;
else % get stim from 30k data
    s = cds.analog{2}.StimTrig;
    ds = [0; diff(s)];
    dds = [0; diff(ds)];
    peaks = [dds(1:end-1)>9 & dds(2:end)<-9; false];
end
idx_start = find(diff(stim_ts) > 1);

num_rand_samps = 1000;

%%
% filter out low-firing units and get spike probability
unit_fr = zeros(1,length(cds.units)-1);
for unit = 1:length(cds.units)-1
    ts = cds.units(unit).spikes.ts;
    idx = zeros(size(ts));
    for i = 1:length(stim_ts)
        idx = idx | ts >= stim_ts(i) & ts <= stim_ts(i)+0.005;
    end
    ts = ts(~idx);
    unit_fr(unit) = length(ts)/cds.meta.duration;
end

%%
% 
for chan = 2:length(idx_start)
    idx = idx_start(chan);
    for unit = 1:length(cds.units)-1
        if unit_fr(unit) > 40
            disp(['stim ' num2str(chan) '; unit ' num2str(unit) '...']);
            stim_N = zeros(num_pulses,length(bin_size:bin_size:win_size(2)));
            rand_N = zeros(num_pulses,length(bin_size:bin_size:win_size(2)));
            
            % get spikes in each bin BEFORE stimulus
            [rand_wave,stim_wave,rand_wave_bin,stim_wave_bin] = deal([]);
            count = 0;
            for i = 1:num_pulses
                ts_idx = find(cds.units(unit).spikes.ts > stim_ts(idx+i)-win_size(2) & cds.units(unit).spikes.ts < stim_ts(idx+i)-0.0001);
                ts = cds.units(unit).spikes.ts(ts_idx);
                % put in small bins
                rand_N(i,:) = histcounts(ts,stim_ts(idx+i)-win_size(2):bin_size:stim_ts(idx+i));
                if sum(ts_idx) > 0
                    for j = 1:length(ts_idx)
                        count = count + 1;
                        rand_wave(count,:) = mean(cds.units(unit).spikes{ts_idx(j),2},1);
                        rand_wave_bin(count) = find(histcounts(ts(j),stim_ts(idx+i)-win_size(2):bin_size:stim_ts(idx+i)));
                    end
                end
            end
            
            % get spikes in each bin AFTER stimulus
            count = 0;
            stim_delay = [];
            for i = 1:num_pulses
                ts_idx = find(cds.units(unit).spikes.ts > stim_ts(idx+i)+win_size(1) & cds.units(unit).spikes.ts < stim_ts(idx+i)+win_size(2)-0.0001);
                ts = cds.units(unit).spikes.ts(ts_idx);
                % put in small bins
                stim_N(i,:) = histcounts(ts,stim_ts(idx+i):bin_size:stim_ts(idx+i)+win_size(2));
                
                if ~isempty(ts)
                    stim_delay = [stim_delay, ts(1)-stim_ts(idx+i)];
                end
                
                if ~isempty(ts_idx)
                    for j = 1:length(ts_idx)
                        count = count + 1;
                        stim_wave(count,:) = mean(cds.units(unit).spikes{ts_idx(j),2},1);
                        stim_wave_bin(count) = find(histcounts(ts(j),stim_ts(idx+i):bin_size:stim_ts(idx+i)+win_size(2)));
                    end
                end
                
            end
            
            figure;
            subplot(2,2,[1,3]); hold all;
            

            
            rand_N = reshape(rand_N,numel(rand_N),1);
            boot_idx = randi(numel(rand_N),numel(rand_N),1000);
            rand_err = prctile(sum(rand_N(boot_idx),1)./numel(rand_N),[2.5,97.5])';
            rand_err = repmat(rand_err,1,size(stim_N,2));
            patch(bin_size*[1:size(stim_N,2), fliplr(1:size(stim_N,2))],[rand_err(1,:), fliplr(rand_err(2,:))],'k','FaceAlpha',0.4,'EdgeAlpha',0);
                        
            p = zeros(1,size(stim_N,2));
            for i = 1:size(stim_N,2)
                [~,p(i)] = ttest2(rand_N,stim_N(:,i));
            end
            
            % bootstrap sum to get error bars
            boot_idx = randi(size(stim_N,1),size(stim_N,1),1000);
            stim_err = zeros(2,size(stim_N,2));
            for i = 1:size(stim_N,2)
                temp = stim_N(:,i);
                temp = temp(boot_idx);
                stim_err(:,i) = prctile(sum(temp,1)./size(temp,1),[2.5 97.5]);
            end
            
            if any(p < 0.05) %any(stim_err(1,:) > rand_err(2,:))
                for i = 1:size(stim_N,2)
                    plot(bin_size*[i; i],[stim_err(1,i); stim_err(2,i)],'LineWidth',2,'Color',plot_colors(i,:));
                    plot(bin_size*i, sum(stim_N(:,i),1)./size(stim_N,1),'o', 'LineWidth',2,'Color',plot_colors(i,:));
                end
                set(gca,'Box','off','TickDir','out','FontSize',14);
                axis('tight');
                title(['stim ' num2str(chan) '; unit ' num2str(unit)],'FontSize',14);
                
                subplot(2,2,2); hold all;
                for i = 1:size(rand_wave,1)
                    plot(rand_wave(i,:),'LineWidth',1,'Color','k');
                end
                set(gca,'Box','off','TickDir','out','FontSize',14);
                axis('tight');
                
                subplot(2,2,4); hold all;
                for i = 1:size(stim_wave,1)
                    plot(stim_wave(i,:),'LineWidth',1,'Color',plot_colors(stim_wave_bin(i),:));
                end
                set(gca,'Box','off','TickDir','out','FontSize',14);
                axis('tight');
                
                [~,s] = pca(stim_wave);
                figure;
                hold all;
                for i = 1:size(stim_wave,1)
                    plot3(s(i,1),s(i,2),s(i,3),'o','Color',plot_colors(stim_wave_bin(i),:));
                end
                
                pause;
            end
            close all;
        end
    end
end
