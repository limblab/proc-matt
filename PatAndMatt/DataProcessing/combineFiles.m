function master_data = combineFiles(filepaths)
% filepaths is cell array of file paths
%   Assumes all files came from the same monkey!!!! AND FROM THE SAME DAY!
possibleArrays = {'M1','PMd'};

% load in data for all desired epochs
data = cell(1,length(filepaths));
for iFile = 1:length(filepaths)
    data{iFile} = load(filepaths{iFile});
end

% see what array data is present
useArrays = intersect(possibleArrays,fieldnames(data{1}));

% if multiple arrays are given, pool them together
[allMT,allPos,allVel,allAcc,allForce,allT,fileIDs] = deal([]);
t_end = zeros(1,length(filepaths));
tasks = cell(1,length(filepaths)); %tasks used
perts = cell(1,length(filepaths)); %conditions
epochs = cell(1,length(filepaths)); %task epochs
for iFile = 1:length(filepaths)
    % stitch together kinematics and trial information
    movement_table = data{iFile}.movement_table;
    movement_table(:,2:6) = movement_table(:,2:6)+t_end(iFile);
    
    allT = [allT; data{iFile}.cont.t+t_end(iFile)];
    allMT = [allMT; movement_table];
    allPos = [allPos; data{iFile}.cont.pos];
    allVel = [allVel; data{iFile}.cont.vel];
    allAcc = [allAcc; data{iFile}.cont.acc];
    allForce = [allForce; data{iFile}.cont.force];
    
    t_end(iFile+1) = allT(end);
    
    fileIDs = [fileIDs;iFile*ones(size(movement_table,1),1)];
    
    tasks{iFile} = data{iFile}.meta.task;
    perts{iFile} = data{iFile}.meta.perturbation;
    epochs{iFile} = data{iFile}.meta.epoch;
end
master_data.cont.t = allT;
master_data.cont.pos = allPos;
master_data.cont.vel = allVel;
master_data.cont.acc = allAcc;
master_data.cont.force = allForce;
master_data.movement_table = allMT;
master_data.meta.file_ids = fileIDs;
master_data.meta.t_ends = t_end;
master_data.meta.filepaths = filepaths;
master_data.meta.task = tasks;
master_data.meta.perturbation = perts;
master_data.meta.epoch = epochs;
master_data.meta.monkey = data{1}.meta.monkey;
master_data.meta.recording_date = data{1}.meta.recording_date;
master_data.params = data{1}.params;
clear allMT allPos allVel allAcc allForce allT;

% now loop along arrays and add spikes
for iArray = 1:length(useArrays)
    useArray = useArrays{iArray};
    
    % build master list of present neurons
    bad_units = checkUnitGuides(cellfun(@(x) x.(useArray).sg,data,'UniformOutput',false));
    master_sg = setdiff(data{1}.(useArray).sg, bad_units, 'rows');
    
    master_data.(useArray).sg = master_sg;
    master_data.(useArray).units = repmat(struct('ts',[],'wf',[],'id',[]),1,size(master_sg,1)); % initialize
    
    wf = cell(1,size(master_sg,1));
    for iFile = 1:length(filepaths)
        [~,idx] = intersect( master_sg,data{iFile}.(useArray).sg,'rows');
        
        % loop along neurons and tack on the spikes
        for i = 1:length(idx)
            % should be a single row
            if size(data{iFile}.(useArray).units(idx(i)).ts,1) > 1
                data{iFile}.(useArray).units(idx(i)).ts = data{iFile}.(useArray).units(idx(i)).ts';
            end
            % add end times from before
            master_data.(useArray).units(i).ts = [master_data.(useArray).units(i).ts, data{iFile}.(useArray).units(idx(i)).ts + t_end(iFile)];
            master_data.(useArray).units(i).id = data{iFile}.(useArray).sg(idx(i),:);
            wf{i} = [wf{i}, data{iFile}.(useArray).units(idx(i)).wf];
        end
    end
    
    for i = 1:size(master_sg,1)
        master_data.(useArray).units(i).wf = wf{i};
    end
    clear wf;
end
clear data;
end

