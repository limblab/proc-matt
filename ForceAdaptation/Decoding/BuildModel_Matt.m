function [filter, varargout]=BuildModel_Matt(binnedData, options)
%    [filter, varargout] = BuildModel(binnedData, options)
%
%
%       binnedData          : data structure to build model from
%       options             : structure with fields:
%           fillen              : filter length in seconds (tipically 0.5)
%           UseAllInputs        : 1 to use all inputs, 0 to specify a neuronID file, or a NeuronIDs array
%           PolynomialOrder     : order of the Weiner non-linearity (0=no Polynomial)
%           PredEMG, PredForce, PredCursPos, PredVeloc, numPCs :
%                               flags to include EMG, Force, Cursor Position
%                               and Velocity in the prediction model
%                               (0=no,1=yes), if numPCs is present, will
%                               use numPCs components as inputs instead of
%                               spikeratedata
%           Use_Thresh,Use_EMGs,Use_Ridge:
%                               options to fit only data above a certain
%                               threshold, use EMGs as inputs instead of
%                               spikes, or use a ridge regression to fit model
%           plotflag            : plot predictions after xval
%
%       Note on options: not all the fields have to be present in the
%       'option' structure provided in arguments. Those that are not will
%       be filled with the values from 'ModelBuildingDefault.m'
%
%       filter: structure of filter data (neuronIDs,H,P,emgguide,fillen,binsize)
%       varargout = {PredData}
%           [PredData]      : structure with EMG prediction data (fit)% 
%
%% Argument handling
   
    if ~isstruct(binnedData)
        binnedData = LoadDataStruct(binnedData);
    end

    binsize = double(binnedData.timeframe(2)-binnedData.timeframe(1));
    
    if nargout > 2
        disp('Wrong number of output arguments');
        return;
    end    
    
    % default values for options:
    default_options = ModelBuildingDefault();
    % fill other options as provided
    all_option_names = fieldnames(default_options);
    for i=1:numel(all_option_names)
        if ~isfield(options,all_option_names(i))
            options.(all_option_names{i}) = default_options.(all_option_names{i});
        end
    end
    
    if ~(options.PredEMGs || options.PredForce || options.PredCursPos || options.PredVeloc || options.PredTarg || options.PredCompVeloc || options.PredMoveDir)
        disp('No Outputs are Selected, Model Building Cancelled');
        filter = [];
        varargout = {};
        return;
    end


%% Inputs    
    
    %Need to be able to find which column(s) is the requested input(s) and only
    %use those to build the models.
    %Default is to use all the available inputs, otherwise ask for a list of
    %the ones you want to use.
    %desiredInputs are the columns in the firing rate matrix that are to be
    %used as inputs for the models
    
    if size(options.UseAllInputs,1)>1
        NeuronIDs = options.UseAllInputs;
        desiredInputs = get_desired_inputs(binnedData.spikeguide, neuronIDs);
    elseif options.UseAllInputs
%        disp('Using all available inputs')
        if ~isfield(binnedData,'neuronIDs')
            neuronIDs=spikeguide2neuronIDs(binnedData.spikeguide);
        else
            neuronIDs=binnedData.neuronIDs;
        end
        desiredInputs=1:size(neuronIDs,1);
    else
        if ~exist('NeuronIDsFile','var')
            [FileName, PathName] =uigetfile([dataPath '\NeuronIDfiles\' '*.mat'],'Filename of desired inputs? ');
            NeuronIDsFile = [PathName FileName];
        end
        neuronIDs = load(NeuronIDsFile);
        field_name = fieldnames(neuronIDs);
        neuronIDs = getfield(neuronIDs, field_name{:});
        desiredInputs = get_desired_inputs(binnedData.spikeguide, neuronIDs);
    end
    if isempty(desiredInputs)
        disp('Incompatible Data; Model Building Aborted');
        filter = [];
        if nargout > 1
            varargout(1) = {[]};
        end
        return;
    end

    numlags= round(options.fillen/binsize);%Designate the length of the filters/number of time lags
                                   % round helps getting rid of floating point error but care should
                                   % be taken in making sure fillen is a multiple of binsize.
    
    numsides=1;    %For a one-sided or causal filter

    %Select decoder inputs:
    if options.Use_EMGs
        Inputs = binnedData.emgdatabin;
        input_type = 'EMG';
    elseif options.numPCs
        disp('Using PCs...')
        Inputs = binnedData.spikeratedata(:,desiredInputs);
        [PCoeffs,Inputs] = princomp(zscore(Inputs));
        Inputs = Inputs(:,1:options.numPCs);
        input_type = 'princomp';
    else
        try
        Inputs = binnedData.spikeratedata(:,desiredInputs);
        input_type = 'spike';
        catch
            keyboard
        end
    end

%%
% DUPLICATE AND SHIFT FOR DISCONTINUITIES
%    Inputs = DuplicateAndShift(binnedData.spikeratedata(:,desiredInputs),numlags); numlags = 1;
    
%% Outputs
    
    Outputs = [];
    OutNames = [];
    
    %Decoder Outputs:
    if options.PredEMGs
       Outputs= [Outputs binnedData.emgdatabin];
       OutNames = [OutNames binnedData.emgguide];
    end
    if options.PredForce
        Outputs = [Outputs binnedData.forcedatabin];
        OutNames = [OutNames; binnedData.forcelabels];
    end
    if options.PredCursPos
        Outputs = [Outputs binnedData.cursorposbin];
        OutNames = [OutNames;  binnedData.cursorposlabels];
    end
    if options.PredVeloc
        Outputs = [Outputs binnedData.velocbin];
        OutNames = [OutNames;  binnedData.veloclabels];
    end
    if options.PredTarg
        Outputs = [Outputs binnedData.targetanglebin];
        OutNames = [OutNames; binnedData.targetanglelabels];
    end
    if options.PredCompVeloc
        % a kind of "compensated" velocity where force is accounted for
        Outputs = [Outputs binnedData.compvelocbin];
        OutNames = [OutNames; binnedData.compveloclabels];
    end
    if options.PredMoveDir
        % a kind of "compensated" velocity where force is accounted for
        Outputs = [Outputs binnedData.movedirbin];
        OutNames = [OutNames; binnedData.movedirlabels];
    end   
%% Calculate Filter

    %The following calculates the linear filters (H) that relate the inputs and outputs
    if options.Use_Ridge
        % Specify condition desired
        condition_desired = 10^4;
        % Duplicate and shift
        Inputs = DuplicateAndShift(Inputs,numlags); numlags = 1;
        % Train ridge model
        H = train_ridge(Inputs',Outputs',condition_desired);
    else
        [H,v,mcc]=filMIMO4(Inputs,Outputs,numlags,numsides,1);
%     H = MIMOCE1(Inputs,Outputs,numlags);
%     H = Inputs\Outputs;
%     Inputs = DuplicateAndShift(Inputs,numlags); numlags = 1;
    end
    
%% Then, add non-linearity if applicable

    fs=1; numsides=1;
    P=[]; T=[];
    patch = [];

    [PredictedData,spikeDataNew,ActualDataNew]=predMIMO4(Inputs,H,numsides,fs,Outputs);
        
    if options.PolynomialOrder
        %%%Find a Wiener Cascade Nonlinearity
        for z=1:size(PredictedData,2)
            if options.Use_Thresh
                %Find Threshold
                T_default = 1.25*std(PredictedData(:,z));
                [T(z,1), T(z,2), patch(z)] = findThresh(ActualDataNew(:,z),PredictedData(:,z),T_default);
                IncludedDataPoints = or(PredictedData(:,z)>=T(z,2),PredictedData(:,z)<=T(z,1));

                %Apply Threshold to linear predictions and Actual Data
                PredictedData_Thresh = PredictedData(IncludedDataPoints,z);
                ActualData_Thresh = ActualDataNew(IncludedDataPoints,z);

                %Replace thresholded data with patches consisting of 1/3 of the data to find the polynomial 
                Pred_patches = [ (patch(z)+(T(z,2)-T(z,1))/4)*ones(1,round(length(nonzeros(IncludedDataPoints))*4)) ...
                                 (patch(z)-(T(z,2)-T(z,1))/4)*ones(1,round(length(nonzeros(IncludedDataPoints))*4)) ];
                Act_patches = mean(ActualDataNew(~IncludedDataPoints,z)) * ones(1,length(Pred_patches));

                %Find Polynomial to Thresholded Data
                [P(:,z)] = WienerNonlinearity([PredictedData_Thresh; Pred_patches'], [ActualData_Thresh; Act_patches'], options.PolynomialOrder,'plot');
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%% Use only one of the following 2 lines:
                %
                %   1-Use the threshold only to find polynomial, but not in the model data
                T=[]; patch=[];                
                %
                %   2-Use the threshold both for the polynomial and to replace low predictions by the predefined value
%                 PredictedData(~IncludedDataPoints,z)= patch(z);
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                %Find and apply polynomial
                [P(:,z)] = WienerNonlinearity(PredictedData(:,z), ActualDataNew(:,z), options.PolynomialOrder);
            end
            PredictedData(:,z) = polyval(P(:,z),PredictedData(:,z));
        end
    end
    
%% Outputs

    filter = struct('neuronIDs', neuronIDs,...
                    'H', H,...
                    'P', P,...
                    'T',T,...
                    'patch',patch,...
                    'outnames', OutNames,...
                    'fillen',options.fillen,...
                    'binsize', binsize,...
                    'input_type',input_type);

    if options.numPCs
        filter.PC = PCoeffs(:,1:options.numPCs);
    end
    
    if nargout > 1
         PredData = struct('preddatabin', PredictedData, 'timeframe', ...
			 binnedData.timeframe(numlags:end),'spikeratedata',spikeDataNew, ...
			 'outnames',OutNames,'spikeguide',binnedData.spikeguide, ...
			 'vaf',RcoeffDet(PredictedData,ActualDataNew),'actualData',ActualDataNew);
        varargout{1} = PredData;
    end
    
    
end

function [Tinf, Tsup, patch] = findThresh(ActualData,LinPred,T)

    thresholding = 1;
    h = figure;
    xT = [0 length(LinPred)];
    offset = mean(LinPred)-mean(ActualData);
    LinPred = LinPred-offset;
    Tsup=mean(LinPred)+T;
    Tinf=mean(LinPred)-T;
    patch = mean(ActualData);
    
        while thresholding
            hold off; axis('auto');
            plot(ActualData,'b');
            hold on;
            plot(LinPred,'r');
            plot(xT,[Tsup Tsup],'k--',xT,[Tinf Tinf],'k--');
            legend('Actual Data', 'Linear Fit','Threshold');
            axis('manual');
            reply = input(sprintf('Redefine High Threshold? [%g] : ',Tsup));
            if ~isempty(reply)
                Tsup = reply;
            else
                thresholding=0;
            end
        end
        thresholding=1;
        while thresholding
            axis('auto');
            hold off;
            plot(ActualData,'b');
            hold on;
            plot(LinPred,'r');
            plot(xT,[Tsup Tsup],'k--',xT,[Tinf Tinf],'k--');
            legend('Actual Data', 'Linear Fit','Threshold');
            axis('manual');
            reply = input(sprintf('Redefine Low Threshold? [%g] : ',Tinf));
            if ~isempty(reply)
                Tinf = reply;
            else
                thresholding=0;
            end
        end
        thresholding=1;
        while thresholding
            axis('auto');
            hold off;
            plot(ActualData,'b');
            hold on;
            plot(LinPred,'r');
            plot(xT,[Tsup Tsup],'k--',xT,[Tinf Tinf],'k--', xT,[patch patch],'g');
            legend('Actual Data', 'Linear Fit','Threshold');
            axis('manual');
            reply = input(sprintf('Redefine Threshold Value? [%g] : ',patch));
            if ~isempty(reply)
                patch = reply;
            else
                thresholding=0;
            end
        end
        Tsup = Tsup+offset;
        Tinf = Tinf+offset;
        patch = patch+offset;
        
    close(h);
end
