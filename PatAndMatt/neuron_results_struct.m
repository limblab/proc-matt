function nr = neuron_results_struct(which_action,nr,analysis_info)
% Need a data structure to store analysis results for units
%
% It's an array of structs where each element is a neuron. Each field of
% that struct entry must be specified in analysis_info
%   General/meta fields:
%       date          : date neuron came from
%       monkey        : monkey neuron came from
%       task          : task neuron came from
%       perturbation  : perturbation neuron came from
%       array         : which array this neuron came from (M1, PMd, etc)
%       neuron_id     : neuron ID number for the session it's from
%       analysis_type : what type ('tuning' currently supported)
% 
%   Analysis_type specific fields:
%       tuning        : another array with an entry for any tuning methods.  Sub-fields:
%           name      : an arbitrary name for this analysis (e.g. 'onpeak' for tuning in a window from movement onset to peak speed)
%           notes     : string with info (can be used for documenting anything)
%           method    : 'regression','glm',etc
%           bootstrap : num_parameters x num_iterations
%           r_squared : num_neurons x num_itrations matrix of r-squareds
%           fr        : firing rates of that cell used for this tuning
%           theta     : directions used for this tuning
%
% INPUTS:
%   which_action      : string specifying how this function should behave
%       'add'         : if a cell exists in nr, appends analysis result appropriately
%       'replace'     : if a cell exists in nr and an analysis with the same name exists, it will replace it
%       'remove'      : if a cell exists in nr and an analysis with the same name exists, it will remove it
%   nr: an existing neuron results struct, or [] to initialize a new one
%   analysis_info: struct of parameters/info/results
%
% OUTPUTS:
%   nr: updated (or initialized) neuron results struct
%
% NOTES:
%   The way it's implemented now, when nr is [] or when a cell doesn't exist in nr yet,
% the value of which_action is basically meaningless. It always adds the entry.

date = analysis_info.date;
monkey = analysis_info.monkey;
task = analysis_info.task;
perturbation = analysis_info.perturbation;
array = analysis_info.array;
neuron_id = analysis_info.neuron_id;
analysis_type = analysis_info.analysis_type;
name = analysis_info.name;

%%
% partition out analysis information
switch lower(analysis_type)
    case 'tuning'
        if ~isfield(analysis_info,'notes')
            notes = [];
        else
            notes = analysis_info.notes;
        end
        sub_struct = struct('name',name, ...
            'notes',notes, ...
            'method',analysis_info.method, ...
            'bootstrap',analysis_info.bootstrap, ...
            'r_squared',analysis_info.r_squared, ...
            'fr',analysis_info.fr, ...
            'theta',analysis_info.theta);
        
    otherwise
        error(['Analysis type ' analysis_type ' is not recognized! Supported types are currently: 1)''tuning''.']);
end

%%
if ~isempty(nr) % if an existing structure was provided
    idx = find(strcmpi({nr.date},date) & strcmpi({nr.monkey},monkey) & strcmpi({nr.task},task) & strcmpi({nr.perturbation},perturbation) & strcmpi({nr.array},array) & ([nr.neuron_id] == neuron_id));
    
    if isempty(idx) % no neurons found
        if strcmpi(which_action,'remove')
            warning('No cell matched the information provided. Could not remove.');
        else
            disp('Cell not identified. Creating new entry.');
            nr_add = struct('date',date, ...
                'monkey',monkey, ...
                'task',task, ...
                'perturbation',perturbation, ...
                'array',array, ....
                'neuron_id',neuron_id, ...
                analysis_type,sub_struct);
            nr = [nr nr_add];
        end
    elseif length(idx) > 1 % multiple neurons found
        error('More than one neuron matched the information provided... Something is fishy here!');
    else % if one neuron matches
        analysis_idx = find(strcmpi({nr(idx).(analysis_type).name},sub_struct.name));
        
        switch lower(which_action)
            case 'remove'
                % check if the analysis name exists
                if isempty(analysis_idx) % if an analysis entry is not found
                    warning('Requested analysis name not found. Could not remove.');
                else
                    disp('Removing analysis entry for requested cell');
                    nr(idx).(analysis_type)(analysis_idx) = [];
                end
            case 'replace'
                if isempty(analysis_idx) % if an analysis entry is not found
                    warning('Requested analysis name not found. Adding entry instead of replacing.');
                    nr(idx).(analysis_type) = [nr(idx).(analysis_type) sub_struct];
                else
                    disp('Replacing analysis entry...');
                    nr(idx).(analysis_type)(analysis_idx) = sub_struct;
                end
            case 'add'
                if ~isempty(analysis_idx) % if an analysis entry is found
                    error('Requested analysis name already exists.');
                else
                    disp('Adding new analysis entry...');
                    nr(idx).(analysis_type) = [nr(idx).(analysis_type) sub_struct];
                end
        end
    end
else % if an empty structure was provided
    nr = struct('date',date, ...
                'monkey',monkey, ...
                'task',task, ...
                'perturbation',perturbation, ...
                'array',array, ....
                'neuron_id',neuron_id, ...
                analysis_type,sub_struct);
end

