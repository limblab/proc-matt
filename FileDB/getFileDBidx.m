function idx = getFileDBidx(filedb,inclusions,exclusions)
% inclusions is cell like {'DBField','value', ...}
% exclusions is cell like {'string_of_logic', ...}
%
% Example exclusions:
%   exclusions = { ...
%       '~(ismember(filedb.Monkey,'Mihili') & datenum(filedb.Date) > datenum('2015-01-01'))',
%       'cellfun(@(x) all(ismember({'M1','PMd'},x)),filedb.Arrays)'};
rem_empty_fn = true;

%%% PROCESS THINGS TO INCLUDE
if rem(length(inclusions),2) ~= 0
    error('Inputs must be provided in pairs stating ...'' VARIABLE '',''VALUE'',...');
end

fn = inclusions(1:2:length(inclusions));
fv = inclusions(2:2:length(inclusions));

idx = ones(size(filedb,1),1);

% Check for the trials that match each criterion
for i = 1:length(fn)
        if ischar(fv{i})
            idx = idx & strcmpi(filedb.(fn{i}),fv{i});
        elseif iscell(fv{i})
            idx = idx & ismember(filedb.(fn{i}),fv{i});
        else
            idx = idx & ismember(filedb.(fn{i}),fv{i});
        end
end

%%% PROCESS THINGS TO EXCLUDE
if rem_empty_fn
    idx = idx & ~cellfun(@isempty,filedb.FileNames);
end

% this whole thingsis kinda hacky, but it makes the code very flexible
for i = 1:length(exclusions)
    idx = idx & eval(exclusions{i});
end

idx = find(idx);


