%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the root directory for this database
rootDir = '/Users/mgp046/Data/';
cerebusDataDir = 'CerebusData';
CDSDir = 'CDS';
TDDir = 'TrialDataFiles';
resultsDir = 'results';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(fullfile(rootDir,'filedb.mat'),'file') % initialize file database
    filedb = filedb_add();
    save(fullfile(rootDir,'filedb.mat'),'filedb');
else % load the file database
    load(fullfile(rootDir,'filedb.mat'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify which arrays each monkey had
%   {'MONKEY','END-DATE',{'ARRAYS'},{'paths to array maps'}}
array_list = {'MrT','2016-01-01',{'M1','PMd'}, ...
              {[rootDir 'MrT\MrT Right M1 SN 6250-0896.cmp'], ...
               [rootDir 'MrT\MrT Right PMd SN 6251-0880.cmp']}; ...
              'Chewie','2016-01-01',{'M1'}, ...
              {[rootDir 'Chewie\Chewie Right M1 SN 1025-0394.cmp']}; ...
              'Chewie','2020-01-01',{'M1','PMd'}, ...
              {[rootDir 'Chewie\Chewie Left M1 SN 6250-001474.cmp'], ...
               [rootDir 'Chewie\Chewie Left PMd SN 6251-001469.cmp']}; ...
              'Mihili','2016-01-01',{'M1','PMd'}, ...
              {[rootDir 'Mihili\Mihili Right M1 SN  6250-000989.cmp'], ...
               [rootDir 'Mihili\Mihili Right PMd SN  6251-000987.cmp']}; ...
              'Mihili','2020-01-01',{'M1','PMd'}, ...
              {[rootDir 'Mihili\Mihili Left M1 SN 1025-001452.cmp'], ...
               [rootDir 'Mihili\Mihili Left PMd SN 6251-001460.cmp']}; ...
             };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO: center out day
% RT: random target day
% FF: force field perturbation
% VR: visual rotation perturbation
% CS: control session
% BC: brain control session
% S(M/P): file has been sorted, M is M1, P is PMd
% ?: consider re-sorting...
sessionList = { ...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Mr T
    'MrT','2013-08-19','FF','CO','CCW',1.48,0.15, ...
        {}; ...   %1  S(M+P)
    'MrT','2013-08-20','FF','RT','CCW',1.48,0.15, ...
        {}; ...   %2  S(M+P)
    'MrT','2013-08-21','FF','CO','CCW',1.48,0.15, ...
        {'AD is split in two'}; ... %3
    'MrT','2013-08-22','FF','RT','CCW',1.48,0.15, ...
        {}; ...   %4  S(M+P)
    'MrT','2013-08-23','FF','CO','CCW',1.48,0.15, ...
        {}; ...   %5  S(M+P)
    'MrT','2013-08-30','FF','RT','CCW',1.48,0.15, ...
        {}; ...   %6  S(M+P)
    'MrT','2013-09-03','VR','CO','CCW',0.52,NaN, ...
        {}; ...   %7  S(M+P)
    'MrT','2013-09-04','VR','RT','CCW',0.52,NaN, ...
        {}; ...   %8  S(M+P)
    'MrT','2013-09-05','VR','CO','CCW',0.52,NaN, ...
        {}; ...   %9  S(M+P)
    'MrT','2013-09-06','VR','RT','CCW',0.52,NaN, ...
        {}; ...   %10 S(M+P)
    'MrT','2013-09-09','VR','CO','CCW',0.52,NaN, ...
        {}; ...   %11 S(M+P)
    'MrT','2013-09-10','VR','RT','CCW',0.52,NaN, ...
        {}; ...   %12 S(M+P)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Chewie
    'Chewie','2013-10-03','VR','CO','CCW',0.52,NaN, ...
        {}; ... %1  S(M)
    'Chewie','2013-10-09','VR','RT','CCW',0.52,NaN, ...
        {}; ... %2  S(M)
    'Chewie','2013-10-10','VR','RT','CCW',0.52,NaN, ...
        {}; ... %3  S(M) ?
    'Chewie','2013-10-11','VR','RT','CCW',0.52,NaN, ...
        {}; ... %4  S(M) ?
    'Chewie','2013-10-22','FF','CO','CCW',1.48,0.15, ...
        {'did second FF adaptation after washout'}; ... %5  S(M)
    'Chewie','2013-10-23','FF','CO','CCW',1.48,0.15, ...
        {}; ... %6  S(M)
    'Chewie','2013-10-28','FF','RT','CCW',1.48,0.15, ...
        {}; ... %7  S(M) ?
    'Chewie','2013-10-29','FF','RT','CCW',1.48,0.15, ...
        {}; ... %8  S(M) ?
    'Chewie','2013-10-31','FF','CO','CCW',1.48,0.15, ...
        {}; ... %9  S(M)
    'Chewie','2013-11-01','FF','CO','CCW',1.48,0.15, ...
        {}; ... %10 S(M)
    'Chewie','2013-12-03','FF','CO','CCW',1.48,0.15, ...
        {}; ... %11 S(M)
    'Chewie','2013-12-04','FF','CO','CCW',1.48,0.15, ...
        {}; ... %12 S(M)
    'Chewie','2013-12-09','FF','RT','CCW',1.48,0.15, ...
        {}; ... %13 S(M) ?
    'Chewie','2013-12-10','FF','RT','CCW',1.48,0.15, ...
        {}; ... %14 S(M) ?
    'Chewie','2013-12-12','VR','RT','CCW',0.52,NaN, ...
        {}; ... %15 S(M) ?
    'Chewie','2013-12-13','VR','RT','CCW',0.52,NaN, ...
        {}; ... %16 S(M) ?
    'Chewie','2013-12-17','FF','RT','CCW',1.48,0.15, ...
        {}; ... %17 S(M) ?
    'Chewie','2013-12-18','FF','RT','CCW',1.48,0.15, ...
        {}; ... %18 S(M) ?
    'Chewie','2013-12-19','VR','CO','CCW',0.52,NaN, ...
        {}; ... %19 S(M) ?
    'Chewie','2013-12-20','VR','CO','CCW',0.52,NaN, ...
        {}; ... %20 S(M) ?
    'Chewie','2015-03-09','CS','CO','',NaN,NaN, ...
        {}; ... %21 S(M) ?
    'Chewie','2015-03-11','CS','CO','',NaN,NaN, ...
        {}; ... %22 S(M) ?
    'Chewie','2015-03-12','CS','CO','',NaN,NaN, ...
        {}; ... %23 S(M) ?
    'Chewie','2015-03-13','CS','CO','',NaN,NaN, ...
        {}; ... %24 S(M) ?
    'Chewie','2015-03-16','CS','RT','',NaN,NaN, ...
        {}; ... %25 S(M) ?
    'Chewie','2015-03-17','CS','RT','',NaN,NaN, ...
        {}; ... %26 S(M) ?
    'Chewie','2015-03-18','CS','RT','',NaN,NaN, ...
        {}; ... %27 S(M) ?
    'Chewie','2015-03-19','CS','CO','',NaN,NaN, ...
        {}; ... %28 S(M) ?
    'Chewie','2015-03-20','CS','RT','',NaN,NaN, ...
        {}; ... %29 S(M) ?
    'Chewie','2015-06-29','FF','CO','CW',1.48,0.15, ...
        {}; ... %30 S(M) SHORT WASHOUT
    'Chewie','2015-06-30','FF','CO','CW',1.48,0.15, ...
        {}; ... %31 S(M)
    'Chewie','2015-07-01','FF','CO','CCW',1.48,0.15, ...
        {}; ... %32 S(M)
    'Chewie','2015-07-03','FF','CO','CCW',1.48,0.15, ...
        {}; ... %33 S(M)
    'Chewie','2015-07-06','FF','CO','CW',1.48,0.15, ...
        {}; ... %34 S(M)
    'Chewie','2015-07-07','FF','CO','CW',1.48,0.15, ...
        {}; ... %35 S(M)
    'Chewie','2015-07-08','FF','CO','CW',1.48,0.15, ...
        {}; ... %36 S(M)
    'Chewie','2015-07-09','VR','CO','CW',0.52,NaN, ...
        {}; ... %37 S(M) ?
    'Chewie','2015-07-10','VR','CO','CW',0.52,NaN, ...
        {}; ... %38 S(M) ?
    'Chewie','2015-07-13','VR','CO','CW',0.52,NaN, ...
        {}; ... %39 S(M) ?
    'Chewie','2015-07-14','VR','CO','CCW',0.52,NaN, ...
        {}; ... %40 S(M) ?
    'Chewie','2015-07-15','VR','CO','CCW',0.52,NaN, ...
        {}; ... %41 S(M) ?
    'Chewie','2015-07-16','VR','CO','CCW',0.52,NaN, ...
        {}; ... %42 S(M) ?
    'Chewie','2015-11-09','GR','CO','CCW',0.52,120, ...
        {}; ... %43 S(M)
    'Chewie','2015-11-10','GR','CO','CCW',0.52,120, ...
        {}; ... %44 S(M)
    'Chewie','2015-11-12','GR','CO','CW',0.52,120, ...
        {}; ... %45 S(M)
    'Chewie','2015-11-13','GR','CO','CCW',0.52,240, ...
        {}; ... %46 S(M)
    'Chewie','2015-11-16','GR','CO','CCW',0.52,240, ...
        {}; ... %47 S(M)
    'Chewie','2015-11-17','GR','CO','CCW',0.52,240, ...
        {}; ... %48 S(M)
    'Chewie','2015-11-18','GR','CO','CCW',0.52,240, ...
        {'smaller target size'}; ... %49 S(M)
    'Chewie','2015-11-19','VR','CO','CCW',1.05,NaN, ...
        {}; ... %50 S(M)
    'Chewie','2015-11-20','GR','CO','CCW',0.52,400, ...
        {}; ... %51 S(M)
    'Chewie','2015-12-01','VR','CO','CCW',0.785,NaN, ...
        {}; ... %52 S(M)
    'Chewie','2015-12-03','VR','CO','CCW',1.05,NaN, ...
        {}; ... %53 S(M)
    'Chewie','2015-12-04','VR','CO','CCW',0.785,NaN, ...
        {}; ... %54 S(M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Chewie (M1 and PMd)
    'Chewie','2016-09-09','VR','CO','CW',0.52,NaN, ...
        {}; ... %1 S(M+P)
    'Chewie','2016-09-12','VR','CO','CCW',0.52,NaN, ...
        {}; ... %2 S(M+P)
    'Chewie','2016-09-14','VR','CO','CW',0.52,NaN, ...
        {}; ... %4 S(M+P)
    'Chewie','2016-09-15','FF','CO','CW',1.48,0.05, ...
        {}; ... %5 S(M+P)
    'Chewie','2016-09-16','VR','CO','CCW',0.52,NaN, ...
        {}; ... %6
    'Chewie','2016-09-19','FF','CO','CCW',1.48,0.05, ...
        {}; ... %7 S(M+P)
    'Chewie','2016-09-20','VR','CO','CW',0.52,NaN, ...
        {}; ... %8
    'Chewie','2016-09-21','FF','CO','CW',1.48,0.05, ...
        {}; ... %9
    'Chewie','2016-09-22','VR','CO','CCW',0.52,NaN, ...
        {}; ... %10
    'Chewie','2016-09-23','FF','CO','CCW',1.48,0.05, ...
        {}; ... %11
    'Chewie','2016-09-27','VR','CO','CW',0.52,NaN, ...
        {'First day of fixed Lab 6 Jacobian'}; ... %12
    'Chewie','2016-09-28','FF','CO','CCW',1.48,0.15, ...
        {'BSOD, Washout is two files'}; ... %13
    'Chewie','2016-09-29','VR','CO','CW',0.52,NaN, ...
        {}; ... %14
    'Chewie','2016-09-30','FF','CO','CCW',1.48,0.15, ...
        {'M1 paused during recording'}; ... %15
    'Chewie','2016-10-05','FF','CO','CCW',1.48,0.15, ...
        {}; ... %16 S(M+P)
    'Chewie','2016-10-06','VR','CO','CW',0.52,NaN, ...
        {}; ... %17 S(M+P)
    'Chewie','2016-10-07','FF','CO','CW',1.48,0.15, ...
        {}; ... %18 S(M+P)
    'Chewie','2016-10-11','FF','CO','CW',1.48,0.15, ...
        {}; ... %19 S(M+P)
    'Chewie','2016-10-12','VR','CO','CCW',0.52,NaN, ...
        {}; ... %20 
    'Chewie','2016-10-13','FF','CO','CW',1.48,0.15, ...
        {}; ... %21  GREAT SESSION
    'Chewie','2016-10-14','CS','CO','',NaN,NaN, ...
        {}; ... %22
    'Chewie','2016-10-17','CS','CO','',NaN,NaN, ...
        {}; ... %23
    'Chewie','2016-10-18','CS','CO','',NaN,NaN, ...
        {}; ... %24
    'Chewie','2016-10-19','CS','CO','',NaN,NaN, ...
        {}; ... %25
    'Chewie','2016-10-20','CS','CO','',NaN,NaN, ...
        {'CO and RT day'}; ... %26
    'Chewie','2016-10-20','CS','RT','',NaN,NaN, ...
        {'CO and RT day'}; ... %26
    'Chewie','2016-10-21','CS','CO','',NaN,NaN, ...
        {'CO and RT day'}; ... %27
    'Chewie','2016-10-21','CS','RT','',NaN,NaN, ...
        {'CO and RT day'}; ... %27
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Mihili
    'Mihili','2014-01-14','VR','RT','CCW',0.52,NaN, ...
        {}; ...    %1  S(M+P) ?
    'Mihili','2014-01-15','VR','RT','CCW',0.52,NaN, ...
        {}; ...    %2  S(M+P) ?
    'Mihili','2014-01-16','VR','RT','CCW',0.52,NaN, ...
        {}; ...    %3  S(M+P) ?
    'Mihili','2014-02-03','FF','CO','CCW',1.48,0.15, ...
        {'didnt complete reaches to a target'}; ...    %4  S(M+P) ?
    'Mihili','2014-02-14','FF','RT','CCW',1.48,0.15, ...
        {}; ...    %5  S(M+P) ?
    'Mihili','2014-02-17','FF','CO','CCW',1.48,0.15, ...
        {}; ...    %6  S(M+P) P?
    'Mihili','2014-02-18','FF','CO','CCW',1.48,0.15, ...
        {'did both perturbations'}; ...    %7  S(M+P) P?
    'Mihili','2014-02-21','FF','RT','CCW',1.48,0.15, ...
        {}; ...    %9  S(M+P) ?
    'Mihili','2014-02-24','FF','RT','CCW',1.48,0.15, ...
        {'did both perturbations'}; ...    %10 S(M+P) ?
    'Mihili','2014-03-03','VR','CO','CCW',0.52,NaN, ...
        {}; ...    %12 S(M+P) ?
    'Mihili','2014-03-04','VR','CO','CCW',0.52,NaN, ...
        {}; ...    %13 S(M+P) ?
    'Mihili','2014-03-06','VR','CO','CCW',0.52,NaN, ...
        {}; ...    %14 S(M+P) ?
    'Mihili','2014-03-07','FF','CO','CCW',1.48,0.15, ...
        {}; ...    %15 S(M+P) P?
    'Mihili','2014-06-26','CS','CO','',NaN,NaN, ...
        {}; ...    %16 S(M+P) ?
    'Mihili','2014-06-27','CS','CO','',NaN,NaN, ...
        {}; ...    %17 S(M+P) ?
    'Mihili','2014-09-29','CS','CO','',NaN,NaN, ...
        {}; ...    %18 S(M+P) ?
%     'Mihili','2014-12-03','CS','CO','',NaN,NaN, ...
%         {}; ...    %19 S(M)   ?
    'Mihili','2015-05-11','CS','CO','',NaN,NaN, ...
        {}; ...    %21 S(M+P) ?
    'Mihili','2015-05-12','CS','CO','',NaN,NaN, ...
        {}; ...    %22 S(M+P) ?
    'Mihili','2015-06-10','FF','CO','CW',1.48,0.15, ...
        {}; ...    %23 S(M+P) P? SHORT WASHOUT
    'Mihili','2015-06-11','FF','CO','CCW',1.48,0.15, ...
        {}; ...    %24 S(M+P) P? SHORT WASHOUT
    'Mihili','2015-06-15','FF','CO','CW',1.48,0.15, ...
        {'Something may be weird in PMd population?'}; ...  %26 S(M+P) P?
    'Mihili','2015-06-16','FF','CO','CW',1.48,0.15, ...
        {}; ...    %27 S(M+P) P?
    'Mihili','2015-06-17','FF','CO','CCW',1.48,0.15, ...
        {'Didnt complete reaches to a target'}; ...  %28 S(M+P) ?
    'Mihili','2015-06-23','VR','CO','CW',0.52,NaN, ...
        {}; ...    %29 S(M+P) ?
    'Mihili','2015-06-25','VR','CO','CW',0.52,NaN, ...
        {}; ...    %30 S(M+P) ?
    'Mihili','2015-06-26','VR','CO','CW',0.52,NaN, ...
        {}; ...    %31 S(M+P) ?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    }; %end

% SORT MIHILI 2014-12-03 PMd file
% RECONSIDER ANALYZING
%   'Mihili','2014-01-17','VR','RT'; ... % Poor work ethic in washout, so it's really long... might be useable if needed
%   'Mihili','2015-06-12','FF','CO'; ...    %25 S(M-P) ? - SHORT WASHOUT, and a fairly garbage session because he didn't complete reaches to a target


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad and questionable data
% %%% Brain control from first implant
%     'Mihili','2015-07-03','BC','CO'; ...    %32 S(M-P) not the best session
%     'Mihili','2015-07-06','BC','CO'; ...    %33 S(M-P) not the best session
%     'Mihili','2015-07-08','BC','CO'; ...    %34 S(-P)
%     'Mihili','2015-07-10','BC','CO'; ...    %35 S(-P)
%     'Mihili','2015-07-13','BC','CO'; ...    %36
%     'Mihili','2015-07-14','BC','CO'; ...    %37
%     'Mihili','2015-07-15','BC','CO'; ...    %38
%     'Mihili','2015-07-16','BC','CO'; ...    %39
% 
% %%% Either bad adaptation or incomplete data. Might be useful for BL though.
% 'Mihili','2013-12-02','FF','RT'; ... S(M-P)
% 'Mihili','2013-12-07','FF','RT'; ... S(M-P)
% 'Mihili','2013-12-08','FF','RT'; ... S(M-P)
% 'Mihili','2013-12-12','FF','CO'; ...
% 'Mihili','2015-12-11','FF','CO'; ...    %20 S(M-P) ? - SOMETHING SEEMS STREANGE IN BEHAVIOR. NO WASHOUT
% 
% 'Mihili','2014-01-17','VR','RT'; ... % Poor work ethic in washout, so it's really long... might be useable if needed
% 'Mihili','2014-01-20','VR','CO'; ... % Adaptation period is half as short as it should be. Maybe not a total waste though
% 'Mihili','2015-06-12','FF','CO'; ... % short washout
% 
% 'Mihili','2014-01-21','VR','CO'; ... % good behavior but busted PMd headstage. M1 should be good
% 'Mihili','2014-01-22','VR','CO'; ... % good behavior but busted PMd headstage. M1 should be good
% 'Mihili','2014-01-28','FF','RT'; ... % good behavior but busted PMd headstage. M1 should be good
% 'Mihili','2014-01-29','FF','RT'; ... % good behavior but busted PMd headstage. M1 should be good
% %%% ALSO A FEW PARTIAL DAYS WITH NO WASHOUT
% mihili_partialDays = {'Mihili','2014-02-11','FF','CO'; ...
%     'Mihili','2014-02-22','FF','RT'; ...
%     'Mihili','2014-02-25','VR','CO'; ...
%     'Mihili','2014-02-27','VR','CO'; ...
%     'Mihili','2014-02-28','VR','CO'; ...
%     'Mihili','2015-06-24','VR','CO'; ...
% note: 12/12 and 12/13 were the first two FF days for
% Mihili... might be worth looking at, even though they
% don't have complete data
% 
% mrt_ff_iffyDays = {'2013-08-27', ... %   RT - poor work ethic, but might be good?
%                    '2013-08-28'};    %   CO - poor work ethic, adaptation period is a bit short
% mrt_ff_badDays = {'2013-08-13', ...  % S CO - force field at 90 degrees
%                    '2013-08-14'};     % S RT - force field changed partway through
% mrt_different_but_good_days = 'MrT','2013-09-24','VRFF','RT'; ... %13 S(P) ? x - 30 CCW
%     'MrT','2013-09-25','VRFF','RT'; ... %14 S(P) ? x - 30 CCW
%     'MrT','2013-09-27','VRFF','RT'; ... %15 S(P) ? x - 30 CCW
%     'MrT','2013-10-11','VR','RT'; ...   %16 S(P) ? x - 45 CCW
%
% chewie_iffyDays = {'2013-10-17', ... % S ? RT 0.13 force mag
%                    '2013-10-18'};    %   RT 0.1 force mag
%
% chewie_m1pmd_baddays
%     'Chewie','2016-09-07','VR','CO'; ... %1  S(M+P) - 30 CCW - MAY HAVE A LOT OF PMD SHUNTING DO NOT USE 
%     'Chewie','2016-09-13','FF','CO','CCW',1.48,0.15; ... %3 S(M+P) - 0.05 1.48 CCW - Had amplifier problems in PMd so don't trust PMd neurons
%
% jaco_stim_recording_days
%     'Jaco','2016-01-27','CS','CO'; ...
%     'Jaco','2016-01-28','CS','CO'; ...
%     'Jaco','2016-01-29','CS','CO'; ...
%     'Jaco','2016-02-02','CS','CO'; ...
%     'Jaco','2016-02-03','CS','CO'; ...
%     'Jaco','2016-02-04','CS','CO'; ...
%       'Jaco','2016-02-15','CS','CO'; ...
%       'Jaco','2016-02-17','CS','CO'; ...
% 'Jaco','2016-02-18','CS','CO'; ...
%    %%% Jaco CF Days
%     'Jaco','2016-04-05','FF','CO'; ... %1 S(M) - 0.15 1.48 CCW
%     'Jaco','2016-04-06','FF','CO'; ... %2 S(M) - 0.15 1.48 CCW
%     'Jaco','2016-04-07','FF','CO'; ... %3 S(M) - 0.15 1.48 CW - short-ish washout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%