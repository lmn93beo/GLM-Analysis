%% Load the data
clear all;

full_load = 0;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\Session 1\';
file = [source '2017-01-23_15-36-21_STRFMeanSpike.mat'];


if full_load
    load(file);
else
    load(file, 'num_seq', 'Lstim', 'numUnits',...
        'unit1', 'unit2', 'unit3', 'unit4', ...
        'unit5', 'unit6', 'unit7', 'unit8', ...
        'unit9', 'unit10', 'unit11', 'unit12', ...
        'unit13', 'unit14', 'unit15', 'unit16', ...
        'unit17', 'unit18', 'unit19', 'unit20', ...
        'unit21', 'unit22', 'unit23', 'unit24', ...
        'unit25', 'unit26', 'unit27', 'unit28', ...
        'unit29');
end

%% Figuring out which unit is MGB/A1
% num_seq: number of the channels
% unit1, unit2,... encode spike times
% Lstim.start_time encodes start times, duration is always 1s
% numUnits: auditory-responsive units

MGB_channels = 1:9;
A1_channels = 10:18;

% Define which unit is MGB and which is A1
responsive_only = 1;

if responsive_only
    % unit is MGB or A1 && is auditory-responsive
    MGB_units = numUnits(numUnits <= 19);
    A1_units = numUnits(numUnits > 19);
else
    % TODO: Hard coded for now...
    MGB_units = 1:19;
    A1_units = 20:29;
end

% Define start and end of 'trials'
end_id = 594; % decrease in time from t = 594 to 595...
ntrials = end_id - 1;
start_times = Lstim.start_time(1:ntrials) + 2;
end_times = Lstim.start_time(1:ntrials) + 3;
deadline = Lstim.start_time(2:end_id);

% Add random jitter
jitter = rand(1, ntrials) * 0;
start_times = start_times - jitter;
end_times = end_times - jitter;

% Define edges
edges = [start_times, end_times];
edges = sort(edges);

%% Start compiling units
glmtrial = struct;

num_units = size(num_seq, 1);
unit_names = {};
var_names = {};

for unit_id = 1 : num_units
    % Add new unit if unit_id is in MGB_units or A1_units
    if ismember(unit_id, MGB_units) %MGB unit
        var_names = [var_names {['unit' num2str(unit_id)]}];
        unit_names = [unit_names {['MGBUnit' num2str(unit_id)]}];
    elseif ismember(unit_id, A1_units) %A1 unit
        var_names = [var_names {['unit' num2str(unit_id)]}];
        unit_names = [unit_names {['A1Unit' num2str(unit_id)]}];
    end
end

num_responsive = numel(MGB_units) + numel(A1_units);

% Define trial structures
for trial = 1 : end_id - 1
    
    % Store duration in ms
    glmtrial(trial).duration = 1000 * (end_times(trial) - start_times(trial));
    glmtrial(trial).cueon = 1000;
    glmtrial(trial).cueoff = 2000;
    
 
    % Store laser ON/OFF type
    glmtrial(trial).laser = mod(trial+1, 2);
    
    % Store spike timings
    for unit_id = 1 : num_responsive
        % Define variable name and unit name
        var_name = var_names{unit_id};
        unit_name = unit_names{unit_id};
    
        % Extract spikes of unit in the trial and put into glmtrial
        eval(['trial_vals = ' var_name '(' var_name '> start_times(trial) & '...
            var_name '< end_times(trial));']);
        
        % Some assertions just to make sure
        assert(all(trial_vals > start_times(trial)) && ...
            all(trial_vals < deadline(trial)));
        assert(all(trial_vals - start_times(trial) < glmtrial(trial).duration));
        
        glmtrial(trial).(unit_name) = ...
            1000 * (trial_vals - start_times(trial));
    end
end

%% Build a GLM object
% Initialize the experiment with appropriate bin size
unitOfTime = 'ms';
binSize = 1;

expt = buildGLM.initExperiment(unitOfTime, binSize);

% Add stimulus timing information
expt = buildGLM.registerTiming(expt, 'cueon', 'Cue Onset');
expt = buildGLM.registerTiming(expt, 'cueoff','Cue Offset');

% Laser ON/OFF condition
expt = buildGLM.registerValue(expt, 'laser', 'Laser ON/OFF');

% Spike trains
% MGB units
for i = 1:numel(MGB_units)
    disp(['MGBUnit' num2str(MGB_units(i))]);
    expt = buildGLM.registerSpikeTrain(expt, ...
        ['MGBUnit' num2str(MGB_units(i))], 'MGB Spike Train');
end

% A1 units
for i = 1:numel(A1_units)
    disp(['A1Unit' num2str(A1_units(i))]);
    expt = buildGLM.registerSpikeTrain(expt, ...
        ['A1Unit' num2str(A1_units(i))], 'A1 Spike Train');
end

expt.trial = glmtrial;

%% Design specification
dspec = buildGLM.initDesignSpec(expt);
binfun = @(t)(t==0)+ceil(t/binSize);

% Add coupling variable from MGB units
for i = 1:numel(MGB_units)
    disp(['Coupling from MGBUnit' num2str(MGB_units(i))]);
    dspec = buildGLM.addCovariateSpiketrain(dspec, ...
        ['MGBUnit' num2str(MGB_units(i))], ['MGBUnit' num2str(MGB_units(i))], ...
        ['Coupling from MGB' num2str(MGB_units(i))]);
end

% Design matrix
trialIndices = 1:2:ntrials;
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);

% Some processing on the design matrix
dm = buildGLM.removeConstantCols(dm);
dm = buildGLM.addBiasColumn(dm); % DO NOT ADD THE BIAS TERM IF USING GLMFIT



%% Perform the regrssion for each
ws_all = cell(1, numel(A1_units));
wvar_all = cell(1, numel(A1_units));


for A1Num = 1:numel(A1_units)
    A1Name = sprintf('A1Unit%d', A1_units(A1Num));
    fprintf('Working on A1 Unit %s...\n', A1Name)
    
    %% Get dependent variable
    y = buildGLM.getBinnedSpikeTrain(expt, A1Name, dm.trialIndices);

    %% Least squares for initialization
    tic
    wInit = dm.X' * dm.X \ dm.X' * y;
    toc

    %% Use matRegress for Poisson regression
    % it requires `fminunc` from MATLAB's optimization toolbox
    fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
    lfunc = @(w)(glms.neglog.poisson(w, dm.X, y, fnlin)); % cost/loss function

    opts = optimset('Algorithm', 'trust-region-reflective', ...
        'GradObj', 'on', 'Hessian','on');

    tic
    [wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
    wvar = diag(inv(hessian));
    toc

    %% Visualize the weights
    ws = buildGLM.combineWeights(dm, wml);
    wvar = buildGLM.combineWeights(dm, wvar);

    %%
    %figure;

    coupled_units = [1, 2, 3, 4, 7, 10, 14, 17, 18];
    ws_all{A1Num} = ws;
    wvar_all{A1Num} = wvar;
    

    % Plot the result
%     for i = 1:numel(dspec.covar)
%         unit = dspec.covar(i).label;
%         unitid = str2double(unit(8:end));
%         weight = ws.(unit).data;
%         time = ws.(unit).tr;
%         ws_all(i, A1Num) = ws.(unit).data;
%         hold on;
% 
%         if ismember(unitid, coupled_units)
%             plot(time,exp(ws.(unit).data), 'r');
%         else
%             plot(time, exp(ws.(unit).data), 'k');
%         end
        %errorbar(ws.(unit).tr, ws.(unit).data, sqrt(wvar.(unit).data));
%     end

%     field_names = {dspec.covar.label};
%     legend(field_names);
end

%% Save!
filename = 'ws_A1_MGB_all_baseline.mat';
if ~exist(filename, 'file')
    save(filename, 'ws_all', 'wvar_all', 'A1_units', 'MGB_units');
else
    error('Error: File exists');
end

