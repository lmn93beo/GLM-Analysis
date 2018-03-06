%% Load the data
clear all;

full_load = 0;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\STRFSessionsCombined';
file = [source '\2017-01-23_15-36-21_STRFMeanSpike.mat'];
file_short = [source '\2017-01-23_15-36-21_STRFMeanSpike_short.mat'];


if full_load
    load(file);

    % Save only the variables we need
%     save(file_short, 'num_seq', 'Lstim', 'unit1', 'unit2', 'unit3', 'unit4', ...
%         'unit5', 'unit6', 'unit7', 'unit8', ...
%         'unit9', 'unit10', 'unit11', 'unit12', ...
%         'unit13', 'unit14', 'unit15', 'unit16', ...
%         'unit17', 'unit18', 'unit19', 'unit20', ...
%         'unit21', 'unit22', 'unit23', 'unit24', ...
%         'unit25', 'unit26', 'unit27', 'unit28', ...
%         'unit29');
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
MGB_units = [];
A1_units = [];

for i = 1:size(num_seq, 1)
    % unit is MGB && is auditory-responsive
    MGB_units = numUnits(numUnits <= 19);
    A1_units = numUnits(numUnits > 19);
end

% Define start and end of 'trials'
end_id = 594; % decrease in time from t = 594 to 595...
ntrials = end_id - 1;
start_times = Lstim.start_time(1:ntrials);
end_times = Lstim.start_time(1:ntrials)+1;
deadline = Lstim.start_time(2:end_id);

% Add random jitter
jitter = 0;
start_times = start_times + jitter * rand(1, ntrials);
end_times = end_times - jitter * rand(1, ntrials);

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

num_responsive = numel(numUnits);

% Define trial structures
for trial = 1 : end_id - 1
    
    % Store duration in ms
    glmtrial(trial).duration = 1000 * (end_times(trial) - start_times(trial));
    
    % Store spike timings
    for unit_id = 1 : num_responsive
        if trial == 4 && unit_id == 1
            disp('here');
        end
        % Define vaSriable name and unit name
        var_name = var_names{unit_id};
        unit_name = unit_names{unit_id};
    
        % Extract spikes of unit in the trial and put into glmtrial
        eval(['trial_vals = ' var_name '(' var_name '> start_times(trial) & '...
            var_name '< end_times(trial));']);
        
%         eval(['groups = discretize(', var_name, ', edges);']);
% 
%         eval(['trial_vals = ', var_name, '(groups == trial * 2 - 1);']);
        %trial_vals = unit1(groups == trial * 2 - 1)';
        
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

% Spike trains
% spike trains
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

% Add coupling variable from MGB units
for i = 1:numel(MGB_units)
    disp(['Coupling from MGBUnit' num2str(MGB_units(i))]);
    dspec = buildGLM.addCovariateSpiketrain(dspec, ...
        ['MGBUnit' num2str(MGB_units(i))], ['MGBUnit' num2str(MGB_units(i))], ...
        ['Coupling from MGB' num2str(MGB_units(i))]);
end

% Design matrix
trialIndices = 1:end_id - 1;
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);

%% Get dependent variable
y = buildGLM.getBinnedSpikeTrain(expt, 'A1Unit20', dm.trialIndices);

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

figure;

% Plot the result
for i = 1:numel(MGB_units)
    unit = ['MGBUnit' num2str(MGB_units(i))];
    weight = ws.(unit).data;
    time = ws.(unit).tr;
    plot(time,exp(ws.(unit).data));
    hold on;
end

legend('MGBUnit1', 'MGBUnit2', 'MGBUnit3', 'MGBUnit4', ...
    'MGBUnit5', 'MGBUnit6', 'MGBUnit7', 'MGBUnit8',...
    'MGBUnit9', 'MGBUnit10', 'MGBUnit11', 'MGBUnit12',...
    'MGBUnit13', 'MGBUnit14', 'MGBUnit15', 'MGBUnit16',...
    'MGBUnit17', 'MGBUnit18', 'MGBUnit19');

