%% Load the data
clear all;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\STRFSessionsCombined';
file = [source '\2017-01-23_15-36-21_STRFMeanSpike.mat'];
load(file);

%% Figuring out which unit is MGB/A1
% num_seq: number of the channels
% unit1, unit2,... encode spike times
% Lstim.start_time encodes start times, duration is always 1s

MGB_channels = 1:9;
A1_channels = 10:18;

MGB_units = [];
A1_units = [];


for i = 1:size(num_seq, 1)
    if ~isempty(find(MGB_channels == num_seq(i, 1), 1))
        MGB_units = [MGB_units, i];
    elseif ~isempty(find(A1_channels == num_seq(i, 1), 1))
        A1_units = [A1_units, i];
    end
end

% Define start and end of 'trials'
end_id = 594; % decrease in time from t = 594 to 595...
start_times = Lstim.start_time(1:end_id-1);
end_times = Lstim.start_time(2:end_id);

% Add random jitter
jitter = 1;
start_times = start_times + 1.1 + rand;
end_times = end_times - 0.1 - rand;

% Group channels into trials
edges = [start_times, end_times];
edges = sort(edges);
groups = discretize(unit1, edges);

%% Start compiling units
glmtrial = struct;

num_units = size(num_seq, 1);
unit_names = cell(1, num_units);
var_names = cell(1, num_units);

for unit_id = 1 : num_units
    var_names{unit_id} = ['unit' num2str(unit_id)];
    if ~isempty(find(MGB_channels == num_seq(unit_id, 1), 1)) %MGB unit
        unit_names{unit_id} = ['MGBUnit' num2str(unit_id)];
    elseif ~isempty(find(A1_channels == num_seq(unit_id, 1), 1)) %A1 unit
        unit_names{unit_id} = ['A1Unit' num2str(unit_id)];
    end
end

% Define trial structures
for trial = 1 : end_id - 1
    % Store duration in ms
    glmtrial(trial).duration = 1000 * (end_times(trial) - start_times(trial));
    
    % Store spike timings
    for unit_id = 1 : num_units
        % Define variable name and unit name
        var_name = var_names{unit_id};
        unit_name = unit_names{unit_id};
    
        % Extract spikes of unit in the trial and put into glmtrial
        eval(['groups = discretize(', var_name, ', edges);']);

        eval(['trial_vals = ', var_name, '(groups == trial * 2 - 1);']);
        %trial_vals = unit1(groups == trial * 2 - 1)';
        assert(all(trial_vals > start_times(trial)) && ...
            all(trial_vals < end_times(trial)));
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
for i = 1:num_units
    if ~isempty(find(MGB_channels == num_seq(i, 1), 1)) %MGB unit
        expt = buildGLM.registerSpikeTrain(expt, ...
            ['MGBUnit' num2str(i)], 'MGB Spike Train');
    elseif ~isempty(find(A1_channels == num_seq(i, 1), 1)) %A1 unit
        expt = buildGLM.registerSpikeTrain(expt, ...
            ['A1Unit' num2str(i)], 'A1 Spike Train');
    end
end

expt.trial = glmtrial;

%% Design specification
dspec = buildGLM.initDesignSpec(expt);

% Add coupling variable from A1 units
for i = 1:num_units
    if ~isempty(find(MGB_channels == num_seq(i, 1), 1))
        dspec = buildGLM.addCovariateSpiketrain(dspec, ...
            ['MGBUnit' num2str(i)], ['MGBUnit' num2str(i)], ...
            ['Coupling from MGB' num2str(i)]);
    end
end

% Design matrix
trialIndices = 1:end_id - 1;
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);




%% Get dependent variable
y = buildGLM.getBinnedSpikeTrain(expt, 'A1Unit21', dm.trialIndices);

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

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
wvar = diag(inv(hessian));

%% Visualize the weights
ws = buildGLM.combineWeights(dm, wml);

figure;

% Plot the result
for i = 1:19
    unit = ['MGBUnit' num2str(i)];
    weight = ws.(unit).data;
    time = ws.(unit).tr;
    plot(time,exp(ws.(unit).data));
    hold on;
end

