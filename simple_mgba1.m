%% Load the data
clear all;

full_load = 0;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\STRFSessionsCombined';
file = [source '\2017-01-24_11-30-27_STRFMeanSpike.mat'];


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
        'unit1', 'unit12');
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
    MGB_units = numUnits(numUnits <= 9);
    A1_units = numUnits(numUnits > 9);
else
    % TODO: Hard coded for now...
    MGB_units = 1:19;
    A1_units = 20:29;
end

% Define start and end of 'trials'
end_id = 594; % decrease in time from t = 594 to 595...
ntrials = end_id - 1;
start_times = Lstim.start_time(1:ntrials) - 1;
end_times = Lstim.start_time(1:ntrials) + 2;
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

% Define trial structures
for trial = 1:end_id-1
    
    % Store duration in ms
    glmtrial(trial).duration = 1000 * (end_times(trial) - start_times(trial));
 
    % Store spike timings
    trial_vals1 = unit1(unit1 > start_times(trial) & unit1 < end_times(trial));
    trial_vals12 = unit12(unit12 > start_times(trial) & unit12 < end_times(trial));
    glmtrial(trial).MGBUnit1 = (trial_vals1 - start_times(trial)) * 1000;
    glmtrial(trial).A1Unit12 = (trial_vals12 - start_times(trial)) * 1000;
end

%% Build a GLM object
% Initialize the experiment with appropriate bin size
unitOfTime = 'ms';
binSize = 1;

param.samplingFreq = 1;
param.mouse = 'MouseX';

expt = buildGLM.initExperiment(unitOfTime, binSize, [], param);


expt = buildGLM.registerSpikeTrain(expt, 'A1Unit12', 'Our Neuron'); % Spike train!!!
expt = buildGLM.registerSpikeTrain(expt, 'MGBUnit1', 'Neighbor Neuron');



expt.trial = glmtrial;

%% Design specification
dspec = buildGLM.initDesignSpec(expt);
%binfun = @(t)(t==0)+ceil(t/binSize);

% Add covariate timing
%bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 25, binfun);
% bs = basisFactory.makeSmoothTemporalBasis('boxcar', 600, 16, binfun);
%offset = 100;
%dspec = buildGLM.addCovariateTiming(dspec, 'cueon', [], [], bs, offset);
%dspec = buildGLM.addCovariateBoxcar(dspec, 'sound', 'cueon', 'cueoff', 'Sound stim');


% Add coupling variable from MGB units
% for i = 1:1
%     disp(['Coupling from MGBUnit' num2str(MGB_units(i))]);
%     dspec = buildGLM.addCovariateSpiketrain(dspec, ...
%         ['MGBUnit' num2str(MGB_units(i))], ['MGBUnit' num2str(MGB_units(i))], ...
%         ['Coupling from MGB' num2str(MGB_units(i))]);
% end
dspec = buildGLM.addCovariateSpiketrain(dspec, 'coupling', 'MGBUnit1', 'Coupling from neuron 2');

% Self history
dspec = buildGLM.addCovariateSpiketrain(dspec, 'hist', 'A1Unit12', 'History filter');

% Design matrix
trialIndices = 1:end_id-1;
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);

%% Get dependent variable
y = buildGLM.getBinnedSpikeTrain(expt, 'A1Unit12', dm.trialIndices);

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
figure;

% Plot the result
for i = 1
    unit = dspec.covar(i).label;
    weight = ws.(unit).data;
    time = ws.(unit).tr;
    %plot(time,(ws.(unit).data));
    errorbar(ws.(unit).tr, ws.(unit).data, sqrt(wvar.(unit).data));
    hold on;
end

field_names = {dspec.covar.label};
legend(field_names(2:end));
% legend('MGBUnit1', 'MGBUnit2', 'MGBUnit3', 'MGBUnit4', ...
%     'MGBUnit5', 'MGBUnit6', 'MGBUnit7', 'MGBUnit8',...
%     'MGBUnit9', 'MGBUnit10', 'MGBUnit11', 'MGBUnit12',...
%     'MGBUnit13', 'MGBUnit14', 'MGBUnit15', 'MGBUnit16',...
%     'MGBUnit17', 'MGBUnit18', 'MGBUnit19');



