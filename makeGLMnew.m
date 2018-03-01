warning off;

dataName = 'WT19/2017-04-29/WT19_SpikeTimes.mat';
behavName = 'WT19/2017-04-29/WT19_SingleSession_Behavior.mat';
[glmtrial, unitOfTime, binSize, nTrials, binfun, numMD, numPFC, RejMD, RejPFC] = GLMpreprocess(dataName, behavName);

%%
for p = 25
    
    clearvars  -except AUC p glmtrial unitOfTime nTrials binSize binfun numMD numPFC
    
    Name = ['PFCUnit' num2str(p)];
    
    binSize = 1;
    
    binfun = @(t)(t==0)+ceil(t/binSize);
%     correct = arrayfun(@(x) glmtrial(x).reward,1:127);
%     glmtrial = glmtrial(find(correct==1));
    
    %% Specify the fields to load
    
    expt = buildGLM.initExperiment(unitOfTime, binSize);
    
    % values = context
    expt = buildGLM.registerValue(expt, 'context', 'Context');
    
    % timings = cue on, cue off
    expt = buildGLM.registerTiming(expt, 'cueon', 'Cue Onset');
    expt = buildGLM.registerTiming(expt, 'cueoff','Cue Offset');
    
    % timings = vision, rule
    % expt = buildGLM.registerTiming(expt, 'vision', 'Attend to Vision Rule');
    % expt = buildGLM.registerTiming(expt, 'audition', 'Attend to Audition Rule');
    
    % rule is value
    expt = buildGLM.registerTiming(expt, 'vision', 'Attend to Vision Rule');
    expt = buildGLM.registerTiming(expt, 'audition', 'Attend to Audition Rule');
    
    % timings = lpf, hpf
    expt = buildGLM.registerTiming(expt, 'lowpasscue', 'Low Pass Filter Cue');
    expt = buildGLM.registerTiming(expt, 'highpasscue', 'High Pass Filter Cue');
    
    % timing = decision
    expt = buildGLM.registerTiming(expt, 'choice', 'Choice Time');
    
    % type = rule / context pair
    expt = buildGLM.registerValue(expt, 'type', 'Trial type');
    
    % type v2
    expt = buildGLM.registerTiming(expt, 'R1C1', 'Rule 1, Context 1');
    expt = buildGLM.registerTiming(expt, 'R2C1', 'Rule 2, Context 1');
    expt = buildGLM.registerTiming(expt, 'R1C2', 'Rule 1, Context 2');
    expt = buildGLM.registerTiming(expt, 'R2C2', 'Rule 2, Context 2');
    
    % spike trains
    for i = 1:numPFC
        expt = buildGLM.registerSpikeTrain(expt, ['PFCUnit' num2str(i)], 'PFC Spike Train');
    end
    for i = 1:numMD
        expt = buildGLM.registerSpikeTrain(expt, ['MDUnit' num2str(i)], 'MD Spike Train');
    end
    
    %% build design spec object and put data in expt
    
    expt.trial = glmtrial;
    dspec = buildGLM.initDesignSpec(expt);
    
    %% trial type
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 600, 15, binfun);
    % bs = basisFactory.makeSmoothTemporalBasis('boxcar', 600, 16, binfun);
    offset = 0;
    dspec = buildGLM.addCovariateTiming(dspec, 'R1C1', [], [], bs, offset);
    dspec = buildGLM.addCovariateTiming(dspec, 'R2C1', [], [], bs, offset);
    dspec = buildGLM.addCovariateTiming(dspec, 'R1C2', [], [], bs, offset);
    dspec = buildGLM.addCovariateTiming(dspec, 'R2C2', [], [], bs, offset);
    
    %% spike history
    % bs = basisFactory.makeSmoothTemporalBasis('boxcar', 24, 12, expt.binfun);
    % dspec = buildGLM.addCovariateSpiketrain(dspec, 'hist', 'sptrain', 'History filter', bs);
    
    %% add MD coupling
    
    for coupleIdx = 1:numMD
        dspec = buildGLM.addCovariateSpiketrain(dspec, ['MDUnit' num2str(coupleIdx)], ['MDUnit' num2str(coupleIdx)], ['Coupling from MD Unit' num2str(coupleIdx)]);
    end
    
    %% build design matrix
    
    % trialIndices = 1:nTrials;
    trialIndices = 1:length(glmtrial);
    dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);
    
    endTrialIndices = cumsum(binfun([expt.trial(trialIndices).duration]));
    X = dm.X(1:endTrialIndices(3),:);
    mv = max(abs(X), [], 1); mv(isnan(mv)) = 1;
    X = bsxfun(@times, X, 1 ./ mv);
    
    %%
    
    y = buildGLM.getBinnedSpikeTrain(expt, Name,  dm.trialIndices);
    
    %% Do some processing on the design matrix
    dm = buildGLM.removeConstantCols(dm);
    % dm = buildGLM.zscoreDesignMatrix(dm, [colIndices{:}]);
    
    dm = buildGLM.addBiasColumn(dm); % DO NOT ADD THE BIAS TERM IF USING GLMFIT
    
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
    
    %% Visualize
    ws = buildGLM.combineWeights(dm, wml);
    wvar = buildGLM.combineWeights(dm, wvar);
    
    ct = 0;
    for kCov = 5 : 5+numMD-1;
        ct = ct+1;
        
        label = dspec.covar(kCov).label;
        figure(1);
        plot(ws.(label).tr, exp(ws.(label).data), 'k'); hold on
        xlim([ 0 50])
        
        drawnow;
        
        xax = find( ws.(label).tr < 50 );
        gain = (ws.(label).data);
        CouplingFilter{p}{ct} = gain;
        AUC{p}(ct) = sum( gain(xax) );
        Range{p}(ct) = range(gain(xax));
        
    end;
end;
