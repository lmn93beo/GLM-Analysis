%% Find props of all filters
load Filters\ws_A1_MGB_all_laser.mat

numMGB = numel(fieldnames(ws_all{1}));
numA1 = numel(ws_all);

MGBnames = fieldnames(ws_all{1});
id = 1;

props = cell(numA1, numMGB);
PTRs = zeros(numA1, numMGB);
Inhibs = zeros(numA1, numMGB);
Excs = zeros(numA1, numMGB);
Ratios = zeros(numA1, numMGB);
Diffs = zeros(numA1, numMGB);
PTTs = zeros(numA1, numMGB);

for i = 1 : numA1
    for j = 1 : numMGB
        % Get properties of the filter
        coupling = ws_all{i}.(MGBnames{j}).data;
        prop = extractProps(coupling);
        props{i, j} = prop;
        PTRs(i, j) = prop.PTR;
        Inhibs(i, j) = prop.InhibSubfield;
        Excs(i, j) = prop.ExcSubfield;
        Ratios(i, j) = prop.Ratio;
        Diffs(i, j) = prop.Diff;
        PTTs(i, j) = prop.PTT;

    end
end

%% Find frequency tuning curves
full_load = 1;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\Session 1\';
file = [source '2017-01-23_15-36-21_STRFMeanSpike.mat'];


if full_load
    load(file);
else
    load(file, 'num_seq', 'Lstim', 'numUnits',...
        'unit1', 'unit12');
end

strf_type = 'STRFSmo';

max_freq_lst = zeros(1, numel(numUnits));
%%
freq_tunings = cell(1, numel(numUnits));

for i = 1 : numel(numUnits)
    cellNum = numUnits(i);
    
    cmd = sprintf('strf1 = %sCell%dlaser1Post;', strf_type, cellNum);
    eval(cmd);
    
    cmd = sprintf('strf2 = %sCell%dlaser2Post;', strf_type, cellNum);
    eval(cmd);
    
    cmd = sprintf('strf3 = %sCell%dlaser3Post;', strf_type, cellNum);
    eval(cmd);

    % Plot the frequency tuning at this maximum
    [U1,~,~] = svd(strf1);
    [U2,~,~] = svd(strf2);
    [U3,~,~] = svd(strf3);
    
    freq_tuning1 = U1(:,1);
    freq_tuning2 = U2(:,1);
    freq_tuning3 = U3(:,1);
    freq_tuning_ave = freq_tuning1.^2 + freq_tuning2.^2 + freq_tuning3.^2;
    [~,freq_max] = max(freq_tuning_ave);
    
    freq_tunings{i} = freq_tuning_ave;
    
    subplot(9, 2, i);
    plot(freq_tuning_ave, 'LineWidth', 2);
    if cellNum <= 19
        title(sprintf('MGBUnit%d, freq = %d', cellNum, freq_max));
    else
        title(sprintf('A1Unit%d, freq = %d', cellNum, freq_max));
    end
    %hold on;
    
    h = refline(0, 0);
    set(h, 'LineStyle', '--', 'Color', 'k');
    
    max_freq_lst(i) = freq_max;
end

%% Find difference |MGB - A1|
MGB_freqs = max_freq_lst(1:13);
A1_freqs = max_freq_lst(14:end);
diff_freqs = abs(MGB_freqs - A1_freqs');

%% Plot diff against the statistics
figure;
subplot(2,3,1)
scatter(diff_freqs(:), PTRs(:));
title('PTR')

subplot(2,3,2)
scatter(diff_freqs(:), PTTs(:));
title('PTT')

subplot(2,3,3)
scatter(diff_freqs(:), Inhibs(:));
title('Inhib')

subplot(2,3,4)
scatter(diff_freqs(:), Excs(:));
title('Exc')

subplot(2,3,5)
scatter(diff_freqs(:), Ratios(:));
title('Ratio')

subplot(2,3,6)
scatter(diff_freqs(:), Diffs(:));
title('Diff')




%% Overlay plots
figure;
plot(freq_tunings{2});
hold on;
plot(freq_tunings{14});













