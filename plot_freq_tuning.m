%% Load the data
clear all;

full_load = 1;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\Session 1\';
file = [source '2017-01-23_15-36-21_STRFMeanSpike.mat'];


if full_load
    load(file);
else
    load(file, 'num_seq', 'Lstim', 'numUnits',...
        'unit1', 'unit12');
end

%% Get the right STRF (TODO: talk to Ian)
strfs = {STRFSmoCell1NL1, STRFSmoCell2NL1, STRFSmoCell3NL1, STRFSmoCell4NL1,...
    STRFSmoCell7NL1, STRFSmoCell8NL1,...
    STRFSmoCell9NL1, STRFSmoCell10NL1, STRFSmoCell12NL1,...
    STRFSmoCell13NL1, STRFSmoCell14NL1, ...
    STRFSmoCell17NL1, STRFSmoCell18NL1, STRFSmoCell20NL1,...
    STRFSmoCell22NL1, STRFSmoCell24NL1,...
    STRFSmoCell27NL1, STRFSmoCell28NL1};

figure;

strf_type = 'STRFSmo';

for i = 1 : numel(numUnits)
    cellNum = numUnits(i);
    
    cmd = sprintf('strf1 = %sCell%dNL1;', strf_type, cellNum);
    eval(cmd);
    
    cmd = sprintf('strf2 = %sCell%dNL2;', strf_type, cellNum);
    eval(cmd);
    
    cmd = sprintf('strf3 = %sCell%dNL3;', strf_type, cellNum);
    eval(cmd);
    
    % Collapse over columns by summing
    strfCols1 = sum(strf1, 1);
    [~,posMax1] = max(strfCols1);
    
    strfCols2 = sum(strf2, 1);
    [~,posMax2] = max(strfCols2);
    
    strfCols3 = sum(strf3, 1);
    [~,posMax3] = max(strfCols3);

    % Plot the frequency tuning at this maximum
%     freq_tuning1 = strf1(:, posMax1);
%     freq_tuning2 = strf2(:, posMax2);
%     freq_tuning3 = strf3(:, posMax3);
    [U1,~,~] = svd(strf1);
    [U2,~,~] = svd(strf2);
    [U3,~,~] = svd(strf3);
    
    freq_tuning1 = U1(:,1);
    freq_tuning2 = U2(:,1);
    freq_tuning3 = U3(:,1);


    
    subplot(9, 2, i);
    plot(freq_tuning1);
    hold on;
    plot(freq_tuning2);
    plot(freq_tuning3);
    %ylim([-5, 5]);
    
    if cellNum <= 19
        title(sprintf('MGBUnit%d', cellNum));
    else
        title(sprintf('A1Unit%d', cellNum));
    end
    %hold on;
    
    h = refline(0, 0);
    set(h, 'LineStyle', '--', 'Color', 'k');
    
end

