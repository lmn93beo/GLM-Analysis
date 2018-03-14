%% Load the data
clear all;

full_load = 1;
source = 'Z:\CheetahData\Pavarrotti mouse\Pavarrotti Mouse\Session 1\';
file = [source '2017-01-23_15-36-21_STRFMeanSpike.mat'];


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

%% Get the right STRF (TODO: talk to Ian)
strfs = {STRFSmoCell1NL1, STRFSmoCell2NL1, STRFSmoCell3NL1, STRFSmoCell4NL1,...
    STRFSmoCell7NL1, STRFSmoCell8NL1,...
    STRFSmoCell9NL1, STRFSmoCell10NL1, STRFSmoCell12NL1,...
    STRFSmoCell13NL1, STRFSmoCell14NL1, ...
    STRFSmoCell17NL1, STRFSmoCell18NL1, STRFSmoCell20NL1,...
    STRFSmoCell22NL1, STRFSmoCell24NL1,...
    STRFSmoCell27NL1, STRFSmoCell28NL1};

figure;
for i = 1 : numel(numUnits)
    cellNum = numUnits(i);
    
    cmd = sprintf('strf1 = STRFRawCell%dNL1;', cellNum);
    eval(cmd);
    
    cmd = sprintf('strf2 = STRFRawCell%dNL2;', cellNum);
    eval(cmd);
    
    cmd = sprintf('strf3 = STRFRawCell%dNL3;', cellNum);
    eval(cmd);
    
    % Collapse over columns by summing
    strfCols1 = sum(strf1, 1);
    [~,posMax1] = max(strfCols1);
    
    strfCols2 = sum(strf2, 1);
    [~,posMax2] = max(strfCols2);
    
    strfCols3 = sum(strf3, 1);
    [~,posMax3] = max(strfCols3);

    % Plot the frequency tuning at this maximum
    freq_tuning1 = strf1(:, posMax1);
    freq_tuning2 = strf2(:, posMax2);
    freq_tuning3 = strf3(:, posMax3);
    
    subplot(9, 2, i);
    plot(freq_tuning1);
    hold on;
    plot(freq_tuning2);
    plot(freq_tuning3);
    
    if cellNum <= 19
        title(sprintf('MGBUnit%d', cellNum));
    else
        title(sprintf('A1Unit%d', cellNum));
    end
    %hold on;
end

