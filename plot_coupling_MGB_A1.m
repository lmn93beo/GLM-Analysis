load ws_A1_MGB_all_stim_onlyB.mat

numMGB = numel(fieldnames(ws_all{1}));
numA1 = numel(ws_all);

MGBnames = fieldnames(ws_all{1});
id = 1;

for i = 1 : numA1
    for j = 1 : numMGB
        coupling = ws_all{i}.(MGBnames{j}).data;
        times = ws_all{i}.(MGBnames{j}).tr;
        subplot(numMGB, numA1, id);
        plot(times, exp(coupling), 'LineWidth', 2)
        
        if id ~= 65
            set(gca,'xtick',[])
            set(gca,'ytick',[])
        else
            xlabel('Time (ms)');
            ylabel('Exp (ws)');
            set(gca,'YAxisLocation','right');
        end
        
        xlim([0, 50]);
        ylim([0.5 2]);
        
        h = refline(0, 1);
        set(h, 'LineStyle', '--', 'Color', 'k');
        
        id = id + 1;
    end
end

%% Label A1 units
for id = 1:numA1
    subplot(numMGB, numA1, id);
    title(sprintf('A1Unit%d', A1_units(id)));
end

%% Label MGB units
for i = 1:numMGB
    id = numA1 * (i - 1) + 1;
    subplot(numMGB, numA1, id);
    ylabel(sprintf('MGB%d', MGB_units(i)));
end