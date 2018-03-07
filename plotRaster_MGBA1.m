function raster_cell = plotRaster_MGBA1(glmstruct, trial_num, start_time, end_time)
% trial_num: number indicating trial

trial_struct = glmstruct(trial_num);

% Make a cell that contains all units
field_names = fieldnames(trial_struct);
raster_cell = cell(1, numel(field_names) - 1);

for i = 2:numel(field_names)
    raster_cell{i - 1} = trial_struct.(field_names{i});
end

plotRaster(raster_cell);
yticks((1:29) + 0.5);
yticklabels(field_names(2:end));

% Highlight stimulus period
patch([start_time end_time end_time start_time], [0 0 numel(field_names) numel(field_names)], [0.8, 0.8, 0.8], 'FaceAlpha',.3);

end