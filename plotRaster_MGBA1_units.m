function raster_cell = plotRaster_MGBA1_units(glmstruct, unit_num, start_time, end_time)
% trial_num: number indicating trial

% Make a cell that contains all units
field_names = fieldnames(glmstruct);
raster_cell = {glmstruct.(field_names{unit_num})};
raster_cell = raster_cell(1:2:200);

plotRaster(raster_cell);
%yticks((1:29) + 0.5);
%yticklabels(field_names(2:end));

% Highlight stimulus period
patch([start_time end_time end_time start_time], [0 0 numel(raster_cell) numel(raster_cell)], [0.8, 0.8, 0.8], 'FaceAlpha',.3);

end