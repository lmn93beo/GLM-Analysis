function raster_cell = plotRaster_MGBA1_units(glmstruct, unit_name, start_time, end_time)
% unit_name: string indicating name of unit
% start_time and end_time: indices of start and end timestamps

% Make a cell that contains all units
raster_cell = {glmstruct.(unit_name)};
raster_cell = raster_cell(1:2:200);

plotRaster(raster_cell);

% Highlight stimulus period
patch([start_time end_time end_time start_time], [0 0 numel(raster_cell) numel(raster_cell)], [0.8, 0.8, 0.8], 'FaceAlpha',.3);

end