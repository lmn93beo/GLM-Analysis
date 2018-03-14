function make_psth(glmstruct, unit_name, binSize)

figure;
% Make a cell that contains all units
raster_cell = {glmstruct.(unit_name)};

raster_times = cell2mat(raster_cell');

% Define edges for histogram
edges = 1:binSize:3000;
counts = histcounts(raster_times, edges);
plot(counts)

end