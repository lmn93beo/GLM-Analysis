trial_start = 3000;
trial_end = 4000;

for i = 1:numel(glmtrial)
    fprintf('Plotting trial %d\n', i);
    clf; c = plotRaster_MGBA1(glmtrial, i , trial_start, trial_end);
    savefig(['trial' num2str(i) '_all_responsive_units_raster.fig']);
end