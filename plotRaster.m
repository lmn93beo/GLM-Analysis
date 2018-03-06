
function plotRaster(SpikeTimesCell, clr);

if nargin < 2; clr = 'k'; end;

set(gcf,'color','w');

mm = max(size(SpikeTimesCell));

if mm > 1

for t = 1 : max(size(SpikeTimesCell));
    times = SpikeTimesCell{t};
    
    for i = 1:length(times);
        line( [times(i) times(i)], [t, t+1], 'color',clr,'linewidth',1.3);
        hold on;
    end;
    
end;

ylim( [1,t] ); xlim( [-0.5, 1.5] );
set(gca,'tickdir','out', 'yaxislocation','left'); box off;

end;