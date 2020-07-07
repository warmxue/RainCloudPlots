% Use like: h = rm_raincloud(data, colours, plot_top_to_bottom, density_type, bandwidth)
% Where 'data' is an M x N cell array of N data series and M measurements
% And 'colours' is an N x 3 array defining the colour to plot each data series
% plot_top_to_bottom: Default plots left-to-right, set to 1 to rotate.
% density_type: 'ks' (default) or 'RASH'. 'ks' uses matlab's inbuilt 'ksdensity' to
% determine the shape of the rainclouds. 'RASH' will use the 'rst_RASH'
% method from Cyril Pernet's Robust Stats toolbox, if that function is on
% your matlab path.
% bandwidth: If density_type == 'ks', determines bandwidth of density estimate
% h is a cell array containing handles of the various figure parts:
% h.p{i,j} is the handle to the density plot from data{i,j}
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}
% h.m(i,j) is the handle to the single, large dot that represents mean(data{i,j})
% h.l(i,j) is the handle for the line connecting h.m(i,j) and h.m(i+1,j)

%% TO-DO:
% Patch can create colour gradients using the 'interp' option to 'FaceColor'. Allow this?

function h = rm_raincloud(data, colours, plot_top_to_bottom, density_type, bandwidth)
%% check dimensions of data

[n_plots_per_series, n_series] = size(data);

% make sure we have correct number of colours
assert(all(size(colours) == [n_series 3]), 'number of colors does not match number of plot series');

%% default arguments

if nargin < 3
    plot_top_to_bottom  = 0;    % left-to-right plotting by default
end

if nargin < 4
    density_type        = 'ks'; % use 'ksdensity' to create cloud shapes
end

if nargin < 5
    bandwidth           = [];   % let the function specify the bandwidth
end

%% Calculate properties of density plots

% Probably okay to hard-code this as it just determines the granularity of
% the density estimate
density_granularity = 200;

n_bins = repmat(density_granularity, n_plots_per_series, n_series);

% calculate kernel densities
for i = 1:n_plots_per_series
    for j = 1:n_series
       
        switch density_type
            
            case 'ks'
                
                % compute density using 'ksdensity'
                [ks{i, j}, x{i, j}] = ksdensity(data{i, j}, 'NumPoints', n_bins(i, j), 'bandwidth', bandwidth);
                
            case 'rash'
                
                % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found 
                assert(exist('rst_RASH', 'file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');
                
                % compute density using RASH
                [x{i, j}, ks{i, j}] = rst_RASH(data{i, j});
                
                % override default 'n_bins' as rst_RASH determines number of bins
                n_bins(i, j) = size(ks{i, j}, 2);
        end
        
        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
        q{i, j}     = (1:n_bins(i, j) - 1)';
        faces{i, j} = [q{i, j}, q{i, j} + 1, q{i, j} + n_bins(i, j) + 1, q{i, j} + n_bins(i, j)];
        
    end
end

% determine spacing between plots, in y direction
spacing     = 4 * mean(mean(cellfun(@max, ks)));% default: 2
ks_offsets  = [0:n_plots_per_series-1] .* spacing;

% flip so first plot in series is plotted on the *top*
ks_offsets  = fliplr(ks_offsets);

% calculate patch vertices from kernel density
for i = 1:n_plots_per_series
    for j = 1:n_series
        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
    end
end


%% jitter for the raindrops

jit_width = spacing / (2*(n_series + 1));% default: 8

%% linewidth
line_width = 2;
%% transparency
alpha = 0.5;
%% black color
bcolor = [0 0 0];

% TO-DO: This should probably not be hardcoded either...
% Although it can be overridden later
raindrop_size = 50;%default: 100;

for i = 1:n_plots_per_series % timepoints
    for j = 1:n_series % group
        jit{i,j} = j * jit_width + rand(1, length(data{i,j})) * jit_width;
    end
    
end

%% means (for mean dots)

cell_means = cellfun(@mean, data);
cell_quartiles = cellfun(@(x) quantile(x, [0.25 0.75 0.5]), data, 'UniformOutput', false);
cell_iqr = cellfun(@(x) x(2) - x(1), cell_quartiles, 'UniformOutput', false);
cell_Xs = cellfun(@(x) sort(x), data, 'UniformOutput', false);
cell_whiskers = cellfun(@(x, y, z) [min(x(x > (y(1) - (1.5 * z)))), max(x(x < (y(2) + (1.5 * z))))], cell_Xs, cell_quartiles, cell_iqr, 'UniformOutput', false);
% 0.25, 0.75, median, min, max
cell_Y = cellfun(@(x, y) [x, y], cell_quartiles, cell_whiskers, 'UniformOutput', false);
%% plot
% note - we *could* plot everything here in one big loop, but then
% different figure parts would overlay each other in a silly way.
% figure,
hold on

% patches
for i = 1:n_plots_per_series
    for j = 1:n_series
        
        % plot patches
        h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', colours(j, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', alpha);
        
        % scatter rainclouds
        h.s{i, j} = scatter(data{i, j}, -jit{i, j} + ks_offsets(i), 'MarkerFaceColor', colours(j, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha, 'SizeData', raindrop_size);
        
    end
end

% plot mean lines
for i = 1:n_plots_per_series - 1 % We have n_plots_per_series-1 lines because lines connect pairs of points
    for j = 1:n_series
        
        h.l(i, j) = line(cell_means([i i+1], j), ks_offsets([i i+1]), 'LineStyle', ':', 'LineWidth', line_width, 'Color', colours(j, :));
        
    end
end

% plot mean dots
for i = 1:n_plots_per_series
    for j = 1:n_series
        
        h.m(i, j) = scatter(cell_means(i, j), ks_offsets(i), 'MarkerFaceColor', colours(j, :), 'MarkerEdgeColor', bcolor, 'MarkerFaceAlpha', 1, 'SizeData', raindrop_size * 2, 'LineWidth', line_width);
    
    end
end

% width of boxplot
% make some space under the density plot for the boxplot and raindrops
yl = get(gca, 'YLim');
% set(gca, 'YLim', [-yl(2)*lwr_bnd yl(2)]);
wdth = 0.1;%yl(2) * 0.25;
% plot whiskers
for i = 1:n_plots_per_series
    for j = 1:n_series
        Y = cell_Y{i, j};
        box_pos = [Y(1), mean(-jit{i, j}) + ks_offsets(i)-0.5*wdth, Y(2)-Y(1), wdth]; %[x, y, w, h]
        % mean line
        h.ml(i, j) = line([Y(3) Y(3)], [(mean(-jit{i, j}) + ks_offsets(i)-(wdth * 0.5)), (mean(-jit{i, j}) +ks_offsets(i) + (wdth * 0.5))], 'col', bcolor, 'LineWidth', line_width);
        
        % whiskers
        h.w1(i, j) = line([Y(2) Y(5)], [(mean(-jit{i, j}) +ks_offsets(i)), (mean(-jit{i, j}) +ks_offsets(i))], 'col', [bcolor, alpha], 'LineWidth', line_width);
        h.w2(i, j) = line([Y(1) Y(4)], [(mean(-jit{i, j}) +ks_offsets(i)), (mean(-jit{i, j}) +ks_offsets(i))], 'col', [bcolor, alpha], 'LineWidth', line_width);
        
        % 'box' of 'boxplot'
        h.b(i,j) = rectangle('Position', box_pos);
        set(h.b(i,j), 'EdgeColor', [bcolor, alpha])
%         set(h.b(i, j), 'FaceColor', [colours(j, :), alpha])
        set(h.b(i,j), 'LineWidth', line_width);
        %set(h{3}, 'FaceColor', bxfacecl);
        % could also set 'FaceColor' here, etc
    end
end


%% clear up axis labels

% 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
% raincloud at the top. So flip the vector around
set(gca, 'YTick', fliplr(ks_offsets));

set(gca, 'YTickLabel', n_plots_per_series:-1:1);

%% determine plot rotation
% default option is left-to-right
% plot_top_to_bottom can be set to 1 
% NOTE: Because it's easier, we actually prepare everything plotted
% top-to-bottom, then - by default - we rotate it here. That's why the
% logical is constructed the way it is.

% rotate and flip
if ~plot_top_to_bottom
    view([90 -90]);
    axis ij
end
