function [im_h, axis_h, fig_h] = prettyPlot(f, labels, xStart, xStep, varargin)
% [im_h, axis_h, fig_h] = prettyPlot(f, labels, xStart, xStep, varargin)
%
% attractive plotting for 2d and 3d images
% 
% f - data to plot, can be Nx1, NxM, or NxMxP
% labels - (optional) cell array of the names for the axes
% xStart - (optional) numDims x 1 vector of location of f(1,1,1)
% xStep - (optional) numDims x 1 vector of the step length in each
% dimension
% varargin - any additional ... 'key', 'value',... pairs are passed to
% matlabs plot() or imagesc() functions
%
% [im_h, axis_h, fig_h] - handles to the plot/image(s), axis(es), and
% figure(s) returned
%
% Note: assumes that the rows of f correspond to x_0 and the columns to x_1,
% meaning that you do _not_ transpose f before passing it
%
% example 1, simple usage:
% [x,y,z] = ndgrid(linspace(-1,1,100));
% prettyPlot(x.^2 + y.^2 + z.^2 < 1);
%
% example 2, setting units and labels
% [x,y,z] = ndgrid(-2:.1:2);
% [im_h, axis_h, fig_h] = prettyPlot( (z+1) .*double(x.^2 + y.^2<1), {'x', 'y', 'z'}, [-2 -2 -2], [.1 .1 .1]);
%
% Michael McCann michael.mccann@epfl.ch

colors = lines(3);

if ~exist('xStart', 'var') || isempty(xStart)
	xStart = ones(1, length(size(f)));
	xStep = ones(1, length(size(f)));
end

xSize = size(f);

if isstruct(xStart) % hack to allow passing params in as xStart
	xStep = xStart.xStep;
	xSize = xStart.xSize;
	xStart = xStart.xStart;
end

if ~exist('labels', 'var') || isempty(labels)
	labels = {'x_0', 'x_1', 'x_2'};
end

if isvector(f) 
	%% 1d plotting code
	if ~isempty(get(gcf, 'CurrentAxes'))
		fig_h = figure;
	end
	im_h = plot( (0:xSize(1)-1)*xStep(1) + xStart(1), f, varargin{:});
	xlabel(labels{1}, 'rot', 0)
	
elseif ismatrix(f) 
	%% 2d plotting code
	if ~isempty(get(gcf, 'CurrentAxes'))
		fig_h = figure;
	end
	im_h = imagesc( (0:xSize(1)-1)*xStep(1) + xStart(1), ...
		(0:xSize(2)-1)*xStep(2) + xStart(2), f', varargin{:});
	
	colorbar
	
	axis xy
	%truesize
	axis equal
	
	xlabel(labels{1}, 'rot', 0)
	ylabel(labels{2}, 'rot', 0)
	
elseif ndims(f) == 3
	%% 3d plotting code
	xEnd = xStart + xStep.*(xSize-1);
	
	
	sliceInd = ceil(size(f)/2);
	

	centerPoint = xStart + (sliceInd-1).*xStep;
	l1_h = zeros(1,3);
	l2_h = zeros(1,3);
	
	im_h = zeros(1, 3);
	fig_h = zeros(1,3);
	axis_h = zeros(1,3);
	
	fig_h(1) = figure;
	im_h(1) = prettyPlot(squeeze(f(sliceInd(1), :, :)), labels(2:3), xStart(2:3), xStep(2:3), varargin{:});
	axis_h(1) = gca;
	update_lines(im_h(1))
	
	fig_h(2) = figure;
	im_h(2) = prettyPlot(squeeze(f(:,sliceInd(2), :)), labels([1 3]), xStart([1 3]), xStep([1 3]), varargin{:} );
	axis_h(2) = gca;
	update_lines(im_h(2))
	
	fig_h(3) = figure;
	im_h(3) = prettyPlot(squeeze(f(:,:,sliceInd(3))), labels(1:2), xStart(1:2), xStep(1:2), varargin{:});
	axis_h(3) = gca;
	update_lines(im_h(3))
	
	% link color scales
	clims = get(axis_h, 'clim');
	clims = cell2mat(clims);

	if ~strcmp(get(gca, 'CLimMode'), 'manual')
		set(axis_h(1), 'clim', [min(f(:)), max(f(:)) + eps]); % + eps avoids limits being equal
	end
	linkprop(axis_h, 'clim');
	
	linkprop(fig_h, 'colormap');
	
	
	% setup callbacks
	set(im_h, 'ButtonDownFcn', @start_updating);
	set(axis_h, 'ButtonDownFcn', @start_updating);
	set(fig_h, 'WindowButtonUpFcn', @stop_updating);
	
	% update the plots when zoomed or panned
	for i = 1:3
		set(zoom(axis_h(i)), 'ActionPostCallback',@update_lines);
		set(pan(fig_h(i)), 'ActionPostCallback',@update_lines);
		set(fig_h(i),'ResizeFcn',@update_lines);
	end
	


	
	
	
else
	error('can only plot 2d and 3d data');
end



	function start_updating(src, ~)
		dim = (src==im_h | src==fig_h);
		set(fig_h(dim), 'WindowButtonMotionFcn', @update_center);
		update_center(src)
	end

	function stop_updating(~, ~)
		set(fig_h(ishandle(fig_h)), 'WindowButtonMotionFcn', '');

	end

	function update_center(src, ~)
		dim = (src==im_h | src==fig_h | src==axis_h);
		
		point = get(axis_h(dim), 'currentpoint');
		
		centerPoint(~dim) = point(1,1:2);
		centerPoint = round((centerPoint-xStart) ./ xStep) .* xStep + xStart;
		
		% fix out of bounds
		small = centerPoint < xStart;
		centerPoint(small) = xStart(small);
		big = centerPoint > xEnd;
		centerPoint(big) = xEnd(big);
		
		
		for i = 1:3
			if ishandle(im_h(i)) % deal with closed figures
				update_lines(im_h(i));
			end
		end
		
		
	end

	function update_lines(src, ~)
		dim = (src==im_h | src==fig_h | src==axis_h);
		
		xlim = get(axis_h(dim), 'xlim');
		ylim = get(axis_h(dim), 'ylim');
		
		if ~l1_h(dim)
			c = colors(~dim, :);
			l1_h(dim) = line('hittest', 'off', 'clipping', 'off', 'color', c(1,:)); % don't let them intercept clicks
			l2_h(dim) = line('hittest', 'off', 'clipping', 'off', 'color', c(2,:));
		end
		
		p = centerPoint(~dim);
		set(l1_h(dim), 'xdata', xlim);
		set(l1_h(dim), 'ydata', p(2) * [1 1]);
		set(l2_h(dim), 'xdata', p(1) * [1 1]);
		set(l2_h(dim), 'ydata', ylim);
		
		inds = {1:xSize(1), 1:xSize(2), 1:xSize(3)};
		inds{dim} =round( (centerPoint(dim)-xStart(dim))/xStep(dim) + 1 );
		set(im_h(dim), 'CData', squeeze(f(inds{:}))');
		
	end

end

