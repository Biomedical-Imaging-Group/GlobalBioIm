function [grid] = makeGrid(varargin)

% utility function for explicitly listing the points in the grid specified
% by (start,step,size). Useful with the KaiserBesselWindow object's
% methods and elsewhere.

% [grid] = makeGrid(xStart, xStep, xSize)
%  [grid] = makeGrid(xGridParams)

if nargin == 1
	xParams = varargin{1};
	xStart = xParams( 1, : );
	xStep  = xParams( 2, : );
	xSize  = xParams( 3, : );
elseif nargin == 3
		xStart = varargin{1};
		xStep = varargin{2};
		xSize = varargin{3};
		
else
	error('bad input')
end


D = length(xStart);


vects = make_grid_vectors(xStart, xStep, xSize);

switch 2
	case 1 % no-mex way, takes >2X the memory of the eventual output
		grid = cell(D,1);
		[grid{:}] = ndgrid( vects{:} );
		
		grid = cellfun(@(a)(a(:)), grid, 'uniformoutput', false)';
		grid = [grid{:}]';
	case 2
		grid = ndgridMatrix( vects{:} );
end
