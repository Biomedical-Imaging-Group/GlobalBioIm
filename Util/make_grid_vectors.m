function vects = make_grid_vectors(xStart, xStep, xSize)

D = length(xStart);


vects = cell(D, 1);
for d = 1:D
	vects{d} = (0:xSize(d)-1)*xStep(d) + xStart(d);
end