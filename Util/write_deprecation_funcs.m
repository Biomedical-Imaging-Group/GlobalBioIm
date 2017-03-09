function write_deprecation_funcs()
% automatically generate functions to point from objects using the old
% naming convention to objects using the new naming convention, with
% warnings.

names = {
	'LinOp/Convolution' 'LinOp/LinOpConv'
 	'LinOp/Diag' 'LinOp/LinOpDiag'
	'LinOp/Grad' 'LinOp/LinOpGrad'
	'LinOp/Hess' 'LinOp/LinOpHess'
	'LinOp/Identity' 'LinOp/LinOpIdentity'
	};

for i = 1:length(names)
	oldName = names{i,1};
	newName = names{i,2};
	
	[p, nNew, e] = fileparts(newName);
	
	[p, n, e] = fileparts(oldName);
	if isempty(e), e = '.m'; end;
	autoName = fullfile(p, 'Auto', [n e]);
	
	fid = fopen(autoName, 'w');
	fprintf(fid, 'function obj = %s(varargin) \n', n);
	fprintf(fid, 'warning(''%s is deprecated, please use %s''); \n', n, nNew);
	fprintf(fid, 'obj = %s(varargin{:}); \n', nNew);
	fclose(fid);
end