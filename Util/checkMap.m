function s = checkMap(H)
% function to check the consistency of the Map H, including specialized
% checks for LinOps
% returns a stats object with information about which tests passed

fprintf('-- Checking Map with name %s--\n', H.name);

% check in
if ~isnumeric(H.sizein) || ~isnumeric(H.sizeout)
	fprintf('H dimensions are not set, cannot do automatic testing\n');
	return;
end

% create random input/output-sized vectors
x = randn(H.sizein);
if H.isComplexIn
	x = x + 1i * randn(H.sizein);
end

y = randn(H.sizeout);
if H.isComplexOut
	y = y + 1i * randn(H.sizein);
end


% apply
try
	Hx = H.apply(x);
	s.applyOK = true;
	fprintf('apply OK\n');
catch ME
	fprintf('apply FAILs:\n\t%s\n', ME.message');
	s.applyOK = false;
end

% applyJacobianT
if H.isDifferentiable
	try
		H.applyJacobianT(y, x);
		s.applyJacobianTOK = true;
		fprintf('applyJacobianT OK\n');
	catch ME
		fprintf('H.isDifferentiable, but applyJacobianT FAILs:\n\t%s\n', ME.message);
		s.applyJacobianTOK = false;
	end
end

% applyInverse
if H.isInvertible
	try
		xhat = H.applyInverse(Hx);
		s.inverseOK = true;
		fprintf('applyInverse OK\n');
	catch ME
		s.inverseOK = false;
		fprintf('H.isInvertible, but applyInverse FAILs:\n\t%s\n', ME.message);
	end
end

if s.applyOK && s.inverseOK
	curSNR = snr(x, x-xhat);
	if curSNR > 70
		okString = 'OK';
	else
		okString = 'FAIL';
	end
	fprintf('\taccurate to %d dB, %s\n', curSNR, okString);
else
	fprintf('\tcannot assess accuracy\n');
end


if isa(H, 'LinOp')
	fprintf('-- LinOp-specific checks --\n')
	
	% adjoint
	try
		HTy = H.applyAdjoint(y);
		s.adjointOK = true;
		fprintf('applyAdjoint OK\n');
	catch ME
		fprintf('applyAdjoint fails:\n\t%s\n', ME.message');
		s.adjointOK = false;
	end
	
	if s.applyOK && s.adjointOK
		lhs = x(:).' * HTy(:);
		rhs = y(:).' * Hx(:);
		curSNR = snr(lhs, lhs-rhs);
		if curSNR > 70
			okString = 'OK';
		else
			okString = 'FAIL';
		end
		fprintf('\taccurate to %d dB, %s\n', curSNR, okString);
	else
		fprintf('\tcannot assess accuracy\n');
	end
	
	% HtH
	try
		HTHx = H.applyHtH(x);
		s.applyHtHOK = true;
		fprintf('applyHtH OK\n');
	catch ME
		fprintf('applyHtH fails:\n\t%s\n', ME.message');
		s.applyHtHOK = false;
	end
	
	if s.applyOK && s.adjointOK && s.applyHtHOK
		lhs = HTHx;
		rhs = H.applyAdjoint( Hx );
		curSNR = snr(lhs, lhs-rhs);
		if curSNR > 70
			okString = 'OK';
		else
			okString = 'FAIL';
		end
		fprintf('\taccurate to %d dB, %s\n', curSNR, okString);
	else
		fprintf('\tcannot assess accuracy\n');
	end
	
	% HHt
	try
		HHty = H.applyHHt(y);
		s.applyHHtOK = true;
		fprintf('applyHHt OK\n');
	catch ME
		fprintf('applyHHt fails:\n\t%s\n', ME.message');
		s.applyHHtOK = false;
	end
	
	if s.applyOK && s.adjointOK && s.applyHHtOK
		lhs = HHty;
		rhs = H.apply( HTy );
		curSNR = snr(lhs, lhs-rhs);
		if curSNR > 70
			okString = 'OK';
		else
			okString = 'FAIL';
		end
		fprintf('\taccurate to %d dB, %s\n', curSNR, okString);
	else
		fprintf('\tcannot assess accuracy\n');
	end

end
