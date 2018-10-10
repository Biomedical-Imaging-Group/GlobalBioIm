function [normA, v] = estimateNorm(A, maxNumApply, targetSNR, v0, verboseFlag)
% Estimates the operator norm of a map, A, using the power iteration
% method, which involves computing (A'*A)*...*(A'*A)*v.
%
% Simple usage: A.norm = estimateNorm(A)
% 
% Complicated usage: [Anorm, v] = estimateNorm(A, maxNumApply, targetSNR, v0, verboseFlag)
%
% Details:
% Anorm - lower bound on the operator norm of A.
%
% v - the unit vector that maximizes \| Av \|. Should be close to the first
% eigenvector of A.
%
% maxNumApply - maximum number of times to apply A. Higher is more
% accurate, but takes longer. Default: 200.
%
% targetSNR - stopping criterion. If snr(A*v, A*v - lambda*v) > targetSNR,
% we assume we have found the largest eigenvector and stop iterating.
% Higher is more accurate, but takes longer. Default 50.
%
% v0 - initial value for c. If you can approximate the largest eigenvector
% of A, using it for v0 will speed up the process, but leaving it out is
% fine. Default: uniform random vector on [-1, 1].
%
% WARNING: if you want expect complex inputs, you should provide a complex
% v0: v0 = 2*(rand(A.sizein)-.5) + 2i*(rand(A.sizein)-.5);
%
% verboseFlag - set to true for more output

% set defaults
if nargin < 2 || isempty(maxNumApply)
	maxNumApply = 200;
end
if nargin < 3 || isempty(targetSNR)
	targetSNR = 50;
end
if nargin < 4 || isempty(v0)
	s = rng;
	rng(42); % make initialization deterministic 
	v0 = 2*(rand(A.sizein)-.5);
	rng(s); % restore RNG state
end
if nargin < 5 || isempty(verboseFlag)
	verboseFlag = false;
end

Av = v0;

A = A' * A;


for iter = 1:maxNumApply
	% the old Av is the new v
	v = Av;	
	
	%% normalize v
	vNorm = sqrt(sum(abs(v(:)).^2));
	v = v / vNorm;

	% compute next power
	Av = A*v;
	
	% check accuracy 
	lambda = v(:)' * Av(:); % order matters if v is complex
	curSNR = snr(Av, Av - lambda*v); 
	
	if verboseFlag
		fprintf('%g, ', curSNR);
	end
	if curSNR >= targetSNR
		break
	end
	

end

if iter == maxNumApply
	warning('I did not find an eigenvector at the %g dB level after %d iterations', targetSNR, maxNumApply);
end


normA = sqrt(abs(lambda));


end
