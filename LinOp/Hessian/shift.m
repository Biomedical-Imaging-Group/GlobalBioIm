function As=shift(A,e,bc)

%SHIFT Shifting operator with reflexive (mirror), periodic (circular) or
%      zero boundary conditions
%   B = SHIFT(A,SHIFTSIZE,BC) shifts the values in the array A by SHIFTSIZE
%   elements. SHIFTSIZE is a vector of integer scalars where the N-th
%   element specifies the shift amount for the N-th dimension of array A.
%   If an element in SHIFTSIZE is positive, the values of A are shifted
%   down (or to the right). If it is negative, the values of A are shifted
%   up (or to the left).
%   BC= 'reflexive' |'circular'|'zero'.

if nargin < 3
  bc='reflexive';
end

if nargin < 2
  error('shift:NoInputs',['No input arguments specified. ' ...
    'There should be at least two input arguments.'])
end

[e, sizeA, numDims, msg] = ParseInputs(A,e);

if (~isempty(msg))
  error('shift:InvalidShiftType','%s',msg);
end

if any(sizeA-abs(e) < 0)
  msg=[ 'The step of the shift at each dimension must be '...
    'less or equal than the size of the dimension itself.'];
  error('shift:InvalidShiftType','%s',msg);
end

idx=cell(numDims,1);

switch bc
  case 'reflexive'
    for k=1:numDims
      if e(k) > 0
        idx{k}=[e(k):-1:1,1:sizeA(k)-e(k)];
      else
        idx{k}=[1-e(k):sizeA(k),sizeA(k):-1:sizeA(k)+e(k)+1];
      end
    end
    As=A(idx{:});
  case 'circular'
    % Loop through each dimension of the input matrix to calculate shifted indices
    for k=1:numDims
      m      = sizeA(k);
      idx{k} = mod((0:m-1)-e(k), m)+1;
    end
    As=A(idx{:});
    case 'zero'
      As=zeros(sizeA);
      idx2=idx;
      for k=1:numDims
        if e(k) > 0	% shift right
          idx{k}=e(k)+1:sizeA(k);
          idx2{k}=1:sizeA(k)-e(k);
        else		% shift left
          idx{k}=1:sizeA(k)+e(k);
          idx2{k}=1-e(k):sizeA(k);
        end
      end
      As(idx{:})=A(idx2{:});
  otherwise
    error('shift:InvalidShiftType','%s','Unknown boundary conditions.');
end


function [e, sizeA, numDimsA, msg] = ParseInputs(A,e)

% default values
sizeA    = size(A);
numDimsA = ndims(A);
msg      = '';

% Make sure that SHIFTSIZE input is a finite, real integer vector
sh        = e(:);
isFinite  = all(isfinite(sh));
nonSparse = all(~issparse(sh));
isInteger = all(isa(sh,'double') & (imag(sh)==0) & (sh==round(sh)));
isVector  = ((ndims(e) == 2) && ((size(e,1) == 1) || (size(e,2) == 1)));

if ~(isFinite && isInteger && isVector && nonSparse)
  msg = ['Invalid shift type: ' ...
    'must be a finite, nonsparse, real integer vector.'];
  return;
end

% Make sure the shift vector has the same length as numDimsA.
% The missing shift values are assumed to be 0. The extra
% shift values are ignored when the shift vector is longer
% than numDimsA.
if (numel(e) < numDimsA)
  e(numDimsA) = 0;
elseif (numel(e) > numDimsA)
  e(numDimsA+1:numel(e))=[];
end

