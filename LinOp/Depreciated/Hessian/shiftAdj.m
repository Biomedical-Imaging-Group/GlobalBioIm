function As=shiftAdj(A,e,bc)

%SHIFTADJ adjoint of the Shifting operator with reflexive (mirror), 
%         periodic (circular) or zero boundary conditions
%   B = SHIFT(A,SHIFTSIZE) shifts the values in the array A by SHIFTSIZE
%   elements. SHIFTSIZE is a vector of integer scalars where the N-th
%   element specifies the shift amount for the N-th dimension of array A.
%   If an element in SHIFTSIZE is positive, the values of A are shifted
%   down (or to the right). If it is negative, the values of A are shifted
%   up (or to the left).
%   BC= 'reflexive' |'circular'|'zero'.

% X=randn(60,60,45,15);Y=randn(60,60,45,15);
% e=[randsrc(1,1,-60:60) randsrc(1,1,-60:60) randsrc(1,1,-45:45) randsrc(1,1,-15:15)]
% bc='reflexive';
% Xs=shift(X,e,bc);Ys_adj=shiftAdj(Y,e,bc);
% % <S(X),Y>=<X,S*(Y)>
% r=Xs(:).'*Y(:)-X(:).'*Ys_adj(:)

if nargin < 3
  bc='reflexive';
end

if nargin < 2
  error('shiftAdj:NoInputs',['No input arguments specified. ' ...
    'There should be at least two input arguments.'])
end

[e, sizeA, numDims, msg] = ParseInputs(A,e);

if (~isempty(msg))
  error('shiftAdj:InvalidShiftType','%s',msg);
end

if any(sizeA-abs(e) < 0)
  msg=[ 'The step of the shift at each dimension must be '...
    'less or equal than the size of the dimension itself.'];
  error('shiftAdj:InvalidShiftType','%s',msg);
end

switch bc
  case 'reflexive'
    idx=cell(numDims,2);
    idx2=idx;%cell(numDims,2);
    nidx=idx;%cell(numDims,2);
    
    for k=1:numDims
      if e(k) > 0
        idx(k,:)={1:sizeA(k)-e(k),1:e(k)};
        idx2(k,:)={1+e(k):sizeA(k),e(k):-1:1};
      else
        idx(k,:)={-e(k)+1:sizeA(k),sizeA(k)+e(k)+1:sizeA(k)};
        idx2(k,:)={1:sizeA(k)+e(k),sizeA(k):-1:sizeA(k)+e(k)+1};
      end
      nidx(k,:)={1:sizeA(k),1:sizeA(k)};
    end
    
    for k=1:numDims
      if k > 1
        A=As;
      end
      As=zeros(sizeA);
      midx=nidx;
      midx(k,:)=idx(k,:);
      midx2=nidx;
      midx2(k,:)=idx2(k,:);
      As(midx{:,1})=A(midx2{:,1});
      As(midx{:,2})=As(midx{:,2})+A(midx2{:,2});
    end
  case 'circular'
    idx=cell(numDims,1);
    for k = 1:numDims
      m      = sizeA(k);
      idx{k} = mod((0:m-1)+e(k), m)+1;
    end
    As=A(idx{:});
  case 'zero'
    As=zeros(sizeA);
    idx=cell(numDims,1);
    idx2=idx;
      for k=1:numDims
        if e(k) > 0	% shift right
          idx{k}=1:sizeA(k)-e(k);
          idx2{k}=1+e(k):sizeA(k);
        else		% shift left
          idx{k}=1-e(k):sizeA(k);
          idx2{k}=1:sizeA(k)+e(k);
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
