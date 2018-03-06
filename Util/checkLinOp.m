function s = checkLinOp(H,checkComplex)
% see checkMap
if ~exist('checkComplex', 'var') || isempty(checkComplex)
    checkComplex = false;
end
s = checkMap(H,checkComplex);