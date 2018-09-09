function H = LinOpDFT(sz,unitary,pad)
% see LinOpSDFT
if nargin <2, unitary=false; end;
if nargin <3, pad=[]; end;
H = LinOpSDFT(sz,[],unitary,pad);
end
