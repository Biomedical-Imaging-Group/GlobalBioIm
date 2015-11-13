function y = applyKernel(kHat,x)
% this function applies the kernel, kHat, on an image, x
% the output has the same size as x, but the multiplication in freq occurs
% padded to the size of kHat
%
% BEWARE: the center of k should have been k(1)


y  =  ifftn( fftn(x,size(kHat)) .* kHat, 'symmetric' );
y  =  y(1:size(x,1), 1:size(x,2), 1:size(x,3)); 