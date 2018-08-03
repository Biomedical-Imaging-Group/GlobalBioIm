% out=rft3d(in) : performs a 3d rft for each 3d element in a dataset
function out=rft3d(in)
numdims=ndims(in);
if numdims > 4
    error('rift3d is only defined for up to 4 dimensions');
else
    in=expanddim(in,4);
end
if (size(in,4) > 1)
    out=newim(size(in,1),floor(size(in,2)/2)+1,size(in,3),size(in,4),'scomplex');
    for e=0:size(in,4)-1
        out(:,:,:,e)=rft(squeeze(in(:,:,:,e)));  % performs individual rift for each element
    end
else
        out=rft(squeeze(in));  % performs individual rift for each element
end

if numdims < 4
    mysize=size(out);
    out=reshape(out,mysize(1:numdims));
end
