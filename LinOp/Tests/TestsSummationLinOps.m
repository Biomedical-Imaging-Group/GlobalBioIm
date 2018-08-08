%% Test LinOp Summations

% Define some operators
sz = [128,128];
D=LinOpDiag(sz,rand(sz));
H=LinOpConv(rand(sz));
x=rand(sz);
fprintf('\nOperators used : \n \t - D : LinOpDiag  \n \t - H : LinOpConv \n \n');
tol=1e-13;
isOk=true;
me=[];

%% Tests
fprintf('Test #1 : \n ');
T=D+H+D;
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t D+H+D            --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (D+H+D)x = Dx+Hx+Dx')
if norm(T*x-(D*x+H*x+D*x)) < tol
    fprintf('  --> OK \n');
else   
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (D+H+D)x ~= Dx+Hx+Dx ');
end
T=D+D+H;
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t D+D+H            --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (D+D+H)x = Dx+Hx+Dx')
if norm(T*x-(D*x+H*x+D*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (D+D+H)x ~= Dx+Hx+Dx ');
end
T=H+D+D;
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t H+D+D            --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (H+D+D)x = Dx+Hx+Dx')
if norm(T*x-(D*x+H*x+D*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (H+D+D)x ~= Dx+Hx+Dx ');
end
T=D+(D+H);
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t D+(D+H)          --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (D+(D+H))x = Dx+Hx+Dx')
if norm(T*x-(D*x+H*x+D*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (D+(D+H))x ~= Dx+Hx+Dx ');
end
T=D+H+H+D;
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t D+H+H+D          --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (D+H+H+D)x = Dx+Hx+Hx+Dx')
if norm(T*x-(D*x+H*x+D*x+H*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (D+H+H+D)x ~= Dx+Hx+Hx+Dx ');
end
T=(D+H)+(H+D);
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t (D+H)+(H+D)      --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if ((D+H)+(H+D))x = Dx+Hx+Hx+Dx')
if norm(T*x-(D*x+H*x+D*x+H*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' ((D+H)+(H+D))x ~= Dx+Hx+Hx+Dx ');
end
T=D-(H+0.5*D);
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t D-(H+0.5*D)      --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (D-(H+0.5*D))x = Dx-Hx-0.5*Dx')
if norm(T*x-(D*x-H*x-0.5*D*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (D-(H+0.5*D))x ~= Dx-Hx-0.5*Dx ');
end
T=D+H-0.5*D;
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t D+H-0.5*D        --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if (D+H-0.5*D)x = Dx+Hx-0.5*Dx')
if norm(T*x-(D*x+H*x-0.5*D*x)) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (D+H-0.5*D)x ~= Dx+Hx-0.5*Dx ');
end
T=(D+H)-0.5*(H+D);
mapCell=cellfun(@(x) x.name,T.mapsCell,'UniformOutput',false);
for ii=1:length(mapCell)-1,mapCell{ii}=[mapCell{ii},' + '];end
fprintf(['\t (D+H)-0.5*(H+D)  --> ',T.name,': ',mapCell{:},'\n']);
fprintf( '\t                      Check if ((D+H)-0.5*(H+D))x = Dx+Hx-0.5*(Hx+Dx)')
if norm(T*x-(D*x+H*x-0.5*(D*x+H*x))) < tol
    fprintf('  --> OK \n');
else
    fprintf('  --> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' ((D+H)-0.5*(H+D))x ~= Dx+Hx-0.5*(Hx+Dx) ');
end

global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    stateTest=isOk;
    message=me;
end