%% Test Listener mechanism
%
% Before to run this script uncomment in Map.m the two lines
%   disp(['     In ',this.name,' property ',src.Name,' has been updated']);
%   disp(['     In ',this.name,' property ',sourc.Name,' has been modified']);
% which are in methods handlePostSetEvents and handleModifiedEvents,
% respectively.
%
clear; clc;
% Define some operators
disp('###### Define operators');
disp('sz = [256,256];h1=rand(sz);h1=h1/sum(h1(:));');
sz = [1024,1024];h1=rand(sz);h1=h1/sum(h1(:));
disp('D=LinOpDiag(sz,rand(sz));')
D=LinOpDiag(sz,rand(sz));
disp('H=LinOpConv(fftn(h1));')
H=LinOpConv(fftn(h1));
disp('L=LinOpDownsample(H.sizein,[1 1]);')
L=LinOpDownsample(H.sizein,[1 1]);
x=rand(sz);

disp('###### Case 1');
disp('T=H+D;')
T=H+D;
disp('CostL2(T.sizeout,rand(sz))*T;')
F=CostL2(T.sizeout,rand(sz))*T;
disp('F.H2.mapsCell{1}=(L*H);')
F.H2.mapsCell{1}=(L*H);
disp(['F.H2.mapsCell{1}.norm = ',num2str(F.H2.mapsCell{1}.norm)])
disp('F.doPrecomputation=1;');
F.doPrecomputation=1;
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),')']);
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- effect of precomputation)']);
disp(['0.5*norm(L*(H*x)+D*x-F.H1.y,''fro'')^2 = ',num2str(0.5*norm(F.H2.apply(x)-F.H1.y,'fro')^2)]);
disp('F.H2.mapsCell{1}.H2.mtf=fftn(3*h1);')
F.H2.mapsCell{1}.H2.mtf=fftn(3*h1);
disp(['F.H2.mapsCell{1}.norm = ',num2str(F.H2.mapsCell{1}.norm)])
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- redo precomputation)']);
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- effect of precomputation)']);
disp(['0.5*norm(L*(H*x)+D*x-F.H1.y,''fro'')^2 = ',num2str(0.5*norm(F.H2.apply(x)-F.H1.y,'fro')^2)]);
disp('F.memoizeOpts.apply=1;');
F.memoizeOpts.apply=1;
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- effect of precomputation)']);
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- effect of memoize)']);
disp('F.H2.mapsCell{2}.diag=rand(sz);')
F.H2.mapsCell{2}.diag=rand(sz);
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- redo precomputation + memoize)']);
t=tic;disp(['F*x = ',num2str(F*x),' (elapsed time ',num2str(toc(t)),' <- effect of memoize)']);