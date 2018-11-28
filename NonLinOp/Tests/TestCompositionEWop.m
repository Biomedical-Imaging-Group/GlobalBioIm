
% Parameters 
sz=[2 2];

%% Build the CostGoodRoughness by compibing operators
% GR(x) sum ||grad x||^2./sqrt(abs(x).^2 bet)
% with bet=0 here
G=LinOpGrad(sz);                     % Gradient operator
sq=OpEWSquaredMagnitude(G.sizeout);  % Elem-wise square magnitude for numerator
S=LinOpSum(G.sizeout,3);             % To sum the two squared gradient components
iv=OpEWInverse(sz);                  % Elem-wise inverse 
sqr=OpEWSqrt(sz);                    % Elem-wise square root
sq2=OpEWSquaredMagnitude(sz);        % Elem-wise square magnitude for denominator (different size)
SS=LinOpSum(sz);                     % Final sum of all pixels

op=SS*((S*sq*G).*(iv*sqr*sq2));      % Combinaison of Maps
GR=CostGoodRoughness(G,0);           % Good Roughness instantiation

% Test if apply and gradient methods give identical results
% By the way, one can note that combining maps is less efficient than a one block implementation  
x=rand(sz);
tic;op*x,toc
tic;GR*x,toc
tic;g=GR.applyGrad(x),toc
tic;gg=op.applyJacobianT(1,x),toc

%% Build the smoothed TV regularizer by compibing operators
G=LinOpGrad(sz);                     % Gradient operator
sq=OpEWSquaredMagnitude(G.sizeout);  % Elem-wise square magnitude for numerator
S=LinOpSum(G.sizeout,3);             % To sum the two squared gradient components
sqr=OpEWSqrt(sz);                    % Elem-wise square root
SS=LinOpSum(sz);                     % Final sum of all pixels

op=SS*sqr*S*sq*G;                           % Combinaison of Maps
TV=CostHyperBolic([sz,2],0,3)*G;  % smoothed TV

% Test if apply and gradient methods give identical results
% By the way, one can note that combining maps is less efficient than a one block implementation  
x=rand(sz);
tic;op*x,toc
tic;TV*x,toc
tic;g=TV.applyGrad(x),toc
tic;gg=op.applyJacobianT(1,x),toc