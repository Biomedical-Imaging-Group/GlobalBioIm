%% GUI-Header
% GUInotation-GR-
% GUIcallOpReg-LinOpIdentity(InputSize)-
% GUIparam-index-vecInt-[]-Dimensions along which the gradient is computed (all by default)
% GUIparam-BC-dropDown/circular/zeros/mirror-circular-Boundary condition (default 'circular'), 'zeros', 'mirror'
% GUIparam-res-vecInt-[]-Vector containing the resolution along each dimension (default all 1)
% GUIcallCostReg-CostGoodRoughness(LinOpGrad(InputSize,index,BC,res),bet)-
% GUIparam-bet-double-0.1-Smoothing parameter (default 1e-1)