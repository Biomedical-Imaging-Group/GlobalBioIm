%% GUI-Header
% GUInotation-Tikh-
% GUIcallOpReg-LinOpGrad(InputSize,index,BC,res)-
% GUIparam-InputSize-vecInt-[]-Input size of the regularizer (e.g. [512 512]). 
% GUIparam-Index-vecInt-[]-Dimensions along which the gradient is computed (all by default)
% GUIparam-BC-dropDown/circular/zeros/mirror-circular-Boundary condition (default 'circular'), 'zeros', 'mirror'
% GUIparam-res-vecInt-[]-Vector containing the resolution along each dimension (default all 1)
% GUIcallCostReg-CostL2(OpReg.sizeout,Weight)-
% GUIparam-Weight-vecInt-[]-Weighting matrix (default 1)