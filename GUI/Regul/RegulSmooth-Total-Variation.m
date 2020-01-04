%% GUI-Header
% GUInotation-sTV-
% GUIcallOpReg-LinOpGrad(InputSize,index,BC,res)-
% GUIparam-index-vecInt-[]-Dimensions along which the gradient is computed (all by default)
% GUIparam-BC-dropDown/circular/zeros/mirror-circular-Boundary condition (default 'circular'), 'zeros', 'mirror'
% GUIparam-res-vecInt-[]-Vector containing the resolution along each dimension (default all 1)
% GUIparam-epsilon-double-0.001-Smoothing parameter (default 1e-3)
% GUIcallCostReg-CostHyperBolic(OpReg.sizeout,epsilon,length(OpReg.sizeout))