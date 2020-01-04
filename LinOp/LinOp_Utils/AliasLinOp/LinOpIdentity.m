function H = LinOpIdentity(sz)

    %% GUI-Header
    % GUInotation-I-
    % GUIcall-LinOpIdentity(InputSize)-
    % GUIparam-InputSize-vecInt-[]-Input size of the gradient operator (e.g. [512 512]).
    
% see LinOpDiag
H = LinOpDiag(sz);
end