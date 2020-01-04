function [listOpti,nameInLib] = GetOptiList(Costs,Ops,NameCosts,NamesOps,isPos,doPrecom)
global sizein
warning off
% Initialize listOpti
nOpti=1;
listOpti={};

sizein=Ops{1}{end}.sizein;

% Precomputations
CompOps={};
CompCosts={};
for ii=1:length(Ops)
    H=Ops{ii}{1};
    for jj=2:length(Ops{ii})
        H = H * Ops{ii}{jj};
    end
    CompOps{ii} = H;
    CompCost{ii} = Costs{ii}*H;
    isdiff(ii) = CompCost{ii}.isDifferentiable;
    isproxComp(ii) = TestProx(CompCost{ii});
    isprox(ii) = TestProx(Costs{ii});
    isL2(ii) = strcmp(Costs{ii}.name,'CostL2') || (strcmp(Costs{ii}.name,'CostMultiplication') && strcmp(Costs{ii}.cost2.name,'CostL2'));
end
nbCosts=length(Costs);


%% Gradient Descent
if ~isPos && all(isdiff)
    % -- General
    listOpti{nOpti}.name = 'Gradient Descent';
    nameInLib{nOpti} = 'GradDsct';
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
    if doPrecom,  listOpti{nOpti}.call{1} = [ listOpti{nOpti}.call{1},' cf.doPrecomputation =1;']; end
    listOpti{nOpti}.call{2} = 'Opt = OptiGradDsct(cf);';
    % - Gam parameter
    listOpti{nOpti}.parameters{1}.name = 'gam';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1e-1';
    listOpti{nOpti}.parameters{1}.default = '0';
    listOpti{nOpti}.parameters{1}.info = 'Descent step (Put 0 to let gam be computed automatically, see Doc to known when this is applicable)';
    listOpti{nOpti}.parameters{1}.toSet =1;
    nOpti=nOpti+1;
end

%% Forward-Backward Splitting
if all(isdiff) && isPos
    % If all terms are differentiables
    listOpti{nOpti}.call{1} = ['F = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
    if doPrecom,  listOpti{nOpti}.call{1} = [ listOpti{nOpti}.call{1},' F.doPrecomputation =1;']; end
    listOpti{nOpti}.call{2} = 'P = CostNonNeg(F.sizein);';
    listOpti{nOpti}.call{3} = 'Opt = OptiFBS(F,P);';
    isFBS=1;
elseif sum(isdiff)>0
    % Otherwise: chck that the sum of the non-diff terms is prox
    F=BuildSumCosts(CompCost,~isdiff,isPos);   % Sum all non-differentiable costs 
    isFBS=TestProx(F);                         % Test if prox
    if isFBS
        listOpti{nOpti}.call{1} = ['F = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
        if doPrecom,  listOpti{nOpti}.call{1} = [ listOpti{nOpti}.call{1},' F.doPrecomputation =1;']; end
        if isPos
            listOpti{nOpti}.call{2} = 'P = CostNonNeg(F.sizein);';
            listOpti{nOpti}.call{3} = ['G = ',GetExprCompCost(NameCosts,NamesOps,find(~isdiff)),' + P;'];
            listOpti{nOpti}.call{4} = 'Opt = OptiFBS(F,G);';
        else
            listOpti{nOpti}.call{2} = ['G = ',GetExprCompCost(NameCosts,NamesOps,find(~isdiff)),';'];
            listOpti{nOpti}.call{3} = 'Opt = OptiFBS(F,G);';
        end
    end
else
    isFBS=0;
end
% If FBS ok set the specific parameters
if isFBS
    % -- General
    listOpti{nOpti}.name = 'Fwd-Bkwd Splitting';
    nameInLib{nOpti} = 'FBS';
    % - Gam parameter
    listOpti{nOpti}.parameters{1}.name = 'gam';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '0';
    listOpti{nOpti}.parameters{1}.default = '0';
    listOpti{nOpti}.parameters{1}.info = 'Descent step (Put 0 to let gam be computed automatically, see Doc to known when this is applicable)';
    listOpti{nOpti}.parameters{1}.toSet =1;
    % -- FISTA parameter
    listOpti{nOpti}.parameters{2}.name = 'fista';
    listOpti{nOpti}.parameters{2}.type = 'boolean';
    listOpti{nOpti}.parameters{2}.val = '0';
    listOpti{nOpti}.parameters{2}.default = '0';
    listOpti{nOpti}.parameters{2}.info = 'Boolean true if the accelerated version FISTA is used (default false)';
    listOpti{nOpti}.parameters{2}.toSet =1;
    % -- updateGam parameter
    listOpti{nOpti}.parameters{3}.name = 'updateGam';
    listOpti{nOpti}.parameters{3}.type = 'dropDown';
    listOpti{nOpti}.parameters{3}.items = {'none','reduced','backtracking'};
    listOpti{nOpti}.parameters{3}.val = 'none';
    listOpti{nOpti}.parameters{3}.default = 'none';
    listOpti{nOpti}.parameters{3}.info = 'Rule for updating gamma (''none'' : default, ''reduced'' : the parameter gam is decreased according to gam/sqrt(iter), ''backtracking'' : backtracking rule). Hint: a simple practical setting is to fix gam arbitrarily and to set the update rule to ''backtracking''';
    listOpti{nOpti}.parameters{3}.toSet =1;
    nOpti=nOpti+1;
end

%% VMLMB
% Requires that all costs are differentiables
if all(isdiff)
    listOpti{nOpti}.name = 'VMLMB';
    nameInLib{nOpti} = 'VMLMB';
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
    if doPrecom,  listOpti{nOpti}.call{1} = [ listOpti{nOpti}.call{1},' cf.doPrecomputation =1;']; end
    if isPos
        % if positivity
        listOpti{nOpti}.call{2} = 'Opt=OptiVMLMB(cf,0,Inf);';
        listOpti{nOpti}.parameters={};
    else
        % if not positivity
        listOpti{nOpti}.call{2} = 'Opt=OptiVMLMB(cf);';
        listOpti{nOpti}.parameters={};
    end
    nOpti=nOpti+1;
end


%% Richardson-Lucy
if (length(Costs)==1 && strcmp(Costs{1}.name,'CostKullLeib') && isPos)
    % Without regularization
    listOpti{nOpti}.name = 'Richardson-Lucy';
    nameInLib{nOpti} = 'RichLucy';
    listOpti{nOpti}.call{1} = ['Opt = OptiRichLucy(',GetExprCompCost(NameCosts,NamesOps,1),');'];
    listOpti{nOpti}.parameters={};
    nOpti=nOpti+1;
elseif (length(Costs)==2 && strcmp(Costs{1}.name,'CostKullLeib') && strcmp(CompCost{2}.cost2.H1.name,'CostHyperBolic') && strcmp(CompCost{2}.cost2.H2.name,'LinOpGrad')  && isPos)
    % With TV regularization
    listOpti{nOpti}.name = 'Richardson-Lucy';
    nameInLib{nOpti} = 'RichLucy';
    listOpti{nOpti}.call{1} = ['Opt = OptiRichLucy(',GetExprCompCost(NameCosts,NamesOps,1),',1,CostReg1.cost2.epsilon);'];
    listOpti{nOpti}.parameters={};
    nOpti=nOpti+1;
end

%% Chambolle-Pock
% Requires that no more than one cost is not prox
isCP=0;
if sum(isproxComp) >= nbCosts-1 
    if all(isproxComp)
        % If all costs are prox, search if the costs can be grouped in the
        % sum of two terms, one beeing prox and the second one like F(Kx)
        % with F prox
        for ii=1:nbCosts
            idx=isproxComp;
            idx(ii)=0;
            F=BuildSumCosts(CompCost,idx,isPos);
            isCP=(TestProx(F) && TestProx(Costs{ii}));
            if isCP, break; end                
        end
        % If ok, the Chambolle-Pock can be used
        if isCP
            listOpti{nOpti}.name = 'Chambolle-Pock';
            nameInLib{nOpti} = 'ChambPock';
            listOpti{nOpti}.call{1} = ['F = ',NameCosts{ii},';'];
            listOpti{nOpti}.call{2} = ['K = ',GetExprCompOp(NamesOps,ii),';'];
            if isPos
                listOpti{nOpti}.call{3} = ['G = ',GetExprCompCost(NameCosts,NamesOps,find(idx)),' + CostNonNeg(K.sizein);'];
            else
                listOpti{nOpti}.call{3} = ['G = ',GetExprCompCost(NameCosts,NamesOps,find(idx)),';'];
            end
            listOpti{nOpti}.call{4} = 'Opt = OptiChambPock(F,K,G);';
        end
    else
        % If one cost is not prox, try if there is the prox when decoupled
        % from the associated LinOp
        isCP = TestProx(Costs{~isproxComp});
        % If ok test if the sum of the remaining costs is prox
        if isCP
            F=BuildSumCosts(CompCost,isproxComp,isPos);
            if ~isempty(F)
                isCP=TestProx(F);
            else
                isCP=0;
            end
        end
        % If ok, the Chambolle-Pock can be used
        if isCP
            listOpti{nOpti}.name = 'Chambolle-Pock';
            nameInLib{nOpti} = 'ChambPock';
            listOpti{nOpti}.call{1} = ['F = ',NameCosts{~isproxComp},';'];
            listOpti{nOpti}.call{2} = ['K = ',GetExprCompOp(NamesOps,find(~isproxComp)),';'];
            tmp=GetExprCompCost(NameCosts,NamesOps,find(isproxComp)); if~isempty(tmp), tmp=[tmp,' + ']; end
            if isPos
                listOpti{nOpti}.call{3} = ['G = ',tmp,'CostNonNeg(K.sizein);'];
            else
                listOpti{nOpti}.call{3} = ['G = ',tmp,';'];
            end
            listOpti{nOpti}.call{4} = 'Opt = OptiChambPock(F,K,G);';
        end
    end
end
% If CP ok set the specific parameters
if isCP
    listOpti{nOpti}.parameters{1}.name='tau';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1';
    listOpti{nOpti}.parameters{1}.default = '1';
    listOpti{nOpti}.parameters{1}.info = 'Parameter of the algorithm (default 1). See the doc for convergence guarantees.';
    listOpti{nOpti}.parameters{1}.toSet =1;
    listOpti{nOpti}.parameters{2}.name='sig';
    listOpti{nOpti}.parameters{2}.type = 'double';
    listOpti{nOpti}.parameters{2}.val = '0';
    listOpti{nOpti}.parameters{2}.default = '0';
    listOpti{nOpti}.parameters{2}.info = 'Parameter of the algorithm which is computed automatically if the norm of the forward model is accessible. See the doc for convergence guarantees.';
    listOpti{nOpti}.parameters{2}.toSet =1;
    listOpti{nOpti}.parameters{3}.name='var';
    listOpti{nOpti}.parameters{3}.type = 'double';
    listOpti{nOpti}.parameters{3}.val = '1';
    listOpti{nOpti}.parameters{3}.default = '1';
    listOpti{nOpti}.parameters{3}.info = 'Select the "bar" variable of the algorithm (see [1])';
    listOpti{nOpti}.parameters{3}.toSet =1;
    nOpti=nOpti+1;
end

%% ADMM
if nbCosts+isPos>1
% Full splitting case
if all(isprox)
    listOpti{nOpti}.name = 'ADMM-FullSplit';
    nameInLib{nOpti} = 'ADMM';
    listOpti{nOpti}.call{2} = 'Fn = {';
    listOpti{nOpti}.call{1} = 'Hn = {';
    for ii=1:nbCosts
        % Combine Costs and Ops while prox
        tmp=Costs{ii};
        idx=0;
        for jj=1:length(Ops{ii})
            tmp=tmp*Ops{ii}{jj};
            if TestProx(tmp) && ~(strcmp(tmp.name,'CostTV') || (strcmp(tmp.name,'CostMultiplication') && strcmp(tmp.cost2.name,'CostTV')))
                idx=jj;                
            else
                break
            end
        end
        listOpti{nOpti}.call{2} =[listOpti{nOpti}.call{2},GetExprCompCost(NameCosts,NamesOps,ii,idx*ones(size(NamesOps)))];
        listOpti{nOpti}.call{1} =[listOpti{nOpti}.call{1},GetExprCompOp(NamesOps,ii,(idx+1)*ones(size(NamesOps)))];
        if ii<nbCosts
             listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},','];
             listOpti{nOpti}.call{1}=[ listOpti{nOpti}.call{1},','];
        end
    end
    if isPos
        listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},'}; Fn = [Fn,{CostNonNeg(Hn{1}.sizein)}];'];
        listOpti{nOpti}.call{1}=[ listOpti{nOpti}.call{1},'}; Hn = [Hn,{LinOpIdentity(Hn{1}.sizein)}];'];
        listOpti{nOpti}.costIndex=1:nbCosts;
    else
        listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},'};'];
        listOpti{nOpti}.call{1}=[ listOpti{nOpti}.call{1},'};'];
    end
    listOpti{nOpti}.call{3}='Opt = OptiADMM([],Fn,Hn,rho);';    
    listOpti{nOpti}.parameters{1}.name='rho';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1';
    listOpti{nOpti}.parameters{1}.default = '-1';
    listOpti{nOpti}.parameters{1}.toSet =0;
    listOpti{nOpti}.parameters{1}.info = 'Lagrangian Multiplier (positive real).';
    listOpti{nOpti}.parameters{2}.name='maxiterCG';
    listOpti{nOpti}.parameters{2}.type = 'double';
    listOpti{nOpti}.parameters{2}.val = '20';
    listOpti{nOpti}.parameters{2}.default = '20';
    listOpti{nOpti}.parameters{2}.info = 'Maximal number of conjugate gradient iterations (when required)';
    listOpti{nOpti}.parameters{2}.toSet =1;
    nOpti=nOpti+1;
end
% Do not split L2 terms
if any(isL2) && all(isL2+isprox)
    listOpti{nOpti}.name = 'ADMM-NoFullSplit';
    nameInLib{nOpti} = 'ADMM';
    listOpti{nOpti}.call{1} = ['F0 = ',GetExprCompCost(NameCosts,NamesOps,find(isL2)),';'];
    listOpti{nOpti}.call{2} = 'Fn = {';
    listOpti{nOpti}.call{3} = 'Hn = {';
    for ii=find(~isL2)
        % Combine Costs and Ops while prox
        tmp=Costs{ii};
        idx=0;
        for jj=1:length(Ops{ii})
            tmp=tmp*Ops{ii}{jj};
            if TestProx(tmp) && ~(strcmp(tmp.name,'CostTV') || (strcmp(tmp.name,'CostMultiplication') && strcmp(tmp.cost2.name,'CostTV')))
                idx=jj;
            else
                break
            end
        end
        listOpti{nOpti}.call{2} =[listOpti{nOpti}.call{2},GetExprCompCost(NameCosts,NamesOps,ii,idx*ones(size(NamesOps)))];
        listOpti{nOpti}.call{3} =[listOpti{nOpti}.call{3},GetExprCompOp(NamesOps,ii,(idx+1)*ones(size(NamesOps)))];
        if ii<find(~isL2,1,'last')
            listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},','];
            listOpti{nOpti}.call{3}=[ listOpti{nOpti}.call{3},','];
        end
    end
    if isPos
        listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},'}; Fn = [Fn,{CostNonNeg(F0.sizein)}];'];
        listOpti{nOpti}.call{3}=[ listOpti{nOpti}.call{3},'}; Hn = [Hn,{LinOpIdentity(F0.sizein)}];'];
        listOpti{nOpti}.costIndex=1:nbCosts;
    else
        listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},'};'];
        listOpti{nOpti}.call{3}=[ listOpti{nOpti}.call{3},'};'];
    end    
    listOpti{nOpti}.call{4}='Opt = OptiADMM(F0,Fn,Hn,rho);';
    listOpti{nOpti}.parameters{1}.name='rho';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1';
    listOpti{nOpti}.parameters{1}.default = '-1';
    listOpti{nOpti}.parameters{1}.toSet =0;
    listOpti{nOpti}.parameters{1}.info = 'Lagrangian Multiplier (positive real).';
    listOpti{nOpti}.parameters{2}.name='maxiterCG';
    listOpti{nOpti}.parameters{2}.type = 'double';
    listOpti{nOpti}.parameters{2}.val = '20';
    listOpti{nOpti}.parameters{2}.default = '20';
    listOpti{nOpti}.parameters{2}.info = 'Maximal number of conjugate gradient iterations (when required)';
    listOpti{nOpti}.parameters{2}.toSet =1;
    nOpti=nOpti+1;
end
end

%% Primal-Dual Condat
% Full splitting case
if all(isprox)
    listOpti{nOpti}.name = 'PrimalDualCondat-FullSplit';
    nameInLib{nOpti} = 'PrimalDualCondat';
    listOpti{nOpti}.call{2} = 'Fn = {';
    listOpti{nOpti}.call{1} = 'Hn = {';
    for ii=1:nbCosts
        % Combine Costs and Ops while prox
        tmp=Costs{ii};
        idx=0;
        for jj=1:length(Ops{ii})
            tmp=tmp*Ops{ii}{jj};
            if TestProx(tmp) && ~(strcmp(tmp.name,'CostTV') || (strcmp(tmp.name,'CostMultiplication') && strcmp(tmp.cost2.name,'CostTV')))
                idx=jj;                
            else
                break
            end
        end
        listOpti{nOpti}.call{2} =[listOpti{nOpti}.call{2},GetExprCompCost(NameCosts,NamesOps,ii,idx*ones(size(NamesOps)))];
        listOpti{nOpti}.call{1} =[listOpti{nOpti}.call{1},GetExprCompOp(NamesOps,ii,(idx+1)*ones(size(NamesOps)))];
        if ii<nbCosts
            listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},','];
            listOpti{nOpti}.call{1}=[ listOpti{nOpti}.call{1},','];
        end
    end
    listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},'};'];
    listOpti{nOpti}.call{1}=[ listOpti{nOpti}.call{1},'};'];
    if isPos
        listOpti{nOpti}.call{3}='P = CostNonNeg(Fn{1}.sizein);';
        listOpti{nOpti}.call{4}='Opt = OptiPrimalDualCondat([],P,Fn,Hn);';    
        listOpti{nOpti}.costIndex=2:nbCosts+1;
    else
        listOpti{nOpti}.call{3}='Opt = OptiPrimalDualCondat([],[],Fn,Hn);';    
    end  
    listOpti{nOpti}.parameters{1}.name='tau';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1e-1';
    listOpti{nOpti}.parameters{1}.default = '-1';
    listOpti{nOpti}.parameters{1}.toSet =1;
    listOpti{nOpti}.parameters{1}.info = 'Parameters tau, sig, and rho have to satisfy an inequality to ensure convergence (see Doc).';
    listOpti{nOpti}.parameters{2}.name='sig';
    listOpti{nOpti}.parameters{2}.type = 'double';
    listOpti{nOpti}.parameters{2}.val = '1e-1';
    listOpti{nOpti}.parameters{2}.default = '-1';
    listOpti{nOpti}.parameters{2}.toSet =1;
    listOpti{nOpti}.parameters{2}.info = 'Parameters tau, sig, and rho have to satisfy an inequality to ensure convergence (see Doc).';
    listOpti{nOpti}.parameters{3}.name='rho';
    listOpti{nOpti}.parameters{3}.type = 'double';
    listOpti{nOpti}.parameters{3}.val = '1.95';
    listOpti{nOpti}.parameters{3}.default = '1.95';
    listOpti{nOpti}.parameters{3}.toSet =1;
    listOpti{nOpti}.parameters{3}.info = 'Parameters tau, sig, and rho have to satisfy an inequality to ensure convergence (see Doc).';
    nOpti=nOpti+1;
end
% Do not split differentiable terms
if any(isdiff) && all(isdiff+isprox)
    listOpti{nOpti}.name = 'PrimalDualCondat-NoFullSplit';
    nameInLib{nOpti} = 'PrimalDualCondat';
    listOpti{nOpti}.call{1} = ['F0 = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
    listOpti{nOpti}.call{2} = 'Fn = {';
    listOpti{nOpti}.call{3} = 'Hn = {';
    for ii=find(~isdiff)
        % Combine Costs and Ops while prox
        tmp=Costs{ii};
        idx=0;
        for jj=1:length(Ops{ii})
            tmp=tmp*Ops{ii}{jj};
            if TestProx(tmp) && ~(strcmp(tmp.name,'CostTV') || (strcmp(tmp.name,'CostMultiplication') && strcmp(tmp.cost2.name,'CostTV')))
                idx=jj;
            else
                break
            end
        end
        listOpti{nOpti}.call{2} =[listOpti{nOpti}.call{2},GetExprCompCost(NameCosts,NamesOps,ii,idx*ones(size(NamesOps)))];
        listOpti{nOpti}.call{3} =[listOpti{nOpti}.call{3},GetExprCompOp(NamesOps,ii,(idx+1)*ones(size(NamesOps)))];
        if ii<find(~isdiff,1,'last')
            listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},','];
            listOpti{nOpti}.call{3}=[ listOpti{nOpti}.call{3},','];
        end
    end
    listOpti{nOpti}.call{2}=[ listOpti{nOpti}.call{2},'};'];
    listOpti{nOpti}.call{3}=[ listOpti{nOpti}.call{3},'};'];
    if isPos
        listOpti{nOpti}.call{4}='P = CostNonNeg(F0.sizein);';
        listOpti{nOpti}.call{5}='Opt = OptiPrimalDualCondat(F0,P,Fn,Hn);';
        listOpti{nOpti}.costIndex=[1:sum(isdiff),sum(isdiff)+2:nbCosts+1];
    else
        listOpti{nOpti}.call{4}='Opt = OptiPrimalDualCondat(F0,[],Fn,Hn);';
    end
    listOpti{nOpti}.parameters{1}.name='tau';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1e-1';
    listOpti{nOpti}.parameters{1}.default = '-1';
    listOpti{nOpti}.parameters{1}.toSet =1;
    listOpti{nOpti}.parameters{1}.info = 'Parameters tau, sig, and rho have to satisfy an inequality to ensure convergence (see Doc).';
    listOpti{nOpti}.parameters{2}.name='sig';
    listOpti{nOpti}.parameters{2}.type = 'double';
    listOpti{nOpti}.parameters{2}.val = '1e-1';
    listOpti{nOpti}.parameters{2}.default = '-1';
    listOpti{nOpti}.parameters{2}.toSet =1;
    listOpti{nOpti}.parameters{2}.info = 'Parameters tau, sig, and rho have to satisfy an inequality to ensure convergence (see Doc).';
    listOpti{nOpti}.parameters{3}.name='rho';
    listOpti{nOpti}.parameters{3}.type = 'double';
    listOpti{nOpti}.parameters{3}.val = '1.95';
    listOpti{nOpti}.parameters{3}.default = '1.95';
    listOpti{nOpti}.parameters{3}.toSet =1;
    listOpti{nOpti}.parameters{3}.info = 'Parameters tau, sig, and rho have to satisfy an inequality to ensure convergence (see Doc).';
    nOpti=nOpti+1;
end

%% FGP
if nbCosts==2 && isa(Costs{1},'CostL2') && isa(Ops{1}{1},'LinOpDiag') && isscalar(Ops{1}{1}.diag) && isa(CompCost{2}.cost2,'CostTV')
    listOpti{nOpti}.name = 'Fast Gradient Proximal';
    nameInLib{nOpti} = 'FGP';
    listOpti{nOpti}.call{1} = ['TV = ',NameCosts{2},'*',NamesOps{2}{1},';'];
    if isPos
        listOpti{nOpti}.call{2} = ['Opt = OptiFGP(',NameCosts{1},',TV,[0,Inf]);'];
    else
        listOpti{nOpti}.call{2} = ['Opt = OptiFGP(',NameCosts{1},',TV);'];
    end
    listOpti{nOpti}.parameters={};
    nOpti=nOpti+1;
end

%% Common parameters
for ii=1:length(listOpti)
    ll=length(listOpti{ii}.parameters);
    listOpti{ii}.parameters{ll+1}.name='maxiter';
    listOpti{ii}.parameters{ll+1}.type = 'double';
    listOpti{ii}.parameters{ll+1}.val = '50';
    listOpti{ii}.parameters{ll+1}.default = '50';
    listOpti{ii}.parameters{ll+1}.info = 'Stopping Criterion: Max number of algorithm iterations (default 50)';
    listOpti{ii}.parameters{ll+1}.toSet =1;
    listOpti{ii}.parameters{ll+2}.name='TolCost';
    listOpti{ii}.parameters{ll+2}.type = 'double';
    listOpti{ii}.parameters{ll+2}.val = '1e-4';
    listOpti{ii}.parameters{ll+2}.default = '-1';
    listOpti{ii}.parameters{ll+2}.info = 'Stopping Criterion: Tolerance on the difference between the cost value at two successive iterates.';
    listOpti{ii}.parameters{ll+2}.toSet =0;
    listOpti{ii}.parameters{ll+3}.name='TolStep';
    listOpti{ii}.parameters{ll+3}.type = 'double';
    listOpti{ii}.parameters{ll+3}.val = '1e-4';
    listOpti{ii}.parameters{ll+3}.default = '-1';
    listOpti{ii}.parameters{ll+3}.info = 'Stopping Criterion: Tolerance on the norm of the difference of two successive iterates.';
    listOpti{ii}.parameters{ll+3}.toSet =0;
end

warning on
end

function expr = GetExprCompCost(NameCosts,NamesOps,idxCosts,nbOps)
if nargin < 4 || isempty(nbOps)
    nbOps=cellfun(@(x) length(x),NamesOps);
end
expr = [];
for ii=idxCosts
    if nbOps(ii)==0
        expr = [expr,NameCosts{ii},' + '];
    else
        expr = [expr,NameCosts{ii},'*(',NamesOps{ii}{1}];
        for jj=2:nbOps(ii)
            expr = [expr,'*',NamesOps{ii}{jj}];
        end
        expr = [expr,') + '];
    end
end
expr=expr(1:end-3);
end

function expr = GetExprCompOp(NamesOps,index,firstOps)
if nargin < 3 || isempty(firstOps)
    firstOps=ones(size(NamesOps));
end
expr = [];
for ii=index
    if firstOps(ii) <= length(NamesOps{ii})
        expr = [expr,NamesOps{ii}{firstOps(ii)}];
        for jj=firstOps(ii)+1:length(NamesOps{ii})
            expr = [expr,'*',NamesOps{ii}{jj}];
        end
    else
        expr = ['LinOpIdentity(',NamesOps{ii}{end},'.sizein)'];
    end
end
end

function F=BuildSumCosts(Costs,idx,isPos)
global sizein
    F=[];
    for ii=find(idx)
        if isempty(F)
            F=Costs{ii};
        else
            F=F + Costs{ii};
        end
    end
    if isPos
        if isempty(F)
            F=CostNonNeg(sizein);
        else
            F=F+CostNonNeg(sizein);
        end
    end
end

function isproxComp=TestProx(F)
    try
        F.applyProx(rand(F.sizein),1);
        isproxComp=1;
    catch ME
        isproxComp=0;
    end
end

