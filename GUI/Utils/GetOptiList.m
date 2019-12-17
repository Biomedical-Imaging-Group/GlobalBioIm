function listOpti = GetOptiList(Costs,Ops,NameCosts,NamesOps,isPos)
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
end
nbCosts=length(Costs);


%% Gradien Descent
if ~isPos && all(isdiff)
    % -- General
    listOpti{nOpti}.name = 'Gradient Descent';
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
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
    listOpti{nOpti}.call{2} = 'P = CostNonNeg(F.sizein);';
    listOpti{nOpti}.call{3} = 'Opt = OptiFBS(F,P);';
    isFBS=1;
elseif sum(isdiff)>0
    % Otherwise: chck that the sum of the non-diff terms is prox
    F=BuildSumCosts(CompCost,~isdiff,isPos);   % Sum all non-differentiable costs 
    isFBS=TestProx(F);                         % Test if prox
    if isFBS
        listOpti{nOpti}.call{1} = ['F = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
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
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
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
    listOpti{nOpti}.call{1} = ['Opt = OptiRichLucy(',GetExprCompCost(NameCosts,NamesOps,1),');'];
    listOpti{nOpti}.parameters={};
    nOpti=nOpti+1;
elseif (length(Costs)==2 && strcmp(Costs{1}.name,'CostKullLeib') && strcmp(CompCost{2}.cost2.H1.name,'CostHyperBolic') && strcmp(CompCost{2}.cost2.H2.name,'LinOpGrad')  && isPos)
    % With TV regularization
    listOpti{nOpti}.name = 'Richardson-Lucy';
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
% Full splitting case
if all(isprox)
    listOpti{nOpti}.name = 'ADMM-FullSplit';
    listOpti{nOpti}.call{2} = 'Fn = {';
    listOpti{nOpti}.call{1} = 'Hn = {';
    for ii=1:nbCosts
        listOpti{nOpti}.call{2} =[listOpti{nOpti}.call{2},NameCosts{ii}];
        listOpti{nOpti}.call{1} =[listOpti{nOpti}.call{1},GetExprCompOp(NamesOps,ii)];
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
%     listOpti{nOpti}.parameters{2}.name='maxiterCG';
%     listOpti{nOpti}.parameters{2}.type = 'double';
%     listOpti{nOpti}.parameters{2}.val = '20';
%     listOpti{nOpti}.parameters{2}.default = '0';
%     listOpti{nOpti}.parameters{2}.info = '';
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

function expr = GetExprCompCost(NameCosts,NamesOps,index)
expr = [];
for ii=index
    expr = [expr,NameCosts{ii},'*(',NamesOps{ii}{1}];
    for jj=2:length(NamesOps{ii})
        expr = [expr,'*',NamesOps{ii}{jj}];
    end
    expr = [expr,') + '];
end
expr=expr(1:end-3);
end

function expr = GetExprCompOp(NamesOps,index)
expr = [];
for ii=index
    expr = [expr,NamesOps{ii}{1}];
    for jj=2:length(NamesOps{ii})
        expr = [expr,'*',NamesOps{ii}{jj}];
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

