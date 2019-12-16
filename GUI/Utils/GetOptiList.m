function listOpti = GetOptiList(Costs,Ops,NameCosts,NamesOps,isPos)
    
% Initialize listOpti
nOpti=1;
listOpti={};

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
    isprox(ii) = 1;
    try
        CompCost{ii}.applyProx(rand(CompCost{ii}.sizein),1);
    catch ME
        isprox(ii) = 0;
    end
end


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
    nOpti=nOpti+1;
end

%% Forward-Backward Splitting
isFBS=0;
if all(isdiff) && isPos
    listOpti{nOpti}.call{1} = ['F = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
    listOpti{nOpti}.call{2} = 'P = CostNonNeg(F.sizein);';
    listOpti{nOpti}.call{3} = 'Opt = OptiFBS(F,P);';
    isFBS=1;
else
    F=[];
    for ii=find(~isdiff)
        if isempty(F)
            F=CompCost{ii};
        else
            F=F + CompCost{ii};
        end
    end
    if isPos, F=F+CostNonNeg(F.sizein); end
    % Test if prox
    isFBS=1;
    try
        F.applyProx(rand(F.sizein),1);
    catch ME
        isFBS=0;
    end
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
end
if isFBS
    % -- General
    listOpti{nOpti}.name = 'Fwd-Bkwd Splitting';
    % - Gam parameter
    listOpti{nOpti}.parameters{1}.name = 'gam';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '0';
    listOpti{nOpti}.parameters{1}.default = '0';
    listOpti{nOpti}.parameters{1}.info = 'Descent step (Put 0 to let gam be computed automatically, see Doc to known when this is applicable)';
    % -- FISTA parameter
    listOpti{nOpti}.parameters{2}.name = 'fista';
    listOpti{nOpti}.parameters{2}.type = 'boolean';
    listOpti{nOpti}.parameters{2}.val = '0';
    listOpti{nOpti}.parameters{2}.default = '0';
    listOpti{nOpti}.parameters{2}.info = 'Boolean true if the accelerated version FISTA is used (default false)';
    % -- updateGam parameter
    listOpti{nOpti}.parameters{3}.name = 'updateGam';
    listOpti{nOpti}.parameters{3}.type = 'dropDown';
    listOpti{nOpti}.parameters{3}.items = {'none','reduced','backtracking'};
    listOpti{nOpti}.parameters{3}.val = 'none';
    listOpti{nOpti}.parameters{3}.default = 'none';
    listOpti{nOpti}.parameters{3}.info = 'Rule for updating gamma (''none'' : default, ''reduced'' : the parameter gam is decreased according to gam/sqrt(iter), ''backtracking'' : backtracking rule). Hint: a simple practical setting is to fix gam arbitrarily and to set the update rule to ''backtracking''';
    nOpti=nOpti+1;
end

%% VMLMB
if all(isdiff)
    listOpti{nOpti}.name = 'VMLMB';
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff)),';'];
    if isPos
        listOpti{nOpti}.call{2} = 'Opt=OptiVMLMB(cf,0,Inf);';
        listOpti{nOpti}.parameters={};
    else
        listOpti{nOpti}.call{2} = 'Opt=OptiVMLMB(cf);';
        listOpti{nOpti}.parameters={};
    end
    nOpti=nOpti+1;
end


%% Richardson-Lucy
if (length(Costs)==1 && strcmp(Costs{1}.name,'CostKullLeib') && isPos)
    listOpti{nOpti}.name = 'Richardson-Lucy';
    listOpti{nOpti}.call{1} = ['Opt = OptiRichLucy(',GetExprCompCost(NameCosts,NamesOps,1),');'];
    listOpti{nOpti}.parameters={};
    nOpti=nOpti+1;
elseif (length(Costs)==2 && strcmp(Costs{1}.name,'CostKullLeib') && strcmp(CompCost{2}.cost2.H1.name,'CostHyperBolic') && strcmp(CompCost{2}.cost2.H2.name,'LinOpGrad')  && isPos)
    listOpti{nOpti}.name = 'Richardson-Lucy';
    listOpti{nOpti}.call{1} = ['Opt = OptiRichLucy(',GetExprCompCost(NameCosts,NamesOps,1),',1,CostReg1.cost2.epsilon);'];
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
    listOpti{ii}.parameters{ll+2}.name='TolCost';
    listOpti{ii}.parameters{ll+2}.type = 'double';
    listOpti{ii}.parameters{ll+2}.val = '1e-4';
    listOpti{ii}.parameters{ll+2}.default = '-1';
    listOpti{ii}.parameters{ll+2}.info = 'Stopping Criterion: Tolerance on the difference between the cost value at two successive iterates.';
    listOpti{ii}.parameters{ll+3}.name='TolStep';
    listOpti{ii}.parameters{ll+3}.type = 'double';
    listOpti{ii}.parameters{ll+3}.val = '1e-4';
    listOpti{ii}.parameters{ll+3}.default = '-1';
    listOpti{ii}.parameters{ll+3}.info = 'Stopping Criterion: Tolerance on the norm of the difference of two successive iterates.';
end

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

