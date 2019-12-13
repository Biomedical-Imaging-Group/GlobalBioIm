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
end


%% Gradien Descent
if ~isPos && all(isdiff)
    listOpti{nOpti}.name = 'Gradient Descent';
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff))];
    listOpti{nOpti}.call{2} = 'Opt = OptiGradDsct(cf);';
    listOpti{nOpti}.parameters{1}.name = 'gam';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1e-1';
    listOpti{nOpti}.parameters{1}.default = '0';
    listOpti{nOpti}.parameters{1}.info = 'Descent step';
    nOpti=nOpti+1;
end

%% Forward-Backward Splitting
if all(isdiff)
    listOpti{nOpti}.name = 'Fwd-Bkwd Splitting';
    listOpti{nOpti}.call{1} = ['F = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff))];
    listOpti{nOpti}.call{2} = 'P = CostNonNeg(F.sizein);';
    listOpti{nOpti}.call{3} = 'Opt = OptiFBS(F,P);';
    listOpti{nOpti}.parameters{1}.name = 'gam';
    listOpti{nOpti}.parameters{1}.type = 'double';
    listOpti{nOpti}.parameters{1}.val = '1e-1';
    listOpti{nOpti}.parameters{1}.default = '0';
    listOpti{nOpti}.parameters{1}.info = 'Descent step';
    listOpti{nOpti}.parameters{2}.name = 'fista';
    listOpti{nOpti}.parameters{2}.type = 'boolean';
    listOpti{nOpti}.parameters{2}.val = '0';
    listOpti{nOpti}.parameters{2}.default = '0';
    listOpti{nOpti}.parameters{2}.info = 'Boolean true if the accelerated version FISTA is used (default false)';
    nOpti=nOpti+1;
end

%% VMLMB
if all(isdiff)
    listOpti{nOpti}.name = 'VMLMB';
    listOpti{nOpti}.call{1} = ['cf = ',GetExprCompCost(NameCosts,NamesOps,find(isdiff))];
    if isPos
        listOpti{nOpti}.call{2} = 'Opt=OptiVMLMB(cf,0,Inf);';
        listOpti{nOpti}.parameters={};
    else
        listOpti{nOpti}.call{2} = 'Opt=OptiVMLMB(cf);';
        listOpti{nOpti}.parameters={};
    end
end


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
    expr = [expr,');'];
end
end

