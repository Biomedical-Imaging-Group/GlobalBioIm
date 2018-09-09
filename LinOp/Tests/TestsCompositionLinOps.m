%% Test LinOp Compositions

isOk=true;
% Define some operators
sz = [128,128];
S=LinOpSelectorPatch(sz,[1 1],[64 64]);
I64=LinOpIdentity(sz/2);
I128=LinOpIdentity(sz);
H=LinOpConv(rand(sz/2));
fprintf('\nOperators used : \n \t - S : LinOpSelectorPatch \n \t - I : LinOpIdentity (LinOpDiag with a diagonal of 1) \n \t - H : LinOpConv \n \n');
me=[];

%% Test #1 (Scalar multiplication)
fprintf('Test #1 (Scalar multiplication) : \n ');
fprintf('\t 2*S*S\''  =2I ');
if isequal(2*S*S',2*I64)
    fprintf('   ---> OK \n');
else
    fprintf('   ---> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 2*S*S\''  ~= 2I');
end
fprintf('\t 2*(S*S\'')=2I ');
if isequal(2*(S*S'),2*I64)
    fprintf('   ---> OK \n');
else
    fprintf('   ---> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 2*(S*S\'') ~= 2I');
end
fprintf('\t S*(2*S\'')=2I ');
if isequal(S*(2*S'),2*I64)
    fprintf('   ---> OK \n');
else
    fprintf('   ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' S*(2*S\'') ~= 2I');
end
fprintf('\t 4*(2*S)*S\''=8I ');
if isequal(4*(2*S)*S',8*I64)
    fprintf(' ---> OK \n');
else
    fprintf(' ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 4*(2*S)*S\'' ~= 8I');
end
fprintf('\t 4*S*(2*S\'')=8I ');
if isequal(4*S*(2*S'),8*I64)
    fprintf(' ---> OK \n');
else
    fprintf(' ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 4*S*(2*S\'') ~= 8I');
end

%% Test #2 (More complex compositions)
fprintf('\nTest #2 (More complex compositions) : \n ');
fprintf('\t S\''*(S*S\'')=S\'' ');
if isequal(S'*(S*S'),S')
    fprintf('                    ---> OK \n');
else
    fprintf('                    ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' S\''*(S*S\'') ~= S\'' ');
end
fprintf('\t (H\''*S)*(S\''*H)=H\''*H ');
if isequal((H'*S)*(S'*H),H'*H)
    fprintf('              ---> OK \n');
else
    fprintf('              ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (H\''*S)*(S\''*H) ~= H\''*H ');
end
fprintf('\t H\''*S*S\''*H=H\''*H ');
if isequal(H'*S*S'*H,H'*H)
    fprintf('                  ---> OK \n');
else
    fprintf('                  ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' H\''*S*S\''*H ~= H\''*H ');
end
fprintf('\t 2*H\''*(2*S)*(3*S\'')*(1*H)=12*H\''*H ');
if isequal(2*H'*(2*S)*(3*S')*(1*H),12*H'*H)
    fprintf(' ---> OK \n');
else
    fprintf(' ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 2*H\''*(2*S)*(3*S\'')*(1*H) ~= 12*H\''*H ');
end
fprintf('\t H\''*(S*S\'')*H=H\''*H ');
if isequal(H'*(S*S')*H,H'*H)
    fprintf('                ---> OK \n');
else
    fprintf('                ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' H\''*(S*S\'')*H ~= H\''*H ');
end

global generateDataUnitTests stateTest message
if ~isempty(generateDataUnitTests)
    stateTest=isOk;
    message=me;
end