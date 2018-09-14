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
Op1=2*S*S';Op1.clearListenerList();
Op2=2*I64;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf('   ---> OK \n');
else
    fprintf('   ---> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 2*S*S\''  ~= 2I');
end
fprintf('\t 2*(S*S\'')=2I ');
Op1=2*(S*S');Op1.clearListenerList();
Op2=2*I64;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf('   ---> OK \n');
else
    fprintf('   ---> FAILED \n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 2*(S*S\'') ~= 2I');
end
fprintf('\t S*(2*S\'')=2I ');
Op1=S*(2*S');Op1.clearListenerList();
Op2=2*I64;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf('   ---> OK \n');
else
    fprintf('   ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' S*(2*S\'') ~= 2I');
end
fprintf('\t 4*(2*S)*S\''=8I ');
Op1=4*(2*S)*S';Op1.clearListenerList();
Op2=8*I64;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf(' ---> OK \n');
else
    fprintf(' ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 4*(2*S)*S\'' ~= 8I');
end
fprintf('\t 4*S*(2*S\'')=8I ');
Op1=4*S*(2*S');Op1.clearListenerList();
Op2=8*I64;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf(' ---> OK \n');
else
    fprintf(' ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 4*S*(2*S\'') ~= 8I');
end

%% Test #2 (More complex compositions)
fprintf('\nTest #2 (More complex compositions) : \n ');
fprintf('\t S\''*(S*S\'')=S\'' ');
Op1=S'*(S*S');Op1.clearListenerList();
Op2=S';Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf('                    ---> OK \n');
else
    fprintf('                    ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' S\''*(S*S\'') ~= S\'' ');
end
fprintf('\t (H\''*S)*(S\''*H)=H\''*H ');
Op1=(H'*S)*(S'*H);Op1.clearListenerList();
Op2=H'*H;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf('              ---> OK \n');
else
    fprintf('              ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' (H\''*S)*(S\''*H) ~= H\''*H ');
end
fprintf('\t H\''*S*S\''*H=H\''*H ');
Op1=H'*S*S'*H;Op1.clearListenerList();
Op2=H'*H;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf('                  ---> OK \n');
else
    fprintf('                  ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' H\''*S*S\''*H ~= H\''*H ');
end
fprintf('\t 2*H\''*(2*S)*(3*S\'')*(1*H)=12*H\''*H ');
Op1=2*H'*(2*S)*(3*S')*(1*H);Op1.clearListenerList();
Op2=12*H'*H;Op2.clearListenerList();
if isequal(Op1,Op2)
    fprintf(' ---> OK \n');
else
    fprintf(' ---> FAILED\n');isOk=false;
    me=sprintf('%s %s \n \t ',me,' 2*H\''*(2*S)*(3*S\'')*(1*H) ~= 12*H\''*H ');
end
fprintf('\t H\''*(S*S\'')*H=H\''*H ');
Op1=H'*(S*S')*H;Op1.clearListenerList();
Op2=H'*H;Op2.clearListenerList();
if isequal(Op1,Op2)
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