sz = [35, 75];
H = LinOpIdentity(sz);

%% apply
x = rand(sz);
y = H*x;
assert(isequal(x,y))

