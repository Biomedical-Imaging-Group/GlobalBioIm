%% Sums should flatten
A = Diagonal( 3 * ones(5) );
B = Identity([5, 5]);
C = Scaling(2, [5, 5]);

h = zeros(5,5);
h(1:2) = [.5 .5];
D = Convolution(h);

sumABCD = A + B + C + D;
assert( length(sumABCD.ALinOp) == 4 );
assert( all(all( sumABCD * ones(5) == 7 )) );

%% we should call HtH and HHt in MulLinOp

A = Diagonal(2 * ones(10) );
AtA = A' * A;
AAt = A * A';

profile on
AtA * rand(10);
AAt * rand(10);
stats = profile('info');

assert( any(ismember({stats.FunctionTable.FunctionName}, 'LinOp>LinOp.HtH')))
assert( any(ismember({stats.FunctionTable.FunctionName}, 'LinOp>LinOp.HHt')))

%% but only where appropriate

A = Diagonal(2 * ones(10) );
B = Identity( [10, 10] );
AtB = A' * B;
ABt = A * B';

profile on
ABt * rand(10);
ABt * rand(10);
stats = profile('info');

assert( ~any(ismember({stats.FunctionTable.FunctionName}, 'LinOp>LinOp.HtH')))
assert( ~any(ismember({stats.FunctionTable.FunctionName}, 'LinOp>LinOp.HHt')))

