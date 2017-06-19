function sample = getsample(current, num)

global N;
psig = 3*eye(N);
R = chol(psig);
sample = repmat(current, num,1) + randn(num,N)*R;