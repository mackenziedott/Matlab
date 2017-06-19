clear all; clc; close all;

N = 1000;
X = linspace(-15,15,N);
al = -3;
be = .15;

%%% Integrand
gX = X.^2; 

myu = 0;
mu1 = -sqrt(-al/be);
mu2 = sqrt(-al/be);
Q = 10;
fX = exp(-(0.5*al*X.^2  + 0.25*be*X.^4)/2/Q);

%%% Find normalization constant and true integral
[xq wq] = lgwt(N, X(1), X(end));
gXq = xq.^2;
fXq = exp(-(0.5*al*xq.^2  + 0.25*be*xq.^4)/2/Q);
normz = fXq'*wq;
%gX1 = exp(-0.5*((X - mu1).^2/sig1^2));
%gX2 = exp(-0.5*((X - mu2).^2/sig1^2));
It = (gXq.*fXq/normz)'*wq;  %true value of integral

%%%Aceptance Rejection  Sampling
sig = 4;
c = 5; 
qX = exp(-0.5*((X - myu).^2/sig^2));

%c1 = 2.5;
plot(X, fX);
hold on
plot(X, c*qX, 'k');
plot(xq, fXq.*gXq, 'g');
set(gca, 'fontsize', 12, 'fontweight', 'bold')
legend('Target pdf (f(x))', 'Proposal (q(x))')
title(['Acceptance Rejection Sampling. c = ' num2str(c)])
% plot(X, c1*gX1, 'm');
% plot(X, c1*gX2, 'm');

ctr = 1;
num = 0;
arsam = [];
fS = [];
while num < N
    x = sig*randn(1) + myu;
    u = rand(1);
    qx = exp(-0.5*((x - myu).^2/sig^2));
    fx = exp(-(0.5*al*x.^2  + 0.25*be*x.^4)/2/Q);
    if u <= fx/c/qx;
        arsam = [arsam; x];
        fS = [fS; fx];
        num = num + 1;
    end
    ctr = ctr + 1;
end
fS = fS/normz;
fprintf('Total number of samples drawn: %d (acceptance rate = %1.3e.\n', ctr, (num/ctr));
figure(2)
[nn, xout] = hist(arsam, linspace(X(1), X(end), 25));
n = nn/sum(nn);
bar(xout, n/(xout(2) - xout(1)));
hold on
plot(X, fX/normz, 'r', 'MarkerFaceColor', 'r', 'linewidth', 2);
set(gca, 'fontsize', 12, 'fontweight', 'bold')
title('Acceptance Rejection Sampling: Comparison after normalization');
%hist(arsam, linspace(-15, 15, 20));
%%% Evaluation of integral:
gXar = arsam.^2;
Iar = sum(gXar)/N;
fprintf('\nTrue Value of Integral: %f\n', It);
fprintf('AR Approximation: %f\n', Iar);

%%%Importance Sampling
issig = 4;
ismu = -4;
isq = exp(-0.5*((X - ismu).^2/issig^2));
isnormz = sqrt(2*pi)*issig;
isq = isq/isnormz;

issam = issig*randn(N,1) + ismu;
fissam = exp(-(0.5*al*issam.^2  + 0.25*be*issam.^4)/2/Q)/normz;
qissam = exp(-0.5*((issam - ismu).^2/issig^2))/isnormz;
iswt = fissam./qissam;
%%% Integral Approximation
gXis = issam.^2;
Iis = gXis'*iswt/sum(iswt);
fprintf('IS Approximation: %f\n', Iis);

figure(3)
%plot(issam, fissam, 'ko', 'markerfacecolor', 'k');
plot(X, fX/normz);
hold on
plot(X, isq, 'k');
%plot(issam, qissam, 'mo', 'markerfacecolor', 'm');
set(gca, 'fontsize', 12, 'fontweight', 'bold')
legend('Target pdf (f(x))', 'Importance Density (q(x))')
title(['Importance Sampling: \sigma = ' num2str(issig)])
plot(issam, 0, 'ko', 'markersize', 4, 'markerfacecolor', 'k');
figure(4)
semilogy(issam, iswt, 'ko', 'markersize', 4, 'markerfacecolor', 'k');
set(gca, 'fontsize', 12, 'fontweight', 'bold')
title(['Importance Sampling: \sigma = ' num2str(issig)])
xlabel('Importance Samples (x_i)');
ylabel('Importance Weights (w_i)');


%%%%alternate ar
%     x1 = sig1*randn(1) + mu1;
%     x2 = sig1*randn(1) + mu2;
%     u1 = rand(1);
%     if u1 < 0.5
%         x = x1; 
%         qx = exp(-0.5*((x - mu1).^2/sig1^2));
%         c = c1;
%     else
%         x = x2;
%         qx = exp(-0.5*((x - mu2).^2/sig1^2));
%         c = c1;
%     end
