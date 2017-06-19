clear all

clc
close all

show_simulation = 0 ;

%read the data from file
load('BallData_forclass.mat','x_hist','z_hist','-mat');
%z_hist is the measurements
%choosing x0

for i = 1:length(z_hist)
    z(:,i) = [z_hist(1,i); z_hist(2,i)];
end
T = 1; %arbitrary unit of time
time = [0:length(z_hist)-1]'*T
max_time_index = length(time)

%Initializa System
%v(t) = Ac*(x(t)+Dc*w(t)
%zk = C*Xk+thethak)
%%%%DESIGN VARIABLES%%%%%
m0 = [0;0;0;0] %Stated to be initialized as a 0 vectors
W = .001*[.0585 .005; .005 .046]
%   W = .05*eye(2)
%  W = .1*[.0025 .01*1/2; .01*1/2 1]*T %Model Noice
Q =  [1 .5; .5 2] %Meausrement Noise
%  Q =  .010*[1 .5; .5 2] 
%defined as the cov(x0,x0)
Sig = 80 * eye(4) 
%%%MODEL CONSTANTS
Ac = [0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 0]
%FOR THIS SYSTEM THE FOLLOWING ARE CONSIDERED 0 FOR CONSTANT VELOCITY
B = 0
uk = 0
%
Dc = [0 0;
    1 0;
     0 0;
     0 1]
C = [1 0 0 0;
    0 0 1 0]
A = eye(4) + Ac*T
D = A*Dc
% P = dlyap(A,D*W*D')
%Assume State noise is zero R

% W = eye(2)
%Assume Measurement Noise is std 1
%Q = eye(2)
% Theta = 1
%Can be ajusted
%GENERATE WHITE NOISE
% w = randn(2,10)*W
% Q = randn(2,10)*Theta
% [X,L,G] = dare(A,C', D*W*D', Q)
% X
% L
% G
% H = W*T
%Sig_gk = (A*Dc*W*Dc'*A')*T
%Initialize k=0
% x(:,1) = x0


% xhat(:,1) = [z_hist(1,1); 0; z_hist(2,1); 0]
% xbar(:,1) = x-xhat
% Test1 = C*Sig*C'
M0 = Sig*C'*(C*Sig*C'+Q)^-1
% z(:,1)
% C*x0
% Test2 = z(1)-C*x0
% xhat_k|k-1 = m0
xhat_kk(:,1) = m0+M0*(z(:,1)-C*m0)
%Execution of FIlter
%k = 1 represenets k = 0 in code
for k = 1:max_time_index-1
%     t(k)= (k-1)*T
    

    %Not needed for this assignment
    %x_new = A*x(k)+B*uk+D*w(k)
%     x_new = A*x(k)+ D*w(k)
%     %Since Z is provided, the following is not needed
%     %Z_new = C*x(k)+Theta(k)
%     L = A*Sig*C'*(C*Sig*C'+Q)^-1
%COMPUTE P
% L = A*Sig*C'*inv(C*Sig*C'+Q)
%NOT NEEDED FOR STATE ESTIMATION - Below are used for 1 Step-Predictors
% xbar = x(:,k)-xhat_k|k1
% xhat_k1|k = A*xhat_k|k1 + L*(z(:,k)-C*xhat_k|k1
%     P_new = A*P*A'+Q
%     K=P_new*C'*inv(C*P*C'

     Sig_k1 = A*(Sig-Sig*C'*((C*Sig*C'+Q)^-1)*C*Sig)*A'+D*W*D'
     %below is not needed for 
%     xhat_new = A*xhat+L*(z_hist(:,k)-C*xhat)
%     xbar_new = x_new-xhat_new
    M_new = Sig_k1*C'*(C*Sig_k1*C'+Q)^-1
    xhat_kk(:,k+1) = A*xhat_kk(:,k)+M_new*(z(:,k+1)-C*(A*xhat_kk(:,k)))
    %This is Most Important
%     xhat_kk(:,k) = x(k-1)+M*(Eta-C*x(k))
    Sig = Sig_k1
    
%     x(:,k+1) = x_new
%     xbar(:,k+1) = xbar_new
%     Sig(k+1) = Sig_new
%     xhat(:,k+1) = xhat_new
%     xhat_kk_plot(:,k) = xhat_kk
end
%Truncate xhat_kk_velocity
xhat_kk_v=xhat_kk
% xhat_kk_v(:,11) = []
% figure;
% ph = plot(time, x_hist(1,:),'r.-',time,z_hist(1,:),'bo-',time,xhat_kk(1,:)); hold on;
% set(gca,'fontsize',18);
% legend('true','measured x-pos', 'xhat k|k');
% xlabel('time');
% ylabel(' x (distance unit)');
% axis([0 10 0 10]);
% ss = sprintf('actual vs. observed x-pos \n (we have actual only because it''s a simulation)');
% title(ss)
% % 
% figure;
% ph = plot(time, x_hist(3,:),'r.-',time,z_hist(2,:),'bo-',time,xhat_kk(3,:)); hold on;
% set(gca,'fontsize',18);
% legend('true','measured y-pos', 'xhat k|k');
% xlabel('time');
% ylabel(' x (distance unit)');
% axis([0 10 0 10]);
% ss = sprintf('actual vs. observed y-pos \n (we have actual only because it''s a simulation)');
% title(ss)
figure;
ph = plot(x_hist(1,:),x_hist(3,:),'ro:', z_hist(1,:),z_hist(2,:),'bo-', xhat_kk(1,:), xhat_kk(3,:),'k*-.', 'MarkerSize', 15)
set(gca,'fontsize',18);
legend('true position','measured position', 'xhat k|k');
xlabel( 'x (distance unit)');
ylabel(' y (distance unit)');
axis([0 10 0 10]);
figure;
subplot(2,1,1)
ph = plot(time,x_hist(2,:),'r.-',time,xhat_kk_v(2,:),'bo-')
set(gca,'fontsize',18);
legend('true','measured x-vel');
xlabel( 'time');
ylabel(' (distance unit/time)');
axis([0 10 0 1.5]);
subplot(2,1,2)
ph2 = plot(time,x_hist(4,:),'r.-',time,xhat_kk_v(4,:),'bo-')
set(gca,'fontsize',18);
legend('true','measured y-vel');
xlabel( 'time');
ylabel(' (distance unit/time)');
axis([0 10 0 1]);