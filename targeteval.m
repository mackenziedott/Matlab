function target = targeteval(x)
global N Dom;


%% Rayleigh
% sig = 2;
% target = x.*exp(-x.^2/2/sig^2)/sig^2;

%% May-West
target = 10./((x(:,2) - x(:,1).^2).^2 + (1 - x(:,1)).^2 + 1) + 5./((x(:,2) - 8).^2 + (5 - x(:,1)).^2 + 1) + 5./((x(:,2) - 8).^2 + (8 - x(:,1)).^2 + 1);

%% Rational Poly- simple
% target = 1./(1 + x(:,1).^2 + x(:,2).^2);

%% Loss function
% xlen = size(x,1);
% indom = inbox(Dom, x);
% target = zeros(xlen,1);
% if isempty(indom)
%     return;
% end
% loss1 = zeros(xlen,1);
% loss2 = ones(xlen,1);

% for i = 1:N
%     loss1 = loss1 + 0.05*x(:,i).^2;
%     loss2 = loss2.*cos(x(:,i));
% end
% target = abs(loss1 - 40*loss2);
% for i = 1:N
%     loss1 = loss1 + 0.05*x(:,i).^2;
%     loss2 = loss2.*cos(x(:,i));
% %     loss1(indom) = loss1(indom) + 0.5*(x(indom,i).^4 - 16*x(indom,i).^2 + 5*x(indom,i));
% end
% loss2(indom) = loss2(indom).*(cos(4*(x(indom,1) + 2.093534)).*cos(2*(x(indom,2) + 2.903534)));
% target(indom) = abs(loss1(indom) - 40*loss2(indom));