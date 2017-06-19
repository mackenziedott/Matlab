function target = targetevalmesh(xm, ym)
% global N;


%% Rayleigh
% sig = 2;
% target = x.*exp(-x.^2/2/sig^2)/sig^2;

%% May-West
target = 10./((ym - xm.^2).^2 + (1 - xm).^2 + 1) + 5./((ym - 8).^2 + (5 - xm).^2 + 1) + 5./((ym - 8).^2 + (8 - xm).^2 + 1);

%% Rational Poly- simple
% target = 1./(1 + x(:,1).^2 + x(:,2).^2);

%% Loss function
% target = abs(0.05*(xm.^2 + ym.^2) - 40*cos(xm).*cos(ym));
% target = 0.5*((xm.^4 - 16*xm.^2 + 5*xm) + (ym.^4 - 16*ym.^2 + 5*ym)) - 10*(cos(4*(xm + 2.093534)).*cos(2*(ym + 2.903534)));