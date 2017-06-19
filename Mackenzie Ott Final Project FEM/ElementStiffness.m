function ke = ElementStiffness(invJ, detJ, k)
si = [+0.77459 -0.77459 0];
ti = [+0.77459 -0.77459 0];
w = [.55555 .55555 .88888];
ke = zeros(4);

for i = 1:3
    for j = 1:3
        B = invJ*(1/4)*[1+ti(j) -1-ti(j) -1+ti(j) 1-ti(j); 1+si(i) 1-si(i) -1+si(i) -1-si(i)];
        x = B.'*k* B;
        ke = ke + w(i)*w(j)* x*detJ;
    end
end
%Test FUnction for second order Quadrature
% si=[.577 -.577]
% ti =[.577 -.577]
% B1 = invJ * (1/4)*[1+ti(1) -1-ti(1) -1+ti(1) 1-ti(1); 1+si(1) 1-si(1) -1+si(1) -1-si(1)]
% B2 = invJ * (1/4)*[1+ti(2) -1-ti(2) -1+ti(2) 1-ti(2); 1+si(1) 1-si(1) -1+si(1) -1-si(1)]
% B3 = invJ * (1/4)*[1+ti(1) -1-ti(1) -1+ti(1) 1-ti(1); 1+si(1) 1-si(1) -1+si(1) -1-si(1)]
% B4 = invJ * (1/4)*[1+ti(2) -1-ti(2) -1+ti(2) 1-ti(2); 1+si(2) 1-si(2) -1+si(2) -1-si(2)]
% ke2 = B1.'*k*B1*detJ + B2.'*k*B2*detJ + B3.'*k*B3*detJ +B4.'*k*B4*detJ

%B1test = B1.'*B1
%B2test = B2.'*B2
%B3test = B3.'*B3
%B4test = B4.'*B4
%k1test = B1test *detJ
%k2test = B2test *detJ
%k3test = B3test *detJ
%k4test = B4test *detJ