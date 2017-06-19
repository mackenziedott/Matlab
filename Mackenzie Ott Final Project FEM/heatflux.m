function q = heatflux(T, si, ti, invJ, k)
B = invJ*(1/4)*[1+ti -1-ti -1+ti 1-ti; 1+si 1-si -1+si -1-si]; %used to check work
invJ; % used to check work
x=[T(3);T(1); T(2); T(4)]; % used to check work
z = -k*B;
q = z*x;
end