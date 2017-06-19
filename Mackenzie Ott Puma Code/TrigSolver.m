function theta = TrigSolver(A,B,D)
ys(1) = B/sqrt(A^2+B^2);

yc(1) = A/sqrt(A^2+B^2);
y = atan2(ys,yc);

theta(1) = y+acos(-D/(sqrt(A^2+B^2)));
theta(2) = y-acos(-D/(sqrt(A^2+B^2)));
end