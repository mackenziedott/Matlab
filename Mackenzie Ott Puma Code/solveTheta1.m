function theta = solveTheta1(X7,Y7,S2, S6, S7,alpha71,a71)
A1 = S6*Y7-S7*sin(alpha71);
B1 = S6*X7+a71;
D1 = S2;
theta = TrigSolver(A1,B1,D1);
end
