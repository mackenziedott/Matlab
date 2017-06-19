function x = Map(A, B, C, D, s, t)
x = (1/4)*((A+B+C+D)+(A-B-C+D)*s+(A+B-C-D)*t + (A-B+C-D)*s*t);
end