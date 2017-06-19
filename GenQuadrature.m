function quadpts = GenQuadrature(numQ)
%%%
% This function generates Quadrature points for numerical integration along
% each dimension, for interior and boundary integration separately. The
% output is scaled to itegration domain of [0,1] and respective weights.
% Inorder to use for the required domain size, new_quadratures =
% (QP(0,1)*domain_length + intitial_point), new_weights =
% (W(0,1)*domain_length)
%%%
global dim;

for i = 1:dim
    [quadpts(i).q(:,1) quadpts(i).q(:,2)] = lgwt(numQ(i), 0, 1);
end