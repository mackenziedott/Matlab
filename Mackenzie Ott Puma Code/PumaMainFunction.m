% Robot Geometry
% Puma

%%Inputs
S6 = input('S6 = ')
%test
% S6 = 4;
% tool_6 = [5; 3; 7];
% tool_F = [24.1120; 20.1130; 18.1670];
% tool_S6_F = [.0790; -.7870; .6120];
% tool_a67_F = [.9970; .0640; -.0470];

tool_6(1) = input('Tool_6 *996+X = ')
tool_6(2) = input('Tool_6 Y = ')
tool_6(3) = input('Tool_6 Z = ')
tool_F(1) = input('Tool_F X = ')
tool_F(2) = input('Tool_F Y = ')
tool_F(3) = input('Tool_F Z = ')
tool_S6_F(1) = input('Tool_S6 X = ')
tool_S6_F(2) = input('Tool_S6 Y = ')
tool_S6_F(3) = input('Tool_S6 Z = ')
tool_a67_F(1) = input('Tool_a67 X = ')
tool_a67_F(2) = input('Tool_a67 Y = ')
tool_a67_F(3) = input('Tool_a67 Z = ')
[phi1,theta2,theta3,theta4,theta5,theta6] = Puma(S6,tool_6,tool_F,tool_S6_F,tool_a67_F);
for i = 1:2
    for x = 1:4
        Soln(1,x+(4*(i-1))) = phi1(i);
    end
    for j = 1:2
        for x = 1:2
            Soln(3,x+2*(j-1)+4*(i-1)) = theta3(j+2*(i-1))*180/pi;
            Soln(2,x+2*(j-1)+4*(i-1)) = theta2(j+2*(i-1))*180/pi;
        end
   
    	for n = 1:2
            Soln(5, n+2*(j-1)+4*(i-1)) = theta5(n+2*(j-1)+4*(i-1))*180/pi;
            Soln(4, n+2*(j-1)+4*(i-1)) = theta4(n+2*(j-1)+4*(i-1))*180/pi;
            Soln(6, n+2*(j-1)+4*(i-1)) = theta6(n+2*(j-1)+4*(i-1))*180/pi;
        end
    end
end
%Determine Good Solutions
for i = 1:8
    for j = 1:6
        if isreal(Soln(j,i)) ~= 1
            Soln(j,:) = [];
            
        end
    end
end
g = size(Soln);
numGoodSoln = g(:,2)
Soln
% Soln(1,:)
% Soln(3,:)
