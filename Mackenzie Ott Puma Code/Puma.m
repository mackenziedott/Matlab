function [phi1,theta2,theta3,theta4,theta5,theta6] = Puma(S6, tool_6, tool_F, tool_S6_F, tool_a67_F)
%Constant Mechanism Parameters
a12 = 0;
a23 = 17;
a34 = 0.8;
a45 =0;
a56 =0;
alpha12 = pi/2;
alpha23 = 0;
alpha34 = 3*pi/2;
alpha45 = pi/2;
alpha56 = pi/2;
S2 = 5.9;
S3 = 0;
S4 = 17;
S5 = 0;
%Arbiterally defnie seventh joint axis using close the loop
a67 = 0;
alpha67 = pi/2; %degrees
%close-the-loop
[theta7, y1, alpha71, S7, S1, a71] = closeTheLoop(tool_6,tool_F, tool_S6_F, tool_a67_F);
%Find Theta 1++++++++++

Y7 = -(sin(alpha71)*cos(alpha67)+cos(alpha71)*sin(alpha67)*cos(theta7));
X7 = sin(alpha67)*sin(theta7)
Z7 = cos(alpha71)*cos(alpha67)-sin(alpha67)*sin(alpha71)*cos(theta7);
theta1 = solveTheta1(X7,Y7,S2, S6,S7, alpha71,a71)
phi1 = (180)/pi*[theta1(1)-y1; theta1(2)-y1];

%Find Theta 3
for i = 1:2
     X1(i) =sin(alpha71)*sin(theta1(i));
%      X1
     Y1(i) = -(sin(alpha12)*cos(alpha71)+cos(alpha12)*sin(alpha71)*cos(theta1(i)));
%      Y1
     Z1(i) = cos(alpha12)*cos(alpha71)-sin(alpha12)*sin(alpha71)*cos(theta1(i));
     X71(i) =X7*cos(theta1(i))-Y7*sin(theta1(i));
%      X71
     Y71(i) = cos(alpha12)*(X7*sin(theta1(i))+Y7*cos(theta1(i)))-sin(alpha12)*Z7;
%      Y71
     Z71(i) = sin(alpha12)*(X7*sin(theta1(i))+Y7*cos(theta1(i)))+cos(alpha12)*Z7;
%      Z71
     Q3(i) = -S6*X71(i)-S7*X1(i)-a71*cos(theta1(i));
%      Q3
     W3(i) = S1-S6*Y71(i)-S7*Y1(i);
%      W3
     A3(i) = 2*a23*a34;
%      A3
     B3(i) = -2*a23*S4;
%      B3
     D3(i) = a23^2 + a34^2 +S4^2 - Q3(i)^2-W3(i)^2;
%      D3
     thetahold3 = TrigSolver(A3(i),B3(i),D3(i));
    theta3(1+2*(i-1)) = thetahold3(1)
      theta3(2+2*(i-1)) = thetahold3(2)
     for j = 1:2
         %Find theta 2 
        A2a(j+2*(i-1)) = a23+a34*cos(theta3(j+2*(i-1)))-S4*sin(theta3(j+2*(i-1)));
        B2a(j+2*(i-1)) = -a34*sin(theta3(j+2*(i-1))) - S4*cos(theta3(j+2*(i-1)))
        A2b(j+2*(i-1)) = -a34*sin(theta3(j+2*(i-1))) - S4*cos(theta3(j+2*(i-1)))
        B2b(j+2*(i-1)) = -a23-a34*cos(theta3(j+2*(i-1)))+S4*sin(theta3(j+2*(i-1)));
        R = [A2a(j+2*(i-1)) B2a(j+2*(i-1));
            A2b(j+2*(i-1)) B2b(j+2*(i-1))]
        Q = [Q3(i); W3(i)]
        % E = [cos theta2; sin theta]
        E = R\Q
        theta2(j+2*(i-1))=atan2(E(2),E(1))
        %Find Theta 5
        X712(j+2*(i-1)) = X71(i)*cos(theta2(j+2*(i-1)))-Y71(i)*sin(theta2((j+2*(i-1))));
        
        Y712(j+2*(i-1)) = cos(alpha23)*(X71(i)*sin(theta2(j+2*(i-1)))+Y71(i)*cos(theta2(j+2*(i-1)))) - sin(alpha23)*Z71(i);
        
        Z712(j+2*(i-1)) = sin(alpha23)*(X71(i)*sin(theta2(j+2*(i-1)))+Y71(i)*cos(theta2(j+2*(i-1)))) + cos(alpha23)*Z71(i);
        
        X7123(j+2*(i-1)) = X712(j+2*(i-1))*cos(theta3(j+2*(i-1)))-Y712(j+2*(i-1))*sin(theta3(j+2*(i-1)));
        
        Y7123(j+2*(i-1)) = cos(alpha34)*(X712(j+2*(i-1))*sin(theta3(j+2*(i-1)))+Y712(j+2*(i-1))*cos(theta3(j+2*(i-1))))-sin(alpha34)*Z712(j+2*(i-1));
        
        Z7123(j+2*(i-1)) = sin(alpha34)*(X712(j+2*(i-1))*sin(theta3(j+2*(i-1)))+Y712(j+2*(i-1))*cos(theta3(j+2*(i-1))))+cos(alpha34)*Z712(j+2*(i-1));
        
        thetamap5(1) = acos(-Z7123(j+2*(i-1)));
        
        thetamap5(2) = -thetamap5(1);
        theta5(1+2*(j-1)+4*(i-1) ) = thetamap5(1);
        theta5(2+2*(j-1)+4*(i-1) ) = thetamap5(2);
    
         for n=1:2
            %find theta 4
            theta4c(1) = X7123(j+2*(i-1))/sin(theta5(n+2*(j-1)+4*(i-1) ));
        
            theta4s(1) = -Y7123(j+2*(i-1))/sin(theta5(n+2*(j-1)+4*(i-1) ));
            theta4(n+2*(j-1)+4*(i-1) ) = atan2(theta4s,theta4c)
            
            %Find Theta 6
            %8 Theta 4
            Xbar4 =   sin(alpha45)*sin(theta4(n+2*(j-1)+4*(i-1)))
            Ybar4 = -(sin(alpha34)*cos(alpha45)+cos(alpha34)*sin(alpha45)*cos(theta4(n+2*(j-1)+4*(i-1))))
            Zbar4 =   cos(alpha34)*cos(alpha45)-sin(alpha34)*sin(alpha45)*cos(theta4(n+2*(j-1)+4*(i-1)))
            %4 Theta 3
            X43 = Xbar4*cos(theta3(j+2*(i-1)))-Ybar4*sin(theta3(j+2*(i-1)))
            Y43 = cos(alpha23)*(Xbar4*sin(theta3(j+2*(i-1)))+Ybar4*cos(theta3(j+2*(i-1))))-sin(alpha23)*Zbar4
            Z43 = sin(alpha23)*(Xbar4*sin(theta3(j+2*(i-1)))+Ybar4*cos(theta3(j+2*(i-1))))+cos(alpha23)*Zbar4
            %4 Theta2
            X432  = X43*cos(theta2(j+2*(i-1)))-Y43*sin(theta2(j+2*(i-1)));
            Y432 = cos(alpha12)*(X43*sin(theta2(j+2*(i-1)))+Y43*cos(theta2(j+2*(i-1))))-sin(alpha12)*Z43;
            Z432 = sin(alpha12)*(X43*sin(theta2(j+2*(i-1)))+Y43*cos(theta2(j+2*(i-1))))+cos(alpha12)*Z43;
            %2 Theta 1
            X4321=X432*cos(theta1(i))-Y432*sin(theta1(i));
            Y4321=cos(alpha71)*(X432*sin(theta1(i))+Y432*cos(theta1(i)))-sin(alpha71)*Z432;
            Z4321=sin(alpha71)*(X432*sin(theta1(i))+Y432*cos(theta1(i)))+cos(alpha71)*Z432;
            %1 theta7
            X43217 = X4321*cos(theta7)-Y4321*sin(theta7);
            Y43217 = cos(alpha67)*(X4321*sin(theta7)+Y4321*cos(theta7))-sin(alpha67)*Z4321;
            theta6s(1) = X43217;

            theta6c(1) = Y43217;
            theta6(n+2*(j-1)+4*(i-1)) = atan2(theta6s, theta6c)

         end
     end
end
end

        