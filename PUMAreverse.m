function [ solutionmatrix ] = PUMAreverse( s6,A67_F,S6_F,Ptool_F,Ptool_6 )
%lower case are scalars, upper case are vectors

% constant mechanism parameters:
a12=0; a23=17; a34=0.8; a45=0; a56=0; 
s2=5.9; s3=0; s4=17; s5=0;
alpha12=pi/2; alpha23=0; alpha34=3*pi/2; alpha45=pi/2; alpha56=pi/2; 


% CLOSE THE LOOP
% outputs: a71,s7,s1,alpha71,theta7,gamma1

i=[1;0;0]; j=[0;1;0]; k=[0;0;1];
P6o_F=Ptool_F-dot(Ptool_6,i)*A67_F-dot(Ptool_6,j)*cross(S6_F,A67_F)-dot(Ptool_6,k)*S6_F;

S1_F=[0;0;1];
S7_F=cross(A67_F,S6_F);
A71_F=cross(S7_F,S1_F)/norm(cross(S7_F,S1_F));

cos71=dot(S7_F,S1_F);
sin71=dot(cross(S7_F,S1_F),A71_F);
alpha71=atan2(sin71,cos71)    
%check if S1 and S7 are parallel
if abs(cos71)>0.999 
    s7=0;
    s1=-dot(P6o_F,S1_F);
    a71=norm(-(P6o_F+s1*S1_F));
    
    %check if S1 and S7 are collinear
    if abs(a71)<0.001
        theta7=0;
        
        cosgamma1=dot(A71_F,[1;0;0]);
        singamma1=dot(cross(A71_F,[1;0;0]),S1_F);
        gamma1=atan2(singamma1,cosgamma1);
        
    else
        A71_F=-(P6o_F+s1*S1_F)/a71;
        
        cos7=dot(A67_F,A71_F);
        sin7=dot(cross(A67_F,A71_F),S7_F);
        theta7=atan2(sin7,cos7);
    
        cosgamma1=dot(A71_F,[1;0;0]);
        singamma1=dot(cross(A71_F,[1;0;0]),S1_F);
        gamma1=atan2(singamma1,cosgamma1);
    
    end

else 

    cos7=dot(A67_F,A71_F);
    sin7=dot(cross(A67_F,A71_F),S7_F);
    theta7=atan2(sin7,cos7);
    
    cosgamma1=dot(A71_F,[1;0;0]);
    singamma1=dot(cross(A71_F,[1;0;0]),S1_F);
    gamma1=atan2(singamma1,cosgamma1);
    
    s7=dot(cross(S1_F,P6o_F),A71_F)/sin71;
    a71=dot(cross(P6o_F,S1_F),S7_F)/sin71
    s1=dot(cross(P6o_F,S7_F),A71_F)/sin71;

end

% DEFINE a67 AND alpha67

a67=0;
alpha67=pi/2;


% SOLVE FOR theta1 AND phi1

X7=sin(alpha67)*sin(theta7)
Y7=-(sin(alpha71)*cos(alpha67)+cos(alpha71)*sin(alpha67)*cos(theta7));
Z7=cos(alpha71)*cos(alpha67)-sin(alpha71)*sin(alpha67)*cos(theta7);
Ath1=s6*Y7-s7*sin(alpha71);
Bth1=s6*X7+a71;
Dth1=s2;
singamma=Bth1/sqrt((Ath1^2)+(Bth1^2));
cosgamma=Ath1/sqrt((Ath1^2)+(Bth1^2));
gamma_1=atan2(singamma,cosgamma);

theta1_a=gamma_1+acos(-Dth1/sqrt((Ath1^2)+(Bth1^2)))
theta1_b=2*pi-theta1_a+2*gamma_1
theta1=[theta1_a theta1_a theta1_a theta1_a theta1_b theta1_b theta1_b theta1_b];

phi1_a=theta1_a-gamma1;
phi1_b=theta1_b-gamma1;
phi1=[phi1_a phi1_a phi1_a phi1_a phi1_b phi1_b phi1_b phi1_b];


% SOLVE FOR theta3

theta3=ones(1,8);
X1=ones(1,8);
Y1=ones(1,8);
Z1=ones(1,8);
X71=ones(1,8);
Y71=ones(1,8);
Z71=ones(1,8);
A3=ones(1,8);
B3=ones(1,8);
gamma3=ones(1,8);
for i=1:8
    X1(i)=sin(alpha71)*sin(theta1(i))
    Y1(i)=-(sin(alpha12)*cos(alpha71)+cos(alpha12)*sin(alpha71)*cos(theta1(i)))
    Z1(i)=cos(alpha12)*cos(alpha71)-sin(alpha12)*sin(alpha71)*cos(theta1(i))
    X71(i)=X7*cos(theta1(i))-Y7*sin(theta1(i))
    Y71(i)=cos(alpha12)*(X7*sin(theta1(i))+Y7*cos(theta1(i)))-Z7*sin(alpha12)
    Z71(i)=sin(alpha12)*(X7*sin(theta1(i))+Y7*cos(theta1(i)))+Z7*cos(alpha12)

    A3(i)=-s6*X71(i)-s7*X1(i)-a71*cos(theta1(i))
    B3(i)=s1-s6*Y71(i)-s7*Y1(i)
    Ath3=2*a23*a34
    Bth3=-2*a23*s4
    Dth3=(a23^2)+(a34^2)+(s4^2)-(A3(i)^2)-(B3(i)^2)
    singamma=Bth3/sqrt((Ath3^2)+(Bth3^2));
    cosgamma=Ath3/sqrt((Ath3^2)+(Bth3^2));
    gamma3(i)=atan2(singamma,cosgamma);
    theta3(i)=gamma3(i)+acos(-Dth3/sqrt((Ath3^2)+(Bth3^2)));
end
for i=3:4
    theta3(i)=2*pi-theta3(i-2)+2*gamma3(i); 
end
for i=7:8
    theta3(i)=2*pi-theta3(i-2)+2*gamma3(i); 
end


% SOLVE FOR theta2

theta2=ones(1,8);
for i=1:8
    A1_2=a23+a34*cos(theta3(i))-s4*sin(theta3(i));
    B1_2=-a34*sin(theta3(i))-s4*cos(theta3(i));
    D1_2=A3(i);
    A2_2=-a34*sin(theta3(i))-s4*cos(theta3(i));
    B2_2=-a23-a34*cos(theta3(i))+s4*sin(theta3(i));
    D2_2=B3(i);

    coeff=[A1_2,B1_2;A2_2,B2_2];
    D=[D1_2;D2_2];
    X=(coeff^-1)*D;
    x=X(1);
    y=X(2);
    theta2(i)=atan2(y,x);
end


% SOLVE FOR theta5

theta5=ones(1,8);
X712=ones(1,8);
Y712=ones(1,8);
Z712=ones(1,8);
X7123=ones(1,8);
Y7123=ones(1,8);
Z7123=ones(1,8);
for i=1:8
    X712(i)=X71(i)*cos(theta2(i))-Y71(i)*sin(theta2(i));
    Y712(i)=cos(alpha23)*(X71(i)*sin(theta2(i))+Y71(i)*cos(theta2(i)))-sin(alpha23)*Z71(i);
    Z712(i)=sin(alpha23)*(X71(i)*sin(theta2(i))+Y71(i)*cos(theta2(i)))+cos(alpha23)*Z71(i);
    X7123(i)=X712(i)*cos(theta3(i))-Y712(i)*sin(theta3(i));
    Y7123(i)=cos(alpha34)*(X712(i)*sin(theta3(i))+Y712(i)*cos(theta3(i)))-sin(alpha34)*Z712(i);
    Z7123(i)=sin(alpha34)*(X712(i)*sin(theta3(i))+Y712(i)*cos(theta3(i)))+cos(alpha34)*Z712(i);
    theta5(i)=acos(-Z7123(i));
end

for i=1:4
    theta5(2*i)=2*pi-theta5(2*i-1);
end


% SOLVE FOR theta4

theta4=ones(1,8);
costheta4=ones(1,8);
sintheta4=ones(1,8);
for i=1:8
    costheta4(i)=X7123(i)/sin(theta5(i));
    sintheta4(i)=-Y7123(i)/sin(theta5(i));
    theta4(i)=atan2(sintheta4(i),costheta4(i));
end


% SOLVE FOR theta6

theta6=ones(1,8);
X_4=ones(1,8);
Y_4=ones(1,8);
Z_4=ones(1,8);
X43=ones(1,8);
Y43=ones(1,8);
Z43=ones(1,8);
X432=ones(1,8);
Y432=ones(1,8);
Z432=ones(1,8);
X4321=ones(1,8);
Y4321=ones(1,8);
Z4321=ones(1,8);
X43217=ones(1,8);
Y43217=ones(1,8);
sintheta6=ones(1,8);
costheta6=ones(1,8);
for i=1:8
    X_4(i)=sin(alpha45)*sin(theta4(i))
    Y_4(i)= -(sin(alpha34)*cos(alpha45)+cos(alpha34)*sin(alpha45)*cos(theta4(i)))
    Z_4(i)= cos(alpha34)*cos(alpha45)-sin(alpha34)*sin(alpha45)*cos(theta4(i))
    
    X43(i)=X_4(i)*cos(theta3(i))-Y_4(i)*sin(theta3(i))
    Y43(i)=cos(alpha23)*(X_4(i)*sin(theta3(i))+Y_4(i)*cos(theta3(i)))-sin(alpha23)*Z_4(i)
    Z43(i)=sin(alpha23)*(X_4(i)*sin(theta3(i))+Y_4(i)*cos(theta3(i)))+cos(alpha23)*Z_4(i)
    
    X432(i)=X43(i)*cos(theta2(i))-Y43(i)*sin(theta2(i))
    Y432(i)=cos(alpha12)*(X43(i)*sin(theta2(i))+Y43(i)*cos(theta2(i)))-sin(alpha12)*Z43(i)
    Z432(i)=sin(alpha12)*(X43(i)*sin(theta2(i))+Y43(i)*cos(theta2(i)))+cos(alpha12)*Z43(i)
    
    X4321(i)=X432(i)*cos(theta1(i))-Y432(i)*sin(theta1(i))
    Y4321(i)=cos(alpha71)*(X432(i)*sin(theta1(i))+Y432(i)*cos(theta1(i)))-sin(alpha71)*Z432(i)
    Z4321(i)=sin(alpha71)*(X432(i)*sin(theta1(i))+Y432(i)*cos(theta1(i)))+cos(alpha71)*Z432(i)
    X43217(i)=X4321(i)*cos(theta7)-Y4321(i)*sin(theta7)
    Y43217(i)=cos(alpha67)*(X4321(i)*sin(theta7)+Y4321(i)*cos(theta7))-sin(alpha67)*Z4321(i)
    sintheta6(i)=X43217(i);
    costheta6(i)=Y43217(i);
    theta6(i)=atan2(sintheta6(i),costheta6(i));
end


% OUTPUT SOLUTION
% each column is a solution set for phi1 through theta6

solutionmatrix=ones(6,8);
solutionmatrix(1,:)=phi1;
solutionmatrix(2,:)=theta2;
solutionmatrix(3,:)=theta3;
solutionmatrix(4,:)=theta4;
solutionmatrix(5,:)=theta5;
solutionmatrix(6,:)=theta6;

solutionmatrix=180*solutionmatrix./pi;

for i=1:8
    for j=1:6
        if solutionmatrix(j,i) > 180
            solutionmatrix(j,i)=solutionmatrix(j,i)-360;
        end
    end
end        