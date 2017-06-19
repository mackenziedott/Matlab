function [theta7, y1, alpha71, S7, S1, a71] = closeTheLoop(tool_6,tool_F, tool_S6_F, tool_a67_F)

%Find coordinates of orgin point of 6th coord system
Origin_F_6 = tool_F - dot(tool_6,[1; 0; 0])*tool_a67_F - dot(tool_6,[0;1;0])*cross(tool_S6_F, tool_a67_F) - dot(tool_6, [0;0;1])*tool_S6_F;
%close-the-loop
S7_F = cross(tool_a67_F, tool_S6_F);
S1_F = [0;0;1];

a71_F = (cross(S7_F, S1_F))/(norm(cross(S7_F,S1_F)));

alpha71c(1) = dot(S7_F, S1_F);
alpha71s(1) = dot(cross(S7_F,S1_F),a71_F);

alpha71 = atan2(alpha71s,alpha71c);

theta7c(1) = dot(tool_a67_F,a71_F);
theta7s(1) = dot(cross(tool_a67_F,a71_F),S7_F);
theta7 = atan2(theta7s,theta7c);
y1c(1) = dot(a71_F, [1;0;0]);

y1s(1) = dot(cross(a71_F,[1;0;0]),S1_F);
y1 = atan2(y1s,y1c);
S7 = dot(cross(S1_F,Origin_F_6),a71_F)/sin(alpha71);
a71 = dot(cross(Origin_F_6,S1_F),S7_F)/sin(alpha71);
S1 = dot(cross(Origin_F_6,S7_F),a71_F)/sin(alpha71);