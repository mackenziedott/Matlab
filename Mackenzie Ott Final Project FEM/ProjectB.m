clear all
clc
% Final Project FEM
% Project B
% FEA Tool - 2D Heat Transfer FEA
% Mackenzie Ott
%Map of Nodes
% 1-3
% | |
% 2-4

%set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%pull from file
num = xlsread('FEMINPUT');
%determine number of Elements
numberOfelements = size(num); 

%pull values
for i = 1:numberOfelements(1)
    
    %Convection Flags
    ConvectFlag12(i) = num(i,13);
    ConvectFlag23(i) = num(i,14);
    ConvectFlag34(i) = num(i,15);
    ConvectFlag41(i) = num(i,16);
    h(i) = num(i,12);
    %Heat Flux Flag and constants
    FluxFlag12(i) = num(i,17);
    FluxFlag23(i) = num(i,18);
    FluxFlag34(i) = num(i,19);
    FluxFlag41(i) = num(i,20);
    q0(i) = num(i,21); 
    %heat generation
    f(i) = num(i,22);
    Ta(i) = num(i,23);
    %RealCoordinates
    RealCoordN1(i,1) = num(i,1);
    RealCoordN1(i,2) = num(i,2);
    RealCoordN2(i,1) = num(i,3);
    RealCoordN2(i,2) = num(i,4);
    RealCoordN3(i,1) = num(i,5);
    RealCoordN3(i,2) = num(i,6);
    RealCoordN4(i,1) = num(i,7);
    RealCoordN4(i,2) = num(i,8);
    %conductivity Constant
    k(i,1) = num(i,9);
    k(i,2) = num(i,10);
end
Le = RealCoordN1(2,1) - RealCoordN1(1,1);
%%% pull k values
%%%%%%%%%%%%%%%%%
%Create Jacobian
%Map Function
gK = zeros(numberOfelements(1)*2 +2);
%Create Global Stiffness Matrix
%Inital Conditions
s0 = 0;
t0 = 0 ;
%Create Global Stiffness Matrix
for i = 1:numberOfelements(1)

%x(i) = Map(RealCoordN1(i,1), RealCoordN2(i,1), RealCoordN3(i,1), RealCoordN4(i,1), s0, t0);
%y(i) = Map(RealCoordN1(i,2), RealCoordN2(i,2), RealCoordN3(i,2), RealCoordN4(i,2), s0, t0);
J = Jacobian(RealCoordN1,RealCoordN2,RealCoordN3,RealCoordN4,s0,t0);
invJ = inv(J);
detJ = det(J);
%Create Element Stiffness Matrix
Ke = ElementStiffness(invJ,detJ, k);
    %set Convection in Stiffness matrix
    if ConvectFlag12(i) == 1
        Ke(1,1) = Ke(1,1) +h(i);
        Ke(3,1) = Ke(3,1) +h(i)/2;
        Ke(1,3) = Ke(1,3) +h(i)/2;
        Ke(3,3) = Ke(3,3) +h(i);
    end
    if ConvectFlag23(i) == 1
       
        Ke(1,1) = Ke(1,1) +h(i);
        Ke(2,1) = Ke(2,1) +h(i)/2;
        Ke(1,2) = Ke(1,2) +h(i)/2;
        Ke(2,2) = Ke(2,2) +h(i);
    end
    if ConvectFlag34(i) == 1
        Ke(2,2) = Ke(2,2) + h(i);
        Ke(4,2) = Ke(4,2) + h(i)/2;
        Ke(2,4) = Ke(2,4) + h(i)/2;
        Ke(4,4) = Ke(4,4) + h(i);
    end
    if ConvectFlag41(i) == 1
        Ke(3,3) = Ke(3,3) + h(i);
        Ke(3,4) = Ke(3,4) + h(i)/2;
        Ke(4,3) = Ke(4,3) + h(i)/2;
        Ke(4,4) = Ke(4,4) + h(i);
    end
    Ke
sizegK=size(gK(:,1));
sizeKe = size(Ke(:,1));
%add global stiffness matrix
    for o = 1:sizegK(1)
     %creating Q side
        Q(o) = 0;
        if o <= sizeKe(1)
            for p = 1:sizegK(1)
                if p <= sizeKe(1)
                    gK(o+2*(i-1),p+2*(i-1)) = gK(o+2*(i-1),p+2*(i-1))+Ke(o,p);
                end
            end
        end
    end
end
gK;
%Compute Boundary Conditions
%Create Q
i = 0;
o = 0;
for i = 1:numberOfelements(1)
    qe = [0;0;0;0];
    %convection
    if ConvectFlag12(i) == 1
        qe(1) = h(i)*Ta(i)*Le/2;
        qe(2) = h(i)*Ta(i)*Le/2;
    end
    if ConvectFlag23(i) == 1
        qe(2) = h(i)*Ta(i)*Le/2;
        qe(3) = h(i)*Ta(i)*Le/2;
    end
    if ConvectFlag34(i) == 1
        qe(3) = h(i)*Ta(i)*Le/2;
        qe(4) = h(i)*Ta(i)*Le/2;
    end
    if ConvectFlag41(i) == 1
        qe(1) = h(i)*Ta(i)*Le/2;
        qe(4) = h(i)*Ta(i)*Le/2;
    end
    %heat Generation
    if f(i) ~= 0
        for w = 1:4
            qe(w) = qe(w) + f(i)*detJ;
        end
    end
    %Heatflux
    if FluxFlag12(i) == 1
        qe(1) = qe(1) + q0(i)*Le/2;
        qe(2) = qe(2) + q0(i)*Le/2;
    end
    if FluxFlag23(i) == 1
        qe(2) = qe(2) + q0(i)*Le/2;
        qe(3) = qe(3) + q0(i)*Le/2;
    end
    if FluxFlag34(i) == 1
        qe(3) = qe(3) + q0(i)*Le/2;
        qe(4) = qe(4) + q0(i)*Le/2;
    end
    if FluxFlag41(i) == 1
        qe(1) = qe(1) + q0(i)*Le/2;
        qe(4) = qe(4) + q0(i)*Le/2;
    end
    %Rearrange qe to match nodes
    qeN(1)=qe(2);
    qeN(2)=qe(3);
    qeN(3)=qe(1);
    qeN(4)=qe(4);
    qeN;
    %Add to global Matrix
    %Node 1 on left matrix is equal to Node 2 on Right matrix
        
    for z = 1:4
        Q(z+2*(i-1)) = Q(z+2*(i-1)) + qeN(z);
    end
end
%Heat Generation
%heatflux

%Compute Temperature Nodes
Q
gK
T = Q/gK
n = 200; %number of divisions
for i= 1:numberOfelements(1)
    %Map Temperature
    Tmap = [T(1+(2*(i-1))); T(2+(2*(i-1)));T(3+(2*(i-1)));T(4+(2*(i-1)))];
    %Create Tempereature Mesh
    %Create X and Y Coordinates
    for smatrix = 1:(n+1)
        s= (smatrix-(n/2)-1)/(n/2);
        for tmatrix = 1:(n+1)
            t = (tmatrix-(n/2)-1)/(n/2);
            TElement(tmatrix + (smatrix-1)*(n+1)) = TemperatureLocation(Tmap,s,t);
            TElementMap(smatrix,tmatrix) = TElement(tmatrix + (smatrix-1)*(n+1));
            xMapElement(tmatrix+ (smatrix-1)*(n+1)) = Map(RealCoordN1(i,1), RealCoordN2(i,1), RealCoordN3(i,1), RealCoordN4(i,1), s, t);
            yMapElement(tmatrix+ (smatrix-1)*(n+1))= Map(RealCoordN1(i,2), RealCoordN2(i,2), RealCoordN3(i,2), RealCoordN4(i,2), s, t);
            %produce heatflux for an element
            qf(:,i) = heatflux(Tmap,s,t,invJ,k);
           
        end
    end
    %Add to Global Meshes
    for y = 1:(n+1)^2
        TTotal(y + (i-1)*(n+1)^2) = TElement(y);
        XTotal(y + (i-1)*(n+1)^2) = xMapElement(y);
        YTotal(y + (i-1)*(n+1)^2) = yMapElement(y);
        ZTotal(y+ (i-1)*(n+1)^2) = 0;
    end
%     for y1 = 1:(n+1)
%         for y2 = (n+1)
%             TTotalMap(y1 + (n+1)*(i-1),y2+(n+1)*(i-1)) = TElementMap(y1,y2);
%         end
%     end
end
%Display heatflux
heatf = qf
% Data is temperature at a given x,y Coordinate
Data(:,1) = XTotal.';
Data(:,2) = YTotal.';
Data(:,3) = TTotal.';
%[Xg, Yg] = meshgrid(XTotal, YTotal);
%figure(3)
%surf(Xg,Yg, TTotal)
%figure(1)
%plot3(XTotal, YTotal, TTotal)

%xlabel('This is an x label','fontsize',14,'fontweight','bold','color',[1 0 0]) 
%ylabel('This is a y label','fontsize',14,'fontweight','bold','color',[0 0 0]) 
%zlabel('This is a z label','fontsize',14,'fontweight','bold','color',[0 0 1]) 
%Using Color Line Function, produces 2-D color Map
figure(2)
color_line3(XTotal, YTotal, ZTotal, TTotal,'o')
xlabel('x','fontsize',14,'fontweight','bold','color',[1 0 0]) 
ylabel('y','fontsize',14,'fontweight','bold','color',[0 0 0])
c = colorbar;
c.Label.String = 'Temperature (Farenheit)';
%Saves as a Figure - Redundency
savefig('TempMapOutput')
%Saves as a JPEG
hgsave('TempMapOutput.jpg')
%write Data to .txt file
header1 = 'X';
header2 = 'Y';
header3 = 'Temperature';
header4 = 'qx ';
header5 = 'qy ';
%header = [header1 header2 header3];
fid=fopen('FEMOutput.txt','w');
fprintf(fid,'%6s %6s %6s\r\n',header1,header2, header3);
fprintf(fid,'%6.3f %6.3f %6.3f\r\n', Data.');
fprintf(fid,'%6s %6s\r\n', header4, header5);
fprintf(fid,'%6.3f %6.3f\r\n', heatf(:,1), heatf(:,2));
fclose(fid);
