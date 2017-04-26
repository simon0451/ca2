%% Simon Popecki
% Computing Assignment 2
% ME603
clear all;
close all;
%% Loading Parameters
%when modifying parameters, do not delete the white space.
%Replace the existing number(s) with a new number(s), make no other changes.
parameters = textread('parameters.txt','%s'); %reading variables in from a text file
qdp = str2double(parameters(4)); %importing the data as a double
h = str2double(parameters(6));
Tinf = str2double(parameters(8));
k = str2double(parameters(10));
L = str2double(parameters(12));
D = str2double(parameters(14));

%% Non-dimensionalize Variables

numX = 10; %the number of nodes in the x direction
numY = numX; %since the matrix needs to be square, the Y and X directions have the same number of nodes, the dimensions of the nodes are different however, if the part is not square 
numdx = numX/L; %(nodes/m) = number of nodes per meter
numdy = D*numY; %(nodes/m) = number of nodes per meter
dx = 1/numdx; %meters per node
dy = 1/numdy; %meters per node


%% Generating The Nodal Type Matrix
%this matrix tells the switch statement which case to execute for which
%parts of the solid - i.e. a corner must be treated differently from an
%internal node
nodetype = zeros(numX,numY); %initializing the matrix
%nodetype(row,column)
nodetype(1,1) = 1; %the top left corner is case 1
nodetype(1,2:(end-1)) = 2; %the top edge is case 2
nodetype(1,end) = 3; %the top right corner is case 3
nodetype(2:(end-1),1) = 4; %the left edge is case 4
nodetype(2:(end-1),2:(end-1)) = 5; %the central nodes are case 5
nodetype(2:(end-1),end) = 6; %the right edge is case 6
nodetype(end,1) = 7; %the lower left corner is case 7
nodetype(end,(2:(end-1))) = 8; %the bottom edge is case 8
nodetype(end,end) = 9; %the lower right corner is case 9
%%
for i = numY %this is the number of rows
    for j = numX %col is the number of columns
        n = j;
        
        switch nodetype(i,j)
            
            case 1 %Top left corner (left blank)
                A(n,n) = -(k*dy/(2*dx)+k*dx/(2*dy)+h(dx/2+dy/2));
                A(n,n+1) = k*dy/(2*dx); %half area conduction from the right side
                A(n,numX) = k*dx/(2*dy); %coming up form the bottom?
                B(n) = -h*(dx/2)*Tinf; %convection, but half value since it is at the corner
            case 2 %Top edge exposed to convection
                A(n,n) = -(k*dx/dx+k+dy/dx-h*dx*Tinf);
                A(n,n-1) = k*dx/(2*dx); %half area conduction from the left side
                A(n,n+1) = k*dx/(2*dx); %half area conduction from the right side
                A(n,n+numX) = k*dy/dx;
                B(n) = -h*dx*Tinf; %there is always convection at the top of the array
                
            case 3 %Top right corner, exposed to convection on two sides
                B(n) = -h*((dx+dy)/2)*Tinf; %coming in from both y and x directions, only getting half area (div/2)
                
            case 4 %Left side, insulated on one side, conduction on the other side
                B(n) = 0; %insulated, no convection
            case 5 %Central node - conduction on all sides
                B(n) = 0; %there is no convection in the central nodes
                
            case 6 %Right side - conduction on one side and convection on the other
                B(n) = -h*dy*Tinf; %coming from y direction
            case 7 %Lower left corner, exposed to a heat flux in on the bottom and insulation on the right
                B(n) = -h*(dx/2)*Tinf; %coming perpendicular to x, half since it's a corner
            case 8 %Bottom side, exposed to heat flux on one side, and conduction on the other
                B(n) = -h*dx*Tinf;
            case 9 %Bottom right corner exposed to to convection on side
                A(n,n) = -(k*dx/(2*dx)+k*dx/(2*dy)+h*(dy/2));
                A(n,n-1) = k*dy/(2*dx);
                A(n,n+1) = k*dx/(2*dy);
                A(n,n-numX) = k*dx/dy;
                B(n) = -h*(dy/2)*Tinf;
        end
    end
end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                