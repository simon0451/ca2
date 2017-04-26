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

row = 10; %the number of nodes in the x direction
col = row; %since the matrix needs to be square, the Y and X directions have the same number of nodes, the dimensions of the nodes are different however, if the part is not square 
numdx = row/L; %(nodes/m) = number of nodes per meter
numdy = D*col; %(nodes/m) = number of nodes per meter
dx = 1/numdx; %meters per node
dy = 1/numdy; %meters per node


%% Generating The Nodal Type Matrix
%this matrix tells the switch statement which case to execute for which
%parts of the solid - i.e. a corner must be treated differently from an
%internal node
nodetype = zeros(row,col); %initializing the matrix
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
for i = 1:1:col %this is the number of rows
    for j = i:1:row %number of columns
        n = i; %?
        
        switch nodetype(i,j)
            
            case 1 %Top left corner (left blank)
                A(n,n) = -((k*.5*dy/(dx/2)+(k*dx/2/(dy/2)))+(h*(dx/2)*Tinf));
                A(n,n+1) = k*.5*dy/(dx/2); %node to the right
                %insulated to the left
                %no node above - convection
                A(n,n+row) = k*dx/2/(dy/2); %node below
                B(n) = -h*(dx/2)*Tinf; %convection, but half value since it is at the corner
                
            case 2 %Top edge exposed to convection
                A(n,n) = -((h*dx)+(k*dx/dy)+2*(k*.5*dy/dx));
                A(n,n+1) = k*.5*dy/dx; %node to the right
                A(n,n-1) = k*.5*dy/dx; %node to the left
                %no node above - convection, this is the top edge!
                A(n,n+row) = k*dx/dy; %this is the node below
                B(n) = -h*dx*Tinf; %there is convection at the top of the solid perpendicular to the x axis
                
            case 3 %Top right corner, exposed to convection on two sides
                
                B(n) = -h*((dx+dy)/2)*Tinf; %coming in from both y and x directions, both only getting half area (div/2)
                
            case 4 %Left side, insulated on one side, conduction on the other side
                B(n) = 0; %insulated, no convection
                
            case 5 %Central node - conduction on all sides
                A(n,n) = -(2*(k*dx/dy)+2*(k*dy/dx)); %we are summing the values equal to the Tm,n
                A(n,n+1) = k*dx/dy;%one on the right
                A(n,n-1) = k*dx/dy;%one on the left
                A(n,n+row) = k*dy/dx;%one below (+numX means we are going up or down a ROW)
                A(n,n-row)= k*dy/dx;%one above
                B(n) = 0; %there is no convection in the central nodes
                
            case 6 %Right side - conduction on one side and convection on the other
                B(n) = -h*dy*Tinf; %coming from y direction
                
            case 7 %Lower left corner, exposed to a heat flux in on the bottom and insulation on the right
                B(n) = -h*(dx/2)*Tinf; %coming perpendicular to x, half since it's a corner
                
            case 8 %Bottom side, exposed to heat flux on one side, and conduction on the other
                A(n,n) = -((k*dx/dy)+2*(k*.5*dy/dx))+qdp*dx;
                A(n,n+1) = k*.5*dy/dx; %node to the right
                A(n,n-1) = k*.5*dy/dx; %node to the left
                A(n,n+row) = qdp*dx; %this is the node below
                A(n,n-row) = k*dx/dy; %this is the node above
                B(n) = 0; %exposed to heat flux only, no convection
                
            case 9 %Bottom right corner exposed to to convection on side

                B(n) = -h*(dy/2)*Tinf;
        end
    end
end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                