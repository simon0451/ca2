%% Simon Popecki
% Computing Assignment 2
% ME603
clear all;
close all;
%% Loading Parameters
%when modifying parameters, do not delete the white space.
%Replace the existing number(s) with a new number(s), make no other changes.
parameters = textread('parameters.txt','%s'); %reading variables in from a text file
global qdp h Tinf k L D dx dy row col nodetype dxnd dynd a b
qdp = str2double(parameters(4)); %importing the data as a double
h = str2double(parameters(6));
Tinf = str2double(parameters(8));
k = str2double(parameters(10));
L = str2double(parameters(12));
D = str2double(parameters(14));

%% MESH GENERATION AND NODE LENGTH

row = 10; %the number of nodes in the x direction %%%%%%%%%%%%%ROWS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = 15; %the number of nodes in the y direction %%%%%%%%%%%%%COLUMNS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = L/(col-1); 
dy = D/(row-1); 
dxnd = dx/L; 
dynd = dy/D; 

xlab = linspace(0,L,col);
ylab = linspace(0,D,row);
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

%% Creating a non-dimensionalized matrix set

%This function computes the solution using nondimensional values, then
%outputs both the non-dimensional solution AND a dimensioned solution
[a,b] = nondim(); %produces a nondimensionalized solution (a,b)

%% Calculating the Matrix Values, Non-Dimensionalized
bp = b';
%C = gauss2(a,b,row,col);


%%
C = a\b; %using this one for faster computations - see gauss2.m
c = vec2mat(C,col); %turning the results of the matrix math into a matrix

%% Dimensioning the Values
cc = dimensionalize(c);

%% Find the position and value of the maximum temperature
[temp,location] = max(cc(:)); %temperature is in Kelvin
[rowmax,colmax] = ind2sub(size(cc),location); %row and column of the maximum temperature node
%column is how far along the L axis it is
disp('Maximum Temperature (K):')
disp(temp)
disp('Maximum Temperature, Position in Matrix [row; column]')
disp(rowmax)
disp(colmax)
disp('Dimensioned Location of Max. Temp. (meters), origin is top left corner (see plot, left-hand rule)')
disp('X-axis:')
xaxis = (colmax-1)*dx;
disp(xaxis)
disp('Y-axis:')
yaxis = (rowmax-1)*dy;
disp(yaxis)

%% Plotting Resulting Values
imagesc(xlab,ylab,cc)
title('Heat Distribution Within Material (left hand coordinates)')
xlabel('X-Position (m)')
ylabel('Y-Position (m)')
colo = colorbar;
xlabel(colo,'Temperature (K)')

%% Heat Flux

for i = 1:1:length(cc(1,:))-1
    tb(i) = abs((cc(1,i+1))-cc(1,i));
end
tmind = find(tb==max(tb)); %this is more efficient than the previous method of indexing
qmaxt = xlab(tmind);
for i = 1:1:(length(cc(:,end))-1)
    rb(i) = abs(cc(i+1,end)-cc(i,end));
end
rmind = find(rb==max(rb));
qmaxr = ylab(rmind);
disp('Position of Maximum Heat Flux (m)')
disp('Top:')
disp(qmaxt)
disp('Right:')
disp(qmaxr)







