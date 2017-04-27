function [a,b] = nondim()
%this function produces a solution with non-dimensionalized variables
global qdp h Tinf k L D row col nodetype dxnd dynd
n = 1; %node counter variable
A = zeros(row*col,row*col);
dy = dynd*D;
dx = dxnd*L;

for i = 1:1:row %this is iterating through a row
    for j = 1:1:col
        
        switch nodetype(i,j)
            
            case 1 %Top left corner (left blank)
                A(n,n) = -((k*.5*dy/(dx)+(k*(dx/2)/(dy)))+h*(dx/2));
                A(n,n+1) = k*.5*dy/(dx); %node to the right
                %insulated to the left
                %no node above - convection
                A(n,n+col) = k*(dx/2)/(dy); %node below
                B(n) = -h*(dx/2)*Tinf; %convection, but half value since it is at the corner
                n = n+1; %increamenting node counter
                
            case 2 %Top edge exposed to convection
                A(n,n) = -((k*dx/dy)+2*(k*.5*dy/dx)+h*dx);
                A(n,n+1) = k*.5*dy/dx; %node to the right
                A(n,n-1) = k*.5*dy/dx; %node to the left
                %no node above - convection, this is the top edge!
                A(n,n+col) = k*dx/dy; %this is the node below
                B(n) = -h*dx*Tinf; %there is convection at the top of the solid perpendicular to the x axis
                n = n+1; %increamenting node counter
                
            case 3 %Top right corner, exposed to convection on two sides
                A(n,n) = -((k*(dy/2)/(dx))+(k*(dx/2)/(dy))+h*((dx+dy)/2));
                %node on the right is convection
                A(n,n-1) = k*dy/2/(dx); %node on the left
                A(n,n+col) = k*dx/2/(dy); %node below
                %node above is convection
                B(n) = -h*((dx+dy)/2)*Tinf; %coming in from both y and x directions, both only getting half area (div/2)
                n = n+1; %increamenting node counter
                
            case 4 %Left side, insulated on one side, conduction on the other side
                A(n,n) = -(k*dx/dy+k*dy/(dx));
                A(n,n+1) = k*dy/(dx); %node on the right
                %node on the left is insulated
                A(n,n+col) = k*dx/2/dy; %node below
                A(n,n-col) = k*dx/2/dy; %node above
                B(n) = 0; %insulated, no convection
                n = n+1; %increamenting node counter
                
            case 5 %Central node - conduction on all sides
                A(n,n) = -(2*(k*dx/dy)+2*(k*dy/dx)); %we are summing the values equal to the Tm,n
                A(n,n+1) = k*dy/dx;%one on the right
                A(n,n-1) = k*dy/dx;%one on the left
                A(n,n+col) = k*dx/dy;%one below (+numX means we are going up or down a ROW)
                A(n,n-col)= k*dx/dy;%one above
                B(n) = 0; %there is no convection in the central nodes
                n = n+1; %increamenting node counter
                
            case 6 %Right side - conduction on one side and convection on the other
                A(n,n) = -((k*dx/dy)+(k*dy/(dx))+h*dy);
                %node on the right is convection
                A(n,n-1) = k*dy/(dx); %node on the left
                A(n,n+col) = k*dx/2/dy; %node below
                A(n,n-col) = k*dx/2/dy; %node above
                B(n) = -h*dy*Tinf; %coming from y direction
                n = n+1; %increamenting node counter
                
            case 7 %Lower left corner, exposed to a heat flux in on the bottom and insulation on the right
                A(n,n) = -((k*(dy/2)/(dx))+(k*(dx/2)/(dy)));
                A(n,n+1) = k*dy/2/(dx); %node on the right
                %node on the left is insulated
                %node below is heat flux, handled by B matrix
                A(n,n-col) = k*dx/2/(dy); %node above
                B(n) = -qdp*dx/2;
                n = n+1; %increamenting node counter
                
            case 8 %Bottom side, exposed to heat flux on one side, and conduction on the other
                A(n,n) = -((k*dx/dy)+2*(k*.5*dy/dx));
                A(n,n+1) = k*.5*dy/dx; %node to the right
                A(n,n-1) = k*.5*dy/dx; %node to the left
                %no node below, heat flux
                A(n,n-col) = k*dx/dy; %this is the node above
                B(n) = -qdp*dx; %exposed to heat flux only, no convection
                n = n+1; %increamenting node counter
                
            case 9 %Bottom right corner exposed to to convection on side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                A(n,n) = -(k*(dy/2)/(dx)+k*(dx/2)/(dy)+h*(dy/2));
                %node on the right is convection
                A(n,n-1) = k*dy/2/(dx); %node on the left
                %node below is handled by the heat flux in
                A(n,n-col) = k*dx/2/(dy); %node above
                B(n) = -h*(dy/2)*Tinf-qdp*dx/2;
                n = n+1; %increamenting node counter
        end
    end
end
a = A*Tinf;
% b = B*Tinf;
% B = B'; %transposing the B  matrix to convert it from a row to a column
b = B'; %same for the non-dimensionalized matrix
end