function [T] = gauss2(a,b)
%references: university of nevada computer science dept
n = length(b);
X = zeros(n,1);
ereval = ones(n,1);

iteration = 0;
while max(ereval) >.000000001
    iteration = iteration+1;
    Z = X;
    for i = 1:n
        j = 1:n; %create j array
        j(i) = []; %this null array eliminates the sum of previous values
        Xtemp = X;
        Xtemp(i) = []; %same thing here
        X(i) = (b(i)-sum(a(i,j) * Xtemp) / (a(i,i)));
    end
%     Xsolution(:,iteration) = X;
    ereval = sqrt((X-Z).^2); %calculating the error
end
T = X;



end