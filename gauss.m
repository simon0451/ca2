function [T] = gauss(a,b,row,col)
% global a b col row

%initializing array
lengtha = length(a(:,1));
T = ones(lengtha(:,1),1);

%initializing the error such that the while loop will operate
%done via newt
newt = 99*ones((length(a(:,1))),1); %random error larger than increment
absdif = abs((newt-T)./T);
err = max(max(absdif));
nnode = row*col;

iterations = 1;

while (err > .000001)
    newt = T;
    for i = 1:1:nnode
        T(i) = 1/a(i,i)*b(i)-a(i,1:i-1)*(T(1:i-1)-a(i,i+1:nnode)*T(i+1:nnode));
    end
    err = max(max(abs(newt - T)./T));
    iterations = iterations+1;
    if (iterations > 1000000) %runaway catcher
        break
    end
end

end