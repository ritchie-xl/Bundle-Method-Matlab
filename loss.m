function [ output ] = loss(w)
%f = @(x) x'*x/2 + max(0, 1 - label * dot(data,x));
%LOSS Summary of this function goes here
%   Detailed explanation goes here
data = load('pima-indians-diabetes.data');
[m n] =size(data);
y = data(:,n);
x = data(:,1:n-1);
for i = 1:length(y)
    if y(i) == 0
        y(i) = -1;
    end
end
output = 0;
l = 1 - y.*(x*w');
for i = 1:length(l)
    if l(i) >=0
        output = output + l(i);
    end
end
        

end

