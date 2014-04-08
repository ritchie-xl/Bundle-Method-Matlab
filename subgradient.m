function [ output ] = subgradient(w)
%SUBGRADIENT Summary of this function goes here
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

sub = [];
for i = 1:length(x)
    if y(i) * x(i,:)*w' > 1
        sub = [sub zeros(n-1,1)];
    else
        sub = [sub (-y(i)*x(i,:))'];
    end
end
    output = mean(sub,2);
end

