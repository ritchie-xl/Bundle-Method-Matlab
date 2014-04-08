function [optval, Ys] = bundle_method(x0, f, g, gamma,m, delta, epislon, omega, n)
%f: the loss function
%g: subgradient function
%gamma: proximity contrl parameter 
%m: descent coefficient
%delta: the tolerance patemeter to determine whether exit the main
            ...iteration
%epislon: tolerence parameter
%omega: decay coeffient
%n: dimension of the data

%Initialization
X =[];
if isempty(x0) == 1
    x = zeros(1,n);
else
    x = x0;
end
X = [X;x];

%start from feasible point (v,d) = (0,0)
y = x;
%iteration counter
t = 1; 

%Compute Linearization Error
I_plus = [];
I_minus = [];
J_plus = zeros(1,n);
g_plus = [];
g_minus = zeros(n,1);
a_plus = [];
a_minus = [];
[num dem] = size(X);
for i =1:num
    a(i) = f(y) - f(X(i,:)) + g(X(i,:))'*(y -X(i,:))';
    if a(i) >=0
        I_plus = [I_plus;X(i,:)];
        g_plus = [g_plus g(X(i,:))];
        a_plus = [a_plus;a(i)];
        if a(i) <= epislon
            J_plus = [J_plus;g(X(i,:))'];
        end
    else
        I_minus = [I_minus;X(i,:)];
        g_minus = [g_minus g(X(i,:))];
        a_minus = [a_minux;a(i)];
    end
end

%main iteration
while sqrt(g(y)'*g(y)/n) / f(y) > delta
    %solve the quadratic problem
    % v = w2v*w
    % d = w2d*w
    n_hat = n + 1;
    w2v = [1 zeros(1,n)];
    w2d = [zeros(n,1) diag(ones(n,1),0)];
    cvx_begin quiet
        cvx_solver mosek
        variable w(n_hat)
        minimize (0.5 * w' * (w2d)' * w2d * w + gamma * w2v * w )
        subject to
            (w2v) * w >= g_plus'*(w2d)*w - a_plus
            (w2v) * w <= g_minus'*(w2d)*w - a_minus
            (w2v) * w <= 0
    cvx_end
    
    v = (w2v)*w;
    d = (w2d)*w;
    
    %Determing new point
    x_tmp = y + d';
fprintf('t = %d, f(y) = %.4f, f(tmp) = %.4f\n', t, f(y), f(x_tmp));    
    
    %check if update the stability point
    if f(x_tmp) <= f(y) + m * v;
        y = x_tmp;
        X = [X;x_tmp];
        I_plus = [];
        %compute the Linearization error for the new stability point
        I_minus = [];
        J_plus = zeros(1,n);
        g_plus = [];
        g_minus = zeros(n,1);
        a_plus = [];
        a_minus = [];
        [num dem] = size(X);
        for i =1:num
            a(i) = f(y) - f(X(i,:)) + g(X(i,:))'*(y -X(i,:))';
            if a(i) >=0
                I_plus = [I_plus;X(i,:)];
                g_plus = [g_plus g(X(i,:))];
                a_plus = [a_plus;a(i)];
                if a(i) <= epislon
                    J_plus = [J_plus;g(X(i,:))'];
                end
            else
                I_minus = [I_minus;X(i,:)];
                g_minus = [g_minus g(X(i,:))];
                a_minus = [a_minus;a(i)];
            end
        end
        
        %check the first exit criterion
        if sqrt(g(y)'*g(y)/n) / f(y) <= delta
            break;
        end
        
        %check the second exit criterion
        cvx_begin quiet
            cvx_solver mosek
            variable lambda(n)
            minimize (lambda' * (J_plus)' * J_plus * lambda)
            subject to
                ones(1,n) * lambda == 1
                lambda >= 0
        cvx_end
        g_star = J_plus * lambda;
        if sqrt(g_star'*g_star/n)/f(y) <= delta
            break;
        end
    else
        %Decrease the proximity control
        gamma = omega * gamma;
    end
    t = t + 1;
end

% while end, get the optimal value
optval = y;

Ys = [];
for i = 1:size(X,1)
    Ys = [Ys; f(X(i,:))];
end

end