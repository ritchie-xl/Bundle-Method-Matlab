
data = load('pima-indians-diabetes.data');
[m n] =size(data);
label = data(:,n);
data = data(:,1:n-1);
x0 = rand(1,n-1);

%[optval] = bundle_method(f, g, gamma,m, delta, epislon, w, n)w
[w, Ys] = bundle_method(x0, @loss,@subgradient, 1, 0.1,  0.01, 0.1, 0.9,n-1); 
