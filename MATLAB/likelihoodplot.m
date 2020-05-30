% Set up parameter space
mu = 2;
alpha = 2:0.1:7;
beta = 2:0.1:7;
P = 1;

% Create a simulation
Y = [3];
dist = ["Constant"];
delta = [5];
N = [0];
params = {};
params{1} = [3];

[t,n,r,intensity] = Hawkessim(mu, Y, dist, delta, N, params, 10);

% create vector of likelihood values
tmp = [];
Q = length(alpha);
R = length(beta);
for i = 1:Q
    for j = 1:R
        tmp(i,j) =  Hawkesll(t, mu, alpha(i), beta(j), P);
    end
end

% beta values move growing along the right
% alpha values grow moving down the column

alin = linspace(min(alpha), max(alpha), 33);
blin = linspace(min(beta), max(beta), 33);
[A,B] = meshgrid(alin, blin);

Z = griddata(alpha, beta, tmp, A, B);
mesh(A,B,Z);
shading interp