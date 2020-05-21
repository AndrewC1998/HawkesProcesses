function [t,n,r,intensity] = Hawkessim(mu, Y, dist, delta, N, params, tmax)
%
% Hawkessim - Simulate a multivariate Hawkes process
%             (asymmetric matrices allowed)
%
% Inputs - mu is the background intensity vector
%        - Y is matrix of initial sizes of self-excited jumps
%        - dist is the distribution function that Y follows
%        - delta is a matrix of the rate of exponential decay
%        - N is the vector of the number of events attributed to process 
%          i observed at and before time 0
%        - params is for the specific parameter choices of a given dist
%        - tmax is the maturity time
%
% Outputs - intensity is the intensity matrix. Each column is intensities
%           for process m
%         - r is the event times for all point processes
%         - N is the M corresponding count processes. Each column is one
%           process
%         - t is the event times
 
 M = length(mu);
 t = [0];
 lambda = {};
 lambda{1} = Y;
 a = {};
 j = 0;
 Nfull = [N];
 ttmp = [];
 
 % Start the while loop till maturity time
 while t(j+1) < tmax
    j = j + 1;
    a{j} = zeros(M+1,M);
    
    for m = 1:M
       a{j}(1,m) = exprnd(1/mu(m));
       for i = 2:(M+1)
          u = rand;
          tmp = 1 - exp(-(1/delta(i-1,m))*(lambda{1}(i-1,m)));
          if u < tmp
             a{j}(i,m) = -(1/delta(i-1,m))*log(1 + ((delta(i-1,m))/(lambda{j}(i-1,m)))*log(1-u));
          else
             a{j}(i,m) = Inf;
          end
       end
    end
    
    minimum = min(min(a{j}));
    t(j+1) = t(j) + minimum;
    
    [X,Z] = find(a{j}==minimum);
    mstar = Z;
    istar = X;
    lambda{j+1} = zeros(M,M);
    
    for m = 1:M
       if dist(mstar,m) == "Exp"
          ymstar = exprnd(1/(params{1}(mstar,m)));
       elseif dist(mstar,m) == "Gamma"      
          ymstar = gamrnd(params{2}(mstar,m), 1/(params{3}(mstar,m)));
       elseif dist(mstar,m) == "Normal" 
          ymstar = normrnd(params{4}(mstar,m), params{5}(mstar,m));      
       else
          error("Distribution not currently supported.") 
       end
       
       for i = 1:M
          lambda{j+1}(i,m) = (lambda{j}(i,m))*exp(-(delta(i,m))*(a{j}(istar,mstar))) + ymstar*(i==mstar);
       end
    N(m) = N(m) + 1*(m == mstar);
    end
    
    Nfull = [Nfull; N];
    k = Nfull(j+1, mstar);
    ttmp(k) = t(j+1);
 end
 
 % Create matrix of intensities
  itmp = zeros(j, M);
  for m = 1:M
     for tt = 1:j
     itmp(tt,m) = mu(m) + sum(lambda{tt}(:,m));
     end
  end
  
 % Create summary statistic to return
 t(end)=[];
 Nfull(end,:) = [];
 n = Nfull;
 r = ttmp(ttmp<tmax);
 intensity = itmp;
end

% To plot using this use stairs()