function [t, Y, X, Z] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed)
%SimulateMarkedHawkesMD Simulate Multivariate Hawkes Process
%  Input:
%        T_stop = maturity time
%            mu = background intensity (D * 1 vector)
%            Y0 = intensity from edge effect (D * 1 vector)
%         delta = decay rates (D * D matrix)
%         distY = distribution of Y (intensity jump)
%          parY = parameters for Y (D * D * NY tensor)
%
%  Optional Input:
%         distX = distribution of X (marked)
%          parX = parameters for X (D * 1 * NX tensor)
%          seed = seed for random number generator
%
%  Output:
%             t = event times (D * 1 cell -> Nm * 1 vector)
%             X = marked      (D * 1 cell -> Nm * 1 vector)
%             Y = jump sizes  (D * 1 cell -> Nm * D matrix) 
%             Z = origins     (D * 1 cell -> Nm * 1 vector)
%
%
%  Note: 
%       1) Supported distY: 'const', 'gamma', 'exponential'
%       2) Supported distX: 'const', 'gamma', 'exponential', 'normal'
%       3) parX has additional column for marked from immigrant
%

if exist('seed','var')
    rng(seed);
else
    seed = int32((now - 7.361e5) * 1e6);
%     disp(['generated seed: ' num2str(seed)]);
    rng(seed);
end


%% Parameter Setting

% infer dimension from mu
D = length(mu);

N = 1000000; % pre-allocate array size


if D <= 0  %%% work for D = 1
    error('D must be at least 1');
end


%% Cache Variables

inv_mu    = 1 ./ mu;
inv_delta = 1 ./ delta;

% boolean on the sampling dist. of Y
Yconst = strcmp(distY, 'const');
Ygam   = strcmp(distY, 'gamma');
Yexp   = strcmp(distY, 'exponential');
Ygbm   = strcmp(distY, 'iidGBM');

% boolean on the sampling dist. of X
if exist('distX','var')
    Xconst = strcmp(distX, 'const');
    Xgam   = strcmp(distX, 'gamma');
    Xexp   = strcmp(distX, 'exponential');
    Xnorm  = strcmp(distX, 'normal');
else
    Xconst = false;
    Xgam   = false;
    Xexp   = false;
    Xnorm  = false;
end

% boolean on whether to store variables
saveZ = nargout >= 4;
saveX = nargout >= 3;
saveY = (nargout >= 2) || Ygbm;


% GBM for Y not implemented
if Ygbm
    error('GBM not implemented.');
end



%% Simulation

% record jump times (inc. both types)
r    = nan(N,1);
r(1) = 0;

% indicator variable showing which process a jump time belongs to
Ztype    = nan(N,1);

% indicator showing which process a jump time comes from
if saveZ
    Zfrom    = nan(N,1);
end

% record jump sizes (intensity jump)
if saveY
    Y_size = nan(N,D);
end

% record marked
if saveX
    X_size = nan(N,1);
end

% initialise lambda
lambda  = Y0;

% inter arrival time (cache)
s  = nan(D,D+1);


%% Simulation Loop for One Path
for j = 2:N
    
    % Sample jump times from background intensity (Dx1)
    % (Will give immigrant)
    s(:,1) = - inv_mu .* log(rand(D,1));
    
    
    % Sample jump times for mutually-exciting process (DxD)
    % (Will give offspring)
    Smi = 1 + delta ./ lambda .* log(rand(D));
    Smi(Smi < 0) = 0;  % remove invalid Smi
    
    Smi = - inv_delta .* log(Smi);
    
    % Combine jump times
    s(:,2:D+1) = Smi;
    
    
    % Find minimum
    [s_trunc, I1] = min(s, [], 1);
    [s_min, I2]   = min(s_trunc);
    
    % record event time
    r(j) = r(j-1) + s_min;
    
    i_star = I2 - 1;
    m_star = I1(I2);
    
    Ztype(j) = m_star;
    if saveZ
        Zfrom(j) = i_star;
    end
    
    
    % sample Y and update intensity
    YfromMstar = nan(D,1);
    if Yconst
        YfromMstar = parY(:,m_star);
    elseif Ygam
        YfromMstar = gamrnd(parY(:,m_star,1), parY(:,m_star,2));
    elseif Yexp
        YfromMstar = exprnd(parY(:,m_star));
    end
        
    if saveY
        Y_size(j,:) = YfromMstar;       % matlab auto transpose
    end
    
    lambda           = lambda .* exp(- delta * s_min);
    lambda(:,m_star) = lambda(:,m_star) + YfromMstar;
    
    
    % sample X
    if saveX
        if Xconst
            X_size(j) = parX(m_star,1);
        elseif Xgam
            X_size(j) = gamrnd(parX(m_star,1,1), parX(m_star,1,2));
        elseif Xexp
            X_size(j) = exprnd(parX(m_star,1));
        elseif Xnorm
            X_size(j) = normrnd(parX(m_star,1,1), parX(m_star,1,2));
        end
    end
    
    % terminate loop if reach maturity time
    if r(j) > T_stop
        r(j) = nan;  % set r_j = nan, thus removed later
        break
    end
    
    
    % give error if pre-allocated vector is too short
    if j == N
        msg = 'Pre-allocated vector is too short!';
        error(msg);
    end
    
end


% trim nan from the data
ind_nan = isnan(r);  % indices of nan

r(ind_nan) = [];
Ztype(ind_nan) = []; % needed
if saveZ
    Zfrom(ind_nan) = [];
end
if saveX
    X_size(ind_nan) = [];
end
if saveY
    Y_size(ind_nan,:) = [];
end



%% Retrieve event times for each Hawkes process

% event times
t = cell(D,1);
if saveZ
    Z = cell(D,1);
end
if saveX
    X = cell(D,1);
end
if saveY
    Y = cell(D,1);
end

for m = 1:D
    ind = Ztype == m;
    
    t{m} = r(ind);
    if saveZ
        Z{m} = Zfrom(ind,:);
    end
    if saveX
        X{m} = X_size(ind);
    end
    if saveY
        Y{m} = Y_size(ind,:);
    end
    
end

