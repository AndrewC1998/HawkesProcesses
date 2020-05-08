function [t, Y, X, Z] = SimulateMarkedHawkes(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed)
%SimulateMarkedHawkes Simulate Multivariate Hawkes Process
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
%          parX = parameters for X (D * D+1 * NX tensor)
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

% infer dimension from mu
D = length(mu);

if nargout == 4
    if D == 0
        error('Dimension of mu is 0');
    elseif D == 1
        if exist('seed','var')
            [t, Y, X, Z] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t, Y, X, Z] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t, Y, X, Z] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY);
        end
    else
        if exist('seed','var')
            [t, Y, X, Z] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t, Y, X, Z] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t, Y, X, Z] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY);
        end
    end
end

if nargout == 3
    if D == 0
        error('Dimension of mu is 0');
    elseif D == 1
        if exist('seed','var')
            [t, Y, X] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t, Y, X] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t, Y, X] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY);
        end
    else
        if exist('seed','var')
            [t, Y, X] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t, Y, X] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t, Y, X] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY);
        end
    end
end

if nargout == 2
    if D == 0
        error('Dimension of mu is 0');
    elseif D == 1
        if exist('seed','var')
            [t, Y] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t, Y] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t, Y] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY);
        end
    else
        if exist('seed','var')
            [t, Y] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t, Y] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t, Y] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY);
        end
    end
end

if nargout == 1
    if D == 0
        error('Dimension of mu is 0');
    elseif D == 1
        if exist('seed','var')
            [t] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t] = SimulateMarkedHawkes1D(T_stop, mu, Y0, delta, distY, parY);
        end
    else
        if exist('seed','var')
            [t] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX, seed);
        elseif exist('distX','var')
            [t] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY, distX, parX);
        else
            [t] = SimulateMarkedHawkesMD(T_stop, mu, Y0, delta, distY, parY);
        end
    end
end
