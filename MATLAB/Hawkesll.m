function [ll] = Hawkesll(t, mu, alpha, beta, P)
%
% Hawkesll - log-likelihood of Hawkes process
%            currently works for P = 1 model and not others
%
% Inputs - t is the maturity time/event times
%        - mu is the background intensity vector
%        - alpha is the alpha matrix parameters
%        - beta is the beta matrix parameters
%        - P is the order of the kernel
%
% Outputs - log-likelihood of the Hawkes process
 
 M = length(mu);
 N = length(t);
 if P == 1
    if M == 1
       R = [];
       R(1) = 0;
       % Find first summation in likelihood
       tmp1 = log(mu + alpha*R(1));
       for k = 2:N
          R(k) = (exp(-beta*(t(k) - t(k - 1))))*(1 + R(k - 1));
          tmp1 = tmp1 + log(mu + alpha*R(k));
       end
     
       % Find second summation in likelihood
       T = max(t);
       tmp2 = 0;
       for k = 1:N
          tmp2 = tmp2 + 1 - exp(-beta*(T - t(k)));
       end
       tmp2 = (alpha/beta)*tmp2;
     
       ll = tmp1 - mu*T - tmp2;
    else
        lltmp = [];
        T = max(t, [], 'all');
        R = zeros(M, M);
        for m = 1:M
            tmp1 = 0;
            for k = 1:N
                if k == 1
                    tmp1 = log(mu(m));
                else
                    tmp2 = 0;
                    for i = 1:M
                        dtmp = 0;
                        for d = 1:N
                            if t(k-1,m) <= t(d,i) && t(d,i) < t(k,m)
                                dtmp = exp(-beta(i,m)*(t(k,m) - t(d,i)));
                            end
                        end
                        R(i,m) = exp(-beta(i,m)*(t(k,m) - t(k-1,m)))*R(i,m) + dtmp;
                        tmp2 = tmp2 + alpha(i,m)*R(i,m);
                    end
                    tmp1 = tmp1 + log(mu(m) + tmp2);
                end
            end
            
            tmp3 = 0;
            for i = 1:M
                for k = 1:N
                    tmp4 = alpha(i,m)/beta(i,m);
                    tmp5 = 1 - exp(-beta(i,m)*(T - t(k,i)));
                    tmp3 = tmp3 + tmp4*tmp5;
                end
            end
            lltmp(m) = tmp1 - mu(m)*T - tmp3;
        end
        ll = sum(lltmp, 'all');
    end
 else
    if M == 1
        T = max(t);
        R = zeros(1,P);
        tmp1 = 0;
        for k = 1:N
            if k == 1
                 tmp1 = log(mu);
            else
                tmp2 = 0;
                for j = 1:P
                   R(j) = (exp(-beta(j)*(t(k) - t(k-1))))*(1 + R(j));
                   tmp2 = tmp2 + alpha(j)*R(j);  
                end
                tmp1 = tmp1 + log(mu + tmp2);
            end
        end
        
        tmp3 = 0;
        for j = 1:P
            for k = 1:N
                tmp4 = alpha(j)/beta(j);
                tmp5 = 1 - exp(-beta(j)*(T - t(k)));
                tmp3 = tmp3 + tmp4*tmp5;
            end
        end
        ll = tmp1 - mu*T - tmp3;
    else
        lltmp = [];
        T = max(t, [], 'all');
        R = {};
            for j = 1:P
                R{j} = zeros(M,M);
            end
        for m = 1:M
            tmp1 = 0;
            for k = 1:N
                tmp2 = 0;
                if k == 1
                    tmp1 = log(mu(m));
                else
                    for i = 1:M
                        for j = 1:P
                            if i == m
                                R{j}(i,m) = (exp(-beta{j}(i,m)*(t(k,m) - t(k-1,m))))*(1 + R{j}(i,m));
                            else
                                dtmp = 0;
                                for d = 1:N
                                    if t(k-1,m) <= t(d,i) && t(d,i) < t(k,m)
                                        dtmp = exp(-beta{j}(i,m)*(t(k,m) - t(d,i)))*(1 + R{j}(i,m));
                                    end
                                end
                                R{j}(i,m) = (exp(-beta{j}(i,m)*(t(k,m) - t(k-1,m))))*(1 + R{j}(i,m)) + dtmp;
                            end
                            tmp2 = tmp2 + alpha{j}(i,m)*R{j}(i,m);
                        end
                    end
                    tmp1 = tmp1 + log(mu(m) + tmp2);
                end
            end
            
            tmp3 = 0;
            for i = 1:M
                for j = 1:P
                	for k = 1:N
                    	tmp4 = alpha{j}(i,m)/beta{j}(i,m);
                        tmp5 = 1 - exp(-beta{j}(i,m)*(T - t(k,i)));
                        tmp3 = tmp3 + tmp4*tmp5;
                	end
                end
            end 
            lltmp(m) = tmp1 - mu(m)*T - tmp3;
        end
        ll = sum(lltmp, 'all');
    end
 end
end
