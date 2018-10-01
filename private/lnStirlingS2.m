function lnS = lnStirlingS2(n,k)
% This function calcultates the log of sterling number, which is the number
% of all possible partitions of assigning n samples into k nonempty clusters
% borrowed from Jie Xiong's Code
%
% Lu Cheng
% 23.01.2012

lnS = 0;
if k > n;
    error(' k < n in lnStirlingS2');
elseif n==0 && k == 0
    lnS = 0;
elseif k == 0
    lnS = -Inf;
elseif k == 1 || k == n 
    lnS = 0;
elseif k==2 
      if n < 50
          lnS = log(2^(n-1) - 1);
      else
          lnS = (n-1)*0.693147;
      end
elseif k == n-1 
    lnS = log(StirlingS2(n,k));
else
    nk = n/k;
    x0 = 1;
    for i = 0:10000
        x = (x0 + nk - nk*exp(-x0))/2;
        if (x-x0)/(1+x) < 1e-10
            break
        end
        x0 = x;
    end
    t0 = n/k -1;
    if (x < 100)
        A = -n * log(x) + k*log(exp(x)-1);
    else
        A = -n * log(x) + k*x;
    end
    A = A + -k*t0 +(n-k)*log(t0);
    lnS = A + (n-k)*log(k) + 0.5*log(t0/((1 + t0)*(x-t0)));
    [c,scalec] = Binomial(n,k);
    lnS = lnS + log(c) + scalec;
end
        

function StirS2 = StirlingS2(n,k)
S = zeros(1,k);
S(1:2) = 1;
if (n == 0 && k == 0)||k == 1||k == n
    StirS2 = 1;
elseif k == 0 || k > n
    StirS2 = 0;
elseif k == 2
    StirS2 = 2^(n-1) - 1;
elseif k == n -1
    StirS2 = n*(n-1)/2;
else
    for i = 3:n
        for j = min(k,i):-1:2
            S(j) = S(j-1) + j*S(j);
        end
    end
    StirS2 = S(k);
end


function [c scale] = Binomial(n,k)
c = 1;
large = 1e99;
scale = 0;
if floor(k) ~= k
    disp('k is not a interger in Binomial.');    
end
if n < 0 && mod(k,2) == 1
    c = -1;
end
if k == 0
    c = 1;
end
if n > 0 && floor(n) == n
    k = min(k,n-k);
end
for i = 1:k
    c = c*(n-k+i)/i;
    if c > large
        scale = scale + log(c);
        c = 1;
    end
end

        
        
