function cout = ChebDerCoeffs(c,rigorous)
% NOTE: remind dividing by the rescaling factor (tmax/2)!
%       If rigous > 0, then this uses the interval arithmetic.
    [n, m] = size(c);
    if rigorous>0
      cout = intval(zeros(n-1, m));
    else
      cout = (zeros(n-1, m));
    end
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*c(2:end,:);
    cout(n-1:-2:1,:) = vcumsum(v(n-1:-2:1,:));
    cout(n-2:-2:1,:) = vcumsum(v(n-2:-2:1,:));
    cout(1,:) = .5*cout(1,:);
end

function csum = vcumsum(V)
n=size(V,1);
csum=tril(ones(n))*V;
end

