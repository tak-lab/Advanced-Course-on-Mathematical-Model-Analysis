function s = exp_taylor(x)
  s = 0;
  t = 1;
  i = 1;
  while 1
    s = s + t;
    if abs(t) < abs(s) * 1e-15, break; end
    t = (x / i)*t;
    i = i+1;
  end
end
