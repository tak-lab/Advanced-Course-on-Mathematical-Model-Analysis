function dy = LotkaVolterra(t,y,a,b,c,d)
dy = zeros(2,1);
dy(1) = y(1).*(a-b*y(2));
dy(2) = y(2).*(-c+d*y(1));