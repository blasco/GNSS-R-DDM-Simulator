freq  = 1575.42e6;
h     = 540e3; 
h0    = 20e6; 
gamma = 70*pi/180; 
Vt    = [-1e3,1e3,0]; 
Vr    = [500,1000,1e2]; %Vr_x ~= 0

lambda = 3e8/freq;
Re     = 6371e3;

r      = [100, -57];
r(3)   = sqrt(Re^2 - r(1)^2 - r(2)^2) - Re;

Rt = [0, h0/tan(gamma), h0]; 
Rr = [0, -h/tan(gamma), h]; 

ni = -(Rt-r)/norm(Rt-r); 
ns = (Rr-r)/norm(Rr-r); 

tau = norm(Rr-r) + norm(r - Rt) - (norm(Rr) + norm(Rt)); 
fd = 1/lambda*(dot(Vt,ni) - dot(Vr,ns)); 

[x1, x2, y1, y2] = myDD2SRF_v1(tau, fd, h, gamma, Vt, Vr, lambda);