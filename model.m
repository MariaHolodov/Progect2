function dzdt = model(t,z,N,In1,In2,PARAMS)
 if ~exist('In1','var')
     %4 parameter does not exist, so default it to 0
      In1 = zeros(N,1);
 end
 if ~exist('In2','var')
     %5 parameter does not exist, so default it to 0
      In2 = zeros(N,1);
 end 
 ni   = PARAMS.ni;
 wn   = PARAMS.wn;
 zeta = PARAMS.zeta;
 c    = PARAMS.c;
 fi  = @(x) tanh(x); %non-linier function on x.
 
 
 %input
 Input = In1.*cos(t*wn )+In2*sin(t*wn );

dzdt= zeros(size(z));
x   = z(1:N,1);
y   = z((N+1):2*N,1);
dydt= z((2*N+1):3*N,1);
W   = reshape(z((3*N+1):end,1),N,N);%CHECK THIS ORDER!@$

fi_a= arrayfun(fi,x);

dzdt(1:N) = (-x+W*fi_a+Input); % dx/dt
dzdt((N+1):2*N) = dydt; % dy/dt
dzdt((2*N+1):3*N) = wn.*fi_a-(wn)^2.*y-2*zeta*wn.*dydt; % d(dy/dt)/dt
DWDT  = ni*(c.*eye(size(W))-fi_a*fi_a'+(fi_a*y'-y*fi_a')); % dw/dt
dzdt((3*N+1):end) = DWDT(:); 

