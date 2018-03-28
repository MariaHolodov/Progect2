%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maria and Ron project
% merch 2018
% Model simple simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%MODEL PARAMS
 PARAMS.ni=0.01;
 PARAMS.wn=1;%2*pi*fc;%???
 PARAMS.zeta=1;%0.7071 %have to be 0-2
 PARAMS.c=0.5;
 
 
N    = 20; %number of neurons in the net
%W   = ;%weight of the synapses
%x   = ;%neurons activity state

%create data input
H = hadamard(N)./sqrt(N);
randN = randi([1,N],1,4);
U = H(:,randN(1));
V = H(:,randN(2));


fi  = @(x) tanh(x); %non-linier function on x.
fi_a= @(X) arrayfun(fi,X);

tspan= [0 5000];
x0   = H(:,5);%ones(N,1);
W0   = randn(N,N);
%W0   = W0*W0';%to make the start matrix symetric
u1 = H(:,randN(3));
v1 = H(:,randN(4));
W0   = (u1*v1'-v1*u1')+4.*W0;%start from leared data+noise
y0   = fi_a(x0);
dydt0= zeros(N,1);
zStart = [x0 , y0 , dydt0, W0];

%no input- let the system some time to stabolize
[t,z] = ode45(@(t,z) model(t,z,N,U,V,PARAMS),tspan,zStart);

%whit input
tspan2= [5000 10000];
[t2,z2] = ode45(@(t2,z2) model(t2,z2,N,U,V,PARAMS),tspan,z(end,:));

%read results
realWEig = zeros(size(z,1),N);
imgWEig = zeros(size(z,1),N);
x   = z(:,1:N);
y   = z(:,(N+1):2*N);
dydt= z(:,(2*N+1):3*N);
for time = 1:size(z,1)
	W   = reshape(z(time,(3*N+1):end),N,N);%CHECK THIS ORDER!@$
    D = eig(W,'matrix');
    realWEig(time,:) = sort(real(diag(D)), 'descend');
    imgWEig(time,:)  = sort(imag(diag(D)), 'descend');
end

realWEig2 = zeros(size(z2,1),N);
imgWEig2 = zeros(size(z2,1),N);

x2   = z2(:,1:N);
y2   = z2(:,(N+1):2*N);
dydt2= z2(:,(2*N+1):3*N);

for time = 1:size(z2,1)
	W2   = reshape(z2(time,(3*N+1):end),N,N);%CHECK THIS ORDER!@$
    D2 = eig(W2,'matrix');
    realWEig2(time,:) = sort(real(diag(D2)), 'descend');
    imgWEig2(time,:)  = sort(imag(diag(D2)), 'descend');
end


%presentation of results ?
figure(1)
for eig = 1:N
    semilogx(t, realWEig(:,eig));
    hold on
end
xlabel('time')
ylabel('real(Eig(W))')
title('Befor input')

figure(2)
for eig = 1:N
    semilogx(t2, realWEig2(:,eig));
    hold on
end
xlabel('time')
ylabel('real(Eig(W))')
title('After input')

figure(3)
for eig = 1:N
    semilogx(t2, imgWEig2(:,eig));
    hold on
end
xlabel('time')
ylabel('imag(Eig(W))')
title('After input')

figure(4)
for eig = 1:N
    semilogx(t, imgWEig(:,eig));
    hold on
end
xlabel('time')
ylabel('imag(Eig(W))')
title('Before input')