function runSimulationAndPresent(N,U,V,PARAMS,maxtime,zStart,print_data)


tspan= [0 5000];%maxtime
[t,z] = ode45(@(t,z) model(t,z,N,U,V,PARAMS),tspan,zStart);

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

%presentation of results ?
figure()
for eig_n = 1:N
    semilogx(t, realWEig(:,eig_n));
    hold on
end
xlabel('time')
ylabel('real(Eig(W))')
title(sprintf('%s input',print_data))

figure()
for eig_n = 1:N
    semilogx(t, imgWEig(:,eig_n));
    hold on
end
xlabel('time')
ylabel('imag(Eig(W))')
title(sprintf('%s input',print_data))

end