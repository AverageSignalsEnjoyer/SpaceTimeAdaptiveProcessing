%% 

clear all
close all

N = 16;
d = 1;
lambda = 2;
AoA = 0; 

angVec = -90:0.01:90;
N = 16;
d = 0.5;  
lambda = 1; 

[iGrid, angGrid] = ndgrid(0:N-1, angVec);

xn = exp(1i*2*pi*iGrid.*d/lambda.*sin(deg2rad(angGrid)));

sigmaJ = db2mag(50);
sigma  = 1;

% Build covariance matrix 
AoAJ = [-65 -40 -25 30 45 60];

R = eye(N);

for i = 1:length(AoAJ)

    sj = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoAJ(i))));
    R = R + sj'*sj*sigmaJ^2;

end

% Training sample covariance matrix 
sn = mvnrnd(zeros(1, N), R, N*2);
Rtrain = 1/(N*2) * sn'*sn;

sn = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoA)));

y = (inv(Rtrain)*sn')'*xn;
yCheb = (inv(Rtrain)*(sn'.*chebwin(N, 30)))'*xn;
yTrue = (inv(R)*sn')'*xn;

figure
hold on  
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(y)/(2*N)))
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yCheb)/(2*N)))
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yTrue)/N))

ylim([-100 0])

for i = 1:length(AoAJ)

    line([d/lambda*sin(deg2rad(AoAJ(i))) d/lambda*sin(deg2rad(AoAJ(i)))], [-100 0], 'Color', 'k')

end

%% Eigendecomposition

[V, D] = eig(Rtrain);

temp = diag(D);
temp(temp < 1e5) = 1;
D = diag(temp);

Rdecomp = V*D/V;

% Rdecomp = V(:, iMax) * D(iMax, iMax) / V(:, iMax);
% Rdecomp = Rdecomp + eye(N);

yDecomp = (inv(Rdecomp)*sn'.*chebwin(N, 30))'*xn;
yCheb = (inv(R)*(sn'.*chebwin(N, 30)))'*xn;

figure
hold on  
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yTrue)/N/2))
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yCheb)/N))