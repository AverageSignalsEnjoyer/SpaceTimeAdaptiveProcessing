%%

clear all
close all

N = 16;
M = 16;
d = 1;
lambda = 2;
AoA = deg2rad(0);
v = 0.8;
T = 1;
fd = 2*v*T/lambda*sin(AoA);

angVec = deg2rad(-90:0.5:90);
% 

fdVec = -0.5:0.01:0.5;
fdVec  = 2*v*T/lambda*sin(angVec);

[iGrid, angGrid] = ndgrid(0:N-1, angVec);

A = exp(1i*2*pi.*iGrid.*d/lambda.*sin(angGrid));

[iGrid, fdGrid] = ndgrid(0:M-1, fdVec);

B = exp(1i*2*pi.*iGrid.*fdGrid);

y = zeros(length(angVec), length(fdVec));

a = exp(1i*2*pi*(0:N-1)*d/lambda*sin(AoA));
b = exp(1i*2*pi*(0:M-1)*0.25);

for i = 1:length(angVec)
    for j = 1:length(fdVec)

        temp = kron(A(:,i), B(:,j));
        y(i, j) = kron(a,b) * temp;

    end
end

imagesc(d/lambda*sin(angVec), fdVec, mag2db(abs(y')/(N*M)))
axis square
colorbar
caxis([-60 0])

%% Build clutter and jamming covariance matrix 

AoAJ = deg2rad([-65 -40]);

R = eye(N*M);
sigmaJ = db2mag(50);

for ang = AoAJ

    a = exp(1i*2*pi*(0:N-1)*d/lambda*sin(ang));
    b = exp(1i*2*pi*(0:M-1)*0);

    temp = kron(a,b);
    R = R + sigmaJ*temp*temp';

end

fdVec = -0.5:0.01:0.5;
% fdVec  = 2*v*T/lambda*sin(angVec);

[iGrid, angGrid] = ndgrid(0:N-1, angVec);

A = exp(1i*2*pi.*iGrid.*d/lambda.*sin(angGrid));

[iGrid, fdGrid] = ndgrid(0:M-1, fdVec);

B = exp(1i*2*pi.*iGrid.*fdGrid);

y = zeros(length(angVec), length(fdVec));

a = exp(1i*2*pi*(0:N-1)*d/lambda*sin(AoA));
b = exp(1i*2*pi*(0:M-1)*0.25);

Rinv = inv(R);

for i = 1:length(angVec)
    for j = 1:length(fdVec)

        temp = kron(A(:,i), B(:,j));
        y(i, j) = kron(a,b) * Rinv * temp;

    end
end

figure
imagesc(d/lambda*sin(angVec), fdVec, mag2db(abs(y')/(N*M)))
axis square
colorbar
caxis([-60 0])
