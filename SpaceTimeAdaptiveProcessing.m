%%

clear all
close all

N = 16; % Number of elements in ULA
M = 16; % Number of pulses? 
lambda = 2; % Operating wavelength
d = lambda/2; % Distance of elements in ULA 
v = d;   % Speed of the aircraft? Example in textbook didn't specify
T = v/2; % Duration of the pulse? 
Nclutter = 200; % How many clutter patches

% Target doppler and AoA 
fd = 0.25;
AoA = deg2rad(0);

% Resolution of clutter patches 
% The response of the beamformer is plotted against these angles and
% dopplers
angVec = d/lambda*sin(deg2rad(linspace(-90, 90, Nclutter))); % norm angle
beta   = 2*v*T/d;
fdVec  = beta*angVec;

% Spatial steering vectors 
[iGrid, angGrid] = ndgrid(0:N-1, angVec);
A = exp(1i*2*pi.*iGrid.*angGrid);

% Doppler steering vectors 
[iGrid, fdGrid] = ndgrid(0:M-1, fdVec);
B = exp(1i*2*pi.*iGrid.*fdGrid);

% Space-time steering vector 
V = zeros(N*M, Nclutter);
for iClutter = 1:Nclutter
    V(:, iClutter) = kron(A(:, iClutter), B(:, iClutter));
end

% Response to a target 
[~, iMin] = min(abs(angVec - AoA));
a = A(:, iMin);
[~, iMin] = min(abs(fdVec - fd));
b = B(:, iMin);

% Covariance matrix 
R = eye(N*M);

% Clutter covariance matrix
R_clutter = V*V';

% Jammer covariance matrix 
AoA_jammer = [32 126];

R_jammer = zeros(N*M);

for ang = AoA_jammer
    
    [~, iMin] = min(abs(angVec - d/lambda*sin(deg2rad(ang))));
    
    bj = B(:, iMin);
    aj = A;
    
    s = kron(bj, aj);
    
    R_jammer = R_jammer + db2mag(50)^2*(s*s');

end

Rinv = inv(R);

% An inefficient, but perhaps more informative way to do the code below 
% Rinv = eye(N*M);

% Y = zeros(length(angVec));
% for i = 1:length(angVec)
%     for j = 1:length(fdVec)
% 
%         x = kron(A(:, i), B(:, j));
%         s = kron(a, b);
%         w = s;
%         % w = Rinv*s;
% 
%         Y(i, j) = w'*x;
%     end
% end

W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*X;
Y = reshape(Y, [Nclutter Nclutter]);

% No jammer or clutter 
subplot(2,2,1)

Rinv = inv(R);
W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*X;
Y = reshape(Y, [Nclutter Nclutter]);

s = surf(angVec, fdVec, mag2db(abs(Y')/(N*M)));
s.EdgeColor = 'none';
view(-90, 90)
colorbar
clim([-60 0])
ylabel('Angle')
xlabel('Doppler')

% Clutter
subplot(2,2,2)

Rinv = inv(R + R_clutter);
W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*X;
Y = reshape(Y, [Nclutter Nclutter]);

s = surf(angVec, fdVec, mag2db(abs(Y')/(N*M)));
s.EdgeColor = 'none';
view(-90, 90)
colorbar
clim([-60 0])
ylabel('Angle')
xlabel('Doppler')
title('Clutter Suppression')

% Clutter and jammer
subplot(2,2,3)

Rinv = inv(R + R_clutter + R_jammer);
W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*X;
Y = reshape(Y, [Nclutter Nclutter]);

s = surf(angVec, fdVec, mag2db(abs(Y')/(N*M)));
s.EdgeColor = 'none';
view(-90, 90)
colorbar
clim([-60 0])
ylabel('Angle')
xlabel('Doppler')
title('Clutter and Jammer Suppression')

% Jammer
subplot(2,2,4)

Rinv = inv(R + R_jammer);
W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*X;
Y = reshape(Y, [Nclutter Nclutter]);

s = surf(angVec, fdVec, mag2db(abs(Y')/(N*M)));
s.EdgeColor = 'none';
view(-90, 90)
colorbar
clim([-60 0])
ylabel('Angle')
xlabel('Doppler')
title('Jammer Suppression')

saveas(gcf, 'stap.png')

%% Space-time Chebyshev windowing 

figure 
hold on

subplot(1,2,1)

Rinv = inv(R);
W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*X;
Y = reshape(Y, [Nclutter Nclutter]);

s = surf(angVec, fdVec, mag2db(abs(Y')/(N*M)));
s.EdgeColor = 'none';
view(-90, 90)
colorbar
clim([-60 0])
ylabel('Angle')
xlabel('Doppler')

% Chebyshev windowing
subplot(1,2,2)

Rinv = inv(R);
W = kron(A, B);
X = kron(a, b);

Y = W'*Rinv*(X .* kron(chebwin(N, 30), chebwin(N, 30)));
Y = reshape(Y, [Nclutter Nclutter]);

s = surf(angVec, fdVec, mag2db(abs(Y')/(N*M)));
s.EdgeColor = 'none';
view(-90, 90)
colorbar
clim([-60 0])
ylabel('Angle')
xlabel('Doppler')

