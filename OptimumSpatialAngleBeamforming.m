%% 

clear all
close all

%% Simple beamforming example

% Beamforming example 
% Signal hits the array at the angle of arrival. The optimal beamformer
% weights derived to maximize response at the angle of arrival is the
% simulated signal itself (see textbook for details, J.R. Guerci 2015)

N      = 16; % Number of elements in uniform linear array (ULA)
d      = 1;  % Spacing between elements in ULA
lambda = 2;  % Opearting wavelength
AoA    = 30; % Angle of arrival

sn = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoA))); % Signal received
wn = sn; % Optimal beamformer weights

angVec = -90:0.01:90;
[iGrid, angGrid] = ndgrid(0:N-1, angVec);

% Output of the N receive channels 
xn = exp(1i*2*pi*iGrid.*d/lambda.*sin(deg2rad(angGrid))); 

y = wn*xn; % Linear combination of the outputs to do beamforming

figure
box on
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(y)/N))
ylim([-50 0])
xlabel('Normalized angle (d/\lambda sin \theta)')
ylabel('Normalized amplitude (dB)')

%% Chebyshev windowing 

% Apply Chebyshev window to demonstrate sidelobe suppression at the cost of
% widening the mainlobe 

figure
hold on
box on
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(y)/N))

yCheb40 = wn.*chebwin(N, 40)'*xn;
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yCheb40)/N))

yCheb60 = wn.*chebwin(N, 60)'*xn;
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yCheb60)/N))

ylim([-100 0])

xlabel('Normalized angle (d/\lambda sin \theta)')
ylabel('Normalized amplitude (dB)')

legend('No window', '40 dB Cheb', '60 dB Cheb', 'Location', 'Northwest')

%% Simulate jamming 

% Simulate jamming at various angles of arrivals. Beamformer nulls can be
% steered to suppress jamming signal

sigmaJ = db2mag(50); % Signal of the jammer

AoA = 0; 
sn = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoA)));
wn = sn; 

AoAJ = [-65 -40 -25 30 45 60]; % Jammer arrives at these angles

% Calculate noise covariance matrix 
R = eye(N);
for i = 1:length(AoAJ)
    sj = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoAJ(i))));
    R = R + sj'*sj*sigmaJ^2;
end

y = (inv(R)*wn')'*xn;
figure
hold on
box on
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(y)/N))

for i = 1:length(AoAJ)

    line([d/lambda*sin(deg2rad(AoAJ(i))) d/lambda*sin(deg2rad(AoAJ(i)))], [-70 0], 'Color', 'k')

end

xlabel('Normalized angle (d/\lambda sin \theta)')
ylabel('Normalized amplitude (dB)')
ylim([-70 0])
legend('Response', 'Jamming')

%% Simulate jamming with Chebyshev window 

sigmaJ = db2mag(50);
sigma  = 1;

AoA = 0; 
sn = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoA)));
wn = sn;

AoAJ = [-65 -40 -25 30 45 60];

R = eye(N);

for i = 1:length(AoAJ)
    sj = exp(-1i*2*pi*(0:N-1)*d/lambda*sin(deg2rad(AoAJ(i))));
    R = R + sj'*sj*sigmaJ^2;
end

y = (inv(R)*sn')'*xn;
yCheb = (inv(R)*(sn'.*chebwin(N, 30)))'*xn;

figure
hold on
box on

plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(y)/N))
plot(d/lambda*sin(deg2rad(angVec)), mag2db(abs(yCheb)/N))

for i = 1:length(AoAJ)
    line([d/lambda*sin(deg2rad(AoAJ(i))) d/lambda*sin(deg2rad(AoAJ(i)))], [-100 0], 'Color', 'k')
end

ylim([-100 0])

xlabel('Normalized angle (d/\lambda sin \theta)')
ylabel('Normalized amplitude (dB)')
legend('No window', '60 dB Cheb', 'Location', 'Northwest')