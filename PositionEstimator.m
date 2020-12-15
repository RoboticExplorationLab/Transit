function Err_mag=PositionEstimator(TagVelocity,TagHeading,TagFreqError,TagLat,TagLong)
% This code estimates the position error of a hypothetical Doppler-based
% tracking system.
% Assumptions:
% 1) The satellite clock and pose are known exactly 
% 2) It assumes adequate SNR
% 3) The ERROR SOURCES are only: uncertainty in the velocity vector of the
% transmitter and uncertainty in the transmitter clock frequency
% 4) The error in measuring the transmitted signal at the satellite is not
% modeled
% 
% The inputs:





% This code works by: 
%
% First version by Shawn Swist
% AA 290 - Manchester
% Later modified by Manchester and MacCurdy


%Sim Parameters. If any of the inputs are left empty, these are
%substituted. If ALL inputs are left empty, this function acts like a
%script and generates plots.

%Tag velocity magnitude in meters/sec
Vtag_ms = 1.0; %default value
%Tag heading in radians
tag_heading = 2*pi*rand();  %choose a random heading in radians
%1-sigma white noise added to doppler measurements. This is the 1-sigma
%frequency extent of random frequency noise added to the Doppler
%measurements. This is assumed to come from noise sources in the tag clock
%that cause it to deviate from the 
freq_noise_hz = 1.0; 
% Tag position
GClat = 37.426622;
GClong = -122.173355;

%% Handle Inputs
if logical(exist('TagVelocity', 'var')) && ~isempty(TagVelocity) Vtag_ms=TagVelocity; end
if logical(exist('TagHeading', 'var')) && ~isempty(TagHeading) tag_heading=TagHeading; end
if logical(exist('TagFreqError', 'var')) && ~isempty(TagFreqError) freq_noise_hz=TagFreqError; end
if logical(exist('TagLat', 'var')) && ~isempty(TagLat) GClat=TagLat; end
if logical(exist('TagLong', 'var')) && ~isempty(TagLong) GClong=TagLong; end
arg_cnt=nargin;

%%
% =============== Constants ===============
% Earth Parameters 
rE = 6378;          % [km]       Earth radius
mu = 3.986e5;       % [km^3/s^2] Earth gravational parameter
%seed = rng(1);      % Keep the seed the same for randn

% Frequency of tag
f0 = 400;           % [MHz]      Frequency
c = 299792458;      % [m/s]      Speed of light

% Time
epoch = 58398;      % [MJD] 10/7/2018
tvec = linspace(.0215,.03,100);
dt = 86400*(tvec(2)-tvec(1)); %seconds


% TURN THESE INTO GEODETIC LAT/LONG! ????
GC_phi = 37.6123;
GC_phi = GClat;
GC_lam = GClong;

% Define satellite orbital elements
a = rE + 500;   % [km]  Semi-major axis
e = 0.01;       %  []   Eccentricity
i = 89;         % [deg] Inclination
Om = 90;        % [deg]    RAAN
w = 0;          % [deg] Argument of Periapsis
nu = deg2rad(0);         % [deg->rad] True anomaly
M = E2M(anom2E(nu,e),e);  % [rad] Mean anomaly
M = rad2deg(M);         % [deg] Redifine mean anomaly in degrees

% Note: E2M, anom2E fcns. use radians, SatProp uses degrees.
oe = [a; e; i; Om; w; M];


% Simulate Satellite motion
[r_ecef,v_ecef,r_enu,AZ,EL,X0tag,ENU] = SatellitePropagator(oe,epoch,mu,tvec+epoch,rE,GC_lam,GC_phi);

% Tag motion
Vtag = Vtag_ms*1e-3*(sin(tag_heading)*ENU(:,1) + cos(tag_heading)*ENU(:,2)); %km/s
Xtag = zeros(size(r_ecef));
Xtag(:,1) = X0tag;
for k = 2:length(Xtag)
    Xtag(:,k) = Xtag(:,k-1) + Vtag*dt;
end

% Find When Satellite is visible
viz = find(EL>10);
rdot = zeros(length(r_ecef),1);
for jj = 1:length(r_enu)
    Xsat = r_ecef(:,jj);
    Vsat = v_ecef(:,jj);
    R = Xsat-Xtag(:,jj);
    V = Vsat - Vtag;
    rdot(jj) = R'*V/sqrt(R'*R);
end


% Lets use less points and see how accurate ...
%viz = viz(100:250);

% Measurement function
X = r_ecef(:,viz);
V = v_ecef(:,viz);
H = @(P) diag((X-P)'*V)./sqrt(diag((X-P)'*(X-P)));
H2 = @(P) f0*(1+-(1000*H(P))/c);  % convert velocity from km/s to m/s
% 
% 
% % Numerical solution
% rmag = vecnorm(r_enu);
% drmag = diff(rmag(viz))./(diff(tvec(viz))*86400);
% drmag = drmag'; % make a column vector
% drmag(end+1) = drmag(end); % repeat last data point (numerical derivatives lose a data point)
% 
% % let the numerical solution be the real data to use...
% % Cost function
% data = drmag;
% data = h(viz) + 0.01*randn(length(viz),1); % add noise, 0 mean 0.01 std
% J = @(P) H(P) - data;
% 
% % Nonlinear least squares solver - velocity based
% opts = optimoptions('lsqnonlin','OptimalityTolerance',1e-15,'Display','off');
% pos = lsqnonlin(J,[rE;0;0],[],[],opts);
%     % seems to do decent with a bad guess. The numerical "true data" looks
%     % like it makes it miss based on the residuals.
% Elat = asind(pos(3)/norm(pos));
% Elong = rad2deg(atan2(pos(2),pos(1)));



% Doppler Shift
freq = f0*(1+-(1000*rdot(viz))/c); %in MHz
noise = freq_noise_hz*1e-6*randn(length(rdot(viz)),1); % sampled noise in MHz
fdata = freq + noise;

J2 = @(P) H2(P) - fdata;

% Nonlinear least squares solver - frequency based
opts = optimoptions('lsqnonlin','Display','off','OptimalityTolerance',1e-15,'FunctionTolerance',1e-15,'StepTolerance',1e-15);
pos2 = lsqnonlin(J2,[rE;0;0],[],[],opts);

Elat2 = asind(pos2(3)/norm(pos2));
Elong2 = rad2deg(atan2(pos2(2),pos2(1)));

err = 1000*(pos2-mean(Xtag,2)); %Zac, what's up with this 1000? Is this 
% meters or KM? The "disp" below indicates meters, but this 1000 makes it
% seem like km...
Err_mag = norm(err);

%% PLOTS

if arg_cnt<1
% plot the estimated locataion vs true location
figure; hold all
grid on
% Load and plot MATLAB built-in Earth topography data
load('topo.mat', 'topo');
topoplot = [topo(:, 181:360), topo(:, 1:180)];
contour(-180:179, -90:89, topoplot, [0, 0], 'black');
[T] = plot(GClong, GClat, 'bs');
%[E] = plot(Elong,Elat,'rp');
[E2] = plot(Elong2,Elat2,'g^');
legend([T,E2],'True Position','Estimated Position');
xlabel('Longitude'); ylabel('Latitude')
title('Estimated Location vs. True Location')


%% other plots
% figure;
% hold all
% earthPlot;
% plot3(r_ecef(1,viz),r_ecef(2,viz),r_ecef(3,viz))
% title('Satellite ECEF when Visible');

% figure;
% plot3(r_enu(1,viz),r_enu(2,viz),r_enu(3,viz))
% xlabel('E');ylabel('N');zlabel('U');
% title('View of sat from GS');

% figure;
% plot(tvec*24,EL)
% title('Elevation angle of satellite');

% figure;
%plot(tvec(viz)*24,drmag);
% plot(tvec(viz)*24,h(viz),':r')
% hold on
% plot(tvec(viz)*24,H(pos),'m-.')
% plot(tvec(viz)*24,data,'k:')
% title('Relative Velocity of Satellite');
% legend('h(p)','estimated pos','data used','location','northwest')


figure;
plot((tvec(viz)-tvec(viz(1)))*1440,1000*(freq-400));
hold on
plot((tvec(viz)-tvec(viz(1)))*1440,1000*(fdata-400),'r--');
xlabel('Time (min)');
ylabel('Doppler (kHz)');
title('Frequency obseved by satellite');
legend('True','Observed')

disp(['Error (meters): ',num2str(Err_mag)]);

end
