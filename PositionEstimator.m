function Err_mag=PositionEstimator(TagVelocity,TagHeading,TagFreqError,TagLat,TagLong)
% This code estimates the position error of a transmitter using a 
% hypothetical Doppler-based tracking system.
% 
% Assumptions:
% 1) The satellite clock and pose are known exactly 
% 2) It assumes adequate SNR
% 3) The ERROR SOURCES are only: uncertainty in the velocity vector of the
% transmitter and uncertainty in the transmitter clock frequency
% 4) The error in measuring the transmitted signal at the satellite is not
% modeled
% 
% The inputs:
% TagVelocity - Tag velocity magnitude in meters/sec
% TagHeading - Tag heading in radians
% TagFreqError - 1-sigma value of white noise added to doppler measurements. This is 
% the 1-sigma frequency extent of random frequency noise added to the Doppler
% measurements. This is assumed to come from noise sources in the tag clock
% that cause it to deviate from the set point
% TagLat and TagLong are the true positions of the tag
% 
% This code works by: Simulating the motion of a satellite (SV) with given
% orbital parameters. It uses a prescribed tag position and finds the
% subset of SV locations in earth-centered earth-fixed coordinates
% that are in view of the tag (when SV is at least 10 deg above horizon). 
% For these positions, the range and velocity between the tag and SV are
% computed. These are the "true" ranges and relative velocities. The dot
% product between these tuples yields the rate of approach/retreat, or the
% scalar relative velocity that appears in the Doppler computation. With
% these approach/retreat speeds, we then form a vector of apparent (to the 
% SV) tag frequencies. We add Gaussian noise to these frequency 
% "measurements", to account for tag Tx frequency error and SV frequency
% estimation error. This array of frequency "measurements" form our
% observables. We then use the equations that describe the approach/retreat
% speeds and the Doppler shift to create a measurement function (called H2
% below) for the true Doppler shift, given a single position and the series
% of known SV positions/velocities. The difference between this 
% measurement function and the data is minimized in a least-squares sense,
% yielding an estimate for the tag position.
%
%
% First version by Shawn Swist as part of AA 290 taught by Manchester
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
f0 = 400;           % [MHz]      Center Frequency
c = 299792458;      % [m/s]      Speed of light

% Time
epoch = 58398;      % [MJD] 10/7/2018
tvec = linspace(.0215,.03,100); %choose a set of times (in MJD) over which 
    %we'll simulate the motion of the satellite
dt = 86400*(tvec(2)-tvec(1)); %seconds between simulation steps


% TURN THESE INTO GEODETIC LAT/LONG! ????
%GC_phi = 37.6123;
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
% This finds the ECEF coordinates in km of the satellite given the orbital
% parameters and a vector of times in MJD. It also gives the tag/transmitter
% position in ECEF.
[r_ecef,v_ecef,r_enu,AZ,EL,X0tag,ENU] = SatellitePropagator(oe,epoch,mu,tvec+epoch,rE,GC_lam,GC_phi);

% Tag motion
% Assume the tag starts from its initial position and moves with constant
% velocity "Vtag_ms" in the direction prescribed by "tag_heading"
Vtag = Vtag_ms*1e-3*(sin(tag_heading)*ENU(:,1) + cos(tag_heading)*ENU(:,2)); %km/s
Xtag = zeros(size(r_ecef));
Xtag(:,1) = X0tag;
for k = 2:length(Xtag)
    Xtag(:,k) = Xtag(:,k-1) + Vtag*dt;
end

% Find the entries of the satellite motion in which it is visible to the 
% tag, indicated when the elevation is above 10 degrees
viz = find(EL>10);

% Now compute 
rdot = zeros(length(r_ecef),1);
for jj = 1:length(r_enu)
    Xsat = r_ecef(:,jj);    %pos of sat, km
    Vsat = v_ecef(:,jj);    %velocity of sat, km/sec
    R = Xsat-Xtag(:,jj);    %sat to tag range
    V = Vsat - Vtag;        %sat to tag relative velocity
    rdot(jj) = R'*V/sqrt(R'*R); %compute the dot product of the velocity 
                                %vector with the unit vector pointed along 
                                %the direction of R. This is the relative
                                %velocity or rate of approach/retreat that
                                %appears in the Doppler computations.
                                %km/sec
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
freq = f0*(1+-(1000*rdot(viz))/c); %Apparent tag freq, in MHz, after Doppler shift
noise = freq_noise_hz*1e-6*randn(length(rdot(viz)),1); % Add Gaussian noise 
                        %with a 1-sigma magnitude scaled by "freq_noise_hz"
        
fdata = freq + noise;   %Apparent tag freq plus frequency noise, after the 
                        %Doppler shift. This is our observable "data" with
                        %which the estimator needs to predict the tag
                        %location.

J2 = @(P) H2(P) - fdata;

% Nonlinear least squares solver - frequency based. Find the single
% position that best explains, in a LS sense, the multiple Doppler measurements
opts = optimoptions('lsqnonlin','Display','off','OptimalityTolerance',1e-15,'FunctionTolerance',1e-15,'StepTolerance',1e-15);
pos2 = lsqnonlin(J2,[rE;0;0],[],[],opts);

%pull out the Lat/Long locations of the estimated position, for plotting
Elat2 = asind(pos2(3)/norm(pos2));
Elong2 = rad2deg(atan2(pos2(2),pos2(1)));

err = 1000*(pos2-mean(Xtag,2)); %compare the difference between the 
                                %estimated tag position and the mean of the
                                %prescribed positions (since the tag moves 
                                %at a constant velocity). Convert from km
                                %to m to get the true error of our estimate.
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
