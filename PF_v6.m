% Particle filter algorithm for a 2 dimensional multiple source array
% amplitude estimation using a 2 dimensional array of microphones and
% normal surface velocity for the actual source surface

% main function
function [xg, yg, X_hat, U, u_plot, u_hat_plot, error_amp, error_vel] = PF_v6()

% parameters and initial conditions
fM = 99; % number of frames
M = 7; % number of microphones
L = 7; % number of sources
N = 100; % number of particles
N_thr = N/10; % threshold sample size
N_eff = zeros(1,(L^2)); % initialize effective sample size array
c = 343; % speed of sound in m/s
rho = 1.2754; % density of air in kg/m^3
pCov = 0.1; % initial covariance of particle distribution
v1Cov = 0.1; % process noise covariance
v2Cov = 0.1; % measurement noise covariance
freq = 100:100:1700; % defining frequency vector
gr = 21; % grid points needed for contour plots of surface velocity
x = 100*ones((L^2), 1); % x = [1 2 3 4 5]'.*ones(L, 1) ; % initial state vector (amplitudes of the sources)

% true states initialization
X = zeros((L^2), fM+1); % initialize sequence of true states
X(:,1) = x; % true state vector at frame 1 is just the initial state vector
U = zeros(gr, gr, fM); % initialize sequence of true velocity vectors
sum_vel = 0;

% particles and estimates initialization
xn = zeros((L^2),N);
for q = 1:N
    xn(:,q) = x; % array of the initial states used later for initializing particles
end
p = normrnd(xn,pCov,(L^2),N); % initializing N particles
xp = zeros((L^2),N); % initialize state estimate
x_hat = zeros((L^2),1); % initial estimated state vector (amplitudes of the sources)
X_hat = zeros((L^2), fM+1); % intialize sequence of state estimates
X_hat(:,1) = x; % estimated state vector at the frame 1 is just the initial state vector
U_hat = zeros(gr, gr, fM); % initialize sequence of estimated velocity vectors
sum_vel_hat = 0;

% initialize particle weights
w = zeros((L^2),fM+1,N);
for l = 1:(L^2)
    for k=1:N
        for fm = 1:fM
            if fm == 1,
                w(l,fm,k) = 1/N; % Equally normalized weights for all the particles at frame 1
            else w(l,fm,k) = 0;
            end
        end
    end
end
v = zeros((L^2),N);

% distance calculation between the source and microphones
xsv = linspace(-0.3, 0.3, L); % x coordinates of the L^2 sources
ysv = linspace(-0.3, 0.3, L); % x coordinates of the L^2 sources
% [xsg,ysg] = meshgrid(xsv,ysv); % defining grid points for L^2 sources
zs = 0.4; % z coordinates of the L^2 sources
xmv = linspace(-0.3, 0.3, M); % x coordinates of the M^2 microphones
ymv = linspace(-0.3, 0.3, M); % y coordinates of the M^2 microphones
[xmg,ymg] = meshgrid(xmv,ymv); % defining grid points for M^2 microphones
zm = 0; % z coordinates of the M^2 sources
rsm = zeros(M, M, L, L);
r = zeros((M^2),(L^2));
id = 1;
for i = 1:L
    for j = 1:L
        rsm(:,:,i,j) = (((xsv(i)-xmg).^2)+((ysv(j)-ymg).^2)+((zs-zm)^2)).^(1/2);
        r(:,id) = reshape(rsm(:,:,i,j),(M^2),1);
        id = id+1;
    end
end

% defining grid points
xv = linspace(-0.5, 0.5, gr);
yv = linspace(-0.5, 0.5, gr);
[xg,yg] = meshgrid(xv,yv); % defining grid points for velocity
zg = 0.35;

% defining distance matrix: distances of grid points from the sources
rsg = zeros(gr, gr, L, L);
rg = zeros((gr^2),(L^2));
id = 1;
for i = 1:L
    for j = 1:L
        rsg(:,:,i,j) = (((xsv(i)-xg).^2)+((ysv(j)-yg).^2)+((zs-zg)^2)).^(1/2);
        rg(:,id) = reshape(rsg(:,:,i,j),(gr^2),1);
        id = id+1;
    end
end

% error matrices initialization
error_amp = zeros(length(freq), fM+1); % initializing the amplitude error matrix
error_vel = zeros(gr, gr, length(freq), fM+1); % initializing the velocity error matrix
% error = zeros(length(freq), L, fM+1);

% initializing covariance matrix
Cxx = zeros((L^2),(L^2),N);
for k = 1:N
    Cxx(:,:,k) = pCov*eye(L^2); % covariance matrix of state vector vs. state vector
end

Uc = zeros((gr^2), length(freq), fM);
Uc_hat = zeros((gr^2), length(freq), fM);
fid=1; % frquency index in error matrix

% main loop
for freq = 100:100:1700 % for different frequencies between 100 Hz to 1700 Hz
    G = (exp((-1i*(2*pi*freq*r))/c))./r; % (M^2)*(L^2) Green's matrix used in measurement equation
    for fm=1:fM % for different frames upto fM
        % nonlinear discrete system in frequency domain with additive noise
        % state space equations
        x = x + normrnd(0, v1Cov); % process equation: true state at frame fm
        y = G*x + normrnd(0, v2Cov); % measurement equation: pressure on M microphones at frame fm corresponding to observations
        
        % true normal surface velocity calculation
        ps = 0; % initializing sound pressure on actual sound surface
        id = 1;
        for i=1:L
            for j = 1:L
                ps = ps + x(id)*((exp((-1i*(2*pi*freq*rsg(:,:,i,j)))/c))./rsg(:,:,i,j));
                id=id+1;
            end
        end
        u = 0; % initializing true normal surface velocity on the actual surface
        for i=1:L
            for j=1:L
                u = u + (1/(1i*rho*freq))*(((zs-zg)./rsg(:,:,i,j)).*((1i*2*pi*freq/c)+(1./rsg(:,:,i,j))).*ps);
            end
        end
        U(:,:,fm) = u; % true normal surface velocity matrix at frame fm
        Uc(:,fid,fm) = reshape(u, gr^2, 1);
        
        % particle filtering
        for j = 1:N
            xp(:,j) = p(:,j) + normrnd(0, v1Cov); % particle prediction
            yp = G*xp(:,j) + normrnd(0, v2Cov); % prediction measurement
            %prediction performance
            dx = xp(:,j) - x; % difference between particle prediction and true state vector
            dy = yp - y; % difference between prediction measurement and observation
            % Covariance matrices
            Cxx(:,:,j) = Cxx(:,:,j) + v1Cov*eye(L^2);
            Cyx = (G*Cxx(:,:,j));
            Cxy = Cyx';
            Cyy = G*Cxx(:,:,j)*G' + v2Cov*eye(M^2);
            temp = zeros((L^2),(L^2)); % temporary L*L matrix for inverse calculation later
            for g=1:(L^2)
                for h=1:(L^2)
                    temp(g,h) = Cxx(g,h,j);
                end
            end
            C = Cyy - Cyx*(inv(temp))*Cxy; % final covariance matrix used in likelihood function
            for l=1:(L^2)
                lf = abs((1/((((2*pi)^(M^2))*det(C))^(1/2)))*exp((-1/2)*dy'*(inv(C))*dy)); % likelihood function
                pd = abs((1/(((2*pi)*(pCov^(L^2)))^(1/2)))*exp((-1/2)*dx'*(inv(pCov*eye(L^2)))*dx)); % prior distribution
                w(l,fm+1,j) = w(l,fm,j)*lf*pd; % weight update equation
            end
            % w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-((abs(d_avg))^2)/(2*v2Cov)));
            % page 324 of steven kay's book
            % *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
        end
        
        % Weight normalization and state estimation
        for l=1:(L^2)
            w(l,fm+1,:) = w(l,fm+1,:)./(sum(w(l,fm+1,:))); % normalize the likelihood of each a priori estimate
            for k = 1:N
                v(l,k) = w(l,fm+1,k);
            end
            x_hat(l) = v(l,:)*p(l,:)'; % weighted sum: the state estimate is the weighted mean of the particles
        end
        
        % Estimated normal surface velocity calculation
        ps_hat = 0; % initializing sound pressure estimate on actual sound surface
        id=1;
        for i=1:L
            for j=1:L
                ps_hat = ps_hat + x_hat(id)*((exp((-1i*(2*pi*freq*rsg(:,:,i,j)))/c))./rsg(:,:,i,j));
                id=id+1;
            end
        end
        u_hat = 0; % initializing estimated normal surface velocity on the actual surface
        for i=1:L
            for j=1:L
                u_hat = u_hat + (1/(1i*rho*freq))*(((zs-zg)./rsg(:,:,i,j)).*((1i*2*pi*freq/c)+(1./rsg(:,:,i,j))).*ps_hat);
            end
        end
        U_hat(:,:,fm) = u_hat; % estimated normal surface velocity matrix at frame fm
        Uc_hat(:,fid,fm) = reshape(u_hat, gr^2, 1);
        
        X(:,fm+1) = x; X_hat(:,fm+1) = x_hat; % update sequence of true states and state estimates
        
        % Error matrices for source amplitude and normal surface velocity
        error_amp(fid,fm+1) = (norm(X_hat(:,fm+1)-X(:,fm+1))/norm(X(:,fm+1)))*100; % source amplitude error matrix
       
        % Sum of true and estimated normal surface velocity matrices over all frames 
        sum_vel = sum_vel + sum(U(:,:,fm));
        sum_vel_hat = sum_vel_hat + sum(U_hat(:,:,fm));
        
        % Resampling
        for l=1:(L^2)
            N_eff(l) = 1/(sum((w(l,fm+1,:).^2))); % effective sample size
            if N_eff(l)<N_thr % Resampling only required when effective sample size is less than threshold sample size
                [p(l,:), v(l,:)] = resample(xp(l,:),v(l,:)); % Resampling
                for k = 1:N
                    w(l,fm+1,k) = v(l,k);
                end
            else p=xp;
            end
        end % End of resampling loop
    end % End of frame loop
    U_mean = sum_vel/fM;
    U_hat_mean = sum_vel_hat/fM;
    error_vel = ((U_hat_mean-U_mean)./U_mean).*100;
    fid=fid+1; % frequency index increment
end % end of main (frequency) loop

% IFFT
np=128;
u_t = zeros((gr^2), np);
u_hat_t = zeros((gr^2), np);

u_t_RMS = zeros((gr^2), 1);
u_hat_t_RMS = zeros((gr^2), 1);

% u_plot = zeros(g, g);
% u_hat_plot = zeros(g, g);

for i=1:(gr^2)
    u_t(i,:) = ifft(Uc(i, :, 99), np);
    u_t_RMS(i) = sqrt(sum(abs((u_t(i,:)).^2))/np);
end
u_plot = reshape(u_t_RMS, gr, gr);

for i=1:(gr^2)
    u_hat_t(i, :) = ifft(Uc_hat(i, :, 99),np);
    u_hat_t_RMS(i) = sqrt(sum(abs((u_hat_t(i,:)).^2))/np);
end
u_hat_plot = reshape(u_hat_t_RMS, gr, gr);

% function definitions used in the main loop
% Resampling function
    function [pnew, wnew] = resample(p,w)
        % Resample the particles p according to weights w
        % Particles must be one ROW of p!!!
        pnew = zeros(1,N); % initializing the particles after resampling
        wnew = zeros(1,N); % initializing the weights after resampling
        % Cumulative Sum of Weights (CSW)
        CSW = zeros(1,N); % initializing the CSW vector
        CSW(1) = w(1);
        for f = 2:N
            CSW(f) = CSW(f-1) + w(f);
        end
        u(1) = rand()/N; % drawing a sample from uniform distribution (0, 1/N)
        % main loop for resampling
        for e = 1:N
            u(e) = u(1) + (e-1)/N;
            f = 1;
            while u(e)>CSW(f)
                f = f+1;
            end
            pnew(e) = p(f);
            wnew(e) = 1/N;
        end % End of resampling main loop
    end % End of resampling function
end % End of main function PF_v5