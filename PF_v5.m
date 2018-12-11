% Particle filter algorithm for the multiple sources amplitude estimation using a linear array of microphones

% main function
function [xg, yg, X_hat, U, U_hat, error_amp, error_vel, w] = PF_v5()

% parameters and initial conditions
fM = 99; % number of frames
M = 11; % number of microphones
L = 5; % number of sources
N = 10; % number of particles
N_thr = N/5; % threshold sample size
N_eff = zeros(1,L); % initialize effective sample size array
c = 343; % speed of sound in m/s
rho = 1.2754; % density of air in kg/m^3
pCov = 0.1; % initial covariance of particle distribution
v1Cov = 0.1; % process noise covariance
v2Cov = 0.1; % measurement noise covariance
freq = 100:100:1700; % defining frequency vector
g = 21; % grid points needed for contour plots of surface velocity
x = 100*ones(L, 1); % x = [1 2 3 4 5]'.*ones(L, 1) ; % initial state vector (amplitudes of the sources)

% true states initialization
X = zeros(L, fM+1); % initialize sequence of true states
X(:,1) = x; % true state vector at frame 1 is just the initial state vector
U = zeros(g, g, fM); % initialize sequence of true velocity vectors
sum_vel = 0;

% particles and estimates initialization
xn = zeros(L,N);
for q = 1:N
    xn(:,q) = x; % array of the initial states used later for initializing particles
end
p = normrnd(xn,pCov,L,N); % initializing N particles
xp = zeros(L,N); % initialize state estimate
x_hat = zeros(L,1); % initial estimated state vector (amplitudes of the sources)
X_hat = zeros(L, fM+1); % intialize sequence of state estimates
X_hat(:,1) = x; % estimated state vector at the frame 1 is just the initial state vector
U_hat = zeros(g, g, fM); % initialize sequence of estimated velocity vectors
sum_vel_hat = 0;

% initialize particle weights
w = zeros(L,fM+1,N);
for l = 1:L
    for k=1:N
        for fm = 1:fM
            if fm == 1,
                w(l,fm,k) = 1/N; % Equally normalized weights for all the particles at frame 1
            else w(l,fm,k) = 0;
            end
        end
    end
end
v = zeros(L,N);

% distance calculation between the sources and microphones
rxs = linspace(1, 1.4, L); % x coordinates of the L sources
rys = 1; % y coordinates of the L sources
rxm = linspace(-0.5, 0.5, M); % x coordinates of the M microphones
rym = 0; % y coordinates of the M microphones
r = zeros(M, L);
for s = 1:L
    for m = 1:M
        r(m,s) = ((((rxs(s)-rxm(m))^2)+((rys-rym)^2))^(1/2))';
    end
end

% defining grid points
xv = linspace(0.7, 1.7, g);
yv = linspace(-0.5, 2.5, g);
[xg,yg] = meshgrid(xv,yv); % defining grid points for velocity

% defining distance matrix: distances of grid points from the sources
rxyg = zeros(g, g, L);
for l = 1:L
    rxyg(:,:,l) = (((xg-rxs(l)).^2)+((yg-rys).^2)).^(1/2);
end

% error matrices initialization
error_amp = zeros(length(freq), fM+1); % initializing the amplitude error matrix
error_vel = zeros(g, g, length(freq), fM+1); % initializing the velocity error matrix
% error = zeros(length(freq), L, fM+1);

% initializing covariance matrix
Cxx = zeros(L,L,N);
for k = 1:N
    Cxx(:,:,k) = pCov*eye(L); % covariance matrix of state vector vs. state vector
end

fid=1; % frquency index in error matrix

% main loop
for freq = 100:100:1700 % for different frequencies between 100 Hz to 1700 Hz
    G = (exp((-1i*(2*pi*freq*r))/c))./r; % M*L Green's matrix used in measurement equation
    for fm=1:fM % for different frames upto fM
        % nonlinear discrete system in frequency domain with additive noise
        % state space equations
        x = x + normrnd(0, v1Cov); % process equation: true state at frame fm
        y = G*x + normrnd(0, v2Cov); % measurement equation: pressure on M microphones at frame fm corresponding to observations
        
        % true normal surface velocity calculation
        ps = 0; % initializing sound pressure on actual sound surface
        for l=1:L
            ps = ps + x(l)*((exp((-1i*(2*pi*freq*rxyg(:,:,l)))/c))./rxyg(:,:,l));
        end
        u = 0; % initializing true normal surface velocity on the actual surface
        for l=1:L
            u = u + (1/(1i*rho*freq))*(((yg-rys)./rxyg(:,:,l)).*((1i*2*pi*freq/c)+(1./rxyg(:,:,l))).*ps);
        end
        U(:,:,fm) = u; % true normal surface velocity matrix at frame fm
        
        % particle filtering
        for j = 1:N
            xp(:,j) = p(:,j) + normrnd(0, v1Cov); % particle prediction
            yp = G*xp(:,j) + normrnd(0, v2Cov); % prediction measurement
            %prediction performance
            dx = xp(:,j) - x; % difference between particle prediction and true state vector
            dy = yp - y; % difference between prediction measurement and observation
            % Covariance matrices
            Cxx(:,:,j) = Cxx(:,:,j) + v1Cov*eye(L);
            Cyx = (G*Cxx(:,:,j));
            Cxy = Cyx';
            Cyy = G*Cxx(:,:,j)*G' + v2Cov*eye(M);
            temp = zeros(L,L); % temporary L*L matrix for inverse calculation later
            for g=1:L
                for h=1:L
                    temp(g,h) = Cxx(g,h,j);
                end
            end
            C = Cyy - Cyx*(inv(temp))*Cxy; % final covariance matrix used in likelihood function
            for l=1:L
                lf = abs((1/((((2*pi)^M)*det(C))^(1/2)))*exp((-1/2)*dy'*(inv(C))*dy)); % likelihood function
                pd = abs((1/(((2*pi)*(pCov^L))^(1/2)))*exp((-1/2)*dx'*(inv(pCov*eye(L)))*dx)); % prior distribution
                w(l,fm+1,j) = w(l,fm,j)*lf*pd; % weight update equation
            end
            % w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-((abs(d_avg))^2)/(2*v2Cov)));
            % page 324 of steven kay's book
            % *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
        end
        
        % Weight normalization and state estimation
        for l=1:L
            w(l,fm+1,:) = w(l,fm+1,:)./(sum(w(l,fm+1,:))); % normalize the likelihood of each a priori estimate
            for k = 1:N
                v(l,k) = w(l,fm+1,k);
            end
            x_hat(l) = v(l,:)*p(l,:)'; % weighted sum: the state estimate is the weighted mean of the particles
        end
        
        % Estimated normal surface velocity calculation
        ps_hat = 0; % initializing sound pressure estimate on actual sound surface
        for l=1:L
            ps_hat = ps_hat + x_hat(l)*((exp((-1i*(2*pi*freq*rxyg(:,:,l)))/c))./rxyg(:,:,l));
        end
        u_hat = 0; % initializing estimated normal surface velocity on the actual surface
        for l=1:L
            u_hat = u_hat + (1/(1i*rho*freq))*(((yg-rys)./rxyg(:,:,l)).*((1i*2*pi*freq/c)+(1./rxyg(:,:,l))).*ps_hat);
        end
        U_hat(:,:,fm) = u_hat; % estimated normal surface velocity matrix at frame fm
        
        X(:,fm+1) = x; X_hat(:,fm+1) = x_hat; % update sequence of true states and state estimates
        
        % Error matrices for source amplitude and normal surface velocity
        error_amp(fid,fm+1) = (norm(X_hat(:,fm+1)-X(:,fm+1))/norm(X(:,fm+1)))*100; % source amplitude error matrix
       
        % Sum of true and estimated normal surface velocity matrices over all frames 
        sum_vel = sum_vel + sum(U(:,:,fm));
        sum_vel_hat = sum_vel_hat + sum(U_hat(:,:,fm));
        
        % Resampling
        for l=1:L
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
