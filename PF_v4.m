% Particle filter algorithm for the single source amplitude estimation 
% using a 2 dimensional array of microphones
function [error_amp_MLE, error_amp_PF, error_vel_MLE, error_vel_PF] = PF_v4()

% parameters and initial conditions
fM = 99; % number of frames
M = 3; % number of microphones
N = 1000; % number of particles
N_thr = N/5; % threshold sample size
c = 3*(10^8); % speed of light in m/s
rho = 1.2754; % density of air in kg/m^3
pCov = 0.1; % initial covariance of particle distribution
freq = 100:100:1700; % defining frequency vector
x = 100; % initial state (amplitude of the source)

% true states initialization
X = zeros(1, fM+1); % initialize sequence of true states
X(1) = x; % true state vector at frame 1 is just the initial state
U = zeros((M^2), fM); % initialize sequence of true velocity vectors

% particles and estimates initialization
p = normrnd(x, pCov, [1 N]); % initializing N particles
xp = zeros(length(x),N); % initialize state estimate
X_hat_PF = zeros(1, fM+1); % intialize sequence of state estimates
X_hat_PF(1) = x; % estimated state at the frame 1 is just the initial state
U_hat_MLE = zeros((M^2), fM); % initialize sequence of estimated velocity vectors using MLE
U_hat_PF = zeros((M^2), fM); % initialize sequence of estimated velocity vectors using PF

% initialize particle weights
w = zeros(fM+1,N);
for k=1:N
    for fm = 1:fM
        if fm == 1,
            w(fm, k) = 1/N;
        else w(fm, k) = 0;
        end
    end
end

% distance calculation between the source and microphones
rxs = 1; % x coordinate of the source
rys = 1; % y coordinate of the source
rzs = 1; % z coordinate of the source
rxym = zeros(M, M, 2);
r = zeros(1, (M^2));
% co-ordiantes of the microphones in a 2-dimensional array
id = 1; % index for distance of source from microphones
for i = 1:M % x coordinates of the M microphones
    for j = 1:M % y coordinates of the M microphones
        for k = 1:2
            if k == 1
                rxym(i,j,k) = (i/2)-1;      % (i/(M-1))-((11-M)/8);
            else if k == 2
                    rxym(i,j,k) = (j/2)-1;
                end
            end
        end
        r(id) = ((((rxs-rxym(i,j,1))^2)+((rys-rxym(i,j,2))^2)+((rzs-0)^2))^(1/2)); % distance vector between the source and M^2 microphones
        id = id+1;
    end
end

% error matrices initialization
error_amp_MLE = zeros(length(freq), fM+1); % initializing the amplitude error matrix for MLE
error_amp_PF = zeros(length(freq), fM+1); % initializing the amplitude error matrix for PF
error_vel_MLE = zeros(1, length(freq)); % initializing the velocity error matrix for MLE
error_vel_PF = zeros(1, length(freq)); % initializing the velocity error matrix for PF

Cxx = pCov*ones(1,N); % covariance matrix of state vs. state

fid=1; % frquency index in error matrix

% main loop: state space equations
for freq = 100 : 100: 1700
    G = (exp((-1i*(2*pi*freq*r'))/c))./r';
    for fm=1:fM
        % nonlinear discrete system in frequency domain with additive noise
        x = x + v1(); % true state at frame fm
        y = G*x + v2();
        
        u = (1/(1i*rho*freq))*((rzs./r').*((1i*2*pi*freq/c)+(1./r')).*y);
        U(:, fm) = u;
        
        Cyy_MLE = v2Cov*eye(M^2); 
        x_hat_MLE = abs((inv(G'*(inv(Cyy_MLE))*G))*G'*(inv(Cyy_MLE))*y);
        y_hat_MLE = G*x_hat_MLE + v2();
        u_hat_MLE = (1/(1i*rho*freq))*((rys./r').*((1i*2*pi*freq/c)+(1./r')).*y_hat_MLE);
        U_hat_MLE(:, fm+1) = u_hat_MLE; 
        error_amp_MLE(fid,fm+1) = ((x_hat_MLE-x)/x)*100;
        
        % particle filtering
        for j = 1:N
            xp(j) = p(j) + v1(); % particle prediction
            yp = G*xp(j) + v2(); % prediction measurement
            %prediction performance
            d = yp - y; % difference between prediction measurement and observation
            % Covariance matrices
            Cxx(j) = Cxx(j) + v1Cov;
            Cyx = (G*Cxx(j));
            Cxy = Cyx';
            Cyy = G*Cxx(j)*G' + v2Cov*eye(M^2);
            C = Cyy - Cyx*(inv(Cxx(j)))*Cxy;
            w(fm+1,j) = w(fm,j)*abs((1/((((2*pi)^(M^2))*det(C))^(1/2)))*exp((-1/2)*d'*(inv(C))*d));
            % w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-((abs(d_avg))^2)/(2*v2Cov)));
            % page 324 of steven kay's book
            % *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
        end
        w(fm+1,:) = w(fm+1,:)./(sum(w(fm+1,:))); % normalize the likelihood of each a priori estimate
        % The state estimate is the weighted mean of the particles
        x_hat_PF = w(fm+1,:)*p'; % weighted sum
        y_hat_PF =  G*x_hat_PF + v2();
        u_hat_PF = (1/(1i*rho*freq))*((rzs./r').*((1i*2*pi*freq/c)+(1./r')).*y_hat_PF);
        U_hat_PF(:, fm+1) = u_hat_PF;
        X(fm+1) = x; X_hat_PF(fm+1) = x_hat_PF; % update sequence of true states and state estimates
        error_amp_PF(fid,fm+1) = ((X_hat_PF(fm+1)-X(fm+1))/X(fm+1))*100;
        
        % resampling
        N_eff = 1/(sum((w(fm+1,:).^2))); % effective sample size
        if N_eff<N_thr
            [p, w(fm+1,:)] = resample(xp,w(fm+1,:));
        else p=xp;
        end   
        % w(j) = w(j)*(1/sqrt(2*pi*v2Cov))*exp(-(d(j)^2)/(2*v2Cov))*(1/sqrt(2*pi*v1Cov))*exp(-(d(m)^2)/(2*v1Cov))
    end
    
    sum_vel = 0;
    sum_vel_MLE = 0;
    sum_vel_PF = 0;
    for fm = 1:fM
        sum_vel = sum_vel + sum(U(:,fm));
        sum_vel_MLE = sum_vel_MLE + sum(U_hat_MLE(:,fm)); 
        sum_vel_PF = sum_vel_PF + sum(U_hat_PF(:,fm));
    end
    U_mean = sum_vel/fM;
    U_MLE_mean = sum_vel_MLE/fM;
    U_PF_mean = sum_vel_PF/fM;
    error_vel_MLE(fid) = (norm(U_MLE_mean-U_mean)/norm(U_mean))*100;
    error_vel_PF(fid) = (norm(U_PF_mean-U_mean)/norm(U_mean))*100;
    fid = fid + 1;
end

    function n = v1() % process noise
        v1Mean = 0; % mean
        v1Cov = 0.1; % covariance
        n = v1Mean + sqrt(v1Cov)*randn;
    end

    function n = v2() % measurement noise
        v2Mean = 0; % mean
        v2Cov = 0.1; % covariance
        n = v2Mean + sqrt(v2Cov)*randn;
    end

    function [pnew, wnew] = resample(p,w)
        % Resample the particles p according to weights w
        % Particles must be one COLUMN of p!!!
        pnew = zeros(1,N);
        wnew = zeros(1,N);
        CSW = zeros(1,N);
        CSW(1) = w(1);
        for f = 2:N
            CSW(f) = CSW(f-1) + w(f);
        end
        u(1) = rand()/N;
        
        for e = 1:N
            u(e) = u(1) + (e-1)/N;
            f = 1;
            while u(e)>CSW(f)
                f = f+1;
            end
            pnew(e) = p(f);
            wnew(e) = 1/N;
        end
    end
end
