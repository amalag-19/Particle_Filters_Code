function [error] = PF()
% initializing N particles
M = 101; % number of frames
N = 100; % number of particles
% L = 100; % total number of virtual sources
c = 3*(10^8); % speed of light in m/s
% d = 0.01;
x = 1; % initial state (amplitude of the virtual source)
X = zeros(1, M); % initialize sequence of true states
pCov = 0.001; % initial covariance of particle distribution
% p = repmat(x',1,N)  + sqrt(pCov)*randn(length(x),N); % initialize particles
p = x + sqrt(pCov)*randn(length(x),N); % initialize particles
xp = zeros(length(x),N); % initialize state estimate
X_hat = zeros(1, M); % intialize sequence of state estimates
q = zeros(1,N); % initialize particle weights
error = zeros(17, M);

% state space equations
rxs = 1; % x coordinate of the source
rys = 1; % y coordinate of the source
rx = linspace(-0.5, 0.5, M); % x coordinates of the M microphones
id=1; % error index
for freq = 100 : 100: 1700
    for m = 1:M
        ry = 0; % y coordinates of the M microphones all being equal to zero
        r = (((rxs-rx(m))^2)+((rys-ry)^2))^(1/2); %distance of microphone m from source
        G = ((exp(-1i*(2*pi*freq*r)/c)))/(r); % free space green's functions
        
        % nonlinear discrete system in frequency domain with additive noise
        x = f(x); % true state at frame m
        y = g(x,G); % observation at frame m
        
        % particle filtering
        for j = 1:N
            xp(j) = f(p(j)); % particle prediction
            yp = g(xp(j),G); % prediction measurement
            %prediction performance
            d = y - yp; % difference between prediction measurement and observation
            q(j) = (1/sqrt(2*pi*v2Cov))*exp(-(d^2)/(2*v2Cov)); % assign importance weight to each particle
            % wt(j) = 
        end
        q = q./sum(q); % normalize the likelihood of each a priori estimate
        
        % The state estimate is the weighted mean of the particles
        x_hat = q*p'; % weighted sum
        X(m) = x; X_hat(m) = x_hat; % update sequence of true states and state estimates
        error(id,m) = ((X_hat(m)-X(m))/X(m))*100;
        p = re_sample(xp,q); % resampling
    end
    id=id+1;
end

%function definitions

    function x = f(x) %  prediction equation or state transition equation
        x = x + v1();
    end

    function y = g(x,G) % observation equation
        % syms p x
        % y = symsum(a(p)*s(n-x(p)), p, 1, 2) + tau;
        y = G*x + v2();
    end

    function n = v1() % process noise
        v1Mean = 0; % mean
        v1Cov = 0.001; % covariance
        n = v1Mean + sqrt(v1Cov)*randn;
    end

    function n = v2() % measurement noise
        v2Mean = 0; % mean
        v2Cov = 0.001; % covariance
        n = v2Mean + sqrt(v2Cov)*randn;
    end

    function p = re_sample(xp,pdf) % resampling
        cdf = cumsum(pdf); % cumulative sum
        diff = cdf'*ones(1,length(pdf)) - ones(length(pdf),1)*rand(1,length(pdf));
        diff = (diff <= 0) * 2 + diff;
        [~, idx] = min(diff);
        p = xp(idx);
    end
end