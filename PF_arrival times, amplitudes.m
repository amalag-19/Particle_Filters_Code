function PF()
%initializing N particles
M = 1000; % number of frames
N = 2000; % number of particles
x = 1; % initial state (amplitude of the virtual source)
X = x; % initialize sequence of true states
pCov = 2; % initial covariance of particle distribution
% p = repmat(x',1,N)  + sqrt(pCov)*randn(length(x),N); % initialize particles
p = x + sqrt(pCov)*randn(length(x),N); % initialize particles
xp = zeros(length(x),N); % initialize state estimate
X_hat = []; % intialize sequence of state estimates
q = zeros(1,N); % initialize particle weights

% state space equations
for m = 1:M
    % nonlinear discrete time system with additive noise
    x = f(x,v1); % true state at time t
    y = g(x,v2); % observation at time t
    % particle filtering
    for i = 1:N
        xp(i) = f(p(i),z); % particle prediction
        yp(i) = h(xp(i),tau); % prediction measurement
        %prediction performance
        d = y - yp; % difference between prediction measurement and observation
        for n = 1:Ns
            for p = 1:2
                q(i) = (1/sqrt(tauCov^Ns))*exp((-1/(2*tauCov))*(yk(n)-)); % assign importance weight to each particle
                k(i) = k(i) + q(i);
            end
            l(i) = l(i) + k(i);
        end
    end
    q = q./sum(q); % normalize the likelihood of each a priori estimate
    % The state estimate is the weighted mean of the particles
    x_hat = q*p'; % weighted sum
    X = [X x]; X_hat = [X_hat x_hat]; % update sequence of true states and state estimates
    p = re_sample(xp,q); % resampling
end

%function definitions

function x = f(x,z) %  prediction equation or state transition equation
x = x + z;
end

function y = h(x,tau) % observation equation
syms p x
y = symsum(a(p)*s(n-x(p)), p, 1, 2) + tau;
end

function n = z() % process noise
zMean = 0; % mean
zCov = 1; % covariance
n = zMean + sqrt(zCov)*randn;
end

function n = tau() % measurement noise
tauMean = 0; % mean
tauCov = 1; % covariance
n = tauMean + sqrt(tauCov)*randn;
end

function p = re_sample(xp,pdf) % resampling
cdf = cumsum(pdf); % cumulative sum
diff = cdf'*ones(1,length(pdf)) - ones(length(pdf),1)*rand(1,length(pdf));
diff = (diff <= 0) * 2 + diff;
[~, idx] = min(diff);
p = xp(idx);
end
end