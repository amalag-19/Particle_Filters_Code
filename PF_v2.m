function [error] = PF_v2()
% initializing N particles
fM = 11; % number of frames
M = 11; % number of microphones
N = 100000; % number of particles
% L = 100; % total number of virtual sources
c = 3*(10^8); % speed of light in m/s
% d = 0.1;
x = 1; % initial state (amplitude of the virtual source)
X = zeros(1, fM); % initialize sequence of true states
pCov = 0.001; % initial covariance of particle distribution
p = x + sqrt(pCov)*randn(length(x),N); % initialize particles
xp = p;
% xp = zeros(length(x),N); % initialize state estimate
X_hat = zeros(1, fM); % intialize sequence of state estimates

% initialize particle weights
w = zeros(fM,N);
for k=1:N
    for fm = 1:fM
        if fm == 1,
            w(fm, k) = 1/N;
        else w(fm, k) = 0;
        end
    end
end

error = zeros(17, fM);

% state space equations
rxs = 1; % x coordinate of the source
rys = 1; % y coordinate of the source
rx = linspace(-0.5, 0.5, M); % x coordinates of the M microphones
ry = 0;
r = (((rxs-rx).^2)+((rys-ry)^2)).^(1/2);
id=1; % error index

for freq = 100 : 100: 1700
    G = ((exp(-1i*(2*pi*freq.*r)/c)))./r;
    for fm=1:fM
        % nonlinear discrete system in frequency domain with additive noise
        x = x + v1(); % true state at frame fm
        y = G.*x + v2();
        % particle filtering
        for j = 1:N
            xp(j) = xp(j) + v1(); % particle prediction
            yp = G.*xp(j) + v2(); % prediction measurement
            %prediction performance
            d = y - yp; % difference between prediction measurement and observation
            d_avg = sum(d)/length(d);
            v2Cov = 0.001; % covariance
            w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-(d_avg^2)/(2*v2Cov)));
            % *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
        end
        w(fm+1,:) = w(fm+1,:)./(sum(w(fm+1,:))); % normalize the likelihood of each a priori estimate
        % The state estimate is the weighted mean of the particles
        x_hat = w(fm+1,:)*p'; % weighted sum
        X(fm+1) = x; X_hat(fm+1) = x_hat; % update sequence of true states and state estimates
        error(id,fm+1) = ((X_hat(fm+1)-X(fm+1))/X(fm+1))*100;
        p = re_sample(xp,w(fm+1,:)); % resampling
        % w(j) = w(j)*(1/sqrt(2*pi*v2Cov))*exp(-(d(j)^2)/(2*v2Cov))*(1/sqrt(2*pi*v1Cov))*exp(-(d(m)^2)/(2*v1Cov))
    end
    id=id+1;
end

% function [error] = PF_v2()
% % initializing N particles
% fM = 11; % number of frames
% M = 11; % number of microphones
% N = 100000; % number of particles
% % L = 100; % total number of virtual sources
% c = 3*(10^8); % speed of light in m/s
% % d = 0.1;
% x = 1; % initial state (amplitude of the virtual source)
% X = zeros(1, fM); % initialize sequence of true states
% pCov = 0.001; % initial covariance of particle distribution
% % p = repmat(x',1,N)  + sqrt(pCov)*randn(length(x),N); % initialize particles
% p = x + sqrt(pCov)*randn(length(x),N); % initialize particles
% xp = p;
% xp = zeros(length(x),N); % initialize state estimate
% X_hat = zeros(1, fM); % intialize sequence of state estimates
% 
% % initialize particle weights
% w = zeros(fM,N);
% for k=1:N
%     for fm = 1:fM
%         if fm == 1,
%             w(fm, k) = 1/N;
%         else w(fm, k) = 0;
%         end
%     end
% end
% 
% % r = zeros(1,M);
% % G = zeros(1,M);
% % y = zeros(1,M);
% error = zeros(17, fM);
% 
% % state space equations
% rxs = 1; % x coordinate of the source
% rys = 1; % y coordinate of the source
% rx = linspace(-0.5, 0.5, M); % x coordinates of the M microphones
% ry = 0;
% r = (((rxs-rx).^2)+((rys-ry)^2)).^(1/2);
% id=1; % error index
% 
% for freq = 100 : 100: 1700
%     G = ((exp(-1i*(2*pi*freq.*r)/c)))./r;
%     for fm=1:fM
%         % nonlinear discrete system in frequency domain with additive noise
%         x = x + v1(); % true state at frame fm
%         y = G.*x + v2();
%         %         for m = 1:M
%         %             ry = 0; % y coordinates of the M microphones all being equal to zero
%         %             r(m) = (((rxs-rx(m))^2)+((rys-ry)^2)).^(1/2); %distance of microphone m from source
%         %             G(m) = ((exp(-1i*(2*pi*freq*r(m))/c)))/(r(m)); % free space green's functions
%         %             y(m) = G(m)*x+v2(); % observation at frame fm
%         % particle filtering
%         for j = 1:N
%             xp(j) = xp(j) + v1(); % particle prediction
%             yp = G.*xp(j) + v2(); % prediction measurement
%             %prediction performance
%             d = y - yp; % difference between prediction measurement and observation
%             % w(j) = w(j)*(1/sqrt(2*pi*v2Cov))*exp(-(d(j)^2)/(2*v2Cov))*(1/sqrt(2*pi*v1Cov))*exp(-(d(m)^2)/(2*v1Cov))
%             % w(m,j) = (1/sqrt(2*pi*v2Cov))*exp(-(d^2)/(2*v2Cov)); % assign importance weight to each particle
%             % v1Cov = 0.001;
%             % a = (xp(j)-(3*v1Cov));
%             % b = (xp(j)+(3*v1Cov));
%             d_avg = sum(d)/length(d);
%             v2Cov = 0.001; % covariance
%             w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-(d_avg^2)/(2*v2Cov)));
%             % *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
%         end
%         %         end
%         w(fm+1,:) = w(fm+1,:)./(sum(w(fm+1,:))); % normalize the likelihood of each a priori estimate
%         % The state estimate is the weighted mean of the particles
%         x_hat = w(fm+1,:)*p'; % weighted sum
%         X(fm+1) = x; X_hat(fm+1) = x_hat; % update sequence of true states and state estimates
%         error(id,fm+1) = ((X_hat(fm+1)-X(fm+1))/X(fm+1))*100;
%         % p = re_sample(xp,w(fm+1,:)); % resampling
%         % w(j) = w(j)*(1/sqrt(2*pi*v2Cov))*exp(-(d(j)^2)/(2*v2Cov))*(1/sqrt(2*pi*v1Cov))*exp(-(d(m)^2)/(2*v1Cov))
%     end
%     id=id+1;
% end
% likelihood

% transition

% proposal
% for fm=1:fM
%     for freq = 100 : 100: 1700
%         for m = 1:M
%             ry = 0; % y coordinates of the M microphones all being equal to zero
%             r = (((rxs-rx(m))^2)+((rys-ry)^2))^(1/2); %distance of microphone m from source
%             G = ((exp(-1i*(2*pi*freq*r)/c)))/(r); % free space green's functions
%             
%             % nonlinear discrete system in frequency domain with additive noise
%             x = f(x); % true state at frame fm
%             y = g(x,G); % observation at frame fm
%             
%             % particle filtering
%             for j = 1:N
%                 xp(j) = f(p(j)); % particle prediction
%                 yp = g(xp(j),G); % prediction measurement
%                 %prediction performance
%                 d = y - yp; % difference between prediction measurement and observation
%                 q(j) = (1/sqrt(2*pi*v2Cov))*exp(-(d^2)/(2*v2Cov)); % assign importance weight to each particle
%                 % wt(j) =
%             end
%             q = q./sum(q); % normalize the likelihood of each a priori estimate
%             
%             % The state estimate is the weighted mean of the particles
%             x_hat = q*p'; % weighted sum
%             X(fm) = x; X_hat(fm) = x_hat; % update sequence of true states and state estimates
%             error(id,fm) = ((X_hat(fm)-X(fm))/X(fm))*100;
%             p = re_sample(xp,q); % resampling
%         end
%         id=id+1;
%     end
% end

%function definitions

%     function x = f(x) %  prediction equation or state transition equation
%         x = x + v1();
%     end

%     function y = g(x,G) % observation equation
%         % syms p x
%         % y = symsum(a(p)*s(n-x(p)), p, 1, 2) + tau;
%         y = G*x + v2();
%     end

%     function z = l(d) % likelihood function
%         v2Cov = 0.001; % covariance
%         z = (1/sqrt(2*pi*v2Cov))*exp(-(d^2)/(2*v2Cov));
%     end

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