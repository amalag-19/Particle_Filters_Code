% (incorrect version) Particle filter algorithm for a 2 dimensional multiple source array
% amplitude and location estimation using a 2 dimensional array of microphones and
% desired & reconstructed normal surface velocity at the actual source surface

% main function
function [xag, yag, X_hat, U, u_plot, u_hat_plot, error_amp, error_loc, error_vel, w] = PF_v8()

% parameters and initial conditions
fM = 99; % number of frames
M = 5; % number of microphones
L = 5; % number of sources
N = 100; % number of particles
N_thr = N/5; % threshold sample size
N_eff = zeros(1,(L^2)); % initialize effective sample size array
c = 343; % speed of sound in m/s
rho = 1.2754; % density of air in kg/m^3
pCov = 0.1; % initial covariance of particle distribution
sxCov = 0.0001; % source location x coordinate covariance
syCov = 0.0001; % source location y coordinate covariance
v1Cov = 0.1; % process noise covariance
v2Cov = 0.1; % measurement noise covariance
freq = 100:100:1700; % defining frequency vector
gr = 21; % grid points needed for contour plots of surface velocity
x = 100*ones((L^2), 1); % x = [1 2 3 4 5]'.*ones(L, 1) ; % initial state vector (amplitudes of the sources)
% uv = 100*ones((L^2), 1);

% true states initialization
X = zeros((L^2), fM+1); % initialize sequence of true states
X(:,1) = x; % true state vector at frame 1 is just the initial state vector
U = zeros((gr^2), length(freq), fM); % initialize sequence of true velocity vectors

% particles and estimates initialization for source amplitudes
xn = zeros((L^2),N);
for q = 1:N
    xn(:,q) = x; % array of the initial states used later for initializing particles
end
p = normrnd(xn,pCov,(L^2),N); % initializing N particles
xp = zeros((L^2),N); % initialize state estimate
x_hat = zeros((L^2),1); % initial estimated state vector (amplitudes of the sources)
X_hat = zeros((L^2), fM+1); % intialize sequence of state estimates
X_hat(:,1) = x; % estimated state vector at the frame 1 is just the initial state vector
U_hat = zeros((gr^2), length(freq), fM); % initialize sequence of estimated velocity vectors

% initialize particle weights
w = zeros(fM+1,N,(L^2));
for n=1:N
    for fm = 1:fM
        if fm == 1,
            w(fm,n,:) = (1/N)*ones((L^2),1); % Equally normalized weights for all the particles at frame 1
        else w(fm,n,:) = 0;
        end
    end
end
v = zeros((L^2),N);

% distance calculation between the source and microphones
xsv_mean = linspace(-0.4, 0.4, L);
xsv_p = zeros(N,L);
for i = 1:N
    xsv_p(i,:) = normrnd(xsv_mean,sxCov,1,L); % x coordinates of the L^2 sources
end

ysv_mean= linspace(-0.4, 0.4, L); 
ysv_p = zeros(N,L);
for i = 1:N
    ysv_p(i,:) = normrnd(ysv_mean,syCov,1,L); % y coordinates of the L^2 sources
end

[xsg_mean,ysg_mean] = meshgrid(xsv_mean,ysv_mean); % defining grid points for L^2 sources

xsg_v_p = zeros((L^2),N);
ysg_v_p = zeros((L^2),N);

xsg_v_mean = reshape(xsg_mean, L^2, 1);
ysg_v_mean = reshape(ysg_mean, L^2, 1);
for i = 1:N
    [xsg,ysg] = meshgrid(xsv_p(i,:),ysv_p(i,:)); % defining grid points for L^2 sources
    xsg_v_p(:,i) = reshape(xsg, L^2, 1);
    ysg_v_p(:,i) = reshape(ysg, L^2, 1);
end
% dels=((max(xsv)-min(xsv))*(max(ysv)-min(ysv)))/(L^2);

figure;plot(xsg_v_mean, ysg_v_mean,'x');
zsg_v = -0.02;

dx_v = zeros((L^2),N);
dy_v = zeros((L^2),N);

rs_v_mean = zeros((L^2),3);
rs_v_mean_mod = zeros((L^2), 1);
rs_v_hat = zeros((L^2),3);
rs_v_hat_mod = zeros((L^2), 1);
Rs_v_hat = zeros((L^2),3,fM);
for j=1:(L^2)
    rs_v_mean(j,:) = [xsg_v_mean(j), ysg_v_mean(j), zsg_v];
    rs_v_mean_mod(j) = (sum(rs_v_mean(j,:).^2))^(1/2);
end
rs_v_p = zeros((L^2), 3, N);
for i = 1:N
    dx_v(:,i) = xsg_v_p(:,i) - xsg_v_mean;
    dy_v(:,i) = ysg_v_p(:,i) - ysg_v_mean;
    for j=1:(L^2)
        rs_v_p(j,:,i) = [xsg_v_p(j,i), ysg_v_p(j,i), zsg_v];
    end
end

xmv = linspace(-0.4, 0.4, M); % x coordinates of the M^2 microphones
ymv = linspace(-0.4, 0.4, M); % y coordinates of the M^2 microphones
[xmg,ymg] = meshgrid(xmv,ymv); % defining grid points for M^2 microphones
hold on; plot (xmg, ymg, 's');

xmg_v = reshape(xmg, M^2, 1);
ymg_v = reshape(ymg, M^2, 1);
zmg_v = 0.35; 
rm_v = zeros((M^2), 3);
for i=1:(M^2)
    rm_v(i,:) = [xmg_v(i), ymg_v(i), zmg_v];
end

rsm_mean = zeros(M^2, L^2);
for j=1:(L^2)
    for i=1:(M^2)
        rsm_mean(i,j) = (sum((rm_v(i,:)-rs_v_mean(j,:)).^2))^(1/2);
    end
end
rsm = zeros(M^2, L^2, N);
for k = 1:N
    for j=1:(L^2)
        for i=1:(M^2)
            rsm(i,j,k) = (sum((rm_v(i,:)-rs_v_p(j,:,k)).^2))^(1/2);
        end
    end
end

% defining grid points
xav = linspace(-0.5, 0.5, gr);
yav = linspace(-0.5, 0.5, gr);
[xag,yag] = meshgrid(xav,yav); % defining grid points for velocity
xag_v = reshape(xag, gr^2, 1);
yag_v = reshape(yag, gr^2, 1);
zag_v = 0; 
ra_v = zeros((gr^2), 3);
for i=1:(gr^2)
    ra_v(i,:) = [xag_v(i), yag_v(i), zag_v];
end

rsa_mean = zeros(gr^2, L^2);
for j=1:(L^2)
    for i=1:(gr^2)
        rsa_mean(i,j) = (sum((ra_v(i,:)-rs_v_mean(j,:)).^2))^(1/2);
    end
end
rsa = zeros(gr^2, L^2, N);
for k=1:N
    for j=1:(L^2)
        for i=1:(gr^2)
            rsa(i,j,k) = (sum((ra_v(i,:)-rs_v_p(j,:,k)).^2))^(1/2);
        end
    end
end

% error matrices initialization
error_amp = zeros(length(freq), fM+1); % initializing the amplitude error matrix
error_loc = zeros(length(freq), fM); % initializing the source location error matrix
error_vel = zeros(length(freq)); % initializing the velocity error matrix

% error = zeros(length(freq), L, fM+1);

% initializing covariance matrix
Cxx = zeros((L^2),(L^2),N);
for k = 1:N
    Cxx(:,:,k) = pCov*eye(L^2); % covariance matrix of state vector vs. state vector
end

fid=1; % frquency index in error matrix

% main loop
for freq = 100:100:1700 % for different frequencies between 100 Hz to 1700 Hz
    sum_vel = 0;
    sum_vel_hat = 0;
    G_mean = (exp((-1i*(2*pi*freq*rsm_mean))/c))./rsm_mean;
    G = (exp((-1i*(2*pi*freq*rsm))/c))./rsm; % (M^2)*(L^2) Green's matrix used in measurement equation
    for fm=1:fM % for different frames upto fM
        % nonlinear discrete system in frequency domain with additive noise
        % state space equations
        x = x + normrnd(0, v1Cov, L^2, 1); % process equation: true state at frame fm
        % uv = uv + normrnd(0, v1Cov, L^2, 1);
        y = G_mean*x + normrnd(0, v2Cov, M^2, 1); % measurement equation: pressure on M microphones at frame fm corresponding to observations
        
        % true normal surface velocity calculation
        ps = zeros((gr^2),1); % initializing sound pressure on actual sound surface
%         for i=1:(L^2)
%             ps = ps + (1i*rho*freq/2)*uv(i)*dels*((exp((-1i*(2*pi*freq*rsa(:,i)))/c))./rsa(:,i));
%         end
        for i=1:(L^2)
            ps = ps + x(i)*((exp((-1i*(2*pi*freq*rsa(:,i)))/c))./rsa(:,i));
        end
        u = zeros((gr^2),1); % initializing true normal surface velocity on the actual surface
        for i=1:(L^2)
            u = u + (1/(1i*rho*freq))*(((zsg_v-zag_v)./rsa(:,i)).*((1i*2*pi*freq/c)+(1./rsa(:,i))).*ps);
        end
        U(:,fid,fm) = u; % true normal surface velocity matrix at frame fm
        
        % particle filtering
        for j = 1:N
            xp(:,j) = p(:,j) + normrnd(0, v1Cov, L^2, 1); % particle prediction
            yp = G(:,:,j)*xp(:,j) + normrnd(0, v2Cov, M^2, 1); % prediction measurement
            %prediction performance
            dx = xp(:,j) - x; % difference between particle prediction and true state vector
            % dy = yp - y; % difference between prediction measurement and observation
            % Covariance matrices
            Cxx(:,:,j) = Cxx(:,:,j) + v1Cov*eye(L^2);
            Cyx = (G(:,:,j)*Cxx(:,:,j));
            Cxy = Cyx';
            Cyy = G(:,:,j)*Cxx(:,:,j)*G(:,:,j)' + v2Cov*eye(M^2);
            temp1 = zeros((L^2),(L^2)); % temporary L*L matrix for inverse calculation later
            for g=1:(L^2)
                for h=1:(L^2)
                    temp1(g,h) = Cxx(g,h,j);
                end
            end
            mu = y + Cyx*(inv(temp1))*dx;
            C = Cyy - Cyx*(inv(temp1))*Cxy; % final covariance matrix used in likelihood function
%             t1=inv(C);
%             t2=(1/((2*pi)^(M^2)*det(C)))^(1/2);
%             t3 = abs((y-mu)'*(t1)*(y-mu));
%             t4 = exp((-1/2)*t3);
%             t5 = abs(t2*t4);
%             D=(C+C')/2;
%             lf = abs(mvnpdf(yp',y',D));
%             pd = abs(mvnpdf(xp(:,j)',x',(pCov*eye(L^2))));
            lf1 = abs(mvnpdf(yp',mu',((C+C')/2)));
            lf2 = abs((1/((((2*pi)^(M^2))*det(C))^(1/2)))*exp((-1/2)*(yp-mu)'*(inv(C))*(yp-mu))); % likelihood function
            pd = abs(((1/(((2*pi)*(sxCov))^((L^2)/2)))*exp((-1/2)*dx_v(:,j)'*(inv(sxCov*eye(L^2)))*dx_v(:,j)))*((1/(((2*pi)*(syCov))^((L^2)/2)))*exp((-1/2)*dy_v(:,j)'*(inv(syCov*eye(L^2)))*dy_v(:,j))))*(10^(-80)); % prior distribution
            w(fm+1,j,:) = w(fm,j,:)*lf*pd; % weight update equation
            % w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-((abs(d_avg))^2)/(2*v2Cov)));
            % page 324 of steven kay's book
            % *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
        end
        
        % Weight normalization and state estimation
        temp2 = zeros(3,N);
        for i=1:(L^2)
            w(fm+1,:,i) = w(fm+1,:,i)./(sum(w(fm+1,:,i))); % normalize the likelihood of each a priori estimate
            for k = 1:N
                v(i,k) = w(fm+1,k,i);
            end
            x_hat(i) = v(i,:)*p(i,:)'; % weighted sum: the state estimate is the weighted mean of the particles
            for t1 = 1:3
                for t2 = 1:N
                    temp2(t1, t2) = rs_v_p(i,t1,t2);
                end
            end
            rs_v_hat(i,:) = v(i,:)*temp2.';
        end
        
        rsa_hat = zeros(gr^2, L^2);
        for j=1:(L^2)
            for i=1:(gr^2)
                rsa_hat(i,j) = (sum((ra_v(i,:)-rs_v_hat(j,:)).^2))^(1/2);
            end
        end
        % Estimated normal surface velocity calculation
        ps_hat = zeros((gr^2),1); % initializing sound pressure on actual sound surface
        for i=1:(L^2)
            ps_hat = ps_hat + x_hat(i)*((exp((-1i*(2*pi*freq*rsa_hat(:,i)))/c))./rsa_hat(:,i));
        end
        u_hat = zeros((gr^2),1); % initializing true normal surface velocity on the actual surface
        for i=1:(L^2)
            u_hat = u_hat + (1/(1i*rho*freq))*(((zsg_v-zag_v)./rsa_hat(:,i)).*((1i*2*pi*freq/c)+(1./rsa_hat(:,i))).*ps_hat);
        end
        U_hat(:,fid,fm) = u_hat; % estimated normal surface velocity matrix at frame fm
        
        X(:,fm+1) = x; X_hat(:,fm+1) = x_hat; Rs_v_hat(:,:,fm)=rs_v_hat; % update sequence of true states and state estimates
        
        % Error matrices for source amplitude and normal surface velocity
        error_amp(fid,fm+1) = (norm(X_hat(:,fm+1)-X(:,fm+1))/norm(X(:,fm+1)))*100; % source amplitude error matrix
        for i = 1:(L^2)
            rs_v_hat_mod(i) = (sum(rs_v_hat(i,:).^2))^(1/2);
        end
        error_loc(fid,fm) = (norm(rs_v_hat_mod-rs_v_mean_mod)/norm(rs_v_mean_mod))*100; % source amplitude error matrix
       
        % Sum of true and estimated normal surface velocity matrices over all frames 
        sum_vel = sum_vel + sum(U(:,fid,fm));
        sum_vel_hat = sum_vel_hat + sum(U_hat(:,fid,fm));
        
        % Resampling
        for i=1:(L^2)
            N_eff(i) = 1/(sum((w(fm+1,:,i).^2))); % effective sample size
            if N_eff(i)<N_thr % Resampling only required when effective sample size is less than threshold sample size
                [p(i,:), v(i,:)] = resample(xp(i,:),v(i,:)); % Resampling
                for k = 1:N
                    w(fm+1,k,i) = v(i,k);
                end
            else p=xp;
            end
        end % End of resampling loop
    end % End of frame loop
    U_mean = sum_vel/fM;
    U_hat_mean = sum_vel_hat/fM;
    error_vel(fid) = ((U_hat_mean-U_mean)./U_mean).*100;
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
    u_t(i,:) = ifft(U(i, :, 99), np);
    u_t_RMS(i) = sqrt(sum(abs((u_t(i,:)).^2))/np);
end
u_plot = reshape(u_t_RMS, gr, gr);

for i=1:(gr^2)
    u_hat_t(i, :) = ifft(U_hat(i, :, 99),np);
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