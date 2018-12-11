% (final version) Particle filter algorithm for a 2 dimensional multiple source array
% amplitude and location estimation using a 2 dimensional array of microphones and
% desired & reconstructed normal surface velocity at the actual source surface

% main function
function [xag, yag, Xsa_hat, U, u_plot, u_hat_plot, error_amp, error_amp_MVU, error_loc, error_vel, w] = PF_v10()
% parameters and initial conditions
fM = 99; % number of frames
M = 5; % number of microphones
L = 5; % number of sources
N = 100; % number of particles
N_thr = N/5; % threshold sample size
N_eff = zeros(1,(L^2)); % initialize effective sample size array
c = 343; % speed of sound in m/s
rho = 1.2754; % density of air in kg/m^3
xsa_pCov = 0.1; % initial covariance of particle distribution
xsl_pCov = 0.0001; % source location x coordinate covariance
ysl_pCov = 0.0001; % source location y coordinate covariance
v1Cov = 0.1; % process noise covariance
v2Cov = 0.1; % measurement noise covariance
freq = 100:100:1700; % defining frequency vector
gr = 21; % grid points needed for contour plots of surface velocity
xsa = 100*ones((L^2), 1); % x = [1 2 3 4 5]'.*ones(L, 1) ; % initial state vector (amplitudes of the sources)
% uv = 100*ones((L^2), 1);

% true states initialization
Xsa = zeros((L^2), fM+1); % initialize sequence of true states
Xsa(:,1) = xsa; % true state vector at frame 1 is just the initial state vector
U = zeros((gr^2), length(freq), fM); % initialize sequence of true velocity vectors

% particles and estimates initialization for source amplitudes
xn = repmat(xsa, 1, N);
sa_p = normrnd(xn,xsa_pCov,(L^2),N); % initializing N source amplitude particles
xsa_p = zeros((L^2),N); % initialize state estimate
xsa_hat = zeros((L^2),1); % initial estimated state vector (amplitudes of the sources)
Xsa_hat = zeros((L^2), fM+1); % intialize sequence of state estimates
Xsa_hat(:,1) = xsa; % estimated state vector at the frame 1 is just the initial state vector
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

% defining virtual source location points mean and particles from gaussian
% distribution
xslv_mean = linspace(-0.4, 0.4, L);
yslv_mean= linspace(-0.4, 0.4, L);
[xslg_mean,yslg_mean] = meshgrid(xslv_mean,yslv_mean); % defining grid points for L^2 sources
xslgv_mean = reshape(xslg_mean, L^2, 1); % x coordinate source location grid vector mean
yslgv_mean = reshape(yslg_mean, L^2, 1); % y coordinate source location grid vector mean
zslgv = -0.02; % z coordinate source location constant for all sources and all particles
rslv_mean = zeros((L^2),3);
rslv_mean_mod = zeros((L^2), 1);
for j=1:(L^2)
    rslv_mean(j,:) = [xslgv_mean(j), yslgv_mean(j), zslgv];
    rslv_mean_mod(j) = (sum(rslv_mean(j,:).^2))^(1/2);
end

xslv_p = zeros(N,L);
yslv_p = zeros(N,L);
xslgv_p = zeros((L^2),N); % initialize x coordinate source location grid vector particles
yslgv_p = zeros((L^2),N); % initialize y coordinate source location grid vector particles
rslv_p = zeros((L^2), 3, N);
rslv_hat = zeros((L^2),3);
rslv_hat_mod = zeros((L^2), 1);
Rslv_hat = zeros((L^2),3,length(freq),fM);
for i = 1:N
    xslv_p(i,:) = normrnd(xslv_mean,xsl_pCov,1,L); % x coordinates of the L^2 sources
    yslv_p(i,:) = normrnd(yslv_mean,ysl_pCov,1,L); % y coordinates of the L^2 sources
    [xslg_p,yslg_p] = meshgrid(xslv_p(i,:),yslv_p(i,:)); % defining grid points for L^2 sources
    xslgv_p(:,i) = reshape(xslg_p, L^2, 1); % x coordinate source location grid vector particle
    yslgv_p(:,i) = reshape(yslg_p, L^2, 1); % y coordinate source location grid vector particle
    for j=1:(L^2)
        rslv_p(j,:,i) = [xslgv_p(j,i), yslgv_p(j,i), zslgv];
    end
end

% dels=((max(xsv)-min(xsv))*(max(ysv)-min(ysv)))/(L^2);
figure;plot(xslgv_mean, yslgv_mean,'x');
axis([-0.5 0.5 -0.5 0.5]);

dx_v = zeros((L^2),N);
dy_v = zeros((L^2),N);
for i = 1:N
    dx_v(:,i) = xslgv_p(:,i) - xslgv_mean;
    dy_v(:,i) = yslgv_p(:,i) - yslgv_mean;
end


% defining microphone location points
xmv = linspace(-0.4, 0.4, M); % x coordinates of the M^2 microphones
ymv = linspace(-0.4, 0.4, M); % y coordinates of the M^2 microphones
[xmg,ymg] = meshgrid(xmv,ymv); % defining grid points for M^2 microphones
xmgv = reshape(xmg, M^2, 1);
ymgv = reshape(ymg, M^2, 1);
zmgv = 0.35; 
rmv = zeros((M^2), 3);
for i=1:(M^2)
    rmv(i,:) = [xmgv(i), ymgv(i), zmgv];
end
hold on; plot (xmg, ymg, 's');

% distance calculation between the virtual sources and microphones
rsm_mean = zeros(M^2, L^2);
for j=1:(L^2)
    for i=1:(M^2)
        rsm_mean(i,j) = (sum((rmv(i,:)-rslv_mean(j,:)).^2))^(1/2);
    end
end
rsm_p = zeros(M^2, L^2, N);
for k = 1:N
    for j=1:(L^2)
        for i=1:(M^2)
            rsm_p(i,j,k) = (sum((rmv(i,:)-rslv_p(j,:,k)).^2))^(1/2);
        end
    end
end

% defining actual source points
xav = linspace(-0.5, 0.5, gr);
yav = linspace(-0.5, 0.5, gr);
[xag,yag] = meshgrid(xav,yav); % defining grid points for velocity
xagv = reshape(xag, gr^2, 1);
yagv = reshape(yag, gr^2, 1);
zagv = 0; 
rav = zeros((gr^2), 3);
for i=1:(gr^2)
    rav(i,:) = [xagv(i), yagv(i), zagv];
end

% distance calculation between the virtual source points and actual source
% points
rsa_mean = zeros(gr^2, L^2);
for j=1:(L^2)
    for i=1:(gr^2)
        rsa_mean(i,j) = (sum((rav(i,:)-rslv_mean(j,:)).^2))^(1/2);
    end
end
rsa_p = zeros(gr^2, L^2, N);
for k=1:N
    for j=1:(L^2)
        for i=1:(gr^2)
            rsa_p(i,j,k) = (sum((rav(i,:)-rslv_p(j,:,k)).^2))^(1/2);
        end
    end
end

% error matrices initialization
error_amp = zeros(length(freq), fM+1); % initializing the amplitude error matrix for PF
error_amp_MVU = zeros(length(freq), fM+1); % initializing the amplitude error matrix for MVU
error_loc = zeros(length(freq), fM); % initializing the source location error matrix
error_vel = zeros(length(freq)); % initializing the velocity error matrix

% error = zeros(length(freq), L, fM+1);

% initializing source amplitude covariance matrix
Cxxsa = zeros((L^2),(L^2),N);
for k = 1:N
    Cxxsa(:,:,k) = xsa_pCov*eye(L^2); % covariance matrix of state vector vs. state vector
end

fid=1; % frquency index in error matrix

% main loop
for freq = 100:100:1700 % for different frequencies from 100 Hz to 1700 Hz
    sum_vel = 0;
    sum_vel_hat = 0;
    G_mean = (exp((-1i*(2*pi*freq*rsm_mean))/c))./rsm_mean;
    G = (exp((-1i*(2*pi*freq*rsm_p))/c))./rsm_p; % (M^2)*(L^2) Green's matrix used in measurement equation
    for fm=1:fM % for different frames upto fM
        % nonlinear discrete system in frequency domain with additive noise
        % state space equations
        xsa = xsa + normrnd(0, v1Cov, L^2, 1); % process equation: true state at frame fm
        % uv = uv + normrnd(0, v1Cov, L^2, 1);
        ysa = G_mean*xsa + normrnd(0, v2Cov, M^2, 1); % measurement equation: pressure on M microphones at frame fm corresponding to observations
        
        % true normal surface velocity calculation
        ps = zeros((gr^2),1); % initializing sound pressure on actual sound surface
%         for i=1:(L^2)
%             ps = ps + (1i*rho*freq/2)*uv(i)*dels*((exp((-1i*(2*pi*freq*rsa(:,i)))/c))./rsa(:,i));
%         end
        for i=1:(L^2)
            ps = ps + xsa(i)*((exp((-1i*(2*pi*freq*rsa_mean(:,i)))/c))./rsa_mean(:,i));
        end
        u = zeros((gr^2),1); % initializing true normal surface velocity on the actual surface
        for i=1:(L^2)
            u = u + (1/(1i*rho*freq))*(((zslgv-zagv)./rsa_mean(:,i)).*((1i*2*pi*freq/c)+(1./rsa_mean(:,i))).*ps);
        end
        U(:,fid,fm) = u; % true normal surface velocity matrix at frame fm
        
        Cyy_MVU = v2Cov*eye(M^2); 
        x_hat_MVU = abs((inv(G_mean'*(inv(Cyy_MVU))*G_mean))*G_mean'*(inv(Cyy_MVU))*ysa);
        error_amp_MVU(fid,fm+1) = ((norm(x_hat_MVU-xsa))/norm(xsa))*100;
        
        % particle filtering
        for j = 1:N
            xsa_p(:,j) = sa_p(:,j) + normrnd(0, v1Cov, L^2, 1); % particle prediction
            ysa_p = G(:,:,j)*xsa_p(:,j) + normrnd(0, v2Cov, M^2, 1); % prediction measurement
            %prediction performance
            dxsa = xsa_p(:,j) - xsa; % difference between particle prediction and true state vector
            % dy = yp - y; % difference between prediction measurement and observation
            % Covariance matrices
            Cxxsa(:,:,j) = Cxxsa(:,:,j) + v1Cov*eye(L^2);
            Cyxsa = (G(:,:,j)*Cxxsa(:,:,j));
            Cxysa = Cyxsa';
            Cyysa = G(:,:,j)*Cxxsa(:,:,j)*G(:,:,j)' + v2Cov*eye(M^2);
            temp1 = zeros((L^2),(L^2)); % temporary L*L matrix for inverse calculation later
            for g=1:(L^2)
                for h=1:(L^2)
                    temp1(g,h) = Cxxsa(g,h,j);
                end
            end
            mu = ysa + Cyxsa*(inv(temp1))*dxsa;
            Csa = Cyysa - Cyxsa*(inv(temp1))*Cxysa; % final covariance matrix used in likelihood function
%             t1=inv(C);
%             t2=(1/((2*pi)^(M^2)*det(C)))^(1/2);
%             t3 = abs((y-mu)'*(t1)*(y-mu));
%             t4 = exp((-1/2)*t3);
%             t5 = abs(t2*t4);
%             D=(C+C')/2;
%             lf = abs(mvnpdf(yp',y',D));
%             pd = abs(mvnpdf(xp(:,j)',x',(pCov*eye(L^2))));

%             lf1 = abs(mvnpdf(ysa_p',mu',((Csa+Csa')/2))); % likelihood function
%             % (mvnpdf(xsa_p(:,j)',xsa',(xsa_pCov*eye(L^2))))*
%             pd = abs((mvnpdf(xslgv_p(:,j)',xslgv_mean',(xsl_pCov*eye(L^2))))*(mvnpdf(yslgv_p(:,j)',yslgv_mean',(ysl_pCov*eye(L^2))))); % prior distribution
% %             lf = abs((1/((((2*pi)^(M^2))*det(C))^(1/2)))*exp((-1/2)*(yp-mu)'*(inv(C))*(yp-mu))); % likelihood function
% %             pd = abs(((1/(((2*pi)*(sxCov))^((L^2)/2)))*exp((-1/2)*dx_v(:,j)'*(inv(sxCov*eye(L^2)))*dx_v(:,j)))*((1/(((2*pi)*(syCov))^((L^2)/2)))*exp((-1/2)*dy_v(:,j)'*(inv(syCov*eye(L^2)))*dy_v(:,j))))*(10^(-80)); % prior distribution
%             w(fm+1,j,:) = w(fm,j,:)*lf*pd; % weight update equation

            lf = abs((1/((((2*pi)^(M^2))*det(Csa))^(1/2)))*exp((-1/2)*(ysa_p-mu)'*(inv(Csa))*(ysa_p-mu))); % likelihood function
            %(xsa_pCov*eye(L^2))
            Cj = [xsa_pCov*eye(L^2), zeros((L^2),2*(L^2)); zeros((2*(L^2)), (L^2)), (xsl_pCov*eye(2*(L^2)))];
            pd = abs((1/((((2*pi)*(xsl_pCov))^((L^2)/2))*(((2*pi)*(ysl_pCov))^((L^2)/2))*(((2*pi)*(xsa_pCov))^((L^2)/2))))*exp((-1/2)*[dxsa', dx_v(:,j)', dy_v(:,j)']*(inv(Cj))*[dxsa; dx_v(:,j); dy_v(:,j)])); % prior distribution
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
            xsa_hat(i) = v(i,:)*sa_p(i,:)'; % weighted sum: the state estimate is the weighted mean of the particles
            for t1 = 1:3
                for t2 = 1:N
                    temp2(t1, t2) = rslv_p(i,t1,t2);
                end
            end
            rslv_hat(i,:) = v(i,:)*temp2.';
        end
        
        rsa_hat = zeros(gr^2, L^2);
        for j=1:(L^2)
            for i=1:(gr^2)
                rsa_hat(i,j) = (sum((rav(i,:)-rslv_hat(j,:)).^2))^(1/2);
            end
        end
        % Estimated normal surface velocity calculation
        ps_hat = zeros((gr^2),1); % initializing sound pressure on actual sound surface
        for i=1:(L^2)
            ps_hat = ps_hat + xsa_hat(i)*((exp((-1i*(2*pi*freq*rsa_hat(:,i)))/c))./rsa_hat(:,i));
        end
        u_hat = zeros((gr^2),1); % initializing true normal surface velocity on the actual surface
        for i=1:(L^2)
            u_hat = u_hat + (1/(1i*rho*freq))*(((zslgv-zagv)./rsa_hat(:,i)).*((1i*2*pi*freq/c)+(1./rsa_hat(:,i))).*ps_hat);
        end
        U_hat(:,fid,fm) = u_hat; % estimated normal surface velocity matrix at frame fm
        
        Xsa(:,fm+1) = xsa; Xsa_hat(:,fm+1) = xsa_hat; Rslv_hat(:,:,fid,fm)=rslv_hat; % update sequence of true states and state estimates
        
        % Error matrices for source amplitude and normal surface velocity
        error_amp(fid,fm+1) = (norm(Xsa_hat(:,fm+1)-Xsa(:,fm+1))/norm(Xsa(:,fm+1)))*100; % source amplitude error matrix
        for i = 1:(L^2)
            rslv_hat_mod(i) = (sum(rslv_hat(i,:).^2))^(1/2);
        end
        error_loc(fid,fm) = (norm(rslv_hat_mod-rslv_mean_mod)/norm(rslv_mean_mod))*100; % source amplitude error matrix
       
        % Sum of true and estimated normal surface velocity matrices over all frames 
        sum_vel = sum_vel + sum(U(:,fid,fm));
        sum_vel_hat = sum_vel_hat + sum(U_hat(:,fid,fm));
        
        % Resampling
        for i=1:(L^2)
            N_eff(i) = 1/(sum((w(fm+1,:,i).^2))); % effective sample size
            if N_eff(i)<N_thr % Resampling only required when effective sample size is less than threshold sample size
                [sa_p(i,:), v(i,:)] = resample(xsa_p(i,:),v(i,:)); % Resampling
                for k = 1:N
                    w(fm+1,k,i) = v(i,k);
                end
            else sa_p=xsa_p;
            end
        end % End of resampling loop
    end % End of frame loop
    U_mean = sum_vel/fM;
    U_hat_mean = sum_vel_hat/fM;
    error_vel(fid) = ((U_hat_mean-U_mean)./U_mean).*100;
    fid=fid+1; % frequency index increment
end % end of main (frequency) loop

figure;plot(Rslv_hat(:,1,10,fM), Rslv_hat(:,2,10,fM),'x');
axis([-0.5 0.5 -0.5 0.5]);
hold on; plot (xmg, ymg, 's');

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
end % End of main function PF_v10