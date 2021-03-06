% (correct version) Particle filter algorithm for a 2 dimensional multiple source array
% amplitude and location estimation using a 2 dimensional array of microphones and
% desired & reconstructed normal surface velocity at the actual source surface

% main function
function [xag, yag, X_hat, U, u_plot, u_hat_plot, error_ampm, error_loc, error_vel, w] = PF_v9()

% parameters and initial conditions
fM = 99; % number of frames
M = 5; % number of microphones
L = 5; % number of sources
Nsl = 20; % number of particles for source location
Nsa = 20; % number of particles for source amplitude
N_thr = Nsa/5; % threshold sample size
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
x = zeros((L^2), Nsl);
for i=1:Nsl
    x(:,i) = 100*ones((L^2), 1); % x = [1 2 3 4 5]'.*ones(L, 1) ; % initial state vector (amplitudes of the sources)
end
% uv = 100*ones((L^2), 1);

% true states initialization
X = zeros((L^2), Nsl, fM+1); % initialize sequence of true states
X(:,:,1) = x; % true state vector at frame 1 is just the initial state vector
U = zeros((gr^2), length(freq), fM, Nsl); % initialize sequence of true velocity vectors
Um = zeros((gr^2), length(freq), fM);

% particles and estimates initialization for source amplitudes
xn = zeros((L^2),Nsl,Nsa);
for r = 1:Nsl
    for q = 1:Nsa
        xn(:,r,q) = x(:,r); % array of the initial states used later for initializing particles
    end
end
p = normrnd(xn,pCov,(L^2),Nsl,Nsa); % initializing N particles
xp = zeros((L^2),Nsl,Nsa); % initialize state estimate
x_hat = zeros((L^2), Nsl); % initial estimated state vector (amplitudes of the sources)
X_hat = zeros((L^2), fM+1, Nsl); % intialize sequence of state estimates
X_hatm = zeros((L^2), fM+1);
for i=1:Nsl
    X_hat(:,1,i) = x(:,i); % estimated state vector at the frame 1 is just the initial state vector
end
U_hat = zeros((gr^2), length(freq), fM, Nsl); % initialize sequence of estimated velocity vectors

% initialize particle weights
w = zeros(fM+1,Nsl,Nsa,(L^2));
for m=1:Nsl
    for n=1:Nsa
        for fm = 1:fM
            if fm == 1,
                w(fm,m,n,:) = (1/Nsa)*ones((L^2),1); % Equally normalized weights for all the particles at frame 1
            else w(fm,n,:) = 0;
            end
        end
    end
end
v = zeros((L^2),Nsa);
wm = zeros((L^2), Nsl);

% distance calculation between the source and microphones
xsv_mean = linspace(-0.4, 0.4, L);
xsv_p = zeros(Nsl,L);
for i = 1:Nsl
    xsv_p(i,:) = normrnd(xsv_mean,sxCov,1,L); % x coordinates of the L^2 sources
end

ysv_mean= linspace(-0.4, 0.4, L); 
ysv_p = zeros(Nsl,L);
for i = 1:Nsl
    ysv_p(i,:) = normrnd(ysv_mean,syCov,1,L); % y coordinates of the L^2 sources
end

[xsg_mean,ysg_mean] = meshgrid(xsv_mean,ysv_mean); % defining grid points for L^2 sources

xsg_v_p = zeros((L^2),Nsl);
ysg_v_p = zeros((L^2),Nsl);

xsg_v_mean = reshape(xsg_mean, L^2, 1);
ysg_v_mean = reshape(ysg_mean, L^2, 1);
for i = 1:Nsl
    [xsg,ysg] = meshgrid(xsv_p(i,:),ysv_p(i,:)); % defining grid points for L^2 sources
    xsg_v_p(:,i) = reshape(xsg, L^2, 1);
    ysg_v_p(:,i) = reshape(ysg, L^2, 1);
end
% dels=((max(xsv)-min(xsv))*(max(ysv)-min(ysv)))/(L^2);

figure;plot(xsg_v_mean, ysg_v_mean,'x');
zsg_v = -0.02;

dx_v = zeros((L^2),Nsl);
dy_v = zeros((L^2),Nsl);

rs_v_mean = zeros((L^2),3);
rs_v_mean_mod = zeros((L^2), 1);
rs_v_hat = zeros((L^2),3);
rs_v_hat_mod = zeros((L^2), 1);
Rs_v_hat = zeros((L^2),3,fM);
for j=1:(L^2)
    rs_v_mean(j,:) = [xsg_v_mean(j), ysg_v_mean(j), zsg_v];
    rs_v_mean_mod(j) = (sum(rs_v_mean(j,:).^2))^(1/2);
end
rs_v_p = zeros((L^2), 3, Nsl);
for i = 1:Nsl
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

% rsm_mean = zeros(M^2, L^2);
% for j=1:(L^2)
%     for i=1:(M^2)
%         rsm_mean(i,j) = (sum((rm_v(i,:)-rs_v_mean(j,:)).^2))^(1/2);
%     end
% end
rsm = zeros(M^2, L^2, Nsl);
for k = 1:Nsl
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

% rsa_mean = zeros(gr^2, L^2);
% for j=1:(L^2)
%     for i=1:(gr^2)
%         rsa_mean(i,j) = (sum((ra_v(i,:)-rs_v_mean(j,:)).^2))^(1/2);
%     end
% end
rsa = zeros(gr^2, L^2, Nsl);
for k=1:Nsl
    for j=1:(L^2)
        for i=1:(gr^2)
            rsa(i,j,k) = (sum((ra_v(i,:)-rs_v_p(j,:,k)).^2))^(1/2);
        end
    end
end

% error matrices initialization
error_amp = zeros(length(freq), fM+1, Nsl); % initializing the amplitude error matrix
error_ampm = zeros(length(freq), fM+1);
error_loc = zeros(length(freq), fM); % initializing the source location error matrix
error_vel = zeros(length(freq)); % initializing the velocity error matrix

% error = zeros(length(freq), L, fM+1);

% initializing covariance matrix
Cxx = zeros((L^2),(L^2),Nsa);
for k = 1:Nsa
    Cxx(:,:,k) = pCov*eye(L^2); % covariance matrix of state vector vs. state vector
end

fid=1; % frquency index in error matrix

% main loop
for freq = 100:100:1700 % for different frequencies between 100 Hz to 1700 Hz
    sum_vel = 0;
    sum_vel_hat = 0;
    G = (exp((-1i*(2*pi*freq*rsm))/c))./rsm; % (M^2)*(L^2) Green's matrix used in measurement equation
    for fm=1:fM % for different frames upto fM
        % nonlinear discrete system in frequency domain with additive noise
        % state space equations
        for m = 1:Nsl % particle filtering for source location
            x(:,m) = x(:,m) + normrnd(0, v1Cov, L^2, 1); % process equation: true state at frame fm
            % uv = uv + normrnd(0, v1Cov, L^2, 1);
            y = G(:,:,m)*x(:,m) + normrnd(0, v2Cov, M^2, 1); % measurement equation: pressure on M microphones at frame fm corresponding to observations
        
            % true normal surface velocity calculation
            ps = zeros((gr^2),1); % initializing sound pressure on actual sound surface
%             for i=1:(L^2)
%                 ps = ps + (1i*rho*freq/2)*uv(i)*dels*((exp((-1i*(2*pi*freq*rsa(:,i)))/c))./rsa(:,i));
%             end
            for i=1:(L^2)
                ps = ps + x(i)*((exp((-1i*(2*pi*freq*rsa(:,i,m)))/c))./rsa(:,i,m));
            end
            u = zeros((gr^2),1); % initializing true normal surface velocity on the actual surface
            for i=1:(L^2)
                u = u + (1/(1i*rho*freq))*(((zsg_v-zag_v)./rsa(:,i,m)).*((1i*2*pi*freq/c)+(1./rsa(:,i,m))).*ps);
            end
            U(:,fid,fm,m) = u; % true normal surface velocity matrix at frame fm
        
            % particle filtering for source amplitudes
            for j = 1:Nsa
                xp(:,m,j) = p(:,m,j) + normrnd(0, v1Cov, L^2, 1); % particle prediction
                yp = G(:,:,m)*xp(:,m,j) + normrnd(0, v2Cov, M^2, 1); % prediction measurement
                %prediction performance
                dx = xp(:,m,j) - x(:,m); % difference between particle prediction and true state vector
                % dy = yp - y; % difference between prediction measurement and observation
                % Covariance matrices
                Cxx(:,:,j) = Cxx(:,:,j) + v1Cov*eye(L^2);
                Cyx = (G(:,:,m)*Cxx(:,:,j));
                Cxy = Cyx';
                Cyy = G(:,:,m)*Cxx(:,:,j)*G(:,:,m)' + v2Cov*eye(M^2);
                temp1 = zeros((L^2),(L^2)); % temporary L*L matrix for inverse calculation later
                for g=1:(L^2)
                    for h=1:(L^2)
                        temp1(g,h) = Cxx(g,h,j);
                    end
                end
                mu = y + Cyx*(inv(temp1))*dx;
                C = Cyy - Cyx*(inv(temp1))*Cxy; % final covariance matrix used in likelihood function
%                 t1=inv(C);
%                 t2=(1/((2*pi)^(M^2)*det(C)))^(1/2);
%                 t3 = abs((y-mu)'*(t1)*(y-mu));
%                 t4 = exp((-1/2)*t3);
%                 t5 = abs(t2*t4);
%                 D=(C+C')/2;
%                 lf = abs(mvnpdf(yp',y',D));
%                 pd = abs(mvnpdf(xp(:,j)',x',(pCov*eye(L^2))));
                lf = abs(mvnpdf(yp',mu',((C+C')/2))); % likelihood function
                pd = abs((mvnpdf(xsg_v_p(:,m)',xsg_v_mean',(sxCov*eye(L^2))))*(mvnpdf(ysg_v_p(:,m)',ysg_v_mean',(syCov*eye(L^2))))); % prior distribution
                w(fm+1,m,j,:) = w(fm,m,j,:)*lf*pd; % weight update equation
%                 w(fm+1,j) = w(fm,j)*((1/sqrt(2*pi*v2Cov))*exp(-((abs(d_avg))^2)/(2*v2Cov)));
%                 page 324 of steven kay's book
%                 *v1())/(a+(b-a)*rand(1)); % weight = likelihood*transition/proposal
            end % End of particle filtering for source amplitudes
        
            % Weight normalization and state estimation
            temp2 = zeros((L^2),Nsa);
            for i=1:(L^2)
                w(fm+1,m,:,i) = w(fm+1,m,:,i)./(sum(w(fm+1,m,:,i))); % normalize the likelihood of each a priori estimate
                for k = 1:Nsa
                    v(i,k) = w(fm+1,m,k,i);
                    temp2(i,k) = p(i,m,k);
                end
                wm(i,m) = mean(w(fm+1, m, :, i));
                x_hat(i,m) = v(i,:)*temp2(i,:).'; % weighted sum: the state estimate is the weighted mean of the particles
            end
       
            X(:,fm+1,m) = x(:,m); X_hat(:,fm+1,m) = x_hat(:,m); Rs_v_hat(:,:,fm)=rs_v_hat; % update sequence of true states and state estimates
        
            % Error matrices for source amplitude and normal surface velocity
            error_amp(fid,fm+1,m) = (norm(X_hat(:,fm+1,m)-X(:,fm+1,m))/norm(X(:,fm+1,m)))*100; % source amplitude error matrix
%             error_loc(fid,fm) = (norm(rs_v_hat_mod-rs_v_mean_mod)/norm(rs_v_mean_mod))*100; % source amplitude error matrix
       
            % Sum of true and estimated normal surface velocity matrices over all frames 
            sum_vel = sum_vel + sum(U(:,fid,fm));
            sum_vel_hat = sum_vel_hat + sum(U_hat(:,fid,fm));
        
            % Resampling
            for i=1:(L^2)
                N_eff(i) = 1/(sum((w(fm+1,m,:,i).^2))); % effective sample size
                if N_eff(i)<N_thr % Resampling only required when effective sample size is less than threshold sample size
                    [temp2(i,:), v(i,:)] = resample(xp(i,:),v(i,:)); % Resampling
                    for k = 1:Nsa
                        w(fm+1,m,k,i) = v(i,k);
                        p(i,m,k) = temp2(i,k);
                    end
                else p=xp;
                end
            end % End of resampling loop
        end % End of particle filtering loop for source location
        
        temp3 = zeros(3,Nsl);
        for i = 1:(L^2)
            X_hatm(i,fm+1) = mean(X_hat(i, fm+1, :));
            for t1 = 1:3
                for t2 = 1:Nsl
                    temp3(t1, t2) = rs_v_p(i,t1,t2);
                end
            end
            wm(i,:) = wm(i,:)./(sum(wm(i,:))); % normalization of the mean weights
            rs_v_hat(i,:) = wm(i,:)*temp3.';
            rs_v_hat_mod(i) = (sum(rs_v_hat(i,:).^2))^(1/2);
        end
        rsa_hat = zeros(gr^2, L^2);
        for j=1:(L^2)
            for i=1:(gr^2)
                rsa_hat(i,j) = (sum((ra_v(i,:)-rs_v_hat(j,:)).^2))^(1/2);
            end
        end
     
        % mean (over source location particles) of error in source amplitude
        error_ampm(fid,fm+1) = mean(error_amp(fid,fm+1,:));
        
        % error in source location
        error_loc(fid,fm) = (norm(rs_v_hat_mod-rs_v_mean_mod)/norm(rs_v_mean_mod))*100;
        
        % mean (over source location particles) of true normal surface velocity
        for i=1:(gr^2)
            Um(i,fid,fm) = mean(U(i,fid,fm,:));
        end
        % Estimated normal surface velocity calculation
        ps_hat = zeros((gr^2),1); % initializing sound pressure on actual sound surface
        for i=1:(L^2)
            ps_hat = ps_hat + X_hatm(i,fm+1)*((exp((-1i*(2*pi*freq*rsa_hat(:,i)))/c))./rsa_hat(:,i));
        end
        u_hat = zeros((gr^2),1); % initializing true normal surface velocity on the actual surface
        for i=1:(L^2)
            u_hat = u_hat + (1/(1i*rho*freq))*(((zsg_v-zag_v)./rsa_hat(:,i)).*((1i*2*pi*freq/c)+(1./rsa_hat(:,i))).*ps_hat);
        end
        U_hat(:,fid,fm) = u_hat; % estimated normal surface velocity matrix at frame fm
        
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
    u_t(i,:) = ifft(Um(i, :, 99), np);
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
        pnew = zeros(1,Nsa); % initializing the particles after resampling
        wnew = zeros(1,Nsa); % initializing the weights after resampling
        % Cumulative Sum of Weights (CSW)
        CSW = zeros(1,Nsa); % initializing the CSW vector
        CSW(1) = w(1);
        for f = 2:Nsa
            CSW(f) = CSW(f-1) + w(f);
        end
        u(1) = rand()/Nsa; % drawing a sample from uniform distribution (0, 1/N)
        % main loop for resampling
        for e = 1:Nsa
            u(e) = u(1) + (e-1)/Nsa;
            f = 1;
            while u(e)>CSW(f)
                f = f+1;
            end
            pnew(e) = p(f);
            wnew(e) = 1/Nsa;
        end % End of resampling main loop
    end % End of resampling function
end % End of main function PF_v5