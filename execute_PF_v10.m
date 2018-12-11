clc;
clear all;
freq = 100:100:1700;
[xag, yag, Xsa_hat, U, u_plot, u_hat_plot, error_amp, error_amp_MVU, error_loc, error_vel, w] = PF_v10();

figure, plot(freq,error_amp(:,100), '-k');
axis([0 1800 0 2.5]);
%title ('Comparison between PF and MVU estimators: Percentage error in norm of the estimated virtual source amplitude vector at the last frame vs. frequency');
xlabel('Frequency (Hz)');
ylabel('Percentage error');
hold on; plot(freq,error_amp_MVU(:,100), '--k');

figure, plot(freq,error_loc(:,99), '-k');
%title ('Percentage error in norm of the estimated virtual source location vector by PF method at the last frame vs. frequency');
xlabel('Frequency (Hz)');
ylabel('Percentage error');

figure;contourf(xag,yag,u_plot),colorbar,shading flat,xlabel('x-coordinate (m)');ylabel('y-coordinate (m)');title('Desired Normal Surface velocity (m/s, Linear)');
figure;contourf(xag,yag,u_hat_plot),colorbar,shading flat,xlabel('x-coordinate (m)');ylabel('y-coordinate (m)');title('Reconstructed Normal Surface velocity (m/s, Linear)');

% figure;contourf(xg,yg,U(:,:,99)),colorbar,shading flat,xlabel('x-coordinate (m)');ylabel('y-coordinate (m)');title('Particle velocity (m/s, Linear)');
% figure;contourf(xg,yg,U_hat(:,:,99)),colorbar,shading flat,xlabel('x-coordinate (m)');ylabel('y-coordinate (m)');title('Particle velocity (m/s, Linear)');

% figure, plot(freq,error_vel, 'r-');
% title ('percentage error in velocity vs. frequency');
% xlabel('frequency (Hz)');
% ylabel('%age error');

% gtext('maximum error');
% gtext('minimum error');
% frame = 1:101;
% surf(frame,freq,error)