clc;
clear all;
freq = 100:100:1700;
[xag, yag, X_hat, U, u_plot, u_hat_plot, error_amp, error_loc, error_vel, w] = PF_v8();

figure, plot(freq,error_amp(:,100), 'r-');
title ('percentage error in source amplitude of the last frame vs. frequency');
xlabel('frequency (Hz)');
ylabel('%age error');

figure, plot(freq,error_loc(:,99), 'r-');
title ('percentage error in source location of the last frame vs. frequency');
xlabel('frequency (Hz)');
ylabel('%age error');

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