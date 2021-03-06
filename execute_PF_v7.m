clc;
clear all;
freq = 100:100:1700;
[xag, yag, X_hat, U, u_plot, u_hat_plot, error_amp1, error_amp_MVU, error_vel, w] = PF_v7(10);
%[xag, yag, X_hat, U, u_plot, u_hat_plot, error_amp2, error_amp_MVU, error_vel, w] = PF_v7(10);

figure, plot(freq,error_amp1(:,100), '-k');
axis([0 1800 0 0.3001]);
%title ('Percentage error in norm of the estimated virtual source amplitude vector at the last frame vs. frequency using the particle filter method with 10 and 100 particles');
% title ('Comparison between PF and MVU estimators: Percentage error in norm of the estimated virtual source amplitude vector at the last frame vs. frequency');
xlabel('Frequency (Hz)');
ylabel('Percentage error');
hold on; plot(freq,error_amp_MVU(:,100), '--k');
%hold on; plot(freq,error_amp2(:,100), '--k');
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