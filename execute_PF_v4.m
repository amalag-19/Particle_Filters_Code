clear all
freq = 100:100:1700;
[error_amp_MLE, error_amp_PF, error_vel_MLE, error_vel_PF] = PF_v4();

figure, plot(freq,error_amp_MLE(:,100), 'r-');
title ('percentage error in source amplitude of the last frame vs. frequency using MVU estimator');
xlabel('frequency (Hz)');
ylabel('%age error');

figure, plot(freq,error_amp_PF(:,100), 'r-');
title ('percentage error in source amplitude of the last frame vs. frequency using PF');
xlabel('frequency (Hz)');
ylabel('%age error');

figure, plot(freq,error_vel_MLE, 'r-');
title ('percentage error in velocity vs. frequency using MVU estimator');
xlabel('frequency (Hz)');
ylabel('%age error');

figure, plot(freq,error_vel_PF, 'r-');
title ('percentage error in velocity vs. frequency using PF');
xlabel('frequency (Hz)');
ylabel('%age error');

% gtext('maximum error');
% gtext('minimum error');
% frame = 1:101;
% surf(frame,freq,error)