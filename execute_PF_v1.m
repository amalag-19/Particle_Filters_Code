clear all
freq = 100:100:1700;
[error]=PF_v2();
figure;
plot(freq,error(:,101), 'r-')
title ('percentage error in source amplitude of the last frame vs. frequency')
xlabel('frequency (Hz)')
ylabel('%age error')
gtext('maximum error')
gtext('minimum error')
% frame = 1:101;
% surf(frame,freq,error)