[r,lags]=xcorr(filtered_C3, filtered_C4);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);ylim([-2*1e8 2*1e8]);
title('delta')