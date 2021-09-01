function orf_int = overlapIntegrand(det1, det2, beta)
%
% calculate (4\pi/5) * root-mean-absolute square of the DC
% antenna pattern for a pair of detectors, as a function of 
% declination
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458;   % speed of light (m/s)

% get detector geometry information
[r1, u1, v1, T1] = getdetectorNew(det1);
[r2, u2, v2, T2] = getdetectorNew(det2);
deltaX = r1-r2;

% (theta, phi) values
eps = 0.001; % small number
Ntheta = 181;
Nphi = 361;
theta = linspace(eps, pi-eps, Ntheta);
phi = linspace(0, 2*pi, Nphi);

% construct overlap reduction function
orf_int = zeros(Ntheta, Nphi);
rms = zeros(1, Ntheta);

% DC component 
f = 0;  % Hz
method = 'lw';

for ii = 1:1:Ntheta
  
  for jj = 1:1:Nphi
    fprintf('working on %d of %d\n', ii, Ntheta);

    [F1p, F1c] = FpFc(f, theta(ii), phi(jj), 0, u1, v1, T1, method, 'id');
    [F2p, F2c] = FpFc(f, theta(ii), phi(jj), 0, u2, v2, T2, method, 'id');

    H = (1/(sin(beta)^2)) * abs( (1/2) * (F1p.*conj(F2p) + F1c.*conj(F2c)) )^2;
    orf_int(ii,jj) = H;

  end

  % average over phi
  rms(ii) = sqrt(sum(orf_int(ii,:)/Nphi)); 

end

% plot rms as a function of dec
figure(1)
dec = pi/2 - theta;
plot(dec*180/pi, rms, 'b*');
xlabel('declination (degrees)');
ylabel('(4 pi/5) RMS of d2 gamma/d2 Omega at f=0');
grid on
xlim([0 90])
ylim([0 0.3])
filename = ['rms_dec_' det1 '_' det2 '.pdf'];
titlestr = [det1 ' - ' det2];
title(titlestr, 'fontsize', 20);
print('-dpdf', filename);

return
