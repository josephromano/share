% simulate the effect of a gw on light in a single arm of a GW interferometer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% laser and arm parameters
c = 1;
L = 1;
L_laser = L/2;
lambda_laser = L/2;
f_laser = c/lambda_laser;
T_laser = 1/f_laser;

% discrete times
t_tot = 8*L;
num_t = 1+36*floor(t_tot/T_laser); % 36 pt/cycle resolution
t = transpose(linspace(0, t_tot, num_t));
delta_t = t(2)-t(1);

% default GW amplitude
h0 = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GWtype = input('Input type of GW: step, ramp, square, triangle, sine\n','s');

switch GWtype
  case 'step'
    % step function
    tau = t_tot/2;
    ind = find(t>=tau);
    h = zeros(num_t,1);
    h(ind) = h0*ones(length(ind),1);

    % analytic expression for round-trip time
    analytic = zeros(num_t,1);
    analytic = 2*L/c + 0.5*h0*(t-t_tot/2).*(t_tot/2<t).*(t<=t_tot/2+2*L/c) ...
                     + h0*(L/c)*(t_tot/2+2*L/c<t).*(t<=t_tot);

  case 'ramp'
    % ramp gw
    h = h0*t/t(end);

    % analytic expression for round-trip time
    analytic = zeros(num_t,1);
    analytic = 2*L/c + 0.25*(h0/t_tot).*(t.^2).*(t<=2*L/c) ...
                     + (h0/t_tot)*(L/c).*(t-L/c).*(2*L/c<t).*(t<=t_tot);

  case 'square'
    % top-hat function gw
    tau1 = t_tot/4;
    tau2 = 3*t_tot/4;
    ind1 = find(t>=tau1);
    ind2 = find(t<=tau2);
    ind = intersect(ind1, ind2);
    h = zeros(num_t,1);
    h(ind) = h0*ones(length(ind),1);

    % analytic expression for round-trip time
    analytic = zeros(num_t,1);
    analytic = 2*L/c + 0.5*h0*(t-t_tot/4).*(t_tot/4<t).*(t<=t_tot/4+2*L/c) ...
                     + 0.5*h0*(2*L/c).*(t_tot/4+2*L/c<t).*(t<=3*t_tot/4) ...
                     + 0.5*h0*(3*t_tot/4-t+2*L/c) ...
                             .*(3*t_tot/4<t).*(t<=3*t_tot/4+2*L/c);

  case 'triangle'
    % triangular gw
    ind1 = find(t<t_tot/2);
    ind2 = setdiff([1:num_t], ind1);
    h = h0*[t(ind1)/(t_tot/2); (t(end)-t(ind2))/(t_tot/2)];

    % analytic expression for round-trip time
    analytic = zeros(num_t,1);
    analytic = 2*L/c + 0.25*(2*h0/t_tot)*(t.^2).*(t<=2*L/c) ...
                     + 0.25*(2*h0/t_tot)*(4*L/c)*(t-L/c)...
                          .*(2*L/c<t).*(t<=t_tot/2)...
                     + (0.25*t_tot*h0 -0.25*(2*h0/t_tot)...
                           *(2*t.^2 -2*t*(2*L/c+t_tot) +t_tot^2 +4*(L/c)^2))...
                          .*(t_tot/2<t).*(t<=t_tot/2+2*L/c)...
                     + 0.25*(2*h0/t_tot)*(4*L/c)*(L/c+t_tot-t)...
                          .*(t_tot/2+2*L/c<t);

 
  case 'sine'
    % sinusoidal gw
    lambda_gw = input('Input wavelength of GW in units of the arm length:\n');
    %lambda_gw = 1*L; 2, 4, 16
    f_gw = c/lambda_gw;
    T_gw = 1/f_gw;
    h = h0*sin(2*pi*f_gw*t);

    % analytic expression for round-trip time
    analytic = zeros(num_t,1);
    analytic = 2*L/c + 0.5*h0/(2*pi*f_gw)*(1-cos(2*pi*f_gw*t)).*(t<=2*L/c) ...
                     + 0.5*h0/(2*pi*f_gw)...
                          *(cos(2*pi*f_gw*(t-2*L/c)) - cos(2*pi*f_gw*t))...
                         .*(2*L/c<t);

  otherwise
    disp('unknown type');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stretch factor (small signal approx)
stretch_factor = 1+h/2;
stretch_factor0 = 1+h0/2;

% initialize discrete x's (proper distance)
delta_x = c*delta_t; % proper speed of light is always = c 
num_x = 1+round(L/delta_x);
x = transpose(linspace(0, L, num_x)); 
num_x_laser = 1+round(L_laser/delta_x);
x_laser = transpose(linspace(-L_laser, 0, num_x_laser));

% intialize EM waves in laser 
y_laser = sin((2*pi/lambda_laser)*x_laser);

% intialize EM waves in interferometer arm
y_l = zeros(num_x, 1);
y_r = zeros(num_x, 1);
for k=1:floor((2*L)/(c*delta_t))
  y_l = [y_l(2:end); -y_r(end)];
  y_r = [sin(-2*pi*f_laser*t(k)); y_r(1:end-1)];
end

% initialize some other variables
counter = 0;
start_times = [];
stop_times  = [];

% define figure
fig=figure(1);

% plot analytic expression for round-trip time
subplot(3,1,3)
ndx = find(t>=2);
plot(t(ndx), analytic(ndx), 'b-');

% loop over discrete times
for k=1:num_t

  % create new wave profiles
  y_l = [y_l(2:end); -y_r(end)];
  y_r = [sin(-2*pi*f_laser*t(k)); y_r(1:end-1)];
  y_laser = sin((2*pi/lambda_laser)*x_laser-2*pi*f_laser*t(k));

  % set start time of each new wave front
  if round(t(k)/T_laser)==t(k)/T_laser
    start_times = [start_times; t(k)];
  end

  % stretch (or compress) the light waves in the arm
  num_x = 1+round(L*stretch_factor(k)/delta_x);
  if k==1
    x_i = transpose(linspace(0, L, num_x));
  else
    x_i = transpose(linspace(0, L*stretch_factor(k-1), num_x));
  end

  % save old x's and construct new ones
  x_old = x;
  x = transpose(linspace(0, L*stretch_factor(k), num_x));

  % interpolate
  y_r = interp1(x_old, y_r, x_i, 'spline');
  y_l = interp1(x_old, y_l, x_i, 'spline');

  % determine location of wave fronts from zero crossings
  [xc_r, indc_r] = findzerocrossings(x, y_r, 'decreasing');
  [xc_l, indc_l] = findzerocrossings(x, y_l, 'increasing');

  % determine stop time of each returning wave front
  if indc_l(1)==1
    tolerance = 0.1;
    if t(k)<(1-tolerance)*2*L/c
      % ignore first few returning wave fronts 
    else
      stop_times = [stop_times; t(k)];
      wf_num = length(stop_times);
      roundtrip_time = stop_times(wf_num)-start_times(wf_num);
      subplot(3,1,3)
      plot(stop_times(wf_num), roundtrip_time, '*r');
      %timelabel=['roundtrip = ' num2str(roundtrip_time,'%1.2f') ' sec'];
      %legend(timelabel, 'Location', 'South');
    end
  end

  % make plots
  subplot(3,1,2)
  hold off
  plot(x_laser, y_laser, 'b', ...
       x, y_r, 'b', x, y_l, 'r', ... 
       xc_r, zeros(size(xc_r)), '*b', xc_l, zeros(size(xc_l)), '*r');
  hold on
  vline(L*stretch_factor(k),'k-');
  vline(0,'k-');
  xlim([-0.5 1.1*max(L,L*stretch_factor0)]);
  ylim([-1.2 1.2]);
  xlabel('Arm length (proper distance)');
  ylabel('EM wave');
  grid on;

  subplot(3,1,1);
  hold on
  plot(t(k), h(k), 'b.');
  xlim([0 t_tot]);
  ylim([-1.2*abs(h0) 1.2*abs(h0)]);
  xlabel('Time (sec)');
  ylabel('GW strain');
  grid on

  subplot(3,1,3)
  hold on
  xlim([0 t_tot]);
  %ylim([0 1.2*2*max(L,L*stretch_factor0)/c]);
  ylim([1.75*L/c 2.25*L/c]);
  xlabel('EM wave arrival time (sec)');
  ylabel('Roundtrip (sec)');
  grid on

  M(k) = getframe(fig);
  %pause(.001);

end

fname = ['interferometerArm_' GWtype '.avi']; 
v = VideoWriter(fname);
open(v)
writeVideo(v, M)
print -depsc2 -f1 interferometerArm
close(v)

return

