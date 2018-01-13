% This is code for simulation of propagation of 2-D SH-waves with CPML ABC
% employing memory variables. It was written by Daniel Grzech as part of
% MSc Individual Project supervised by Dr Panagiotis Kosmas.

% 24th August 2016

close all;
clear all;
running_time = tic;

%=============================basic variables=============================%

%x_length = 0.15;
%big_delta_x = 0.0015;
%z_length = 0.15;
%big_delta_z = 0.0015;

%i_length = x_length / big_delta_x;
%k_length = z_length / big_delta_z;

i_length = 200;
k_length = 200;
big_delta_x = 0.00075;
big_delta_z = big_delta_x;

no_time_steps = 1000;
delta_t = big_delta_x / 2;
t = no_time_steps * delta_t;

fprintf('no. of time steps %d \n', no_time_steps);
disp('');

rho = 1000;
mu = 1000;

%====================CPML thickness in each direction=====================%

cpml_thickness = 10;
x_cpml_1 = cpml_thickness;
x_cpml_2 = cpml_thickness;
z_cpml_1 = cpml_thickness;
z_cpml_2 = cpml_thickness;

i_length_total = i_length + x_cpml_1 + x_cpml_2;
k_length_total = k_length + z_cpml_1 + z_cpml_2;

x_computational_zone = (x_cpml_1+1):(i_length_total-x_cpml_2);
z_computational_zone = (z_cpml_1+1):(k_length_total-z_cpml_2);

%===========================antenna information===========================%

no_antennas = 1;

% position of source--spherical wavefront

source_pos = [round(i_length_total / 2) round(k_length_total /2)];

%============================source parameters============================%

omega_p = 2*pi*100;
f_a = 1;
delay = 50 * delta_t;
source_val = (1 - 1/2 * omega_p^2 * ((0:no_time_steps - 1) * delta_t - delay).^2) .* exp(-1/4 * omega_p^2 * ((0:no_time_steps - 1) * delta_t - delay).^2);

figure
ax1 = gca;
plot(ax1, (0:no_time_steps - 1) * delta_t, source_val);
xlabel('time (s)')
ylabel('velocity (m/s)')
grid on;
drawnow;

%=================================sensor==================================%

sensor_pos = [x_cpml_1+5 round(k_length_total / 2)];
sensor_1 = [];
sensor_2 = [];
sensor_3 = [];

%=============================CPML parameters=============================%

n_1 = 3; n_2 = 0; n_3 = 3;

kappa_x_max = 0;
sigma_x_1_max = 0.8 * (n_1 + 1) / (sqrt(mu/rho) * big_delta_x);
sigma_x_2_max = 0.8 * (n_1 + 1) / (sqrt(mu/rho) * big_delta_x);
alpha_x_max = 2 * pi * f_a;

kappa_z_max = 0;
sigma_z_1_max = 0.8 * (n_1 + 1) / (sqrt(mu/rho) * big_delta_z);
sigma_z_2_max = 0.8 * (n_1 + 1) / (sqrt(mu/rho) * big_delta_z);
alpha_z_max = 2 * pi * f_a;

kappa_x_1_v_y = zeros(x_cpml_1, 1);
kappa_x_2_v_y = zeros(x_cpml_2, 1);
kappa_z_1_v_y = zeros(z_cpml_1, 1);
kappa_z_2_v_y = zeros(z_cpml_2, 1);
sigma_x_1_v_y = zeros(x_cpml_1, 1);
sigma_x_2_v_y = zeros(x_cpml_2, 1);
sigma_z_1_v_y = zeros(z_cpml_1, 1);
sigma_z_2_v_y = zeros(z_cpml_2, 1);
alpha_x_1_v_y = zeros(x_cpml_1, 1);
alpha_x_2_v_y = zeros(x_cpml_2, 1);
alpha_z_1_v_y = zeros(z_cpml_1, 1);
alpha_z_2_v_y = zeros(z_cpml_2, 1);
beta_x_1_v_y = zeros(x_cpml_1, 1);
beta_x_2_v_y = zeros(x_cpml_2, 1);
delta_x_1_v_y = zeros(x_cpml_1, 1);
delta_x_2_v_y = zeros(x_cpml_2, 1);
beta_z_1_v_y = zeros(z_cpml_1, 1);
beta_z_2_v_y = zeros(z_cpml_2, 1);
delta_z_1_v_y = zeros(z_cpml_1, 1);
delta_z_2_v_y = zeros(z_cpml_2, 1);
b_x_1_v_y = zeros(x_cpml_1, 1);
c_x_1_v_y = zeros(x_cpml_1, 1);
b_x_2_v_y = zeros(x_cpml_2, 1);
c_x_2_v_y = zeros(x_cpml_2, 1);
b_z_1_v_y = zeros(z_cpml_1, 1);
c_z_1_v_y = zeros(z_cpml_1, 1);
b_z_2_v_y = zeros(z_cpml_2, 1);
c_z_2_v_y = zeros(z_cpml_2, 1);

kappa_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
kappa_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
sigma_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
sigma_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
alpha_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
alpha_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
beta_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
beta_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
delta_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
delta_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
b_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
c_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
b_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
c_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);

kappa_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
kappa_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
sigma_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
sigma_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
alpha_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
alpha_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
beta_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
beta_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
delta_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
delta_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
b_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
c_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
b_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
c_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);

% x-axis--left

for i = 1:x_cpml_1
    
    kappa_x_1_v_y(i) = 1.0 + kappa_x_max * ((x_cpml_1 - i) / (x_cpml_1 - 1.0))^n_1;
    sigma_x_1_v_y(i) = sigma_x_1_max * ((x_cpml_1 - i) / (x_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_x_1_v_y(i) = alpha_x_max * ((i - 1.0) / (x_cpml_1 - 1.0))^(n_3);
    beta_x_1_v_y(i) = sigma_x_1_v_y(i) / kappa_x_1_v_y(i);
    delta_x_1_v_y(i) = sigma_x_1_v_y(i) / kappa_x_1_v_y(i) + alpha_x_1_v_y(i);
    
    b_x_1_v_y(i) = (2.0 - delta_t * beta_x_1_v_y(i)) / (2.0 + delta_t * beta_x_1_v_y(i));
    c_x_1_v_y(i) = (2.0 * delta_t * -delta_x_1_v_y(i)) / (2.0 + delta_t * beta_x_1_v_y(i));
     
end

for i = 1:x_cpml_1 - 1
   
    kappa_x_1_tau_yx(i) = 1.0 + kappa_x_max * ((x_cpml_1 - i - 0.5) / (x_cpml_1 - 1.0))^n_1;
    sigma_x_1_tau_yx(i) = sigma_x_1_max * ((x_cpml_1 - i - 0.5) / (x_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_x_1_tau_yx(i) = alpha_x_max * ((i - 0.5) / (x_cpml_1 - 1.0))^(n_3);
    beta_x_1_tau_yx(i) = sigma_x_1_tau_yx(i) / kappa_x_1_tau_yx(i);
    delta_x_1_tau_yx(i) = sigma_x_1_tau_yx(i) / kappa_x_1_tau_yx(i) + alpha_x_1_tau_yx(i);
    
    b_x_1_tau_yx(i) = (2.0 - delta_t * beta_x_1_tau_yx(i)) / (2.0 + delta_t * beta_x_1_tau_yx(i));
    c_x_1_tau_yx(i) = (2.0 * delta_t * -delta_x_1_tau_yx(i)) / (2.0 + delta_t * beta_x_1_tau_yx(i));
    
end

% x-axis--right

for i = 1:x_cpml_2
    
    kappa_x_2_v_y(i) = 1.0 + kappa_x_max * ((x_cpml_2 - i) / (x_cpml_2 - 1))^n_1;
    sigma_x_2_v_y(i) = sigma_x_2_max * ((x_cpml_2 - i) / (x_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_x_2_v_y(i) = alpha_x_max * ((i - 1.0) / (x_cpml_2 - 1.0))^(n_3);
    beta_x_2_v_y(i) = sigma_x_2_v_y(i) / kappa_x_2_v_y(i);
    delta_x_2_v_y(i) = sigma_x_2_v_y(i) / kappa_x_2_v_y(i) + alpha_x_2_v_y(i);
    
    b_x_2_v_y(i) = (2.0 - delta_t * beta_x_2_v_y(i)) / (2.0 + delta_t * beta_x_2_v_y(i));
    c_x_2_v_y(i) = (2.0 * delta_t * -delta_x_2_v_y(i)) / (2.0 + delta_t * beta_x_2_v_y(i));
    
end

for i = 1:x_cpml_2 - 1
   
    kappa_x_2_tau_yx(i) = 1.0 + kappa_x_max * ((x_cpml_2 - i - 0.5) / (x_cpml_2 - 1.0))^n_1;
    sigma_x_2_tau_yx(i) = sigma_x_2_max * ((x_cpml_2 - i - 0.5) / (x_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_x_2_tau_yx(i) = alpha_x_max * ((i - 0.5) / (x_cpml_2 - 1.0))^(n_3);
    beta_x_2_tau_yx(i) = sigma_x_2_tau_yx(i) / kappa_x_2_tau_yx(i);
    delta_x_2_tau_yx(i) = sigma_x_2_tau_yx(i) / kappa_x_2_tau_yx(i) + alpha_x_2_tau_yx(i);
    
    b_x_2_tau_yx(i) = (2 - delta_t * beta_x_2_tau_yx(i)) / (2 + delta_t * beta_x_2_tau_yx(i));
    c_x_2_tau_yx(i) = (2 * delta_t * -delta_x_2_tau_yx(i)) / (2 + delta_t * beta_x_2_tau_yx(i));
    
end

% z-axis--bottom

for k = 1:z_cpml_1
    
    kappa_z_1_v_y(k) = 1.0 + kappa_z_max * ((z_cpml_1 - k) / (z_cpml_1 - 1.0))^n_1;
    sigma_z_1_v_y(k) = sigma_z_1_max * ((z_cpml_1 - k) / (z_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_z_1_v_y(k) = alpha_z_max * ((k - 1.0) / (z_cpml_1 - 1.0))^(n_3);
    beta_z_1_v_y(k) = sigma_z_1_v_y(k) / kappa_z_1_v_y(k);
    delta_z_1_v_y(k) = sigma_z_1_v_y(k) / kappa_z_1_v_y(k) + alpha_z_1_v_y(k);
    
    b_z_1_v_y(k) = (2 - delta_t * beta_z_1_v_y(k)) / (2 + delta_t * beta_z_1_v_y(k));
    c_z_1_v_y(k) = (2 * delta_t * -delta_z_1_v_y(k)) / (2 + delta_t * beta_z_1_v_y(k));
    
end

for k = 1:z_cpml_1 - 1
   
    kappa_z_1_tau_yz(k) = 1.0 + kappa_z_max * ((z_cpml_1 - k - 0.5) / (z_cpml_1 - 1.0))^n_1;
    sigma_z_1_tau_yz(k) = sigma_z_1_max * ((z_cpml_1 - k - 0.5) / (z_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_z_1_tau_yz(k) = alpha_z_max * ((k - 0.5) / (z_cpml_1 - 1.0))^(n_3);
    beta_z_1_tau_yz(k) = sigma_z_1_tau_yz(k) / kappa_z_1_tau_yz(k);
    delta_z_1_tau_yz(k) = sigma_z_1_tau_yz(k) / kappa_z_1_tau_yz(k) + alpha_z_1_tau_yz(k);
    
    b_z_1_tau_yz(k) = (2 - delta_t * beta_z_1_tau_yz(k)) / (2 + delta_t * beta_z_1_tau_yz(k));
    c_z_1_tau_yz(k) = (2 * delta_t * -delta_z_1_tau_yz(k)) / (2 + delta_t * beta_z_1_tau_yz(k));
    
end

% z-axis--top

for k = 1:z_cpml_2
    
    kappa_z_2_v_y(k) = 1.0 + kappa_z_max * ((z_cpml_2 - k) / (z_cpml_2 - 1.0))^n_1;
    sigma_z_2_v_y(k) = sigma_z_2_max * ((z_cpml_2 - k) / (z_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_z_2_v_y(k) = alpha_z_max * ((k - 1.0) / (z_cpml_1 - 1.0))^(n_3);
    beta_z_2_v_y(k) = sigma_z_2_v_y(k) / kappa_z_2_v_y(k);
    delta_z_2_v_y(k) = sigma_z_2_v_y(k) / kappa_z_2_v_y(k) + alpha_z_2_v_y(k);
    
    b_z_2_v_y(k) = (2 - delta_t * beta_z_2_v_y(k)) / (2 + delta_t * beta_z_2_v_y(k));
    c_z_2_v_y(k) = (2 * delta_t * -delta_z_2_v_y(k)) / (2 + delta_t * beta_z_2_v_y(k));
 
end

for k = 1:z_cpml_2 - 1
   
    kappa_z_2_tau_yz(k) = 1.0 + kappa_z_max * ((z_cpml_2 - k - 0.5) / (z_cpml_2 - 1.0))^n_1;
    sigma_z_2_tau_yz(k) = sigma_z_2_max * ((z_cpml_2 - k - 0.5) / (z_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_z_2_tau_yz(k) = alpha_z_max * ((k - 0.5) / (z_cpml_2 - 1.0))^(n_3);
    beta_z_2_tau_yz(k) = sigma_z_2_tau_yz(k) / kappa_z_2_tau_yz(k);
    delta_z_2_tau_yz(k) = sigma_z_2_tau_yz(k) / kappa_z_2_tau_yz(k) + alpha_z_2_tau_yz(k);
       
    b_z_2_tau_yz(k) = (2 - delta_t * beta_z_2_tau_yz(k)) / (2 + delta_t * beta_z_2_tau_yz(k));
    c_z_2_tau_yz(k) = (2 * delta_t * -delta_z_2_tau_yz(k)) / (2 + delta_t * beta_z_2_tau_yz(k));
    
end

kappa_x_v_y = zeros(i_length_total, k_length_total);
sigma_x_v_y = zeros(i_length_total, k_length_total);
alpha_x_v_y = zeros(i_length_total, k_length_total);
kappa_z_v_y = zeros(i_length_total, k_length_total);
sigma_z_v_y = zeros(i_length_total, k_length_total);
alpha_z_v_y = zeros(i_length_total, k_length_total);
kappa_x_tau_yx = zeros(i_length_total, k_length_total);
sigma_x_tau_yx = zeros(i_length_total, k_length_total);
alpha_x_tau_yx = zeros(i_length_total, k_length_total);
kappa_z_tau_yz = zeros(i_length_total, k_length_total);
sigma_z_tau_yz = zeros(i_length_total, k_length_total);
alpha_z_tau_yz = zeros(i_length_total, k_length_total);

kk = z_cpml_2 - 1;

for k = 1:k_length_total - 1
    
    if k <= z_cpml_1 - 1
        
        kappa_x_tau_yx(:, k) = kappa_x_1_tau_yx(k);
        sigma_x_tau_yx(:, k) = sigma_x_1_tau_yx(k);
        alpha_x_tau_yx(:, k) = alpha_x_1_tau_yx(k);
        
    elseif k >= k_length_total + 1 - z_cpml_2
        
        kappa_x_tau_yx(:, k) = kappa_x_2_tau_yx(kk);
        sigma_x_tau_yx(:, k) = sigma_x_2_tau_yx(kk);
        alpha_x_tau_yx(:, k) = alpha_x_2_tau_yx(kk);
        kk = kk - 1;
        
    else
        
        kappa_x_tau_yx(:, k) = 1.0;
        sigma_x_tau_yx(:, k) = 0;
        alpha_x_tau_yx(:, k) = 0;
        
    end
    
end

ii = x_cpml_2 - 1;

for i = 1:i_length_total - 1
    
    if i <= x_cpml_1 - 1
        
        kappa_z_tau_yz(i, :) = kappa_z_1_tau_yz(i);
        sigma_z_tau_yz(i, :) = sigma_z_1_tau_yz(i);
        alpha_z_tau_yz(i, :) = alpha_z_1_tau_yz(i);
        
    elseif i >= i_length_total + 1 - x_cpml_2
        
        kappa_z_tau_yz(i, :) = kappa_z_2_tau_yz(ii);
        sigma_z_tau_yz(i, :) = sigma_z_2_tau_yz(ii);
        alpha_z_tau_yz(i, :) = alpha_z_2_tau_yz(ii);
        ii = ii - 1;
        
    else
        
        kappa_z_tau_yz(i, :) = 1.0;
        sigma_z_tau_yz(i, :) = 0;
        alpha_z_tau_yz(i, :) = 0;
        
    end
    
end

kk = z_cpml_2;

for k = 2:k_length_total - 1
   
   if k <= z_cpml_1
       
       kappa_x_v_y(:, k) = kappa_x_1_v_y(k);
       sigma_x_v_y(:, k) = sigma_x_1_v_y(k);
       alpha_x_v_y(:, k) = alpha_x_1_v_y(k);
       
   elseif k >= k_length_total + 1 - z_cpml_2
       
       kappa_x_v_y(:, k) = kappa_x_2_v_y(kk);
       sigma_x_v_y(:, k) = sigma_x_2_v_y(kk);
       alpha_x_v_y(:, k) = alpha_x_2_v_y(kk);
       kk = kk - 1;
       
   else
       
       kappa_x_v_y(:, k) = 1.0;
       sigma_x_v_y(:, k) = 0;
       alpha_x_v_y(:, k) = 0;
       
   end
       
end

ii = x_cpml_2;
   
for i = 2:i_length_total - 1
       
       if i <= x_cpml_1
           
           kappa_z_v_y(i, :) = kappa_z_1_v_y(i);
           sigma_z_v_y(i, :) = sigma_z_1_v_y(i);
           alpha_z_v_y(i, :) = alpha_z_1_v_y(i);
           
       elseif i >= i_length_total + 1 - x_cpml_2
           
           kappa_z_v_y(i, :) = kappa_z_2_v_y(ii);
           sigma_z_v_y(i, :) = sigma_z_2_v_y(ii);
           alpha_z_v_y(i, :) = alpha_z_2_v_y(ii);
           ii = ii - 1;
           
       else
          
           kappa_z_v_y(i, :) = 1.0;   
           sigma_z_v_y(i, :) = 0;
           alpha_z_v_y(i, :) = 0;
           
       end
        
end

%======matrices for u_y, v_y, tau_yx, tau_yz, A_x, B_z, C_x, and D_z======%

u_y = zeros(i_length_total, k_length_total);
v_y = zeros(i_length_total, k_length_total);
tau_yx = zeros(i_length_total, k_length_total);
tau_yz = zeros(i_length_total, k_length_total);

% CPML

A_x_1 = zeros(x_cpml_1, k_length_total);
A_x_2 = zeros(x_cpml_2, k_length_total);
C_x_1 = zeros(x_cpml_1 - 1, k_length_total);
C_x_2 = zeros(x_cpml_2 - 1, k_length_total);
B_z_1 = zeros(i_length_total, z_cpml_1);
B_z_2 = zeros(i_length_total, z_cpml_2);
D_z_1 = zeros(i_length_total, z_cpml_1 - 1); 
D_z_2 = zeros(i_length_total, z_cpml_2 - 1);

%======================propagation in time and space======================%

plot_freq = 100;

u_y_min = sum(source_val(source_val < 0)) * delta_t * 1.25;
u_y_max = sum(source_val(source_val > 0)) * delta_t * 1.25;

% name = '2d_sh_wave_ade_spherical.avi';
% name = '2d_sh_wave_ade_plane.avi';
% v = VideoWriter(name); % create object to write visualisation to file
% open(v);

for n = 1:no_time_steps
    
    % update tau_yx and tau_yz

    i = 1:i_length_total - 1;
    k = 1:k_length_total - 1;
    
    tau_yx(i, k) = tau_yx(i, k) + mu .* 1 ./ kappa_x_tau_yx(i, k) * delta_t .* (v_y(i + 1, k) - v_y(i, k)) ./ big_delta_x;
    
    tau_yz(i, k) = tau_yz(i, k) + mu .* 1 ./ kappa_z_tau_yz(i, k) * delta_t .* (v_y(i, k + 1) - v_y(i, k)) ./ big_delta_z;
        
    % update C_x
    
    % left tau_yx
    
    k = 1:k_length_total - 1;
    
    for i = 1:x_cpml_1 - 1
        
        C_x_1_old = C_x_1;
        C_x_1(i, k) = b_x_1_tau_yx(i) .* C_x_1(i, k) + c_x_1_tau_yx(i) .* (v_y(i + 1, k) - v_y(i, k)) / big_delta_x;
        tau_yx(i, k) = tau_yx(i, k) + delta_t * mu * 1 ./ kappa_x_tau_yx(i, k) .* 1/2 .* (C_x_1(i, k) + C_x_1_old(i, k));
        
    end
        
    % right tau_yx
    
    ii = x_cpml_2 - 1;
    k = 1:k_length_total - 1;
    
    for i = i_length_total + 1 - x_cpml_2:i_length_total - 1
        
        C_x_2_old = C_x_2;
        C_x_2(ii, k) = b_x_2_tau_yx(ii) .* C_x_2(ii, k) + c_x_2_tau_yx(ii) .*  (v_y(i + 1, k) - v_y(i, k)) / big_delta_x;
        tau_yx(i, k) = tau_yx(i, k) + delta_t * mu * 1 ./ kappa_x_tau_yx(i, k) .* 1/2 .* (C_x_2(ii, k) + C_x_2_old(ii, k));
        ii = ii - 1;
        
    end
    
    % update D_z
    
    % bottom tau_yz
    
    i = 1:i_length_total - 1;
    
    for k = 1:z_cpml_1 - 1
       
        D_z_1_old = D_z_1;
        D_z_1(i, k) = b_z_1_tau_yz(k) .* D_z_1(i, k) + c_z_1_tau_yz(k) .* (v_y(i, k + 1) - v_y(i, k)) / big_delta_z;
        tau_yz(i, k) = tau_yz(i, k) + delta_t * mu .* 1 ./ kappa_z_tau_yz(i, k) .* 1/2 .* (D_z_1(i, k) + D_z_1_old(i, k));
        
    end
    
    % top tau_yz
    
    kk = z_cpml_2 - 1;
    i = 1:i_length_total - 1;
    
    for k = k_length_total + 1 - z_cpml_2:k_length_total - 1
        
        D_z_2_old = D_z_2;
        D_z_2(i, kk) = b_z_2_tau_yz(kk) .* D_z_2(i, kk) + c_z_2_tau_yz(kk) .* (v_y(i, k + 1) - v_y(i, k)) / big_delta_z;
        tau_yz(i, k) = tau_yz(i, k) + delta_t * mu .* 1 ./ kappa_z_tau_yz(i, k) .* 1/2 .* (D_z_2(i, kk) + D_z_2_old(i, kk));
        kk = kk - 1;
        
    end

    % update v_y
    
    i = 2:i_length_total - 1;
    k = 2:k_length_total - 1;
        
    v_y(i, k) = v_y(i, k) + delta_t * 1 / rho * (1 ./ kappa_x_v_y(i, k) .* (tau_yx(i, k) - tau_yx(i - 1, k)) ./ big_delta_x + 1 ./ kappa_z_v_y(i, k) .* (tau_yz(i, k) - tau_yz(i, k - 1)) ./ big_delta_z);
        
    % update A_x and B_z
    
    % A_x for left v_y
    
    k = 2:k_length_total - 1;
        
    for i = 2:x_cpml_1
    
        A_x_1_old = A_x_1;
        A_x_1(i, k) = b_x_1_v_y(i) .* A_x_1(i, k) + c_x_1_v_y(i) .* (tau_yx(i, k) - tau_yx(i - 1, k)) ./ big_delta_x;
        v_y(i, k) = v_y(i, k) + delta_t ./ rho .* 1 ./ kappa_x_v_y(i, k) .* 1/2 .* (A_x_1(i, k) + A_x_1_old(i, k));
        
    end
        
    % A_x for right v_y
    
    ii = x_cpml_2;
    k = 2:k_length_total - 1;
    
    for i = i_length_total + 1 - x_cpml_2:i_length_total - 1
       
        A_x_2_old = A_x_2;
        A_x_2(ii, k) = b_x_2_v_y(ii) .* A_x_2(ii, k) + c_x_2_v_y(ii) .* (tau_yx(i, k) - tau_yx(i - 1, k)) / big_delta_x;
        v_y(i, k) = v_y(i, k) + delta_t / rho .* 1 ./ kappa_x_v_y(i, k) .* 1/2 .* (A_x_2(ii,k) + A_x_2_old(ii, k));
        ii = ii - 1;
        
    end
     
    % B_z for bottom v_y
    
    i = 2:i_length_total - 1;
    
    for k = 2:z_cpml_1
        
        B_z_1_old = B_z_1;
        B_z_1(i, k) = b_z_1_v_y(k) .* B_z_1(i, k) + c_z_1_v_y(k) .* (tau_yz(i, k) - tau_yz(i, k - 1)) / big_delta_z;
        v_y(i, k) = v_y(i, k) + delta_t / rho .* 1 ./ kappa_z_v_y(i, k) .* 1/2 .* (B_z_1(i, k) + B_z_1_old(i, k));
        
    end
    
    % B_z for top v_y
    
    kk = z_cpml_2;
    i = 2:i_length_total - 1;
    
    for k = k_length_total + 1 - z_cpml_2:k_length_total - 1
        
        B_z_2_old = B_z_2;
        B_z_2(i, kk) = b_z_2_v_y(kk) .* B_z_2(i, kk) + c_z_2_v_y(kk) .* (tau_yz(i, k) - tau_yz(i, k - 1)) / big_delta_z;
        v_y(i, k) = v_y(i, k) + delta_t / rho .* 1 ./ kappa_z_v_y(i, k) .* 1/2 .* (B_z_2(i, kk) + B_z_2_old(i, kk));
        kk = kk - 1;
    
    end
    
    % source--spherical wavefront
   
%     v_y(source_pos(1), source_pos(2)) = v_y(source_pos(1), source_pos(2)) + source_val(n);
    
    % source--plane wavefront
    
    v_y(round(i_length_total / 2), :) = v_y(round(i_length_total / 2), :) + source_val(n);
    
    u_y = u_y + v_y * delta_t;
    
    % sensor
    
    sensor_1 = [sensor_1 v_y(sensor_pos(1,1), sensor_pos(1,2))];
    
    % plot current image every plot_freq time steps
    
    if mod(n, plot_freq) == 0
          
          disp(n)
          
          figure
          ax2 = gca;
          temp = double(v_y(1:i_length_total, 1:k_length_total));
          imagesc(ax2, [big_delta_x big_delta_x * i_length_total], [big_delta_z big_delta_z * k_length_total], (squeeze(temp).'));
          cmax = max(max(max(temp)), -min(min(temp)));
          caxis([-cmax cmax]);
          title('CPML for 2-D SH-waves');
          xlabel('x (m)');
          ylabel('z (m)');
          colorbar;
          axis image;
          
    end
    
%     % visualise wave propagation
%     
%     figure(2)
%     ax3 = gca;
%     axis tight manual;
%     ax3.NextPlot = 'replaceChildren';
%     surf(ax3, 0:big_delta_x:(i_length_total - 1) * big_delta_x, 0:big_delta_z:(k_length_total - 1) * big_delta_z, u_y);
%     zlim([u_y_min u_y_max]);
%     title('Propagation of a 2-D SH-wave with a spherical wavefront without an ABC');
%     title('Propagation of a 2-D SH-wave with a plane wavefront without an ABC');
%     title('Propagation of a 2-D SH-wave with a spherical wavefront with the CPML ABC');
%     title('Propagation of a 2-D SH-wave with a plane wavefront with the CPML ABC');
%     xlabel('x (m)');
%     ylabel('z (m)');
%     zlabel('displacement (m)');
%     grid on;
%     drawnow;
%     writeVideo(v, getframe(figure(2)));
         
end

% close(v);

% plot sensor values

figure
ax4 = gca;
plot(ax4, (0:no_time_steps - 1) * delta_t, 20*log10(abs(sensor_1./max(abs(sensor_1)))));
ylim([-140, 0]);
xlabel('time (s)')
ylabel('20*log_{10}(v_y/max(v_y))')
grid on;

%==============================running time===============================%

toc(running_time);