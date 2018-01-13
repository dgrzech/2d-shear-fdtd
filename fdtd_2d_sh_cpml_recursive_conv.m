% This is code for simulation of propagation of 2-D SH-waves with CPML ABC
% employing recursive convolution method. It was written by Daniel Grzech
% as part of MSc Individual Project supervised by Dr Panagiotis Kosmas.

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
analysis = -1/4 * omega_p^2 * ((0:no_time_steps - 1) * delta_t - delay).^2;
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

b_x_1_v_y = zeros(x_cpml_1, 1);
b_x_2_v_y = zeros(x_cpml_2, 1);
c_x_1_v_y = zeros(x_cpml_1, 1);
c_x_2_v_y = zeros(x_cpml_2, 1);
b_z_1_v_y = zeros(z_cpml_1, 1);
b_z_2_v_y = zeros(z_cpml_2, 1);
c_z_1_v_y = zeros(z_cpml_1, 1);
c_z_2_v_y = zeros(z_cpml_2, 1);

kappa_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
kappa_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
sigma_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
sigma_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
alpha_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
alpha_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
b_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
b_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);
c_x_1_tau_yx = zeros(x_cpml_1 - 1, 1);
c_x_2_tau_yx = zeros(x_cpml_2 - 1, 1);

kappa_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
kappa_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
sigma_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
sigma_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
alpha_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
alpha_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
b_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
b_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);
c_z_1_tau_yz = zeros(z_cpml_1 - 1, 1);
c_z_2_tau_yz = zeros(z_cpml_2 - 1, 1);

% x-axis--left

for i = 1:x_cpml_1
    
    kappa_x_1_v_y(i) = 1.0 + kappa_x_max * ((x_cpml_1 - i) / (x_cpml_1 - 1.0))^n_1;
    sigma_x_1_v_y(i) = sigma_x_1_max * ((x_cpml_1 - i) / (x_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_x_1_v_y(i) = alpha_x_max * ((i - 1.0) / (x_cpml_1 - 1.0))^(n_3);
    b_x_1_v_y(i) = exp(-(sigma_x_1_v_y(i) / kappa_x_1_v_y(i) + alpha_x_1_v_y(i)) * delta_t);
    
    if((sigma_x_1_v_y(i) == 0.0) && (alpha_x_1_v_y(i) == 0.0) && (i == x_cpml_1))
        
        c_x_1_v_y(i) = 0;
        
    else
        
        c_x_1_v_y(i) = sigma_x_1_v_y(i) * (b_x_1_v_y(i) - 1.0) / (sigma_x_1_v_y(i) * kappa_x_1_v_y(i) + kappa_x_1_v_y(i)^2 * alpha_x_1_v_y(i)); 
        
    end
    
end

for i = 1:x_cpml_1 - 1
   
    kappa_x_1_tau_yx(i) = 1.0 + kappa_x_max * ((x_cpml_1 - i - 0.5) / (x_cpml_1 - 1.0))^n_1;
    sigma_x_1_tau_yx(i) = sigma_x_1_max * ((x_cpml_1 - i - 0.5) / (x_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_x_1_tau_yx(i) = alpha_x_max * ((i - 0.5) / (x_cpml_1 - 1.0))^(n_3);
    b_x_1_tau_yx(i) = exp(-(sigma_x_1_tau_yx(i) / kappa_x_1_tau_yx(i) + alpha_x_1_tau_yx(i)) * delta_t);
    c_x_1_tau_yx(i) = sigma_x_1_tau_yx(i) * (b_x_1_tau_yx(i) - 1.0) / (sigma_x_1_tau_yx(i) * kappa_x_1_tau_yx(i) + kappa_x_1_tau_yx(i)^2 * alpha_x_1_tau_yx(i)); 
    
end

% x-axis--right

for i = 1:x_cpml_2
    
    kappa_x_2_v_y(i) = 1.0 + kappa_x_max * ((x_cpml_2 - i) / (x_cpml_2 - 1))^n_1;
    sigma_x_2_v_y(i) = sigma_x_2_max * ((x_cpml_2 - i) / (x_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_x_2_v_y(i) = alpha_x_max * ((i - 1.0) / (x_cpml_2 - 1.0))^(n_3);
    b_x_2_v_y(i) = exp(-(sigma_x_2_v_y(i) / kappa_x_2_v_y(i) + alpha_x_2_v_y(i)) * delta_t);
    
    if((sigma_x_2_v_y(i) == 0.0) && (alpha_x_2_v_y(i) == 0.0) && (i == x_cpml_2))
        
        c_x_2_v_y(i) = 0.0;
        
    else
        
        c_x_2_v_y(i) = sigma_x_2_v_y(i) * (b_x_2_v_y(i) - 1.0) / (sigma_x_2_v_y(i) * kappa_x_2_v_y(i) + kappa_x_2_v_y(i)^2 * alpha_x_2_v_y(i)); 

    end
    
end

for i = 1:x_cpml_2 - 1
   
    kappa_x_2_tau_yx(i) = 1.0 + kappa_x_max * ((x_cpml_2 - i - 0.5) / (x_cpml_2 - 1.0))^n_1;
    sigma_x_2_tau_yx(i) = sigma_x_2_max * ((x_cpml_2 - i - 0.5) / (x_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_x_2_tau_yx(i) = alpha_x_max * ((i - 0.5) / (x_cpml_2 - 1.0))^(n_3);
    b_x_2_tau_yx(i) = exp(-(sigma_x_2_tau_yx(i) / kappa_x_2_tau_yx(i) + alpha_x_1_tau_yx(i)) * delta_t);
    c_x_2_tau_yx(i) = sigma_x_2_tau_yx(i) * (b_x_2_tau_yx(i) - 1.0) / (sigma_x_2_tau_yx(i) * kappa_x_2_tau_yx(i) + kappa_x_2_tau_yx(i)^2 * alpha_x_2_tau_yx(i)); 
    
end

% z-axis--bottom

for k = 1:z_cpml_1
    
    kappa_z_1_v_y(k) = 1.0 + kappa_z_max * ((z_cpml_1 - k) / (z_cpml_1 - 1.0))^n_1;
    sigma_z_1_v_y(k) = sigma_z_1_max * ((z_cpml_1 - k) / (z_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_z_1_v_y(k) = alpha_z_max * ((k - 1.0) / (z_cpml_1 - 1.0))^(n_3);
    b_z_1_v_y(k) = exp(-(sigma_z_1_v_y(k) / kappa_z_1_v_y(k) + alpha_z_1_v_y(k)) * delta_t);
    
    if((sigma_z_1_v_y(k) == 0.0) && (alpha_z_1_v_y(k) == 0.0) && (k == z_cpml_1))
        
        c_z_1_v_y(k) = 0.0;
        
    else
        
        c_z_1_v_y(k) = sigma_z_1_v_y(k) * (b_z_1_v_y(k) - 1.0) / (sigma_z_1_v_y(k) * kappa_z_1_v_y(k) + kappa_z_1_v_y(k)^2 * alpha_z_1_v_y(k)); 

    end
    
end

for k = 1:z_cpml_1 - 1
   
    kappa_z_1_tau_yz(k) = 1.0 + kappa_z_max * ((z_cpml_1 - k - 0.5) / (z_cpml_1 - 1.0))^n_1;
    sigma_z_1_tau_yz(k) = sigma_z_1_max * ((z_cpml_1 - k - 0.5) / (z_cpml_1 - 1.0))^(n_1 + n_2);
    alpha_z_1_tau_yz(k) = alpha_z_max * ((k - 0.5) / (z_cpml_1 - 1.0))^(n_3);
    b_z_1_tau_yz(k) = exp(-(sigma_z_1_tau_yz(k) / kappa_z_1_tau_yz(k) + alpha_z_1_tau_yz(k)) * delta_t);
    c_z_1_tau_yz(k) = sigma_z_1_tau_yz(k) * (b_z_1_tau_yz(k) - 1.0) / (sigma_z_1_tau_yz(k) * kappa_z_1_tau_yz(k) + kappa_z_1_tau_yz(k)^2 * alpha_z_1_tau_yz(k)); 
 
end

% z-axis--top

for k = 1:z_cpml_2
    
    kappa_z_2_v_y(k) = 1.0 + kappa_z_max * ((z_cpml_2 - k) / (z_cpml_2 - 1.0))^n_1;
    sigma_z_2_v_y(k) = sigma_z_2_max * ((z_cpml_2 - k) / (z_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_z_2_v_y(k) = alpha_z_max * ((k - 1.0) / (z_cpml_1 - 1.0))^(n_3);
    b_z_2_v_y(k) = exp(-(sigma_z_2_v_y(k) / kappa_z_2_v_y(k) + alpha_z_2_v_y(k)) * delta_t);
    
    if((sigma_z_2_v_y(k) == 0.0) && (alpha_z_2_v_y(k) == 0.0) && (k == z_cpml_2))
        
        c_z_2_v_y(k) = 0.0;
        
    else
        
        c_z_2_v_y(k) = sigma_z_2_v_y(k) * (b_z_2_v_y(k) - 1.0) / (sigma_z_2_v_y(k) * kappa_z_2_v_y(k) + kappa_z_2_v_y(k)^2 * alpha_z_2_v_y(k)); 

    end 
 
end

for k = 1:z_cpml_2 - 1
   
    kappa_z_2_tau_yz(k) = 1.0 + kappa_z_max * ((z_cpml_2 - k - 0.5) / (z_cpml_2 - 1.0))^n_1;
    sigma_z_2_tau_yz(k) = sigma_z_2_max * ((z_cpml_2 - k - 0.5) / (z_cpml_2 - 1.0))^(n_1 + n_2);
    alpha_z_2_tau_yz(k) = alpha_z_max * ((k - 0.5) / (z_cpml_2 - 1.0))^(n_3);
    b_z_2_tau_yz(k) = exp(-(sigma_z_2_tau_yz(k) / kappa_z_2_tau_yz(k) + alpha_z_2_tau_yz(k)) * delta_t);
    c_z_2_tau_yz(k) = sigma_z_2_tau_yz(k) * (b_z_2_tau_yz(k) - 1.0) / (sigma_z_2_tau_yz(k) * kappa_z_2_tau_yz(k) + kappa_z_2_tau_yz(k)^2 * alpha_z_2_tau_yz(k)); 
    
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

delta_x_v_y = sigma_x_v_y ./ kappa_x_v_y;
beta_x_v_y = delta_x_v_y + alpha_x_v_y;
delta_z_v_y = sigma_z_v_y ./ kappa_z_v_y;
beta_z_v_y = delta_z_v_y + alpha_z_v_y;
delta_x_tau_yx = sigma_x_tau_yx ./ kappa_x_tau_yx;
beta_x_tau_yx = delta_x_tau_yx + alpha_x_tau_yx;
delta_z_tau_yz = sigma_z_tau_yz ./ kappa_z_tau_yz;
beta_z_tau_yz = delta_z_tau_yz + alpha_z_tau_yz;

%========matrices for u_y, v_y, tau_yx, tau_yz, and CPML variables========%

u_y = zeros(i_length_total, k_length_total);
v_y = zeros(i_length_total, k_length_total);
tau_yx = zeros(i_length_total, k_length_total);
tau_yz = zeros(i_length_total, k_length_total);

% CPML
 
psi_x_1_v_y = zeros(x_cpml_1, k_length_total);
psi_x_2_v_y = zeros(x_cpml_2, k_length_total);
psi_z_1_v_y = zeros(i_length_total, z_cpml_1);
psi_z_2_v_y = zeros(i_length_total, z_cpml_2);
 
psi_z_1_tau_yz = zeros(i_length_total, z_cpml_1 - 1);
psi_z_2_tau_yz = zeros(i_length_total, z_cpml_2 - 1);
psi_x_1_tau_yx = zeros(x_cpml_1 - 1, k_length_total);
psi_x_2_tau_yx = zeros(x_cpml_2 - 1, k_length_total);

%======================propagation in time and space======================%

plot_freq = 100;

u_y_min = sum(source_val(source_val < 0)) * delta_t * 1.25; % * 0.1;
u_y_max = sum(source_val(source_val > 0)) * delta_t * 1.25; % * 0.1;

% name = '2d_sh_wave_recursive_conv_spherical.avi';
name = '2d_sh_wave_recursive_conv_plane.avi';
v = VideoWriter(name); % create object to write visualisation to file
open(v);

for n = 1:no_time_steps
            
    % update tau_yx and tau_yz

    i = 1:i_length_total - 1;
    k = 1:k_length_total - 1;
    
    tau_yx(i, k) = tau_yx(i, k) + mu .* 1 ./ kappa_x_tau_yx(i, k) * delta_t .* (v_y(i + 1, k) - v_y(i, k)) ./ big_delta_x;
    
    tau_yz(i, k) = tau_yz(i, k) + mu .* 1 ./ kappa_z_tau_yz(i, k) * delta_t .* (v_y(i, k + 1) - v_y(i, k)) ./ big_delta_z;
        
    % update CPML for tau_yx
    
    % left tau_yx
    
    k = 1:k_length_total - 1;
    
    for i = 1:x_cpml_1 - 1
        
        psi_x_1_tau_yx(i, k) = b_x_1_tau_yx(i) .* psi_x_1_tau_yx(i, k) + c_x_1_tau_yx(i) .* (v_y(i + 1, k) - v_y(i, k)) / big_delta_x;
        tau_yx(i, k) = tau_yx(i, k) + delta_t * mu * psi_x_1_tau_yx(i, k);
        
    end
    
    % right tau_yx
    
    ii = x_cpml_2 - 1;
    k = 1:k_length_total - 1;
    
    for i = i_length_total + 1 - x_cpml_2:i_length_total - 1
        
        psi_x_2_tau_yx(ii, k) = b_x_2_tau_yx(ii) .* psi_x_2_tau_yx(ii, k) + c_x_2_tau_yx(ii) .* (v_y(i + 1, k) - v_y(i, k)) / big_delta_x;
        tau_yx(i,k) = tau_yx(i,k) + delta_t * mu * psi_x_2_tau_yx(ii, k);
        ii = ii - 1;
        
    end
    
    % update CPML for tau_yz
    
    % bottom tau_yz
    
    i = 1:i_length_total - 1;
    
    for k = 1:z_cpml_1 - 1

        psi_z_1_tau_yz(i, k) = b_z_1_tau_yz(k) .* psi_z_1_tau_yz(i, k) + c_z_1_tau_yz(k) .* (v_y(i, k + 1) - v_y(i, k)) / big_delta_z;
        tau_yz(i, k) = tau_yz(i, k) + delta_t * mu * psi_z_1_tau_yz(i, k);
       
    end
    
    % top tau_yz
    
    kk = z_cpml_2 - 1;
    i = 1:i_length_total - 1;
    
    for k = k_length_total + 1 - z_cpml_2:k_length_total - 1
        
        psi_z_2_tau_yz(i, kk) =	b_z_2_tau_yz(kk) .* psi_z_2_tau_yz(i, kk) + c_z_2_tau_yz(kk) .* (v_y(i, k + 1) - v_y(i, k)) / big_delta_z;
        tau_yz(i, k) = tau_yz(i, k) + delta_t * mu * psi_z_2_tau_yz(i, kk);
        kk = kk - 1;
        
    end
    
    % update v_y
    
    i = 2:i_length_total - 1;
    k = 2:k_length_total - 1;
    
    v_y(i, k) = v_y(i, k) + delta_t * 1 / rho * (1 ./ kappa_x_v_y(i, k) .* (tau_yx(i, k) - tau_yx(i - 1, k)) ./ big_delta_x + 1 ./ kappa_z_v_y(i, k) .* (tau_yz(i, k) - tau_yz(i, k - 1)) ./ big_delta_z);
        
%     % update CPML for v_y
%     
%     % left v_y
%     
%     k = 2:k_length_total - 1;
%         
%     for i = 2:x_cpml_1
% 
%         psi_x_1_v_y(i, k) = b_x_1_v_y(i) .* psi_x_1_v_y(i, k) + c_x_1_v_y(i) .* (tau_yx(i,k) - tau_yx(i - 1, k)) ./ big_delta_z;
%         v_y(i, k) = v_y(i, k) + delta_t * 1 / rho * psi_x_1_v_y(i, k);
%         
%     end
%      
%     % right v_y
%     
%     ii = x_cpml_2;
%     k = 2:k_length_total - 1;
%     
%     for i = i_length_total + 1 - x_cpml_2:i_length_total - 1
%     
%         psi_x_2_v_y(ii, k) = b_x_2_v_y(ii) .* psi_x_2_v_y(ii, k) + c_x_2_v_y(ii) .* (tau_yx(i, k) - tau_yx(i - 1, k)) / big_delta_x;
%         v_y(i, k) = v_y(i, k) + delta_t * 1 / rho * psi_x_2_v_y(ii, k);
%         ii = ii - 1;
%         
%     end
%     
%     % bottom v_y
%     
%     i = 2:i_length_total - 1;
%     
%     for k = 2:z_cpml_1
%         
%         psi_z_1_v_y(i, k) = b_z_1_v_y(k) .* psi_z_1_v_y(i, k) + c_z_1_v_y(k) .* (tau_yz(i, k) - tau_yz(i, k - 1)) / big_delta_z;
%         v_y(i, k) = v_y(i, k) + delta_t * 1 / rho * psi_z_1_v_y(i, k);
%         
%     end
%     
%     % top v_y
%     
%     kk = z_cpml_2;
%     i = 2:i_length_total - 1;
%     
%     for k = k_length_total + 1 - z_cpml_2:k_length_total - 1
%         
%         psi_z_2_v_y(i, kk) = b_z_2_v_y(kk) .* psi_z_2_v_y(i, kk) + c_z_2_v_y(kk) .* (tau_yz(i, k) - tau_yz(i, k - 1)) / big_delta_z;
%         v_y(i, k) = v_y(i, k) + delta_t * 1 / rho * psi_z_2_v_y(i, kk);
%         kk = kk - 1;
%     
%     end
        
    % source--spherical wavefront
   
     v_y(source_pos(1), source_pos(2)) = v_y(source_pos(1), source_pos(2)) + source_val(n);
     %n
     v_y(source_pos(1), source_pos(2))
     %x = [n v_y(source_pos(1), source_pos(2))]
         
%     % source--plane wavefront
    
%    v_y(round(i_length_total / 2), :) = v_y(round(i_length_total / 2), :) + source_val(n);
    
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
% %     title('Propagation of a 2-D SH-wave with a spherical wavefront without an ABC');
%     title('Propagation of a 2-D SH-wave with a plane wavefront without an ABC');
% %     title('Propagation of a 2-D SH-wave with a spherical wavefront with the CPML ABC');
% %     title('Propagation of a 2-D SH-wave with a plane wavefront with the CPML ABC');
%     xlabel('x (m)');
%     ylabel('z (m)');
%     zlabel('displacement (m)');
%     grid on;
%     drawnow;
%     writeVideo(v, getframe(figure(2)));
  
end
 
close(v);
 
% % plot sensor values
% 
% figure
% ax4 = gca;
% plot(ax4, (0:no_time_steps - 1) * delta_t, 20*log10(abs(sensor_1./max(abs(sensor_1)))));
% ylim([-140, 0]);
% xlabel('time (s)')
% ylabel('20*log_{10}(v_y/max(v_y))')
% grid on;

%==============================running time===============================%

toc(running_time);