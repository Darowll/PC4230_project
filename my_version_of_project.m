% Project Janani version

% Using Hermite polynomials to initialise the first state and the excited
% state

%Initializing all the parameters that we are to use

a = -80;
b = +80;
L = b-a;
N = 2048;
X = a+L*(0:N-1)/N;                     % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1];        % Dimensionless momentum
T = 20*pi;                                % Time duration of the evolution, can adjust this to see the difference between fast and slow move
M = 10^3;                              % Total No. of steps in the evolution
dt = T/M;                              % Time step
Binsize=25;                            % This is to control every Binsize steps, we take snapshots

%defining the parameters thae question wants us to investigate
A = 0.01;
omega = 0.9;

%defining the propogators
UT = exp(-1i*(P.^2/2)*dt);             % One-step propagator in momentum space

% Using the Hermite polynomials to define the initial state
% As given in code No7, considering the initial state to be a gaussian
% wavepacket with wave vector k0, centered at X0 and width DEL0

K0 = 0;                                % Wavevector of the Gaussian  %to demonstrate FFT-p, use K0=0, for dynamics, K0=3
X0 = 0;                                % Center of the Gaussian 
DEL0 = 1;                              % Width of the Gaussian

% Defining an initial state using hermite polynomials 
Poly_0 = hermiteH(0,X);                                  % 0 for ground state
V_INI_TEMP = Poly_0.*exp(-(X(1:N)-X0).^2/(2*DEL0^2));    % exp(-(X(1:N)-X0).^2/(2*DEL0^2)) is the gaussian
V_INI = V_INI_TEMP / sqrt(sum(abs(V_INI_TEMP).^2));      


% Defining the excited state using Hermite Polynomials
Poly_ex = hermiteH(1,X);                                      % ex for excited state
psi_ex_prep = Poly_ex .* exp(-(X(1:N)-X0).^2 / 2);            % First excited state wavefunction
psi_ex = psi_ex_prep / sqrt(sum(abs(psi_ex_prep).^2));        % Normalization


Demofftp=fft(V_INI);                                          % FFT of the initial state for momentum space propagation
                                                     

VE_FIN(1:M/Binsize,1:N) = 0;	                              % Store the image of intermediate states every T=100
VE_FIN_store = V_INI;                                         % Store the initial state into a running vector
                           
% Using the fft to evolve the initial state

for m = 1:M
    
    t = m * dt;
    f_t = cos(omega* t);                                      %time dependent perturbation
    V_pert = A * sin(X) * f_t;                                %total perturbation

    UV = exp(-1i*(X.^2/2 + V_pert) * dt/2);                   %split operator

    VE_temp_1 = UV.*VE_FIN_store;
    VE_temp_2 = fft(VE_temp_1);
    VE_temp_3 = UT.*VE_temp_2;
    VE_temp_4 = ifft(VE_temp_3);
    VE_temp_5 = UV.*VE_temp_4;
    VE_FIN_store = VE_temp_5;
    
    
    % Projecting  

  
    proj_0 = trapz(X, conj(V_INI) .* VE_FIN_store);
    prob_0(m) = abs(proj_0)^2;

    proj_ex = trapz(X, conj(psi_ex) .* VE_FIN_store);        %project onto the ground state
    prob_ex(m) = abs(proj_ex)^2;
    
end
V_INI = VE_FIN_store;

% Plotting the transition probability as a function of time
time = (0:M-1) * dt;
plot(time, prob_0, 'LineWidth', 2);
plot(time, prob_ex, 'LineWidth', 2);
xlabel('Time');
ylabel('Transition Probability');
title('Transition Probability to First Excited State');
grid on;
