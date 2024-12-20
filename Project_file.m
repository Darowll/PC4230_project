%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using fast fourier transport (fft) method
%   In matlab, there is a internal function fft(...) which can do the job
%   For more detials on fft, type "help fft" in matlab command window and
%   press enter button to check the matlab help file  
% Unit of energy: hbar*omega, where h_bar is the Planck constant and
%   omega is the frequency of the trap
%   Unit of length: l=sqrt(h_bar/(m*omega)), where sqrt(...) is the square
%   root function and m is the mass of the particle
%   Unit of momentum: hbar/l
%    energy unit: hbar\omega,  Hamiltonian --> dimensionless
%%   time dimensionless: omega*t    i d/dt | >= dimension H |>
%    dimensionless time = 2pi. one classical period
%--------------------------------------------------------------------------
a = -20;                       % Left end point 
b = +20;                       % Right end point 
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T= 100*pi;                         % Time duration of the evolution
M = 10^4;                     % Total No. of steps in the evolution
dt = T/M;                       % Time step
H0 = zeros(M);                 % Hamiltonian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters to vary 
A = 0.01;      %Driving amplitude
omega = 1.98;  %Driving Frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV = exp(-1i*(X.^2/2)*dt/2);   % One-step propagator in position space, only taking diagonal form
UT = exp(-1i*(P.^2/2)*dt);       % One-setp propagator in momentum space
% note, hbar=1 in our dimensionless units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
X0=0.0;
sigma=1.0;  % sigma is the width of the initial wavepacket
%psiprep=exp(-(X(1:N-1)-X0).^2/0.5)  squeezed
psiprep=exp(-(X(1:N)-X0).^2/(2*sigma^2));  %Gaussian state
psi=psiprep/sqrt(sum(abs(psiprep).^2));%normalized state
%plot (X(1:N),abs(psi(1:N)).^2);   % plotting initial state
hold on 
%plot (P(1:N),abs(phi(1:N)).^2) 
psi_0=psi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the first excited state (Hermite polynomial H1(x))
%psi_excited_prep = hermiteH(1,(X(1:N)-X0)).* psiprep; % First excited state
psi_excited_prep = hermiteH(2,(X(1:N)-X0)).* psiprep; %Second excited state
%psi_excited_prep = hermiteH(4,(X(1:N)-X0)).* psiprep; %Fourth excited state
psi_excited = psi_excited_prep/sqrt(sum(abs(psi_excited_prep).^2)); %Normalise Excited state

% Initialize transition probability storage
P_1 = zeros(1, M); % Transition probability at each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:M %time steps
    % Time-dependent perturbation V(t) = A * sin(x) * cos(omega * t)
    f_t = cos(omega* m * dt); %time dependent perturbation
    V_pert = A * sin(X) * f_t; %total perturbation

    UV = exp(-1i * (X.^2/2 + V_pert) * dt /2);
    
    %Split Operator method
    psi_1 = UV.*psi_0;
    phi_2 = fft(psi_1);   %wavefunction in momentum space
    phi_3 = UT.*phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV.*psi_3;
    psi_0 = psi_4; %prepare a new cycle 

    %Project onto first excited state and compute probability
    %C_1 = trapz(X, conj(psi_excited) .* psi_0); % Overlap with first excited state
    C_1 = sum(conj(psi_excited) .* psi_0) * L/N;
    P_1(m) = abs(C_1)^2;                % Transition probability

end
psi=psi_0; %final state updated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot (X(1:N),abs(psi(1:N)).^2)  %plotting the final state profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the transition probability as a function of time
time = (0:M-1) * dt/pi;
plot(time, P_1, 'LineWidth', 2);
xlabel('Time (\pi)','FontSize', 16);
ylabel('Transition Probability','FontSize', 16);
title('Transition Probability to Fourth Excited State','FontSize', 16);
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
