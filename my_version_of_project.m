% Project Janani version

% using Hermite polynomials to initialise the first stae and the excited
% state (?)

%Initializing all the parameters that we are to use

a = -80;
b = +80;
L = b-a;
N = 2048;
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1];   % Dimensionless momentum
T = 8;                            % Time duration of the evolution, can adjust this to see the difference between fast and slow move
M = 10^3;                         % Total No. of steps in the evolution
dt = T/M;                         % Time step
Binsize=25                        % This is to control every Binsize steps, we take snapshots

%defining the parameters thae question wants us to investigate
A = 1.0;
omega = 1.0;

%defining the propogators
%UV = exp(-1i*(X.^2/2)*dt/2);   % One-step propagator in position space, only taking diagonal form
UT = exp(-1i*(P.^2/2)*dt);      % One-step propagator in momentum space

% Using the Hermite polynomials to define the initial state
% As given in code No7, considering the initial state to be a gaussian
% wavepacket with wave vector k0, centered at X0 and width DEL0

K0 = 0;                         % Wavevector of the Gaussian  %to demonstrate FFT-p, use K0=0, for dynamics, K0=3
X0 = 0;                         % Center of the Gaussian 
DEL0 = 1;                       % Width of the Gaussian

% Defining an initial state using hermite polynomials as a ground state
%wavefunctions
Poly=hermiteH(4,X);

VE_INI_temp = Poly.*exp(-(X-X0).^2/(2*DEL0^2));      % 1i means "i"
                                                     % Normalized initial state as ground state or excited of harmonic
                                                     % oscillators
VE_INI = VE_INI_temp/sqrt(VE_INI_temp*VE_INI_temp'); % normalization
Demofftp=fft(VE_INI);                                % FFT of the initial state for momentum space propagation

% Using the fft to evolve the initial state

VE_FIN(1:M/Binsize,1:N) = 0;	 % Store the image of intermediate states every T=100
VE_FIN_store = VE_INI;           % Store the initial state into a running vector
vv=1;                            % moving speed of the potential,

for m = 1:M
    UV = exp(-1i*(X-vv*m*dt).^2/2*dt/2);
    VE_temp_1 = UV.*VE_FIN_store;
    VE_temp_2 = fft(VE_temp_1);
    VE_temp_3 = UT.*VE_temp_2;
    VE_temp_4 = ifft(VE_temp_3);
    VE_temp_5 = UV.*VE_temp_4;
    VE_FIN_store = VE_temp_5;
    
    if (mod(m,Binsize)==0)       %s napshot every Binsize cycles
                                 % Take a picture of the wavepacket and save it in VE_FIN
        VE_FIN(m/Binsize,:) = VE_FIN_store;
        
        toc
    end
end

% Projection of the quantum amplitude onto the eigenstates of the
% unperturbed Harmonic oscillator

n_states = 4;                    % Number of Eigenstates to project onto
amplitudes = zeros(n_states, 1);    % Store projection amplitudes

% creating a loop through the first 4 eigenstates 
for n = 0:n_states-1;
    Poly_n = hermiteH(n, X);                  % Hermite polynomial of degree n
    psi_n = Poly_n .* exp(-X.^2 / 2);          % Multiply by Gaussian
    psi_n = psi_n / sqrt(sum(abs(psi_n).^2));  % Normalize the eigenstate

% Calculate projection amplitude
    amplitudes(n+1) = trapz(X, conj(psi_n) .* VE_FIN_store); % Overlap integral
end

% Display the projection amplitudes
disp('Projection amplitudes onto Hermite polynomial eigenstates:');
disp(amplitudes);

