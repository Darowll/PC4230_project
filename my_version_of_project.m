% Project Janani version

% Using Hermite polynomials to initialise the first state and the excited
% state

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
Binsize=25;                       % This is to control every Binsize steps, we take snapshots

%defining the parameters thae question wants us to investigate
A = 1.0;
omega = 1.0;

%defining the propogators
UT = exp(-1i*(P.^2/2)*dt);      % One-step propagator in momentum space

% Using the Hermite polynomials to define the initial state
% As given in code No7, considering the initial state to be a gaussian
% wavepacket with wave vector k0, centered at X0 and width DEL0

K0 = 0;                         % Wavevector of the Gaussian  %to demonstrate FFT-p, use K0=0, for dynamics, K0=3
X0 = 0;                         % Center of the Gaussian 
DEL0 = 1;                       % Width of the Gaussian

% Defining an initial state using hermite polynomials 
Poly_0 = hermiteH(0,X);                                  % 0 for ground state
V_INI_TEMP = Poly_0.*exp(-(X-X0).^2/(2*DEL0^2));
V_INI = V_INI_TEMP / sqrt(V_INI_TEMP * V_INI_TEMP');


% Defining the excited state using Hermite Polynomials
Poly_ex = hermiteH(1,X);                                 % ex for excited state
psi_ex = Poly_ex .* exp(-(X-X0).^2 / 2);                      % First excited state wavefunction
psi_ex = psi_ex / sqrt(sum(abs(psi_ex).^2));              % Normalization


Demofftp=fft(V_INI);                                % FFT of the initial state for momentum space propagation
                                                      
%   Check the probability distribution of the initial state
%area(X,abs(VE_INI).^2);
%   Verify that the initial state is indeed normalized
%sum(abs(VE_INI).^2);


% Using the fft to evolve the initial state

VE_FIN(1:M/Binsize,1:N) = 0;	 % Store the image of intermediate states every T=100
VE_FIN_store = V_INI;           % Store the initial state into a running vector
vv=1;                            % moving speed of the potential,


for m = 1:M
    
    t = m * dt;
    f_t = cos(omega* t);         %time dependent perturbation
    V_pert = A * sin(X) * f_t;   %total perturbation

    UV = (exp(-1i*(X-vv*m*dt).^2/2 + V_pert)*dt/2);     %split operator
    VE_temp_1 = UV.*VE_FIN_store;
    VE_temp_2 = fft(VE_temp_1);
    VE_temp_3 = UT.*VE_temp_2;
    VE_temp_4 = ifft(VE_temp_3);
    VE_temp_5 = UV.*VE_temp_4;
    VE_FIN_store = VE_temp_5;
    
    if (mod(m,Binsize)==0)       %snapshot every Binsize cycles
                                 % Take a picture of the wavepacket and save it in VE_FIN
        VE_FIN(m/Binsize,:) = VE_FIN_store;

       
        % Projecting

        proj_0 = trapz(X, conj(psi_0) .* VE_FIN_store);      %project onto the ground state
        proj_ex = trapz(X, conj(psi_ex) .* VE_FIN_store);    %project onto the ground state

        % Store the transition probabilities 

        prob_0(m) = abs(proj_0)^2;
        prob_ex(m) = abs(proj_ex)^2;
    end
end

%   Make a movie to show the wavepacket evolution
for k = 1:2*M/Binsize
    area(X,abs(VE_FIN(k,:)).^2);    %Plotting the evolution of the wave packet
    hold on;

    if k*Binsize < M
        plot(X,(1+k*vv*Binsize)*X.^2/30000,'r--');  % illustration for the trap only
    else
        plot(X,(1+M*vv)*X.^2/30000,'r--');          % Final Trap position
    end 

    hold off
    axis([a,b,0,0.1])  % x and y axis range for better view
    xlabel('Position');
    ylabel('Probability Density');
    title(sprintf('Wavepacket Evolution at Step %d', k));


    f = getframe;  % animation command, check help page
    im=frame2im(f); % animation command

    pause(0.1);
end
