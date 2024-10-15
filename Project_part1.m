%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using fast fourier transport (fft) method
%   In matlab, there is a internal function fft(...) which can do the job
%   For more detials on fft, type "help fft" in matlab command window and
%   press enter button to check the matlab help file  
% Unit of energy: hbar*omega, where h_bar is the Planck constant and
%   omega is the frequency of the trap
%   Unit of length: l=sqrt(h_bar/(m*omega)), where sqrt(...) is the square
%  mo
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
T=5*pi;                         % Time duration of the evolution
M = 10^3;                     % Total No. of steps in the evolution
dt = T/M;                       % Time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters to vary 
A = 1;      %Driving amplitude
omega = 5;  %Driving Frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UV = exp(-1i*(X.^2/2)*dt/2);    % One-step propagator in position space, only taking diagonal form
UV = exp(-1i*(X.^2/2+0.1*X.^4)*dt/2);
UT = exp(-1i*(P.^2/2)*dt);       % One-setp propagator in momentum space
% note, hbar=1 in our dimensionless units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
X0=4.0;
sigma=1.0;  % sigma is the width of the initial wavepacket
%psiprep=exp(-(X(1:N-1)-X0).^2/0.5)  squeezed
psiprep=exp(-(X(1:N)-X0).^2/(2*sigma^2));  %Gaussian state
psi=psiprep/sqrt(sum(abs(psiprep).^2));%normalized state
%plot (X(1:N),abs(psi(1:N)).^2);   % plotting initial state
hold on 
%plot (P(1:N),abs(phi(1:N)).^2) 
psi_0=psi;

%Creating figure and initializing gif file
figure;
set(gcf,'position',[500 500 1000 500])
h1 = animatedline;
h1.Color = [1 0 1];
h2 = animatedline;
h2.Color = [0 0 1];
axis([a b 0 0.2]);

% Add a legend to label both lines
legend({'Perturbed (magenta)', 'Unperturbed (blue)'}, 'Location', 'northeast');

% Create a dynamic filename based on A and omega
gifFile = sprintf('projectpart1_A%.2f_omega%.2f.gif', A, omega);

% Initial wavefunction for both perturbed and unperturbed systems
psi_pert = psi_0;    % Perturbed system starts with psi_0
psi_unpert = psi_0;  % Unperturbed system starts with psi_0

for m = 1:M
    % Time-dependent perturbation V(t) = A * sin(x) * cos(omega * t)
    f_t = cos(omega* m * dt); %time dependent perturbation
    V_pert = A * sin(X) * f_t; %total perturbation
    UV_pert = exp(-1i * (X.^2/2) + V_pert * dt /2);

   % Evolve wavefunctions for the perturbed system
    psi_pert = UV_pert .* psi_pert;
    phi_pert = fft(psi_pert);
    phi_pert = UT .* phi_pert;
    psi_pert = ifft(phi_pert);
    psi_pert = UV_pert .* psi_pert;

    % Evolve wavefunctions for the unperturbed system
    psi_unpert = UV .* psi_unpert;
    phi_unpert = fft(psi_unpert);
    phi_unpert = UT .* phi_unpert;
    psi_unpert = ifft(phi_unpert);
    psi_unpert = UV .* psi_unpert;

    clearpoints(h1);
    clearpoints(h2);
    addpoints(h1, X(1:N), abs(psi_pert(1:N)).^2);  % Perturbed
    addpoints(h2, X(1:N), abs(psi_unpert(1:N)).^2);  % Unperturbed

    % Update the title with the values of A and omega
    title(sprintf('With (magenta) and Without (blue) Perturbation | A = %.2f, \\omega = %.2f', A, omega));
    
    % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if m == 1
        imwrite(imind, cm, gifFile, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end

    pause(0.05);
    drawnow;
end
%psi=psi_0; %final state updated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot (X(1:N),abs(psi(1:N)).^2)  %plotting the final state profile