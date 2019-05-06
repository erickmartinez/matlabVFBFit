function [C_t,VFB,time_,depth_um] = FDNP_SiNxDevice1Dirichlet(D,C0,...
        thickness,TempC,voltage,SimulationTime,CSiNx,area_m)
% Finite differences code for solving Sodium Migration-Diffusion with 
% constant source and mass transfer coeffiecients at intermediate
% Parameters:
% -----------
% D: double
%   The Diffusion coefficient (cm2/s)
% C0: double
%   The sodium ion concentration at the Al/SiNx interface (cm-3)
% thicknesses: double
%   The thickness of SiNx layer (um)
% TempC: float
%   The temperature (C)
% voltage: float
%   The stress voltage (V)
% SimulationTime: int, float
%   The simulation time (s)
% CSiNx: double
%   The capacitance of SiNx (F)
% area_m: double
%   The area of the gate (m)    

    C0_glass = C0;
    TempK = TempC + 273.15;
    
    


    % diameter_mm = 1; % The diameter of the gate electrode (mm)
    % % Get the area in m2
    % diameter_m  = diameter_mm * 1e-3;
    % area_m      = pi*(diameter_m/2)^2;

    % Define the layer stack
    % L1: EVA, L2: SiNx, L3: Si

    % Arrhenius D0 for diffusion coefficient
    L1.D0 = 1.30e-11; % Not used for fitting

    % Activation energies for diffusion coefficient
    L1.Ea = 0.45; % eV  % Not used for fitting


    % Get the Diffusion coefficients at current temperature:
    L1.D = D;%ArrheniusVariable(L2.D0,L2.Ea,TempK);

    L1.C0 = 1E10; % The concentration at t=0 for layer 2 (cm-3)

    % Thicknesses in um
    L1.thickness = thickness;

    % Electric field
    L1.EField = voltage/L1.thickness/100;% %MV/cm

    % Boundary conditions

    % Discretization parameters
    % Number of space points per layer
    M = 10000;
    % Number of time steps
    N = 5000;


    % Setup FD for the first layer
    % set mu1 gamma1
    % time_ = logspace(0,log10(SimulationTime*100),N)/100;
    time_ = linspace(0,SimulationTime,N);
    dt = mean(gradient(time_));
    L1.dx  = L1.thickness / (M-1); % in um


    L1.M = round(L1.thickness/L1.dx)+1;

    % Lump all the exponents of the constants at the end
    q_ = 1.60217662;
    k_  = 1.3806485;

    % Change the units of distance of the mu and gamma to um
    % Change D from (cm2/s) to um2/s (factor of 10^8)
    % Change Efield from MV/cm to V/um (factor of 100)

    L1.qEDkT    = q_*L1.EField*L1.D*1e14/(k_*TempK); % 1.60217662E-19*L2.EField*1E2*L2.D*1e8/(1.38064852E-23*TempK)


    % Set the initial concentration
    C1 = L1.C0*ones(L1.M,1)./C0_glass;
    C_t = L1.C0*ones(L1.M,N)./C0_glass;


    % filename = strcat('Simulated_SIMS_',filetag);
    % movieFile  = strcat(filename,'.avi');

    x1 = linspace(0,L1.thickness,L1.M);

    
    depth_um = x1;

    VFB = zeros(N,1);

    
    L1.mu       = (L1.D*1E8).*dt./(2*L1.dx^2);
    L1.gamma    = L1.qEDkT.*dt./(4*L1.dx); 

    dgAe1 = -(L1.mu+L1.gamma);   % The element that repeats mostly over the diagonal 1
    dgAe2 =  (1+2*L1.mu);        % The element that repeats mostly over the diagonal 3
    dgAe3 = -(L1.mu-L1.gamma);   % The element that repeats mostly over the diagonal 3

    diagA1 = dgAe1*ones(1,L1.M-1);
    diagA2 = dgAe2*ones(1,L1.M);
    diagA3 = dgAe3*ones(1,L1.M-1);

    dgBe1 = (L1.mu+L1.gamma);   % The element that repeats mostly over the diagonal 1
    dgBe2 =  (1-2*L1.mu);        % The element that repeats mostly over the diagonal 3
    dgBe3 = (L1.mu-L1.gamma);   % The element that repeats mostly over the diagonal 3

    diagB1 = dgBe1*ones(1,L1.M-1);
    diagB2 = dgBe2*ones(1,L1.M);
    diagB3 = dgBe3*ones(1,L1.M-1);

    % Construct the sparce matrices A1 and B1 from the diagonal elements
    A1_sp = spdiags([[diagA1,0]',diagA2',[0,diagA3]'],[-1,0,1],L1.M,L1.M);
    B1_sp = spdiags([[diagB1,0]',diagB2',[0,diagB3]'],[-1,0,1],L1.M,L1.M);

    d1 = zeros(L1.M,1);
    d1(1) = 2*(L1.mu+L1.gamma);
    d1(L1.M) = 2*(L1.mu-L1.gamma)*L1.C0/C0_glass;
   
    % Get the LU decomposition
    [LL1,U1,P1] = lu(A1_sp);

    % Iterate over time
    for i=1:N 
        % Record C_(t)
        C_t(:,i) = C1(:,1);
        % solve the linear system;
        b1 = (B1_sp*C1+d1);
        y1 = LL1\(P1*b1);
        C1 = abs(U1\y1);
%         C1 = abs(A1_sp\(B1_sp*C1+d1));
        
        qC = q_.*C1*C0_glass;
        VFB(i,1) = -1/((CSiNx/area_m)*L1.thickness)*trapz(x1',x1'.*qC)*1E-19;    
    end
    
%     time_ = downsample(time_,2);
%     VFB   = downsample(VFB,2);
%     C_t = (downsample(C_t',2)');
    C_t = C_t.*C0_glass;
end