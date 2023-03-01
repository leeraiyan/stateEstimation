function [V,U,converged, i] = state_est(branch, Ybus, Yf, Yt, Sbus, V0, ref, pv, pq, U)
%   STATE_EST  Solves a state estimation problem.
%   [V, converged, i] = state_est(branch, Ybus, Yf, Yt, Sbus, V0, ref, pv, pq, mpopt)
%   State estimator with WLS algorithm


[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
 TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
 ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%------------------------------------------------------------------------- 
%% Initial setting & options
t0        = clock;
j         = sqrt(-1);
tol       = 1e-8;                       %% Tolerance set 
max_it    = 500;                        %% Maximum number of Newton iteration
flag      = 1;                          %% Activation of  bad data detection
converged = 0;                          %% Flag of Newton iteration
i         = 0;                          %% Number flag of Newton iteration
nb        = length(V0);                 %% Bus number 
nbr       = size(Yf, 1);                %% Branch number 
nref      = [pv;pq];                    %% The number of PV and PQ bus
f         = branch(:, F_BUS);           %% "From" buses
t         = branch(:, T_BUS);           %% "To" buses


%%-------------------------------------------------------------------------
%% Measurement Jacobian Matrix H
%%-------------------------------------------------------------------------


%% The full matrix is provided below with intention to use all possible
%% measurements. 
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V0);
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V0);

H = [
    real(dSf_dVa)   real(dSf_dVm);      %% Real power flow 
    real(dSt_dVa)   real(dSt_dVm);      
    real(dSbus_dVa) real(dSbus_dVm);    %% Real power injection 
    speye(nb)       sparse(nb,nb);
    imag(dSf_dVa)   imag(dSf_dVm);      %% Reactive power flow 
    imag(dSt_dVa)   imag(dSt_dVm);
    imag(dSbus_dVa) imag(dSbus_dVm);    %% Reactive power injection 
    sparse(nb,nb)   speye(nb);
    ];

%%-------------------------------------------------------------------------
%% The full Measurement Set
%%-------------------------------------------------------------------------



%% The measurements from SCADA system are "simulated" with power flow
%% solutions plus noises following Gaussian Distribution. In below the full
%% measurement set with perfect values are given. Full measurement means
%% all the power flow measurements, all bus voltage phasors and all
%% injection power to the buses. 
%% Note: V0 is the voltage phasor calculated from power flow calculation
%%-------------------------------------------------------------------------
z = [
    real(Sf);
    real(St);
    real(Sbus);
    angle(V0);
    imag(Sf);
    imag(St);
    imag(Sbus);
    abs(V0);
    ];

%%-------------------------------------------------------------------------
%% The full Covariance Matrix (R) inverse (Weight Matrix)
%%-------------------------------------------------------------------------



%% The noise of SCADA measurements are assumed to follow Gaussian
%% Distribution with a standard deviation that consists of two components:
%% one is propotioanl to the full scale of the meter being used and the
%% other component being proportional to the actual measurement.
%% Reference:
%% J. J. Allemong, L. Radu and A. M. Sasson, “A fast and reliable state...
%% estimation algorithm for AEP’s new control center,” IEEE Trans....
%% PAS, vol. 101, no. 4, pp. 933-944, 1982.

%% The noise of angle measurements are assumed to be due to GPS
%% synchronization error which is 1 ms. 

%% Measurement devices' full scale
fullscale_vol = 1.2;
fullscale=30;

sigma = [
    0.02 * abs(Sf)      + 0.0052 * fullscale * ones(nbr,1);
    0.02 * abs(St)      + 0.0052 * fullscale * ones(nbr,1);
    0.02 * abs(Sbus)    + 0.0052 * fullscale * ones(nb,1);
    0.018 * pi/180 *ones(nb,1)*3;                                
    0.02 * abs(Sf)      + 0.0052 * fullscale * ones(nbr,1);
    0.02 * abs(St)      + 0.0052 * fullscale * ones(nbr,1);
    0.02 * abs(Sbus)    + 0.0052 * fullscale * ones(nb,1);
    0.02 * abs(V0)      + 0.0052 * fullscale_vol * ones(nb,1) ;
] ./ 3;

ns = length(sigma);

%% Weight Matrix (the inverse of Matrix R)
WInv = spdiags( 1 ./ sigma .^ 2, 0, ns, ns );          

%%  Modeling the noise allowing Gaussian Distribution

err = normrnd(zeros(size(sigma)), sigma ); 

%% Introduce noise to the measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = z + err;                                 %% Introduce errors to the measurements  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%-------------------------------------------------------------------------    
%% Setting initial state (flat start)
%% All voltage magnitude are 1 p.u and all angles are 0.
V = ones(nb,1);

%% Compute estimates based on intial state
Sfe = V(f) .* conj(Yf * V);
Ste = V(t) .* conj(Yt * V);
Sbuse = V .* conj(Ybus * V);

%% The intial estimates corrsponding to the measurements 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_est = [                           %% Construct estimate set corresponding
    real(Sfe);
    real(Ste);
    real(Sbuse);
    angle(V);
    imag(Sfe);
    imag(Ste);
    imag(Sbuse);
    abs(V);                                    %% to the measurement set
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Measurement residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delz = z - z_est;                                                         %% residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normF =  delz' * WInv * delz;                             %% construct J function Equation 3.10
%chusqu = err' * WInv * err; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% Check tolerance before 
fprintf('\n----------------- State Estimation Records -------------------\n');
fprintf('\n it     norm( F )     Chi-threshold     step size');
fprintf('\n----  --------------  --------------  --------------');
fprintf('\n%3d    %10.4e     %10.3e     %10.3e', i, normF,0,0);

if normF < tol                               
   converged = 1;
    fprintf('\nConverged!\n');
end

%%-------------------------------------------------------------------------
%% The measurements filter through which measuremented are chosen.
%% ------------------------------------------------------------------------

%%-------------------------------------------------------------------------
%% ID of all possible measurements
pf=[1:nbr];                                                       %% all Pf 
pt=[nbr+1:2*nbr];                                                 %% all Pt
pbus=[2*nbr+1:2*nbr+nb];                                 %% all P injection
Va=[2*nbr+nb+1:2*nbr+2*nb];               %% all angle measurement from PMU

qf=[2*nbr+2*nb+1:3*nbr+2*nb];                                     %% all Qf
qt=[3*nbr+2*nb+1:4*nbr+2*nb];                                     %% all Qt
qbus=[4*nbr+2*nb+1:4*nbr+3*nb];                          %% all Q injection 
Vm=[4*nbr+3*nb+1:4*nbr+4*nb];                                     %% all Vm  

%%-------------------------------------------------------------------------
%% "vv" vector represent the ID the selected measurements. 
%% Select measurement set real power flow measurement"From", reactive
%% power flow "To" and voltage magnitude measurement for state estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vv=[pf,...                                                        %% all Pf
    qt,...                                                       %% all Qt
    Vm                                                              %% all Vm
   ]';                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ------------------------------------------------------------------------
zz=z(vv);                                       %% selected measurement set
ddelz = delz(vv);               %% Residual of the selected measurement set


%% "ww" is a vector containning the index of state variables that will be
%% updated each iteration. The voltage angle and magnitude of reference bus
%% is unchanged. 

ww = [ nref; nb+nref ]; 

%% the maximum number of data can be rejected
max_it_bad_data =length(vv);
% 
ibd = 1;                              %% Flag for the bad data modification        
nm = length(vv);                         %% Number of the used measurements
   
while ((~converged) & (ibd <= max_it_bad_data))
     

%%-------------------------------------------------------------------------
%% The reduced H, covariance matrix inverse and measurements sets according
%% to selected measurements. 
    HH = H(vv,ww);                 %% Select the relevent J matrix entities
                                           %% based on selected measurement
                                                           %% configuration
                                          
    WWInv = WInv(vv,vv);         %% Select the relevent inverse matrix of R
                                           %% based on selected measurement
                                                           %% configuration
    
    
    
    
    VVa = angle(V(nref));        %% Angles of the none-reference buses that
                                              %% will be updated after each
                                                               %% iteration
                                          
    VVm = abs(V(nref));           %% Magnitudes of the none-reference buses
                                              %% will be updated after each
                                                               %% iteration
    
%%---------------------------------------------------------------------
%% Newton Iteration
   
%% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normF = delz' * WInv * delz;              %% construct J function acorrding to Equation 3.10 with
                                            %% all matrix of reduced sizes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  




        
    %% when the result does not converge and the iteration is smaller than 
    %% the limitation
    while ((~converged) & (i < max_it))
        %% update iteration counter
        i = i + 1;
        
        %% Compute delta x based on equation 2.12 of Ali Abur's book
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F = HH' * WWInv * ddelz;                             %% the right side of equation 3.24 
        G = HH' * WWInv * HH;                              %% the left side of equation 3.24 
        dx = (G \ F);                     %% compute delta x based on equation 3.24
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        %% "VVa" is firt half of the dx vector where the angle of reference
        %% bus is excluded
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VVa = VVa + dx(1:nb-1);                                                 %% update VVa   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %% "VVm" is the second half of solution where the magnitude
        %% of reference bus is excluded
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VVm = VVm + dx(nb:2*nb-2);                                                 %% update VVm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        
        
        %% Update voltage phasor vector "V"       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        V(nref) = VVm .* exp(1j * VVa);  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        Sfe = V(f) .* conj(Yf * V);
        Ste = V(t) .* conj(Yt * V);
        Sbuse = V .* conj(Ybus * V);
        
        %% Compute estimates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_est = [
            real(Sfe);
            real(Ste);
            real(Sbuse);
            angle(V);
            imag(Sfe);
            imag(Ste);
            imag(Sbuse);
            abs(V);
        ];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        
        zz_est=z_est(vv);
        
        %% update residual and objective function J with the new estimates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delz = z - z_est;
        ddelz = delz(vv);                                       %% update residual 
        normF = ddelz' * WWInv * ddelz;                               %% update objective function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
   
        Chi=chi2inv(0.99,(nm-2*nb));
        
        %% Check for convergence
        %% The square of the step size is used intead of absolute value as an indication
        %% for convergence.

step = dx' * dx;                   %(blanks here) The square of the step size

     
        fprintf('\n%3d    %10.3e      %10.3e     %10.3e', i, normF, Chi,  step);

        if (step < tol)                   %% when the step size square is smaller than the threhold
                                          %% an indication of convergence
                                          %% is generated. 
            converged = 1;
            fprintf('\nState estimator converged in %d iterations.\n', i);
        end
    end
    
        if ~converged
            fprintf('\nState estimator did not converge in %d iterations.\n', i);
        end 
        
        
        
        if (normF>=Chi)                    %% (Blank here)
                                           %% Chi-square detection is used here
                                           %% to generate event "Bad data".
            fprintf('\n Attention: The measurement set may contain bad data.\n', i);
        end 

    
        
    %%---------------------------------------------------------------------

 end
se_t = etime(clock, t0);

sprintf('\t-------------------------------------------------\');
fprintf('\n State Estimation converge in %d seconds:\n',se_t);
sprintf('\t-------------------------------------------------\n');
