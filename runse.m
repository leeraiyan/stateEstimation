%function [returnvalues] = runese (values) 

function [MVAbase, bus, gen, branch, success, et] = runse(casename, mpopt, fname, solvedcase)
%RUNSE  Runs a state estimator.
%   [baseMVA, bus, gen, branch, success, et] = runse(casename, mpopt,
%   fname, solvedcase)
%   Runs a state estimator after a Newton power flow.

%%-----  initialize  -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
 VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
 TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
 ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
 MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
 QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%%-------------------------------------------------------------------------
%% argument setting 

if nargin < 4
    solvedcase = '';                    %% don't save solved case
    if nargin < 3
        fname = '';                     %% don't print results to a file
        if nargin < 2
            mpopt = mpoption;           %% use default options
            if nargin < 1
                casename = 'case14';  %% default data file is 'case30.m'z
            end
        end
    end
end

%%-------------------------------------------------------------------------
%% data reading & convertion

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch] = loadcase(casename);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on
gbus = gen(on, GEN_BUS);                %% what buses are they at

%%-------------------------------------------------------------------------
%% Power Flow
t0 = clock;                             %% time record 

%% initial state with a flat start   
V0  = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
V0(gbus) = gen(on, VG) ./ abs(V0(gbus)).* V0(gbus);
    
%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
    
%% compute complex bus power injections (generation - load)
Sbus = makeSbus(baseMVA, bus, gen);
    
%% run the Newton Ralphson power flow calculation
% alg = mpopt(1);
[V, success, iterations] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
   
%% update data matrices with solution
[bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);

et = etime(clock, t0);

%%-------------------------------------------------------------------------
%% Create perfect measurements based on Load Flow Calculation

%% save some values from load flow solution
Pflf=branch(:,PF);                           %% Real Power Flow From Bus
Qflf=branch(:,QF);                           %% Reactive Power Flow From Bus
Ptlf=branch(:,PT);                           %% Real Power Flow to Bus 
Qtlf=branch(:,QT);                           %% Reactive Power FLow to Bus
Sbuslf = V .* conj(Ybus * V);                %% Complex Power Injection at Bus
Vlf=V;                                       %% Volatge Phasor 

%%-------------------------------------------------------------------------

%% WLS state estimator
k=0;
U=0;
while k<1
[V,U, converged, i] = state_est(branch, Ybus, Yf, Yt, Sbuslf, Vlf, ref, pv, pq,U);
k=k+1;
end
VV=V;


%% update data matrices to match estimator solution ...
%% ... bus injections at PQ buses
Sbus = V .* conj(Ybus * V);
bus(pq, PD) = -real(Sbus(pq)) * baseMVA;
bus(pq, QD) = -imag(Sbus(pq)) * baseMVA;

%% ... gen outputs at PV buses
on = find(gen(:, GEN_STATUS) > 0);           %% "on" generators
gbus = gen(on, GEN_BUS);                     %%  the bus "on" generators are at
refgen = find(gbus == ref);                  %%  the reference gen(s)
gen(on, PG) = real(Sbus(gbus)) * baseMVA + bus(gbus, PD);   %% inj P + local Pd

%% ... line flows, reference bus injections, etc.
[bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
%%-------------------------------------------------------------------------


%% plot differences from load flow solution
Pfe=branch(:,PF);
Qfe=branch(:,QF);
Pte=branch(:,PT);
Qte=branch(:,QT);
nbr = length(Pfe);

subplot(3,2,1), plot(1/pi*(angle(Vlf)-angle(V)),'.'), title('Voltage Angle (p.u.)'), xlabel('Bus Number');
subplot(3,2,2), plot(abs(Vlf)-abs(V),'.'), title('Voltage Magnitude (p.u.)'), xlabel('Bus Number');
subplot(3,2,3), plot([1:nbr],(Pfe-Pflf)/100,'r.',[1:nbr],(Pte-Ptlf)/100,'b.'), title('Real Flow (p.u.)'), xlabel('Line Number');
subplot(3,2,4), plot([1:nbr],(Qfe-Qflf)/100,'r.',[1:nbr],(Qte-Qtlf)/100,'b.'), title('Reactive Flow (p.u.)'), xlabel('Line Number');
subplot(3,2,5), plot(baseMVA*real(Sbuslf-Sbus)/100, '.'), title('Real Injection (p.u.)'), xlabel('Bus Number');
subplot(3,2,6), plot(baseMVA*imag(Sbuslf-Sbus)/100, '.'), title('Reactive Injection (p.u.)'), xlabel('Bus Number');


Pff=sum(abs(Pfe-Pflf))
qff=sum(abs(Qfe-Qflf))
Ptt=sum(abs(Pte-Ptlf))
Qtt=sum(abs(Qte-Qtlf))

PP=sum(abs(baseMVA*real(Sbuslf-Sbus)))
QQ=sum(abs(baseMVA*imag(Sbuslf-Sbus)))
%%-------------------------------------------------------------------------  
%% end state estimator code
%printpf(baseMVA, bus, gen, branch, [], success, et, 1, mpopt);

