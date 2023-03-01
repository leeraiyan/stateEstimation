function Sbus = makeSbus(baseMVA, bus, gen)
%MAKESBUS   Builds the vector of complex bus power injections.
%   Sbus = makeSbus(baseMVA, bus, gen) returns the vector of complex bus
%   power injections, that is, generation minus load. Power is expressed
%   in per unit.


%% constants
j = sqrt(-1);

%% define named indices into bus, gen matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%% form net complex bus power injection vector
nb = size(bus, 1);
ngon = size(on, 1);
Cg = sparse(gbus, [1:ngon]', ones(ngon, 1), nb, ngon);  %% connection matrix
                                                        %% element i, j is 1 if
                                                        %% gen on(j) at bus i is ON
Sbus =  ( Cg * (gen(on, PG) + j * gen(on, QG)) ...  %% power injected by generators
            - (bus(:, PD) + j * bus(:, QD)) ) / ... %% plus power injected by loads
        baseMVA;                                    %% converted to p.u.

return;
