function [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, Ybus, Yf, Yt, V, ref, pv, pq);
%PFSOLN  Updates bus, gen, branch data structures to match power flow soln.
%   [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, ...
%                                   Ybus, Yf, Yt, V, ref, pv, pq)


%% constants
j = sqrt(-1);
nl = size(branch0, 1);      %% number of lines

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% initialize return values
bus     = bus0;
gen     = gen0;
branch  = branch0;

%%----- update bus voltages -----
bus(:, VM) = abs(V);
bus(:, VA) = angle(V) * 180 / pi;

%%----- update Qg for all gens and Pg for swing bus -----
%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?
refgen = find(gbus == ref);             %% which is(are) the reference gen(s)?

%% compute total injected bus powers
%% This is slow in Matlab 5 ...
% Sg = V(gbus) .* conj(Ybus(gbus, :) * V);
%% ... so we do this instead ...
temp = Ybus.';
Sg = V(gbus) .* conj(temp(:, gbus).' * V);

%% update Qg for all generators
gen(:, QG) = zeros(size(gen, 1), 1);                %% zero out all Qg
gen(on, QG) = imag(Sg) * baseMVA + bus(gbus, QD);   %% inj Q + local Qd
%% ... at this point any buses with more than one generator will have
%% the total Q dispatch for the bus assigned to each generator. This
%% must be split between them. We do it first equally, then in proportion
%% to the reactive range of the generator.

if length(on) > 1
    %% build connection matrix, element i, j is 1 if gen on(i) at bus j is ON
    nb = size(bus, 1);
    ngon = size(on, 1);
    Cg = sparse([1:ngon]', gbus, ones(ngon, 1), ngon, nb);

    %% divide Qg by number of generators at the bus to distribute equally
    ngb = sum(Cg)';         %% nb x 1, number of gens at this bus
    ngg = Cg * sum(Cg)';    %% ngon x 1, number of gens at this gen's bus
    gen(on, QG) = gen(on, QG) ./ ngg;
    
    
    %% divide proportionally
    Cmin = sparse([1:ngon]', gbus, gen(on, QMIN), ngon, nb);
    Cmax = sparse([1:ngon]', gbus, gen(on, QMAX), ngon, nb);
    Qg_tot = Cg' * gen(on, QG);     %% nb x 1 vector of total Qg at each bus
    Qg_min = sum(Cmin)';            %% nb x 1 vector of min total Qg at each bus
    Qg_max = sum(Cmax)';            %% nb x 1 vector of max total Qg at each bus
    ig = find(Cg * Qg_min == Cg * Qg_max);  %% gens at buses with Qg range = 0
    Qg_save = gen(on(ig), QG);
    gen(on, QG) = gen(on, QMIN) + ...
        (Cg * ((Qg_tot - Qg_min)./(Qg_max - Qg_min + eps))) .* ...
            (gen(on, QMAX) - gen(on, QMIN));    %%    ^ avoid div by 0
    gen(on(ig), QG) = Qg_save;
end                                             %% (terms are mult by 0 anyway)

%% update Pg for swing bus
gen(on(refgen(1)), PG) = real(Sg(refgen(1))) * baseMVA + bus(ref, PD);  %% inj P + local Pd
if length(refgen) > 1       %% more than one generator at the ref bus
    %% subtract off what is generated by other gens at this bus
    gen(on(refgen(1)), PG) = gen(on(refgen(1)), PG) - sum(gen(on(refgen(2:length(refgen))), PG));
end

%%----- update/compute branch power flows -----
out = find(branch(:, BR_STATUS) == 0);      %% out-of-service branches
br = find(branch(:, BR_STATUS));            %% in-service branches
Sf = V(branch(br, F_BUS)) .* conj(Yf(br, :) * V) * baseMVA; %% complex power at "from" bus
St = V(branch(br, T_BUS)) .* conj(Yt(br, :) * V) * baseMVA; %% complex power injected at "to" bus
branch(br, [PF, QF, PT, QT]) = [real(Sf) imag(Sf) real(St) imag(St)];
branch(out, [PF, QF, PT, QT]) = zeros(length(out), 4);

return;
