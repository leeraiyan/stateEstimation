function TorF = isload(gen)
%ISLOAD  Checks for dispatchable loads.
%   TorF = isload(gen) returns a column vector of 1's and 0's. The 1's
%   correspond to rows of the gen matrix which represent dispatchable loads.
%   The current test is Pmin < 0 AND Pmax == 0.
%   This may need to be revised to allow sensible specification
%   of both elastic demand and pumped storage units.



[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

TorF = gen(:, PMIN) < 0 & gen(:, PMAX) == 0;

return;
