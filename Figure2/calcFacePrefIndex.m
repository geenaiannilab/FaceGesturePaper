function [faceExpPref,depthPref] = calcFacePrefIndex(Ris, nExp)

[Rmaxs, ~] = max(Ris,[],2);
[Rmins, ~] = min(Ris,[],2);

faceExpPref = (nExp - (sum(Ris) ./ Rmaxs)) ./ (nExp -1);
depthPref = (Rmaxs - Rmins ) ./ (Rmaxs + Rmins);


end