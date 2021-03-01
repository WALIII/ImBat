function Angle = ProjAngleV(L1, L2, LP)
% Angle between 2 vectors projected into plane
% Angle = ProjAngleV(L1, L2, LP)
% INPUT:
%   L1, L2, LP: Vectors of dimension [T x 3] or [1 x 3].
%               They need not be normalized.
% OUTPUT:
%   Angle: Angle between L1 and L2 projected into the plane orthogonal to LP.
%      In other words: the angle between L1 and L2 when looking along LP.
%      The sign of the angle is positive, if [L1, L2, LP] is a right handed
%      system, so the sign changes if L1 and L2 are interchanged.
%      Range of output: -180 < Angle <= 180:
%        L1 = L2  => Angle = 0
%        L1 ~= L2, L1 <-> L2 => Angle <-> -Angle
%      Discontinuity at antiparallel lines:
%        L1 = -L2 => Angle = 180
%        L1 = -L2 + [0,eps,0] => Angle = -180
% Author: Jan Simon, Heidelberg, (C) 2005-2013 matlab.2010@n-simon.de
% License: Use, copy, modify as you like, don't blame me for troubles

% General algorithm:
%   A = L1 x LP;
%   B = L2 x LP;
%   Angle = ATAN2(Norm(A x B), Dot(A, B)) * sign(Dot(L1 x L2, LP))

% Project L1 and L2 into plane perpendicular to LP:
A = CrossTx3(L1, LP);
B = CrossTx3(L2, LP);

% Dot product between projected lines:
AdotB = A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3);

AxB = [A(:, 2) .* B(:, 3) - A(:, 3) .* B(:, 2), ...
       A(:, 3) .* B(:, 1) - A(:, 1) .* B(:, 3), ...
       A(:, 1) .* B(:, 2) - A(:, 2) .* B(:, 1)];

Angle = atan2(NormTx3(AxB), AdotB);

% Sign(LX) is zero for LX==0, but this would destroy Inf values.
negS        = (DotTx3(AxB, LP) < 0);
Angle(negS) = -Angle(negS);

% Care for infinite values:
% Angle(isinf(AdotB)) = Inf;    % -Inf -> Inf
Angle(~isfinite(AdotB)) = Inf;  % NaN and -Inf -> Inf

% No valid angle for to0 short input vectors::
Small = 1.490116119384766e-008;   % SQRT(EPS)
Angle(sum(L1 .* L1, 2) < Small) = Inf;
Angle(sum(L2 .* L2, 2) < Small) = Inf;
Angle(sum(LP .* LP, 2) < Small) = Inf;

function R = CrossTx3(X, Y)
R = [ X(:,2) .* Y(:,3) - X(:,3) .* Y(:,2), ...
      X(:,3) .* Y(:,1) - X(:,1) .* Y(:,3), ...
      X(:,1) .* Y(:,2) - X(:,2) .* Y(:,1)];

function R = NormTx3(X)
R = sqrt(X(:, 1) .* X(:, 1) + X(:, 2) .* X(:, 2) + X(:, 3) .* X(:, 3));

function R = DotTx3(X, Y)
R = X(:, 1) .* Y(:, 1) + X(:, 2) .* Y(:, 2) + X(:, 3) .* Y(:, 3);