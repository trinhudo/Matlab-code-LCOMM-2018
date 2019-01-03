function [D_H1k, D_H2l, D_G2p, D_G1q, D_f2_pq_kl, ...
    D_f1ql, positionDownlinkUsers, ...
    positionUplinkUsers] ...
    = CreateLargeScaleFading(K, L, P, Q, Parameters)
%%
RadiusOfCell        = Parameters(1);
InnerZoneRadius     = Parameters(2);
RadiusOfNearestUser = Parameters(3);
StandardDeviation   = 10^(Parameters(4)/10);
ploss               = Parameters(5);

%% Ramdom deployment

% large-fading for K downlink users - inner zone
Zvector = StandardDeviation*randn(1,K);
rvector ...
    = RadiusOfNearestUser*ones(1,K) ...
    + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,K);
anglevector = 2*pi*rand(1,K);

positionDownlinkUsers ...
    = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
distanceBS_To_DownlinkUsers_Zone1 ...
    = sqrt(sum(positionDownlinkUsers'.^2));
PL_downlink_Zone1 ...
    = GetPathloss(30.18, 26, distanceBS_To_DownlinkUsers_Zone1);

D_H1k = (diag(PL_downlink_Zone1)).^(0.5);

% large-fading for L downlink users - outer zone
Zvector = StandardDeviation*randn(1,L);
rvector ...
    = InnerZoneRadius*ones(1,L) ...
    + (RadiusOfCell-InnerZoneRadius)*rand(1,L);
anglevector = 2*pi*rand(1,L);

positionDownlinkUsers = [positionDownlinkUsers; ...
    (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
distanceBS_To_DownlinkUsers_Zone2 ...
    = sqrt(sum(positionDownlinkUsers(K+1:K+L)'.^2));
PL_downlink_Zone2 ...
    = GetPathloss(30.18, 26, distanceBS_To_DownlinkUsers_Zone2);

D_H2l = (diag(PL_downlink_Zone2)).^(0.5);

% large-fading for P uplink users - outer zone
Zvector = StandardDeviation*randn(1,P);
rvector ...
    = InnerZoneRadius*ones(1,P) ...
    + (RadiusOfCell-InnerZoneRadius)*rand(1,P);
anglevector = 2*pi*rand(1,P);

positionUplinkUsers ...
    = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
distanceBS_To_UplinkUsers_Zone2 = sqrt(sum(positionUplinkUsers'.^2));
PL_uplink_Zone2 = GetPathloss(30.18, 26, distanceBS_To_UplinkUsers_Zone2);

D_G2p = (diag(PL_uplink_Zone2)).^(0.5);

% large-fading for Q uplink users - inner zone
Zvector = StandardDeviation*randn(1,K);
rvector ...
    = RadiusOfNearestUser*ones(1,Q) ...
    + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,Q);
anglevector = 2*pi*rand(1,Q);

positionUplinkUsers = [positionUplinkUsers; ...
    (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
distanceBS_To_UplinkUsers_Zone1 ...
    = sqrt(sum(positionUplinkUsers(P+1:P+Q)'.^2));
PL_uplink_Zone1 ...
    = GetPathloss(30.18, 26, distanceBS_To_UplinkUsers_Zone1);

D_G1q = (diag(PL_uplink_Zone1)).^(0.5);

% large-fading for CCI P ULUs to K DLUs
positionDownUsers_PQtimes ...
    = kron(positionDownlinkUsers(1:K+L,:),ones(P+Q,1));
positionUpUsers_KLtimes ...
    = repmat(positionUplinkUsers(1:P+Q,:), K+L, 1);
distance_PQup2KLdown ...
    = sqrt(sum((positionDownUsers_PQtimes-positionUpUsers_KLtimes).^2 ,2));
PL_PQup2KLdown = GetPathloss(145.4, 37.5, distance_PQup2KLdown/1000);

D_f2_pq_kl = (reshape(PL_PQup2KLdown, P+Q, K+L)).^(0.5);

% large-fading for CCI Q ULUs to L DLUs
D_f1ql = D_f2_pq_kl(P+1:P+Q,K+1:K+L);

end

function [Pathloss] = GetPathloss(firstPar, secondPar, distance)
% distance (km)
Pathloss_dB = firstPar + secondPar*log10(distance);
Pathloss = 10.^(-Pathloss_dB/10);
end