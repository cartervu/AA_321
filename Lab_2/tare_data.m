%% Tare data
% The correction is should be subtracted from the curves in order to tare
% the data. For example, for every value of tare_QIndpsf, we should
% subtract the value tare_qcorrect.
%% The naming convention of the variables is (type of test)_(measured dimension)_(correction or variance)
%% first the average zero wind value is calculated, then the variance of the data selection is calculated
%% then corrected curve is calculated.

load('windtunnel_data.mat')

%% tare 
tare_qcorrect = mean(tare_QIndpsf(1:81));
tare_qcorrect_var = sqrt(var(tare_QIndpsf(1:81)));
tare_q = tare_QIndpsf-tare_qcorrect;    % (psf) the corrected dynamic pressure 

tare_Axialforcecorrect = mean(tare_AxialForceLb(1:81));
tare_Axialforcecorrect_var = sqrt(var(tare_AxialForceLb(1:81)));
tare_F = tare_AxialForceLb-tare_Axialforcecorrect;  %(lb) the corrected Drag

%% disk
disk_qcorrect = mean(disk_QIndpsf(1:128));
disk_qcorrect_var = sqrt(var(disk_QIndpsf(1:128)));
disk_q = disk_QIndpsf - disk_qcorrect;

disk_Axialforcecorrect = mean(disk_AxialForceLb(1:128));
disk_Axialforcecorrect_var = sqrt(var(disk_AxialForceLb(1:128)));
disk_F = disk_AxialForceLb - disk_Axialforcecorrect;

disk_DeltaPpsfcorrect = mean(disk_DeltaPpsf(1:128));
disk_DeltaPpsfcorrect_var = sqrt(var(disk_DeltaPpsf(1:128)));
disk_deltaP = disk_DeltaPpsf - disk_DeltaPpsfcorrect;   %(psf) The corrected difference between the front and back pressures

%% sphere, no ring
notrip_qcorrect = mean(notrip_QIndpsf(1:70));
notrip_qcorrect_var = sqrt(var(notrip_QIndpsf(1:70)));
notrip_q = notrip_QIndpsf - notrip_qcorrect;

notrip_Axialforcecorrect = mean(notrip_AxialForceLb(1:70));
notrip_Axialforcecorrect_var = sqrt(var(notrip_AxialForceLb(1:70)));
notrip_F = notrip_AxialForceLb - notrip_Axialforcecorrect ;

notrip_DeltaPpsfcorrect = mean(notrip_DeltaPpsf(1:70));
notrip_DeltaPpsfcorrect_var = sqrt(var(notrip_DeltaPpsf(1:70)));
notrip_deltaP = notrip_DeltaPpsf - notrip_DeltaPpsfcorrect;

%% sphere, small ring
smalltrip_qcorrect = mean(smalltrip_QIndpsf(1:78));
smalltrip_qcorrect_var = sqrt(var(smalltrip_QIndpsf(1:78)));
smalltrip_q = smalltrip_QIndpsf - smalltrip_qcorrect;

smalltrip_Axialforcecorrect = mean(smalltrip_AxialForceLb(1:78));
smalltrip_Axialforcecorrect_var = sqrt(var(smalltrip_AxialForceLb(1:78)));
smalltrip_F = smalltrip_AxialForceLb - smalltrip_Axialforcecorrect ;

smalltrip_DeltaPpsfcorrect = mean(smalltrip_DeltaPpsf(1:78));
smalltrip_DeltaPpsfcorrect_var = sqrt(var(smalltrip_DeltaPpsf(1:78)));
smalltrip_deltaP = smalltrip_DeltaPpsf - smalltrip_DeltaPpsfcorrect;

%% sphere, large ring
largetrip_qcorrect = mean(largetrip_QIndpsf(1:92));
largetrip_qcorrect_var = sqrt(var(largetrip_QIndpsf(1:92)));
largetrip_q = largetrip_QIndpsf - largetrip_qcorrect;

largetrip_Axialforcecorrect = mean(largetrip_AxialForceLb(1:92));
largetrip_Axialforcecorrect_var = sqrt(var(largetrip_AxialForceLb(1:92)));
largetrip_F = largetrip_AxialForceLb - largetrip_Axialforcecorrect ;

largetrip_DeltaPpsfcorrect = mean(largetrip_DeltaPpsf(1:92));
largetrip_DeltaPpsfcorrect_var = sqrt(var(largetrip_DeltaPpsf(1:92)));
largetrip_deltaP = largetrip_DeltaPpsf - largetrip_DeltaPpsfcorrect;

%% Find the reynolds number curve for all of the q curves
q2re = sqrt(2*0.00232662)*0.416667 /(0.3778E-6); % sqrt(2*rho)*d/mu

tare_Re = q2re * sqrt(tare_q);
disk_Re = q2re * sqrt(disk_q);
notrip_Re = q2re * sqrt(notrip_q);
smalltrip_Re = q2re * sqrt(smalltrip_q);
largetrip_Re = q2re * sqrt(largetrip_q);

%% Find Cd from the axial force and the reynolds number
DRe2cd = 8 * 0.0023267 / pi / (0.3778E-6)^2; % 8*rho/pi/mu^2

tare_cd = DRe2cd*( -tare_F./(tare_Re.^2));
disk_cd = DRe2cd*( -disk_F./(disk_Re.^2));
notrip_cd = DRe2cd*( -notrip_F./(notrip_Re.^2));
smalltrip_cd = DRe2cd*( -smalltrip_F./(smalltrip_Re.^2));
largetrip_cd = DRe2cd*( -largetrip_F./(largetrip_Re.^2));

%% Find Cp as a function of the reynolds number
close all
Re2cp = 2*0.416667 ^2*0.0023267 / (0.3778E-6 )^2;

disk_cp = Re2cp*( disk_deltaP./(disk_Re.^2));
notrip_cp = Re2cp*( notrip_deltaP./(notrip_Re.^2));
smalltrip_cp = Re2cp*( smalltrip_deltaP./(smalltrip_Re.^2));
largetrip_cp = Re2cp*( largetrip_deltaP./(largetrip_Re.^2));


% figure(2) % 
% loglog(tare_Re(135:510), tare_cd(135:510))% tare drag vs Re plot
% 
% figure(3) % 495522
% loglog(disk_Re(134:490), disk_cd(134:490))
% 
% figure(4) % 305037
% loglog(notrip_Re(81:462), notrip_cd(81:462))
% 
% figure(5) % 61148
% loglog(smalltrip_Re(88:525), smalltrip_cd(88:525))
% 
% figure(6) % 399000
% loglog(largetrip_Re(104:777), largetrip_cd(104:777))

figure(7)
loglog(disk_Re(134:490), disk_cd(134:490),'linewidth',2)
hold on
loglog(notrip_Re(81:462), notrip_cd(81:462),'--','linewidth',2)
hold on
loglog(smalltrip_Re(88:525), smalltrip_cd(88:525),':','linewidth',2)
hold on
loglog(largetrip_Re(104:777), largetrip_cd(104:777),'-.','linewidth',2)
hold on
loglog(tare_Re(135:512), tare_cd(135:512),':','linewidth',1)
xlabel('R_E')
ylabel('C_D')
legend('disk','sphere','sphere w/smallring','sphere w/largering','tare')


% figure(8)
% loglog(disk_Re(134:490), disk_cp(134:490))
% 
% figure(9)
% loglog(notrip_Re(81:462), notrip_cp(81:462))
% 
% figure(10)
% loglog(smalltrip_Re(88:525), smalltrip_cp(88:525))
% 
% figure(11)
% loglog(largetrip_Re(104:777), largetrip_cp(104:777))

figure(12)
loglog(disk_Re(134:490), disk_cp(134:490),'linewidth',2)
hold on
loglog(notrip_Re(81:462), notrip_cp(81:462),'--','linewidth',2)
hold on
loglog(smalltrip_Re(88:525), smalltrip_cp(88:525),':','linewidth',2)
hold on
loglog(largetrip_Re(104:777), largetrip_cp(104:777),'-.','linewidth',2)

xlabel('R_E')
ylabel('C_{P,f-r}')
legend('disk','sphere','sphere w/smallring','sphere w/largering')


