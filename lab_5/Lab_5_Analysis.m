% AA 321 Lab 5 Analysis Questions 1, 3, 4

clear all; close all; clc;


% Path information
absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/lab_5/";
% Abs changes for each user that pulls from the github
rel = "KWT_Finite_Wing_Data/";

procd = "PROCESSED_DATA";
rawd = "RAW_DATA";

% Chord: 15"
% 
% Length: 87"
% 
% Endcap diameter: 17"
% 
% Trip-dots 1.5" from leading edge


%% Problem 1
clear all; close all; clc; 

% Path information
absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/lab_5/";
% Abs changes for each user that pulls from the github
rel = "KWT_Finite_Wing_Data/";

procd = "PROCESSED_DATA";
rawd = "RAW_DATA";

% Load the data
cd(strcat(absp,rel))

CorrMatrix = readmatrix("BalanceCalibrationMat.csv");
CorrMatrix(:,1) = [];

PMWT_EO = readmatrix("PMWT_Endcaps_Off.csv")


cd(strcat(absp,rel,rawd))

% readmatrix("RUN_0002C.xls")

RAW = cell(1,20);
PROC = cell(1,20);

for i = 1:9
    RAW{1,i} = readmatrix(strcat("RUN_000",num2str(i),"C.xls"));
end

for i = 10:20
    RAW{1,i} = readmatrix(strcat("RUN_00",num2str(i),"C.xls"));
end

% I only need the raw dataset for 2, I don't need any others.

cd(strcat(absp,rel,procd))
for i = [2,3,5,6,9]
    PROC{1,i} = readmatrix(strcat("run_000",num2str(i),"C.csv"));
end

for i = [10,12,13,19,20]
    PROC{1,i} = readmatrix(strcat("run_00",num2str(i),"C.csv"));
end

DataMatrix = RAW{1,2};
% 1: CODE       2: TP	3: AlphaENC	4: Psi	5: QNOM	6: QA       7: LIFTR	
% 8: DRAGR      9: PMR	10: YMR     11: RMR	12: SFR	13% PRESSTS	14: TEMPTS	
% 15: Fouling	16: Time

% We want the Alpha 3, the QNOM 5, the QA 6, the lift 7, the drag 8, the PM 9
% how do we correct for alpha? Not certain...

% I want a column vector of L, D, M for each data point, and I then
% multiply each one by CorrMatrix to get the corrected L, D, M.

% Then, I take D = D - NDRAG * q, M - M - PMWT

% In this section, we create an array that has a LDM vector for each
% data point for the raw datset 2

% row 1 is L
% row 2 is D
% row 3 is M


s = size(DataMatrix);
LDM = zeros(3, s(1));
LDM(1,:) = transpose(DataMatrix(:,7)); % Lift
LDM(2,:) = transpose(DataMatrix(:,8)); % Drag
LDM(3,:) = transpose(DataMatrix(:,9)); % Moment
LDM = CorrMatrix * LDM;
% Dynamically change matrix size
LDM(2:4,:) = LDM(1:3,:);
LDM(1,:) = transpose(DataMatrix(:,3)); % Add in the angle of attack alpha as row 1



NDrag = 0.1;
LDM(3,:) = LDM(3,:) - transpose(NDrag * DataMatrix(:,6)); % NDrag = 0.1 correction


LDM(4,:) = LDM(4,:) - transpose(PMWT_EO(:,2));
% We approximate the PMWT alpha values as the exact same values as the
% alpha values for the wind tunnel tests, there is error in this



% RAW Run 2 has the last 3 columns all at 21 degrees -- to keep all the rows
% consistent, instead of averaging, we'll just delete the last 3 columns.
% LDM(:,end-2:end) = [];

% We use QA, not QCorr, we do not correct the Q-value, this has error

% disp(LDM)
% 10 degrees is column 8 of LDM matrix
% We take CL = L / q S, CD = D/Q S, CM = M / q S c
q = DataMatrix(8,6); % 8th row of DataMatrix, since DataMatrix is rows, not columns
S = 15*87 / 12^2; % We convert S to square feet, since q is in PSF.
c = 15; % We convert c to square feet, since M is in ft-lbs.
% NO it is not! THE TAS HAVE LIED TO US IT IS IN INCH-POUNDS!


% CL = LDM(2,8) / (q*S);
% CD = LDM(3,8) / (q*S);
% CM = LDM(4,8) / (q*S*c);

CLCDCM = LDM;
CLCDCM(2, : ) = CLCDCM(2, : ) / (q*S);
CLCDCM(3, : ) = CLCDCM(3, : ) / (q*S);
CLCDCM(4, : ) = CLCDCM(4, : ) / (q*S*c);


%Compare corrected LDM matrix to corrected KWT LDM Matrix
CLCDCM = transpose(CLCDCM);
disp(CLCDCM(:,2:4))
disp(PROC{1,2}(:,9:11))
% 9, 10, 11 is CL, CD, CM

cd("/Users/carter/Software/MATLAB/321_Matlab_Scripts/lab_5/KWT_Finite_Wing_Data/Figures")

figure()
plot(CLCDCM(:,1),CLCDCM(:,2),'DisplayName',"CLCorr")
hold on
plot(PROC{1,2}(:,7),PROC{1,2}(:,9),'DisplayName',"CLKWT")
hold off
legend()
plotaxes("\alpha","C_L","")
saveas(gcf,('CL_vs_alpha_comparison.jpg'))

figure()
plot(CLCDCM(:,1),CLCDCM(:,3),'DisplayName',"CDCorr")
hold on
plot(PROC{1,2}(:,7),PROC{1,2}(:,10),'DisplayName',"CDKWT")
hold off
legend()
plotaxes("\alpha","C_D","")
saveas(gcf,('CD_vs_alpha_comparison.jpg'))

figure()
plot(CLCDCM(:,1),CLCDCM(:,4),'DisplayName',"CMCorr")
hold on
plot(PROC{1,2}(:,7),PROC{1,2}(:,11),'DisplayName',"CMKWT")
hold off
legend()
plotaxes("\alpha","C_M","")
saveas(gcf,('CM_vs_alpha_comparison.jpg'))

% Note that this procedurewas completed with the approximation of NDRAG =
% 0.1 and without the NLIFT and NPM corrections, which likely account for
% the primary differences

%%
close all; clc;

AR = 5.8;

% one CL plot, one CD plot, one CM plot, one L/D = CL / CD, one CL vs CD
% plot

% Total of 5 plots, with 4 of them vs alpha.

% Include 10 with, 10 without, 35 with, 35 without, each as a curve in each
% plot

% The runs are all the PROC runs, but which ones do we need? 
% 2, 10 no cap
% 3 is back
% 5, 35 no cap
% 6 is back

% 9, 10 with cap
% 10 back
% 12 35 with cap
% 13 back

% 19, 20 are redo of 2, 3, except EOT instead of DFFV

cd(strcat(absp,rel,"Figures/"))
    
% figure()
% for qindex = 1:length(qfolder)
%     q = qfolder(qindex);
%     plot(str2double(alpha),CLArray(qindex,:),'DisplayName',q)
%     hold on
% end

% CL, CD, CM, L/D = CL / CD, LC vs CD
% 1: Run  2: Test  3: TP      4: QA      5: QC      6: ALPHAENC 7: ALPHAC 
% 8: PSI  9: CLWA 10: CDWA   11: CMWA25 12: RE_MAC 13: RE_FT   14: PMTARE

% columns: 9 vs 7, 10 vs 7, 11 vs 7, 9/10 vs 7

% set2 = PROC{1,2}; % 10 no cap increasing
% set3 = PROC{1,3}; % 10 no cap decreasing
% set5 = PROC{1,5}; % 35 no cap increasing
% set6 = PROC{1,6}; % 35 no cap decreasing
% set9 = PROC{1,9}; % 10 cap increasing
% set10 = PROC{1,10}; % 10 cap decreasing
% set12 = PROC{1,12}; % 35 cap increasing
% set13 = PROC{1,13}; % 35 cap decreasing

names = ["P=10, no cap increasing", "P=10, no cap decreasing", ...
    "P=35, no cap increasing", "P=35, no cap decreasing", ...
    "P=10, cap increasing", "P=10, cap decreasing", ...
    "P=35, cap increasing", "P=35, cap decreasing"];
% sets = [set2, set3, set5, set6, set9, set10, set12, set13];
runindices = [2, 3, 5, 6, 9, 10, 12, 13];


shortnames = ["P=10, no cap", ...
    "P=35, no cap", ...
    "P=10, cap",  ...
    "P=35, cap"];
shortrunindices = [2,5,9,12];

figure()
for plotindex = 1:4
    runindex = shortrunindices(plotindex);
    name = shortnames(plotindex);
    plot(PROC{1,runindex}(:,7),PROC{1,runindex}(:,9),'DisplayName',name)
    hold on 
end
plotaxes("\alpha (deg)","C_L","NACA 23012: C_L vs \alpha")
hold off


figure()
for plotindex = 1:4
    runindex = shortrunindices(plotindex);
    name = shortnames(plotindex);
    plot(PROC{1,runindex}(:,7),PROC{1,runindex}(:,10),'DisplayName',name)
    hold on 
end
plotaxes("\alpha (deg)","C_D","NACA 23012: C_D vs \alpha")
hold off

figure()
for plotindex = 1:4
    runindex = shortrunindices(plotindex);
    name = shortnames(plotindex);
    plot(PROC{1,runindex}(:,7),PROC{1,runindex}(:,11),'DisplayName',name)
    hold on 
end
plotaxes("\alpha (deg)","C_M","NACA 23012: C_M vs \alpha")
hold off

figure()
for plotindex = 1:4
    runindex = shortrunindices(plotindex);
    name = shortnames(plotindex);
    plot(PROC{1,runindex}(:,7),PROC{1,runindex}(:,9) ./ PROC{1,runindex}(:,10),'DisplayName',name)
    hold on 
end
plotaxes("\alpha (deg)","C_L/C_D","NACA 23012: C_L/C_D vs \alpha")
hold off



figure()
for plotindex = 1:4 % 8
    runindex = shortrunindices(plotindex);
    name = shortnames(plotindex);
    plot(PROC{1,runindex}(:,10),PROC{1,runindex}(:,9),'DisplayName',name)
    hold on 
end
plotaxes("C_D","C_L","NACA 23012: C_L vs C_D")
hold off

figure()
for plotindex = 1:4 % 8
    runindex = shortrunindices(plotindex);
    name = shortnames(plotindex);
    plot(PROC{1,runindex}(:,10),PROC{1,runindex}(:,9) .* PROC{1,runindex}(:,9),'DisplayName',name)
    hold on 
end
xvals = 0:0.01:0.15;
plot(xvals,16.93*xvals-0.2109,'DisplayName',"Best Fit CL^2 vs CD, no caps")
plot(xvals,21.38*xvals-0.4147,'DisplayName',"Best Fit CL^2 vs CD, caps")
plotaxes("C_D","C_L^2","NACA 23012: C_L^2 vs C_D")
hold off

%        p1 =       21.38  (20.59, 22.17)
%        p2 =     -0.4147  (-0.4568, -0.3726)

%% Problem 4
close all; clc;

% Path information
absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/lab_5/";
% Abs changes for each user that pulls from the github
rel = "KWT_Finite_Wing_Data/";

procd = "PROCESSED_DATA";
rawd = "RAW_DATA";


PROC = cell(1,20);

cd(strcat(absp,rel,procd))
for i = [2,3,5,6,9]
    PROC{1,i} = readmatrix(strcat("run_000",num2str(i),"C.csv"));
end

for i = [10,12,13,19,20]
    PROC{1,i} = readmatrix(strcat("run_00",num2str(i),"C.csv"));
end

% we use the increasing datasets, we take the average of all of them,
% assuming alpha is the same, we use nominal alpha values, we get a line
% OR: we fit all the data simultaneously. Much better. Get a single column
% vector that is all the CL values, square each value .^2, and a single
% column value of all the alpha values, then fit using curve fitting
% toolbox.

% set2 = PROC{1,2}; % 10 no cap increasing
% set3 = PROC{1,3}; % 10 no cap decreasing
% set5 = PROC{1,5}; % 35 no cap increasing
% set6 = PROC{1,6}; % 35 no cap decreasing
% set9 = PROC{1,9}; % 10 cap increasing
% set10 = PROC{1,10}; % 10 cap decreasing
% set12 = PROC{1,12}; % 35 cap increasing
% set13 = PROC{1,13}; % 35 cap decreasing


CL2vecwith = []; % no cap
CDvecwith = [];
alphawith = [];
for plotindex = [1,3]
    runindex = runindices(plotindex);
    CL2vecwith = [CL2vecwith; PROC{1,runindex}(:,9).^2];
    CDvecwith = [CDvecwith; PROC{1,runindex}(:,10)];
    alphawith = [alphawith; PROC{1,runindex}(:,7)];
end

% for CDindex = 1:len(CDvec)
CDindex = 1;
numelements = length(CDvecwith);
while numelements > CDindex
    if CDvecwith(CDindex) > 0.1 || alphawith(CDindex) > 15
       CDvecwith(CDindex) = [];
       CL2vecwith(CDindex) = [];
       alphawith(CDindex) = [];
       numelements = numelements - 1;
    else
       CDindex = CDindex + 1;
    end
end
if CDvecwith(end) > 0.1 || alphawith(end) > 15
       CDvecwith(end) = [];
       CL2vecwith(end) = [];
       alphawith(CDindex) = [];
end

CL2vec = [];
CDvec = [];
alphavec = [];
for plotindex = [5,7]
    runindex = runindices(plotindex);
    CL2vec = [CL2vec; PROC{1,runindex}(:,9).^2];
    CDvec = [CDvec; PROC{1,runindex}(:,10)];
    alphavec = [alphavec; PROC{1,runindex}(:,7)];
end

% for CDindex = 1:len(CDvec)
CDindex = 1;
numelements = length(CDvec);
disp(CDvec)
while numelements > CDindex
    if CDvec(CDindex) > 0.1 || alphavec(CDindex) > 15
       CDvec(CDindex) = [];
       CL2vec(CDindex) = [];
       numelements = numelements - 1;
    else
       CDindex = CDindex + 1;
    end
end

if CDvec(end) > 0.1 || alphavec(end) > 15
       CDvec(end) = [];
       CL2vec(end) = [];
end
disp(CDvec)







figure()
for plotindex = [1,3] % 8
    runindex = runindices(plotindex);
    name = names(plotindex);
    plot(PROC{1,runindex}(:,10),PROC{1,runindex}(:,9) .* PROC{1,runindex}(:,9),'DisplayName',name)
    hold on 
end
xvals = 0:0.01:0.15;
plot(xvals,16.82*xvals-0.1884,'DisplayName',"Best Fit CL^2 vs CD")

plot(CDvecwith,CL2vecwith,'ro')

plotaxes("C_D","C_L^2","NACA 23012: C_L^2 vs C_D")
hold off




% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =       16.93  (16.53, 17.33)
%        p2 =     -0.2109  (-0.2311, -0.1908)
% 
% Goodness of fit:
%   SSE: 0.01411
%   R-square: 0.9971
%   Adjusted R-square: 0.997
%   RMSE: 0.02532



% e = p1 / pi AR
e = 16.93 / pi / AR % This is 1.1547, this is a non-physical value. 
C_DP = 0.2109 / pi / e / AR
% saveas(gcf,strcat('CL_vs_alpha_',section,'.jpg'))



% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =       21.38  (20.59, 22.17)
%        p2 =     -0.4147  (-0.4568, -0.3726)
% 
% Goodness of fit:
%   SSE: 0.009055
%   R-square: 0.9973
%   Adjusted R-square: 0.997
%   RMSE: 0.03009



% e without
ew = 21.38 / pi / AR
C_DPw = 0.4147 / pi / ew / AR

%% Functions

function nooutput = plotaxes(xaxis,yaxis,plotname)
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','black','DisplayName','axes');  %y-axis
    line(xL, [0 0],'Color','black','DisplayName','');  %x-axis
    title(plotname)
    xlabel(xaxis)
    ylabel(yaxis)
    legend('Location','southeast')
    legend
    set(gca,'FontSize',20)
    hold off
end