%% AA 321 Lab 1: Analsis Task 3
clear all; close all; clc;


Path = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_1/QCalData/Part2/";
Filenames = ["qcal_50psf.csv" "qcal_50psf_lower.csv"];

DataMatrix = [];
for name = Filenames
    DataMatrix = [DataMatrix; readmatrix(strcat(Path,name))];
end

DataMatrix = DataMatrix(:,3:5);

% Col DNE/2: time
% Col 1/3: z-coord from top in inches
% Col 2/4: q_ind, psf, rough estimate of wind speed from pressure ports in inlet
% and contraction sections
% Col 3/5: dynamic pressure q, psf, measured by pitot static tube.

% Considerable uncertainty in all measurements due to stagnation pressure 
% averaging over tube inlet diameter and inexact z-values (all are off by
% 0.5 tube diameters, should be 0.5 tube diameters lower than reported).

% Unknown pitot tube diameter, so this effect cannot be quantitatively
% corrected for

% Plot: q vs z, col 5 vs col 3, with error bars.
% To do this: separate values for each z, get the stats on each, and then
% plot those points

% pitot tube measurement is digital, z-value measurement is also digital

% Statistical uncertainties
% Systematic uncertainties could not be quanitified due to lack of
% preciseinformation regarding the potential inaccuracies in the
% experimental setup (such as the diameter of the pitot tube and the true
% z-value)
% Split out into each separate z-value, make a vector of q-values for each
% z-value

% The z-values that are not precisely at 0.1 increments are due to either
% imprecision in both the position of the pitot tube and the measurement of
% the position of the pitot tube. To get a better sense of the plot, we
% apply some smoothing, taking a single average value for every z that is
% at a 0.1 step, including all values with \pm 0.05 in in the averaging.

% Solution: make a dictionary. Keys are the z-values, values are a vector
% of q-values. Dictionaries can't change size, though, so instead we make
% an array with the first column the key and the next columns the average
% and the uncertainty.

% Iterate over each z-value, construct a vector for it from the DataMatrix

% plot(DataMatrix(:,1),DataMatrix(:,3),'.')



Dict = [];
for step = 0:360
    start = 0; %exclude this value
    dz = -0.1;
    stop = 36;
    z = start + dz*step;
    qvalues = [];
    for row = 1:size(DataMatrix,1)
        if abs(DataMatrix(row,1) - z) < 0.05
            qvalues = [qvalues DataMatrix(row,3)];
        end
    end
    if size(qvalues) == [0 0]
        continue
    else
        Dict = [Dict; z stat_digital(qvalues)];
    end
end
Dict2 = Dict


plot(Dict(:,1),Dict(:,2),'k')
hold on


% Dict2(16:2:321,:)=[];
% Dict2(10:2:end-20,:)=[];
% Dict2(3:3:end,:)=[];
% Dict2(2:2:end,:)=[];
DictSize=size(Dict2);
zerror = zeros(DictSize(1),1)+0.05;

errorbar(Dict2(:,1),Dict2(:,2),Dict2(:,3),Dict2(:,3),zerror,zerror,'b.')



% title('True Dynamic Pressure q vs Vertical Position z')
legend('Data','Measurement Uncertainty')
ylabel('True Dynamic Pressure q (psf)')
xlabel('Vertical Position z (in)')
set(gca,'FontSize',20)
axis([-37 0 0 55])
%% Functions

function average = avg(dataset)
    numvals = length(dataset);
    total = 0;
    for valindex = 1:numvals
       total = total + dataset(valindex);
    end
    average = total/numvals;
end

function uncertainty = uncert_digital(dataset)
    uncertainty = (max(dataset)-min(dataset))/(2*sqrt(3));
end

function statistics = stat_digital(dataset)
    statistics = [avg(dataset) uncert_digital(dataset)];
end

function statistics = stat_analog(dataset)
    statistics = [avg(datset) uncert_digital(dataset)/sqrt(2)];
end