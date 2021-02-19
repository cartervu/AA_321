%% Plotting Scripts Q1, Q2


%% Pressure vs port number plots
clear all; close all; clc;

absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
% Abs changes for each user that pulls from the github
rel = "2d_Wing/";
qfolder = ["q20"; "q25"; "q30"; "q35"; "q40"; "q45"];
alpha = ["-4", "-2", "+0", "+2", "+4", "+6", "+8", "+10", "+12", "+14", "+16", "+18", "+20"];


cd(absp)
DatasetIndex = "E" %'A' 'B' 'C' 'D' 'E' 
load(strcat("2D_Wing_Data_Aggregate",DatasetIndex,".mat"))

% -4 0 4 8 12 16 20

%% P vs Portnum (one for each q)
close all; clc;
cd(strcat(absp,"Figures/Fig_1"))

for qindex = 1:length(qfolder)
    figure()
    for aindex = 1:length(alpha)
        a = alpha(aindex);
        q = qfolder(qindex);
        disp(a)
        Dataset = Aggregate{qindex,aindex};
        if size(Dataset) ~= [0 0]
%         Dataset(:,2) = -Dataset(:,2);
%         Dataset(:,2) = Dataset(:,2) - min(Dataset(:,2));
%         Dataset(:,2) = Dataset(:,2) - avg(Dataset(:,2));
    %     Dataset(end+1,:) = Dataset(1,:);
    %     errorbar(Dataset(:,1),Dataset(:,2),Dataset(:,3))
            hold on
            plot(Dataset(:,1),Dataset(:,2),'DisplayName',a)
            title(strcat(q))
            xlabel('portnum')
            ylabel('pressure')
            set(gca,'FontSize',20)
        end
    end
    legend
    hold off
    saveas(gcf,strcat('p_vs_portnum_',q,DatasetIndex,'.jpg'))
end

%% P vs portnum (one for each alpha)
close all; clc;


for aindex = 1:length(alpha)
    figure()
    a = alpha(aindex);
    disp(a)
    for qindex = 1:length(qfolder)
        q = qfolder(qindex);
%         disp(q)
        Dataset = Aggregate{qindex,aindex};
        if size(Dataset) ~= [0 0]
        Dataset(:,2) = -Dataset(:,2); % enable
%         Dataset(:,2) = Dataset(:,2) - min(Dataset(:,2));
        Dataset(:,2) = Dataset(:,2) - avg(Dataset(:,2)); % enable

    %     Dataset(end+1,:) = Dataset(1,:);
    %     errorbar(Dataset(:,1),Dataset(:,2),Dataset(:,3))
            hold on
            plot(Dataset(:,1),Dataset(:,2),'DisplayName',q)
            title(strcat("alpha=",a))
            xlabel('portnum')
            ylabel('pressure')
            set(gca,'FontSize',20)
        end
    end
    legend
    hold off
end

%% Plot all plots for a given q separately
close all; clc;

qindex = 1;
for aindex = 1:length(alpha)
    a = alpha(aindex);
    q = qfolder(qindex);
    disp(a)
    Dataset = Aggregate{qindex,aindex};
%     Dataset(:,2) = -Dataset(:,2);
%     Dataset(:,2) = Dataset(:,2) - min(Dataset(:,2));
%     Dataset(end+1,:) = Dataset(1,:);
    figure()
    errorbar(Dataset(:,1),Dataset(:,2),Dataset(:,3))
    hold on
    plot(Dataset(:,1),Dataset(:,2))
    title(strcat(q,' alpha=',a))
    xlabel('portnum')
    ylabel('pressure')
    set(gca,'FontSize',20)
end

%% Cp vs x/c plots
clear all; close all; clc;

absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
% Abs changes for each user that pulls from the github
rel = "2d_Wing/";
qfolder = ["q20"; "q25"; "q30"; "q35"; "q40"; "q45"];
alpha = ["-4", "-2", "+0", "+2", "+4", "+6", "+8", "+10", "+12", "+14", "+16", "+18", "+20"];

cd(absp)
load("portlocations.mat")
port_loc = [1:44; port_loc(:,:)];

% plot(port_loc(2,:),port_loc(3,:))


% we have q_infty for each point
% we have p for each point
% how do we get p_infty for each point
% Bernoulli's equation
% We use p port1 + q port1 = pinfty + qinfty
% we know q port1 is 0, because that is the stagnation point
% so we have p_port1 = pinfty + q_infty
% so we have pinfty = p_port1 - q_infty
% Actually, assuming no-slip condition, the q at each port is 0
% so we have p_port - p_infty = q_infty everywhere.
% But this is clearly not the case since then cp = 1 eerywhere.
% Let's try with the p_infty from the first port and see how it goes


% First, pick a q/alpha combo
% Then get the p and q_infty from the first port measurement
% Then perform the subtraction
% combined uncertainty sqrt(dp^2 + dq^2)

% iterate over all the datasets initially, but at first we just pick one
%               1    2     3        4      5      6
% qfolder = ["q20"; "q25"; "q30"; "q35"; "q40"; "q45"];
%            1     2     3     4    5      6    7      8       9     10     11
% alpha = ["-4", "-2", "+0", "+2", "+4", "+6", "+8", "+10", "+12", "+14", "+16",
% "+18", "+20"];
%  12      13

% Dataset = Aggregate{1,6}; % q20, alpha -4

% p_infty = p_port - q_infty for port 0 or port 1?


DatasetIndex = "E" %'A' 'B' 'C' 'D' 'E' 
load(strcat("2D_Wing_Data_Aggregate",DatasetIndex,".mat"))

cd("Figures/Fig_2")

for qindex = 1:length(qfolder);
    for aindex = 1:length(alpha)
        a = alpha(aindex);
        q = qfolder(qindex);
        disp(a)
        Dataset = Aggregate{qindex,aindex};
        if size(Dataset) ~= [0,0]
            Dataset(:,2) = -Dataset(:,2);
        %     Dataset(:,2) = Dataset(:,2) - min(Dataset(:,2));
            Dataset(:,2) = Dataset(:,2) - avg(Dataset(:,2));

            P_infty = Dataset(1,2) - Dataset(1,4);
            P_infty_uncert = sqrt(Dataset(1,3)^2+Dataset(1,5)^2);

            [rownum,colnum] = size(Dataset);
            for rowindex = 1:rownum
                Diff = Dataset(rowindex,2) - P_infty;
                Dataset(rowindex,8) = Diff / Dataset(rowindex,4); % 8 Cp and 9 Cp uncert
                Dataset(rowindex,9) = abs(Dataset(rowindex,8))*sqrt( ...
                    (sqrt(Dataset(rowindex,3)^2 + P_infty_uncert^2)/Diff)^2 + ...
                    (Dataset(rowindex,5)/Dataset(rowindex,4))^2 );
            end
            Dataset(end+1,:) = Dataset(1,:);

            figure()
            errorbar(port_loc(2,:),Dataset(:,8),Dataset(:,9))
            hold on
            plot(port_loc(2,:),Dataset(:,8))
            title(strcat(q,' alpha=',a))
            xlabel('x/c')
            ylabel('Cp')
            set(gca,'FontSize',20)
            set(gca,'YDir','reverse')
            xlim([0 1])
            saveas(gcf,strcat("Cp_vs_xc_",DatasetIndex,"_",q,"alpha",a,".jpg"))
        end
    end
end

close all
% pressure coefficient at each point
% unit vector in direction perpendicular to surface at each point
% normal/axial pressure at each point
% integrate using trapezoidal approximation
% total normal/axial for each alpha/q combination
% total Lift/Drag for each alpha/q combination



%% Single Dataset Cp plot

Dataset = Aggregate{1,2};
% Dataset(:,2) = -Dataset(:,2);
% Dataset(:,2) = Dataset(:,2) - min(Dataset(:,2));

P_infty = Dataset(1,2) - Dataset(1,4);
P_infty_uncert = sqrt(Dataset(1,3)^2+Dataset(1,5)^2);

[rownum,colnum] = size(Dataset);
for rowindex = 1:rownum
    Diff = Dataset(rowindex,2) - P_infty;
    Dataset(rowindex,8) = Diff / Dataset(rowindex,4); % 8 Cp and 9 Cp uncert
    Dataset(rowindex,9) = Dataset(rowindex,8)*sqrt( ...
        (sqrt(Dataset(rowindex,3)^2 + P_infty_uncert^2)/Diff)^2 + ...
        (Dataset(rowindex,5)/Dataset(rowindex,4))^2 );
end
Dataset(end+1,:) = Dataset(1,:);

figure()
plot(port_loc(2,:),Dataset(:,8))


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