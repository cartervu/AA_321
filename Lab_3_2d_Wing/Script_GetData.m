%% Lab 3 2D Wing
clear all; close all; clc;

absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
% Abs changes for each user that pulls from the github
rel = "2d_Wing/";
qfolder = ["q20"; "q25"; "q30"; "q35"; "q40"; "q45"];
alpha = ["-4", "-2", "+0", "+2", "+4", "+6", "+8", "+10", "+12", "+14", "+16", "+18", "+20"];

MMatrix = cell(length(qfolder),length(alpha));
Aggregate = cell(length(qfolder),length(alpha));

section = "D" %'A' 'B' 'C' 'D' 'E' 



for qindex = 1:length(qfolder)
    q = qfolder(qindex);
    for aindex = 1:length(alpha)
        a = alpha(aindex);
        % This gives us the dataset for a specific q/alpha combination
        folderpath = strcat(absp,rel,q,"/","alpha",a);
        
%         disp(q)
%         disp(a)
        cd(folderpath)
%         disp('--------')
%         disp(pwd)
%         disp(ls)
        filelist = dir;
        [m,n] = size(filelist);
        Datamatrix = []; % get a Datamatrix for each q/alpha combination
        for i = 3:m
            if filelist(i).name(2) == section % check the section (AA, AB, AC, etc.)
                disp(filelist(i).name)
                Datamatrix = [Datamatrix; readmatrix(filelist(i).name)];
            end
        end
        % Datamatrix(:,1) zeros
        % Datamatrix(:,2) port number
        % Datamatrix(:,3) pressure reading at port (psf)
        % Datamatrix(:,4) q ind for wind tunnel at time of reading (psf)
        % Datamatrix(:,5) temperature (ambient, deg C)
        
        % If two values for the same qind, a, have the same port number, I
        % want to combine the measurement into a single reading for temp and 
        % dynamic pressure with a statistical uncertainty?
        % Yes, if we assume the q_ind is constant, so the only difference
        % between measurements is port number, pressure reading, and
        % temperature reading. we average pressure and temperature
        % readings.
        % error on the q seems to be +- 0.15 psf max.
        
        
        % aggregate the data 
        AggMatrix = [];
        [sized,y] = size(Datamatrix);
        for portnum = 1:43 %if include0, start from 0 instead of 1
            rows = [];
            rownum = 0;
            % Get rows that match with a given port number by checking each
            % row's port number
            for dataindex = 1:sized
                if Datamatrix(dataindex,2) == portnum
                    rownum = rownum + 1;
                    rows(rownum,:)=Datamatrix(dataindex,:);
                end
            end
            if size(rows) ~= 0 % if we had any rows that matched the port number, add a row to the final array
                static = stat_digital(rows(:,3));
                dynamic = stat_digital(rows(:,4));
                temp = stat_digital(rows(:,6));
                AggMatrix(portnum,:) = [portnum, static, dynamic, temp]; % if include0, replace portnum with portnum+1 in first occurence
            end
        end
%         disp(q)
%         disp(qindex)
%         disp(a)
%         disp(aindex)
        MMatrix{qindex,aindex} = Datamatrix; 
        Aggregate{qindex,aindex} = AggMatrix;
    end
end


cd(absp)

Mname = strcat("2D_Wing_Data",section,".mat")
Aname = strcat("2D_Wing_Data_Aggregate",section,".mat")

save(Mname,"MMatrix")
save(Aname,"Aggregate")


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