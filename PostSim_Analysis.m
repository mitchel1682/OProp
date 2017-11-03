clear all; clc; close all;

%Post Sim Analysis Script
%M. McDonald
% load verif_data.mat
load time_vec;

%file read params:
fileID = fopen('statevector.txt','r');
formatSpec = '%f %f %f %f %f %f %f';
sizeA = [7 Inf];
statevector = fscanf(fileID,formatSpec,sizeA)';

%time_vec = [2017 10 29 0 0 0];
Recef_i = dcmeci2ecef('IAU-2000/2006',time_vec);

tout = statevector(:,1);
r_ECI = [statevector(:,2) statevector(:,3) statevector(:,4)];
v_ECI = [statevector(:,5) statevector(:,6) statevector(:,7)];

[size_r size_c] = size(r_ECI);
r_ECEF = zeros(size_r,3);

%time_vec = [2017 10 29 0 0 0];
for i = 1:size_r
    t = datetime(time_vec) + seconds(tout(i,1));
    
    Recef_i = dcmeci2ecef('IAU-2000/2006',datevec(t));
    r_ECEF(i,:) = (Recef_i*r_ECI(i,:)')';
end


%output plots:
plot(tout,r_ECEF);
legend('r_x','r_y','r_z');
xlabel('time, seconds');
ylabel('meters');

% figure(2); 
% subplot(1,3,1);
% plot(tout,r_ECEF(:,1));
% hold on;
% plot(time,r_ECEF_A(:,1));
% xlabel('time, seconds');
% ylabel('position, meters');
% legend('C analysis','NESC verif_A');
% hold off;
% 
% subplot(1,3,2);
% plot(tout,r_ECEF(:,2));
% hold on;
% plot(time,r_ECEF_A(:,2));
% xlabel('time, seconds');
% ylabel('position, meters');
% legend('C analysis','NESC verif_A');
% hold off;
% 
% subplot(1,3,3);
% plot(tout,r_ECEF(:,3));
% hold on;
% plot(time,r_ECEF_A(:,3));
% xlabel('time, seconds');
% ylabel('position, meters');
% legend('C analysis','NESC verif_A');
% hold off;