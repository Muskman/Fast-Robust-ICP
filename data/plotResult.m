%% Plot output from the cpp code

clc
clear

pcFileName = 'happy_vrip';
pcFolder = strcat(pcFileName,'/res');

for i = 0:2
    resFile = strcat(pcFolder,'/m',num2str(i),'output.txt');
    res = textread(resFile);
    figure
    plot(res(:,1),res(:,3),'LineWidth', 3);
    title(strcat('M',num2str(i),' Result'));
    xlabel('Time'); ylabel('Error');
    set(gca, 'YScale', 'log')
end