clear
clc
totalTimes = 15;
totalInstances = 256; 

ii = [1:1:256];
parfor nn = 1:length(ii)
    ins = ii(nn);
    main1(ins,totalTimes);
end