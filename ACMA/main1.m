function [] = main1(ins,totalTimes)

eval(['load instance',num2str(ins)]);

S = zeros(W, T);
for i = 1:W
    for j = 1:T
        a = ceil(Tfind_len(i,j)/deta(i,1));
        b = ceil(Tfind_len(i,j)/Texe_len(i,j));
        temps1 = min(a,b); 
        temps2 = min(temps1, N(i,1)); 
        S(i,j) = min(temps2, lk(1,j));
    end
end
wuqihe = sum(S,2);
min11 = min(wuqihe,N);
min1 = sum(min11);

mubiaohe = sum(S,1);
min22 = min(mubiaohe, lk);
min2 = sum(min22);

nvar = min(min1, min2);
maxnvar = sum(S(:));


nPop = 100; 
evaluation = nPop*100*3;
para = struct('W',W,'T',T,'V',V,'deta',deta,...
    'TF',Tfind_stp,'LEN',Tfind_len, 'TE', Texe_len,...
    'TM',Tmeet_stp,'S',S,'maxnvar',maxnvar,...
    'PL',PL,'PH',PH,'P_a',P_a,'P_c',P_c,'P_type',P_type,'N',N,'lk',lk,...
    'nvar',nvar,'nPop',nPop,'evaluation',evaluation);


aa = zeros(totalTimes,31);
bb = zeros(totalTimes,31);

for times = 1:totalTimes
    [a,b] = main2(para);
    aa(times,1:length(a)) = a;
    bb(times,1:length(b)) = b;

end
disp(ins)
names1 = ['convergefit_ACMA.xlsx'];
writematrix(aa,names1,'Sheet',ins);
names2 = ['convergeiter_ACMA.xlsx'];
writematrix(bb,names2,'Sheet',ins);

end

