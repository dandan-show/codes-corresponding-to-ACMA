function [newoff,evaluationadd] = mutate(sol,para)
W = para.W;
T = para.T;


evaluationadd = 0;

originguanxi = sol.guanxi;
origintimepos = sol.timepos;
off.guanxi = originguanxi;
off.timepos = origintimepos;

len = length(originguanxi);

for e = 1:len


    rr = rand();
    if rr<=0.25
        index = randi([1,W*T]);
        while index == originguanxi(e)
            index = randi([1,W*T]);
        end

       
        off.guanxi(e) = index;
        
    elseif 0.25<rr && rr<=0.5
        index = randi([1,W*T]);
        while index == originguanxi(e)
            index = randi([1,W*T]);
        end
        off.guanxi(e) = index;

        thistimepos = roundn(unifrnd(0,1),-2);
        while thistimepos == origintimepos(e)
            thistimepos = roundn(unifrnd(0,1),-2);
        end
        off.timepos(e) = thistimepos;
        
    elseif rr>0.75
        thistimepos = roundn(unifrnd(0,1),-2);
        while thistimepos == origintimepos(e)
            thistimepos = roundn(unifrnd(0,1),-2);
        end
        off.timepos(e) = thistimepos;
        
    end
end



newoff = decode2(off.guanxi, off.timepos, para);






