function [Archive, localpop, MCR, MF, k, evaluationadd] = shade(Archive, localpop, para, MCR, MF, k, elitesize)

evaluationadd = 0;

lencols = para.W*para.T;

CR  = randn(elitesize,1).*sqrt(0.1) + MCR(randi(end,elitesize,1));
CR  = repmat(max(0,min(1,CR)),1, lencols); 
F   = min(1,trnd(1,elitesize,1).*sqrt(0.1) + MF(randi(end,elitesize,1)));
while any(F<=0)
    F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
end

F = repmat(F,1,lencols);

Site         = rand(size(CR)) < CR; 

[~,rank] = sort([localpop.Cost],'descend'); 
if elitesize == 1
    p = 1;
else
    p = rank(ceil(rand(1,elitesize).*max(2,rand(1,elitesize)*0.2*elitesize))); % 1*elitesize
end

%%%%%%%%%%%%
PP = [localpop
    Archive];

Offspring = localpop;

for n = 1:elitesize
    timechange = 0;
    
    r1 = randi(elitesize);
    r2 = randi(numel(PP));
    while r2==r1
        r2 = randi(numel(PP));
    end
    
    xguanxi = localpop(n).guanxi;
    xtimepos = localpop(n).timepos;
    % pbest
    pbestguanxi = localpop(p(n)).guanxi;
    pbesttimepos = localpop(p(n)).timepos;
    % xr1
    r1guanxi = localpop(r1).guanxi;
    r1timepos = localpop(r1).timepos;
    % xr2
    r2guanxi = PP(r2).guanxi;
    r2timepos = PP(r2).timepos;


    
    newofftimepos = xtimepos;

    remains = zeros(1,lencols);

    
    for e = 1:length(localpop(n).timepos)

        
        guanxi = xguanxi(e);

        timepos = xtimepos(e);

        
        if Site(n,guanxi) == 0
            continue;
        end

        % pbest
        poses = find(pbestguanxi == guanxi);
        if isempty(poses)
            Xpb = timepos;
        else
            pos = poses(randi([1,length(poses)]));
            Xpb = pbesttimepos(pos);
        end

        % xr1
        pos1s = find(r1guanxi == guanxi);
        if isempty(pos1s)
            Xr1 = timepos;
        else
            pos1 = pos1s(randi([1,length(pos1s)]));
            Xr1 = r1timepos(pos1);
        end


        % xr2
        pos2s = find(r2guanxi == guanxi);
        if isempty(pos2s)
            Xr2 = timepos;
        else
            pos2 = pos2s(randi([1,length(pos2s)]));
            Xr2 = r2timepos(pos2);
        end

        
        % F
        f = F(n, guanxi);

        newofftimepos(e) = timepos + f*(Xpb-timepos+Xr1-Xr2);

        
        if newofftimepos(e) < 0
            newofftimepos(e) = 0;
        elseif newofftimepos(e) > 1
            newofftimepos(e) = 1;
        end


        
        if newofftimepos(e)~=timepos
            
            timechange = 1;
        end

        remains(1, guanxi) = remains(1, guanxi) + 1;
    end

    
    if timechange
        Offspring(n).timepos = roundn(newofftimepos,-2);
        Offspring(n) = decode(Offspring(n), para);
        evaluationadd = evaluationadd + 1;
    end
end





% Update the pop and archive
offcost = [Offspring.Cost];
delta   = [localpop.Cost] - offcost;


replace = delta < 0;
localpop(replace) = Offspring(replace);
Archive = [Archive
    localpop(replace)];
Archive = Archive(randperm(end, min(end,para.nPop)));


if any(replace)
    w      = delta(replace)./sum(delta(replace));
    MCR(k) = dot(w,CR(replace,1));
    MF(k)  = dot(w,F(replace,1).^2)./dot(w,F(replace,1));
    k      = mod(k,length(MCR)) + 1;
end

end


