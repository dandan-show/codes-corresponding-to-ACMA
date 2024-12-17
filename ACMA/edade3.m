function [convergefit, convergeiter] = edade3(para)


empty_individual.guanxi = [];
empty_individual.time = [];
empty_individual.timepos = [];
empty_individual.mP = [];
empty_individual.Cost = [];

convergeiter = zeros(1,31);
convergefit = zeros(1,31);


pairindexes = cell(1,para.W*para.T);
inds = 1;
startp = 1;
for jj = 1:para.T
    for ii = 1:para.W
        len = para.S(ii,jj);
        pairindexes{inds} = [startp:startp+len-1];
        startp = startp + len;
        inds = inds+1;
    end
end

originv = zeros(1,para.maxnvar);
inds = 1;
startp = 1;
for jj = 1:para.T
    for ii = 1:para.W
        len = para.S(ii,jj);
        originv(startp:startp+len-1) = para.V(jj);
        startp = startp + len;
        inds = inds+1;
    end
end


crossp = 1;




pop = repmat(empty_individual, para.nPop, 1);

pop = initialize(para, pop);

Costs = [pop.Cost];
[bestvalue, bestind] = max(Costs);

fit = bestvalue;
convergeiter(1) = 1;
convergefit(1) = fit;



it = para.nPop;
iter = 1;



MCR = zeros(para.nPop,1) + 0.5;
MF  = zeros(para.nPop,1) + 0.5;
Archive = [];
convergetimes = 1;
dbtimes = 1;
k   = 1;
notimproveiter = 1;
while it <= para.evaluation

    
    pop = ga(pop, para, crossp);
    it = it + para.nPop;


    rates = Costs./bestvalue;
    xiangjiao = zeros(para.nPop,para.nPop);
    thislens = zeros(1,para.nPop);

    for i = 1:para.nPop
        thispair = pop(i).guanxi;
        thislens(i) = length(thispair);
        for j = i+1:para.nPop
            counts = 0;
            compair = pop(j).guanxi;
            for ee = 1:thislens(i)
                element = thispair(ee);
                record = find(compair==element);
                if ~isempty(record)
                    counts = counts + 1;
                    compair(record(1)) = [];
                end
            end

            xiangjiao(i,j) = counts;
            xiangjiao(j,i) = xiangjiao(i,j);
        end
    end
  



    pairsims = sum(xiangjiao,1)./(para.nPop-1); 
    pairsims = pairsims./thislens; 
    choselocal = rand(1,para.nPop)<pairsims.*rates;


    locals = find(choselocal);
    localpop = pop(locals);
    localsize = length(locals);

    if localsize~=0
        [Archive, localpop, MCR, MF,k,evaluationadd] = shade(Archive, localpop, para, MCR, MF, k, localsize);

        pop(locals) = localpop;
        it = it + evaluationadd;
    end



    lastbestvalue = bestvalue;
    Costs = [pop.Cost];
    [bestvalue, ~] = max(Costs);
    if bestvalue > lastbestvalue
        notimproveiter = 1;
    else
        notimproveiter = notimproveiter + 1;
    end





    if notimproveiter > 5
        transp = zeros(para.nPop,para.maxnvar);
        transv = zeros(para.nPop,para.maxnvar);

        for nn = 1:para.nPop
            temptransv = originv;
            temptransp = zeros(1, para.maxnvar);
            temppairindex = pairindexes;
            for ee = 1:length(pop(nn).guanxi)
                varind = pop(nn).guanxi(ee);
                pos = temppairindex{varind}(1);
                temppairindex{varind}(1) = [];
                temptransp(pos) = pop(nn).mP(ee);
            end
            transp(nn,:) = temptransp;
            zhi0 = find(temptransp==0);
            temptransv(zhi0) = 0;
            transv(nn,:) = temptransv;

        end

        meancost = mean(Costs);

       
        features = [transp transv];
        D=pdist2(features,features);
        mmean = mean(D(:));
        epsilon = 0.8*mmean;   
        MinPts = 5;                               

        IDX=DBSCAN(para,D,epsilon,MinPts);         
        maxIDX = max(IDX);
        dbtimes = dbtimes + 1;


        for ll = 1:maxIDX

            thisleiset = find(IDX== ll);
            thisleicosts = [pop(thisleiset).Cost];
            thismax = max(thisleicosts);
            thismin = min(thisleicosts);
            if thismax==thismin
                [~,deleteinds] = sort(thisleicosts);
                popdeleteinds = thisleiset(deleteinds(1:end-1));
            else
                deleteinds = find(thisleicosts<meancost);
                popdeleteinds = thisleiset(deleteinds);
            end
            % restart
            if ~isempty(popdeleteinds)
                for rr = 1:length(popdeleteinds)
                    popind = popdeleteinds(rr);
                    [pop(popind), evaluationadd2] = mutate(pop,para);
                end
                it = it + evaluationadd2;
            end

        end
    end


    Costs = [pop.Cost];
    [bestvalue, bestind] = max(Costs);

    if it>=convergetimes*1000

        fit = bestvalue;
        convergeiter(convergetimes+1) = it;
        convergefit(convergetimes+1) = fit;
        convergetimes = convergetimes + 1;
    end


    iter = iter + 1;

  
end

Costs = [pop.Cost];
[bestvalue, bestind] = max(Costs);

fit = bestvalue;
convergeiter(convergetimes+1) = it;
convergefit(convergetimes+1) = fit;



end

