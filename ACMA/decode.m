function [solnew] = decode(sol, para)
W = para.W;
T = para.T;
TE = para.TE;

solnew.guanxi = [];
solnew.time = [];
solnew.mP = [];
solnew.timepos = [];
solnew.Cost = [];



N = zeros(W, 1);
%
lk = zeros(1,T);
assign = ones(W, T);


e = 1;


TF = cell(W, T);
LEN = cell(W, T);
for i = 1:W
    for j = 1:T
        TF{i,j} = para.TF(i,j);
        LEN{i,j} = para.LEN(i,j);
    end
end

origintimepos = sol.timepos;
originguanxi = sol.guanxi;
mTF = [];
mP = []; 
mPosition = [];
timepos = [];

while e <= length(originguanxi)


    index = originguanxi(e);
    target = floor((index-1)/W) + 1;
    weapon = index - (target-1)*W;

    if assign(weapon, target) == 0
        e = e + 1;
        continue;
    end

   
    totallen = sum(LEN{weapon,target});
    thispos = origintimepos(e);
    twrate = LEN{weapon,target}./totallen;
    temptl = 0;
    tempth = twrate(1);
    
    for nn = 1:length(TF{weapon,target})
        if nn == length(TF{weapon,target})
            starttime = roundn(TF{weapon,target}(end) + (thispos-temptl)*totallen,-2);
            break;
        end
        if thispos <= tempth
            starttime = roundn(TF{weapon,target}(nn) + (thispos-temptl)*totallen,-2);
            break;
        else
            temptl = tempth;
            tempth = tempth + twrate(nn+1);
            
        end
        
    end

    mPosition = [mPosition index];
    mTF = [mTF starttime];
    endtime = starttime + TE(weapon, target);
    timepos = [timepos thispos];




   
    N(weapon, 1) = N(weapon, 1) + 1;
    Nflag = 1;

    lk(1,target) = lk(1,target) + 1;
    lkflag = 1;


    if N(weapon, 1) >= para.N(weapon,1)
        assign(weapon,:) = 0;
        Nflag = 0;
        for j = 1:T
            TF{weapon,j} = [];
            LEN{weapon,j} = [];
        end
    end


    if lk(1,target) >= para.lk(1,target)
        assign(:,target) = 0;

        lkflag = 0;
        for i = 1:W
            TF{i,target} = [];
            LEN{i,target} = [];
        end
    end

    
    if Nflag
        Lmin = starttime - para.deta(weapon,1);
        Lmax = starttime + para.deta(weapon,1);

        for t = 1:T
            
            if assign(weapon,t) == 0
                continue;
            end

            thisTF = TF{weapon,t};

           
            flagAdd2 = 0;
            flagDelete2 = [];


            TFadd2 = [];
            LENadd2 = [];


            for nn = 1:length(thisTF)
                t1 = thisTF(nn);
                t2 = t1 + LEN{weapon,t}(nn);

                
                if t1<Lmin && t2>Lmin && t2<=Lmax
                    LEN{weapon,t}(nn) = Lmin-t1;

                    if Lmin-t1<0.1
                        flagDelete2(end+1) = nn;
                    end



                elseif t1>=Lmin && t2<=Lmax 
                    flagDelete2(end+1) = nn;
                    
                elseif t1>=Lmin && t1<Lmax && t2>Lmax
                    LEN{weapon,t}(nn) = t2-Lmax;
                    TF{weapon,t}(nn) = Lmax;
                    

                    if t2-Lmax<0.1
                        flagDelete2(end+1) = nn;
                    end

                    
                elseif t1<Lmin && t2>Lmax

                    
                    if Lmin-t1>0.1 && t2-Lmax>0.1

                        
                        LEN{weapon,t}(nn) = Lmin-t1;


                       
                        flagAdd2 = 1;
                        LENadd2 = [LENadd2, t2-Lmax];
                        TFadd2 = [TFadd2, Lmax];


                        
                    elseif Lmin-t1>0.1 && t2-Lmax<=0.1
                        
                        LEN{weapon,t}(nn) = Lmin-t1;



                       
                    elseif Lmin-t1<=0.1 && t2-Lmax>0.1
                        
                        LEN{weapon,t}(nn) =  t2-Lmax;
                        TF{weapon,t}(nn) =  Lmax;

                    else
                        flagDelete2(end+1) = nn;
                    end

                end

            end
            
            if flagAdd2
                [TF{weapon,t}, shunxu] = sort([TF{weapon,t} TFadd2]);
                LEN{weapon,t} = [LEN{weapon,t} LENadd2];
                LEN{weapon,t} = LEN{weapon,t}(shunxu);
            end
            
            if ~isempty(flagDelete2)

                TF{weapon,t}(flagDelete2) = [];
                LEN{weapon,t}(flagDelete2) = [];

                if isempty(TF{weapon,t})
                    assign(weapon,t) = 0;
                end
            end
        end
    end




    if lkflag
        for w = 1:W

            if assign(w,target) == 0
                continue;
            end
            L = starttime;
            H = endtime;
            thisTF = TF{w,target};

            flagAdd = 0;
            flagDelete = [];
            TFadd = [];
            LENadd = [];



            for nn = 1:length(thisTF)
                t1 = thisTF(nn);
                t2 = t1 + LEN{w,target}(nn);

        
                t3 = t1+TE(w,target);
                t4 = t2+TE(w,target);


               
                if t3<L && t4>L && t4<=H
                    LEN{w,target}(nn) = L-t3;

                    if L-t3<0.1
                        flagDelete(end+1) = nn;
                    end

                elseif t3>=L && t4<=H
                    flagDelete(end+1) = nn;
                elseif t3>=L && t3<H && t4>H
                    
                    if t2<=H
                        flagDelete(end+1) = nn;
                        
                    elseif t2>H
                        LEN{w,target}(nn) = t2-H;
                        TF{w,target}(nn) = H;


                        if t2-H<0.1
                            flagDelete(end+1) = nn;
                        end

                    end
                elseif t3<L && t4>H
                    
                    if t2<=H
                        LEN{w,target}(nn) = L-t3;

                        if L-t3<0.1
                            flagDelete(end+1) = nn;
                        end


                        
                    elseif t2>H

                        
                        if  L-t3>0.1 && t2-H>0.1

                            
                            LEN{w,target}(nn) = L-t3;


                            
                            flagAdd = 1;
                            LENadd = [LENadd, t2-H];
                            TFadd = [TFadd, H];

                            
                        elseif  L-t3>0.1 && t2-H<=0.1
                            
                            LEN{w,target}(nn) = L-t3;


                            
                        elseif L-t3<=0.1 && t2-H>0.1
                            
                            LEN{w,target}(nn) = t2-H;
                            TF{w,target}(nn) = H;


                        else
                            flagDelete(end+1) = nn;
                        end
                    end
                elseif t3>=H
                    if t2<=H
                        flagDelete(end+1) = nn;
                        
                    elseif t2>H && t1<H
                        LEN{w,target}(nn) = t2-H;
                        TF{w,target}(nn) = H;

                        if t2-H<0.1
                            flagDelete(end+1) = nn;
                        end
                    end
                end
            end

            if flagAdd
                [TF{w,target}, shunxu] = sort([TF{w,target} TFadd]);
                LEN{w,target} = [LEN{w,target} LENadd];
                LEN{w,target} = LEN{w,target}(shunxu);

            end
            if ~isempty(flagDelete)
                TF{w,target}(flagDelete) = [];
                LEN{w,target}(flagDelete) = [];
                if isempty(TF{w,target})
                    assign(w,target) = 0;
                end
            end
        end
    end
    e = e + 1;
    if assign == 0
        break;
    end

end
mP = calProb(mPosition, mTF, para);

fit = CalFitness(mPosition, mP, para);


solnew.guanxi = mPosition;
solnew.timepos = timepos;
solnew.time = mTF;
solnew.mP = mP;
solnew.Cost = fit;






end
