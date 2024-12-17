function offspring = OperatorGA(p1,p2,crossp,para)
offspring = p1;
if rand()<=crossp
    allguanxi = [p1.guanxi, p2.guanxi];
    alltimepos = [p1.timepos, p2.timepos];
    sorts = randperm(length(allguanxi));
    allguanxi = allguanxi(sorts);
    alltimepos = alltimepos(sorts);
    offspring = decode2(allguanxi, alltimepos, para);
end
   


end