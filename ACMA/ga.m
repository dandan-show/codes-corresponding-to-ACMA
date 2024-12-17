
function pop = ga(pop, para, crossp)


empty_individual.guanxi = [];
empty_individual.time = [];
empty_individual.mP = [];
empty_individual.timepos = [];
empty_individual.Cost = [];
Offspring = repmat(empty_individual, para.nPop, 1);

for n = 1:para.nPop

    chosens = randi([1, para.nPop],1,2); 
    Offspring(n)  = OperatorGA(pop(chosens(1)),pop(chosens(2)),crossp, para);
end

pop = [pop
    Offspring];
newcosts = [pop.Cost];
[~,rank]   = sort(newcosts,'descend');
pop = pop(rank(1:para.nPop));

end