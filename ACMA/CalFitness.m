function [Fitness] = CalFitness(mPosition, mP, para)

p_fail = ones(1,para.T);
for i = 1:length(mPosition)
    ind_target = floor((mPosition(i)-1)/para.W) + 1;
    p_fail(ind_target) = p_fail(ind_target)*(1-mP(i));
end

Fitness = roundn(sum(para.V.*(1 - p_fail)),-3); 


end