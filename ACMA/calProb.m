function P = calProb(guanxi,time, para)

P = zeros(1, length(time));
for e = 1:length(time)

    index = guanxi(e);

    j = floor((index-1)/para.W) + 1;
    i = index - (j-1)*para.W; 



    if para.P_type(i,j) == 1
        x = time(e)-(para.TF(i,j)+para.LEN(i,j));
        P(e) = roundn(para.P_a(i,j)*(x^2)+para.P_c(i,j),-2);


    else
        x = time(e)-para.TF(i,j);
        P(e) = roundn(para.P_a(i,j)*(x^2)+para.P_c(i,j),-2);

    end
end

end

