
function P = calProb2(i,j,time,para)


 
    if para.P_type(i,j) == 1
        x = time-(para.TF(i,j)+para.LEN(i,j));
        P = para.P_a(i,j)*(x^2)+para.P_c(i,j);


    elseif para.P_type(i,j) == 2
        x = time-para.TF(i,j);
        P = para.P_a(i,j)*(x^2)+para.P_c(i,j);

    else
        x = time-(para.TF(i,j)+para.LEN(i,j)/2);
        P = para.P_a(i,j)*(x^2)+para.P_c(i,j);
    end


end

