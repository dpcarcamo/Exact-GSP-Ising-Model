function [Js, connum] =  Jpairs(J, D)

    
    Jcon = J ~= 0;
    Jcon = triu(Jcon);


    Js = zeros(size(Jcon)); % I have no idea why i cant update Jcon
    count = 0;
    for iter = 1:size(D,1)
        
        i = D(iter,1);
        j = D(iter,2);
        k = D(iter,3);

        count = count + 1;
        Js(i, j) = count;
        if k ~=0 
            count = count + 1;
            Js(i, k) = count; 
        end

    end

    Js = Js + Js.';
    connum = count;

end
