function [J, h] = RandomGSPIsing(numSpins)
    % Function to create a randome GSP network defined in an Ising
    % interaction matrix J and external field vector h. 

    J = zeros(numSpins); % Initialize random Spin matrix
    
    h = randn(numSpins,1);
    
    J(1,2) = randn();
%     J(1,3) = randn();
%     J(2,3) = randn();

    Jcons = (J ~= 0 ); % Use this to account for nodes we "dont" add. I.e phantom nodes

    for i = 3:numSpins
        [row, col] = find(Jcons);
        randidx = randi(length(row));

        if randi(2) ~= 3 % Connect two nodes
            J(row(randidx),i) = randn();
            J(col(randidx),i) = randn(); 
            Jcons(row(randidx),i) = 1;
            Jcons(col(randidx),i) = 1; 
        else % Connect to one
            if randi(2) == 2
                J(row(randidx),i) = randn();
                Jcons(row(randidx),i) = 1;
                Jcons(col(randidx),i) = 1; 
            else
                J(col(randidx),i) = randn(); 
                Jcons(row(randidx),i) = 1;
                Jcons(col(randidx),i) = 1; 
            end
        end
    end



    J = J + J.';

end
