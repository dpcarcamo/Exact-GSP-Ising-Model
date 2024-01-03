function [H] = Entropy(mean)
    H = -mean.*log(mean)-(1-mean).*log(1-mean);
end
