%% Lets create a function to calculate breakage probability (P)
% B = Critical breakage ratio
% Cb = Breakage coefficient
% gamma = shear rate

function P = BreakProb(B,Cb,gamma)

P = zeros(1,length(B));

for i=1:length(B)
    if B(i)<1
        P(i) = 0;
    else
        P(i) = Cb.*gamma.*(1-exp(1-B(i)));
%         P(i) = Cb.*(1-exp(1-B(i)));
    end
end



