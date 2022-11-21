%% Lets create a function to calculate child generation rate(R)

% M = number of sections the fiber length is divided into 
% delta_l = the resolution of the fiber length/ or smallest fiber length
% considered
% S = Dimensionless fitting parameter.
% P = breakage probability


function R = ChildGen(M,delta_l,S,P)
R = zeros(M,M);
% Loop to allocate values to child generation rate Rij:
for i = 1:M
    for j = 1:M
        if j<i
            fr = 0;
        else
            fr = 1;
        end
        R(i,j) = fr*normpdf(i.*delta_l,(j.*delta_l)/2,S.*j*delta_l);

    end
end

% Now we need to normalise the matrix R:
R_sum_column  = sum(R,1);

for i=1:M

    x = (2.*P(i))./R_sum_column(i);
    R(:,i)=x.*R(:,i);    
            
end

end