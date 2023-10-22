


%% Difference from above is the for-loop goes from i = 1:(t2-t1). 
% (also dR=cell(t2-t1,2) is different. instead of i = 1:(t2-t1-1)
%Also, when concatenating the data in line ~108  difference in for-loop: i = 1:(t2-t1). instead of i = 1:(t2-t1-1)
%% Calculating lambda and dR
Lambda = cellfun(@(x) abs(diff(x)), R,'UniformOutput',false); 
Lambda(cellfun(@isempty,Lambda)) = {nan}; 

%calculating dR, j=1 is the lhs wavecrests, and j=2 is rhs.
dR = cell(t2-t1,2);
for i = 1:(t2-t1)
    for j = 1:2
        if numel( R{i,j} ) == numel( R{i+1,j})
            dR{i,j} = abs(R{i+1,j} - R{i,j});
           
        else
            dR{i,j} = [];
        end
    end
end
dR(cellfun(@isempty,dR)) = {nan}; 


%% Concatenate the data

for j=1:2
Lambda_cat = Lambda(:,j);
dR_cat = dR(:,j);

%Plotting the experimental values, and collecting the corresponding lambda
%and dR/dt values
% for i = 1:(t2-t1-1)
for i = 1:(t2-t1)
    if numel(Lambda_cat{i})+1 == numel(dR_cat{i})
    lam_i{i} = [Lambda_cat{i}*f];
    c_exp{i} = [dR_cat{i}/dt*f];

    c_exp_lead{i} = c_exp{i}(2:end);
    c_exp_trail{i} = c_exp{i}(1:end-1);


    T(i,:) = i;
    end
end
t_vec = find(T);


Lam_i = cat(1, lam_i{:}).';%[m]
C_exp_trail = cat(1, c_exp_trail{:}).'; %[m/s]
C_exp_lead = cat(1, c_exp_lead{:}).'; %[m/s]
expdata = [Lam_i',C_exp_trail',C_exp_lead'];

ExDat{j} = expdata; %store experimental data in cellarray to use later


end

% concatenate datacells
EXPDATA = cat(1, ExDat{:})*1e2; %multiply by 1e2 to get units in [cm] and [cm/s]



