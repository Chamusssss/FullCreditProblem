function [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, loadFixed)
% N: Number of creditor

tail = 0.20; % cursive l in paper
C = 4; % number of credit states

if(loadFixed)
           filename = strcat(pwd,'/Experiments/S=',num2str(S),'/params.mat');
           load(filename); 
else
    % We will be putting all creditors in credit state C with prob p of
    % moving to D
    p = 0.01*(1 + sin(16*pi*(1:N)/N));
    
    % Credit state matrix N x C. Contains probs of each creditor to move to each
    % credit state for THIS timestep
    CMM = zeros(N,C);
    CMM(:,1) = p;
    CMM(:,2) = 1-p;

    % Homogenerous
    EAD = 0.5 + rand(N, 1);           % exposure of the nth obligor, vector
    EAD = EAD / sum(EAD);             %still vector
    CN = 2 * ones(N, 1);        % initial credit state of the nth obligor
    %all are state 2?
    BETA = 1/sqrt(S)*(-1 + (2)*rand(N,S));    % sensitivity to each individual risk
    % N*S matrix
    %BETA = repelem(0.01,N,S);

    LGC = zeros(N, C);
    LGC(:,1) = floor(5*(1:N)/N).^2';

    cumCMM = cumsum(CMM, 2);
    %H_indc = zeros(C, C);
    H = norminv(cumCMM, 0, 1);
    %H_indc(cumCMM<=0) = -1e100;
    %H_indc(cumCMM>=1) = 1e100;
end

end

