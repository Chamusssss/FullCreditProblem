clear all;
N = 2500; % Number of creditors
NZ = 5600; % Number of samples from MC
nE = 2; % Number of epsilion samples to take PER z sample
NRuns = 1;
S = 20; % Dimension of Z
C = 4; % Number of credit states

%Initialize data
[H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);

disp('BEGIN COMPUTING CREDIT STATE SELECTION MATRIX')
t = cputime;
CSMat = sparse(1:N,CN,ones(N,1));
disp(strcat('FINISH COMPUTING CREDIT STATE SELECTION MATRIX...',num2str(cputime - t),'s'))

disp('BEGIN SAMPLING')
t = cputime;
sampleZ = randn(S,NZ);
disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING PNCZ')
t = cputime;
denom = (1-sum(BETA.^2,2)).^(1/2); 
BZ = BETA*sampleZ;
CH = CSMat*H;
CHZ = repmat(CH,1,1,NZ);
BZ = reshape(BZ,N,1,NZ);
CBZ = repelem(BZ,1,C);
PINV = (CHZ - CBZ) ./ denom;
PHI = normcdf(PINV);
PHI = [zeros(N,1,NZ) PHI];
pncz = diff(PHI,1,2); %column wise diff
clear sampleE;
clear sampleZ;
clear BETA;
clear BZ;
clear CH;
clear CHZ;
clear CBZ;
clear PHI;
disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

disp('BEGIN SAMPLING PNCZ')
t = cputime;
cdf = cumsum(pncz,2);
cdf = repelem(cdf,1,1,nE);
u = rand([N,1,nE*NZ]);
isOne = (cdf >= u) == 1;
ind = isOne & (cumsum(isOne,2) == 1);
disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING LOSS')
t = cputime;
weights = EAD.*LGC;
weights = repelem(weights,1,1,nE*NZ);
LossMat = weights.*ind;
Loss = sum(sum(LossMat,2),1);
Loss = reshape(Loss,1,nE*NZ);
l = double(Loss > tail);
vpa(mean(l))
vpa(var(l))
disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
%clear all;

