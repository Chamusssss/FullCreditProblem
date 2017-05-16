clear all;
N = 2500; % Number of creditors
NZ = 5000; % Number of samples from MoG (pi*) 
nE = 2; % Number of epsilion samples to take PER z sample
NPi = 600; % Number of samples from MCMC of pi
NRuns = 10; % Number of times to recompute integral before averaging results
S = 10; % Dimension of Z
k = 2; % Number of Gaussians in MoG
burninRatio = 0.2;
C = 4;

%Initialize data
[H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);

disp('BEGIN COMPUTING CREDIT STATE SELECTION MATRIX')
t = cputime;
CSMat = sparse(1:N,CN,ones(N,1));
disp(strcat('FINISH COMPUTING CREDIT STATE SELECTION MATRIX...',num2str(cputime - t),'s'))

%Sample from pi
disp('BEGIN MCMC SAMPLING FROM PI')
t = cputime;
B = floor(NPi * burninRatio);
f = @(z) DensityAtZ(z,H,BETA,tail,EAD,LGC,CSMat);
sampleZ = slicesample(rand(1,S), NPi, 'pdf', f, 'thin', 3, 'burnin', B);
%sampleZ = rand(NZ,S); %FOR TESTING PURPOSES
disp(strcat('FINISH MCMC SAMPLING FROM PI...',num2str(cputime - t),'s'))

disp('BEGIN TRAINING MOG')
t = cputime;
[~, model, ~] = Emgm(sampleZ', k);  
MoGWeights = model.weight;
MoGMu = model.mu;
MoGSigma = model.Sigma;
disp(strcat('FINISH TRAINING MOG...',num2str(cputime - t),'s'))

disp('BEGIN SAMPLING')
t = cputime;
sampleZ = SampleMoG(MoGWeights,MoGMu,MoGSigma,NZ)';
MoGDen = arrayfun(@(i) EvalMoG(MoGWeights,MoGMu,MoGSigma,sampleZ(:,i)),1:NZ);
ZDen = arrayfun(@(i) f(sampleZ(:,i)'),1:NZ);
MoGIntegrand = ZDen./MoGDen;
disp('MoG Estimate')
vpa(mean(MoGIntegrand))
vpa(var(MoGIntegrand))
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
clear BETA;
clear BZ;
clear CH;
clear CHZ;
clear CBZ;
clear PHI;
disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING THETA')
t = cputime;
weights = EAD.*LGC;
[pTheta,theta] = GlassermanPTheta(pncz,weights,tail);
disp(strcat('FINISH COMPUTING THETA...',num2str(cputime - t),'s'))
 
disp('BEGIN SAMPLING PNCZ')
t = cputime;
cdf = cumsum(pTheta,2);
cdf = repelem(cdf,1,1,nE);
u = rand([N,1,nE*NZ]);
isOne = (cdf >= u) == 1;
ind = isOne & (cumsum(isOne,2) == 1);
disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

disp('BEGIN COMPUTING LOSS')
t = cputime;
LossMat = repelem(weights,1,1,NZ*nE).*ind;
Loss = sum(sum(LossMat,2),1);
theta = reshape(theta,[1,1,NZ]);
B = zeros([N C NZ]);
for j=1:NZ
    B(:,:,j) = theta(:,:,j)*weights;
end
psi = sum(log(sum(pncz.*exp(B),2)),1);
LRE = reshape(exp(-repelem(theta,1,1,nE).*Loss + repelem(psi,1,1,nE)),1,nE*NZ,1);
Loss = reshape(Loss,1,nE*NZ);
LRZ = repelem(arrayfun(@(i) mvnpdf(sampleZ(:,i))/MoGDen(i),1:NZ),1,nE);
LR = LRE.*LRZ;
l = double(Loss > tail).*LR;
vpa(mean(l))
vpa(var(l))
disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
%clear all;

