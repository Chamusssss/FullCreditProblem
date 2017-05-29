clear all;
N = 2500; % Number of creditors
NZ = 5000; % Number of samples from MoG (pi*) 
nE = 3; % Number of epsilion samples to take PER z sample
NPi = 600; % Number of samples from MCMC of pi
NRuns = 1; % Number of times to recompute integral before averaging results
S = 20; % Dimension of Z
k = 2; % Number of Gaussians in MoG
burninRatio = 0.1;
C = 4;

a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM] = ProblemParams(N, S, true);
    %Sample from pi
    disp('BEGIN MCMC SAMPLING FROM PI')
    t = cputime;
    B = floor(NPi * burninRatio);
    f = @(z) DensityAtZ(z,H,BETA,tail,EAD,LGC);
    sampleZ = slicesample(rand(1,S), NPi, 'pdf', f, 'thin', 3, 'burnin', B);
    %sampleZ = rand(NZ,S); %FOR TESTING PURPOSES
    sampleE = randn(N,nE*NZ);
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
    clear MoGIntegrand;
    clear MoGMu;
    clear MoGSigma;
    clear MoGWeights; 
    clear model;
    disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING Y')
    t = cputime;
    denom = (1-sum(BETA.^2,2)).^(1/2); %Not used as denom but keeping notation consistant
    BZ = BETA*sampleZ;
    BZ = reshape(BZ,N,1,NZ);
    BZ = repelem(BZ,1,1,nE);
    sampleE = reshape(sampleE,N,1,nE*NZ);
    Y = BZ + bsxfun(@times,sampleE,denom);
    clear sampleE;
    clear BETA;
    clear denom;
    clear BZ;
    disp(strcat('FINISH COMPUTING Y...',num2str(cputime - t),'s'))
    
    disp('BEGIN COMPUTING INDICATORS')
    t = cputime;
    CH = H;
    CHZE = repmat(CH,1,1,nE*NZ);
    isOne = ((Y <= CHZE) == 1);
    ind = isOne & (cumsum(isOne,2) == 1);
    disp(strcat('FINISH COMPUTING INDICATORS...',num2str(cputime - t),'s'))

    disp('BEGIN COMPUTING LOSS')
    t = cputime;
    weights = EAD.*LGC;
    LossMat = repelem(weights,1,1,NZ*nE).*ind;
    Loss = sum(sum(LossMat,2),1);
    Loss = reshape(Loss,1,nE*NZ);
    LRZ = repelem(arrayfun(@(i) mvnpdf(sampleZ(:,i))/MoGDen(i),1:NZ),1,nE);
    LR = LRZ;
    l = double(Loss > tail).*LR;
    a(r) = vpa(mean(l));
    v(r) = vpa(var(l));
    clear C;
    clear CMM;
    clear CN;
    clear denom;
    clear EAD;
    clear f;
    clear H;
    clear ind;
    clear l;
    clear LGC;
    clear Loss;
    clear LossMat;
    clear LR;
    clear LRE;
    clear LRZ;
    clear MoGDen;
    clear psi;
    clear pTheta;
    clear sampleZ;
    clear theta;
    clear weights;
    clear ZDen;
    clear pncz;
    disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))%clear all;
end
%[vpa(a); vpa(v)]'
vpa(a)
vpa(mean(a))
vpa(mean(v))
