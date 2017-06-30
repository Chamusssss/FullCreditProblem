clear all;
N = 2500; % Number of creditors
NZ = 5000; % Number of samples from MC
nE = 30; % Number of epsilion samples to take PER z sample
S = 10; % Dimension of Z
NRuns = 5; % Number of times to recompute integral before averaging results
  
a = zeros(1,NRuns);
v = zeros(1,NRuns);

for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    %Initialize data
    [H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);

    disp('BEGIN SAMPLING')
    t = cputime;
    sampleZ = randn(S,NZ);
    sampleE = randn(N,nE*NZ);
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
    l = double(Loss > tail);
    a(r) = vpa(mean(l));
    v(r) = vpa(var(l));
end

%[vpa(a); vpa(v)]'
vpa(a)
vpa(mean(a))
%vpa(mean(v))

