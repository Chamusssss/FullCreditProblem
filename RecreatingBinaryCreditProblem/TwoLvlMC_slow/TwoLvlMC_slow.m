clear all;
N = 2500; % Number of creditors
NZ = 300; % Number of samples from MC
nE = 200; % Number of epsilion samples to take PER z sample
S = 10; % Dimension of Z
NRuns = 5; % Number of times to recompute integral before averaging results
  
a = zeros(1,NRuns);
v = zeros(1,NRuns);
  
az = zeros(1,NZ);
vz = zeros(1,NZ);

%Initialize data
[H, BETA, tail, EAD, CN, LGC, CMM, C] = ProblemParams(N, S, true);
for r=1:NRuns
    totalT = cputime;
    disp(strcat('RUN NUMBER',num2str(r)))
    for zIndex=1:NZ
        %disp(num2str(zIndex))

        %disp('BEGIN SAMPLING')
        %t = cputime;
        sampleZ = randn(S,1);
        %disp(strcat('FINISH SAMPLING...',num2str(cputime - t),'s'))

        %disp('BEGIN COMPUTING PNCZ')
        %t = cputime;
        denom = (1-sum(BETA.^2,2)).^(1/2); 
        BZ = BETA*sampleZ;
        CH = H;
        CHZ = repmat(CH,1,1,1);
        BZ = reshape(BZ,N,1,1);
        CBZ = repelem(BZ,1,C);
        PINV = (CHZ - CBZ) ./ denom;
        PHI = normcdf(PINV);
        PHI = [zeros(N,1,1) PHI];
        pncz = diff(PHI,1,2); %column wise diff
%         clear sampleZ;
%         clear BETA;
%         clear BZ;
%         clear CH;
%         clear CHZ;
%         clear CBZ;
%         clear PHI;
%         clear PINV;
%         clear denom;
        %disp(strcat('FINISH COMPUTING PNCZ...',num2str(cputime - t),'s'))

        %disp('BEGIN SAMPLING PNCZ')
        %t = cputime;
        cdf = cumsum(pncz,2);
        cdf = repelem(cdf,1,1,nE);
        u = rand([N,1,nE*1]);
        isOne = (cdf >= u) == 1;
        ind = isOne & (cumsum(isOne,2) == 1);
        %clear isOne;
        %clear u;
        %clear cdf;
        %disp(strcat('FINISH SAMPLING PNCZ...',num2str(cputime - t),'s'))

        %disp('BEGIN COMPUTING LOSS')
        %t = cputime;
        weights = EAD.*LGC;
        weights = repelem(weights,1,1,nE*1);
        LossMat = weights.*ind;
        Loss = sum(sum(LossMat,2),1);
        Loss = reshape(Loss,1,nE*1);
        l = double(Loss > tail);
        az(zIndex) = vpa(mean(l));
        %vz(zIndex) = vpa(var(l));
        %disp(strcat('FINISH COMPUTING LOSS...',num2str(cputime - t),'s'))
        %clear C;
        %clear CMM;
        %clear CN;
        %clear EAD;
        %clear H;
        %clear l;
        %clear LGC;
        %clear Loss;
        %clear LossMat;
        %clear pncz;
        %clear weights;
        
    end
    a(r) = vpa(mean(az));
    %v(r) = vpa(var(vz));
    disp(strcat('TOTAL RUNTIME...',num2str(cputime - totalT),'s'))
end

%[vpa(a); vpa(v)]'
vpa(a)
vpa(mean(a))
%vpa(mean(v))

