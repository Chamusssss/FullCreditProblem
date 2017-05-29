hold on

NSim = 10;
xMax = 3;
yMax = 0.00001;
trueAns = 0.0000059367146221281147369847920602925;
X = 1:xMax;

TwoLvlMC = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

TwoLvlIS = [ 0, 0, 0, 0, 0, 0, 0.0000051995272514346987973365629698286, 0, 0, 0];

plot([0 X],repmat(trueAns,1,xMax + 1),'r')

scatter(ones(NSim,1),TwoLvlMC,'g')
scatter(2*ones(NSim,1),TwoLvlIS,'b')

legend('True Ans')
%title('S:20, l:0.10, sample num:10600')
xlabel('')
xticks([0 1 2])
xticklabels({'','2LvlMC(6533,13066)','2LvlIS(4600,5000,10000)'})
ylabel('P(L > l)')
axis([0,xMax,0 yMax])

hold off
