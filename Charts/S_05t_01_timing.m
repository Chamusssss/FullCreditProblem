hold on

NSim = 10;
xMax = 3;
yMax = 0.017;
trueAns = 0.0089;
X = 1:xMax;

TwoLvlMCX = [ones(NSim,1)];
TwoLvlMCY = [0.10923076923076923076923076923077];

TwoLvlISX = [2*ones(NSim,1)];
TwoLvlISY = [0.1121395568097153516706754317056147];

scatter(TwoLvlMCX,TwoLvlMCY,'g')
scatter(TwoLvlISX,TwoLvlISY,'b')

legend('True Ans')
%title('S:5, l:0.10, sample num:10600')
xlabel('seconds')
ylabel('P(L > l)')
axis([0,xMax,0 yMax])

hold off
