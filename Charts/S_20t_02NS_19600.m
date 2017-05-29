hold on

NSim = 10;
xMax = 3;
yMax = 0.0001;
trueAns = 5.4000e-05;
X = 1:NSim;

TwoLvlMC = [ 0, 0.00027110066871498283029098138442075, 0.0001807337791433218868606542562805, 0.000090366889571660943430327128140249, 0, 0, 0, 0, 0.00027110066871498283029098138442075, 0];

TwoLvlIS = [ 0.000052982884825502731874526107791468, 0.000054218606872382901188008041959421, 0, 0.000051941165079023722918934863335139, 0.000054824878114741391042265272748324, 0.000047764114608294980343812519141267, 0.000048389623749758529501919374071761, 0.000052599925372850377320362719180125, 0.000050505792356074061982048467589479, 0.000050507977046288157080018016431566];

plot([0 X],repmat(trueAns,1,NSim + 1),'r')

scatter(ones(NSim,1),TwoLvlMC,'g')
scatter(2*ones(NSim,1),TwoLvlIS,'b')

legend('True Ans')
%title('S:20, l:0.20, sample num:10600')
xlabel('')
xticks([0 1 2])
xticklabels({'','2LvlMC(6533,13066)','2LvlIS(4600,5000,10000)'})
ylabel('P(L > l)')
axis([0,xMax,0 yMax])

hold off
