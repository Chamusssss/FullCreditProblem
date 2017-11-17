hold on

xMin = 0;
xMax = 4500;
yMin = 0;
yMax = 0.00008;
trueAns =  0.000054139945393358308464029462188075;

OneLvlISX = [600,600,600,600,600,...
    1200,1200,1200,1200,1200,...
    1800,1800,1800,1800,1800,...
    2400,2400,2400,2400,2400,...
    3000,3000,3000,3000,3000,...
    3600,3600,3600,3600,3600,...
    4200,4200,4200,4200,4200
             ];
        
OneLvlISY = [0.000051042369407049549643456964886923, 0.000067310716026318522810754918506149, 0.000048603165541661117040814210632504, 0, 0.000033236150753846131625842547085981,...
0.000047655366216505102422689582608584, 0.000056301235361382295221618909941697, 0.000051882117541775882188590657273508, 0, 0.000046726055612218626064257448460282,...
0.000061205184114872240993751595361516, 0.000042037343331109728466771918720113, 0.000040446124775129788661912066949355, 0.000053190961538799554542047282046369, 0.000040591321622507866900744083471508,...
0.000052790713817333679864357565136856, 0.000058332107526688423138061523331999, 0, 0.000049859305701841280652683585650209, 0,...
0.000065969783425809283387841741586044, 0.00003750128988565839789973144213775, 0.000060335171910682647126739769349868, 0.000052151010597798655823217200966369, 0.00004932678034863185139021965475159,...
0.000048080901046087013653990149819606, 0.000041266716981021253380246954778343, 0.000046798099636901527677897844048616, 0.000052638639862224287824417345849071, 0.000049784892262948987856456539713434,...
0.000062591861996224003946050207236595, 0.000059649834749873700020898759088794, 0, 0.000052133088324953556042754693189423, 0.000055343795033012367027419609044614
];
   
OneLvlISYMean = arrayfun(@(i) mean(OneLvlISY((5*i-4):(5*i))), 1:(length(OneLvlISY)/5));

s1 = scatter(OneLvlISX,OneLvlISY,'g', 'filled');

r1 = refline(0,trueAns);
r1.Color = 'r';
r1.LineStyle = '--';

legend('1LvlIS','True Ans');
xlabel('Number of samples used to train pi*');
ylabel('P(L > l)');
axis([xMin,xMax,yMin,yMax]);

hold off
