hold on

xMin = 0;
xMax = 350;
yMin = -20;
yMax = 0;
trueAns =  0.0053089481291651767420404487985914;
X = 1:xMax;

TwoLvlMCX = [38.3438,38.3438,38.3438,38.3438,38.3438,...
    85.5469,85.5469,85.5469,85.5469,85.5469,...
    222.9375,222.9375,222.9375,222.9375,222.9375,...
             ];
        
TwoLvlMCY = [ 0.015399866912103214039442278249226, 0.0060738329191282244140515089725341, 0.0055027766975200695029335484775856, 0.0079470215101107091748566801925335, 0.014764386147504356822102167257071,...
0.0084339464231081712197735100744467, 0.001285850380314723997998704874135, 0.0094036774839297062461573872838017, 0.011392245576521612196452259979651, 0.0094036774839322805757957368655298,...
 0.010182197784780425142692195095151, 0.0091742890971545365541572891743272, 0.0024300756243040462500903942100194, 0.0095667226668890471830142274711761, 0.0064025004740006601078317061137568,...
];

TwoLvlMCYMean = arrayfun(@(i) mean(TwoLvlMCY((5*i-4):(5*i))), 1:(length(TwoLvlMCY)/5));

OneLvlISX = [56.6875,56.6875,56.6875,56.6875,56.6875,...
    103.6094,103.6094,103.6094,103.6094,103.6094,...
    240.375,240.375,240.375,240.375,240.375,...
            ];
        
OneLvlISY = [0.00083185249450099820432702468764319, 0.00046589650523166231937260972628678, 0.00076208376537766195270989788568272, 0.00074658683742487828589179388316666, 0.0011194766644655738171282299830978,...
             0.00056878031689221999686600916845691, 0.00065325289585038530133093148677403, 0.0006741402966588267830355674092857, 0.00063378892559478889473512097652019, 0.00080759380842911499681635589809048,...
             0.00065714618215382571402710665253721, 0.00077099350956053166409459587171682, 0.00092370737111763115759210940325374, 0.00065132217879650990133288956940305, 0.00067734593438821457528353775501273,...
             ];

OneLvlISYMean = arrayfun(@(i) mean(OneLvlISY((5*i-4):(5*i))), 1:(length(OneLvlISY)/5));

TwoLvlISX = [87.2656,87.2656,87.2656,87.2656,87.2656,...
    126.3594,126.3594,126.3594,126.3594,126.3594,...
    248.2188,248.2188,248.2188,248.2188,248.2188,...
            ];
        
TwoLvlISY = [0.00049020493339274899460139556239824, 0.00053364781111672721872901670181477, 0.00040068716502519142281732444921261, 0.00036373501838045233826302027324573, 0.00054452560712132036177535177046138,...
             0.00043480040468063649219260247136276, 0.00039074855065322349227316389708164, 0.00050685980530886648714322673114907, 0.00037766770555116283740806903956866, 0.00046175369927030926811561961642383,...
             0.00032155730274799766265497735773238, 0.00039846594404843044007896724068019, 0.00047765766150163882280618721232202, 0.00027892365359956979714572966599917, 0.00039653781697697123740861679941361,...
         ];

TwoLvlISYMean = arrayfun(@(i) mean(TwoLvlISY((5*i-4):(5*i))), 1:(length(TwoLvlISY)/5));

OneLvlISEX = [34.9375,34.9375,34.9375,34.9375,34.9375,...
    70.25,70.25,70.25,70.25,70.25,...
    134.875,134.875,134.875,134.875,134.875,...
            ];
        
OneLvlISEY = [0.015073749645676882621958547758823, 0.003176306771052253129922204877289, 0.009395921857079858260930471658412, 0.0051293046359196011602099574133717, 0.0073510452362445530619106115466366,...
    0.011399672606093860446896393057159, 0.0075053761397202841884612745104732, 0.0042838293966263934695049009349077, 0.0083820174263332467706755224412518, 0.0063918949908048717026276541730567,...
    0.0037098761012412232764001718265945, 0.0064244110156403038194850019237947, 0.0072140612678426975848400282131934, 0.011520820751136952309057193133413, 0.011705649551243289124835555981008,...
             ];

OneLvlISEYMean = arrayfun(@(i) mean(OneLvlISEY((5*i-4):(5*i))), 1:(length(OneLvlISEY)/5));

% GISX = [54.9375,54.9375,54.9375,54.9375,54.9375,...
%     102.8125,102.8125,102.8125,102.8125,102.8125,...
%     193.2344,193.2344,193.2344,193.2344,193.2344,...
%             ];
%         
% GISY = [0.00073436109230224603117848092281861, 0.0012372192757783761693501922707128, 0.00038127225843657660056540414927895, 0.00034426423955001577868403894733262, 0.0011115512096425321855730494746695,...
%     0.0010503936016067324729639231861711, 0.000031929577307915751871719017840334, 0.00018855787983682566004949021643, 0.00068597344259690366634341085472215, 0.00018283859989662480725718107787969,...
%     0.0016237153984460057573602220770681, 0.000034404740168452812253390565855327, 0.0001323420461603407147395394805045, 0.0001887586098645013452362179107169, 0.00012036504452725141020314608697461,...
%              ];
% 
% GISYMean = arrayfun(@(i) mean(GISY((5*i-4):(5*i))), 1:(length(GISY)/5));

OneLvlISRatio = (OneLvlISYMean ./ TwoLvlMCYMean) .* 100;

TwoLvlISRatio = (TwoLvlISYMean ./ TwoLvlMCYMean) .* 100;

OneLvlISERatio = (OneLvlISEYMean ./ TwoLvlMCYMean) .* 100;

% GISRatio = (GISYMean ./ TwoLvlMCYMean) .* 100;

plot(unique(TwoLvlMCX),log10(TwoLvlMCYMean),'r','LineWidth',3)
plot(unique(OneLvlISX),log10(OneLvlISYMean),'g','LineWidth',3)
plot(unique(OneLvlISEX),log10(OneLvlISEYMean),'black','LineWidth',3)
plot(unique(TwoLvlISX),log10(TwoLvlISYMean),'b','LineWidth',3)
% plot(unique(GISX),log(GISYMean),'y','LineWidth',3)

axis([xMin,xMax,yMin,yMax]);

% r1 = refline(0,log(trueAns));
% r1.Color = 'r';
% r1.LineStyle = '--';

% legend('2LvlMC','1LvlISZ','1LvlISE','2LvlIS','GlassermanIS');
legend('2LvlMC','1LvlISZ','1LvlISE','2LvlIS');
xlabel('Runtime (Seconds)');
ylabel('log(Varience)');


hold off