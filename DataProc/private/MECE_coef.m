function [ C ] = MECE_coef( shotnum )


%MECE system��choose proper coefficients to input into 'Teprofile.m'
if  shotnum<9923&&shotnum>=9651       %��Χ����
    C=0.75/0.35*[1.0000 0.9821 9.5127 0.4560 0.2151 0.4238 0.1715 0.1779 0.2465 0 0 0.1196 0.1076 4.3707 0 0.0703];
    %C9699_r=ECE_calib_update(9699,760,350,[ 10  11 15]);
    %to be abosolutely calibrated
elseif shotnum<=11000&&shotnum>=9923
    C= 1.4375*[1.0000    0.9642    9.4594    0.3890    0.5669    0.4004    0.3838    0.4188    0.4309    0    0    0.1637    0.1412 1.5731 0.5163   0.0137];
% ECE_calib_2shots
%shotnum=10037;shotnum2=10017;
%t1=[600 700];t2=[400 500];
%t1=[270 320];t2=[400 450];
%��10009�ھ��Ա궨
elseif shotnum<=11836&&shotnum>=11000           %10000�ķ�Χ���޴���

    %load C11785;
    C=[9.0074 1.5236 0.4359 3.0112 0.2376 0.8400 0.3504 0.3849 4.8450 1.9716 0.1958 0.2370 1.6063 4.6308  0  0]*1.5;
    %C11791_r=ECE_calib_update(11791,680,300,[9 15 16]);
    %abosolutely calibrated according to C11785;C(9) is modified from 0 to C11785(9)
    %load C11785;
    
    %C=9.0074*1.7*[1.0000    0.1848    0.0534    0.4227    0.0329    0.1299    0.0522    0.0653    0.0341    0.3111    0.0303    0.0373    0.2192    0.3692 0 0];
    C=9.0074*1.7*[1.0000    0.1805    0.0507    0.3807    0.0282    0.1058    0.0400    0.0475    0.0220    0.1972    0.0178    0.0207    0.1089    0.1621   0   0];
    %obtained by using ECE_calib_2shots with {shotnum=11523;shotnum2=11524; t1=456;t2=370;}
    C=9.0074*1.7*[1.0000    0.1848    0.0531    0.4121    0.0316    0.1230    0.0482    0.0593    0.0289    0.2809    0.0247    0.0256    0.0799   -0.0654  0   0];
    
    C=9.0074*1.7*[1.0000    0.1883    0.0551    0.4373    0.0343*1.05    0.1362    0.0545    0.0683    0.0341    0.3399    1*0.0308    1*0.0334    0.1222   abs(-0.0417) 0   0];
    %obtained by using ECE_calib_2shots with
    %shotnum=11523;shotnum2=11524;t1=[358 401]-23;t2=[324 378];

elseif shotnum<13057&&shotnum>=12970            %accurate range

    C=[1.0000  1.9793  5.9115 10.9757 3.5170 6.9255 3.3592 5.6018 4.7188 0 16.9165  4.0818  4.1576  1.8111 51.6639  44.0200/1.0584];
    %C12983_r=ECE_calib_update(12983,550,950,[10 16])
    %�����Ա궨to be abosolutely calibrated
    C=C* 1.0584  ;%according to C13215;C(16) is modified to C13215(16)
elseif shotnum<=13530&&shotnum>=13057       %accurate range
    %load C13215
    C=[1.0584 2.2217 6.1153 12.2164 3.7318 7.0575 3.1839 5.5135 4.8235 41.1065 18.8986 4.2042 3.9323 1.7597 50.4087 44.0200]
elseif shotnum<=14453&&shotnum>13562       %��Χ��14453����
%     load C13566   
    C=[1.0000   74.3055   14.4386   41.9891   9.5093   23.9137   8.5086   20.9507   5.3525   193.7153   39.2483   9.2947   29.6868   3.4224   -246.5518   110.5279];
    % t1=450, t2=990;
    C=[1.0000   76.4764  16.6848   58.0817   16.7049   48.1207   20.8911   56.9053   16.0375   673.1508   139.7688  34.4801   117.5921  14.9643  -1060.3  1063.1];
    %��εľ�ֵ��t1=300:10��450����;t2=t1+500��;
    
                            %���Ա궨��shot13724����Ҫ�¶�Խ��Խ��
elseif shotnum<=14787&&shotnum>=14733
    C=[6.0851    3.0842   48.7562  109.7346   22.3857   56.1159   22.2072   44.8355   12.9899    1.2780    0.9008   24.9655    1.9316   11.5632         0      0];
    %C=ECE_calib_update( 14787,1050,650,[15 16]);C14787=abscalib_Tetls(14787,C),(300,1100)֮����Ա궨
elseif shotnum<=14951&&shotnum>=14788
    C=[6.4981    2.0386   42.8300  107.3148   23.0967   53.3470   21.0595   47.1874   16.5161    0.9904    0.9193   30.8174    2.4607   17.7351         0         0];
    %C=ECE_calib_update( 14789,1100,780,[15 16]);C14789=abscalib_Tetls(14792,C),(300,1100)֮����Ա궨
elseif shotnum<=15164&&shotnum>=14952
    C=[ 4.7487    2.1264   35.7688  108.3792   17.9673   51.5158   18.0502   38.8569   12.9439    0.8832    0.7168   20.1985    1.6741   10.6497         0         0];
    %C=ECE_calib_update( 15106,700,400,[15 16]);C15106=abscalib_Tetls(15106,C),(300,1100)֮����Ա궨
elseif shotnum<=15343&&shotnum>=15165
    C=[66.4996   79.5898  143.8472  251.8683   20.4639  112.1768   58.2710  114.3032   71.0851  169.9785   76.1698   57.4063   56.9518   23.3788  299.7238  166.5289];
    %C=ECE_calib_update( 15182,440,860,[]);C15182=abscalib_Tetls(15189,C),(300,1100)֮����Ա궨
elseif shotnum>=15344&&shotnum<15407
    C=[4.1911/1.9531     1.3182    0.6765    3.7516    0.4763    1.2636    0.4838    0.6518    1.3255    1.5390    1.0876    0.7706    2.9234   32.3759         0         0];
    %C=ECE_calib_update( 15351,815,455,[15 16]);C15351= abscalib_Tetls( 15364,C ),(600,950)֮����Ա궨
    %��һ����������������C=ECE_calib_update( 15350,530,260,[15 16])
elseif shotnum>=15407&&shotnum<=15781
    %C=[66.4996   79.5898  143.8472  251.8683   20.4639  112.1768   58.2710  114.3032   71.0851  169.9785   76.1698   57.4063   56.9518   23.3788  299.7238  166.5289];
%     C=[25.4847   34.4137   73.3687  124.8211   10.6927   63.8083   31.7957   73.5372   48.9432  107.5586   45.9703   32.6097   37.5285   13.8152  192.6249  106.3536];
%     C=[21.7517   29.3728   62.6216  106.5373    9.1264   54.4616   27.1383   62.7655   41.7740   91.8034   39.2366   27.8330   32.0313   11.7915  164.4092   90.7749];
    C=[ 26.3620   35.5984   75.8944  129.1181   11.0608   66.0049   32.8903   76.0687   50.6281  111.2613   47.5528   33.7323   38.8204   14.2908  199.2560  110.0148];
    %C1=ECE_calib_update(15609,400,850,[]);C2=ECE_calib_update(15638,400,850)%���ڱ궨�����Դ���
    %C=[C1(1:7)/C1(7) C2(8:end)/C2(7)];Cabs = abscalib_Tetls( 15614,CC,[450 750] )
    C=[26.3620   35.1430   73.8369  124.0612   10.4731   61.7706   30.5568   69.6315   46.8313   99.9354   45.8063   32.3171   35.2975   12.6538  178.7851  100.7589];
elseif shotnum>=15881&&shotnum<=15923
    C=[    0.4689    0.3284    0.2258    0.9618    3.7427    0.3324    3.7304    0.2010    0.3421    0.2392    0.2390    0.1897    0.3898    2.0831         0         0];
    %C=ECE_calib_update(15902,1405,1150,[8 15 16]);֮�����adjustC.m//[trunc profile rr]=Teprofile(15903,[900:100:1300],[],C,[]);%%%%%%%%%%%%%%����2
    %Cabs=abscalib_Tetls(15911,C,[600 1500])
elseif shotnum>=15982&&shotnum<16020
    C=[1.0165*1.5   19.1026    6.7880    1.2291    4.9182   13.8751    5.1011    9.9899   11.1123    0.6917    0.9512    9.3183    0.7413    5.2465         0         0];
    %C=ECE_calib_update(15991,885,635,[15 16]);%C(1)ƫС������MECE1�����ƺ����ԣ����Ի��ɾ���λ�ÿ����������ֻ��ߵ�������ˮƽλ��Dh�����ľ���ֵ����
    %Cabs=abscalib_Tetls(15996,C,[200 800]);
elseif shotnum>=16020&&shotnum<16416
    C=[   1.2*1.3636   24.9446    7.3155    1.1565    4.9456   12.9891    4.4357    6.8621    6.5118    0.6047    0.6274    4.8255    0.3625    2.2905         0         0];
    C=1.8*0.65*0.62/0.6998*[ 0.89*0.8984*1.2*1.3636   24.9446    7.3155    1.1565    4.9456   12.9891    4.4357    6.8621    6.5118    1.0*0.6047    0.6274    4.8255    0.3625    2.2905         0         0]; 
    %C=ECE_calib_update(16179,1135,910,[ 15 16]);
    %[ Cabs ] = abscalib_Tetls( 16183,C,[700 1300] )
elseif shotnum>=16416&&shotnum<16616
    C=[1.1446   22.6206    5.3439    0.9341    4.8181   13.0070    4.3196    6.2034    5.7958    0.6375    0.6688    6.1906    0.5351    2.5030         0         0];
    %C=ECE_calib_update(16592,600,300,[15 16]);C(2)=1.05*C(2);%΢����C(2)
    %Cabs=C*(Cprevious/C) ����һ��궨ȷ��
elseif shotnum>16800&&shotnum<=17077
    C=[1.0000  0.4582 0.3100 0.3217 0.2536 0.1694 1.2791 0.3469 1.0911  0.0764 0.4938  0.2792  0.9074 0.0818  0.1208   0.0324];
    C=10*[1.0000  0.4582 0.3100 0.3217 0.2536 0.1694 1.2791 0.3469 1.0911  0.0764 0.4938  0.2792  0.9074 0.0818  0.1208   0.0324];
    C=10*[1.2000  0.4571 0.3146 0.3265 0.2666 0.1850 1.4524 0.4092 1.3880  0.1061 0.7352  0.4536  1.5827 0.1537  0.2451  0.0758];
    %δ���о��Ա궨
elseif shotnum>17077&&shotnum<17300
    
    C=[1.0000  0.85*0.4582 0.85*0.3100 0.95*0.3217 0.2536 0.1694 1.2791 0.3469 1.0911  0.0764 0.4938  0.2792  0.0245*0.9074 0.0818  0.1208  0.6*0.0324];
    C=10*[1.0000  0.4582 0.3100 0.3217 0.2536 0.1694 1.2791 0.3469 1.0911  0.0764 0.4938  0.2792  0.0245*0.9074 0.0818  0.1208   0.6*0.0324];
   % C=[1.0000 0.4162  0.2644  0.2665 0.2062  0.1361  1.0177  0.2741 0.8575  0.0593  0.3799  0.2130  0.6910  0.0615  0.0909  0.0239];�����ӵ�
    %ECE_calib_2shots
%shotnum=16897;shotnum2=16899;
%t1=[1300 1400];t2=[800 900];
%δ���о��Ա궨
    %C=10*[ 1.0000    0.5593*.7    0.3900*.75    0.3463*1.1    0.3058    0.2124    1.6425    0.4590    1.6557    0.1227    0.8527    0.5622    0.0569    0.2143    0.3658    0.0762];
elseif shotnum>=17300&&shotnum<=17871
    C=10*[  1.0000    0.5387    0.3590    0.4059   0.3184    0.2108    1.8172    0.4030    3.7351    0.0960    0.7441    0.6393    0.0522    0.1870    2.9122/10    0.0464];
    C=10*[  1.0000    0.5387    0.3590    0.4059   0.3184    0.2108    1.8172    0.4030    3.7351    0.0960    0.7441    0.6393    0.0522    0.1870    2.9122    0.0464];
    C=[1.0000    0.5493    0.3670    0.4162    0.3280    0.2116    1.8262    0.3978    3.6002    0.0899    0.6748    0.5601    0.0438    0.1618    1.7512    0.0233];
    C=[1.0000    0.5498    0.3690    0.4219    0.3297    0.2116    1.8361    0.3912    3.4838    0.0862    0.6010    0.4317    0.0213    0.0301    0.9733   -0.0296];
 
    %C=ECE_calib_update(17445,1000,700)֮���Ϊ�������궨ϵ�� C=ECE_calib_update(17445,1000��840)֮���Ϊ���ĸ��궨ϵ�� C=ECE_calib_update(17445,840,1000)
    %δ���о��Ա궨
    %SDD����У׼
elseif shotnum>=17872&&shotnum<18285
%     C=[1.1000    0.6527    0.4497    0.4767    0.3875    0.2610    2.2747    0.5033    4.9682    0.1227    0.9030    0.8662    0.0831    0.2995    6.2918    0.2527];
%     
%     %C=[1.4807    0.6527    0.4418    0.4965    0.3875    0.2541    2.2747    0.5033    5.5829    0.1227    0.9030    0.8662    0.0831    0.2995    6.2918    0.2527];
%     % 1     3     4     6     9�޸ĺ��ϵ��
% C=ECE_calib_update(17875,920,720);
    C=10*[ 1.1279^2    0.5485    0.3532    0.3442    0.2630    0.1737    1.4396    0.3120    3.4361    0.0801    0.5635    0.5237    0.0441    0.1285    2.1505    0.0653];
%SDD����У׼

elseif shotnum>=18285&&shotnum<18790
    C=12.8*[ 1.0000    0.4107    0.2505    0.2435    0.2082    0.1392    1.0490    0.2357    2.2635    0.0585    0.3797    0.3191    0.0259    0.0621    1.1399    0.0279];
%     shotnum=18300;shotnum2=18299;
%     t1=[700 800];t2=t1;
%SDD����У׼
elseif shotnum>=18790&&shotnum<19135
    C=12*[1.0000    0.5337    0.3525    0.3134    0.2481    0.1718    1.4603    0.3403    3.5700    0.0929    0.6414    0.6207    0.0582    0.1531    2.1857    0.1036];
    % C=ECE_calib_update(18984,840,650);    
    
     %C=12*[1.0000    0.4344    0.2965    0.3412    0.2481    0.1718    1.4603    0.3403    3.5700    0.0929    0.6414    0.6207    0.0582    0.1531    2.1857    0.1036];
   
     %C=12*[1.0000    0.4344    0.2965    0.3412    0.2483    0.1709    1.4784    0.3169    3.6268    0.0837    0.6874    0.6941    0.0618    0.1999    3.1350*1    0.0544*1];
%     shotnum=18998;shotnum2=18997;
%     t1=[320 360];t2=t1;
%���Ը���ǰ���ϵ�����Ա궨
elseif shotnum>=19135&&shotnum<19829
C=12* 1.173*[  1.0000    0.4471    0.3026    0.3074    0.2222    0.1492    1.2337    0.2371    2.1887    0.0523    0.3867    0.3054       0.0582*.46    0.0809    0.7478    0.0222];
elseif shotnum>=19829&&shotnum<19865
C=12* 1.173*[  1.0000    0.4471    0.3026    0.3074    0.2222    0.1492    1.2337    0.2371    2.1887    0.0523    0.3867    0.3054       0.0582*.46*4.2    0.0809    0.7478    0.0222];
%������ͬ������13��
elseif shotnum>=19865&&shotnum<20248
C=13.6*[ 1.0000    0.4659    0.2025    0.3664    0.2461    0.1583    1.4974    0.3445    2.9668    0.0764    0.6250    0.6063    0.2558    0.2006    4.2159    0.1638];
% shotnum=19870;shotnum2=19876;%ƫ����λ�ηŵ�
% t1=[350 480];t2=t1;

%C=13*[1.0000    0.4899    0.2231/1.04    0.3814    0.2919    0.2019    1.9197    0.4508    4.1196    0.1141    0.9538    0.9081    0.1895    0.2948    6.5066    0.2017];
%C=13*[1.0000    0.4899*.8    0.2231*.8    0.3814    0.2919    0.2019    1.9197    0.4508    4.1196    0.1141    0.9538    0.9081    0.1895    0.2948    6.5066    0.2017];
%C=ECE_calib_update(19975,1000,750);�����ŵ�
elseif shotnum>=20248&&shotnum<20937    
C=1.33*2*13.6*[1.0000    0.2546    0.2003    0.2994    0.1824    0.1081    0.9094    0.2058    1.3194    0.0342    0.2941    0.2469    0.0541    0.0701    0.6950    0.1018];
%shotnum=20255;shotnum2=20256;
%t1=[345 350];t2=t1;
elseif shotnum>=20937&&shotnum<21235
    %��һ��ʵ��(2013)
 C1=[ 1.0000    0.6030    0.1719    0.5546    0.2550    0.4011    0.4202    2.1770    0.7034    0.1244    1.0259    0.5859   16.7046    0.1412    3.6053    0.1584];
 % shotnum1=20937, shotnum2=20938  t1=[1000:1200] t2=[1000:1200]
 % Bt1=1.28;Bt2=1.3066
 C2=[1.0000  0.4479  0.1303  0.6281  0.2949  0.5123  0.6029  2.9116  1.0029   0.2109  1.4204  0.6945  58.1903  0.1869  3.9603  0.1085];
  % shotnum1=20943, shotnum2=20945  t1=[1000:1200] t2=[1000:1200]
 % Bt1=1.526;Bt2=1.545
 C=(C1+C2)/2;
elseif shotnum>=21235&&shotnum<21425
 C=[1.0000    0.2626    0.0825    0.4336    0.1580    0.2781    0.3322    1.5711    0.6889    0.0912    0.7141    0.3881    5.0165    0.4955    5.5553    0.4530];
%shotnum1=21237 shotnum2=21238 t1=[1200 1400]  t2=[400 600];  Bt1=1.3365
%Bt2=1.3141
elseif shotnum>=21425 && shotnum<21839
    C=[1.0000    0.5039    0.1425    0.5157    0.2070    0.3508    0.3810    1.7412    0.7653    0.1077    1.0330    0.5541    0.0291    0.1134    2.9170    0.0894];
    %shotnum1=21684 shotnum2=21685 t1=[1300 1400]  t2=[1300 1400];
    %Bt1=1.3555,Bt2=1.3814;
elseif shotnum>=21839 && shotnum< 22049
%     C1=[1.0000    0.5343    0.1502    0.5371    0.2147    0.3554    0.3443    1.7081    0.7256    0.0975    1.0489    0.6184    0.0333    0.1561    3.4596    0.2427];
   %shotnum1=21840 shotnum=21841; t1=t2=900:1000;
   C=[ 1.0000    0.5378    0.1539    0.5652    0.2366    0.4094    0.4059    2.0973    0.9258    0.1245    1.3815    0.8429    0.0465    0.2260    5.5484    0.3796];
 
   %shotnum1=21840 shotnum=21841; t1=t2=750:950;
    C=[ 1.0000    0.5467    0.1551    0.5460    0.2295    0.4041    0.4067    1.9080    0.8698    0.1187    1.1709    0.6324    0.0285    0.1185    2.1475    0.0761];
    %shotnum1=21840 shotnum=21841; t1=t2=300:400;
%    C=(C1+C2)/2;
elseif shotnum>=22049 &&shotnum<22200
    C=[ 1.0000    0.3862    0.0855    0.2861    0.0960    0.1313    0.1124    0.4822    0.1815    0.0218    0.1690    0.0680    0.0025    0.0097    0.2475    0.0066];
    %shot1=22055,shot2= 22054,t1=[290:310],t2= [690:710]
elseif shotnum>=22200  
%     &&  shotnum<22435
     C=[ 1.0000    0.5305    0.1434    0.5743    0.2476    0.4156    0.4354    2.2444    0.9467    0.1460    1.4078    0.7950    0.0447    0.2000    2.3003   0.0060];
    %shotnum1=22206;shotnum1=22207  t1=t2=300:400;
elseif shotnum>22435
    %C=[ 1.0000    0.4787    0.1104    0.2973    0.0879    0.1138    0.1129    0.3916    0.1216    0.0133    0.0927    0.0292    0.0007    0.0010    0.0067         0];
%     C=[ 1.0000    0.4230    0.1120    0.4706    0.1739    0.2700    0.2751    1.2982    0.5206    0.0726    0.6649    0.4053    0.0216    0.0908    0.4532    0.0244];
%      %shotnum1=22474,shotnum2=22475,t1=t2=420:440
    C=[ 1.0000    0.3510    0.0891    0.3807    0.1324    0.1987    0.1943    0.8774    0.3399    0.0441    0.3546    0.1935    0.0093    0.0339    0.3279    0.0163];
     %shotnum1=22474,shotnum2=22475,t1=t2=500:600
   C1=[ 1.0000    0.3756    0.0923    0.3645    0.1354    0.1991    0.2047    0.8714    0.3486    0.0501    0.4568    0.2220    0.0116    0.0402    0.9866    0.0254];
%     shotnum1=22508,shotnum2=22509 t1=t2=250:300
   C2=[  1.0000    0.3781    0.0931    0.3663    0.1340    0.1919    0.1917    0.7869    0.3045    0.0430    0.3755    0.1801    0.0097    0.0337    0.8393    0.0251];
%     shotnum1=22508,shotnum2=22509 t1=t2=300:350   
   C=(C1+C2)/2;
%   C=[ 1.0000    0.3837    0.0963    0.4086    0.1573    0.2390    0.2543    1.0958    0.4541    0.0665    0.5996    0.2954    0.0162    0.0552    1.3484    0.0358];
 %     shotnum1=22508,shotnum2=22509 t1=t2=1000:1100
else
    fprintf('No calibration coefficients!\n')
    C=[];
    return
end

end
