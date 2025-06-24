% Arm Force Field Method
% copright Nicholas La Delfa & Jim Potvin (2025) CC-BY-4.0
% last updated 2025-06-24

% If you use this code or method, please cite the original publication:
% Nicholas J. La Delfa, Jim R. Potvin, 
% The ‘Arm Force Field’ method to predict manual arm strength based on only hand location and force direction, 
% Applied Ergonomics, Volume 59, Part A, 2017, Pages 410-421
% ISSN 0003-6870, https://doi.org/10.1016/j.apergo.2016.09.012. 
% https://www.sciencedirect.com/science/article/pii/S0003687016302095

clear
clc
%clf('reset')                   %clears all graphics
format short

% Coordinate system: +x axis is Lateral Left, +y is Superior, +z is Anterior

% Constants for 50th Female _______________________________________

        CV = 0.277;             % coefficient of variation for strength
        bm = 73.5;              % body mass (kg)
        ht = 1.619;              % stature (m)
        pc = 75;                % percent capable
        GravityG = [0, -1, 0];  % gravity in Global Axis System

    % CofG distance ratios 
        UAcogR = 0.436;         % Upper Arm CofG % of length, from shoudler   
        FAcogR = 0.430;         % Forearm CofG % of length, from elbow
        HcogR = 0.740;          % Hand CofG % of length, from wrist

    % Body mass ratios 
        UAmassR = 0.0280;       % Upper Arm mass as % of total body mass
        FAmassR = 0.0170;       % Forearm
        HmassR = 0.0060;        % Hand

% Variable Inputs _______________________________________
    % coordinates for each joint and the direction of the force vectors
        %   1 = Left Hand/Knuckle     
        %   2 = Left Wrist     
        %   3 = Left Elbow        
        %   4 = Left Shoulder  
        
        %   5 = Right Hand/Knuckle
        %   6 = Right Wrist 
        %   7 = Right Elbow 
        %   8 = Right Shoulder 

        %   9 = C7T1 
        %   10 = L5S1 

        %   11 = Left Hand Force Vector (in Global)
        %   12 = Right Hand Force Vector (in Global)

        % Example coordinates (m) and hand loads (N)
            xyz(1,:) = [0.29034,0.91934, 0.29134];
            xyz(2,:) = [0.26170, 0.96619, 0.23538];
            xyz(3,:) = [0.19100, 1.11068, 0.09171];
            xyz(4,:) = [0.09191, 1.30461, -0.06136];
            xyz(5,:) = [-0.37368, 1.20097, 0.360862];
            xyz(6,:) = [-0.36905, 1.20173, 0.28260];
            xyz(7,:) = [-0.42676, 1.28636, 0.09280];
            xyz(8,:) = [-0.21538, 1.32291, -0.06481];
            xyz(9,:) = [-0.06385, 1.35457, -0.06182];
            xyz(10,:) = [-0.06300, 0.95253, -0.07984];
            
            xyz(11,:) = [0.000, 1.000, 0.000];   % unit vector for force applied at the left knuckle (not the reaction force)
            xyz(12,:) = [0.000, 1.000, 0.000];   %    at the right knuckle
            LActualLoad = 50                     % actual force at the left knuckle (used to determine % capable)
            RActualLoad = 50                     %    at the right knuckle


% Calculations ____________________________________________________________________________________________    

    % hand forces in Global Axis System (GAS) (x: +Lleft, y: +Up, z: +Anterior)
        LFG = xyz(11,:);
        RFG = xyz(12,:);

    % establish Shoulder Axis System (SAS)
        SAS(3,:) = xyz(8,:)-xyz(4,:);
            temp = (SAS(3,1)^2 + SAS(3,2)^2 + SAS(3,3)^2)^0.5;
        SAS(3,:) = SAS(3,:)/temp;               % Lateral unit vector in Shoulder Axis System (z)

        Trunk = xyz(9,:)-xyz(10,:);
            temp = (Trunk(1)^2 + Trunk(2)^2 + Trunk(3)^2)^0.5;
        Trunk = Trunk/temp; 

        SASt(1,:) = cross(Trunk,SAS(3,:));       % Anterior unit vector in Shoulder Axis System (x)
            temp = (SASt(1,1)^2 + SASt(1,2)^2 + SASt(1,3)^2)^0.5;
        SAS(1,:) = SASt(1,:)/temp;

        SAS(2,:) = cross(SAS(3,:), SAS(1,:)) ;   % Superior unit vector in Shoulder Axis System (x)

    % joint locations in GAS
        LHG = xyz(1,:)-xyz(4,:);                                % Left Hand wrt gravity
        LWG = xyz(2,:)-xyz(4,:);                                % Left Wrist wrt gravity
        LEG = xyz(3,:)-xyz(4,:);                                % Left Elbow wrt gravity

        RHG = xyz(5,:)-xyz(8,:);                                % Right Hand wrt gravity
        RWG = xyz(6,:)-xyz(8,:);                                % Right Wrist wrt gravity
        REG = xyz(7,:)-xyz(8,:);                                % Right Elbow wrt gravity

    % joint locations in SAS
        LHSAS = mtimes(LHG,SAS.');                               % Left Hand wrt SAS
            LHSAS(3) = LHSAS(3)*-1;    % note: wwitches lateral polarity so left is lateral for the left arm
        LWSAS = mtimes(LWG,SAS.');                               % Left Wrist wrt SAS
            LWSAS(3) = LWSAS(3)*-1;
        LESAS = mtimes(LEG,SAS.');                               % Left Elbow wrt SAS
            LESAS(3) = LESAS(3)*-1;
            LSAS = [LHSAS; LWSAS; LESAS]; % created for display

        RHSAS = mtimes(RHG,SAS.');                               % Right remains lateral for the right arm
        RWSAS = mtimes(RWG,SAS.');                               % Right Wrist wrt SAS
        RESAS = mtimes(REG,SAS.');                               % Right Elbow wrt SAS
            RSAS = [RHSAS; RWSAS; RESAS]; % created for display - supressed for now

    % gravity in SAS
        GravitySAS = mtimes(GravityG, SAS.');     % Gravity vector wrt SAS

    % hand force unit vectors in Global
        LFG = LFG / ((LFG(1)^2 + LFG(2)^2 + LFG(3)^2)^0.5);     % Left hand force unit vector wrt Gravity
        RFG = RFG / ((RFG(1)^2 + RFG(2)^2 + RFG(3)^2)^0.5);     % Right hand force unit vector wrt Gravity
            FG = [LFG; RFG]; % created for display - supressed for now

    % hand force unit vectors in SAS
        LFSAS = mtimes(LFG, SAS.');                              % Left hand force wrt SAS
        RFSAS = mtimes(RFG, SAS.');                              % Right hand force wrt SAS
            FSAS = [LFSAS; RFSAS]; % created for display

    % segment CofG coordinates in SAS
        LUAcog = LESAS * UAcogR;                                % Left Upper Arm CofG
        LFAcog = LESAS + (LWSAS - LESAS) * FAcogR;              % Left Forearm CofG
        LHcog = LWSAS + (LHSAS - LWSAS) * HcogR;                % Left Hand CofG
            Lcog = [LUAcog; LFAcog; LHcog]; % created for display

        RUAcog = RESAS * UAcogR;                                % Right Upper Arm CofG
        RFAcog = RESAS + (RWSAS - RESAS) * FAcogR;              % Right Forearm CofG
        RHcog = RWSAS + (RHSAS - RWSAS) * HcogR;                % Right Hand CofG
            Rcog = [RUAcog; RFAcog; RHcog]; % created for display

    % Semgment weights (N)
        UAwt = bm * UAmassR * 9.81;                             % Upper arm weight (right or left) 
        FAwt = bm * FAmassR * 9.81;                             % Forearm weight
        Hwt = bm * HmassR * 9.81;                               % Hand weigth

    % Segment moments caused at shoulder
        % Left Arm
        LUAmom = cross(LUAcog,GravitySAS) * UAwt;                           % shoulder moment caused by Left Upper Arm
        LFAmom = cross(LFAcog,GravitySAS) * FAwt;                           % shoulder moment caused by Left Forearm
        LHmom = cross(LHcog,GravitySAS) * Hwt;                              % shouler moment caused by Left Hand
        LTotmom = LUAmom + LFAmom + LHmom;                                  % Left: total moment casued by gravity
        LTotmomRes = (LTotmom(1)^2 + LTotmom(2)^2 + LTotmom(3)^2)^0.5;      % Total shoulder moment caused by Left segments
            Lmom = [LUAmom, 0; LFAmom, 0; LHmom, 0; LTotmom, LTotmomRes];    % created for display

        LTotmomUV = LTotmom / LTotmomRes;
        LReach = (LHSAS(1)^2 + LHSAS(2)^2 + LHSAS(3)^2)^0.5;                 % Left reach distance
        LReachUV = LHSAS / LReach;                                          % Left reach unit vector
        Lassist = cross(LTotmomUV, LReachUV);                               % Left direction of gravity contribution to MAS
        LGFEres = LTotmomRes / LReach;                                      % Left Gravity Force Effect resultant
        LGFE = Lassist * LGFEres / norm(Lassist);                           % Left Gravity Force Effect vector
            Lvectors = [LTotmomUV, 0; LReachUV, 0; Lassist, 0; LGFE, LGFEres]; % created for display

        % Right Arm
        RUAmom = cross(RUAcog,GravitySAS) * UAwt;                           % shoulder moment caused by Right Upper Arm
        RFAmom = cross(RFAcog,GravitySAS) * FAwt;                           % shoulder moment caused by Right Forearm
        RHmom = cross(RHcog,GravitySAS) * Hwt;                              % shouler moment caused by Right Hand
        RTotmom = RUAmom + RFAmom + RHmom;                                  % Total shoulder moment caused by Right segments
        RTotmomRes = (RTotmom(1)^2 + RTotmom(2)^2 + RTotmom(3)^2)^0.5;
            Rmom = [RUAmom, 0; RFAmom, 0; RHmom, 0; RTotmom, RTotmomRes];    % created for display

        RTotmomUV = RTotmom / RTotmomRes;
        RReach = (RHSAS(1)^2 + RHSAS(2)^2 + RHSAS(3)^2)^0.5;                 % Right reach distance
        RReachUV = RHSAS / RReach;                                          % Right reach unit vector
        Rassist = cross(RTotmomUV, RReachUV);                               % Right direction of gravity contributon to MAS
        RGFEres = RTotmomRes / RReach;                                      % Right Gravity Force Effect resultant        
        RGFE = Rassist * RGFEres / norm(Rassist);                           % Right Gravity Force Effect vector
            Rvectors = [RTotmomUV, 0; RReachUV, 0; Rassist, 0; RGFE, RGFEres]; % created for display

% ___________________________________________________________________________

% Artificial Neural Network

    % ANN coefficients (18 inputs, 13 nodes)

        in1offset = [-0.4145454545454540,0.0000000000000000,-0.2000000000000000,0.0211450020362057,0.0000000000000000,0.0000000000000000,-1.0000000000000000,-1.0000000000000000,-1.0000000000000000, ...
            -0.4720000000000000,-0.5072230027811940,-0.5021474062240250,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.4720000000000000,-0.5135528785197460,-0.5326927207332780];

        in1gain = [2.2437989556135800,3.8331616396348900,2.8838355012166200,3.9808719139870100,3.7210886709687600,3.5496613795916400,1.0000000000000000,1.0000000000000000,1.0000000000000000, ...
            2.1041188125756200,2.0185852492594400,1.9796783114071000,3.7360574357315700,6.9790625815425600,13.0370893261040000,2.1186440677966100,1.9379162887927500,1.8726304261237400];

        in1min = -1;

     % Layer 1  
        Layer1b = [-0.249390763050715,-0.100250151663260,0.095639766968571,0.101680098004717,0.245807318675936,1.196620804832640,-0.640172565850211,0.504404303005015,-0.934511666215270,-0.088019829568919, ...
            -0.159442240916475,0.250629528273471,0.365963496484053];

        Layer1c = [0.0396468197036280,0.5135914913409210,-0.1615158131460990,0.5335782521154120,0.3926468149028430,0.7784365394317120,-0.0902124024580994,-0.0258725059064285,-0.3186976135026210,0.1408568009578810, ...
            -0.1015030485416650,0.2236865262109910,0.0184515483855661,0.2848686155088860,0.0966284464548426,0.2746562193587130,0.0259812747144168,0.0879899598845416
        -0.5746118508217600,-0.5027904960725710,-0.0346106776749378,0.3400709517981290,0.2055955545287030,0.3016142159084500,0.3007020715811980,-0.0492985615840271,-0.0303543877659774,0.8771017058082360, ...
        -0.7025433938175290,-0.7009810281763620,-0.4055899728486460,0.3227272624391930,0.4020874234982020,0.3610266765209100,-0.1992278230825480,-0.2430789659264560
        -0.0557183046775558,0.0490560616244855,0.0732281213380645,-0.3808356272955540,0.2394342636225010,-0.4721690370427150,0.4100268921804240,0.2704500232626520,-0.2550484876609980,0.1449649272503230, ...
        -0.4809743876368540,-0.2546406454871900,0.3465404330579830,0.0708754262465684,0.4788433598444420,0.2215891388376550,-0.4773030365151250,-0.2050673971388950
        0.3610611665751980,0.2245754371636250,0.2554306518605670,-0.4629606992678790,-0.3229704656758240,-0.5781072481730050,0.2132757087698990,0.1465432624230660,0.0376969558157140,0.4972601070092010, ...
        -0.2098604304412510,0.6795406140198610,0.1627211073321460,0.0219398458443392,-0.0658912388596506,0.4337293482561090,0.1039589318217870,0.7510560601092920
        0.2705797866315510,-0.6402561185960470,0.0848149984877266,0.3125525727147880,-0.4783292551509150,0.4418097040331320,-0.0565511811718371,-0.5084718686313590,0.1351866169817260,0.5187860156051750, ...
        0.0646553732286881,-0.2940020003199560,-0.6337273006332680,-0.2026063794675500,0.5239938084219410,-0.1340013610736720,-0.4554200176011580,0.4528727227080690
        -0.3349960804891440,-0.4357978550461050,0.4465358482261200,-0.0375157672679093,-0.7067229365556190,0.2756001172307310,0.2248372901459630,0.3299628845887790,0.1822251533236540,0.3533438912919030, ...
        0.5186303608981890,-0.9299265033661340,0.3321282113954300,0.4971312002270840,-0.2924830392757020,-0.1625864627104710,0.2983953130234120,-0.0127645885137380
        -0.2164840626998560,0.4989879882631440,0.0217402997176494,0.1149004951516220,0.5062698137551010,0.2317911665459110,0.6100899160989360,0.2567140735923170,0.0144884567379115,0.5710515402598770, ...
        -0.2946443936700260,0.0265906078074100,0.5048313837962140,-0.3792877826710210,-0.0454203902615683,0.1801269116220240,0.1372308233697510,0.0979411270189525
        0.3913441374502880,-0.0073456205366371,0.1755164415474560,0.0523500610091085,0.1243779626016780,-0.3666722086134470,0.9308840320885010,0.6956057633373600,-0.1477020919957110,-0.1508098949617220, ...
        0.2167013349329770,0.1728450799063930,-0.8902712375254020,-0.0679969677511183,0.4471591656305260,-0.0896015899235026,0.3044721158016960,0.4884305865702950
        -0.4946760291059190,0.1939477889212450,-0.1364066060231560,0.9642073332147800,-0.3823007930166790,-0.8581877087935020,0.5678240113666970,0.0183348911745601,0.0734474674712959,-0.3038824379926810,...
        -0.2105586924908500,-0.6436023360344460,-1.7207455397837500,-0.1053675074478800,0.4871803180932090,-0.4073649576265590,0.0790913947610594,0.0685909517456038
        -0.3072776347547340,0.1253616866532380,0.3129490393046960,-0.6842875264248760,0.7812594038954920,0.3974967146848980,0.1525550347153040,0.1361716306263450,0.0096179114168465,0.0172555667366513, ...
        -0.4714076889981510,0.2259879519205360,-0.2578286175424030,0.3522799993269420,-0.1745347203626280,-0.0405847095874450,-0.5426321171031580,-0.1430194565848930
        -0.0791367309675499,-0.4438012585473360,-0.4787994868168900,0.5757443432838650,0.0303014125868337,-0.6913000111042270,0.2852278064428020,-0.0843136504085405,-0.2629878149433620,-0.3190129683251240, ...
        0.0875081476021028,-0.5423870522968970,0.2141963456722130,0.4311146825841830,0.0024140834105431,0.0722952121344132,-0.3238967625044890,-0.4829474364743040
        0.1418259625408540,0.0760784150445575,-0.5283565483560720,-0.0011468411281056,-0.2474033651048630,-0.1732445962983420,-0.3189666135670600,-0.4959499926905650,0.1190199212636580,0.5403429274270760, ...
        0.0638667163600597,-0.2305280571308560,-0.0812576791473448,0.2206287611054260,-0.0940649002307790,-0.3564957901562640,-0.6196230305451320,0.1155897334030620
        0.0132468525879893,-0.8765998536722000,-0.0504484159893567,-0.0875696222338695,0.1449172648350600,0.8244006592929960,-0.0540479238588384,0.1144408906733490,-0.0399914698933143,0.4346560301984630, ...
        -0.0813723658678813,0.5068478617188480,1.0401589725829900,-0.4421187830849380,-0.2862908924579780,0.2308997829770930,-0.0823541762175492,-0.5335282381151980]; 

     % Layer 2  
        Layer2b = 0.319619509557245;

        Layer2c = [-0.726008728459171,-0.896119899552364,1.303898287490790,-0.963599520641712,0.892304171343111,-1.109496084791920,1.320456762815650,-0.881510511623027,1.197329700532090,-1.197641852318800, ...
            -0.836429746366561,-0.853341101811818,1.697291156353840];

    % Output Layer  
        OUTmin = -1;

        OUTgain = 0.0111758422535829;

        OUToffset = 44.2451801385399;

    %inputs to the ANN    (SI (y) , AP (z), ML (x))
        % Left Hand Inputs
        LANNin(1) = LHSAS(2);                                               % Hand Location wrt Shoudler (r)
        LANNin(2) = LHSAS(1);
        LANNin(3) = LHSAS(3);
        LANNin(7) = LFSAS(2);                                               % direction cosine (DC) of Force unit vector (F)
        LANNin(8) = LFSAS(1);
        LANNin(9) = -LFSAS(3);       % reverse for Left arm
        LANNin(4) = (LANNin(2)^2 + LANNin(3)^2)^0.5;                        % 2D Projection or r on plane
        LANNin(5) = (LANNin(1)^2 + LANNin(3)^2)^0.5;
        LANNin(6) = (LANNin(1)^2 + LANNin(2)^2)^0.5;
        LANNin(10) = (LANNin(2) * LANNin(9)) - (LANNin(3) * LANNin(8));     % DC of 3D moment arm (DC or r x F)
        LANNin(11) = (LANNin(3) * LANNin(7)) - (LANNin(1) * LANNin(9));
        LANNin(12) = (LANNin(1) * LANNin(8)) - (LANNin(2) * LANNin(7));
        LANNin(13) = (LANNin(10)^2 + LANNin(11)^2 + LANNin(12)^2)^0.5;      % resultant of 3D moment arm (3DMA)
        LANNin(14) = LANNin(13)^2;                                          % 3DMA^2
        LANNin(15) = LANNin(13)^3;                                          % 3DMA^3
        LANNin(16) = LANNin(4) * LANNin(7);                                 % DC of F x 2D projection of r
        LANNin(17) = LANNin(5) * LANNin(8);
        LANNin(18) = LANNin(6) * LANNin(9);

        % Right Hand Inputs (same equations as for Left above)
        RANNin(1) = RHSAS(2);      
        RANNin(2) = RHSAS(1);
        RANNin(3) = RHSAS(3);
        RANNin(7) = RFSAS(2);       
        RANNin(8) = RFSAS(1);
        RANNin(9) = RFSAS(3);      
        RANNin(4) = (RANNin(2)^2 + RANNin(3)^2)^0.5;
        RANNin(5) = (RANNin(1)^2 + RANNin(3)^2)^0.5;
        RANNin(6) = (RANNin(1)^2 + RANNin(2)^2)^0.5;
        RANNin(10) = (RANNin(2) * RANNin(9)) - (RANNin(3) * RANNin(8));
        RANNin(11) = (RANNin(3) * RANNin(7)) - (RANNin(1) * RANNin(9));
        RANNin(12) = (RANNin(1) * RANNin(8)) - (RANNin(2) * RANNin(7));
        RANNin(13) = (RANNin(10)^2 + RANNin(11)^2 + RANNin(12)^2)^0.5;
        RANNin(14) = RANNin(13)^2;
        RANNin(15) = RANNin(13)^3;
        RANNin(16) = RANNin(4) * RANNin(7);
        RANNin(17) = RANNin(5) * RANNin(8);
        RANNin(18) = RANNin(6) * RANNin(9);

        ANNin = [LANNin.', RANNin.'];        % for display

    % MaxMin Function to Modify the Original Input
        for i = 1:18        % 18 Inputs
            Lp(i) = in1gain(i) * (LANNin(i) - in1offset(i)) - 1;    % Left p values
            Rp(i) = in1gain(i) * (RANNin(i) - in1offset(i)) - 1;    % Right
        end
        p = [Lp.',Rp.'];     % for display
    
    % Layer 1 operations
        for n = 1:13        % 13 nodes
            Lsum(n) = Layer1b(n);                       % Layer 1 bias values 
            Rsum(n) = Layer1b(n);
            for i = 1:18    % 18 inputs
                Latemp(n,i) = Lp(i) * Layer1c(n,i);     % summing (p * Layer 1 coeficients)
                Lsum(n) = Lsum(n) + Latemp(n,i);        % for each Node

                Ratemp(n,i) = Rp(i) * Layer1c(n,i);     % right side
                Rsum(n) = Rsum(n) + Ratemp(n,i);
            end
            La(n) = tanh(Lsum(n));                      % TanH of sum from Layer 1
            Ra(n) = tanh(Rsum(n));    
        end
        a = [La.', Ra.'];    % for display

    % Later 2 operations
        Lsum2 = Layer2b;                                % Layer 2 bias values
        Rsum2 = Layer2b;
        for n = 1:13
            Lsum2 = Lsum2 + (La(n) * Layer2c(n));       % summing values for Layer 2
            Rsum2 = Rsum2 + (Ra(n) * Layer2c(n));
        end
       Lsum2;       % for display                        
       Rsum2;

    % ANN Output
        Lmas = (Lsum2 + 1)/OUTgain + OUToffset;          % ANN raw maximum arm strength estimates
        Rmas = (Rsum2 + 1)/OUTgain + OUToffset;
    
% __________________________________________________________________________________________________
% Bounding Strength to Min and Max values observed in our studies
    % values are bounded based on the height of the hand wrt shoulder, and the direction of the force vector
    % (eg. if the force was anterior (+1) and inferior (-1) with no
    % medial/lateral component, the codes would be 1, -1, 0)

    % arrays of minimums and maximums for each code
    MinMAS = [51.3,49.9,49.9,51.3,51.3,49.9,51.3,49.9,49.9,51.3,64.3,49.9,72.3,9999.0,49.9,51.3,72.3,49.9,52.9,49.9,49.9,52.9,52.9,49.9,52.9,49.9,49.9,47.8,47.1,47.1,47.8,47.8,47.1,47.8,47.1, ...
        47.1,47.8,64.3,47.1,68.9,9999.0,53.1,47.8,72.3,47.1,52.9,51.4,51.4,52.9,52.9,51.4,52.9,51.4,51.4,44.2,44.2,44.2,44.2,44.2,44.2,44.2,44.2,44.2,44.2,64.3,44.2,65.5,9999.0,56.3,44.2,72.3, ...
        44.2,52.9,52.9,52.9,52.9,52.9,52.9,52.9,52.9,52.9];

    MaxMAS = [223.2,223.2,223.2,223.2,223.2,223.2,223.2,223.2,223.2,223.2,120.0,223.2,133.4,-9999.0,99.3,223.2,125.7,223.2,184.1,184.1,184.1,184.1,184.1,184.1,184.1,184.1,184.1,220.5,220.5,220.5,220.5, ...
        199.6,220.5,210.4,210.4,199.6,220.5,149.8,220.5,165.5,-9999.0,134.7,210.4,128.0,199.6,190.9,190.9,181.9,190.9,175.5,181.9,190.9,190.9,177.1,217.7,217.7,217.7,217.7,176.0,217.7,197.7,197.7, ...
        176.0,217.7,179.7,217.7,197.7,-9999.0,170.2,197.7,130.2,176.0,197.7,197.7,179.7,197.7,166.9,179.7,197.7,197.7,170.2];

    % Height Code   (-1 for < - 0.01, 0 for between -0.01 & 0.01, 1 for > 0.01)
        Htband = 0.01;
        Lht = 0;
        if LHSAS(2) < -Htband
            Lht = -1;
        elseif LHSAS(2) > Htband
            Lht = 1;
        end
        Rht = 0;
        if RHSAS(2) < -Htband
            Rht = -1;
        elseif RHSAS(2) > Htband
            Rht = 1;
        end

    % Ant/Post Force Code (-1 for negative, 0 for 0, +1 for positive)
        if LFSAS(1) == 0
            Lap = 0;
        else
            Lap = LFSAS(1)/abs(LFSAS(1));
        end
        if RFSAS(1) == 0
            Rap = 0;
        else
            Rap = RFSAS(1)/abs(RFSAS(1));
        end

    % Sup/Inferior Force Code  (-1 for negative, 0 for 0, +1 for positive)
        if LFSAS(2) == 0
            Lsi = 0;
        else
            Lsi = LFSAS(2)/abs(LFSAS(2));
        end
        if RFSAS(2) == 0
            Rsi = 0;
        else
            Rsi = RFSAS(2)/abs(RFSAS(2));
        end

    % Med/Lateral Force Code  (-1 for negative, 0 for 0, +1 for positive)
        if -LFSAS(3) == 0           % must switch polarity for left side
            Lml = 0;
        else
            Lml = -LFSAS(3)/abs(LFSAS(3));  %do these need to be switched?
        end
        if RFSAS(3) == 0
            Rml = 0;
        else
            Rml = RFSAS(3)/abs(RFSAS(3));  %do these need to be switched?
        end
        
        Lcode = [Lht, Lap, Lsi, Lml];
        Rcode = [Rht, Rap, Rsi, Rml];

    % creating 4D arrays of Minimums and Maximums from MinMax array above
        c = 0;
        for ht = 1:3
            for ap = 1:3
                for si = 1:3
                    for ml = 1:3
                        c = c + 1;
                        Min(ht, ap, si, ml) = MinMAS(c);
                        Max(ht, ap, si, ml) = MaxMAS(c);
                        temp = [c, ht-2, ap-2, si-2, ml-2, Min(ht, ap, si, ml), Max(ht, ap, si, ml)];
                    end
                end
            end
        end

    % add 2 to each code for look-up in array
        Lmin = Min(Lht+2, Lap+2, Lsi+2, Lml+2);
        Lmax = Max(Lht+2, Lap+2, Lsi+2, Lml+2);
        Rmin = Min(Rht+2, Rap+2, Rsi+2, Rml+2);
        Rmax = Max(Rht+2, Rap+2, Rsi+2, Rml+2);
        
    % bounding to Mininum and Maximum MAS values
        if Lmas < Lmin
            Lmas = Lmin;
        end
        if Lmas > Lmax
            Lmas = Lmax;
        end
        if Rmas < Rmin
            Rmas = Rmin;
        end
        if Rmas > Rmax
            Rmas = Rmax;
        end
        
      %  t = 'Mean 0G MAS values';
      %  Lmas   % for display;
      %  Rmas;
% _____________________________________________________________________________________________        
% Estimate of Zero-Gravity maximum arm strength (MAS) for selected population
    Lsd = Lmas * CV;                         % estimate of standard deviation based on mean and global CV value
    Rsd = Rmas * CV;

    prob = 1 - pc/100;
   % t = 'MAS - including gravity and percentile';
    L0gMAS = norminv(prob, Lmas, Lsd);       % Zero G MAS values based on percent capable selected
    R0gMAS = norminv(prob, Rmas, Rsd);
    
% component of gravity acting along the force vector (gravity assist)
    Lga = (LGFE(1) * LFSAS(1)) + (LGFE(2) * LFSAS(2)) + (LGFE(3) * -LFSAS(3));       % must reverse lateral side for Left
    Rga = (RGFE(1) * RFSAS(1)) + (RGFE(2) * RFSAS(2)) + (RGFE(3) * RFSAS(3));
    
% Final MAS value (with Gravity)
    LMASwG = L0gMAS + Lga;
    RMASwG = R0gMAS + Rga;
    
% Percentage capable of the actual loads
    LActualLoad0G = LActualLoad - Lga;
    RActualLoad0G = RActualLoad - Rga;
    % t = 'Percent Capable - including gravity and percentile';
    LMASprob = (1 - normcdf(LActualLoad0G, Lmas, Lsd)) * 100;
    RMASprob = (1 - normcdf(RActualLoad0G, Rmas, Rsd)) * 100;
    
    showStrength = [LMASwG, RMASwG]
    showPerCap = [LMASprob, RMASprob]
               
