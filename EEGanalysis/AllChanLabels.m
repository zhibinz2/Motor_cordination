%% Left hemisphere
LeftLabels=[1 
4 
6 7 8 9 10 
17 18 19 20 
26 27 28 29 
35 36 37 38 
44 
46 47 48 49 
55 56 
60 
63 

65 
67 68 
72 73 74 
78 
80 81 82 83 
88 89 
92 93 94 95 
100 101 102 103 
108 109 110 111 
116 
118 119 
122 123 124];

'FP1',
'AF3',
'F11','F7','F5','F3','F1',
'FT11','FC5','FC3','FC1',
'T7','C5','C3','C1',
'TP7','CP5','CP3','CP1',
'M1',
'P7','P5','P3','P1',
'PO7','PO3',
'O1',
'CB1',

'AFP1',
'AF7','AF5',
'AFF5H','AFF3H','AFF1H',
'F9',
'FFT7H','FFC5H','FFC3H','FFC1H',
'FT9','FT7',
'FTT7H','FCC5H','FCC3H','FCC1H',
'TTP7H','CCP5H','CCP3H','CCP1H',
'TPP7H','CPP5H','CPP3H','CPP1H',
'P9',
'PPO3H','PPO1H',
'PO9','PO5','PO1',

%% Right hemisphere
RightLabels=[3 
5 
12 13 14 15 16 
22 23 24 25 
31 32 33 34 
40 41 42 43 
45 
51 52 53 54 
58 59 
62 
64 

66 
70 71 
75 76 77 
79 
84 85 86 87 
90 91 
96 97 98 99 
104 105 106 107 
112 113 114 115 
117 
120 121 
125 126 127];

'FP2',
'AF4',
'F2','F4','F6','F8','F12',
'FC2','FC4','FC6','FT12',
'C2','C4','C6','T8',
'CP2','CP4','CP6','TP8',
'M2',
'P2','P4','P6','P8',
'PO4','PO8',
'O2',
'CB2',

'AFP2',
'AF6','AF8',
'AFF2H','AFF4H','AFF6H',
'F10',
'FFC2H','FFC4H','FFC6H','FFT8H',
'FT8','FT10',
'FCC2H','FCC4H','FCC6H','FTT8H',
'CCP2H','CCP4H','CCP6H','TTP8H',
'CPP2H','CPP4H','CPP6H','TPP8H',
'P10',
'PPO2H','PPO4H',
'PO2','PO6','PO10',

%% Central
CentralLabels=[2 11 21 30 39 50 57 61 ...
    69 128];

'FPZ','FZ','FCZ','CZ','CPZ','PZ','POZ','OZ',
'AFZ','CBZ',

%% Additional channels
'VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'

%% All channels
'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'

'FP1-1','FPZ-2','FP2-3','AF3-4','AF4-5','F11-6','F7-7','F5-8','F3-9','F1-10','FZ-11','F2-12','F4-13','F6-14','F8-15','F12-16','FT11-17','FC5-18','FC3-19','FC1-20','FCZ-21','FC2-22','FC4-23','FC6-24','FT12-25','T7-26','C5-27','C3-28','C1-29','CZ-30','C2-31','C4-32','C6-33','T8-34','TP7-35','CP5-36','CP3-37','CP1-38','CPZ-39','CP2-40','CP4-41','CP6-42','TP8-43','M1-44','M2-45','P7-46','P5-47','P3-48','P1-49','PZ-50','P2-51','P4-52','P6-53','P8-54','PO7-55','PO3-56','POZ-57','PO4-58','PO8-59','O1-60','OZ-61','O2-62','CB1-63','CB2-64','AFP1-65','AFP2-66','AF7-67','AF5-68','AFZ-69','AF6-70','AF8-71','AFF5H-72','AFF3H-73','AFF1H-74','AFF2H-75','AFF4H-76','AFF6H-77','F9-78','F10-79','FFT7H-80','FFC5H-81','FFC3H-82','FFC1H-83','FFC2H-84','FFC4H-85','FFC6H-86','FFT8H-87','FT9-88','FT7-89','FT8-90','FT10-91','FTT7H-92','FCC5H-93','FCC3H-94','FCC1H-95','FCC2H-96','FCC4H-97','FCC6H-98','FTT8H-99','TTP7H-100','CCP5H-101','CCP3H-102','CCP1H-103','CCP2H-104','CCP4H-105','CCP6H-106','TTP8H-107','TPP7H-108','CPP5H-109','CPP3H-110','CPP1H-111','CPP2H-112','CPP4H-113','CPP6H-114','TPP8H-115','P9-116','P10-117','PPO3H-118','PPO1H-119','PPO2H-120','PPO4H-121','PO9-122','PO5-123','PO1-124','PO2-125','PO6-126','PO10-127','CBZ-128','VEOG-129','HEOG-130','EMG1-131','EMG2-132','HL 1-133','HL 2-134','EMG3-135','EMG4-136','EMG5-137','EMG6-138','TRIGGER-139'


%%
load('Neuroscan_spherical.mat')
plot3(x,y,z,'ro');
labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
text(x,y,z,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');

plot(x,y,'ro');
labels = num2cell(string([1:128]));
text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');


plot3(x,y,z,'ro');
labels = {'FP1-1','FPZ-2','FP2-3','AF3-4','AF4-5','F11-6','F7-7','F5-8','F3-9','F1-10','FZ-11','F2-12','F4-13','F6-14','F8-15','F12-16','FT11-17','FC5-18','FC3-19','FC1-20','FCZ-21','FC2-22','FC4-23','FC6-24','FT12-25','T7-26','C5-27','C3-28','C1-29','CZ-30','C2-31','C4-32','C6-33','T8-34','TP7-35','CP5-36','CP3-37','CP1-38','CPZ-39','CP2-40','CP4-41','CP6-42','TP8-43','M1-44','M2-45','P7-46','P5-47','P3-48','P1-49','PZ-50','P2-51','P4-52','P6-53','P8-54','PO7-55','PO3-56','POZ-57','PO4-58','PO8-59','O1-60','OZ-61','O2-62','CB1-63','CB2-64','AFP1-65','AFP2-66','AF7-67','AF5-68','AFZ-69','AF6-70','AF8-71','AFF5H-72','AFF3H-73','AFF1H-74','AFF2H-75','AFF4H-76','AFF6H-77','F9-78','F10-79','FFT7H-80','FFC5H-81','FFC3H-82','FFC1H-83','FFC2H-84','FFC4H-85','FFC6H-86','FFT8H-87','FT9-88','FT7-89','FT8-90','FT10-91','FTT7H-92','FCC5H-93','FCC3H-94','FCC1H-95','FCC2H-96','FCC4H-97','FCC6H-98','FTT8H-99','TTP7H-100','CCP5H-101','CCP3H-102','CCP1H-103','CCP2H-104','CCP4H-105','CCP6H-106','TTP8H-107','TPP7H-108','CPP5H-109','CPP3H-110','CPP1H-111','CPP2H-112','CPP4H-113','CPP6H-114','TPP8H-115','P9-116','P10-117','PPO3H-118','PPO1H-119','PPO2H-120','PPO4H-121','PO9-122','PO5-123','PO1-124','PO2-125','PO6-126','PO10-127','CBZ-128'};
text(x,y,z,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');


%%
% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};
% AllchanNames{1}

%% Calculate subplot coordinates for scalp map
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load('coordinates.mat');
X = [chanlocs128.X];
Y = [chanlocs128.Y];
Z = [chanlocs128.Z];
% labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};
labels = {'1-FP1', '2-FPZ', '3-FP2', '4-AF3', '5-AF4', '6-F11', '7-F7', '8-F5', '9-F3', '10-F1', '11-FZ', '12-F2', '13-F4', '14-F6', '15-F8', '16-F12', '17-FT11', '18-FC5', '19-FC3', '20-FC1', '21-FCZ', '22-FC2', '23-FC4', '24-FC6', '25-FT12', '26-T7', '27-C5', '28-C3', '29-C1', '30-CZ', '31-C2', '32-C4', '33-C6', '34-T8', '35-TP7', '36-CP5', '37-CP3', '38-CP1', '39-CPZ', '40-CP2', '41-CP4', '42-CP6', '43-TP8', '44-M1', '45-M2', '46-P7', '47-P5', '48-P3', '49-P1', '50-PZ', '51-P2', '52-P4', '53-P6', '54-P8', '55-PO7', '56-PO3', '57-POZ', '58-PO4', '59-PO8', '60-O1', '61-OZ', '62-O2', '63-CB1', '64-CB2', '65-AFP1', '66-AFP2', '67-AF7', '68-AF5', '69-AFZ', '70-AF6', '71-AF8', '72-AFF5H', '73-AFF3H', '74-AFF1H', '75-AFF2H', '76-AFF4H', '77-AFF6H', '78-F9', '79-F10', '80-FFT7H', '81-FFC5H', '82-FFC3H', '83-FFC1H', '84-FFC2H', '85-FFC4H', '86-FFC6H', '87-FFT8H', '88-FT9', '89-FT7', '90-FT8', '91-FT10', '92-FTT7H', '93-FCC5H', '94-FCC3H', '95-FCC1H', '96-FCC2H', '97-FCC4H', '98-FCC6H', '99-FTT8H', '100-TTP7H', '101-CCP5H', '102-CCP3H', '103-CCP1H', '104-CCP2H', '105-CCP4H', '106-CCP6H', '107-TTP8H', '108-TPP7H', '109-CPP5H', '110-CPP3H', '111-CPP1H', '112-CPP2H', '113-CPP4H', '114-CPP6H', '115-TPP8H', '116-P9', '117-P10', '118-PPO3H', '119-PPO1H', '120-PPO2H', '121-PPO4H', '122-PO9', '123-PO5', '124-PO1', '125-PO2', '126-PO6', '127-PO10', '128-CBZ'};

sort(Z) %just to examine
ZZ=((max(Z)-min(Z))*1.5-(Z-min(Z)))./((max(Z)-min(Z))*1.5); %Make higher Z smaller fraction
sort(ZZ) %just to examine

XXPLOT=(10*Y)*0.95;
YYPLOT=(10*X)*0.8;

XXPLOT=XXPLOT.*ZZ;
YYPLOT=YYPLOT.*ZZ;

XXPLOT=-1*XXPLOT;
YYPLOT=YYPLOT;

XXPLOT=XXPLOT*0.6;
YYPLOT=YYPLOT*0.7;

XXPLOT=XXPLOT+0.5;
YYPLOT=YYPLOT+0.5;

figure('units','normalized','outerposition',[0 0 1 1]);
plot(XXPLOT,YYPLOT,'r.');

% labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
text(XXPLOT,YYPLOT,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');