%% 读取文档
clc;clear
%------AF--------
%load AF folder里的所有mat文件
AF_folder=dir(fullfile( '.\ECGAF\ECGAF\*.mat'));  
% %只保留文档的名字,用于后面的数据的load
AF_folder_Name={AF_folder.name};
% AF folder里的实验者人数
AF_folder_num=size(AF_folder,1);

%创建零cell array存放所有实验者的AF数据
AF_sample=cell([AF_folder_num 1]);
AF_data=cell([AF_folder_num 1]);
%年龄 sex
AF_age=cell(AF_folder_num ,1);
AF_sex=cell(AF_folder_num ,1);

for i=1:AF_folder_num
    AF_sample{i,1}=load(AF_folder_Name{1,i});
    AF_data{i,1}=AF_sample{i,1}.ECG.data;
     AF_age{i,1}=AF_sample{i,1}.ECG.age;
     AF_sex{i,1}=AF_sample{i,1}.ECG.sex;
end

%-----------Normal-------
Normal_folder=dir(fullfile( '.\ECGnormal\ECGnormal\*.mat'));  
% %只保留文档的名字,用于后面的数据的load
Normal_folder_Name={Normal_folder.name};
% Normal folder里的实验者人数
Normal_folder_num=size(Normal_folder,1);

%创建零cell array存放所有实验者的数据
Normal_sample=cell([Normal_folder_num 1]);
Normal_data=cell([Normal_folder_num 1]);

%创建零矩阵存放实验者的年龄和性别
Normal_age=cell(Normal_folder_num ,1);
Normal_sex=cell(Normal_folder_num ,1);

for i=1:Normal_folder_num
    Normal_sample{i,1}=load(Normal_folder_Name{1,i});
    Normal_data{i,1}=Normal_sample{i,1}.ECG.data;
    Normal_age{i,1}= Normal_sample{i,1}.ECG.age;
    Normal_sex{i,1}= Normal_sample{i,1}.ECG.sex;
end



%% 滤波 计算BPM
kernel=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
leads_num=12;

AF_BPM=zeros(AF_folder_num,leads_num);
Normal_BPM=zeros(Normal_folder_num,leads_num);

AF_RR_Mean=zeros(AF_folder_num,leads_num);
Normal_RR_Mean=zeros(Normal_folder_num,leads_num);

AF_RR_Var=zeros(AF_folder_num,leads_num);
Normal_RR_Var=zeros(Normal_folder_num,leads_num);

AF_SS_Mean=zeros(AF_folder_num,leads_num);
Normal_SS_Mean=zeros(Normal_folder_num,leads_num);

AF_SS_Var=zeros(AF_folder_num,leads_num);
Normal_SS_Var=zeros(Normal_folder_num,leads_num);

%--------AF----------
for i=1:AF_folder_num    %第i个实验者
    for j=1:leads_num    %第j个lead
        %滤波 
        AF_filter_1=conv(AF_data{i,1}(j,:),kernel);
        fs=500; %采样率
        fmaxd_1=2; %截至频率
        fmaxn_1=fmaxd_1/(fs/2);
        [B,A]=butter(1,fmaxn_1,'low');
        ecg_low=filtfilt(B,A,AF_filter_1);%通过2Hz低通滤波器的信号
        AF_filter=AF_filter_1-ecg_low; %去除这一段信号，得到去基线漂移的信号`在这里插入代码片
        %Detect the R wave - signal with min amplitude of 0.1 
        %and distance between peaks >=200
        %[Y,X]
        [~,Rwave] = findpeaks(AF_filter,'MinPeakHeight',0.1,'MinPeakDistance',200);
        
        %find the averaged beat per minute
        % calculate the R-R interval time period
        RRinterval=(1:length(Rwave)-1);
        for k=1:length(Rwave)-1
         RRinterval(k)=Rwave(k+1)-Rwave(k);
        end
        Mean=mean(RRinterval);% find the average R-R interval
        AF_RR_Mean(i,j)=Mean;
        secperbeat=Mean*1/(500);%second per beat (sampling at 500Hz)
        bpm=60/secperbeat;%beat per minute
        AF_BPM(i,j)=bpm;  
        % R wave 方差
        Var=var(RRinterval);
        AF_RR_Var(i,j)=Var;
        
        
        %---------------S-S--------------
        %Detect the S wave - the inverted ECG signal with min amplitude of 0.1 
        %and distance between peaks >=200
        AF_filter_inverted = -AF_filter;
        [~,Swave] = findpeaks(AF_filter_inverted,'MinPeakHeight',0.1,'MinPeakDistance',200);
        
        SSinterval=(1:length(Swave)-1);
        for k=1:length(Swave)-1
         SSinterval(k)=Swave(k+1)-Swave(k);
        end
        MeanS=mean(SSinterval);
        AF_SS_Mean(i,j)= Mean;  
        % R wave 方差
        Var=var(SSinterval);
        AF_SS_Var(i,j)=Var;
        
    end
end

%% 
%-------------Normal-----------

for i=1:Normal_folder_num    %第i个实验者
    for j=1:leads_num    %第j个lead
        %滤波 
        Normal_filter_1=conv(Normal_data{i,1}(j,:),kernel);
        fs=500; %采样率
        fmaxd_1=2; %截至频率
        fmaxn_1=fmaxd_1/(fs/2);
        [B,A]=butter(1,fmaxn_1,'low');
        ecg_low=filtfilt(B,A,Normal_filter_1);%通过5Hz低通滤波器的信号
        Normal_filter=Normal_filter_1-ecg_low; %去除这一段信号，得到去基线漂移的信号`在这里插入代码片
        %Detect the R wave - signal with min amplitude of 0.1 
        %and distance between peaks >=200
        %[Y,X]
        [~,Rwave] = findpeaks(Normal_filter,'MinPeakHeight',0.1,'MinPeakDistance',200);
       
        %find the averaged beat per minute
        % calculate the R-R interval time period
        RRinterval=(1:length(Rwave)-1);
        for k=1:length(Rwave)-1
         RRinterval(k)=Rwave(k+1)-Rwave(k);
        end
        Mean=mean(RRinterval);% find the average R-R interval
        Normal_RR_Mean(i,j)=Mean;
        secperbeat=Mean*1/(500);%second per beat (sampling at 500Hz)
        bpm=60/secperbeat;%beat per minute
        Normal_BPM(i,j)=bpm;   
        
        %R wave 方差
        Var=var(RRinterval);
        Normal_RR_Var(i,j)=Var;
        
        %-------S-S---------
         %Detect the S wave - the inverted ECG signal with min amplitude of 0.1 
        %and distance between peaks >=200
        Normal_filter_inverted = -Normal_filter;
        [~,Swave] = findpeaks(Normal_filter_inverted,'MinPeakHeight',0.1,'MinPeakDistance',200);
        
        SSinterval=(1:length(Swave)-1);
        for k=1:length(Swave)-1
        SSinterval(k)=Swave(k+1)-Swave(k);
        end
        MeanS=mean(SSinterval);% find the average SS interval
        Normal_SS_Mean(i,j)=MeanS;
        VarS=var(SSinterval);
        Normal_SS_Var(i,j)=VarS;
        
    end
end





%% Lable 特征表格

%---AF Lable------
AF_label = strings([size(AF_BPM,1) 1]);
AF_label(:,1) = 'AF';

%---------Normal Lable----------
Normal_label = strings([size(Normal_BPM,1) 1]);
Normal_label(:,1) = 'Normal';

%---------Lable to Table
Lable_Array=[AF_label;Normal_label];
Table_Lable=array2table(Lable_Array,'VariableNames',{'Label'});


%% Age  Sex Table
Age_Array=cell2mat([AF_age;Normal_age]);
Age_Table=array2table(Age_Array,'VariableNames',{'Age'});


Sex_Array=[AF_sex;Normal_sex];
Sex_Truth=zeros(length(Lable_Array),1);
for i=1:2016
    if Sex_Array(i,1)== "Female"
        Sex_Truth(i,1)=0;
    else
        Sex_Truth(i,1)=1;     
    end
end
Sex_Table=array2table(Sex_Truth,'VariableNames',{'Sex'});

%% 1个患者对应1个lable，12个leads。1个leads有5个特征，12*5=60个特征
%[BPM_lead_1  RR_Mean_lead_1  RR_Var_lead_1  SS_Mean_lead_1  SS_Var_lead_1] 
feature_num=5;
feature_total=leads_num*feature_num;
Total_Lable_num=length(Lable_Array);
Feature_Array=zeros(Total_Lable_num,feature_total);

%% 

%-----------combine AF Normal----------
BPM=[AF_BPM;Normal_BPM];
RR_Mean=[AF_RR_Mean;Normal_RR_Mean];
RR_Var=[AF_RR_Var;Normal_RR_Var];
SS_Mean=[AF_SS_Mean;Normal_SS_Mean];
SS_Var=[AF_SS_Var;Normal_SS_Var];


for i=1:Total_Lable_num
    for j=1:5:feature_total
        Feature_Array(i,j)=BPM(i,floor(j/5)+1);
        Feature_Array(i,j+1)=RR_Mean(i,floor(j/5)+1);
        Feature_Array(i,j+2)=RR_Var(i,floor(j/5)+1);
        Feature_Array(i,j+3)=SS_Mean(i,floor(j/5)+1); 
        Feature_Array(i,j+4)=SS_Var(i,floor(j/5)+1); 
    end
end


Feature_Table=array2table(Feature_Array);
Feature_Table_total=[Feature_Table,Age_Table,Sex_Table,Table_Lable];



%% 朴素贝叶斯
Data=Feature_Table_total(:,1:62);
Lable=Feature_Table_total(:,63);
Mdl=fitcnb(Data,Lable);

% K fold
CVMdl = crossval(Mdl);
kfoldLoss(CVMdl)
