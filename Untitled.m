%load the AF data from file
clear;
clc;
AF_output = dir(fullfile('.\ECGAF\ECGAF\*.mat')); %list the files in AF_folder
AF_file_names = {AF_output.name}; % creat cell arrays for AF files names
AF_num = size(AF_output,1); %number of AF

AF_file_dir_names = strings(AF_num,1); %return 1098 by 1 array

%Create cells array to store AF data for all experimenters
AF_ECG = cell([AF_num 1]);  %create 1098*1 empty cells array
AF_data = cell([AF_num 1]);
AF_age=cell(AF_num ,1);
AF_sex=cell(AF_num ,1);

for i=1:AF_num
    AF_file_dir_names(i,1) = ['.\ECGAF\ECGAF\',AF_file_names{1,i}];
    AF_ECG{i,1} = load(AF_file_dir_names(i,1)); %assign values i row 1 coloum
    AF_data{i,1} = AF_ECG{i,1}.ECG.data;
    AF_age{i,1}= AF_ECG{i,1}.ECG.age;
    AF_sex{i,1}= AF_ECG{i,1}.ECG.sex;
end
%--------------------------------------------------------------
%load the normal data from file
Normal_output = dir(fullfile('.\ECGnormal\ECGnormal\*.mat'));
Normal_file_names = {Normal_output.name}; %normal files names
Normal_num = size(Normal_output,1); %number of normal

Normal_file_dir_names = strings(Normal_num,1); %return 918 by 1 arra

%Create cells array to store Normal data for all experimenters
Normal_ECG = cell([Normal_num 1]);
Normal_data = cell([Normal_num 1]);
Normal_age=cell(Normal_num ,1);
Normal_sex=cell(Normal_num ,1);

for i=1:Normal_num
    Normal_file_dir_names(i,1) = ['.\ECGnormal\ECGnormal\',Normal_file_names{1,i}];
    Normal_ECG{i,1} = load(Normal_file_dir_names(i,1));
    Normal_data{i,1} = Normal_ECG{i,1}.ECG.data;
    Normal_age{i,1}= Normal_ECG{i,1}.ECG.age;
    Normal_sex{i,1}= Normal_ECG{i,1}.ECG.sex;
end

%----------------------------------------------------------------------
% Filter the images and calculate variance
kernel=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %kernel for the low pass filter
leads_num = 12;
%----------------------------------------------------
% calculate AF patient variance
AF_i_lead_j_RR_Mean = zeros([AF_num 12]);
AF_i_lead_j_SS_Mean = zeros([AF_num 12]);
AF_i_lead_j_RR_variance = zeros([AF_num 12]);  %return 1098 by 12 arrays repreasents R waves 12 leads of 1 patient
AF_i_lead_j_SS_variance = zeros([AF_num 12]);  %S ways
AF_i_lead_j_BPM = zeros([AF_num 12]);% return 1098 by 12 arrays about BPM

for i=1:AF_num % the i experimenter
    for j=1:leads_num  % the j leads
        % Assign lead_j data
        AF_i_lead_j = AF_data{i,1}(j,:); %get the (j,:) in the contents of {i,1}
        AF_i_lead_j_filtered = conv(AF_i_lead_j,kernel);%apply convolution for filtering out the signals       
    %K2= medfilt2(AF_i_lead_j,[5 5]); %5 × 5 
    %K3= medfilt2(AF_i_lead_j,[7 7]); %7 × 7
    %K4= medfilt2(AF_i_lead_j,[9 9]); %9 × 9

        fs=500; 
        fmaxd_1=2;
        fmaxn_1=fmaxd_1/(fs/2);
        [B,A]=butter(1,fmaxn_1,'low');
        ecg_low=filtfilt(B,A,AF_i_lead_j_filtered);%Signal passing through a 2Hz low-pass filter
        AF_filter=AF_i_lead_j_filtered-ecg_low; %Remove this segment of signal to get a signal that removes baseline drift
        
        %--------------------------------------
        %Detect the R wave - signal with min amplitude of 0.1 and distance between
        %peaks >=200
        [~,AF_i_lead_j_Rwave] = findpeaks(AF_i_lead_j_filtered,'MinPeakHeight',0.1,'MinPeakDistance',200);
        %calculate the RR interval time period
        RRinterval = (1:length(AF_i_lead_j_Rwave)-1);
        for k=1:length(AF_i_lead_j_Rwave)-1
            RRinterval(k) = AF_i_lead_j_Rwave(k+1) - AF_i_lead_j_Rwave(k);
        end
        MeanR = mean(RRinterval);% find the average R-R interval
        AF_i_lead_j_RR_Mean(i,j)=MeanR;
        %Caculate the BPM for jth lead of the patient
        AF_i_lead_j_secperbeat =MeanR*1/(500);
        AF_i_lead_j_BPM(i,j) = 60/AF_i_lead_j_secperbeat;
        %Caculate the variance of R wave for jth lead of ith patient
        AF_i_lead_j_RR_variance(i,j) = var(RRinterval);
       
        %--------------------------------------------------------
        %Detect the S wave - the inverted ECG signal with min amplitude of 0.1 and distance between
        %peaks >=200
        AF_i_lead_j_ECG_interted = -AF_i_lead_j_filtered;
        [~,AF_i_lead_j_Swave] = findpeaks(AF_i_lead_j_ECG_interted,'MinPeakHeight',0.1,'MinPeakDistance',200);
        % Caculate the distance of peak of S wave
        SSinterval = (1:length(AF_i_lead_j_Swave)-1);
        for k=1:length(AF_i_lead_j_Swave)-1
            SSinterval(k) = AF_i_lead_j_Swave(k+1) - AF_i_lead_j_Swave(k);
        end
        MeanS = mean(SSinterval);% find the average SS interval
        AF_i_lead_j_SS_Mean(i,j)=MeanS;
        %Caculate the variance of S wave for jth lead of ith patient
        AF_i_lead_j_SS_variance(i,j) = var(SSinterval);
    end
end   
%------------------------------------------------------
% calculate Normal people variance
Normal_i_lead_j_RR_Mean = zeros([Normal_num 12]);
Normal_i_lead_j_SS_Mean = zeros([Normal_num 12]);
Normal_i_lead_j_RR_variance = zeros([Normal_num 12]);  %return 1098 by 12 arrays repreasents R waves 12 leads of 1 normal
Normal_i_lead_j_SS_variance = zeros([Normal_num 12]);  %S ways
Normal_i_lead_j_BPM = zeros([Normal_num 12]);% return 1098 by 12 arrays about BPM

for i=1:Normal_num %i person
    for j=1:leads_num %j leads
        Normal_i_lead_j = Normal_data{i,1}(j,:);
        Normal_i_lead_j_filtered = conv(Normal_i_lead_j,kernel);%apply convolution for filtering out the signals
        
        fs=500; %
        fmaxd_1=2; %
        fmaxn_1=fmaxd_1/(fs/2);
        [B,A]=butter(1,fmaxn_1,'low');
        ecg_low=filtfilt(B,A,Normal_i_lead_j_filtered);% 5Hz low pass
        Normal_filter=Normal_i_lead_j_filtered-ecg_low; 
        %--------------------------------------
        %Detect the R wave - signal with min amplitude of 0.1 and distance between
        %peaks >=200
        [~,Normal_i_lead_j_Rwave] = findpeaks(Normal_i_lead_j_filtered,'MinPeakHeight',0.1,'MinPeakDistance',200);
         %calculate the RR interval time period
        RRinterval = (1:length(Normal_i_lead_j_Rwave)-1);
        for k=1:length(Normal_i_lead_j_Rwave)-1
            RRinterval(k) = Normal_i_lead_j_Rwave(k+1) - Normal_i_lead_j_Rwave(k);
        end
        MeanR = mean(RRinterval);% find the average R-R interval
        Normal_i_lead_j_RR_Mean(i,j)=MeanR;
        %Caculate the BPM for jth lead of the patient
        Normal_i_lead_j_secperbeat =MeanR*1/(500);
        Normal_i_lead_j_BPM(i,j) = 60/Normal_i_lead_j_secperbeat;
        %Caculate the variance of R wave for jth lead of ith patient
        AF_i_lead_j_RR_variance(i,j) = var(RRinterval);
        
        %--------------------------------------------------------------
        %Detect the S wave - the inverted ECG signal with min amplitude of 0.1 and distance between
        %peaks >=200
        Normal_i_lead_j_ECG_interted = -Normal_i_lead_j_filtered;
        [~,Normal_i_lead_j_Swave] = findpeaks(Normal_i_lead_j_ECG_interted,'MinPeakHeight',0.1,'MinPeakDistance',200);
        % Caculate the distance of peak of S wave
        SSinterval = (1:length(Normal_i_lead_j_Swave)-1);
        for k=1:length(Normal_i_lead_j_Swave)-1
            SSinterval(k) = Normal_i_lead_j_Swave(k+1) - Normal_i_lead_j_Swave(k);
        end
        MeanS = mean(SSinterval);% find the average SS interval
        Normal_i_lead_j_SS_Mean(i,j)=MeanS;
        %Caculate the variance of S wave for jth lead of ith normal
        Normal_i_lead_j_SS_variance(i,j) = var(SSinterval);
    end
end
    
%------------------------------------------------
   %Create table for classifier
   AF_label = strings([AF_num 1]);
   AF_label(:,1) = 'AF';
   
   Normal_label = strings([Normal_num 1]);
   Normal_label(:,1) = 'Normal'; 
   
   %Each column of label becomes a variable in table
   label_combine = [AF_label;Normal_label];
   %label_combine = [AF_label.',Normal_label.'];
   label_table = array2table(label_combine,'VariableNames',{'Label'}); 
   %age table
   Age_label=cell2mat([AF_age;Normal_age]);
   Age_Table=array2table(Age_label,'VariableNames',{'Age'});
   
   %sex table
   Sex_label=[AF_sex;Normal_sex];
   Sex_Truth=zeros(length(label_combine),1);
for i=1:2016
    if Sex_label(i,1)== "Female"
        Sex_Truth(i,1)=0;
    else
        Sex_Truth(i,1)=1;     
    end
end
Sex_Table=array2table(Sex_Truth,'VariableNames',{'Sex'});
%[BPM_lead_1  RR_Mean_lead_1  RR_Var_lead_1  SS_Mean_lead_1  SS_Var_lead_1] 
feature_num=5;
feature_total=leads_num*feature_num;
Total_Lable_num=length(label_combine);
Feature_Array=zeros(Total_Lable_num,feature_total);

BPM=[AF_i_lead_j_BPM;Normal_i_lead_j_BPM];
RR_Mean=[AF_i_lead_j_RR_Mean;Normal_i_lead_j_RR_Mean];
RR_Var=[AF_i_lead_j_RR_variance;Normal_i_lead_j_RR_variance];
SS_Mean=[AF_i_lead_j_SS_Mean;Normal_i_lead_j_SS_Mean];
SS_Var=[AF_i_lead_j_SS_variance;Normal_i_lead_j_SS_variance];
  
for i=1:Total_Lable_num
    for j=1:5:feature_total
        Feature_Array(i,j)= BPM(i,floor(j/5)+1);
        Feature_Array(i,j+1)=RR_Mean(i,floor(j/5)+1);
        Feature_Array(i,j+2)=RR_Var(i,floor(j/5)+1);
        Feature_Array(i,j+3)=SS_Mean(i,floor(j/5)+1); 
        Feature_Array(i,j+4)=SS_Var(i,floor(j/5)+1); 
    end
end
Feature_Table=array2table(Feature_Array);
Feature_Table_total=[Feature_Table,Age_Table,Sex_Table,label_table];
