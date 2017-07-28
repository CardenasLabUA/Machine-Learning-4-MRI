function [Img, Bruker_Info] = Bruker_reader32(experiment_dir)
%Reads the 'acqp' and 'reco' and 'visu_pars' and 'method' files written out by Bruker Biospin system and
%returns scan parameters. This function also loads in the image file.

%Inputs
%experiment_dir: Directory containing the 'method', 'acqp', and 'reco'
%                files as well as the 2dseq images

%Outputs
%Img:            Image from the 2dseq file
%Bruker_Info:    A structure containing the imaging parameters from the scan

%Written by Adam Bernstein and Julio Cardenas 7/29/2014
%University of Arizona
%Email: asb2@email.arizona.edu

%Set Info Field "MagneticField"
Bruker_Info.MagneticFieldStrength = 7;
mycurrentdirectory=pwd;
cd(experiment_dir);

% Start with ACQP file
fid = fopen('acqp');
while ~feof(fid)
    tline = fgetl(fid);
    if strncmp(tline, '##$NSLICES=', 11)
        Bruker_Info.num_slices = strread(tline, '##$NSLICES=%f');
    end
    
    if strncmp(tline, '##$NR=', 50)
        Bruker_Info.NumberRepititions = strread(tline, '##$NR=%f');
    end
end
fclose(fid);

%Read METHOD file
index = 1;
fid = fopen('method');

while ~feof(fid)
    tline = fgetl(fid);
    if strncmp(tline, '##$NumCestExps=', 14)
        Bruker_Info.NumberCESTExperiments = strread(tline, '##$NumCestExps=%f');
    end
    if strncmp(tline, '##$CestArray=', 13)
        tline = fgetl(fid);
        while ~strncmp(tline, '$$', 2)
            values = cell2mat(textscan(tline, '%f'))';
            str = sprintf('cest_array%s = values;', num2str(index));
            eval(str)
            str = sprintf('array_size = size(cest_array%s, 2);', num2str(index));
            eval(str)
            tline = fgetl(fid);
            index = index + 1;
        end
        Bruker_Info.cest_array = cest_array1;
        for i = 2:(index-1)
            str = sprintf('Bruker_Info.cest_array = horzcat(Bruker_Info.cest_array, cest_array%s);', num2str(i));
            eval(str)
        end
    end
end
fclose(fid);

cd([experiment_dir,'/pdata/1'])

%Read VISU_PARS file
fid = fopen('visu_pars');

while ~feof(fid)
    tline = fgetl(fid);

    if strncmp(tline, '##OWNER=', 8); % time
        tline = fgetl(fid);
        time_str = tline(15:22);
        Bruker_Info.hour_str = time_str(1:2);
        Bruker_Info.minute_str = time_str(4:5);
        Bruker_Info.second_str = time_str(7:8);
        Bruker_Info.scan_time = str2double(Bruker_Info.hour_str)*3600+str2double(Bruker_Info.minute_str)*60+str2double(Bruker_Info.second_str);
    end
    
    if strncmp(tline, '##$VisuCoreSize',15); % img_size
        tline=fgetl(fid);
        Bruker_Info.ImageDimensions = strread(tline);
        if size(Bruker_Info.ImageDimensions, 2) == 2
            Bruker_Info.ImageDimensions(3) = 1;
        end
    end
    
    if strncmp(tline, '##$VisuAcqEchoTime',16); % TE
        tline = fgetl(fid);
        TE=strread(tline);
        Bruker_Info.TE = TE/1000;
    end
    
    if strncmp(tline, '##$VisuAcqRepetitionTime',20); % TR
        tline = fgetl(fid);
        TR = strread(tline);
        Bruker_Info.TR = TR/1000;
    end
    
    if strncmp(tline, '##$VisuCoreDataMax=',18); % Number of imgs
        Bruker_Info.Num_imgs = strread(tline,'%*s %f');
    end
    
    if strncmp(tline, '##$VisuAcqNumberOfAverages=',25); % NA
        Bruker_Info.NumberAverages = strread(tline,'##$VisuAcqNumberOfAverages=%f');
    end
end
fclose(fid);

%Read RECO file
fid = fopen( 'reco');
while ~feof(fid)
    tline = fgetl(fid); 
    if strncmp(tline, '##$RECO_map_slope', 15); %slope
        tline = fgetl(fid);
        Bruker_Info.slope = strread(tline,'%f',1);    
    end
end
fclose(fid);    

fid = fopen( '2dseq');
Img = fread(fid, 'int32');
dim4 = length(Img)/(Bruker_Info.ImageDimensions(1)*Bruker_Info.ImageDimensions(2)*Bruker_Info.ImageDimensions(3));
fclose(fid);

Img = reshape(Img, Bruker_Info.ImageDimensions(1), Bruker_Info.ImageDimensions(2), Bruker_Info.ImageDimensions(3), dim4);
Img = Img/Bruker_Info.NumberAverages;
Img = Img/Bruker_Info.slope;

for n = 1:size(Img, 4)
    for m=1:size(Img, 3)
        Img(:,:,m,n)=flipud(Img(:,:,m,n));
        Img(:,:,m,n)=imrotate(Img(:,:,m,n),-90);
    end

    
cd(mycurrentdirectory);
    
end

