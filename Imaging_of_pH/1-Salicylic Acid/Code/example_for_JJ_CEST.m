clear all; close all
cd ../data/pH6.77/CEST/

filenames=dir('*mM.txt');

number_of_files=length(filenames);
ppm=importdata('ppm.txt');
ppm=ppm(5:end);

for q=1:number_of_files
        S=importdata(filenames(q).name);
            Z=( S(5:end)./S(4) );
                [~,i]=max(1-Z);
                    plot(ppm-ppm(i),Z); hold all;
end

