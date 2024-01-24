clc

cd(evalin('base','exp_path'));

Filescsv = dir(strcat('*.csv'));
Filesxlsx = dir(strcat('*.xlsx'));

Files = [Filesxlsx;Filescsv];

if length(Files) >1
disp(['Single run started, but more than one datasheet (',num2str(length(Files)),' datasheets) found. Codes shifted to BATCH run with CONSTANT inputs. Good luck!'])
else
disp('Single run started - one datasheet found. Good luck!')
end
assignin('base','Files',Files);

cd ..


for tttttttt = 1:length(Files)
    files_ = evalin('base','Files');
    disp('###################################################################################################################')
    disp(['Datasheet ',num2str(tttttttt),' of ',num2str(length(files_))])
    file_name = files_(tttttttt).name; 
    assignin('base','datasheet_name',file_name);
    run ECM_Easy.m
end

disp('###################################################################################################################')
disp('Congratulations!!! Mission accomplished !!! Please go to the folder \Results\ to check all the results.')
disp('(For debugging (or donation) purpose, please contact tao.zhu@imperial.ac.uk :D)')
