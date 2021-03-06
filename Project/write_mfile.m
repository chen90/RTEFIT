function write_mfile(mfile,path_mfile)
%% read the EFIT input data from mfile
% -inputs:
% mfile             struct, the data of the mfile
% path_mfile    string, the path of the mfile
% -outputs:

% Edited by Shuying SUN in 2019/05/27
% Modified by Bin Chen in 2019/05/31

listname_mfile=fieldnames(mfile); %fields name of the structure of mfile
fid=fopen(path_mfile,'w');
fprintf(fid,['!This file is generated by MATLAB script, ',...
    'please do not change it unless you understand all the meaning of it!\n']);
for i=1:length(listname_mfile)
    listvalue=getfield(mfile,listname_mfile{i}); %eval(['listvalue=mfile.',listname_mfile{i}])
    fprintf(fid,['&',listvalue.namelist,'\n']);
    sprintf(['&',listvalue.namelist,'\n']);
    for j=1:length(listvalue.varname)
        varname =listvalue.varname{j};
        varvalue=listvalue.varvalue{j};
        vartype =listvalue.vartype{j};
        varlength=length(varvalue);
        if(j==length(listvalue.varname)) % for the last variable
            comma='';
        else % for others
            comma=',';
        end
        if(varlength==1)
            if(strcmpi(vartype,'double'))
                valuestr=num2str(varvalue,'%12e');
            elseif(strcmpi(vartype,'integer'))
                valuestr=num2str(varvalue,'%d');
            elseif(strcmpi(vartype,'logical'))
                if(varvalue)
                    valuestr='T';
                elseif(~varvalue)
                    valuestr='F';
                end
            end
            fprintf(fid,[' ',varname,' = ',valuestr,comma,'\n']);
        else % varlength>1
            if(strcmpi(vartype,'double'))
                fprintf(fid,[' ',varname,' = ',num2str(varvalue(1),'%12e'),',\n']);
                for k=2:(varlength-1)
                    fprintf(fid,['            ',num2str(varvalue(k),'%12e'),',\n']);
                end
                fprintf(fid,['            ',num2str(varvalue(end),'%12e'),comma,'\n']);
            elseif(strcmpi(vartype,'integer'))
                fprintf(fid,[' ',varname,' = ',num2str(varvalue(1),'%d'),',\n']);
                for k=2:(varlength-1)
                    fprintf(fid,['            ',num2str(varvalue(k),'%d'),',\n']);
                end
                fprintf(fid,['            ',num2str(varvalue(end),'%d'),comma,'\n']);
            end
        end % if(varlength)
    end
    fprintf(fid,'/\n');
end % for length(listname_mfile)
fclose(fid);

end
