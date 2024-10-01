%%Wrapper for SILEX ans Cavity-Plate test case
classdef wrapperSILEX < handle
    properties
        shellExecute='Main_porous_cavity7.sh';  %shell script for executing computation
        resultFile='results/cavity7_with_porous.mat';  %mat-file provided by python
        nbProc=20; %number of processors for MPI
        %parameters
        paraValFull=[];
        paraValHist=[];
        %results
        resultStore=[];
        resultSave=[];
        varResult={};
        varResultFinal={};
        %log and save file
        FileBased='porous_cavity';
        logFileFull='';
        saveFileFull='';
    end
    methods
        %constructor
        function obj=wrapperSILEX()
            dt=datetime;
            fprintf('============================\n');
            fprintf(' >>> Initialize Wrapper <<<\n');
            fprintf(' %s\n',datestr(dt));
            %generate logfile and savefile
            obj.logFileFull=[datestr(dt,'YYYY-mm-DD_HH-MM-SS_') obj.FileBased '.log'];
            obj.saveFileFull=[datestr(dt,'YYYY-mm-DD_HH-MM-SS_') obj.FileBased '.mat'];
            %show files
            fprintf(' > Logfile: %s\n',obj.logFileFull);
            fprintf(' > Savefile: %s\n',obj.saveFileFull);
            %show information
            fprintf(' > Number of processors: %i\n',obj.nbProc);
            fprintf(' > Shell script: %s\n',obj.shellExecute);
            fprintf(' > Results from SILEX: %s\n',obj.resultFile);
            %write on logFile
            fid = fopen(obj.logFileFull,'w');fclose(fid);
            writeLog(obj,'Initialize \n','date');
        end
        %setter for paraValHist
        function set.paraValHist(obj,paraVal)
            obj.paraValHist(end+1,:)=paraVal;
        end
        %setter for paraValHist
        function set.varResult(obj,resIn)
            obj.varResult{end+1}=resIn;
        end
        %execute one run
        function execOne(obj,paraVal)
            fprintf(' > Run for parameters: ');
            fprintf('%d ',paraVal);fprintf('\n');
            obj.paraValHist=paraVal;
            %write in logfile
            writeLog(obj,['Run Parameters: ' sprintf('%g ',paraVal) '\n'],'date');
            %create command
            cmdSH=['./' obj.shellExecute  ' ' sprintf('%g ',paraVal) ' ' sprintf('%g',obj.nbProc) '>>' obj.logFileFull];
            %execute
            system(cmdSH);
            %write in logfile
            writeLog(obj,'STOP\n','date');
        end
        %save result
        function saveCurrent(obj)
            %write in MatlabFile
            S.varResult=obj.varResult;
            S.paraValFull=obj.paraValFull;
            S.paraValHist=obj.paraValHist;
            save(obj.saveFileFull,'-struct','S');
            writeLog(obj,'Save current results\n','date');
        end
        %save final result
        function saveFinal(obj)
            %write in MatlabFile
            S.varResult=obj.varResultFinal;
            S.paraValFull=obj.paraValHist;
            save(obj.saveFileFull,'-struct','S');
            writeLog(obj,'Save final results\n','date');
        end
        %save current result
        function S=readResult(obj)
            %write in logfile
            writeLog(obj,'Read results','date');
            %load matlab file
            S=load(obj.resultFile);
            obj.varResult=S;
            %delete file
            delete(obj.resultFile);
        end
        %write on logfile
        function writeLog(obj,txt,type)
            if nargin<3;type='simple';end
            fileId=fopen(obj.logFileFull,'a');
            switch type
                case 'simple'
                case 'date'
                    txt=['[' datestr(datetime,'YYYY/mm/DD HH-MM-SS') '] ' txt];
            end
            fprintf(fileId,txt);
            fclose(fileId);
        end
        %execute and read result(s) of one or many runs
        function outS=compute(obj,paraIn)
            %ordering data
            sP=[size(paraIn,1) size(paraIn,2) size(paraIn,3)];
            if size(paraIn,3)~=1
                obj.paraValFull=reshape(paraIn,[sP(1)*sP(2) sP(3)]);
            else
                obj.paraValFull=paraIn;
            end
            
            fprintf(' >> Computation for %i sets of parameters\n',size(obj.paraValFull,2))
            writeLog(obj,['Run Parameters: ' sprintf('Computation for %i sets of parameters\n',size(paraV,2)) '\n'],'date');
            
            %along the parameters
            for itV=1:size(obj.paraValFull,2)
                %execute computation
                obj.execOne(obj.paraValFull(itV,:));
                %read result
                obj.readResult;
                %save current result
                obj.saveCurrent;
            end
            
            %reordering data
            obj.varResultFinal=cell(sP(1),sP(2));
            obj.varResultFinal=obj.varResult;
            outS=obj.varResultFinal;
            obj.paraValFull=paraIn;
            %save final
            obj.saveFinal;
        end
    end
end