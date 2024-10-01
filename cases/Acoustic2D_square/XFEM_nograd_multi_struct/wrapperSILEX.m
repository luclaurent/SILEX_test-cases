%% Wrapper for SILEX and Cavity-Plate test case
%% Luc Laurent -- luc.laurent@lecnam.net


classdef wrapperSILEX < handle
    properties
        debugMod=false;     % display additional information
        %scripts and files for running
        folderWork='exportSILEX';
        shellExecute='Main.sh';  %shell script for executing computation
        pythonExecute='python3';
        mpiExecute='mpirun';
        caseDefine='';
        %files for running
        resultFile='results/xfem_3_results.mat';
        pythonMakeMesh={};
        pythonCompute={'Main_xfem_3_mod_complex.py'};%{'Main_xfem_3_mod_complex.py'};
        %MPI options
        mpiOptNbCores='-np';
        varOpenblas='export OPENBLAS_NUM_THREADS=1';    %number of thread(s) per execution (=1 when running MUMPS sequential)
        varMKL='export MKL_NUM_THREADS=1';              %avoid running MKL parallel for classical scipy functions
        multiThread=false;                              %no multithreads per processors when running using mumps sequential
        forceMPI=false;                                 %force MPI even if 1 processor is requested
        %running cases
        numCurrent=0;   %current number of computations
        numTotal=0;     %total number of computations (for running case from a last one)
        statusFlag=false; %flag for status of computation
        initFlag=false; %flag for initialization (true=done)
        daemon; %declaring a situation with a OFFLINE/ONLINE mode (depending of the size of pythonCompute)
        %parameters
        paraValFull=[];
        paraValFullOrder=[];
        paraValHist=[];
        %results
        resultStore=[];
        resultSave=[];
        varResult={};
        varResultFinal={};
        %log and save file
        FileBased='xfem_3';
        logFileFull='';
        saveFileFull='';
        %time of the execution
        elapsedInitMeshCPUTime=[];
        elapsedInitMeshTime=[];
        elapsedMeshTime=[];
        elapsedMeshCPUTime=[];
        elapsedInitTime=[];
        elapsedInitCPUTime=[];
        elapsedComputeTime=[];
        elapsedComputeCPUTime=[];
        elapsedTotalTime=[];
        elapsedTotalCPUTime=[];
        %parameters for running a case
        nbProc=1;      %number of processors for MPI
        nbSteps=NaN;    %number of steps in the frequency range
        freqMax=NaN;    %maximum frequency
        freqMin=NaN;    %minimum frequency
        nbModesFluid=NaN;   %number of modes used for fluid (Craig-Bampton)
        nbModesSolid=NaN;   %number of modes used for solid (Craig-Bampton)
    end
    methods
        %constructor
        function obj=wrapperSILEX(paraVal,nbProc)
            dt=datetime;
            fprintf('============================\n');
            fprintf(' >>> Initialize Wrapper <<<\n');
            fprintf(' %s\n',datestr(dt));
            %generate logfile and savefile
            obj.logFileFull=[obj.folderWork '/' datestr(dt,'YYYY-mm-DD_HH-MM-SS_') obj.FileBased '.log'];
            obj.saveFileFull=[obj.folderWork '/' datestr(dt,'YYYY-mm-DD_HH-MM-SS_') obj.FileBased '.mat'];
            %create folder if it does not exist
            if ~exist(obj.folderWork,'dir');mkdir(obj.folderWork);end
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
            %create last log file
            obj.lastLog;
            %
            if nargin>1
                obj.nbProc=nbProc;
            end
            if nargin>0
                obj.compute(paraVal);
            end
        end
        %setter for paraValHist
        function set.paraValHist(obj,paraVal)
            obj.paraValHist(end+1,:)=paraVal;
        end
        %setter for varResult
        function set.varResult(obj,resIn)
            obj.varResult{end+1}=resIn;
        end
        %setter for elapsedInitMeshCPUTime
        function set.elapsedInitMeshCPUTime(obj,resIn)
            obj.elapsedInitMeshCPUTime(end+1)=resIn;
        end
        %setter for elapsedMeshCPUTime
        function set.elapsedMeshCPUTime(obj,resIn)
            obj.elapsedMeshCPUTime(end+1)=resIn;
        end
        %setter for elapsedComputeCPUTime
        function set.elapsedComputeCPUTime(obj,resIn)
            obj.elapsedComputeCPUTime(end+1)=resIn;
        end
        %setter for elapsedInitCPUTime
        function set.elapsedInitCPUTime(obj,resIn)
            obj.elapsedInitCPUTime(end+1)=resIn;
        end
        %setter for elapsedInitMeshTime
        function set.elapsedInitMeshTime(obj,resIn)
            obj.elapsedInitMeshTime(end+1)=resIn;
        end
        %setter for elapsedMeshTime
        function set.elapsedMeshTime(obj,resIn)
            obj.elapsedMeshTime(end+1)=resIn;
        end
        %setter for elapsedComputeTime
        function set.elapsedComputeTime(obj,resIn)
            obj.elapsedComputeTime(end+1)=resIn;
        end
        %setter for elapsedInitTime
        function set.elapsedInitTime(obj,paraVal)
            obj.elapsedInitTime(end+1,:)=paraVal;
        end
        %getter for elapsedTotalTime
        function TT=get.elapsedTotalTime(obj)
            TT=sum(obj.elapsedComputeTime)+sum(obj.elapsedMeshTime)+sum(obj.elapsedInitTime);
        end
        %getter for elapsedTotalCPUTime
        function TT=get.elapsedTotalCPUTime(obj)
            TT=sum(obj.elapsedComputeCPUTime)+sum(obj.elapsedMeshCPUTime)+sum(obj.elapsedInitCPUTime);
        end
        %getter for daemon or not
        function TT=get.daemon(obj)
            TT=(iscell(obj.pythonCompute)&&numel(obj.pythonCompute)>1);
        end
        %setter for freqMin
        function obj=set.freqMin(obj,value)
            if value>0
                obj.changeValue('minimal frequency',obj.freqMin,value,true);
                obj.freqMin=value;
            else
                fprintf('Minimal frequency must be positive (value: %d)\n',value);
            end
        end
        %setter for freqMaxtxt
        function obj=set.freqMax(obj,value)
            if value>0
                obj.changeValue('maximum frequency',obj.freqMax,value,true);
                obj.freqMax=value;
            else
                fprintf('Maximum frequency must be positive (value: %d)\n',value);
            end
        end
        %setter for nbProc
        function obj=set.nbProc(obj,value)
            %obtain number of processors
            [~,s]=system('lscpu|grep "^CPU(s):"|tr -dc ''0-9''');
            nbProcAvail=str2double(s);
            if ~(value<=nbProcAvail);value=nbProcAvail;end
            if value>0&&floor(value)==value
                obj.changeValue('number of processors',obj.nbProc,value,true);
                obj.nbProc=value;
            else
                fprintf('Number of processors must be a positiv integer value (value: %d)\n',value);
            end
        end
        %setter for nbSteps
        function obj=set.nbSteps(obj,value)
            if value>1&&floor(value)==value
                obj.changeValue('number of frequency steps',obj.nbSteps,value,true);
                obj.nbSteps=value;
            else
                fprintf('Number of processors must be a positive integer value larger than 1 (value: %d)\n',value);
            end
        end
        %create last.log file (linked to the last log)
        function lastLog(obj)
            lastLogFile='last.log';
            if exist(lastLogFile,'file')==2
                delete(lastLogFile);
            end
            %create link
            system(['ln -s ' obj.logFileFull ' ' lastLogFile]);
        end
        %customize system command (for dealing with errors)
        function systemRUN(obj,cmd)
            obj.statusFlag=false;
            %run command
            [statFlag,r]=system(cmd);
            if statFlag~=0
                obj.printWindows('################\n################\nExecution error\n################\n################\n');
                obj.printWindows(['Executed command: ' cmd '\n']);
                obj.printWindows(['Log:\n################\n' r '\n\n################\n\n\n']);
            else
                obj.statusFlag=true;
            end
        end
        %execute meshing
        function meshExecute(obj,paraVal)
            if numel(obj.pythonMakeMesh)~=0
                %command
                printWindows(obj,['MESHING parameters: ' sprintf('%g ',paraVal) '\n']);
                writeLog(obj,['Start MESHING (parameters: ' sprintf('%g ',paraVal) ')\n'],'date');
                %build final command
                if obj.daemon
                    if obj.initFlag
                        pythonScript= obj.pythonMakeMesh{2};
                    else
                        pythonScript= obj.pythonMakeMesh{1};
                    end
                else
                    pythonScript= obj.pythonCompute;
                end
                cmdMesh=[ obj.pythonExecute  ' ' pythonScript ' ' sprintf('%g ',paraVal) '>> ' obj.logFileFull ' 2>&1'];
                obj.systemRUN(cmdMesh);
                writeLog(obj,'End MESHING\n','date');
            else
                printWindows(obj,'no MESH\n');
                writeLog(obj,'no MESH\n','date');
            end
        end

        %build command for running computation w/- or w/o mpi
        function cmd=buildCmdRun(obj,paraVal)
            if nargin==1;paraVal=NaN;end
            cmd=[];
            cmdOptPython=[];
            %specific options for ONLINE/OFFLINE mode
            modROM = obj.daemon&&obj.initFlag;
            %mpi part
            if ~obj.multiThread
                cmd=[cmd obj.varMKL ';'];
                if obj.numCurrent==1
                    obj.printWindows(' >> No MultiThreading on MKL\n')
                end
            end
            if (obj.nbProc>1||obj.forceMPI)&&~isnan(obj.nbProc)   
                cmd=[cmd obj.varOpenblas ';' obj.mpiExecute ' ' obj.mpiOptNbCores ' ' num2str(obj.nbProc)];
                if obj.numCurrent==1
                    obj.printWindows([' >> Number of processors :' num2str(obj.nbProc) '\n'])
                end
            end
            %opt for mechanical definition
            if ~isnan(obj.freqMax)&&(modROM||~obj.daemon)
                cmdOptPython=[cmdOptPython ' -F ' num2str(obj.freqMax)];
                if obj.numCurrent==1
                    obj.printWindows([' >> Maximum frequency :' num2str(obj.freqMax) '\n'])
                end
            end
            if ~isnan(obj.freqMin)&&(modROM||~obj.daemon)
                cmdOptPython=[cmdOptPython ' -f ' num2str(obj.freqMin)];
                if obj.numCurrent==1
                    obj.printWindows([' >> Minimum frequency :' num2str(obj.freqMin) '\n'])
                end
            end
            if ~isnan(obj.nbSteps)&&(modROM||~obj.daemon)
                cmdOptPython=[cmdOptPython ' -s ' num2str(obj.nbSteps)];
                if obj.numCurrent==1
                    obj.printWindows([' >> Number of frequency steps :' num2str(obj.nbSteps) '\n'])
                end
            end
            if ~any(isnan(paraVal))&&(modROM||~obj.daemon)
                if numel(paraVal)>1
                    txtTmp=sprintf('%g,',paraVal(1:end-1));
                else
                    txtTmp='';
                end
                txtTmp=[txtTmp sprintf('%g',paraVal(end))];
                cmdOptPython=[cmdOptPython ' -p ' txtTmp];                
                if obj.numCurrent==1
                    obj.printWindows([' >> Number of frequency steps :' sprintf('%g',paraVal) '\n'])
                end
            end
            
            if ~isempty(obj.caseDefine)&&(modROM||~obj.daemon)
                cmdOptPython=[cmdOptPython ' -c ' obj.caseDefine];
                if obj.numCurrent==1
                    obj.printWindows([' >> Case defined :' num2str(obj.caseDefine) '\n'])
                end
            end
            %build final command
            if obj.daemon
                if obj.initFlag
                    pythonScript= obj.pythonCompute{2};
                else
                    pythonScript= obj.pythonCompute{1};
                end
            else
                pythonScript= obj.pythonCompute{1};
            end
            cmd=[cmd ' ' obj.pythonExecute ' ' pythonScript ' ' cmdOptPython ' >> ' obj.logFileFull ' 2>&1'];
        end
        %initialization (ONLINE/OFFLINE mode)
        function init(obj,paraVal)
            if obj.daemon&&~obj.initFlag
                %build mesh
                obj.printWindows('Start INITIALIZATION MESH\n');
                writeLog(obj,'Start INITIALIZATION MESH\n','date');
                countTime=mesuTime;
                obj.meshExecute(paraVal);
                countTime.stop();
                writeLog(obj,'End INITIALIZATION MESH\n','date');
                obj.printWindows('End INITIALIZATION MESH\n');
                %store the times for meshing
                obj.elapsedInitMeshTime=countTime.tElapsed;
                obj.elapsedInitMeshCPUTime=countTime.tCPUElapsed;
                %create command
                cmdRun=obj.buildCmdRun;
                %execute
                obj.printWindows('Start INITIALIZATION COMPUTE\n');
                writeLog(obj,'Start INITIALIZATION COMPUTE\n','date');
                countTime=mesuTime;
                obj.systemRUN(cmdRun);
                countTime.stop();
                writeLog(obj,'End INITIALIZATION COMPUTE\n','date');
                obj.printWindows('End INITIALIZATION COMPUTE\n');
                %store the times
                obj.elapsedInitTime=countTime.tElapsed;
                obj.elapsedInitCPUTime=countTime.tCPUElapsed;
                %
                obj.initFlag=true;
            end
        end
        %execute one run
        function res=execOne(obj,paraVal)
            % increment of the coutner of computations
            obj.numCurrent=obj.numCurrent+1;
            %
            obj.printWindows(['>> Run for parameters: ' sprintf('%g ',paraVal) '\n']);
            obj.paraValHist=paraVal;
            %write in logfile
            writeLog(obj,['Run Parameters: ' sprintf('%g ',paraVal) '\n'],'date');
            %build mesh
            countTime=mesuTime;
            obj.meshExecute(paraVal);
            countTime.stop();
            %store the times for meshing
            obj.elapsedMeshTime=countTime.tElapsed;
            obj.elapsedMeshCPUTime=countTime.tCPUElapsed;
            %initialization (OFFLINE/ONLINE mode)
            obj.init(paraVal);
            %build mesh
            obj.printWindows('Start MESH\n');
            writeLog(obj,'Start MESH\n','date');
            countTime=mesuTime;
            obj.meshExecute(paraVal);
            countTime.stop();
            writeLog(obj,'End MESH\n','date');
            obj.printWindows('End MESH\n');
            %store the times for meshing
            obj.elapsedMeshTime=countTime.tElapsed;
            obj.elapsedMeshCPUTime=countTime.tCPUElapsed;
            %create command
            cmdRun=obj.buildCmdRun(paraVal);
            %execute
            obj.printWindows('Start COMPUTE\n');
            writeLog(obj,'Start COMPUTE\n','date');
            countTime=mesuTime;
            obj.systemRUN(cmdRun);
            countTime.stop();
            writeLog(obj,'End COMPUTE\n','date');
            obj.printWindows('End COMPUTE\n');
            %store the times
            obj.elapsedComputeTime=countTime.tElapsed;
            obj.elapsedComputeCPUTime=countTime.tCPUElapsed;
            %write in logfile
            %writeLog(obj,['STOP ' sprintf('%g',countTime.tElapsed) ' ' sprintf('%g',countTime.tCPUElapsed) ' \n'],'date');
            writeLog(obj,'STOP\n');
            %read result
            res=obj.readResult;
            %save current result
            obj.saveCurrent;
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
            S.paraValFull=obj.paraValFull;
            save(obj.saveFileFull,'-struct','S');
            writeLog(obj,'Save final results\n','date');
        end
        %save current result
        function S=readResult(obj)
            %write in logfile
            writeLog(obj,'Read results\n','date');
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
            txtOut='wrapperSILEX';
            switch type
                case 'simple'
                case 'date'
                    txtOut=['[' txtOut ' - ' datestr(datetime,'YYYY/mm/DD HH:MM:SS.FFF') '] '];
            end
            txtOut=[txtOut '(' sprintf('%i',obj.numCurrent) '/' sprintf('%i',obj.numTotal) ') ' txt];
            fprintf(fileId,txtOut);
            fclose(fileId);
        end
        %write on log in the command window
        function printWindows(obj,txt)
            txtPre='wrapperSILEX||';
            txtPre=[txtPre '(' sprintf('%i',obj.numCurrent) '/' sprintf('%i',obj.numTotal) ') '];
            %add prefix in every new lines except for if the string finish
            %with newline
            if numel(txt)>1
                if strcmp(txt(end-1:end),'\n')
                    txtOut=strrep(txt(1:end-2),'\n',['\n' txtPre]);
                    txtOut=[txtOut '\n'];
                else
                    txtOut=strrep(txt,'\n',['\n' txtPre]);
                end
            else
                txtOut=txt;
            end
            fprintf([txtPre txtOut]);
        end
        %execute and read result(s) of one or many runs
        function outS=compute(obj,paraIn)
            %ordering data
            sP=[size(paraIn,1) size(paraIn,2) size(paraIn,3)];
            if sP(3)~=1
                obj.paraValFull=reshape(paraIn,[sP(1)*sP(2) sP(3)]);
            else
                obj.paraValFull=paraIn;
            end
            obj.paraValFullOrder=obj.paraValFull;

            fprintf(' >> Computation for %i sets of parameters\n',size(obj.paraValFull,1))
            writeLog(obj,['Run Parameters: ' sprintf('Computation for %i sets of parameters\n',size(obj.paraValFull,1)) '\n'],'date');

            %along the parameters
            nbRun=size(obj.paraValFull,1);
            obj.numTotal=nbRun;
            for itV=1:nbRun
                %execute computation
                obj.execOne(obj.paraValFull(itV,:));
            end
            %reordering data
            if sP(3)~=1
                obj.varResultFinal=reshape(obj.varResult,[sP(1) sP(2)]);
            else
                obj.varResultFinal=reshape(obj.varResult,[sP(1) 1]);
            end
            outS=obj.varResultFinal;
            obj.paraValFull=paraIn;
            %save final
            obj.saveFinal;
        end
        %function for checking changes on value
        function flag=changeValue(obj,txt,oldValue,newValue,flagLog)
            flag=false;
            if nargin<5;flagLog=false;end
            if newValue~=oldValue
                flag=true;
                txtOk=['Change of ' txt ' to ' sprintf('%g ',newValue) '(old value:' sprintf('%g ',oldValue) ')\n'];
                obj.printWindows(txtOk);
                if flagLog||obj.debugMod
                    obj.writeLog(txtOk,'date');
                end
            end
        end
    end
end
