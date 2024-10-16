%% Script dedicated to the initialisation, the measure and the display of the spent time
%% L.LAURENT -- 15/05/2012 -- luc.laurent@lecnam.net

classdef mesuTime < handle
    properties (Access = private)
        counterTicToc=[];
        tInit=[];
        tFinal=[];
    end
    properties
        tElapsed=[];
        tCPUElapsed=[];
    end
    methods
        %constructor
        function obj=mesuTime()
            %initialization tic-toc count
            obj.counterTicToc=tic;
            %Measure initial CPU time
            obj.tInit=cputime;
        end
        %start the measure
        function start(obj)
            %Measure initial CPU time
            obj.tInit=cputime;
        end
        %stop the measure
        function stop(obj,flagDisp)
            if nargin==1;flagDisp=true;end
            %Measure initial CPU time
            obj.tElapsed=toc(obj.counterTicToc);
            %Elapsed time
            obj.tCPUElapsed=cputime-obj.tInit;
            %show the result
            if flagDisp
                show(obj);
            end
        end
        %get CPU time elasped
        function timeElapsed=cpuElapsed(obj)
            if isempty(obj.tCPUElapsed);obj.stop(false);end
            timeElapsed=obj.tCPUElapsed;
        end
        %get time elasped
        function timeElapsed=Elapsed(obj)
            if isempty(obj.tCPUElapsed);obj.stop(false);end
            timeElapsed=obj.tElapsed;
        end
        %display the elapsed time
        function show(obj)
            fprintf('  #### Time/CPU time (s): %4.2f s / %4.2f s\n',obj.tElapsed,obj.tCPUElapsed);
        end
    end
end
