classdef Test < HandleObject
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        activeTesters;
        % for benchmarks
        tester;
        testNo;
%        testFields;
        testFields={};
        groups=[];
        reportFields;
        IsTextReport=true;
        toExcel=true;
        fileName='c:\temp\report.xls';
        excel;
        ExcelReport;
        sheetName = 'TestsResult';
        currentConfig=0;
        configs;
    end % properties

    methods
        function vec=getReportVec(obj, fieldName, vecSize)
            vec=cell(1,vecSize);
            for i=1:vecSize
                vec{i}=[fieldName '(' num2str(i) ')'];
            end % for
        end


        function result=getReportFields(obj)
            result={};
        end

        % old version. see setTestField
        %function obj=setTestFields (obj,varargin)
        %    obj.testFields={obj.testFields  varargin};
        %end % setFieldValue

        function obj=setTestField(obj, fieldName, values, varargin)
            if nargin>=4 && ischar('values')
                varargin=[values varargin];
                clear values;
            end
            n=length(obj.reportFields);

            p= inputParser;
            p.addParamValue('order', n+1);
            p.addParamValue('groupNo', 0);
            p.parse(varargin{:});
            p=p.Results;
            groupNo=p.groupNo;
            order=p.order;

            putTestField=true;
            if groupNo~=0
                if length(obj.groups)<groupNo
                    obj.groups{groupNo}=length(obj.testFields)/2+1;
                else
                    putTestField = false;
                    obj.groups{groupNo}=[obj.groups{groupNo}  {fieldName, values}];
                end
            end
            if exist('values','var') && putTestField
                obj.testFields=[obj.testFields  {fieldName values}];
            end
            if order<=n
                obj.reportFields{n+1}=obj.reportFields{order};
            end
            %fieldName
            obj.reportFields{order}=eval('fieldName');
       end

       function textReport=runCartesianProduct(obj)
                    textReport=obj.runTests();
        end % runCartesianProduct

         function obj=initTest(obj)
             obj.currentConfig=0;
             obj.configs=obj.makeCartesianProduct();
         end
         % Check which testers were changed
         function updateChangedTesters(obj, c)
            a=obj.getActiveTesters;
            if c>1
                for i=1:length(a)
                    obj.tester{a(i)}.changed=false;
                end
                current=obj.configs.m(:,c);
                previous=obj.configs.m(:,c-1);
                % indexes of changes fields
                fieldIdxs  = find(current-previous);
                fieldNames = obj.testFields(2*fieldIdxs-1);
                for j=1:length(fieldNames)
                    fieldName=fieldNames{j};
                    if strfind(fieldName, 'tester{')==1
                        strcmp(fieldName(1:7) , 'tester{');
                        e=strfind(fieldName,'}')-1;
                        testerNo=str2num(fieldName(8:e));
                        obj.tester{testerNo}.changed=true;
                    end
                end
            else % c==1
                for i=1:length(a)
                    obj.tester{a(i)}.changed=true;
                end
            end % if
         end
         function [obj test]=nextTest(obj)
             current =obj.currentConfig;
             if current < obj.configs.nCols
                current=current+1;
                obj.currentConfig=current;
                c=obj.configs.m(:,current);
                obj.updateChangedTesters(current);
                test=obj;
                %debug=arrayfun(@(x) {x  obj.testFields{x} obj.testFields{x+1}}, 1:2:length(obj.testFields),'UniformOutput',false)
                for fieldNo=1:length(obj.testFields)/2
                    fieldName=obj.testFields{2*fieldNo-1};
                    fieldVec=obj.testFields{2*fieldNo};
                    if iscell(fieldVec)
                        fieldValue=fieldVec{c(fieldNo)};
                    else
                        fieldValue=fieldVec(c(fieldNo));
                    end
                    test=Utils.setFieldValue(test, fieldName, fieldValue);
                end % for fieldName
                for groupNo=1:length(obj.groups)
                    g=obj.groups{groupNo};
                    mainFieldNo=g{1};
                    g=g(2:end);
                    for fieldNo=1:length(g)/2
                        fieldName=g{2*fieldNo-1};
                        fieldVec=g{2*fieldNo};
                        fieldValue = fieldVec(c(mainFieldNo));
                        test=Utils.setFieldValue(test, fieldName, fieldValue);
                    end
                end
             else
                 test=[];
             end % if
          end % function nextTest

          function a=getActiveTesters(obj)
            if isempty(obj.activeTesters)
                a=1:length(obj.tester);
            else
                a=obj.activeTesters;
            end
          end
          function obj=init(obj)
            a=obj.getActiveTesters;
            for t=1:length(a)
                obj=obj.tester{a(t)}.init(obj, a(t));
            end
          end % init
          function configs=makeCartesianProduct(obj)
            nFields=length(obj.testFields)/2;
            sizeCell=cell(1,nFields);
            for fieldNo=1:nFields
               valuesVec=obj.testFields{2*fieldNo};
               sizeCell{fieldNo}=1:length(valuesVec);
            end
            configs=Utils.setProd(sizeCell)';
          end % makeCartesian
        function obj=runAll(obj)
            a=obj.getActiveTesters;
            if isempty(a)
                obj.run();
            else
                for t=1:length(a)
                        tr=obj.tester{a(t)};
                        if tr.changed || tr.runIfUnchanged
                            tr.run();
                        end
                end
            end
        end
        % return set of points of another tester
        function P=previousP(obj, testerNo)
            P=obj.tester{testerNo-1}.P;
        end
        function obj=setTestFieldFromTester(obj, testerNo, fieldname, varargin)
            str=['tester{' num2str(testerNo) '}.'];
            obj=obj.setTestField([str fieldname], varargin{:});
        end
        function result=globalTimeRatio(obj)
            result=obj.tester{1}.totalTime/obj.tester(2).totalTime;
        end
        function result=globalSizeRatio(obj)
            result=obj.tester(1).size/obj.tester(2).size;
        end

        function textReport=runTests(obj)
            tic;
            obj=obj.open();
            obj=initTest(obj);
            obj.print(obj.getTitle);
            [obj test]=obj.nextTest();
            while not (isempty(test))
                 test=test.runAll();
                 reportLine=test.report2Line();
                 obj.print(reportLine);
                 %disp(['test no' num2str(obj.currentConfig) ' over']);
                 [obj test]=obj.nextTest();
             end % while
            obj.close();
            textReport=[];
            toc;
        end % runTests
        function print(obj, reportLine)
            if obj.toExcel
                obj.excel.addRows (reportLine);
            else
                cell2csv( obj.fileName, reportLine, ';');
                %disp(reportLine)
            end % if
        end % print(reportLine)
        function obj=open(obj)
            if obj.toExcel
                 obj.excel=Excel(obj.fileName, obj.sheetName); % use current directory
            end % if
        end
        function close(obj)
            if obj.toExcel
                obj.excel.closeFile();
            end
        end
        function reportLine=makeReportLine(obj)
           reportLine=cell(1,length(obj.reportFields));
           for i=1:length(obj.reportFields)
               if ~isempty(obj.reportFields{i})
                fieldValue = Utils.getFieldValue(obj, obj.reportFields{i});
                fieldValue=Utils.asString(fieldValue);
                reportLine{i}=fieldValue;
               else
                   reportLine{i}='';
               end
           end % for
        end

        function reportLine = report2Line (obj)
            reportLine=obj.makeReportLine();
        end
        function text=getTitle(obj)
            obj=obj.setReportFields();
            text=obj.reportFields;
        end % getTitle
    % for back compatability only
    function obj=setReportFields(obj)
    end % setReportFields
    end % methods
end

