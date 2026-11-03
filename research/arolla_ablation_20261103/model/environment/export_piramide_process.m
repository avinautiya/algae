function export_piramide_process(source_mat, output_dir)
% Export ONLY process observations, never load ablation or calibration files.
% Requires MATLAB timetable support. Preserve source clock, do not infer UTC.
s=load(source_mat,'debtemp','sm','trh','wind');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
families={'debtemp','sm','trh','wind'}; meta=struct();
for f=1:numel(families)
    family=families{f}; assert(isfield(s,family));
    for j=1:3
        tower=sprintf('Tower%d',j); tt=s.(family).(tower);
        assert(istimetable(tt),'Source object is not a timetable');
        rt=tt.Properties.RowTimes;
        assert(isdatetime(rt),'Expected native datetime row times');
        assert(~any(isnat(rt)) && issorted(rt) && numel(unique(rt))==numel(rt),'Invalid time axis');
        info=struct('row_count',height(tt),'variable_names',{tt.Properties.VariableNames}, ...
          'variable_units',{tt.Properties.VariableUnits},'source_time_zone',rt.TimeZone, ...
          'clock_authority_verified',~isempty(rt.TimeZone));
        rt.Format='yyyy-MM-dd''T''HH:mm:ss.SSS';
        tab=timetable2table(tt,'ConvertRowTimes',false);
        tab=addvars(tab,string(rt),'Before',1,'NewVariableNames','SOURCE_LOCAL_LABEL');
        filename=sprintf('%s_%s.csv',family,tower);
        writetable(tab,fullfile(output_dir,filename));
        info.file=filename; meta.(family).(tower)=info;
    end
end
fid=fopen(fullfile(output_dir,'SOURCE_CLOCK_AND_COLUMNS.json'),'w');
assert(fid~=-1); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(meta,'PrettyPrint',true));
end
