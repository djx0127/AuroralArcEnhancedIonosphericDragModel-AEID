function [Time_day,time1,sma1,ecc1,inc1,OMG1,aop1,M1,r1] = TLEelments_Mutation(DeleteTime,time,julian_time,sma,ecc,inc,OMG,aop,M,r)
% Remove the time periods based on the time point of orbit maintenance control.
stage_delete = size(DeleteTime,1);
Time_day = julian_time-julian_time(1);
% 创建逻辑索引，标记需要保留的数据点
keep_mask = true(length(julian_time), 1);
if ~isempty(DeleteTime)
    for idx = 1:stage_delete
        % 标记需要移除的数据点
        remove_mask = (Time_day >= DeleteTime(idx,1)) & (Time_day <= DeleteTime(idx,2));
        keep_mask(remove_mask) = false;
    end
end

Time_day = Time_day(keep_mask);
time1 = time(keep_mask);
sma1 = sma(keep_mask);
ecc1 = ecc(keep_mask);
inc1 = inc(keep_mask);
OMG1 = OMG(keep_mask);
aop1 = aop(keep_mask);
M1   = M(keep_mask);
r1   = r(keep_mask);

end