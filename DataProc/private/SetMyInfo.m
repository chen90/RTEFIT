function y=SetMyInfo

DasInfo.FileType='swip_das01';        %10 文件类型,字符串
DasInfo.ChnlId=0;         %2  通道号(整型)
mystr=zeros(1,12);
DasInfo.ChnlName=char(mystr);       %12 通道名,字符串
DasInfo.Addr=0;           %4  数据地址(长整型)
DasInfo.Freq=0;         %4  采数频率(单精度)
DasInfo.Len=0;            %4  数据长度(长整型)
DasInfo.Post=0;           %4  触发后长度(长整型)
DasInfo.MaxDat=0;         %2  满量程时的A/D转换(整型)
DasInfo.LowRang=0;      %4  量程下限(单精度)
DasInfo.HighRang=0;     %4  量程上限(单精度)
DasInfo.Factor=1;       %4  系数因子(单精度)
DasInfo.Offset=0;       %4  信号偏移量(单精度)
mystr=zeros(1,8);
DasInfo.Unit=char(mystr);  %8  物理量单位,字符串
DasInfo.Dly=0;          %4  延迟(单精度)
DasInfo.AttribDt=0;       %2  数据属性(整型)
DasInfo.DatWth=0;         %2  数据宽度(整型)
                                            %Sum=74
                                            %sparing for later use
DasInfo.SparI1=0;         %2  备用字段1 数字 (整型)
DasInfo.SparI2=0;         %2  备用字段2 数字 (整型)
DasInfo.SparI3=0;         %2  备用字段3 数字 (整型)
DasInfo.SparF1=0;       %4  备用字段1 数字 (单精度)
DasInfo.SparF2=0;       %4  备用字段2 数字 (单精度)
mystr=zeros(1,8);
DasInfo.SparC1=char(mystr);         %8  备用字段1 数字 (字符串)
%mystr=zeros(1,16);
%DasInfo.SparC2=char(mystr);         %16 备用字段2 数字 (字符串)
DasInfo.SparC2='Developed by SXM';         %16 备用字段2 数字 (字符串)
mystr=zeros(1,10);
DasInfo.SparC3=char(mystr);         %10 备用字段3 数字 (字符串)
                                            %Total  Sum=122
y=DasInfo;
