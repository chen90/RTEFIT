function A = read_afile(filenamestr)
%% read the EFIT output data from afile
% -inputs:
% filenamestr   string, the full path of the afile
% -outputs:
% A             struct, the data passed from afile

% Edited by Shuying SUN on 2019/05/28
% Contact: bangrsun@163.com, sunshuyingc@enn.cn
% ENN Group 1989-2019, all rights reserved.

% Updated by Shuying SUN on 2019/05/29, comments after reading each line
% are from EFIT source file weqdud129.f90

fid=fopen(filenamestr);
if(fid<=0)
    error(['Can not open file:', filenamestr]);
end
A={};
% read first 3 lines, but do not handle the information
fgetl(fid);
%% read 2nd line
A.ishot =fscanf(fid,'%d',1);
A.ktime1=fscanf(fid,'%d\n',1);
% if (ishot.le.99999) then
% write (neqdsk,1050) ishot,ktime1
% else
% write (neqdsk,1053) ishot,ktime1
% endif
%% read 3th line
A.time=fscanf(fid,'%g\n',1);
%write (neqdsk,1040) (time(j),j=i,i)
%% read 4th line
temp=fscanf(fid,'%c',1);
if(temp~='*')
    error('Wrong afile format!');
end
A.time=fscanf(fid,'%g',1);
fscanf(fid,'%c',9); % 9 space
A.jflag=fscanf(fid,'%5d',1);
fscanf(fid,'%c',11); % 11 space
A.lflag=fscanf(fid,'%5d',1);
fscanf(fid,'%c',1);  % 1 space
A.limloc=fscanf(fid,'%c',3);
fscanf(fid,'%c',1);  % 1 space
A.mco2v=fscanf(fid,'%3d',1);
fscanf(fid,'%c',1);  % 1 space
A.mco2r=fscanf(fid,'%3d',1);
fscanf(fid,'%c',1);  % 1 space
A.qmflag=fscanf(fid,'%c',3);
fscanf(fid,'%c',1);  % 1 space
A.nlold=fscanf(fid,'%5d',1);
A.nlnew=fscanf(fid,'%5d',1);
% end of 4th line
%write (neqdsk,1060) time(jj),jflag(jj),lflag,limloc(jj),mco2v,mco2r,qmflag,nlold,nlnew
%1060 format (1h*,f8.3,9x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5)
%% read 5th line
temp=fscanf(fid,'%e',4);
A.tsaisq=temp(1); A.rcencm=temp(2); A.bcentr=temp(3); A.pasmat=temp(4);
%write (neqdsk,1040) tsaisq(jj),rcencm,bcentr(jj),pasmat(jj)
%% read 6th line
temp=fscanf(fid,'%e',4);
A.cpasma=temp(1); A.rout=temp(2); A.zout=temp(3); A.aout=temp(4);
%write (neqdsk,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)
%% read 7th line
temp=fscanf(fid,'%e',4);
A.eout=temp(1); A.doutu=temp(2); A.doutl=temp(3); A.vout=temp(4);
%write (neqdsk,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)
%% read 8th line
temp=fscanf(fid,'%e',4);
A.rcurrt=temp(1); A.zcurrt=temp(2); A.qsta=temp(3); A.betat=temp(4);
%write (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
%% read 9th line
temp=fscanf(fid,'%e',4);
A.betap=temp(1); A.ali=temp(2); A.oleft=temp(3); A.oright=temp(4);
%write (neqdsk,1040) betap(jj),ali(jj),oleft(jj),oright(jj)
%% read 10th line
temp=fscanf(fid,'%e',4);
A.otop=temp(1); A.obott=temp(2); A.qpsib=temp(3); A.vertn=temp(4);
%write (neqdsk,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)
%% read 11th line
A.rco2v=fscanf(fid,'%e',A.mco2v);
%write (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
%% read 12th line
A.dco2v=fscanf(fid,'%e',A.mco2v);
%write (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
%% read 13th line
A.rco2r=fscanf(fid,'%e',A.mco2r);
%write (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
%% read 14th line
A.dco2r=fscanf(fid,'%e',A.mco2r);
%write (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
%% read 15th line
temp=fscanf(fid,'%e',4);
A.shearb=temp(1); A.bpolav=temp(2); A.s1=temp(3); A.s2=temp(4);
%write (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
%% read 16th line
temp=fscanf(fid,'%e',4);
A.s3=temp(1); A.qout=temp(2); A.olefs=temp(3); A.orighs=temp(4);
%write (neqdsk,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)
%% read 17th line
temp=fscanf(fid,'%e',4);
A.otops=temp(1); A.sibdry=temp(2); A.areao=temp(3); A.wplasm=temp(4);
%write (neqdsk,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
%% read 18th line
temp=fscanf(fid,'%e',4);
A.terror=temp(1); A.elongm=temp(2); A.qqmagx=temp(3); A.cdflux=temp(4);
%write (neqdsk,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj)
%% read 19th line
temp=fscanf(fid,'%e',4);
A.alpha=temp(1); A.rttt=temp(2); A.psiref=temp(3); A.xndnt=temp(4);
%write (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
%% read 20th line
temp=fscanf(fid,'%e',4);
A.rseps(1)=temp(1); A.zseps(1)=temp(2); A.rseps(2)=temp(3); A.zseps(2)=temp(4);
%write (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj),zseps(2,jj)
%% read 21st line
temp=fscanf(fid,'%e',4);
A.sepexp=temp(1); A.obots=temp(2); A.btaxp=temp(3); A.btaxv=temp(4);
%write (neqdsk,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
%% read 22nd line
temp=fscanf(fid,'%e',4);
A.aaq1=temp(1); A.aaq2=temp(2); A.aaq3=temp(3); A.seplim=temp(4);
%write (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
%% read 23th line
temp=fscanf(fid,'%e',4);
A.rmagx=temp(1); A.zmagx=temp(2); A.simagx=temp(3); A.taumhd=temp(4);
%write (neqdsk,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
%% read 24th line
temp=fscanf(fid,'%e',4);
A.betapd=temp(1); A.betatd=temp(2); A.wplasmd=temp(3); A.fluxx=temp(4);
%write (neqdsk,1040) betapd(jj),betatd(jj),wplasmd(jj),fluxx
%% read 25th line
temp=fscanf(fid,'%e',4);
A.vloopt=temp(1); A.taudia=temp(2); A.qmerci=temp(3); A.tavem=temp(4);
%write (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem
%% read 26th line
temp=fscanf(fid,'%d',4);
A.nsilop=temp(1); A.magpri=temp(2); A.nfcoil=temp(3); A.nesum=temp(4);
%write (neqdsk,1041) nsilop0,magpri0,nfcoil0,nesum0
%%
A.csilop=fscanf(fid,'%e',A.nsilop);
A.cmpr2 =fscanf(fid,'%e',A.magpri);
A.ccbrsp=fscanf(fid,'%e',A.nfcoil);
A.eccurt=fscanf(fid,'%e',A.nesum);
%write (neqdsk,1040) (csilop(k,jj),k=1,nsilop),(cmpr2(k,jj),k=1,magpri)
%write (neqdsk,1040) (ccbrsp(k,jj),k=1,nfcoil)
%write (neqdsk,1040) (eccurt(jj,k),k=1,nesum)
%%
temp=fscanf(fid,'%e',4);
A.pbinj=temp(1); A.rvsin=temp(2); A.zvsin=temp(3); A.rvsout=temp(4);
%write (neqdsk,1040) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
%%
temp=fscanf(fid,'%e',4);
A.zvsout=temp(1); A.vsurfa=temp(2); A.wpdot=temp(3); A.wbdot=temp(4);
%write (neqdsk,1040) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
%%
temp=fscanf(fid,'%e',4);
A.slantu=temp(1); A.slantl=temp(2); A.zuperts=temp(3); A.chipre=temp(4);
%write (neqdsk,1040) slantu(jj),slantl(jj),zuperts(jj),chipre
%%
temp=fscanf(fid,'%e',4);
A.cjor95=temp(1); A.pp95=temp(2); A.ssep=temp(3); A.yyy2=temp(4);
%write (neqdsk,1040) cjor95(jj),pp95(jj),ssep(jj),yyy2(jj)
%%
temp=fscanf(fid,'%e',4);
A.xnnc=temp(1); A.cprof=temp(2); A.oring=temp(3); A.cjor0=temp(4);
%write (neqdsk,1040) xnnc(jj),cprof,oring(jj),cjor0(jj)
%%
temp=fscanf(fid,'%e',4);
A.fexpan=temp(1); A.qqmin=temp(2); A.chigamt=temp(3); A.ssi01=temp(4);
%write (neqdsk,1040) fexpan,qqmin,chigamt,ssi01
%%
temp=fscanf(fid,'%e',4);
A.fexpvs=temp(1); A.sepnose=temp(2); A.ssi95=temp(3); A.rqqmin=temp(4);
%write (neqdsk,1040) fexpvs,sepnose,ssi95(jj),rqqmin
%%
temp=fscanf(fid,'%e',4);
A.cjor99=temp(1); A.cj1ave=temp(2); A.rmidin=temp(3); A.rmidout=temp(4);
%write (neqdsk,1040) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
%%
temp=fscanf(fid,'%e',4);
A.psurfa=temp(1); A.peak=temp(2); A.dminux=temp(3); A.dminlx=temp(4);
%write (neqdsk,1040) psurfa(jj), peak(jj),dminux(jj),dminlx(jj)
%%
temp=fscanf(fid,'%e',4);
A.dolubaf=temp(1); A.dolubafm=temp(2); A.diludom=temp(3); A.diludomm=temp(4);
%write (neqdsk,1040) dolubaf(jj),dolubafm(jj),diludom(jj),diludomm(jj)
%%
temp=fscanf(fid,'%e',4);
A.ratsol=temp(1); A.rvsiu=temp(2); A.zvsiu=temp(3); A.rvsid=temp(4);
%write (neqdsk,1040) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
%%
temp=fscanf(fid,'%e',4);
A.zvsid=temp(1); A.rvsou=temp(2); A.zvsou=temp(3); A.rvsod=temp(4);
%write (neqdsk,1040) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
%%
temp=fscanf(fid,'%e',4);
A.zvsod=temp(1); A.condno=temp(2); A.psin32=temp(3); A.psin21=temp(4);
%write (neqdsk,1040) zvsod(jj),condno,psin32(jj),psin21(jj)
%%
temp=fscanf(fid,'%e',4);
A.rq32in=temp(1); A.rq21top=temp(2); A.chilibt=temp(3); A.ali3=temp(4);
%write (neqdsk,1040) rq32in(jj),rq21top(jj),chilibt,ali3(jj)
%%
temp=fscanf(fid,'%e',4);
A.xbetapr=temp(1); A.tflux=temp(2); A.xdum=temp(3); A.xdum=temp(4);
%write (neqdsk,1040) xbetapr,tflux(jj),xdum,xdum
%%
fscanf(fid,'%c',2); % include a LINE-FEED character before this space
A.header=fscanf(fid,'%c',42);
fscanf(fid,'%c',1);
A.fit_type=fscanf(fid,'%c',3);
%write (neqdsk,1042) header,fit_type
%1042 format (1x,a42,1x,a3)
%%
fclose(fid);

end