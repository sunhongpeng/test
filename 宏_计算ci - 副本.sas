/*bootstrap���� CI ����95%��������*/
/*%let indata=a;*/
/*%let xvar=chronic1;*/
/*%let wealth=prin1;*/
/**/
/*%CI(indata=final2013,xvar=chronic1,wealth=prin1,outdata=out);*/

%macro CI1(indata=,xvar=,wealth=,outdata=);
/* ��һ��*/
data _null_;
  set &indata nobs=nobs end=end;
  if end;
  call symputx("number",nobs);
run;

/* �ڶ���:�Լ�ͥ�Ƹ��������*/
proc rank data =&indata out = &xvar.ranka; /* ������ݼ�Ϊ&xvar.ranka*/
var &wealth; /* ���ݳ����ŷ���ֱ����ͥ�Ƹ�����*/
ranks wealthrank; /* ��ͥ�Ƹ����ȱ�ʾΪwealthrank*/
run;
proc means data = &xvar.ranka noprint; /*�����ݼ�&xvar.ranka����ͳ������*/
var wealthrank;
/* ���ݳ����ŷ���ֱ����ͥ�Ƹ��������wealthrank_max*/
output out = wealthrank_max max = wealthrank_max;
run;
/* ��ͥ�Ƹ��������wealthrank _ max��������ݼ�wealthrank_max*/
data &xvar.rank; /* �������ݼ�&xvar.rank */
set &xvar.ranka ;
if _n_ = 1 then set wealthrank_max; /* �ϲ����ݼ�&xvar.ranka��wealthrank_max*/
wealthrank_frac = wealthrank/wealthrank_max; /* �����ͥ�Ƹ��ķ�����*/
run;

/* ������: ���ֱ���������ȵ�Э����*/
proc corr data = &xvar.rank COV outp=cov(type=cov) noprint;
var &xvar wealthrank_frac; /* �����ͥ�Ƹ����������ֱ���֮���Э����*/
run;
data cova; set cov; /* �������ݼ�cova*/
if _name_ ="wealthrank_frac" and _type_="COV"; cov = &xvar; /* ��ͥ�Ƹ���������&xvar. ֮���Э�����ʾΪcov*/
keep  cov; /* ���ݼ�cova ֻ��������replicate cov*/
run;

/* ���Ĳ�: ����ָ��*/
proc means data = &indata noprint;
var &xvar; output out = mu mean = mu; /* ���������&xvar. �ľ���mu ������Ϊ���ݼ�mu*/
run;
data mua; set mu; keep mu; /* �������� mu �������ݼ�mua*/
run;
data ci; merge cova mua; /* �ϲ�Э�������ݼ�cova �����ݼ�mua���������ݼ�ci*/
CI = 2* COV/mu; /* ���ݹ�ʽ�������ļ���ָ��*/
run;

data &outdata;
  length xvar $ 15;
  xvar="&xvar";
  n=&number;
  set ci(keep=ci);
run;

proc datasets lib=work noprint;
  delete &xvar.: ci cov cova mu mua wealthrank_max ;
run;
quit;
%mend;





%macro decilogit(data=,y=,w=,x=,out=);
/*1.Probitģ��*/
ods listing close;
ods output  ParameterEstimates=param;
proc genmod   data=&data  DESCENDING ;
   model &y=&x;
run;
data param1;
  set param;
  if  _n_^=1;
  if Parameter="�߶�" then delete;
  keep Parameter Estimate StdErr ProbChiSq;
run;

/*2.�Ա���CI*/
data xci;
  set _null_;
run;
  %let i=1;
  %do %while(%scan(&x,&i,%str( ))  ne %str());
    %let xx = %scan(&x,&i,%str( )); 
	%let i = %eval(&i+1);
    %CI1(indata=&data,xvar=&xx,wealth=&w,outdata=xci1)
	data xci;
	  set xci xci1;
	run;
   %end;
/*3.Ӧ����CI*/
%CI1(indata=all,xvar=&y,wealth=&w,outdata=yci)
data _null_;
  set yci;
  call symput("ciy",ci);
run;

/*4.�����Ա�����ֵ*/
proc means data=&data;
  var &x;
  output out=meanx mean=;
run;
proc transpose data=meanx(drop=_TYPE_ _FREQ_) out=meanx1(rename=COL1=mux keep=col1);
run;
/*5.����Ӧ������ֵ*/
proc means data=&data;
  var &y;
  output out=meany mean=;
run;
data _null_;
  set meany;
  call symput("muy",&y);
run;
/*5.�����Ա������׺Ͱٷֱ�*/
data outa;
  merge param1 meanx1 xci(keep=ci) end=end;
  muy=symgetn('muy');
  Elasti=Estimate*mux/muy;
  Contri=Elasti*ci;
  sumci+Contri;
  if end=1 then call symput("sumci",sumci);
  ciy=symgetn('ciy');
  Percent=Contri/ciy*100;
/*  drop sumci muy mux ciy;*/
run;
/*5.����в�׺Ͱٷֱ�*/
data resi;
  Parameter="Residual";
  Contri=symgetn('ciy')-symgetn('sumci');
  Percent=Contri/symgetn('ciy')*100;
  output;
  Parameter="Total";
  Contri=symgetn('ciy');
  Percent=Contri/symgetn('ciy')*100;
  output;
run;

data &out;
  set outa resi;
run;

proc datasets lib=work noprint;
  delete Meanx Meanx1 Meany Outa param param1 resi xci xci1 yci;
run;
quit;
%mend;

data all;
  set diabetes.all2011_2;
  array var(8) chronic1 class_income agegroup sex education urban_nbs Uninsured BMI;
  do i=1 to 8;
   if var(i)=. then delete;
   end;
  keep chronic1-chronic14 income agegroup sex education urban_nbs Uninsured BMI  class_income;
run;
%decilogit(data=all,y=chronic8,w=income,
x=sex education urban_nbs Uninsured BMI,
out=www);

