/*bootstrap法求 CI 及其95%置信区间*/
/*%let indata=a;*/
/*%let xvar=chronic1;*/
/*%let wealth=prin1;*/
/**/
/*%CI(indata=final2013,xvar=chronic1,wealth=prin1,outdata=out);*/

%macro CI1(indata=,xvar=,wealth=,outdata=);
/* 第一步*/
data _null_;
  set &indata nobs=nobs end=end;
  if end;
  call symputx("number",nobs);
run;

/* 第二步:对家庭财富求分数秩*/
proc rank data =&indata out = &xvar.ranka; /* 输出数据集为&xvar.ranka*/
var &wealth; /* 根据抽样号分组分别求家庭财富的秩*/
ranks wealthrank; /* 家庭财富的秩表示为wealthrank*/
run;
proc means data = &xvar.ranka noprint; /*对数据集&xvar.ranka进行统计描述*/
var wealthrank;
/* 根据抽样号分组分别求家庭财富的最大秩wealthrank_max*/
output out = wealthrank_max max = wealthrank_max;
run;
/* 家庭财富的最大秩wealthrank _ max，输出数据集wealthrank_max*/
data &xvar.rank; /* 建立数据集&xvar.rank */
set &xvar.ranka ;
if _n_ = 1 then set wealthrank_max; /* 合并数据集&xvar.ranka与wealthrank_max*/
wealthrank_frac = wealthrank/wealthrank_max; /* 求出家庭财富的分数秩*/
run;

/* 第三步: 求结局变量与分数秩的协方差*/
proc corr data = &xvar.rank COV outp=cov(type=cov) noprint;
var &xvar wealthrank_frac; /* 求出家庭财富分数秩与结局变量之间的协方差*/
run;
data cova; set cov; /* 建立数据集cova*/
if _name_ ="wealthrank_frac" and _type_="COV"; cov = &xvar; /* 家庭财富分数秩与&xvar. 之间的协方差表示为cov*/
keep  cov; /* 数据集cova 只保留变量replicate cov*/
run;

/* 第四步: 求集中指数*/
proc means data = &indata noprint;
var &xvar; output out = mu mean = mu; /* 求出各组中&xvar. 的均数mu 并保存为数据集mu*/
run;
data mua; set mu; keep mu; /* 保留变量 mu 建立数据集mua*/
run;
data ci; merge cova mua; /* 合并协方差数据集cova 与数据集mua，建立数据集ci*/
CI = 2* COV/mu; /* 根据公式求出各组的集中指数*/
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
/*1.Probit模型*/
ods listing close;
ods output  ParameterEstimates=param;
proc genmod   data=&data  DESCENDING ;
   model &y=&x;
run;
data param1;
  set param;
  if  _n_^=1;
  if Parameter="尺度" then delete;
  keep Parameter Estimate StdErr ProbChiSq;
run;

/*2.自变量CI*/
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
/*3.应变量CI*/
%CI1(indata=all,xvar=&y,wealth=&w,outdata=yci)
data _null_;
  set yci;
  call symput("ciy",ci);
run;

/*4.计算自变量均值*/
proc means data=&data;
  var &x;
  output out=meanx mean=;
run;
proc transpose data=meanx(drop=_TYPE_ _FREQ_) out=meanx1(rename=COL1=mux keep=col1);
run;
/*5.计算应变量均值*/
proc means data=&data;
  var &y;
  output out=meany mean=;
run;
data _null_;
  set meany;
  call symput("muy",&y);
run;
/*5.计算自变量贡献和百分比*/
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
/*5.计算残差贡献和百分比*/
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

