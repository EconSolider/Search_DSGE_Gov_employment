cd "C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_est\data\unemployment benifit"
do 失業給付clearing.do
gen n_time=quarterly(time,"YQ")
cap drop time
gen zsum=ub_mean
tempfile ub
save `ub',replace

cd C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_est\data
do 実効税率\実効税率clearing.do
tempfile efficiency_tax
save `efficiency_tax',replace

do 毎月勤労統計調査\毎勤clearing.do
tempfile hours
save `hours',replace

import excel "data.xlsx", sheet("全データ") firstrow clear
gen year=year(time)
gen quarter=quarter(time)
gen n_time=yq(year,quarter)
tsset n_time,quarterly

merge 1:1 n_time using `ub'
cap drop if _merge!=3
cap drop _merge

merge 1:1 year quarter using `hours'
cap drop if _merge!=3
cap drop _merge


*データ整理
local varlist1 "r国内総生産10億 r民間最終消費支出10億 r民間住宅投資10億 r民間企業設備投資10億 r公的固定資本形成10億 r政府最終消費支出公務員人件費含む10億 r雇用者報酬10億 "
local varlist2 "y c1 c2 I GI GC wN"
local length_varlist:word count `varlist1'
forvalue i=1(1)`length_varlist'{
local old_var:word `i' of `varlist1'
local new_var:word `i' of `varlist2'
gen `new_var'=`old_var'
drop `old_var'
local old_var2=subinstr("`old_var'","10億","",.)
label variable `new_var' "`old_var2'"
}
sort n_time

gen pi=GDPデフレーター/L.GDPデフレーター
label variable pi インフレ率デフレーター

rename 名目利子率コールレート i
replace i=1+i/400

rename 完全失業率 u
replace u=u/100

rename 雇用者万人 N
replace N=N*10^4
label variable N 雇用者数

gen w=wN/N
label variable w 賃金

rename 労働力人口万人 pop
replace pop=pop
label variable pop 労働力人口

rename 名目国債残高百万 b
replace b=b/10^3/(GDPデフレーター/100)
label variable b 国債残高

rename 有効求人数原数 V

rename 完全失業者万人 U
replace U=U*10^4
label variable U 完全失業者数

gen z=zsum/GDPデフレーター*100
label variable z 失業給付

*求人部門割合近似
gen xi=V/N
label variable xi 求人部門割合
*就業率
gen n=1-u
label variable n 就業率 

*民間消費
gen c=c1+c2
label variable c 住宅投資_最終消費支出

*per capita
gen c_pc=c/pop
gen I_pc=I/pop
gen GC_pc=GC/pop
gen GI_pc=GI/pop
gen b_pc=b/pop
gen y_pc=c_pc+I_pc+GC_pc+GI_pc

*政府雇用
rename 政府雇用千人 gov_emp_thou
gen gov_emp=gov_emp_thou*1000
reg gov_emp N U V
predict NG,xb
gen ng=NG/(pop*10000)

*実効税率
merge m:1 year using `efficiency_tax'
drop if _merge!=3
cap drop _merge
sort n_time
tsset n_time

gen theta=V/U

export excel using "data_after_clear", firstrow(variables) replace

