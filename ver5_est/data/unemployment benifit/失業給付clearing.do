cd "C:\Users\lzkzz\Desktop\matlearn\Dynare_train\Search_fiscal_bond_correct_generaluti_modularized_est\data\unemployment benifit"
import excel "C:\Users\lzkzz\Desktop\matlearn\Dynare_train\Search_fiscal_bond_correct_generaluti_modularized_est\data\unemployment benifit\基本手当受給人員数.xlsx", sheet("受給者とデータ") firstrow clear

format time %td
gen n_time=ym(year(time),month(time))
gen quarter=quarter(time)
gen year=year(time)
gen month=month(time)
format n_time %tm
tsset n_time
reg 名目失業給付千円 失業給付受給人員 L.失業給付受給人員 L2.失業給付受給人員 L3.失業給付受給人員
predict ub,xb
cap drop if year==.
cap drop if ub==.
gen ub_mean=ub/失業給付受給人員
keep quarter year month ub ub_mean
bysort year quarter:egen ub_qua=sum(ub) 
duplicates drop year quarter,force

gen time=""
local N=_N
forvalue i=1(1)`N'{
local a=year[`i']
local b=quarter[`i']
replace time="`a'"+"Q"+"`b'" in `i'
}
keep time ub ub_mean
order time ub
