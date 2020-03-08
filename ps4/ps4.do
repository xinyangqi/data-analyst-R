// STATS 506 Fall 2019
// This script contains solutions for problem set 4
//
// The question 1 create four models to calculate the relative effect of condition
// for four variables and 95% CI.
// The question 2 calculates the census divisions' disparity between urban and rural
// areas in terms of the proportion of homes with internet access.
// The question 3 creates two models.
// The first model using only day 1 data , fits a logistic regression investigating how the 
// odds of drinking water are change for weekdays relative to weekends, while controlling 
// for "winter", age and age square, gender, and pir.
// The second model fit a mixed logistic model, using data of both days and including a 
// random intercept for each respondent. 

// Author : Xinyang Qi 
// Updated : Nov 18 , 2019


// Document the process
version 16.0
log using ps4.txt , text replace

// Question 1
// import data
import delimited M:/506/mouse.csv
egen group_con=group(condition)

// create model log(tot_dist) ~ condition
generate log_tot_dist=log(tot_dist)
mixed log_tot_dist group_con || _all: R.subject_nr || _all: R.exemplar
matrix model1=r(table)

// create model log(max_abs_dev) ~ condition
generate log_max_abs_dev=log(max_abs_dev)
mixed log_max_abs_dev group_con || _all: R.subject_nr || _all: R.exemplar
matrix model2=r(table)

// create model log(avg_abs_dev) ~ condition
generate log_avg_abs_dev=log(avg_abs_dev)
mixed log_avg_abs_dev group_con || _all: R.subject_nr || _all: R.exemplar
matrix model3=r(table)

// create model log(auc) ~ condition
generate log_auc=log(auc)
mixed log_auc group_con || _all: R.subject_nr || _all: R.exemplar
matrix model4=r(table)

// export relative effect and CI(95%)
mata
model1=st_matrix("model1")
model2=st_matrix("model2")
model3=st_matrix("model3")
model4=st_matrix("model4")
model1=round(model1[(1,5,6),1],0.001)
model2=round(model2[(1,5,6),1],0.001)
model3=round(model3[(1,5,6),1],0.001)
model4=round(model4[(1,5,6),1],0.001)
st_matrix("s_model1",model1)
st_matrix("s_model2",model2)
st_matrix("s_model3",model3)
st_matrix("s_model4",model4)
end

putexcel set ps4_1.xls , replace
putexcel A2="relative effect"
putexcel A3="lwr"
putexcel A4="upr"
putexcel B1="tot_dist"
putexcel C1="max_abs_dev"
putexcel D1="avg_abs_dev"
putexcel E1="AUC"
putexcel B2=matrix(s_model1)
putexcel C2=matrix(s_model2)
putexcel D2=matrix(s_model3)
putexcel E2=matrix(s_model4)
capture clear

// Question2
// import data
import delimited M:/506/recs2015_public_v4.csv
generate area=uatyp10=="R"
// select variables
keep doeid area nweight internet division brrwt*
preserve
// calculate proportion
generate wei_inter=internet*nweight
collapse (sum) wei_inter nweight , by(area division)
generate propor=wei_inter/nweight
drop wei_inter nweight
reshape wide propor , i(division) j(area)
generate disparsity=abs(propor1-propor0)
drop propor0 propor1
save q2_1.dta , replace
restore
drop nweight
// caculate replicate proportion
reshape long brrwt , i(doeid) j(replicate)
generate wei_inter_repl=internet*brrwt
collapse (sum) wei_inter_repl brrwt , by(replicate division area)
generate propor_repl=wei_inter_repl/brrwt
drop wei_inter_repl brrwt
reshape wide propor_repl , i(division replicate) j(area)
generate disparity_repl=abs(propor_repl1-propor_repl0)
drop propor_repl1 propor_repl0
// join data to calculate CI
merge m:1 division using q2_1.dta
drop _merge
generate square=(disparity_repl-disparsity)^2
collapse (first) disparsity (sum) square , by(division)
generate se=sqrt(4/96*square)
drop square
generate upr=disparsity+1.96*se
generate lwr=disparsity-1.96*se
gsort -disparsity
// export result
export delimited ps4_q2.csv , replace

// Question 3
// (a)
// Handle the information data
fdause M:/506/DEMO_D.XPT , clear
keep seqn riagendr ridageyr indfmpir ridexmon
save MEMO.dta , replace

// Handle the day1 data
fdause M:/506/DR1TOT_D.XPT
keep seqn dr1day dr1_320 dr1_330 dr1bwatz
generate any_prewater1=0
replace any_prewater1=1 if dr1_320>0 | dr1_330>0 | dr1bwatz>0
replace any_prewater1=. if dr1_320==. & dr1_330==. & dr1bwatz==.
drop dr1_320 dr1_330 dr1bwatz
generate is_weekday1=dr1day
replace is_weekday1=0 if is_weekday1==1 | is_weekday1==7
replace is_weekday1=1 if is_weekday1!=0 & is_weekday1!=.
rename dr1day drday1
save day1.dta , replace

// Handle the day2 data
fdause M:/506/DR2TOT_D.XPT , clear
keep seqn dr2day dr2_320 dr2_330 dr2bwatz
generate any_prewater2=0
replace any_prewater2=1 if dr2_320>0 | dr2_330>0 | dr2bwatz>0
replace any_prewater2=. if dr2_320==. & dr2_330==. & dr2bwatz==.
drop dr2_320 dr2_330 dr2bwatz
generate is_weekday2=dr2day
replace is_weekday2=0 if is_weekday2==1 | is_weekday2==7
replace is_weekday2=1 if is_weekday2!=0 & is_weekday2!=.
rename dr2day drday2
save day2.dta , replace

// join three data frames
merge 1:1 seqn using day1.dta 
drop _merge
// transformate to a long table
reshape long drday any_prewater is_weekday , i(seqn) j(survey_day)
merge m:1 seqn using MEMO.dta 
replace ridexmon=0 if ridexmon==2
rename ridexmon is_winter
drop _merge
save q3_a.dta , replace

// (b)

generate missing=0
foreach var of varlist * {
replace missing=1 if `var'==.
}
keep if missing==0
summarize(ridageyr) 
// The mean of ridageyr is 27.97578
generate cen_age=ridageyr-27.97578
summarize(indfmpir)
// The mean of indfmpir is 2.395235
generate cen_pir=indfmpir-2.395235
drop ridageyr indfmpir missing
replace cen_age=cen_age/10
save q3_b.dta , replace

// (c)
//create logistic model
logistic any_prewater i.is_weekday i.is_winter i.riagendr c.cen_pir c.cen_age c.cen_age#c.cen_age if survey_day==1
// extract useful information of the model
matrix mod=r(table)
margins , dydx(*)
matrix mar=r(table)
mata
mod=st_matrix("mod")
mod=round(mod[(1,5,6),(2,4,6,7,8,9)],0.001)
mar=st_matrix("mar")
mar=round(mar[(1,5,6),(2,4,6,7,8)],0.001)
st_matrix("s_mod",mod)
st_matrix("s_mar",mar)
end

// export the result to an excel
putexcel set ps4_3c.xls , replace
putexcel A2="odds ratio"
putexcel A3="odds lwr"
putexcel A4="odds upr"
putexcel A5="margin effect"
putexcel A6="ME lwr"
putexcel A7="ME upr"
putexcel B1="is_weekday"
putexcel C1="is_winter"
putexcel D1="gender"
putexcel E1="cen_pir"
putexcel F1="cen_age"
putexcel G1="cen_age^2"
putexcel B2=matrix(s_mod)
putexcel B5=matrix(s_mar)

// (d)
// create a glm model
meglm any_prewater i.is_weekday i.is_winter i.riagendr c.cen_pir c.cen_age c.cen_age#c.cen_age || seqn: , family(binomial) link(logit) or
// extract useful information of the model
matrix gmod=r(table)
margins , dydx(*)
matrix gmar=r(table)
mata
gmod=st_matrix("gmod")
gmod=round(gmod[(1,5,6),(2,4,6,7,8,9)],0.001)
gmar=st_matrix("gmar")
gmar=round(gmar[(1,5,6),(2,4,6,7,8)],0.001)
st_matrix("s_gmod",gmod)
st_matrix("s_gmar",gmar)
end

// export the result to an excel
putexcel set ps4_3d.xls , replace
putexcel A2="odds ratio"
putexcel A3="odds lwr"
putexcel A4="odds upr"
putexcel A5="margin effect"
putexcel A6="ME lwr"
putexcel A7="ME upr"
putexcel B1="is_weekday"
putexcel C1="is_winter"
putexcel D1="gender"
putexcel E1="cen_pir"
putexcel F1="cen_age"
putexcel G1="cen_age^2"
putexcel B2=matrix(s_gmod)
putexcel B5=matrix(s_gmar)

// From the output above, we can find that the two models are almost the same.











