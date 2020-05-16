/* Paper to study causal effect of Republic of India's federal health scheme 
called "Ayushman Bharat Yojana" which allows citizens below a certain income level 
to get free health insurance and avail health care at private facilities
at no cost */
//Scheme implemented in September 2018
/*Out of 36 states, four did not choose to adopt it, presenting a natural experimental setting
to study causal effect of scheme on health outcomes in treated states. Health outcomes studied are
MMR and IMR using Diff-in-diff and Triple differencing estimators on three year 
unbalanced panel data using twoway state and year fixed effects.*/
/*Wild cluster bootstrap is used to adjust p values and t stats as 
the number of groups is just 36. Clustered standard errors used.*/

//Mrinalini Darswal, md42877

clear

use "C:\Users\Mrinalini\Desktop\Project Data\CausalFinalPaper.dta", replace

//Labelling variables ex. label variable income "Gross income in 2008, in Euro"

label variable state "One of Republic of India’s 36 states/Union Territories"

label variable stateid "State unique identification number"

label variable district "One of India’s 720 administrative units at sub-state level like counties"

label variable distictid "District ID unique identification number"

label variable year "2017 2018 2019"

label variable aby "Treatment Var: Indicator for whether a state adopted the Central Health Benefit "

label variable aby_18 "Interaction term: 1 If a state had the scheme in year 2018”]"

label variable aby_19dd "Interaction term, DD estimator: If a state had the scheme in year 2019”]"

label variable instdel "% Institutional deliveries to Total Reported Deliveries in a year"

label variable livebirths "Total Number of reported live births in a year"

label variable v15 "% live births of Total Reported Births"

label variable lbw "% Newborns having weight less than 2.5 kg to Newborns weighed at birth"

label variable opv "% Newborns given Oral Polio Vaccine (OPV) at birth to Reported live birth"

label variable bcg "% Newborns given BCG (Tuberculosis vaccine) to Reported live birth"

label variable hepb "“% Newborns given Hepatitis-B vaccine (Birth Dose) at birth to Reported live bir"

label variable hepb "% Newborns given Hepatitis-B vaccine (Birth Dose) at birth to Reported live bir"

label variable mmrvac "“% Infants 0 to 11 months old who received Measles, Mumps and Rubella vaccine to"

label variable mmrvac "% Infants 0 to 11 months old who received Measles, Mumps and Rubella vaccine to"

label variable mrboost "No of children 16-24 months age given Measles Rubella Vaccine 2nd dose"

label variable infantdeaths "Total Number of Infant Deaths reported in the year"

label variable motherdeaths "Total no. maternal deaths in the year in total reported births"

label variable totpregs "Total number of pregnant women Registered for Ante Natal Checkups"

label variable mmr "Maternal Mortality Rate: mother deaths per 100,000 live births"

label variable imr "Infant Mortality Rate: number of deaths per 1,000 live births under one year of "

label variable doctor "Total doctors in public heath system at state level averaged for districts"

label variable subcentre "most peripheral health facility in Indian Public Health System"

label variable phc "Primary Health Centre: state-owned rural single doctor facilities"

label variable chc "First referral health units, accepts patients referred from phc’s"

label variable subdivhosp "hospital at sub-district level"

label variable disthosp "final referral centers for the primary and secondary levels of the public health"

label variable totalhlthinfra "Total public health facilities, does not include private clinics and hospitals"

label variable gdp_pc "per capita gdp at current prices"

label variable poor "per capita gdp below median value of 110605 INR (USD 1455, 76 INR=1USD)"

label variable ddd "Triple Difference estimator: “Interaction term:1 If a state had the scheme in ye"

//desc

summarize imr mmr lbw gdp_pc, detail

// installing software

ssc install estout
//ssc install cmogram

/* install program for wild cluster bootstrap : http://cameron.econ.ucdavis.edu/research/papers.html
Stata command for One-way Wild Cluster Bootstrap Robust Standard Errors (with asymptotic refinement)
- Stata user-written command boottest written by the following authors.
- David Roodman, James MacKinnon, Morten Nielsen, Matthew Webb (2018), "Fast and Wild Bootstrap Inference in Stata using boottest".
[https://ideas.repec.org/p/qed/wpaper/1406.html#download]*/

ssc install boottest

/* DD estimate of Infant Mortality rate 
Two-way state time fixed effect, clustered se, without covariates*/

/*#delimit ;
global xlist aby aby_19dd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost mmr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year
;*/

xtset stateid
xi: xtreg imr aby aby_19dd i.year,fe vce (cluster stateid)
eststo IMRFEwoutCov

boottest aby_19dd, reps(999999) seed(42424)

/* DD estimate of Infant Mortality rate 
Two-way state time fixed effect, clustered se, with covariates*/
#delimit ;
xi: xtreg imr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost mmr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe vce (cluster stateid)
;
eststo IMRFEwithCov

boottest aby_19dd, reps(999999) seed(42424)

/* DD estimate of Infant Mortality rate 
Two-way state time random effects, clustered se, with covariates*/

#delimit ;
xi: xtreg imr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost mmr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,re vce (cluster stateid)
;
eststo IMRREwithCov

//Hausman Test for FE vs. Re
qui xtreg imr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb mmrvac mrboost mmr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe 
estimates store fixed

xtreg imr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb mmrvac mrboost mmr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,re 
estimates store random

hausman fixed random

////*****Comparison Table*****\\\\\

#delimit ;
esttab IMRFEwoutCov IMRFEwithCov IMRREwithCov fixed random using "Table 1.rtf",p
title("{\b Table 1.} Two-way state time panel data regressions of Infant Mortality rate (IMR)on Health Scheme (ABY) treatment and other variables")
compress label onecell
mgroups("FEwoutCov" "FEwithCov" "REwithCov" "fsimple" "rsimple" pattern(1 1 1 1 1)) 
addnotes("Notes: This table compares Diff-in-Diff coefficients to study Average Treatment Effect of implementation of 
Central Government Health Insurance Scheme ABY on Reproductive and Child Health outcomes measured by 
IMR in different states of Republic of India. The treatment year is 2018.
The estimates measures DD coefficients using one year each pre(2017) and post treatment(2019),
and their significance across five different regressions. It uses primary official data from Government of India.
IMRFEwoutCov : DD estimate using two-way state time Fixed Effects with clustered se's, without covariates
IMRFEwithCov : DD estimate using Fixed Effects with clustered se's, with covariates
IMRREwithCov : DD estimate using two-way state time Random Effects with clustered se's, with covariates
fsimple: DD estimate using two-way state time Fixed Effects without clustered se's, with covariates
rsimple: DD estimate using two-way state time Random Effects without clustered se's, with covariates
stars ***, **, * denote statistical significance at the 1%, 5% and 10% respectively.") replace
;

//////////*********MMR************\\\\\\\\\\\


/* DD estimate of Maternal Mortality rate 
Two-way state time fixed effect, clustered se, without covariates*/

xtset stateid
xi: xtreg mmr aby aby_19dd i.year,fe vce (cluster stateid)
eststo MMRFEwoutCov

boottest aby_19dd, reps(999999) seed(42424)
/* DD estimate of Maternal Mortality rate 
Two-way state time fixed effect, clustered se, with covariates*/
#delimit ;
xi: xtreg mmr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost imr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe vce (cluster stateid)
;
eststo MMRFEwithCov

boottest aby_19dd, reps(999999) seed(42424)

/* DD estimate of Maternal Mortality rate 
Two-way state time random effects, clustered se, with covariates*/

#delimit ;
xi: xtreg mmr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost imr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,re vce (cluster stateid)
;
eststo MMRREwithCov

//Hausman Test for FE vs. Re
qui xtreg mmr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb mmrvac mrboost imr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe 
estimates store fixedmmr

xtreg mmr aby aby_19dd instdel percentlivebirths lbw opv bcg hepb mmrvac mrboost imr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,re 
estimates store randommmr

hausman fixedmmr randommmr

////****COMPARISON TABLE DD MMR ****\\\\

#delimit ;
esttab MMRFEwoutCov MMRFEwithCov MMRREwithCov fixedmmr randommmr using "Table 2.rtf",p
title("{\b Table 2.} Two-way state time panel data regressions of Infant Mortality rate (IMR)on Health Scheme (ABY) treatment and other variables")
compress label onecell
mgroups("FEwoutCov" "FEwithCov" "REwithCov" "fsimple" "rsimple" pattern(1 1 1 1 1)) 
addnotes("Notes: This table compares Diff-in-Diff coefficients to study Average Treatment Effect of implementation of 
Central Government Health Insurance Scheme ABY on Reproductive and Child Health outcomes measured by 
MMR in different states of Republic of India. The treatment year is 2018.
The estimates measures DD coefficients using one year each pre(2017) and post treatment(2019),
and their significance across five different regressions. It uses primary official data from Government of India.
MMRFEwoutCov : DD estimate using two-way state time Fixed Effects with clustered se's, without covariates
MMRFEwithCov : DD estimate using Fixed Effects with clustered se's, with covariates
MMRREwithCov : DD estimate using two-way state time Random Effects with clustered se's, with covariates
fsimple: DD estimate using two-way state time Fixed Effects without clustered se's, with covariates
rsimple: DD estimate using two-way state time Random Effects without clustered se's, with covariates
stars ***, **, * denote statistical significance at the 1%, 5% and 10% respectively.") replace
;


//////////////////////////******************TRIPLE DIFFERENCING DDD ESTIMATION*******************\\\\\\\\\\\\\\\\\\\\\\\\\\

/* DDD estimate of Infant Mortality rate 
Two-way state time fixed effect, clustered se, without covariates*/

xtset stateid
xi: xtreg imr aby poor ddd i.year,fe vce (cluster stateid)
eststo DDDIMRFEwoutCov

boottest ddd, reps(999999) seed(42424)

/* DDD estimate of Infant Mortality rate 
Two-way state time fixed effect, clustered se, with covariates*/
#delimit ;
xi: xtreg imr aby poor ddd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost mmr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe vce (cluster stateid)
;
eststo DDDIMRFEwithCov

boottest ddd, reps(999999) seed(42424)

////****comparison table DDD IMR****\\\\

#delimit ;
esttab  DDDIMRFEwoutCov DDDIMRFEwithCov using "Table 3.rtf",p
title("{\b Table 3.} Two-way state time panel data DDD regressions of Infant Mortality rate (IMR)on Health Scheme (ABY) treatment and other variables")
compress label onecell
mgroups("FEwoutCov" "FEwithCov" pattern(1 1)) 
addnotes("Notes: This table compares Diff-in-Diff-in-Diff coefficients to study Average Treatment Effect of implementation of 
Central Government Health Insurance Scheme ABY on Reproductive and Child Health outcomes measured by 
IMR in different states of Republic of India. The treatment year is 2018.
The estimates measures DDD coefficients using one year each pre(2017) and post treatment(2019),
and their significance across two different regressions. It uses primary official data from Government of India.
DDDIMRFEwoutCov : DDD estimate using two-way state time Fixed Effects with clustered se's, without covariates
DDDIMRFEwithCov : DDD estimate using Fixed Effects with clustered se's, with covariates
stars ***, **, * denote statistical significance at the 1%, 5% and 10% respectively.") replace
;


/////////////****MMR***\\\\\\\\\\\\\\\

/* DDD estimate of Maternal Mortality rate 
Two-way state time fixed effect, clustered se, without covariates*/

xtset stateid
xi: xtreg mmr aby poor ddd i.year,fe vce (cluster stateid)
eststo DDDMMRFEwoutCov

boottest ddd, reps(999999) seed(42424)

/* DDD estimate of Maternal Mortality rate 
Two-way state time fixed effect, clustered se, with covariates*/
#delimit ;
xi: xtreg mmr aby poor ddd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost imr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe vce (cluster stateid)
;
eststo DDDMMRFEwithCov

boottest ddd, reps(999999) seed(42424)

////****comparison table DDD MMR****\\\\

#delimit ;
esttab  DDDMMRFEwoutCov DDDMMRFEwithCov using "Table 4.rtf",p
title("{\b Table 4.} Two-way state time panel data DDD regressions of Infant Mortality rate (IMR)on Health Scheme (ABY) treatment and other variables")
compress label onecell
mgroups("FEwoutCov" "FEwithCov" pattern(1 1)) 
addnotes("Notes: This table compares Diff-in-Diff-in-Diff coefficients to study Average Treatment Effect of implementation of 
Central Government Health Insurance Scheme ABY on Reproductive and Child Health outcomes measured by 
MMR in different states of Republic of India. The treatment year is 2018.
The estimates measures DDD coefficients using one year each pre(2017) and post treatment(2019),
and their significance across two different regressions. It uses primary official data from Government of India.
DDDMMRFEwoutCov : DDD estimate using two-way state time Fixed Effects with clustered se's, without covariates
DDDMMRFEwithCov : DDD estimate using Fixed Effects with clustered se's, with covariates
stars ***, **, * denote statistical significance at the 1%, 5% and 10% respectively.") replace
;

///////////////////////**********ROBUSTNESS TESTS************\\\\\\\\\\\\\\\\\\\\

//placebo tests//

xtset stateid
xi: xtreg gdp_pc aby aby_19dd i.year,fe vce (cluster stateid)
eststo DDGDPwoutCov

boottest aby_19dd, reps(999999) seed(42424)

xtset stateid
#delimit ;
xi: xtreg gdp_pc mmr aby ddd instdel percentlivebirths lbw opv bcg hepb 
mmrvac mrboost imr doctor subcentre phc chc subdivhosp disthosp gdp_pc i.year,fe vce (cluster stateid)
;
eststo DDGDPwithCov

boottest ddd, reps(999999) seed(42424)

////****comparison table DD GDP per Capita****\\\\

#delimit ;
esttab  DDGDPwoutCov DDGDPFEwithCov using "Table 5.rtf",p
title("{\b Table 5.} Two-way state time panel data Placebo DD regressions of GDP per capita on Health Scheme (ABY) treatment and other variables")
compress label onecell
mgroups("FEwoutCov" "FEwithCov" pattern(1 1)) 
addnotes("Notes: This table compares Placebo Diff-in-Diff coefficients to study Average Treatment Effect of implementation of 
Central Government Health Insurance Scheme ABY on  
GDP per capita in different states of Republic of India. The treatment year is 2018.
The estimates measures DD coefficients using one year each pre(2017) and post treatment(2019),
and their significance across two different regressions. THESE ARE EXPECTED TO BE STATISTICALLY INSIGNIFICANT AT ALL LEVELS OF SIGNIFICANCE.It uses primary official data from Government of India.
DDGDPwoutCov : DD estimate using two-way state time Fixed Effects with clustered se's, without covariates
DDGDPwithCov : DD estimate using Fixed Effects with clustered se's, with covariates
stars ***, **, * denote statistical significance at the 1%, 5% and 10% respectively.") replace
;


















