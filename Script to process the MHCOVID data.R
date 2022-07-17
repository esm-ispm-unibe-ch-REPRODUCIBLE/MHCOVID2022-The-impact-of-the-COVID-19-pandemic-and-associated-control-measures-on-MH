### PROCESS THE ORIGINAL DATA ---- 
outcomes<-outcomes %>% rename(
  ESI=EconomicSupportIndex,
  CHIndex=ContainmentHealthIndex,
  School_closing=C1,
  Workplace_closing=C2,
  Stay_home_req=C6,
  Restrictions_movement=C7,
  Facial_cover=H6
)


## format the time variables
outcomes$timepoint<- as.Date(outcomes$timepoint,format = "%m/%d/%y")

## enumerate study population vs subgroups-
#the variable stpop tells me whether the study population is the same as the analysed population (as opposed to subgroup)
outcomes$stpop<-as.numeric(outcomes$study_population==outcomes$population)

## clean study and subgroup population names--
outcomes$study_population[outcomes$study_population=="Women,Adults"]<-"Adults"
outcomes$population[outcomes$population=="Women,Adults"]<-"Women"

## create **out_pop**, which is an outcome population to analyse --
outcomes$stpop<-as.numeric(outcomes$study_population==outcomes$population)
out_pop<-outcomes$population
stpop<-outcomes$stpop
out_pop[stpop==0 & (out_pop=="Men"|out_pop=="Women")]<-paste(
  out_pop[stpop==0 & (out_pop=="Men"|out_pop=="Women")],
  outcomes$study_population[stpop==0 & (out_pop=="Men"|out_pop=="Women")],
  sep=" ")
outcomes$out_pop<-out_pop
rm(out_pop,stpop)

## create sex variable--
outcomes$sex<-"Mixed"
outcomes$sex[grepl("Women",outcomes$out_pop)]<-"Women"
outcomes$sex[grepl("Men",outcomes$out_pop)]<-"Men"


## CREATE A DATASET WITH LONGITUDINAL STUDIES AND SMDs and match with the meta-data.
#Create a dataset with only longitudinal studies
longidataset=outcomes[outcomes$is_longitudinal==1,]#keep only longitudinal data
longidataset$mainanalysis<-
  longidataset$population==longidataset$study_population # mainanalysis variable takes 1 for the rows indicating the study population so we will not perform subgroup analyses, but the main analyses.
longidataset$continuousout<-!is.na(longidataset$sd) #continuousout variable takes 1 for rows with continuous outcomes
trulylongi<-names(table(longidataset$record_id[longidataset$mainanalysis]))[table(longidataset$record_id[longidataset$mainanalysis])>=2]#keep studies with trully at least 2 timepoints
longidataset<-longidataset[longidataset$record_id %in% trulylongi,]

#  match the dataset with the characteristics --
metadata<-subset(metadata, select = -c(country) )
inlongi<-unique(longidataset$record_id)
longidmetadata<-metadata[metadata$record_id %in% inlongi,]
longi<-merge(longidataset,longidmetadata, by = "record_id")
longi<-select(longi,!ends_with(".y"))
rm(longidataset,longidmetadata,inlongi, trulylongi)


## CREATE SMD for CONTINUOUS outcome ----
longicont=longi[!is.na(longi$sd) & longi$mainanalysis, ]
rm(longi)
longicont<- filter(longicont,record_id!="66150") #this study has different instruments per time point
longicont=arrange(longicont,record_id,condition,scale)
#the data we need should have 1. all sd observed 2. no single arms 3. study_design only 1 and 2

###### Calculate SMD by study design and create the "data" dataset ready for analysis--
source("calculate_smd_vectorised.R")
grouped <- longicont %>% group_by(record_id, scale, population)
longicont <- grouped %>%
  filter(mainanalysis)  %>%
  filter(!is.na(sd)) %>%
  filter(n() > 1) %>%
  dplyr::select(everything())

length(table(longicont$record_id))

#calculate SMDs
longicont<-ungroup(longicont)
longicont<-arrange(longicont, record_id, condition, scale, timepoint)
SMDs<-calculate_smd(longicont) 
SMDs<-arrange(SMDs, record_id, condition, scale, timepoint)%>% 
  dplyr::select(record_id, s2,d,var,starts_with("cov"))
SMDs<-ungroup(SMDs)
#keep useful variables from the longicont
sublongicont<-longicont %>% dplyr::select(record_id , doi,author_1,year,population , 
                                          study_population , country , condition , scale , timepoint , 
                                          n_timepoints , is_longitudinal , is_prepandemic , sample_size ,
                                          is_binary , score ,  sd ,  y , sey , stringency , days_after_first ,
                                          days_after_pandemic,confirmed_cumulative , deaths_cumulative ,
                                          confirmed_avg , deaths_avg , cumulative_stringency ,country_population_2019,
                                          confirmed_per_100000 , deaths_per_100000 ,confirmed_avg_per_100000 , deaths_avg_per_100000,
                                          gdp_per_capita_2019 , gini_2019 , effect_direction , rob_info_bias.x , num_invited.x , 
                                          num_assessed.x , response_rate.x , rob_is_target_pop.x , rob_non_bias.x ,
                                          stpop , out_pop , sex , mainanalysis , continuousout ,  pop0_name , 
                                          tag_num_timepoints , study_design , method_recruitment , method_collection ,
                                          multi_country , pop0_central_age, pop0_sample_size , pop0_mean_med ,
                                          pop0_sd , pop0_min_age , pop0_max_age , pop0_percent_female ,
                                          pop_percent_phys_con , pop_percent_psych_con ,  pop_percent_covid19 ,
                                          pop_ethnicity , pop_ethnicity_other, ESI,CHIndex,School_closing,
                                          Workplace_closing,Stay_home_req,Restrictions_movement,Facial_cover)

#create the full data for analysis. Exclude columns with dupicate names in SMDs
if(sum(range(longicont$record_id-SMDs$record_id))==0){
  dataset <- cbind.data.frame(sublongicont, dplyr::select(SMDs,-"record_id",-"population",-"scale"))}
## sort the dataset                                   
dataset=arrange(dataset, record_id, condition, scale, timepoint)
dataset$exact_days_after_first<-dataset$days_after_first
dataset$days_after_first[dataset$days_after_first<0]=0

dataset<-filter(dataset,!(record_id==121277 & condition=="Social anxiety"))# one study has two different measures of anxiety (social anxiety and GAD)
dataset<-filter(dataset,!(record_id==121277 & condition=="Panic/Somatic symptoms"))

###Create covariates---- 
dataset$condition[dataset$condition %in% c("Social anxiety","Generalized Anxiety Disorder", "Panic/Somatic symptoms")]<-"Anxiety"
dataset$effect_direction[is.na(dataset$effect_direction)]<-1
dataset$d<-dataset$effect_direction*dataset$d ## changing the effect direction!
dataset$giniMinusmin<-dataset$gini_2019-min(dataset$gini_2019,na.rm=T)
dataset$gdpDiv10000minusmin<-dataset$gdp_per_capita_2019/10000-min(dataset$gdp_per_capita_2019/10000,na.rm=T)
dataset$study_design01<-dataset$study_design-1
dataset$scaled.score<-(dataset$score)/sqrt(dataset$s2)
dataset$logcumulative_stringency<-log(dataset$cumulative_stringency+1)
dataset$sqrtcumulative_stringency<-sqrt(dataset$cumulative_stringency)
dataset$logconfirmed_cumulative100k=log(dataset$confirmed_cumulative/dataset$country_population_2019*100000+1)
dataset$logdeaths_cumulative100k=log(dataset$deaths_cumulative/dataset$country_population_2019*100000+1)
dataset$sqrtconfirmed_cumulative100k=sqrt(dataset$confirmed_cumulative/dataset$country_population_2019*100000)
dataset$sqrtdeaths_cumulative100k=sqrt(dataset$deaths_cumulative/dataset$country_population_2019*100000)
dataset$days_after_pandemic[is.na(dataset$days_after_pandemic)]<-dataset$days_after_first[is.na(dataset$days_after_pandemic)]
#create RoB variables
repfun<-function(x){
  x[x=="High risk"]<-2
  x[x=="Low risk"]<-0
  x[x=="Unclear risk"]<-1
  as.numeric(x)}
numRoB<-as.data.frame(
  apply(cbind(dataset$rob_is_target_pop.x,dataset$rob_info_bias.x,dataset$rob_non_bias.x),2,repfun))
names(numRoB)<-c("RoBRepres","RoBInfoBias","RoBNonResp")
numRoB$record_id=dataset$record_id
numRoB$condition=dataset$condition
numRoB$sample_size=dataset$sample_size
numRoB=ungroup(numRoB %>% group_by(record_id,condition) 
               %>% group_modify(~ mutate(.x,RoBRepStudy=max(RoBRepres))) 
               %>% group_modify(~ mutate(.x,RoBInfoStudy=max(RoBInfoBias))) 
               %>% group_modify(~ mutate(.x,RoBNonRespStudy=max(RoBNonResp)))
               %>% group_modify(~ mutate(.x,meanStudySS=mean(sample_size))))

numRoBdich<-as.data.frame(t(apply(select(numRoB,RoBRepStudy,RoBInfoStudy,RoBNonRespStudy),1,function(x) as.numeric(x>0))))
names(numRoBdich)<-c("RoBRepresdich","RoBInfoBiasdich","RoBNonRespdich")
dataset<-cbind.data.frame(dataset,numRoB %>% select(RoBRepStudy,RoBInfoStudy,RoBNonRespStudy,meanStudySS),numRoBdich)
dataset$USA<-as.numeric(dataset$country=="United States")
dataset$CHINA<-as.numeric(dataset$country=="China")
#recode funny names
dataset$author_1[dataset$author_1=="Br√§scher"]<-"Braescher"
dataset$author_1[dataset$author_1=="Gim√©nez-Das√≠"]<-"Gimenez-Dasi"
dataset$author_1[dataset$author_1=="Sch√§fer"]<-"Schaefer"
dataset$author_1[dataset$author_1=="Sch√ºtzwohl"]<-"Schuetzwohl"
dataset$author_1[dataset$author_1=="S√∏nderskov"]<-"Sonderskov"

dataset$ESI[is.na(dataset$ESI) & dataset$is_prepandemic==1]<-0
dataset$CHIndex[is.na(dataset$CHIndex)&  dataset$is_prepandemic==1]<-0
dataset$School_closing[is.na(dataset$School_closing)&  dataset$is_prepandemic==1]<-0
dataset$Workplace_closing[is.na(dataset$Workplace_closing)&  dataset$is_prepandemic==1]<-0
dataset$Restrictions_movement[is.na(dataset$Restrictions_movement)&  dataset$is_prepandemic==1]<-0
dataset$Stay_home_req[is.na(dataset$Stay_home_req)&  dataset$is_prepandemic==1]<-0
dataset$Facial_cover[is.na(dataset$Facial_cover)&  dataset$is_prepandemic==1]<-0

dataset$authoryear<-paste(dataset$author_1,dataset$year)

rm(SMDs,sublongicont,grouped,numRoB,numRoBdich, longicont)   


write.csv(dataset,"MHCOVIDdataset.csv" )