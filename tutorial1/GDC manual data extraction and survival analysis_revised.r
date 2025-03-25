


#### GDC - manual extraction of data for survival analysis in R    
# - https://portal.gdc.cancer.gov/
# - bone marrow - NOTCH 1
#  - Here I'm going to manually extract data for individuals who have cancer related to tissue:bone marrow, and any individual within that cohort that has a NOTCH1 mutation
#     - For this example, there are a good number of individuals included (n=7,625)
#     - NOTCH1 was chosen simply as it was one of the most common genes with a mutation in this cohort
# 
#   - GDC manual extraction
#     - Click on 'Cohort Builder' tab, this is where you select individuals for your dataset
#     - From Tissue or Organ of Origin, selected bone marrow
#     - Now click on 'Analysis Center' > 'Mutation Frequency'
#     - Selected NOTCH1 - So there are 7,625 cases, 181 with a NOTCH1 mutation
#     - Save Cohort > Save As...
#     - click on the 7,625 CASES button > Clinical tab > TSV
#       - this will download the data as a Zip file - place it somewhere on your laptop where its easy to point to in R. Below I've placed it in a folder called GDCdata on my desktop

# 1. read in the follow_up.tsv file
F = read.delim("~/Downloads/clinical.cohort.2025-03-25/follow_up.tsv", sep = "\t", header = TRUE, fill = TRUE)

# 2. read in the clinical.tsv file
C = read.delim("~/Downloads/clinical.cohort.2025-03-25/clinical.tsv", sep = "\t", header = TRUE, fill = TRUE)


# 3. extract days to death data from the clinical data file
# Cdtd=C[,c("case_id","days_to_death")] # previously I had thes two columns named like this ("case_id","days_to_death"), but in the download they may be named like this ("cases.case_id","demographic.days_to_death")
Cdtd=C[,c("cases.case_id","demographic.days_to_death")]
names(Cdtd)=c("case_id","days_to_death")
Cdtd=Cdtd[c(which(Cdtd$days_to_death != "'--")),]
head(Cdtd) # have a look at the head of the data
length(which(duplicated(Cdtd$case_id))) # how many duplicated Ids are there
Cdtd=Cdtd[c(which(!duplicated(Cdtd$case_id))),] # remove duplicated Ids
dim(Cdtd)
hist(as.numeric(Cdtd$days_to_death),breaks=100) # lets look at a histogram of days until death


# 4. now extract days to last follow-up from the clinical data file
# this will be used to code follow-up time for individuals who were still alive at the end of the observed follow-up
#Cfu=C[,c("case_id","days_to_last_follow_up")] # just extract the variables that we need
Cfu=C[,c("cases.case_id","diagnoses.days_to_last_follow_up")] # just extract the variables that we need
names(Cfu)=c("case_id","days_to_last_follow_up")
Cfu=Cfu[c(which(Cfu$days_to_last_follow_up != "'--")),] # need to remove rows with this 'filler'
dim(Cfu)
head(Cfu)
Cfu=Cfu[c(which(!duplicated(Cfu$days_to_last_follow_up))),] # this file includes alot of duplicates to help 'pad' for other variables that may have more than one level recorded for one patient, i.e. they will have multiple rows which cause duplicates to appear for this variable that need to be removed. 

# 5. now we can merge the death and follow-up variables together, using the case Id as the merge variable
M=merge(Cdtd,Cfu,by="case_id",all=TRUE)

# because we have 2 columns of data, will write a loop to pull out 'time' (time to death or last followup) and 'status' (dead/alive at that time)
# this loop will first identify those that have died...
# ...then if the individual did not die, look in the days_to_last_follow_up column for a day value

M$time=NA
M$status=0
for(i in 1:nrow(M)) {
	if(!is.na(M$days_to_death[i])) {
		M$time[i]=M$days_to_death[i]
		M$status[i]=1
	} else {
		if(!is.na(M$days_to_last_follow_up[i])) {M$time[i]=M$days_to_last_follow_up[i]} 
	}
}
head(M)
length(which(is.na(M$time))) # checking if there are any NA values in the time column that we just coded


#  - now we need to define who has a NOTCH1 mutation - back to GDC to manually extract this
#    - there are n=181 with a NOTCH1 mutation
#    - Click on '+' to Save a new cohort of these cases
#    - click on 181 Cases tab
#    - click on Clinical > TSV to download
#    - now we will read this NOTCH1 subset in, process and merge it into our dataset


NOTCH1 = read.delim("~/Downloads/clinical.cohort.2025-03-25(1)/clinical.tsv", sep = "\t", header = TRUE, fill = TRUE)
NOTCH1$mutant=1 # indicator
NOTCH1=NOTCH1[,c("cases.case_id","mutant")] # too many columns, lets just pull out the Id and NOTCH1 indicator
names(NOTCH1)[1]="case_id"
NOTCH1=NOTCH1[c(which(!duplicated(NOTCH1$case_id))),]
dim(NOTCH1) # n=181, this should be the same number of individuals with the mutation shown on the GDC portal

# now merge this in with the follow-up data
dim(M); M=merge(M,NOTCH1,by="case_id",all=TRUE); dim(M)
# there were slightly more (n=4) after the merge - this means we have 4 individuals with missing followup data
M=M[c(which(!is.na(M$status))),] # remove NA values that have popped up in the merge
M$mutant[c(which(is.na(M$mutant)))]=0
table(M$mutant) # note there will likely be less individuals with the mutation after the merge - this is simply because not all individuals with the mutation will have a days to last follow up value
M$time=as.numeric(M$time)
M[1:30,]


#  - now we are ready to run the survival analysis
  
library(ggsurvfit)
#dev.new(height=4, width=4)
survfit2(Surv(time, status) ~ mutant, data = M) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
)



