# revised version
	# the names of the clinical.tsv and follow_up.tsv were not correct (see revised details below)
	# I've also removed adding the 'days_to_follow_up' variable as we do not need it, we just need 'days_to_death' and 'days_to_last_follow_up' from the clinical.tsv file
	# looks like we only need the 'clinical.tsv' file to run the basic survival analysis

#### GDC - manual extraction of data for survival analysis in R    
# - https://portal.gdc.cancer.gov/
# - bone marrow cohort, NOTCH 1 mutations
#  - Here I'm going to manually extract data for individuals who have cancer related to tissue:bone marrow, and any individual within that cohort that has a NOTCH1 mutation
#     - For this example, there are a good number of individuals included (n=7,625)
#     - NOTCH1 was chosen simply as it was one of the most common genes with a mutation in this cohort
#
# - IMPORTANT NOTE for manual extraction
#     - You will need to download twice, firstly (STEP 1. below) for your main cohort and secondly (STEP 2. below) to define those individuals who had a mutation within your cohort
#     
# STEP 1. define your main cohort/group first
#     - Click on 'Cohort Builder' tab, this is where you select individuals for your cohort/group
#     - For my cohort, I selected 'bone marrow' under 'Tissue or Organ of Origin' - there are 7,625 cases (this will appear in the upper right window)
#     - Click on the CASES tab - in my cohort there are '7,625 CASES' showing - there should be a drop down menu with 'Files', 'Custom Filters', 'Biospecimen' and 'Clinical'
# 		- Click on the arrow next to 'Clinical' and select TSV to download .tsv files for your main cohort
#		- This will download a folder with 5 files (follow_up.tsv, pathology_detail.tsv, exposure.tsv, family_history.tsv, clinical.tsv) - for this simple survival analysis, we'll only be using the clinical.tsv file
# STEP 2. define your mutation next 
#     - with your selected cohort still showing, you now need to select a mutation
#     - click on the 'Analysis Center' tab, then click on the 'Mutation Frequency' play button
#     - Here I'm going to select NOTCH1 from the gene list - out of my selected 7,625 cases, there are 181 with a NOTCH1 mutation
#     		- Click on the '+' next to NOTCH1 where it says 181
#		- It will ask to save those 181 cases as a new cohort - provide a name and save as a new cohort
#		- From your saved cohorts list, now select the new cohort that you just saved
#		- You should see the number of cases change in the upper right, mine now shows '181 CASES'
#     		- Now I click on the 181 CASES tab, which will bring up the drop down menu with 'Files', 'Custom Filters', 'Biospecimen' and 'Clinical'
# 		- Same process as before, click on the arrow next to 'Clinical' and select TSV to download .tsv files that we'll use to define individuals with our chosen mutation
#		- This will again download another folder with 5 files (follow_up.tsv, pathology_detail.tsv, exposure.tsv, family_history.tsv, clinical.tsv) - again we're just using the clinical.tsv file from this second download to define those who had the NOTCH1 mutation in our cohort


# Now we can process the data in R

# 1. read in the first clinical.tsv file, downloaded from STEP 1. above
# C = read.delim("~/Downloads/clinical.cohort.2025-03-25/clinical.tsv", sep = "\t", header = TRUE, fill = TRUE)
# I've saved my main clinical.tsv file online - you can read this directly into your R session if you need to see what your data should look like
C = read.delim("https://bioinformatics.erc.monash.edu/~sbya0003/tmp/TRM5006_6006_2025/tutorial1/main_clinical/clinical.tsv", sep = "\t", header = TRUE, fill = TRUE)


# 2. extract days to death data from the clinical data file
# Cdtd=C[,c("case_id","days_to_death")] # previously I had thes two columns named like this ("case_id","days_to_death"), but in the download they may be named like this ("cases.case_id","demographic.days_to_death")
Cdtd=C[,c("cases.case_id","demographic.days_to_death")]
names(Cdtd)=c("case_id","days_to_death")
Cdtd=Cdtd[c(which(Cdtd$days_to_death != "'--")),]
head(Cdtd) # have a look at the head of the data
length(which(duplicated(Cdtd$case_id))) # how many duplicated Ids are there - 3866
Cdtd=Cdtd[c(which(!duplicated(Cdtd$case_id))),] # remove duplicated Ids
dim(Cdtd) # 1496 rows
hist(as.numeric(Cdtd$days_to_death),breaks=100) # lets look at a histogram of days until death


# 3. now extract days to last follow-up from the clinical data file
# this will be used to code follow-up time for individuals who were still alive at the end of the observed follow-up
#Cfu=C[,c("case_id","days_to_last_follow_up")] # just extract the variables that we need
Cfu=C[,c("cases.case_id","diagnoses.days_to_last_follow_up")] # just extract the variables that we need
names(Cfu)=c("case_id","days_to_last_follow_up")
Cfu=Cfu[c(which(Cfu$days_to_last_follow_up != "'--")),] # need to remove rows with this 'filler'
dim(Cfu) # 7866 rows
head(Cfu)
Cfu=Cfu[c(which(!duplicated(Cfu$days_to_last_follow_up))),] # this file includes alot of duplicates to help 'pad' for other variables that may have more than one level recorded for one patient, i.e. they will have multiple rows which cause duplicates to appear for this variable that need to be removed. 

# 4. now we can merge the death and follow-up variables together, using the case Id as the merge variable
M=merge(Cdtd,Cfu,by="case_id",all=TRUE)

# because we have 2 columns of data, will write a loop to pull out 'time' (time to death or last followup) and 'status' (dead/alive at that time)
# this loop will first identify those that have died...
# ...then if the individual did not die, look in the days_to_last_follow_up column for a day value

# 5. create the time (time to event or last follow-up) and status (binary, yes=event, no=no event) variables
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
length(which(is.na(M$time))) # checking if there are any NA values in the time column that we just coded - should be 0


# now we need to add the data from STEP 2 above, which we'll use to define those individuals in our cohort who have a NOTCH1 mutation
# 6. extracting and adding the mutation data
#NOTCH1 = read.delim("~/Downloads/clinical.cohort.2025-03-25(1)/clinical.tsv", sep = "\t", header = TRUE, fill = TRUE)
# I've saved my mutation clinical.tsv file online - you can read this directly into your R session if you need to see what your data should look like
NOTCH1 = read.delim("https://bioinformatics.erc.monash.edu/~sbya0003/tmp/TRM5006_6006_2025/tutorial1/NOTCH1_clinical/clinical.tsv", sep = "\t", header = TRUE, fill = TRUE)

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


# 7. now we are ready to run a basic survival analysis
library(ggsurvfit)
#dev.new(height=4, width=4)
survfit2(Surv(time, status) ~ mutant, data = M) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
)
