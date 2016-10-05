###################################################################
# EFDM core
# October 6, 2015
# Authors: Seija Sirkiä, Maria Eronen, Tuula Packalen
# Original Author: Seija Sirkiä, January 25th 2013
#   
#   
#   Copyright 2015 European Union
#   
#   Licensed under the EUPL, Version 1.1 or – as soon they
#   will be approved by the European Commission - subsequent
#   versions of the EUPL (the "Licence");
#   You may not use this work except in compliance with the
#   Licence.
#   You may obtain a copy of the Licence at:
#   
#   http://joinup.ec.europa.eu/software/page/eupl/licence-eupl
#   
#   Unless required by applicable law or agreed to in
#   writing, software distributed under the Licence is
#   distributed on an "AS IS" basis,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
#   express or implied.
#   See the Licence for the specific language governing
#   permissions and limitations under the Licence.
#   
###################################################################
#--------------------------------
# functions needed:

library(abind) #for function abind

readtransprobs<-function(filename)
  #for reading a ready made transition matrix or several from
  #a text file (file extension .txt), or from an R save file 
  #(file extension .RData)
{
  if(tail(strsplit(filename,"[.]")[[1]],1)=="RData")
  {
    tmp<-load(filename) #load should create a single new object
    #in this local environment, it still needs to be returned
    return(get(tmp))
    #if there were several, if it is the wrong shape...
    #who knows what will happen!
  }
  #else
  rawdata<-readLines(filename)
  if(length(rawdata)==basedim) # if so, 
    #there's just one transition matrix used everywhere
  {
    tmp<-matrix(as.numeric(unlist(strsplit(rawdata,split=" "))),
                ncol=basedim,byrow=TRUE)
    transmat<-rep(tmp,times=prod(otherdim))
    dim(transmat)<-c(basedim,basedim,otherdim)
    return(transmat)
  }
  #otherwise, there are several matrices, for different factor levels
  
  #Which factors are involved? This should be mentioned on the first line
  tmp<-unlist(strsplit(rawdata[1],split=" "))
  tmp<-unlist(strsplit(tmp,split="="))
  fs<-tmp[seq(1,length(tmp),by=2)] #first, third etc are the factor names
  fn<-factdims[fs]
  #rest of the file is the matrices, separated by "header" lines
  frows<-seq(from=1,length=prod(fn),by=basedim+1) #these are the header lines
  tmp<-matrix(as.numeric(unlist(strsplit(rawdata[-frows],split=" "))),
              ncol=basedim,byrow=TRUE)
  #Actually, the header lines are thrown away! It doesn't matter what's on them,
  #but they probably help keeping the file organized as assumed
  dim(tmp)<-c(basedim,fn,basedim)
  tmp<-aperm(tmp,c(1,length(fn)+2,1:length(fn)+1))
  os<-setdiff(otherdimnames,fs)
  transmat<-rep(tmp,times=prod(factdims[os]))
  dim(transmat)<-c(basedim,basedim,fn,factdims[os])
  nr<-1:length(otherdimnames)
  names(nr)<-c(fs,os)
  transmat<-aperm(transmat,c(1,2,nr[otherdimnames]+2))
  return(transmat)
}

mmf<-function(M)
  #simply for doing the necessary "matrix times a vector" multiplications
  #somewhat efficiently
{
  m<-dim(M)[1]
  M[,1:m]%*%M[,m+1]
}

dividebyA<-function(oldstate)
{
  result<-vector("list",nact)
  #if dynamic updating of A, use a function to get a new A:
  #actproblist<<-do.call(updateA,actproblist)
  
  #while(TRUE) { #if iterating A
  for(i in 1:nact) {
    nr<-1:nfact
    names(nr)<-factnames    
    tmp<-oldstate*actproblist[[i]]
    #tmp is the subpopulation receiving the activity i
    result[[i]]<-tmp
  } 
  #if(!iterating) break #if no iteration, we're done here
  #else {
  #tmp<-do.call(iterateA,args=list(result,actproblist))
  #if(is.null(tmp)) break 
  #else if(usethesameAonnextstep) actproblist<<-tmp
  #else actproblist<-tmp #only change the local A
  #} #else, we iterated
  #} while
  names(result)<-actnames
  result
}


newstate<-function(oldstatediv) 
  #for producing the next state from the current one, using the obtained
  #activity and transition probabilities
{
  result<-oldstatediv
  for(i in 1:nact) {
    nr<-1:nfact
    names(nr)<-factnames
    acti<-activities[[i]]
    factnames.here<-c(acti,setdiff(factnames,acti))
    tmp<-aperm(oldstatediv[[actnames[i]]],nr[factnames.here])
    #tmp is the subpopulation receiving the activity i, with appropriate dims
    dim(tmp)<-c(prod(dim(tmp)[1:length(acti)]),1,prod(dim(tmp)[-(1:length(acti))]))
    transmati<-transmats[[i]]
    #transmati is the transition matrix corresponding to activity i
    dim(transmati)<-c(dim(transmati)[1:2],prod(dim(transmati)[-(1:2)]))
    tmp<-apply(abind(transmati,tmp,along=2),3,mmf)
    #tmp now has the new states of this subpopulation
    dim(tmp)<-factdims[factnames.here]
    names(nr)<-factnames.here
    tmp<-aperm(tmp,nr[factnames])
    dimnames(tmp)<-dimnames(tmp)
    result[[actnames[i]]]<-tmp
  }

  result
}




#describing the factors:
factors.tmp<-readLines(STATESPACE_FILENAME)
nfact<-length(factors.tmp)
factors.tmp<-strsplit(factors.tmp," ")
factnames<-character(nfact)
for(i in 1:nfact) {
  factnames[i]<-factors.tmp[[i]][1]
  factors.tmp[[i]]<-factors.tmp[[i]][-1]
}
names(factors.tmp)<-factnames
factlvls<-factors.tmp
factdims<-sapply(factlvls,length)


#describing the activities
activities.tmp<-readLines(ACTIVITIES_FILENAME)
nact<-length(activities.tmp)
activities<-strsplit(activities.tmp," ")
actnames<-character(nact)
actmethod<-character(nact)
actfiles<-character(nact)
rm(activities.tmp)

for(i in 1:nact) {
  actnames[i]<-activities[[i]][1]
  actmethod[i]<-activities[[i]][2]
  actfiles[i]<-activities[[i]][3]
  activities[[i]]<-activities[[i]][-(1:3)]
}
names(activities)<-actnames

#setting up the transition matrices
transmats<-list()
for(i in 1:length(activities))
{
  #the levels of these factors may change under this activity:
  basedimnames<-activities[[i]] 
  #these are fixed under this activity:
  otherdimnames<-setdiff(factnames,basedimnames)
  basedim<-prod(factdims[basedimnames])
  otherdim<-factdims[otherdimnames]
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #11.2.2014: This part here. Apparently the transition probability matrices
  #should come to the core in a file, end of story. The estimation routine would
  #be used beforehand (as an add-on, as we've been calling it) and written in a file
  #so that core wouldn't need to care where it came from.
  #17.4.2014: Nope, not gonna happen
  if(actmethod[i]=="read") 
    transmat.forthisround <- readtransprobs(actfiles[i]) 
  if(actmethod[i]=="estimate") {
    if(!exists("estimatetransprobs")) source("efdmestim.r")
    transmat.forthisround <- estimatetransprobs(actfiles[i])
  }
    
  if(actmethod[i]=="custom") 
    transmat.forthisround<- eval(call(actfiles[i]))
  #the functions that are called here have to be written and loaded by
  #the user!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  
  
  #the rest is making sure dimensions have correct and sensible names
  lvlnames<-factlvls[[basedimnames[1]]]
  if(length(basedimnames>1))
    for(j in 2:length(basedimnames))
    {
      lvlnames.more<-factlvls[[basedimnames[j]]]
      lvlnames<-paste(rep(lvlnames,times=length(lvlnames.more)),
                      rep(lvlnames.more,each=length(lvlnames)),sep=".")
    }
  dimnames(transmat.forthisround)<-c(list(lvlnames,lvlnames),factlvls[otherdimnames])
  transmats<-c(transmats,list(transmat.forthisround))
  names(transmats)[i]<-actnames[i]
}

rm(basedim,basedimnames,otherdim,otherdimnames,transmat.forthisround,lvlnames,
   lvlnames.more)

#setting up the initial state and the probabilities of activities
statespace<-read.table(INITSTATE_FILENAME,header=TRUE)
actprobtable<-read.table(ACTPROBS_FILENAME,header=TRUE)
#this next merge is a bit dumb but the fixing of factor levels and ordering
#needs to be done for both, so I'm doing it at once
statespace<-merge(statespace,actprobtable)
for(fname in factnames)
  statespace[[fname]]<-factor(statespace[[fname]],
                              levels=factlvls[[fname]],ordered=TRUE)
statespace<-statespace[do.call(order,statespace[rev(factnames)]),]

initstate<-statespace[[setdiff(names(statespace),union(factnames,actnames))]]
dim(initstate)<-factdims
dimnames(initstate)<-factlvls


actproblist<-list()
for(i in 1:nact) {
  actproblist[[i]]<-array(statespace[[actnames[i]]],dim=factdims)
  dimnames(actproblist[[i]])<-factlvls
}

#-----------------------------
#finally, the actual simulation

#resultstates<-statespace[factnames]

state<-dividebyA(initstate)
resultstates<-list(step0=state)


#for(i in 1:nact)
#resultstates[[actnames[i]]]<-data.frame(step0=as.vector(tmp[[i]]))

nrofsteps<-as.numeric(NROFSTEPS) #how many steps should be taken?
#a silly compatibility issue there with the names, sorry about that

for (i in 1:nrofsteps) {
  nstate<-newstate(state)  
  state<-rowSums(do.call(cbind,lapply(nstate,as.vector)))
  dim(state)<-factdims
  state<-dividebyA(state)
  resultstates[[paste0("step",i)]]<-state
}
#resultstates is now a list of lists, outer list going over the simulation steps, 
#and inner lists going over the activities. The inmost elements are multiway arrays
#in the form of the statespace (that is, factdims)

#Reformulate resultstates to be data frame whose columns are arrays, with nact columns.
#First, vectorize and cbind the multiway arrays:
resultstates<-lapply(resultstates,function(elem) {do.call(cbind,lapply(elem,as.vector))})
#Second, make the data frame:
resultstates<-data.frame(resultstates)

rawoutput<-statespace[factnames]
rawoutput<-data.frame(rawoutput,resultstates)
#rawoutput has 'nact' columns per step, with names suchs as step0.actname
write.table(rawoutput,file="rawoutput.txt")

tmp<-as.matrix(resultstates)
dim(tmp)<-c(dim(tmp)[1],nact,nrofsteps+1)
tmp<-apply(tmp,c(1,3),FUN=sum)
resultstates<-data.frame(tmp)
names(resultstates)<-paste0("step",0:nrofsteps)
resultstates<-data.frame(statespace[factnames],resultstates)
#resultstates is like rawoutput, but with the different columns with name
#such as step0.something summed together and named just step0
write.table(resultstates,file="resultstates.txt")