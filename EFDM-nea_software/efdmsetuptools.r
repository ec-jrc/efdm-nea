###################################################################
# EFDMsetuptools
# October 6, 2015
# Authors:  Maria Eronen, Seija Sirkiä, Tuula Packalen
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
setupstatespace<-function(filename)
  #This is a function for initializing the state space.
  #(Does the same as core has to do when it starts, 
  #should ultimately be written just once somewhere!
{
  factors.tmp<-readLines(filename)
  nfact<<-length(factors.tmp)
  factors.tmp<-strsplit(factors.tmp," ")
  factnames<<-character(nfact)
  for(i in 1:nfact) {
    factnames[i]<<-factors.tmp[[i]][1]
    factors.tmp[[i]]<-factors.tmp[[i]][-1]
  }
  names(factors.tmp)<-factnames
  factlvls<<-factors.tmp
  factdims<<-sapply(factlvls,length)
}


pre.estimate<-function(inputfilename,factorsfilename,changing,resultfilename)
  #This is a function for performing an estimation of a transition matrix
  #"beforehand", and not during the running of the whole model/simulation
  #Arguments:
  #inputfilename: 
  #-the name of the estimation control file (that would be given 
  #in the activities control file)
  #factorsfilename:
  #-the name of the state space defining control file (factors.txt by default)
  #changing: 
  #-names of the factors that "have" the transitions (that would be given in
  #the activities control file)
  #resultfilename:
  #-name of the file where the result is to be saved (do use the extension .RData
  #if you plan to read this in during model/simulation run!)


  #the call of this function would then look something like this:
  #pre.estimate("estiminput.txt","factors.txt","vol age","nomgmtP.RData")
  #(with the needed data and possible prior available of course)
{
  setupstatespace(factorsfilename)
  basedimnames<<-strsplit("stem vol",split=" ")[[1]]
  otherdimnames<<-setdiff(factnames,basedimnames)
  basedim<<-prod(factdims[basedimnames])
  otherdim<<-factdims[otherdimnames]
  if(!exists("estimatetransprobs")) source("efdmestim.r")
  #create a name for the object based on the time and date
  #(I want to be fairly certain it is unique...)
  oname<-paste0(strsplit(date()," ")[[1]],collapse="")
  oname<-paste0(strsplit(oname,":")[[1]],collapse="")
  oname<-paste0("estimate",oname,collapse="")
  #note, this should not be relevant to the user, the file name
  #in which this object will be saved is a different thing
  assign(x=oname,value=estimatetransprobs(inputfilename))
  save(list=oname,file=resultfilename)
}

neutralprior<-function(nvol,nage,resultfilename)
  #This is a function for creating a simple prior for the
  #estimation procedure, for a "vol-age" model, implementing
  #the idea "grow one age class older, stay in the same vol
  #class, unless already in the highest age class, in which
  #case stay where you are"
  #Arguments:
  #nvol and nage
  #-number of vol and age classes
  #resultfilename
  #-name of the (text) file the output will be written in
{
  tmp<-diag(nage)
  tmp<-tmp[,c(2:nage,nage)]
  diag(nvol)%x%tmp
  write.table(tmp%x%diag(nvol),file=resultfilename,row.names=FALSE,col.names=FALSE)  
}

makeanewP<-function(inputfilename,shiftarray,shiftdims,resultfilename)
  #This is a function that takes a ready-made transition matrix
  #and creates a new one by adding shifts (transitions with probability 1)
  #from a state to another before the given transitions. These shifts
  #can only happen within the subspace of the states given by the 
  #dynamic factors of the original transition matrix
  #Arguments:
  #inputfilename
  #-name of the save file where the R object for the original matrix is
  #shiftarray
  #-an array giving the new state for each old state. The rows must
  #match the subspace, that is, must refer to the same states (level
#combinations) as the first and second dimension of the given
#original transition matrix. The columns give the new new state
#by specifying the rank of the new level of each
#dynamic factor, and must be in the same order as the factors
#appear in the subspace ordering. In other words, if the dynamic
#factors for the original matrix were named as "f1 f2 f3" then
#the levels of f1 run the fastest on the rows, and the columns
#are in the same order.
#shiftdims
#-a vector giving the number of levels of the dynamic factors,
#so that nrow(shiftarray)==prod(shiftdims)
#-resultfilename
#-name of the file where the result is to be saved (do use the extension .RData
#if you plan to read this in during model/simulation run!)
{
  oname<-load(inputfilename)
  origP<-get(oname)
  Pdims<-dim(origP)
  newP<-origP
  dim(newP)<-c(Pdims[1:2],prod(Pdims[-(1:2)]))
  dimcoef<-cumprod(c(1,shiftdims[-length(shiftdims)]))
  ind <- apply(t(shiftarray-1)*dimcoef,2,sum)+1
  newP<-newP[,ind,]
  dim(newP)<-Pdims
  #create a name for the object based on the time and date
  #(I want to be fairly certain it is unique...)
  oname<-paste0(strsplit(date()," ")[[1]],collapse="")
  oname<-paste0(strsplit(oname,":")[[1]],collapse="")
  oname<-paste0("thin",oname,collapse="")
  #note, this should not be relevant to the user, the file name
  #in which this object will be saved is a different thing
  assign(x=oname,value=newP)
  save(list=oname,file=resultfilename)
}

makeathinP.ea<-function(inputfilename,voldrops,resultfilename)
  #Special purpose function for making a transition matrix out of
  #a pre-made "no management" transition matrix, in case of a
  #"vol-age" model. The idea is that thinning means "drop so and
  #so many volume classes, then grow as usual", leading to moving 
  #columns right. 
  
  #In the following nvol refers to the number of volume classes
  #and nage to the number of age classes. These can be deduced
  #from the inputs: the length of voldrops = nvol, nage can 
  #then be calculated from the dimensions of the input matrix.
  
#Assumption 1: the pre-made matrix in the input file is a saved
#R object of type array of dimensions
#(nvol*nage) x (nvol*nage) x something, 
#where something is either nothing (single forest type), or 
#a single number (that many forest types defined by a single factor),
#or several numbers (that many forest type defining factors). This
#will be the case if it was created by "estimatetransprobs" function.

#Assumption 2: The "first" 2-dimensional array is made of blocks
#of size nvol x nvol arranged in nage x nage array. That is,
#the first nvol columns refer to the first age class (and not the
#other way around, first nage columns referring to the first vol
#class.) This will be the case if the matrix in the input file
#was created by the aforementioned function and specifying the 
#changing factors to be "vol age", in that order (the order is
#what matters, not the names used)

#Arguments:
#inputfilename
#-name of the save file where the R object for "no management" matrix is
#voldrops
#-a vector of integers, how many levels to drop at each level (length 
#equal to number of volume levels)
#-resultfilename
#-name of the file where the result is to be saved (do use the extension .RData
#if you plan to read this in during model/simulation run!)

{
  oname<-load(inputfilename)
  tmp<-get(oname)
  nvol<-length(voldrops)
  dims<-dim(tmp)
  nage<-dims[1]/nvol
  ind<-1:nvol-voldrops
  shiftarr<-cbind(rep(ind,times=nage),rep(1:nage,each=nvol))
  makeanewP(inputfilename,shiftarr,c(nvol,nage),resultfilename)
}


