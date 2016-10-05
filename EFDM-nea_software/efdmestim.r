###################################################################
# EFDMestim
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
estimatetransprobs<-function(filename)
  #for estimating a transition matrix from data using the Bayes-like 
  #algorithm explained elsewhere
{
  tmp<-strsplit(readLines(filename)," ")
  #the data:
  pairdata<-read.table(tmp[[1]][2],header=TRUE)
  #the prior, either a file name or an instruction to use uninformative prior 
  #(unit matrix)
  prior<-tmp[[2]][2]
  #the order of factors, "from most influential to least influential", 
  #together with corresponding weights
  factorder<-unlist(strsplit(unlist(strsplit(tmp[[3]],split=" ")),split="="))
  
  #the following ensures that factors are represented in the objects correctly
  #and consistently
  for(j in 1:length(otherdimnames))
    pairdata[[otherdimnames[j]]]<-factor(pairdata[[otherdimnames[j]]],
                                         levels=factlvls[[otherdimnames[j]]],ordered=TRUE)
  
  name0<-paste0(basedimnames,0)
  name1<-paste0(basedimnames,1)
  for(j in 1:length(basedimnames))
  {
    pairdata[[name0[j]]]<-factor(pairdata[[name0[j]]],
                                 levels=factlvls[[basedimnames[j]]],ordered=TRUE)
    pairdata[[name1[j]]]<-factor(pairdata[[name1[j]]],
                                 levels=factlvls[[basedimnames[j]]],ordered=TRUE)
  }
  #prior is either read from another file or set up as a unit matrix:
  if(prior=="uninformative") s<-diag(basedim) else
    s<-matrix(as.numeric(unlist(strsplit(readLines(prior),split=" "))),ncol=basedim,byrow=TRUE)
  
  dim(factorder)<-c(2,length(factorder)/2)
  pairdata<-pairdata[,c(name1,name0,factorder[1,])]
  
  #start of estimations algorithm
  freq<-table(pairdata)
  dim(freq)<-c(basedim,basedim,factdims[factorder[1,]])
  priorweights<-c(1,as.numeric(factorder[2,]))
  #To try and explain what is going on in the following for loop:
  
  #1: On the 1st round, tmp is 2-dimensional, the frequency table of
  #the changing factor levels, pooled over every non-changing factor.
  #The 2-dimensional table of prior "frequencies" in s is summed to 
  #to tmp, and saved again to s.
  
  #2: On the 2nd round, tmp is 3-dimensional, standing for the 2-dimensional
  #frequency tables of the changing factor levels, one for each level of the
  #first non-changing factor.
  #The 2-dimensional table of prior "frequencies" in s obtained on the first
  #round is summed to each 2-dimensional factor1 specific subtable 
  #of tmp, and the 3-dimensional result is saved again to s.
  
  #3: On the 3rd round, tmp is 4-dimensional, standing for the 2-dimensional
  #frequency tables of the changing factor levels, one for each combination 
  #of the levels of the first 2 non-changing factors.
  #The 2-dimensional tables of prior "frequencies" in s obtained on the second
  #round for each level of the first non-changing factor, are summed to each 
  #2-dimensional factor1 and factor2 specific subtable of tmp, and the 
  #4-dimensional result is saved again to s.
  
  #And so on.
  for (i in 1:length(priorweights))
  {
    tmp<-apply(freq,1:(1+i),sum)
    s<-sweep(tmp,1:max(2,i),s*priorweights[i],'+')
  }
  #in the end, s contains all the numerator values for the singular 
  #probability estimates, and the denominators are the column sums of
  #all those 2-dimensional tables, obtained here:
  N<-apply(s,2:length(dim(s)),sum)
  #then, numbers in s are divided by the corresponding number in N:
  transmat<-sweep(s,2:length(dim(s)),N,'/')
  #end of estimation algorithm
  
  nr<-1:length(otherdimnames)
  names(nr)<-factorder[1,]
  transmat<-aperm(transmat,c(1,2,nr[otherdimnames]+2))
  return(transmat) 
}