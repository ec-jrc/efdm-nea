###################################################################
# EFDMutils
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
efdmui<-function() {
  VARS<-character(0)
  
  cat("Hello, I am EFDM.\n") 
  
  ###############################
  VARS<-c(VARS,"STATESPACE_FILENAME")
  cat("I need you to describe the state space for me.\n")
  cat("Choose an existing description file (click cancel if not available)...\n")
  
  tmp<-try(file.choose(),silent=FALSE)
  if(class(tmp)=="try-error") {
    cat("Looks like you cancelled.\n")
    cat("Type in the names of the factors then, separated by whitespaces:\n")
    tmp<-readline()
    tmp<-strsplit(tmp,split=" ")
    factnames<-tmp[[1]]
    nfact<-length(factnames)
    factlvls<-vector("list",nfact)
    factdims<-rep(NA,nfact)
    names(factlvls)<-factnames
    cat(paste("I got", nfact,"factors:", paste0(factnames,collapse=" ","\n")))
    
    for(i in 1:nfact){
      cat(paste0("Type in the levels in the correct order for factor ",factnames[i],":\n"))
      tmp<-readline()
      factlvls[[i]]<-strsplit(tmp,split=" ")[[1]]
      factdims[i]<-length(factlvls[[i]])
      cat(paste("I got",factdims[i],"levels:", paste0(factlvls[[i]],collapse=" "),"\n"))
    }
    
    assign(x=VARS[length(VARS)],value="factors.txt")
    
    cat(paste(factnames,unlist(lapply(factlvls,FUN=function(item) paste(item,collapse=" ")))),
        file=get(VARS[length(VARS)]),sep="\n")
    cat("I have written a file",get(VARS[length(VARS)]), "of these for future use.\n")
  }
  else #tmp should be a filename
  {
    assign(x=VARS[length(VARS)],value=tmp)
  }
  #################################
  VARS<-c(VARS,"INITSTATE_FILENAME")
  cat("Choose the file for initial state...\n")
  assign(x=VARS[length(VARS)],value=file.choose())
    
  #################################
  VARS<-c(VARS,"ACTIVITIES_FILENAME")
  cat("I need you to describe the activities for me.\n")
  cat("Choose an existing description file (click cancel if not available)...\n")
  
  tmp<-try(file.choose(),silent=TRUE)
  
  if(class(tmp)=="try-error") {
    cat("Looks like you cancelled.\n")
    cat("Type in the names of the activities then, separated by whitespaces:\n")
    tmp<-readline()
    tmp<-strsplit(tmp,split=" ")
    actnames<-tmp[[1]]
    nact<-length(actnames)
    cat(paste("I got", nact,"activities:", paste0(actnames,collapse=" "),"\n"))
    
    activities<-vector("list",nact)
    actmethod<-actfiles<-character(nact)
    for(i in 1:nact){
      cat(paste0("Type in the factors that have transitions under ",actnames[i],":\n"))
      tmp<-readline()
      activities[[i]]<-strsplit(tmp,split=" ")[[1]]
      cat(paste("I got",length(activities[[i]]),"factors:", paste0(activities[[i]],collapse=" "),"\n"))
      cat("How do I get the corresponding P-matrix?\n")
      repeat {
        cat("Type 1 for reading it from a file, 2 for estimating it, 3 for a custom function:\n")
        tmp<-readline()
        tmp<-as.numeric(tmp)
        if(tmp%in%1:3) break
      }
      switch(tmp,{actmethod[i]<-"read"
                  cat("Choose the file to read...\n")
                  actfiles[i]<-file.choose()},
             {actmethod[i]<-"estimate"
              cat("Choose the input file for the estimation addon...\n")
              actfiles[i]<-file.choose()},
             {actmethod[i]<-"custom"
              cat("Type in the name of the custom function to call:\n")
              actfiles[i]<-readline()})
    }
    
    assign(x=VARS[length(VARS)],value="activities.txt")
    
    cat(paste(actnames,actmethod,actfiles, unlist(lapply(activities,FUN=function(item) paste(item,collapse=" ")))),
        file=get(VARS[length(VARS)]),sep="\n")
    cat("I have written a file",get(VARS[length(VARS)]), "of these for future use.\n")
  }
  else #tmp should be a filename
  {
    assign(x=VARS[length(VARS)],value=tmp)
  }
  #################################
  VARS<-c(VARS,"ACTPROBS_FILENAME")
  cat("Choose the file for activity probabilities...\n")
  assign(x=VARS[length(VARS)],value=file.choose())
  
  ########################################
  VARS<-c(VARS,"NROFSTEPS")
  cat("How many steps to run?\n")
  assign(x=VARS[length(VARS)],value=readline())
  
  
  ######################################
  VARS<-c(VARS,"OUTPUTREQUEST_FILENAME")
  cat("Do you want any output in addition to the raw results? (0=no 1=yes)\n")
  tmp<-readline()
  #actually anything other than 0 = yes
  if(tmp=="0") {
    assign(x=VARS[length(VARS)],value=tmp)
  } else {
    cat("I need you to describe the outputs you want.\n")
    cat("Choose an existing output request file (click cancel if not available)...\n")
    
    tmp<-try(file.choose(),silent=FALSE)
    if(class(tmp)=="try-error") {
      cat("Looks like you cancelled.\n")
      outputreq<-outputcoef<-character(0)
      repeat {
        cat("Type in the (first/next) output request then (hit enter for none):\n")
        tmp<-readline()
        if(!(tmp=="")) {
          outputreq<-c(outputreq,tmp)
          cat("Choose a file with the coefficients for this output variable...\n")
          outputcoef<-c(outputcoef,file.choose())
        }
        else break
      }
      cat("I got",length(outputreq),"output requests.")
      if(length(outputreq!=0)) {
        assign(x=VARS[length(VARS)],value="outputrequests.txt")      
        cat(c(rbind(outputreq,outputcoef)),file=get(VARS[length(VARS)]),sep="\n")
        cat("I have written a file",get(VARS[length(VARS)]), "of these for future use.\n")
      } else #they didn't want any after all
        assign(x=VARS[length(VARS)],value="0")
    }
    else #tmp should be a filename
    {
      assign(x=VARS[length(VARS)],value=tmp)
    }
    
  }
  
  ####################################
  
  inputfile<-"efdminput.txt"
  
  cat(paste(VARS, unlist(mget(VARS))),file=inputfile,sep="\n")
  cat("Thank you, I have written a file", inputfile, "of these inputs.\n")
  cat("Type", paste0("runefdm(\"",inputfile,"\")"), "when you are ready to run this scenario.\n")
  invisible(0) #just to have a sensible return value
} #end of efdmui function

#------------------------------------
runefdm<-function(inputfile)
{
 tmp<-readLines(inputfile)
 tmp<-matrix(unlist(strsplit(tmp,split=" ")),ncol=2,byrow=TRUE)
 VARS<-tmp[,1]
 apply(tmp,1,FUN=function(x) assign(x[1],x[2],pos=.GlobalEnv))
 source("efdmcore.r")
 if(OUTPUTREQUEST_FILENAME!=0) {
   source("efdmoutput.r")
   outputcalc(OUTPUTREQUEST_FILENAME) 
 }
 cat("I have run through the code. Results are on files.\n") 
}
  