###################################################################
# EFDMoutput
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
#outputrequests <- "outputrequests.txt"

outputcalc <- function(outputrequests){
  
requests<-readLines(outputrequests)
#Reads the file which includes output requests

nline<-length(requests)
#Number of lines in outputrequests file

for(i in seq(1,nline-1,by=2)){ #1 3 5 ...

  outputrequest<-requests[i]
  #Reads the output request from a textfile
  
  strsplit(outputrequest, " ") -> strlist
  #Splits the string from the textfile to pieces

  nvars <- 0
  k <- 1
  while (strlist[[1]][k] != "by"){
    nvars <- nvars + 1
    k = k+1
  }
  nvars <- nvars - 1 #Minus keyword
  #Number of output variables
  
  howvar <- strlist[[1]][1]
  #The 1st word from the string (TOTAL or PERHA)
  
  outputvar <-strlist[[1]][2:(1+nvars)]

  specact<-FALSE #specifically mentioned activities (there is a dot in outputvar(s))
  if(grepl("[.]",outputvar[1])){
  specact<-TRUE  
  outputv <- simplify2array(strsplit(outputvar,"[.]"))
  #E.g. drain.thin1 -> drain thin1 etc. as columns
  outputvar <- outputv[1,1]
  outputa <- outputv[2,]
  }

  byvar <- strlist[[1]][(3+nvars):length(strlist[[1]])]
  #Output/by-variable

  coeffile<-read.table(requests[i+1],header=TRUE) 
  #The name of the coefficient files

  var.acts<-paste(outputvar,actnames,sep=".")
  #all of the potential ones! and in the same order as general activities

  acts <- intersect(names(coeffile),var.acts)
  #Number of activities for which there are coefficients
  
  if(length(acts) != 0){
   #If there are some activities (it is of type "drain")
    
    for(j in setdiff(var.acts,acts)){          
      coeffile[[j]] <- 0
    } #set the ones not mentioned in the coeffile to 0
    if (specact) for(j in acts){
      if(!j %in% paste(outputvar,outputa,sep=".")) 
      coeffile[[j]] <- 0
    #Activities which are not counted are 0
    }

  output<-merge(rawoutput,coeffile)
  #Combines two data together

  
  tmp <- matrix(0,nrow=nrow(output), ncol=nrofsteps+1)
  #Matrix which has the same number of columns as there 
  #are steps (or the size of variable nrofsteps)
  
  for(k in 1:nact){ 
    table <-as.matrix(subset(output,select=paste0("step", 0:nrofsteps,".", actnames[k])))
    #Chooses the steps to another matrix
    coef <- output[[var.acts[k]]]
    #Chooses the coefficients of activities
    tmp <- sweep(table,1,coef,FUN="*") + tmp    
  } #For ends
    
  out <- cbind(output[factnames], tmp)
    #Combines two data together
  names(out) <- c(factnames,paste0("step",1:(nrofsteps+1)))
    #Gives names to some columns
  nsteps <- paste0("step",1:(nrofsteps))
    
  }else{
   #If there are not any activities
  
   output<-merge(resultstates,coeffile)
   #Combines two data together
   
   table <-as.matrix(subset(output,select=paste0("step",0:nrofsteps)))
   #Chooses the steps to another matrix
   coef <- output[[outputvar]]  
   tmp <- sweep(table,1,coef,FUN="*")
  
   out <- cbind(output[factnames], tmp)
   #Combines two data together
   names(out) <- c(factnames,paste0("step",0:nrofsteps))
   #Gives names to some columns
   nsteps <- paste0("step",0:nrofsteps)
  } #Else ends
  
  area <- resultstates
  #To get the amount of area
  
  if(howvar=="TOTAL"){
    output<-aggregate(x=subset(out,select=nsteps),
                      by=subset(out,select=byvar),FUN=sum) 
    #The amount of output variable
  } else if (howvar=="PERHA"){
     aoutvar<-aggregate(x=subset(out,select=nsteps),
                      by=subset(out,select=byvar),FUN=sum)
    #The amount of output variable
	ha<-aggregate(x=subset(area,select=paste0("step",0:(length(nsteps)-1))),
	              by=subset(area,select=byvar),FUN=sum)
     #The amount of area
 
    #output variable / area
	output <- cbind(aoutvar[,1:length(byvar)],aoutvar[,length(byvar)+(1:length(nsteps))]/ha[,length(byvar)+(1:length(nsteps))])
     #Combines output and byvars together
	names(output)<-c(byvar,nsteps)

    
  }

  if(specact) outputvartext<-paste(c(outputvar,outputa),collapse="_") else
    outputvartext<-outputvar
  outfname<-paste0(howvar,"_",outputvartext, "_by_",paste(byvar,collapse="_"),".txt")
  write.table(output, file=outfname)
    #Saves the results to a text file
  efdmplot(output,outputvartext,byvar,howvar)


} #For-statement ends
} #Function ends

efdmplot <- function(dname,outputvar,byvar,howvar){
  #Draws barplots of data frames given as argument or in a file
  if(is.character(dname)){ 
    output <- read.table(dname) #Reads the data for barplot
  }else{ #dname is presumably a data frame
    output<-dname
  }
  
  #M<-as.matrix(output[-(1:length(byvar))])
  M<-as.matrix(output[setdiff(names(output),byvar)])
  #Data to matrix without the names of byvars
  
  #rownames(M) <- as.character(paste(output[[1]],output[[length(byvar)]]))
  rownames(M) <- do.call(paste,output[byvar])
  #Rownames to data
  
  if(colnames(M)[1]=="step0"){
    vartype <- TRUE #output variable is type volume
  }else{
    vartype <- FALSE #output variable is type drain
  } 
  
  names <- c()
  if(vartype==TRUE){
    names<- paste0("step",0:(nrofsteps))
  }else{ 
    for(j in 1:(nrofsteps)){
      step <- paste0("step",j-1,"-step",j)
      names <- c(names,step)
    }
  }
  #Naming the bars of barplots
  #varype=TRUE for the outputs of the type "volume"
  #vartype=FALSE for the outputs of the type "drain"
  
  jpeg(paste0("barplot_of_totals_of_",howvar,"_of_",outputvar,"_by_",
              paste(byvar,collapse="_"),".jpg"))
  barplot(M, legend.text=TRUE, col=c(1:nrow(M)), ylab=outputvar, names.arg=names)
  dev.off()
  #Barplot with totals
  
  jpeg(paste0("barplot_of_",howvar,"_of_",outputvar,"_by_",
              paste(byvar,collapse="_"),".jpg"))
  barplot(M,beside=TRUE, legend.text=TRUE, col=c(1:nrow(M)), ylab=outputvar, names.arg=names)
  dev.off()
  #Barplot without totals
  
} #Function ends

