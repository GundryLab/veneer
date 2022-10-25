library(shiny)
#library(shinyjs)
library(plyr)
library(dplyr)
library(DataCombine)
library(openxlsx)
library(readxl)
library(reticulate)
library(reshape2)
#library(RMariaDB)
library(data.table)
library(ggplot2)
library(ggrepel)


options(shiny.maxRequestSize = 50*1024^2)
source('functions.R')

###############################################################################
#
# The following three lines must be uncommneted in the shinyapps.io version
# in order for it to work there. The .Rprofile must also be in place (see git
# repo).  The name of the virtual environment (below) must be the same as
# the VIRTUALENV_NAME variable in the .Rprofile. 
#
# In theory, with the .Rprofile, we shouldn't need to comment out the three 
# virtual environment lines (below) to work locally, but I get complaints about
# pip on (some of) my local machines when installing pandas to the virtual 
# environment...
###############################################################################

virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "pandas")
use_virtualenv("r-reticulate", required = T)

# this is how we access the python script
source_python('./ref/cScIFTING.py')

# put this here so that it only runs once per instance on shinyapps.io
#annotation <- read.delim("./ref/CIRFESSannotations.tsv", header = TRUE)
annotation <- "a"

#set up data structures for database
protein_fields <- c('MPA','MPAnoIso', 'numPep', 'numPSM', 'psmExclusive', 'pctExclusive', 'PSMwSCM', 'pctPSMwSCM', 'SCMonePSM')
peptide_fields <- c('MPA','MPAnonSplit', 'MPAnoIso', 'pepSeq', 'numMPA', 'pepPSM', 'pepPSMwSCM', 'pctPepPSMwSCM', 'hasSCM')
psm_fields <- c('MPA','MPAnonSplit', 'MPAnoIso', 'pepSeq', 'annSeq', 'numMPA', 'hasSCM')



shinyServer(function(input, output, session) {
  

  data_input <- eventReactive(input$userfile1, {
    
    output$readerror <- renderText("")
    
    d <- list()
    files <- c()
    
    parts <- strsplit(input$userfile1$name, "\\.")[[1]]
    ext <- parts[length(parts)]
    if(ext == "zip"){
      files <- processZip(input$userfile1$datapath)
    } else {
      files <-processSingle(input$userfile1$datapath, input$userfile1$name)
    }
    for (i in seq_along(files)) {
      withProgress(message = paste0('Reading and Processing ', files[i]), value = 0, {
#        print(files[i])
        parts <- strsplit(files[i], "\\.")[[1]]
#        print(parts)
        ext <- parts[length(parts)]
        rn <- parts[length(parts)-1]
#        print(rn)
        moreparts <- strsplit(rn,"\\/")[[1]]
#        print(moreparts)
        fn <- moreparts[length(moreparts)]
        if(ext == "csv"){
          df <- read.csv(files[i], header=TRUE)
        } else if(ext=="txt" || ext=="tsv" || ext=="tab") {
          df <- read.delim(files[i], header=TRUE)
        } else if(ext=="xlsx" || ext=="xls") {
          # use the readxl library instead of openxlsx. openxlsx seems to have
          # a hard time reading xl files output by PD.
          # df <- read.xlsx(files[i], sep.names=" ", na.strings="")
          predf <- read_excel(files[i])
          df <- as.data.frame(predf)
        } else {
          validate(need(ext=="csv" || ext=="tsv" || ext=="xlsx" || ext=="xls" || ext=="tab" || ext=="txt", "You have the wrong file extension.  Please see the instructions for possible extensions and associated file types." ))
        }
        incProgress(0.5)
        #      withProgress(message = 'Reading and Processing Uploaded Data', value = 0, {
        #      print(input$userfile1$datapath)
        #      xl <- read.xlsx(input$userfile1$datapath, sep.names=" ", na.strings="")
        #      print("read file")
        df[is.na(df)] <- "" # R likes NA.  I don't
        
        # call the python script to do the cScIFTING
        # p is a list of data frames -p1, p2 are SCM and NSB proteins
        # p3, p4 are SCM nad NSB peptides. p5, p6 are the SCM and NSB PSMs
        # p7 - MIAPE, p8 - Reagent Analysis, p9 - Motif Analysis p10 - Specificity
        p <- cScIFTING(df)
        
        if(typeof(annotation) == "character"){
#          print("at first annotation reading")
          annotation <<- read.delim("./ref/CIRFESSannotations.tsv", header = TRUE)
        }
        # add the annotation to the proteins (pep and psm, would be straight-forward)
        p[[1]]<-join(p[[1]], annotation, by=c("MPAnoIso"), type="left", match="first")
#        print(colnames(p[[1]]))
        p[[2]]<-join(p[[2]], annotation, by=c("MPAnoIso"), type="left", match="first")
        
        # perform the GO term analysis
        tdf <- as.data.frame(table(unlist(sapply(as.character(p[[1]]$Gene.Ontology..GO...Uniprot.),strsplit,split=';'),use.names=FALSE)))
        colnames(tdf)<-c('Gene Ontology Term', 'Count')
        p[[11]]<-tdf[order(-tdf$Count),]
        
        # perform the keyword analysis
        tdf <- as.data.frame(table(unlist(sapply(as.character(p[[1]]$Keywords..Uniprot.),strsplit,split=';'),use.names=FALSE)))
        colnames(tdf)<-c('Keywords (Uniprot)', 'Count')
        p[[12]]<-tdf[order(-tdf$Count),]
        
        ###################
        # Make the Protter output.
        # Start with Accession and Annotated Sequence
        # Make a new field called PeptideSequence by removing the extra sequence 
        # from Annotated Sequence by using strsplit on the '.'  that indicated the
        # start and end of the extra sequence.  Then use the FindReplace function
        # from the DataCombine library (holy god would it have been hard if that
        # didn't exist - string manipulation in R sucks rocks) to change N to [N+1]
        # and C to C[+57].  FindReplace uses a dataframe to determine what to find
        # and what to replace it with.  that dataframe 'Replaces' is for.
        ###################
        
        tdf<-data.frame(p[[5]]$MPA, p[[5]]$'Annotated Sequence')
        colnames(tdf)<-c('Accession', 'Annotated Sequence')
        Replaces <- data.frame(from = c("n", "c", "m"), to = c("N[+1]", "C[+57]", "M[+16]"))
        tdf<-data.frame(tdf, 'PeptideSequence'=sapply(strsplit(as.character(tdf$'Annotated Sequence'),'.',fixed=TRUE), "[", 2))
        tdf$PeptideModifiedSequence<-tdf$PeptideSequence
        tdf <- FindReplace(data = tdf, Var = "PeptideModifiedSequence", replaceData = Replaces,from = "from", to = "to", exact = FALSE)
        p[[13]]<-data.frame(tdf$Accession, tdf$PeptideSequence, tdf$PeptideModifiedSequence)
        colnames(p[[13]])<-c('ProteinName', 'PeptideSequence', 'PeptideModifiedSequence')
        
        # input for CIRFESS
        p[[14]]<-data.frame(p[[1]]$MPA)
        colnames(p[[14]])<-c('MPA')
        p[[15]]<-p[[1]][which(p[[1]]$cirfessScore!=0|p[[1]]$SPC==3|p[[1]]$SPC==4|p[[1]]$"Signal.Peptide..PredSi."=="Yes"|p[[1]]$"Signal.Peptide..SignalP."=="Yes"|p[[1]]$"Signal.Peptide..Phobius."=="Yes"),]
#        p[[15]]<-p[[1]][which(p[[1]]$"Signal.Peptide..PredSi."=="Yes"),]
      })
      d[[fn]] <- p
    }
    return(d)
  })

  
  ####################################
  # display specificity mini-report
  ####################################
  output$dtSummary <- renderTable({
    req(input$userfile1)
    s<-c()
    for(i in 1:length( data_input() )) {
      p<-data_input()[[i]]
      exp <- names(data_input())[i] 
      scm <- length(p[[1]][,1])
      nsb <- length(p[[2]][,1])
      pctProt <- p[[10]][2,3]
      pctPSM <- p[[10]][5,3]
      pngasef <- p[[8]][1,1]
      strep <- p[[8]][2,1]
      tryp <- p[[8]][3,1]
      s <- c( s, c(exp, scm, nsb, pctProt, pctPSM, pngasef, strep, tryp) )
    }
    df <- data.frame(matrix(s,ncol=8, byrow=TRUE))
    colnames(df) <- c('Experiment', 'SCM Proteins', 'NSB Proteins', '%Protein Specificity', '%PSM Specificity', 'PNGaseF PSMs', 'Streptavidin PSMs', 'Trypsin PSMs' )
    df
  })
  
  
  protterfy <- function(df) {
    prots <- distinct(df['ProteinName'])
    protter <-c()
    for( i in 1:length(prots[,1]) ) {
      prot<-prots[i,]
      output <- paste0(prot, ' ', paste0(distinct(subset(df,ProteinName==prot))[,3], collapse=','))
      protter <- c(protter, output)
    }
    protter<-as.data.frame(protter)
    return(protter)
  }
  

  output$dlAnno <- downloadHandler(
    filename = function() {
      paste0(input$userfile1$name, '.zip')
    },
    content = function(filename) {
      output$readerror <- renderText("")
      withProgress(message = 'Writing Files', max=length(data_input())+1, value = 0, {
        protter <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ProteinName", "PeptideSequence", "PeptideModifiedSequence"))
        d <- data_input()
        for(i in 1:length( d )) {
          p <- d[[i]]
          fn <- names(d)[i]
          cscfile = paste0(fn, '_Veneer.xlsx')
          protterfile = paste0(fn, '_protter.tsv')
          protter <- rbind(protter, p[[13]])
#          l = list("SCM Proteins"=p[[1]], "SCM - Filtered"=p[[15]], "NSB Proteins"=p[[2]], "SCM Peptides"=p[[3]], "NSB Peptides"=p[[4]], "SCM PSMs"=p[[5]], "NSB PSMs"=p[[6]], "MIAPE"=p[[7]], "Reagent Analysis"=p[[8]], "Motif Analysis"=p[[9]], "Specificity"=p[[10]], "GO Terms"=p[[11]], "Keywords"=p[[12]])
          l = list("SCM Proteins"=p[[1]], "SCM - Filtered"=p[[15]], "NSB Proteins"=p[[2]], "SCM Peptides"=p[[3]], "NSB Peptides"=p[[4]], "SCM PSMs"=p[[5]], "NSB PSMs"=p[[6]], "Reagent Analysis"=p[[8]], "Motif Analysis"=p[[9]], "Specificity"=p[[10]], "GO Terms (Uniprot)"=p[[11]], "Keywords (Uniprot)"=p[[12]])
#          write.xlsx(l, cscfile, colNames=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
          write.xlsx(l, cscfile, colNames=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
          protterout <- protterfy( p[[13]] )
          write.table(protterout, protterfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

          zip(zipfile=filename, files=c(cscfile, protterfile))
        }
      }) # end of withProgress
    } # end of downloadHandler content argument
  )  # end of downloadHandler 




}) # end of shinyServer
