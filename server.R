library(shiny)
#library(shinyjs)
library(plyr)
library(dplyr)
library(DataCombine)
library(openxlsx)
library(readxl)
library(reticulate)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggrepel)
#library(airtabler)


options(shiny.maxRequestSize = 150*1024^2)
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
# There is a chance that some high/med/low/zero groups will have no proteins in them.  Low is most likely.
# I set a variable here to set later so that it can be seen in each function.


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
        ext <- parts[length(parts)]
        rn <- parts[length(parts)-1]
        moreparts <- strsplit(rn,"\\/")[[1]]
        fn <- moreparts[length(moreparts)]
        if(ext == "csv"){
          df <- read.csv(files[i], header=TRUE, check.names=FALSE)
        } else if(ext=="txt" || ext=="tsv" || ext=="tab") {
          df <- read.delim(files[i], header=TRUE, check.names=FALSE)
        } else if(ext=="xlsx" || ext=="xls") {
          # use the readxl library instead of openxlsx. openxlsx seems to have
          # a hard time reading xl files output by PD.
          predf <- read_excel(files[i])
          df <- as.data.frame(predf)
        } else {
          validate(need(ext=="csv" || ext=="tsv" || ext=="xlsx" || ext=="xls" || ext=="tab" || ext=="txt", "You have the wrong file extension.  Please see the instructions for possible extensions and associated file types." ))
        }
        incProgress(0.5)
        
        df[is.na(df)] <- "" # R likes NA.  I don't
        # need an id column to allow joining information later. See comment immediately below
        df['ID'] <- seq.int(nrow(df))
        print(fn)
        # call the python script to do the cScIFTING
        # p is a list of data frames -p1 - p4 are the high, med, low, and zero proteinss.
        # p5 -p8 are the H/M/L/Z peptides p9 - p12 are the H/M/L/Z PSMs
        # p13 - Reagent Analysis, p14 - Sequon Analysis and p15 - Specificity
       ## Also, I am only feeding in the two columns cScIFTING uses from the input PD
        # spreadsheet to save memory.  I will need to put the information back into the 
        # PSM tabs in the output spreadsheet.  The ID column is used for this purpose
        p <- cScIFTING( df[c('ID', 'Master Protein Accessions', 'Annotated Sequence')] )

        if( length(p[[3]]) >0 ) {
          lowOK = TRUE
        } else {
          lowOK = FALSE
        }
        
        if(typeof(annotation) == "character"){
#          print("at first annotation reading")
          annotation <<- read.delim("./ref/CIRFESSannotations.tsv", header = TRUE)
        }

        # add the annotation to the proteins (pep and psm, would be straight-forward)
        p[[1]]<-join(p[[1]], annotation, by=c("MPAnoIso"), type="left", match="first")
        p[[2]]<-join(p[[2]], annotation, by=c("MPAnoIso"), type="left", match="first")
        if( lowOK ) {
          p[[3]]<-join(p[[3]], annotation, by=c("MPAnoIso"), type="left", match="first")
        }
        p[[4]]<-join(p[[4]], annotation, by=c("MPAnoIso"), type="left", match="first")
        
        # put the PSM information from the PD file back into the PSM spreadsheet tabs
        p[[9]] <- join(p[[9]], df, type='left', match='first')
        p[[9]] <- p[[9]][, !names(p[[9]]) %in% c("ID")]
        p[[10]] <- join(p[[10]], df, type='left', match='first')
        p[[10]] <- p[[10]][, !names(p[[10]]) %in% c("ID")]
        p[[11]] <- join(p[[11]], df, type='left', match='first')
        p[[11]] <- p[[11]][, !names(p[[11]]) %in% c("ID")]
        p[[12]] <- join(p[[12]], df, type='left', match='first')
        p[[12]] <- p[[12]][, !names(p[[12]]) %in% c("ID")]
        
        # perform the GO term analysis
        if( lowOK ){
          go <- c(p[[1]]$Gene.Ontology..GO...Uniprot., p[[2]]$Gene.Ontology..GO...Uniprot., p[[3]]$Gene.Ontology..GO...Uniprot.)
        } else {
          go <- c(p[[1]]$Gene.Ontology..GO...Uniprot., p[[2]]$Gene.Ontology..GO...Uniprot.)
        }
        tdf <- as.data.frame(table(unlist(sapply(as.character(go),strsplit,split=';'),use.names=FALSE)))
        colnames(tdf)<-c('Gene Ontology Term', 'Count')
        p[[16]]<-tdf[order(-tdf$Count),]

        # perform the keyword analysis
        if( lowOK ) {
          kw <- c(p[[1]]$Keywords..Uniprot., p[[2]]$Keywords..Uniprot., p[[3]]$Keywords..Uniprot.)
        } else {
          kw <- c(p[[1]]$Keywords..Uniprot., p[[2]]$Keywords..Uniprot.)
        }
        tdf <- as.data.frame(table(unlist(sapply(as.character(kw),strsplit,split=';'),use.names=FALSE)))
        colnames(tdf)<-c('Keywords (Uniprot)', 'Count')
        p[[17]]<-tdf[order(-tdf$Count),]

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
        
        tdf <-p[[9]][c('MPA', 'Annotated Sequence')]
        tdf<- rbind(tdf, p[[10]][c('MPA', 'Annotated Sequence')])
        if( lowOK ) {
          tdf<- rbind(tdf, p[[11]][c('MPA', 'Annotated Sequence')])
        }
        colnames(tdf)<-c('Accession', 'Annotated Sequence')
        Replaces <- data.frame(from = c("n", "c", "m"), to = c("N[+1]", "C[+57]", "M[+16]"))
        tdf<-data.frame(tdf, 'PeptideSequence'=sapply(strsplit(as.character(tdf$'Annotated Sequence'),'.',fixed=TRUE), "[", 2))
        tdf$PeptideModifiedSequence<-tdf$PeptideSequence
        tdf <- FindReplace(data = tdf, Var = "PeptideModifiedSequence", replaceData = Replaces,from = "from", to = "to", exact = FALSE)
        p[[18]]<-data.frame(tdf$Accession, tdf$PeptideSequence, tdf$PeptideModifiedSequence)
        colnames(p[[18]])<-c('ProteinName', 'PeptideSequence', 'PeptideModifiedSequence')
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
      highprot <- length(p[[1]][,1])
      medprot <- length(p[[2]][,1])
      if( length(p[[3]]) > 0){
        lowprot <- length(p[[3]][,1])
      } else {
        lowprot<-0
      }
      zeroprot <- length(p[[4]][,1])
      totprot <- highprot + medprot + lowprot + zeroprot
      pcthighprot <- highprot/totprot
      pctmedprot <- medprot/totprot
      pctlowprot <- lowprot/totprot
      pctzeroprot <- zeroprot/totprot

      s <- c( s, c(exp, highprot, round(pcthighprot*100,1), medprot, round(pctmedprot*100,1), lowprot, round(pctlowprot*100,1), zeroprot, round(pctzeroprot*100,1) ))
    }
    df <- data.frame(matrix(s,ncol=9, byrow=TRUE))
    colnames(df) <- c('Experiment', '# High Proteins', '% High Proteins', '# Medium Proteins', '% Medium Proteins', '# Low Proteins', '% Low Proteins', '# Zero Proteins', '% Zero Proteins' )
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
      paste0(input$userfile1$name, '_Veneer.zip')
    },
    content = function(filename) {
      output$readerror <- renderText("")
      size <- length(data_input())+1
      withProgress(message = 'Writing Files', max= size, value = 0, {
        protter <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ProteinName", "PeptideSequence", "PeptideModifiedSequence"))
        d <- data_input()
        for(i in 1:length( d )) {
          p <- d[[i]]
          fn <- names(d)[i]
          cscfile = paste0(fn, '_Veneer.xlsx')
          protterfile = paste0(fn, '_protter.tsv')
          protter <- rbind(protter, p[[18]])
          
          l = list("High Proteins"=p[[1]], "Medium Proteins"=p[[2]], "Low Proteins"=p[[3]], "Zero Proteins"=p[[4]], "High Peptides"=p[[5]], "Medium Peptides"=p[[6]], "Low Peptides"=p[[7]], "Zero Peptides"=p[[8]],  "High PSMs"=p[[9]], "Medium PSMs"=p[[10]], "Low PSMs"=p[[11]], "Zero PSMs"=p[[12]],  "Reagent Analysis"=p[[13]], "Sequon Analysis"=p[[14]], "Specificity"=p[[15]], "GO Terms (Uniprot)"=p[[16]], "Keywords (Uniprot)"=p[[17]])
          write.xlsx(l, cscfile, colNames=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
          
          protterout <- protterfy( p[[18]] )
          write.table(protterout, protterfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
          
          incProgress (1) 
          zip(zipfile=filename, files=c(cscfile, protterfile))
        }
      }) # end of withProgress
    } # end of downloadHandler content argument
  )  # end of downloadHandler 




}) # end of shinyServer
