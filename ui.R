library(shiny)
#library(shinyjs)
#library(RMariaDB)
#library(airtabler)

source('functions.R')


shinyUI(navbarPage("", id="main", theme = "bootstrap.css",
  
  tabPanel( "Veneer",

    tags$head(tags$script(HTML('
      var fakeClick = function(tabName) {
        var dropdownList = document.getElementsByTagName("a");
        for (var i = 0; i < dropdownList.length; i++) {
          var link = dropdownList[i];
          if(link.getAttribute("data-value") == tabName) {
            link.click();
          };
        }
      };
    '))),
    
#    useShinyjs(),        
    tags$img(src="Veneer_OverviewNew.png", width="583px", align="right"),
    h4("Welcome to ", span(class ="text-success", "Veneer"),"!"),
    p(tags$i("Filtering, annotating and curating cell surface protein data")),
    p("Veneer is a web app for analyzing mass spectrometry data from cell surface capture and related methods. 
       Veneer rapidly processes search results to filter out non-specific binders and annotate data with 
       information from >10 resources"),
    p("Veneer was written in Python and R and the web application was developed using the Shiny library. Source 
       code and all reference lookup tables are publicly available at ", tags$a(href="https://github.com/GundryLab/veneer", "GitHub"), "."),

    br(), 
    #     div(style="width:20%;display:block;margin-left:auto;margin-right:auto",
    div(style="width:25%;display:block;margin-right:auto",
    tags$script(type="text/javascript", id="clustrmaps", src="https://clustrmaps.com/map_v2.js?d=6phPuEaX8fEeYstJuTRlTq7TJ2L8rTHD0NRumBE7gQs&cl=ffffff&w=a")    ) 
#   tags$script(type="text/javascript", id="clustrmaps", src="https://cdn.clustrmaps.com/map_v2.js?d=oOpnxVR26WnLiF16Exi-XrxGb3rCX9xwJ4nBPUTKf1E&cl=ffffff&w=a")     )

  ),

  tabPanel( "Instructions",
#    p(tags$i("Before you begin:")),
    p(style="font-size: 17px", tags$i("Before you begin:")),
    p("It is strongly recommended that all users read the ", tags$a(href="Veneer_UserGuide.pdf", "User Guide"), " which contains 
      step-by-step tutorials. The User Guide comprehensively defines the functions and annotations provided in Veneer."
    ),
    p(style="font-size: 17px", tags$i("File Format Requirements:")),
    p("Import file format is critical to success.",  tags$b("Filter and Annotate "), "accepts xlsx files 
      containing a list of protein identifiers (UniProt Accession) and corresponding 
      peptide spectrum matches (PSM). There are four requirements:",
      tags$ol(
        tags$li("One column must be present and contain Uniprot Accessions labeled: ", tags$i("Master Protein Accessions")), 
        tags$li("Another column must be present containing PSM sequence information labeled: ", tags$i("Annotated Sequence")), 
        tags$li("The PSM entries in ", tags$i("Annotated Sequence"),  " must be in the following format: [R].TQDEILFSnSTR.[L], where flanking amino acids are in brackets and deamidation (release of the N-glycan from the asparagine) is denoted as small letter n"),
#        tags$li("Other columns may be present but are not necessary and will be ignored except for being reported back in the output."),
        tags$li("The file cannot exceed 50 MB") 
      ),
      "Optional:",
      tags$ol(
        tags$li("Other columns may be included in the input file, such as other output from Proteome Discoverer.  
                These columns will be ignored except for being reported back in the output"),
        tags$li("Multiple xlsx files can be combined into a zip file.")
      )
    ),
    p(style="font-size: 17px", tags$i("Example files:")),
    p("Examples of correctly formatted files formatted can be downloaded using the 
      links below. For more information, please refer to the User Guide."
    ),
    p("Example File 1 - ", tags$a(href="Veneer Example File 1.xlsx", "Single Experiment")),
    p("Example File 2 - ", tags$a(href="Veneer Example File Zip.zip", "Zip File of Multiple Experiments")),

    p(style="font-size: 17px", tags$i("Conversion to Uniprot Accession IDs:")),
    p("Veneer operates with Uniprot accession identifiers. Bulk conversion of alternate 
       IDs to Uniprot IDs can be performed using the ‘Retrieve/ID mapping tool’ 
       available on the Uniprot website, found ", tags$a(href="https://www.uniprot.org/uploadlists/", "here"), ". Note that conversion between IDs 
       is not always one-to-one. Manual curation of the results from the ID mapping is advisable."
    ),
    p(style="font-size: 17px", tags$i("Species availability:")),
    p("Currently, annotations are available for human. If you 
       have requests for additional species, please contact us.")
    ),
  # Sidebar with a slider input for number of bins
  tabPanel("Filter and Annotate",
    sidebarPanel(
      fileInput("userfile1", "Upload file", multiple =FALSE, buttonLabel = "Browse Files", placeholder = "Select File"),
      p('Upload a single proteomics result file or a zip file of multiple results')
    ),
    mainPanel(
      h4("Retrieve Results"),
      downloadButton("dlAnno", label = "Download"),
      textOutput("readerror"),
      h4("Specificity"),
      tableOutput('dtSummary')
    )
  ),

  
  ##########    References   ##########
  
  tabPanel(
    "References",
    div(
      h4("How to reference ", span(class ="text-success", "Veneer") ),
      p("If you use any of the Veneer tools in your work, please cite the original manuscript:"),
      p("Publication Pending")
      # p("Authors, Title,", 
      #   tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/32053146/", "https://www.ncbi.nlm.nih.gov/pubmed/32053146/"))
    ),
    br(),
    div(
      h4("Publications that cite ", span(class ="text-success", "Veneer") ),
      p("Coming Soon!")
    ),
    br(),
    # div(
    #   h4("Publications that support the ", span(class ="text-success", "SPC Score") ),
    #   tags$ol(
    #     tags$li( "Bausch-Fluck D, et al. (2018) The in silico human surfaceome. Proc Natl Acad Sci U S A 115(46):E10988-E10997."),
    #     tags$li( "da Cunha JP, et al. (2009) Bioinformatics construction of the human cell surfaceome. Proc Natl Acad Sci U S A 106(39):16752-16757"),
    #     tags$li( "Town J, et al. (2016) Exploring the surfaceome of Ewing sarcoma identifies a new and unique therapeutic target. Proc Natl Acad Sci U S A 113(13):3603-3608" ),
    #     tags$li( "Diaz-Ramos MC, Engel P, & Bastos R (2011) Towards a comprehensive human cell-surface immunome database. Immunol Lett 134(2):183-187.")
    #   )
    # ),
    br()
    ),
  
  ##########    Contact   ##########
  
  tabPanel(
    "Contact",
    div(
      h4(span(class ="text-success", "Contact "), "Us!"),
      p("If you have questions or suggestions for additional features, please contact us by email:"),
      p(class="text-info", style="text-indent:1.5em", "rebekah.gundry at unmc.edu"),
      p("Additional cell surface-related information and tools can be found at our growing website:"),
      p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net")
    ),
    br()
  )
  
))
