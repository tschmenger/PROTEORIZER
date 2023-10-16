library(shiny)
library(shinyjs)
library(rmarkdown)
#library(gmailr)
library(NGLVieweR)
library(dplyr)
library(shinycssloaders)
library(htmltools)
library(DT)
library(stringr)
library(shinyFeedback)

#######################################################################################################################################

############# Working on the Input and Result Generation
###### This is to make the input easier
idmapper <- read.delim("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoX/Lookups/Uniprot_GeneNames_Reviewed.tsv", header= T, sep = "\t", stringsAsFactors = FALSE)
name_uniprot <- setNames(idmapper[,"Entry"],idmapper[,"Gene_Names_primary"])

#uniprot_name <- setNames(idmapper[,"Gene_Names"],idmapper[,"Entry"])
col3_split <- strsplit(as.character(idmapper$Gene_Names)," ")
col1_vals <- rep(idmapper$Entry,lengths(col3_split))
col3_vals <- unlist(col3_split)
uniprot_name <- setNames(col3_vals,col1_vals)



server <- function(input, output, session) {
  shinyjs::hide("filter")
  file_ready <- reactiveVal(FALSE)
  structure_ready <- reactiveVal(FALSE)
  resfile_paths <- reactiveVal()
  alignmentpathway <- reactiveVal()
  
  observeEvent(input$example, {
    updateTextInput(session, "request",value="RHOA/Y34C,L191A,K140A")
  })
  
  observeEvent(input$dataset_input, {
    if (input$dataset_input == "Humsavar") {
      showFeedbackWarning(
        inputId = "dataset_input",
        text = "No scoring can be provided since parts of this dataset were used for training."
      )  
    } else if (input$dataset_input == "Both") {
      showFeedbackWarning(
        inputId = "dataset_input",
        text = "No scoring can be provided since parts of this dataset were used for training."
      )} 
    else {
      hideFeedback("dataset_input")
    }})
  
  observeEvent(input$submit, {
    ### reset file_ready to false
    file_ready(FALSE)
    shinyjs::show("spinner1")
    shinyjs::show("spinner2")
    shinyjs::show("spinner3")
    
    ### reset any already-rendered results
    output$mytable1 <- NULL
    output$image_1 <- NULL
    output$image_2 <- NULL
    
    
    # Get user input
    request <- input$request
    option <- input$option_input
    datensatz <- input$dataset_input
    filepathus <- input$file_input$datapath
    if (!is.null(filepathus)){
      timedate <- format(Sys.time(),"%Y%m%d_%H%M%S")
      customfile_name <- paste0(timedate,"_",input$file_input$name)
      customfilepath <- file.path("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/User_alignments/",
                                  customfile_name)
      file.copy(filepathus, customfilepath)
    }
    else{filepathus <- "FALSE"}
    
    
    # Write user input to file
    write_input_to_file <- function(request,option,file) {
      target_dir <- "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/User_Inputs"
      timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
      filename <- paste0("user_input_",timestamp,".txt")
      filepath <- file.path(target_dir,filename)
      write.table(data.frame(request,option,file),file = filepath, sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
    }
    write_input_to_file(request,option,filepathus)
    
    that_is_it <- list.files(path = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/User_Inputs")
    foldernumber <- length(that_is_it)
    if (foldernumber > 100){
      deletecommand <- "/home/bq_tschmenger/anaconda2/bin/python /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/CleanerOperation.py"
      system(deletecommand)
    }
    tryCatch({
      DF <- read.table(text = request, sep = "/", as.is = TRUE)
      raw_id <- toupper(DF[[1,1]])
      mutations <- toupper(DF[[1,2]])
      if (is.na(uniprot_name[raw_id])){
        genename <- raw_id
        uniprot <- name_uniprot[raw_id]
      } else {
        genename <- uniprot_name[raw_id]
        uniprot <- raw_id
      }
      # Display success message
      showModal(modalDialog("Success!", "Your request has been submitted. 
                            Head over to the \"Explore the results\" tab and wait for your request to finish.", 
                            easyClose = TRUE, 
                            size = "l",
                            title = "Request submitted"))
    },error = function(e){
      uniprot <- "none"
      mutations <- "none"
      genename <-"none"   
      # Display error message
      showModal(modalDialog("Error!", "Your request was incomplete or erroneous. Please check the input format. The app will automatically reload.", 
                            easyClose = FALSE,
                            size = "l",
                            title = "Error"))
      Sys.sleep(4)
      session$reload()
    }
    )
    
    # Execute shell command to call external Python script
    if (filepathus != "FALSE"){
      customalignfile = paste("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/User_alignments/",
                              customfile_name,
                              sep = "")
      
      cmd <- paste("/home/bq_tschmenger/anaconda2/bin/python /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/STABLE_Proteorizer_alpha_20231016.py",uniprot,mutations,genename,"R_Submissions FALSE",customalignfile,"FALSE","FALSE",option,datensatz,foldernumber, sep=" ")
    }
    else {
      cmd <- paste("/home/bq_tschmenger/anaconda2/bin/python /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/STABLE_Proteorizer_alpha_20231016.py",uniprot,mutations,genename,"R_Submissions FALSE",filepathus,"FALSE","FALSE",option,datensatz,foldernumber, sep=" ")
    }
    
    #cat(cmd)
    tryCatch({
      system(cmd)},
      error =function(e){
        if (datensatz!="Uniprot"){
          showModal(modalDialog("Not enough data available.", easyClose = F))  
        }})  
    
    stringler <- paste(foldernumber,uniprot,genename,sep="-")
    correct_path <- file.path("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/R_Submissions",
                              stringler)
    
    tryCatch({
      graphical_path <- file.path(correct_path,
                                  "ClusterPlot_unlim.svg")},
      error = function(e){graphical_path=""})
    
    tryCatch({  
      if (length(grep(',', as.character(mutations)))){
        alignment_first <- file.path(correct_path,
                                     "AnnotatedAlignment_30000_15.svg")
      }
      else {
        alignment_first <- file.path(correct_path,
                                     "AnnotatedAlignment_30_15.svg")
      }
    },
    error = function(e){alignment_first=""})
    
    tryCatch({
      uniprotler_raw <- strsplit(stringler,"-")
      uniprotler <- uniprotler_raw[[1]][2]
      realalphafold <- paste("AF-",uniprotler,"-F1-model_v3.pdb",sep="")  
      realstrucpath <- file.path(correct_path,realalphafold)
    },
    error = function(e){realstrucpath=""})
    
    if (option == "RandomWalk"){
      tryCatch({   
        resultfile <- file.path(correct_path,
                                "FinalResults_unlim_conservation_scored.txt")},
        error = function(e){resultfile=""})
      
    }
    else {
      tryCatch({  
        resultfile <- file.path(correct_path,
                                "FinalResults_unlim_conservation_scored_hclust.txt")},
        error = function(e){resultfile=""}) 
    }
    resfile_paths(list(
      uni_identifier = uniprot,
      mutatoos = mutations,
      onepathforall = correct_path,
      tablepath = resultfile,
      clusterspath = graphical_path,
      annotated = alignment_first,
      strupath = realstrucpath))
    alignmentpathway(list(annotated = alignment_first))
    
    file_ready(TRUE)
    
  })
  
  observeEvent(file_ready(),{
    if (file_ready()){
      shinyjs::show("filter")
      kopftext <- c("Cluster number. Filter to show positions close to each other.",
                    "Source of the data, one of [Input, Uniprot] for the protein of interest or [Orthos] when mapped from similar proteins.",
                    "Amino acid residue number based on the protein of interest",
                    "Prediction outcome. The prediction result is meaningful if BAYES and Predictor transformed values combined are >= 0.25 (medium) or > 0.5 (high).",
                    "Score derived from naive bayesian. Score transformation: 10 to 12.5 -> 0.125|12.5 to 15 -> 0.25|>15 -> 0.375.",
                    "Score calculated using a random forest predictor. Score transformation: 0.5 to 0.7 -> 0.125|0.7 to 0.9 -> 0.25|> 0.9 -> 0.375",
                    "Sequence identity based on the alignment.",
                    "Amino acid type identity based on the alignment. Types: hydrophobic(A,I,L,M,F,W,V,C), positive(K,R), negative(E,D), polar(N,Q,S,T), aromatic(H,Y).",
                    "Functional information based on prior knowledge.",
                    "COSMIC count.",
                    "gnomAD counts (heterozygous)",
                    "gnomAD counts (homozygous)",
                    "Mechismo results on protein interaction perturbations. For more visit mechismo3.russelllab.org.")    
      headerCallback <- c("function(thead, data, start, end, display){",
                          sprintf(" var tooltips = [%s];", toString(paste0("'",kopftext,"'"))),
                          " for(var i = 1; i <= tooltips.length; i++){",
                          "   $('th:eq('+i+')',thead).attr('title',tooltips[i-1]);",
                          " }",
                          "}")
      tryCatch({
        ### preparing the table
        theresultsfile <- read.delim(resfile_paths()$tablepath, header = T, sep ="\t", quote="")
        #View(theresultsfile)
        betterorder <- c("Clusternumber","Data_Source","Position","Verdict","BAYES","Predictor","SeqIdent%","TypeIdent%","Functional_Information","COSMIC","gnomAD_Het","gnomAD_Hom","Mechismo_Predictions")
        colnames(theresultsfile)<- c("Clusternumber","Data_Source","Position","Functional_Information","Mechismo_Predictions","SeqIdent%","TypeIdent%","BAYES","COSMIC","gnomAD_Het","gnomAD_Hom","Predictor","Verdict")
        theresultsfile <- theresultsfile[, betterorder]
        theresultsfile[, 'SeqIdent%'] <- as.integer(theresultsfile[, 'SeqIdent%'])
        theresultsfile[, 'TypeIdent%'] <- as.integer(theresultsfile[, 'TypeIdent%'])
        result_dt <- datatable(theresultsfile,
                               options =list(headerCallback=JS(headerCallback)),
                               filter = list(position = 'top', clear = FALSE, plain = TRUE
                               )) %>% formatStyle("Verdict", backgroundColor = styleEqual(c("Impact(high)","Impact(medium)","Impact(low)","No_Impact","-"),c("mediumspringgreen","lightskyblue","#EB7B3B","white","white")))%>%
          formatStyle("SeqIdent%", background = styleColorBar(range(theresultsfile$"SeqIdent%"), 'lightblue'),backgroundSize = '98% 88%',
                      backgroundRepeat = 'no-repeat',backgroundPosition = 'center')%>%
          formatStyle("TypeIdent%", background = styleColorBar(range(theresultsfile$"TypeIdent%"), 'lightblue'),backgroundSize = '98% 88%',
                      backgroundRepeat = 'no-repeat',backgroundPosition = 'center')
        output$mytable1 <- renderDataTable({result_dt}) 
        
      },error =function(e){
        showModal(modalDialog("Table Error!", "Tabular results could not be generated. Click & proceed to the remaining results.", easyClose = TRUE))
      })
      
      ### preparing the images
      tryCatch({
        output$image_1 <- renderImage({list(
          src = resfile_paths()$clusterspath,
          contentType = 'image/svg+xml'
        )
        },deleteFile = F)},
        error =function(e){
          showModal(modalDialog("Plotting Error!", "Clusterplot could not be generated. Click & proceed to the remaining results.", easyClose = TRUE))  
        })
      
      tryCatch({
        output$image_2 <- renderUI({
          HTML(htmltools::includeHTML(resfile_paths()$annotated))
        })},
        error =function(e){
          showModal(modalDialog("Plotting Error!", "Annotated alignment could not be generated. Click & proceed to the remaining results.", easyClose = TRUE))  
        })
      output$download_alignment <- downloadHandler(
        filename = function() {
          paste("Alignment_",Sys.Date(),".svg",sep="")
        },
        content = function(file) {
          cat(alignmentpathway()$annotated)
          svg_file_path <- alignmentpathway()$annotated
          svg_content <- readLines(svg_file_path)
          writeLines(svg_content,file)
        }
      )
      output$download_table <- downloadHandler(
        filename = function() {
          paste("FuncInfo_",Sys.Date(),".csv",sep="")
        },
        content = function(con) {
          write.csv(theresultsfile, con)
        }
      )
      
      ##### Bringing the structure to life here
      ## Initially generating the 3D structure
      if (nchar(resfile_paths()$strupath)>5){
        output$structure <- renderNGLVieweR({
          NGLVieweR(resfile_paths()$strupath) %>%
            addRepresentation("cartoon",
                              param = list(name = "cartoon",
                                           colorScheme = "bfactor")
            ) %>%
            stageParameters(backgroundColor = "black") %>%
            setQuality("high") %>%
            setFocus(0) %>%
            setSpin(FALSE)
        })
        structure_ready(TRUE)}
    }})
  
  ### Changing the alignment 
  observeEvent(input$filter, {
    tryCatch({
      windowgroesse <- input$windo
      sequenzenanzeige <- input$topseqs
      sequenzdatei <- file.path(resfile_paths()$onepathforall,"UsedSequences_unlim.fasta")
      positionsdatei <- file.path(resfile_paths()$onepathforall, "positiondictionary.txt")
      featuredatei <-file.path(resfile_paths()$onepathforall,"featurefile.txt")
      clusterdatei <- file.path(resfile_paths()$onepathforall,"clusterfile.txt")
      translatione <- file.path(resfile_paths()$onepathforall,"translationfile.txt")
      commando <- paste("/home/bq_tschmenger/anaconda2/bin/python /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/STABLE_Annotate_Alignment_V8_Proteorizer.py",
                        resfile_paths()$uni_identifier,
                        resfile_paths()$mutatoos,
                        as.character(windowgroesse),
                        sequenzdatei,
                        positionsdatei,
                        featuredatei,
                        as.character(sequenzenanzeige),
                        clusterdatei,
                        translatione,
                        sep=" ")
      setwd(resfile_paths()$onepathforall)
      system(commando)
      if (length(grep(',', as.character(resfile_paths()$mutatoos)))){
        windowgroesse <- "30000"
      }
      resultalignment <- paste("AnnotatedAlignment_",
                               as.character(windowgroesse),
                               "_",
                               as.character(sequenzenanzeige),
                               ".svg",
                               sep="")
      resultalignmentpath <- file.path(resfile_paths()$onepathforall,
                                       resultalignment)
      alignmentpathway(list(annotated = resultalignmentpath))
      
      output$image_2 <- renderUI({
        HTML(htmltools::includeHTML(alignmentpathway()$annotated))
      })
    },
    error =function(e){
      showModal(modalDialog("Plotting Error!", "Make sure you correctly submitted your input.", easyClose = TRUE))  
    })
  })
  
  
  ### Resetting the structure
  observeEvent(input$reset, {
    if (structure_ready()){
      output$structure <- renderNGLVieweR({
        NGLVieweR(resfile_paths()$strupath) %>%
          addRepresentation("cartoon",
                            param = list(name = "cartoon",
                                         colorScheme = "bfactor")
          ) %>%
          stageParameters(backgroundColor = "black") %>%
          setQuality("high") %>%
          setFocus(0) %>%
          setSpin(FALSE)
      })}})
  ### changing structure background
  observeEvent(input$color, {
    tryCatch({
      NGLVieweR_proxy("structure") %>%
        updateStage(param = list('backgroundColor' = isolate(input$color))
        )})
    
  })
  ### taking a Snapshot of the structure
  observeEvent(input$snapshot, {
    if (structure_ready()){
      NGLVieweR_proxy("structure") %>%
        snapShot("Snapshot.png", param = list(antialias = TRUE,
                                              trim = TRUE,
                                              transparent = TRUE,
                                              scale = 1))
    }})
  # # ### Adding residue highlights to the structure
  observeEvent(input$add, {
    if (structure_ready()){
      tryCatch({
        colorpicker_inputdata = 1
        colorlist_inputdata <- c("red","blue","green","brown","purple","yellow","orange")
        ### Getting the information that needs to be highlighted on the structure
        tryCatch({
          tablepath_inputdata <- file.path(resfile_paths()$tablepath)
          partfilus_inputdata <- read.delim(tablepath_inputdata, header=TRUE, quote="")
          inputclusters_idata <- c()
          inputpositions_idata <- c()
          positions_to_highlight_idata <- c()
          for (row in 1:nrow(partfilus_inputdata))
          {
            clusterinfo_idata <- partfilus_inputdata[row,"Clusternumber"]
            clusterinfo_idata <- gsub(" ", "", clusterinfo_idata)
            positional <- partfilus_inputdata[row,"Position"]
            positional <- gsub(" ","",positional)
            sourcus <- partfilus_inputdata[row,"Data_Source"]
            sourcus <- gsub(" ","",sourcus)
            ### If I know a certain cluster number corresponds to a clustered input position
            if (clusterinfo_idata %in% inputclusters_idata){
              newinfo <- paste(clusterinfo_idata,positional,sep="/") ### X/Y, where X is the clusternumber and Y the residue number
              positions_to_highlight_idata <- c(positions_to_highlight_idata,newinfo)
            }
            else {}
            
            if (sourcus == "Input"){
              if (clusterinfo_idata != "NoClusterMembership"){
                if (nchar(clusterinfo_idata)>=1){
                  goodtoknow <- paste(clusterinfo_idata,positional,sep="/") ### X/Y, where X is the clusternumber and Y the residue number
                  inputpositions_idata <- c(inputpositions_idata, goodtoknow)
                  inputclusters_idata <- c(inputclusters_idata, clusterinfo_idata)}
                else{}
              }
            }
            else {}
          }
        },
        error = function(e) {
          return()})
        
        for (mutation in inputpositions_idata){ ### "3/127" "2/97"
          realmut <- strsplit(mutation,"/")[[1]][2]
          realclust <- strsplit(mutation,"/")[[1]][1]
          NGLVieweR_proxy("structure") %>%
            addSelection("ball+stick",
                         param =
                           list(
                             name = realmut,
                             sele = isolate(as.character(gsub('\\D+','', realmut))),
                             colorValue = isolate(nth(colorlist_inputdata,colorpicker_inputdata))
                           ))
          for (position in positions_to_highlight_idata){ ### "3/128" "2/98"  "3/127" "2/98"  "3/128"
            posclust <- strsplit(position,"/")[[1]][1]
            pospos <- strsplit(position,"/")[[1]][2]
            if (posclust == realclust){
              NGLVieweR_proxy("structure") %>%
                addSelection("ball+stick",
                             param =
                               list(
                                 name = mutation,
                                 sele = isolate(as.character(gsub('\\D+','', pospos))),
                                 colorValue = isolate(nth(colorlist_inputdata,colorpicker_inputdata))
                               )) %>%
                addSelection("label",
                             param = list(
                               sele = realmut,
                               labelType = "format",
                               labelFormat = "[%(resname)s]%(resno)s", # or enter custom text
                               labelGrouping = "residue", # or "atom" (eg. sele = "20:A.CB")
                               color = isolate(nth(colorlist_inputdata,colorpicker_inputdata)),
                               fontFamiliy = "sans-serif",
                               xOffset = 1,
                               yOffset = 0,
                               zOffset = 0,
                               fixedSize = TRUE,
                               radiusType = 1,
                               radiusSize = 1.5, # Label size
                               showBackground = FALSE
                             ))
            }else {}
            
          }
          colorpicker_inputdata = colorpicker_inputdata + 1
        }
        
      })}})
  
  
  
  
  #############   #############   #############   #############   #############   #############   ############# 
  ############# Working on the PreLoaded data
  observeEvent(input$dataset,{
    parentfolder <- "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/"
    #parentfolder <- "/srv/shiny-server/proteorizer/"
    output$mytable1_discover <- NULL
    output$image_1_discover <- NULL
    output$image_2_discover <- NULL
    
    if (nchar(input$dataset)>=2){
      directions <- basename(list.dirs(path=file.path(parentfolder, "www/",input$dataset,"/")))
      choiceis <- c("",directions)
      
      updateSelectInput(session,'proteincase',
                        choices =choiceis,
                        select = "",
                        label = "Select case")
      output[["select"]]<- NULL}
    
  })
  #############
  observeEvent(input$proteincase,{
    parentfolder <- "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/"
    #parentfolder <- "/srv/shiny-server/proteorizer/"
    
    if (nchar(input$dataset)>=2){
      #        print("step1")
      if (as.character(input$proteincase) != as.character(input$dataset)){
        #          print("step2")
        if (input$proteincase != ""){
          
          if (input$dataset == "VUS_Humsavar_RW"){ #if (input$dataset == "Proteorizer_RW_VUS_Humsavar"){
            completeTABLEfilepath <- file.path(parentfolder,"www",input$dataset,input$proteincase,"FinalResults_unlim_conservation_scored.txt")}
          else {completeTABLEfilepath <- file.path(parentfolder, "www",input$dataset,input$proteincase,"FinalResults_unlim_conservation_scored_hclust.txt")}
          
          completeCLUSTERpath <- file.path(parentfolder, "www",input$dataset,input$proteincase,"ClusterPlot_unlim.svg")
          
          completeALIGNMENTpath <- file.path(parentfolder, "www",input$dataset,input$proteincase,"AnnotatedAlignment_30_15.svg")
          if (file.exists(completeALIGNMENTpath)== T){}
          else {
            completeALIGNMENTpath <- file.path(parentfolder, "www",input$dataset,input$proteincase,"AnnotatedAlignment_30000_15.svg")
          }
          
          ### preparing the table
          thediscovertable <- read.delim(completeTABLEfilepath, header = T, sep ="\t", quote="")
          betterorder <- c("Clusternumber","Data_Source","Position","Verdict","BAYES","Predictor","SeqIdent%","TypeIdent%","Functional_Information","COSMIC","gnomAD_Het","gnomAD_Hom","Mechismo_Predictions")
          colnames(thediscovertable)<- c("Clusternumber","Data_Source","Position","Functional_Information","Mechismo_Predictions","SeqIdent%","TypeIdent%","BAYES","COSMIC","gnomAD_Het","gnomAD_Hom","Predictor","Verdict")
          thediscovertable <- thediscovertable[, betterorder]
          thediscovertable[, 'Clusternumber'] <- str_replace(thediscovertable[, 'Clusternumber'], "NoClusterMembership", "-")
          thediscovertable[, 'Clusternumber'] <- as.integer(thediscovertable[, 'Clusternumber'])
          thediscovertable[, 'Verdict'] <- as.character(thediscovertable[, 'Verdict'])
          thediscovertable[, 'SeqIdent%'] <- as.integer(thediscovertable[, 'SeqIdent%'])
          thediscovertable[, 'TypeIdent%'] <- as.integer(thediscovertable[, 'TypeIdent%'])
          result_discover <- datatable(thediscovertable,
                                       filter = list(position = 'top', clear = FALSE, plain = TRUE
                                       )) %>% formatStyle("Verdict", backgroundColor = styleEqual(c("Impact(high)","Impact(medium)","Impact(low)","No_Impact","-"),c("mediumspringgreen","lightskyblue","#EB7B3B","white","white")))%>%
            formatStyle("SeqIdent%", background = styleColorBar(range(0,100), 'lightblue'),backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',backgroundPosition = 'center')%>%
            formatStyle("TypeIdent%", background = styleColorBar(range(0,100), 'lightblue'),backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',backgroundPosition = 'center')
          output$mytable1_discover <- renderDataTable({result_discover}) 
          
          output$download_table_discover <- downloadHandler(
            filename = function() {
              paste("FuncInfo_",Sys.Date(),".csv",sep="")
            },
            content = function(con) {
              write.csv(thediscovertable, con)
            }
          )
          
          output$image_1_discover <- renderImage({list(
            src = completeCLUSTERpath,
            contentType = 'image/svg+xml'
          )
          },deleteFile = F)
          
          output$image_2_discover <- renderImage({list(
            src = completeALIGNMENTpath,
            contentType = 'image/svg+xml',
            height = 700
            
          )
          },deleteFile = F)
          
          
          #########################################
          ### Initially generating the 3D structure
          proteinID_raw <- strsplit(input$proteincase,"-")
          proteinID <- proteinID_raw[[1]][1]
          alphafoldfile_discover <- paste("AF-",proteinID,"-F1-model_v3.pdb",sep="")  
          structurepath_discover <- file.path(parentfolder,"www",input$dataset,input$proteincase,alphafoldfile_discover)
          output$structure_discover <- renderNGLVieweR({
            NGLVieweR(structurepath_discover) %>%
              addRepresentation("cartoon",
                                param = list(name = "cartoon",
                                             colorScheme = "bfactor")
              ) %>%
              stageParameters(backgroundColor = "black") %>%
              setQuality("high") %>%
              setFocus(0) %>%
              setSpin(FALSE)
          })
        }}}})
  
  ### Resetting the structure in "Discover"
  observeEvent(input$reset_discover, {
    parentfolder <- "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/"
    #parentfolder <- "/srv/shiny-server/proteorizer/"
    if (nchar(input$proteincase)>=2){
      proteinID_raw <- strsplit(input$proteincase,"-")
      proteinID <- proteinID_raw[[1]][1]
      alphafoldfile_discover <- paste("AF-",proteinID,"-F1-model_v3.pdb",sep="")  
      structurepath_discover <- file.path(parentfolder, "www",input$dataset,input$proteincase,alphafoldfile_discover)
      output$structure_discover <- renderNGLVieweR({
        NGLVieweR(structurepath_discover) %>%
          addRepresentation("cartoon",
                            param = list(name = "cartoon",
                                         colorScheme = "bfactor")
          ) %>%
          stageParameters(backgroundColor = "black") %>%
          setQuality("high") %>%
          setFocus(0) %>%
          setSpin(FALSE)
      })}})
  ### changing structure background in "Discover"
  observeEvent(input$color_discover, {
    if (nchar(input$proteincase)>=2){
      NGLVieweR_proxy("structure_discover") %>%
        updateStage(param = list('backgroundColor' = isolate(input$color_discover))
        )}
    
  })
  ### taking a Snapshot of the structure in "Discover"
  observeEvent(input$snapshot_discover, {
    if (nchar(input$proteincase)>=2){
      NGLVieweR_proxy("structure_discover") %>%
        snapShot("Snapshot.png", param = list(antialias = TRUE,
                                              trim = TRUE,
                                              transparent = TRUE,
                                              scale = 1))
    }})
  ### Adding residue highlights to the structure in "Discover"
  observeEvent(input$add_discover, {
    parentfolder <- "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/R_Shiny_Interface/"
    #parentfolder <- "/srv/shiny-server/proteorizer/"
    if (nchar(input$proteincase)>=2){
      colorpicker = 1
      colorlist <- c("red","blue","green","brown","purple","yellow","orange")
      ### Getting the information that needs to be highlighted on the structure
      tryCatch({
        if (input$dataset=="VUS_Humsavar_RW"){ #if (input$dataset=="Proteorizer_RW_VUS_Humsavar"){
          tablepath <- file.path(parentfolder,"www",input$dataset,input$proteincase,"FinalResults_unlim_conservation_scored.txt")}
        else {tablepath <- file.path(parentfolder, "www",input$dataset,input$proteincase,"FinalResults_unlim_conservation_scored_hclust.txt")}
        partfilus <- read.delim(tablepath, header=TRUE, quote="")
        inputclusters <- c()
        inputpositions <- c()
        positions_to_highlight <- c()
        for (row in 1:nrow(partfilus))
        {
          clusterinfo <- partfilus[row,"Clusternumber"]
          clusterinfo <- gsub(" ", "", clusterinfo)
          positional <- partfilus[row,"Position"]
          positional <- gsub(" ","",positional)
          sourcus <- partfilus[row,"Data_Source"]
          sourcus <- gsub(" ","",sourcus)
          ### If I know a certain cluster number corresponds to a clustered input position
          if (clusterinfo %in% inputclusters){
            newinfo <- paste(clusterinfo,positional,sep="/") ### X/Y, where X is the clusternumber and Y the residue number
            positions_to_highlight <- c(positions_to_highlight,newinfo)
          }
          else {}
          
          if (sourcus == "Input"){
            if (clusterinfo != "NoClusterMembership"){
              if (nchar(clusterinfo)>=1){
                #print(clusterinfo)
                goodtoknow <- paste(clusterinfo,positional,sep="/") ### X/Y, where X is the clusternumber and Y the residue number
                inputpositions <- c(inputpositions, goodtoknow)
                inputclusters <- c(inputclusters, clusterinfo)}
              else{}
            }
          }
          else {}
        }
      },
      error = function(e) {
        return()})
      
      for (mutation in inputpositions){ ### "3/127" "2/97" 
        realmut <- strsplit(mutation,"/")[[1]][2]
        realclust <- strsplit(mutation,"/")[[1]][1]
        NGLVieweR_proxy("structure_discover") %>%
          addSelection("ball+stick",
                       param =
                         list(
                           name = realmut,
                           sele = isolate(as.character(gsub('\\D+','', realmut))),
                           colorValue = isolate(nth(colorlist,colorpicker))
                         ))
        for (position in positions_to_highlight){ ### "3/128" "2/98"  "3/127" "2/98"  "3/128"
          posclust <- strsplit(position,"/")[[1]][1] 
          pospos <- strsplit(position,"/")[[1]][2] 
          if (posclust == realclust){
            NGLVieweR_proxy("structure_discover") %>%
              addSelection("ball+stick",
                           param =
                             list(
                               name = mutation,
                               sele = isolate(as.character(gsub('\\D+','', pospos))),
                               colorValue = isolate(nth(colorlist,colorpicker))
                             )) %>%
              addSelection("label",
                           param = list(
                             sele = realmut,
                             labelType = "format",
                             labelFormat = "[%(resname)s]%(resno)s", # or enter custom text
                             labelGrouping = "residue", # or "atom" (eg. sele = "20:A.CB")
                             color = isolate(nth(colorlist,colorpicker)),
                             fontFamiliy = "sans-serif",
                             xOffset = 1,
                             yOffset = 0,
                             zOffset = 0,
                             fixedSize = TRUE,
                             radiusType = 1,
                             radiusSize = 1.5, # Label size
                             showBackground = FALSE
                           ))
          }else {}
          
        }
        colorpicker = colorpicker + 1
      }
      
    }})
  
}

