library(shiny)
library(Biobase)
library(data.table)
library(rjson)
library(DT)
library(shinyBS)
library(shinythemes)
library(ggplot2)
source("ggheat.continuous.R")

##make sure to load dat beforehand
##load app using source_HEPG2.R or source_MCF10A.R

uniq_id<-"dose (uM)"

if(dat$title == "MCF10A Portal"){
	uniq_id <- "unique_ID_by_chem"
}

domain<-paste0("https://carcinogenome.org/data/", sub(" .*", "", dat$title))

defaults<-list(landmark_de = FALSE, 
		summarizefunc_de = "median", 
		filterbyinput_de = c("score", "number"),
		range_de = c(-2, 2), 
		numberthresleft_de = 10, 
		numberthresright_de = 10)

tas<-dat[["Profile Annotation"]]$TAS
minProf<-3
maxTAS<-tas[order(tas, decreasing = TRUE)][minProf]
maxTAS<- floor(maxTAS/0.1)*0.1

to.hex<-function(x){
  cols<-col2rgb(x)
  red<-cols[1]
  green<-cols[2]
  blue<-cols[3]
  return(rgb(red, green, blue, maxColorValue = 255))
}


data.table.round<-function(dt, digits = 3){
	cols<-sapply(colnames(dt), function(i) is.numeric(dt[,i]))
	cols<-names(which(cols))
	for(i in cols)
		dt[,i]<-round(dt[,i], digits)
	dt<-data.table(dt)
}

ui = bootstrapPage(
  numericInput('n', 'Number of bins', 100),
  plotOutput('plot')
)

server = function(input, output) {
  output$plot <- renderPlot({ hist(dat, breaks = input$n) })
}

customTable<-function(i, ...){
	DT::renderDataTable(data.table.round(i), 
		escape = FALSE,
		extensions = 'Buttons', 
		server = TRUE,
		options = list(dom = 'T<"clear">Blfrtip', 
			deferRender = FALSE,
			scrollX = TRUE, 
			scrollY = 400,
			scrollCollapse = TRUE,
			pageLength =15, 
			lengthMenu = c(15,50,100,200,1000),
			buttons=c('copy','csv','print')))
}

de_opts<-function(){

	bsCollapse(id = "de_opt_panel", 
	open = "Options",
	bsCollapsePanel("Options", 
	fluidRow(
	column(2, 


		checkboxInput("landmark_de", "Landmark only", value =defaults[["landmark_de"]])),
	column(2, selectInput("summarizefunc_de", "Summarization:",
		sumrules,
	  selected = defaults[["summarizefunc_de"]])),

	column(2,checkboxGroupInput("filterbyinput_de", "Filter by:",
	             c("score" = "score",
	               "number" = "number"),
	             selected = defaults[["filterbyinput_de"]])),
	column(2,sliderInput("range_de", "score threshold", min = -10, max = 10, 
	  value = defaults[["range_de"]], step = 0.01)),
	column(2,sliderInput("numberthresleft_de", "Num +",
	  min = 0, max = 1000, value = defaults[["numberthresleft_de"]], ticks = FALSE, step = 10)),
	column(2,sliderInput("numberthresright_de", "Num -",
	  min = 0, max = 1000, value = defaults[["numberthresleft_de"]], ticks = FALSE, step = 10)) 
	
	),
actionButton("restore", "Restore Defaults", style=mybutton)
	))
}

gs_opts<-function(){
	fluidRow(column(3,
	selectInputWithTooltip(inputId = "gsname",
	            	  				label = "Gene set name", 
	            	  				choices = names(dsmap), 
	            	  				bId = "Bgsname",  
	            	  				helptext = helptextgsname)),
	column(3,selectInputWithTooltip(inputId = "gsmethod", 
		label = "Projection method", 
		choices = gsmethods, 
		bId = "Bgsmethod", 
		helptext =helptextgsmethod)),
	column(3,selectInput("summarize_gs", 
		"Sort by:", 
		sumrules,
		selected = "median")))
}

marker_gs_opts<-function(inputId = "marker_gsname",
	bId = "Bgsname_marker",
	inputId2 = "marker_gsmethod"){
	fluidRow(column(3,
	selectInputWithTooltip(inputId = inputId,
	            	  				label = "Gene set name", 
	            	  				choices = names(dsmap), 
	            	  				bId = bId,  
	            	  				helptext = helptextgsname)),
	column(3,selectInputWithTooltip(inputId = inputId2, 
		label = "Projection method", 
		choices = gsmethods, 
		bId = "Bgsmethod_marker", 
		helptext =helptextgsmethod))
	)
}

conn_opts<-function(){
	fluidRow(
		column(3,
			selectInput(inputId = "conn_name",
				label = "Connectivity Level", 
				choices = names(connmap))
			),
		column(8, 
			selectInput("summarizefunc_conn", 
			"Sort by:", 
			sumrules,
			selected = "median")
		)
	)
}

marker_conn_opts<-function(inputId = "marker_conn_name", hm = FALSE){
	choices = names(connmap)
	if(hm) choices = "Perturbagen Classes"
	fluidRow(
		column(3,
		selectInput(inputId = inputId,
			label = "Connectivity Level", 
			choices = choices)
		)
	)
}

get_ids_pdat<-function(pdat, 
  cols = c("Chemical Name", "CAS", "BUID"), 
  col.unique = "BUID",
  val.ignore = c("", " ", NA, "NA", "NOCAS")){
  tab<-unique(pdat[, cols])
  res<-lapply(cols, function(i){
  x<-as.character(tab[,i])
  x.uniq<- setdiff(unique(x), union(x[duplicated(x)], val.ignore))
  x.uniq<-sort(x.uniq)
  })
  names(res)<-cols
  return(res)
}

get_BUID<-function(input, tab){
  as.character(tab[which(apply(tab, 1, function(i) any(i %in% input)))[1], "BUID"])
}

summarize_eset<-function(mat, 
  summarize.func = c("mean", "median", "max", "min", "Q1", "Q3"),
  do.scorecutoff = TRUE, scorecutoff = c(-0.6, 0.6), 
  do.nmarkers = TRUE, nmarkers = c(100, 100)
  ){

  summarize.func<- match.arg(summarize.func)
  x<-apply(mat, 1, match.fun(summarize.func))
  x<-as.numeric(x)
  n<-length(x)
  
  if(do.nmarkers){
  	
  	ind0<-sum(x > 0)
  	n1<-min(nmarkers[1], ind0)
  	n2<-min(nmarkers[2], n-ind0)
    ord<-order(x, decreasing = TRUE)
    n2ind<-n-n2+1
    if(n1 == 0 & n2 == 0) x.ind.nmarkers<- NULL
    else if(n2 == 0) x.ind.nmarkers<-ord[1:n1]
	else  x.ind.nmarkers<-c(ord[1:n1], ord[n2ind:n])
  } else
    x.ind.nmarkers<-1:n

  if(do.scorecutoff)
    #TODO: rank by score here too
    x.ind.scorecutoff<-which(x > scorecutoff[2] | x < scorecutoff[1])
  else
    x.ind.scorecutoff<-1:n

  inds<-intersect(x.ind.nmarkers, x.ind.scorecutoff)
  inds<-inds[order(x[inds], decreasing = TRUE)]
  return(list(inds = inds, scores = x[inds]))
} 

get_connectivity<-function(input, connlist, tab, 
	conn_name,
	annot_prof,
	match_id = "sig_id",
	col_id = uniq_id,
	header = "Connectivity Score",
	summarize.func = c("mean", "median", "max", "min", "Q1", "Q3")){
	 i<-get_BUID(input, tab)
	 currname<-connmap[[conn_name]]
	 res<-connlist[[currname]]
	 ind2<-match(pData(res)[, match_id], annot_prof[, match_id])
	 pData(res)<-annot_prof[ind2,]
	 res<-res[, res$BUID %in% i]
	 res<-summarize_gsproj(eset = res, 
	 	annot_prof, match_id, col_id, header, summarize.func)
}

get_gsproj<-function(input, gslist, tab, gsname, gsmethod,
	annot_prof, 
	match_id = "sig_id", 
	col_id = uniq_id, 
	header = "GS Score",
	summarize.func = c("mean", "median", "max", "min", "Q1", "Q3")
	){
  i<-get_BUID(input, tab)
  ind<-grep(paste0(".*", dsmap[[gsname]], ".*", gsmethod, ".*"), names(gslist))
  res<-gslist[[ind]]
 
  ind2<-match(pData(res)[,match_id], annot_prof[, match_id])
  pData(res)<-annot_prof[ind2,]

  res<-res[, res$BUID %in% i]

  res<-summarize_gsproj(eset = res, 
  	annot_prof, match_id, col_id, header, summarize.func)

  #return hyperlink to MSigDB genesets
  if(gsname %in% c("Hallmark", "C2"))
    res$genesets<-sapply(as.character(res$genesets), get_geneset_link)

  return(res)
}

get_de_eset<-function(annot_prof, match_id = "sig_id", hm = FALSE){
  res<-dat[["Gene Expression"]]
  ind2<-match(pData(res)[,match_id], annot_prof[, match_id])
  pData(res)<-annot_prof[ind2,]
  if(hm){
  	inds<-fData(res)[, "Landmark Gene"] %in% "Yes"
  	res<-res[inds,]
  }
  return(res)
}

get_gs_eset<-function(gslist, gsname, gsmethod, annot_prof,
	match_id = "sig_id"){
  ind<-grep(paste0(".*", dsmap[[gsname]], ".*", gsmethod, ".*"), names(gslist))
  res<-gslist[[ind]]
  ind2<-match(pData(res)[,match_id], annot_prof[, match_id])
  pData(res)<-annot_prof[ind2,]
  return(res)
}

get_conn_eset<-function(connlist, 
 				conn_name, 
 				annot_prof,
				match_id = "sig_id"){
  ind<-which(names(connmap) %in% conn_name)
  res<-connlist[[ind]]
  ind2<-match(pData(res)[,match_id], annot_prof[, match_id])
  pData(res)<-annot_prof[ind2,]

  return(res)
}

summarize_gsproj<-function(eset, 
	annot_prof, 
	match_id = "sig_id", 
	col_id = uniq_id, 
	header = "GS Score",
	summarize.func = c("mean", "median", "max", "min", "Q1", "Q3")){

	head(pData(eset))
	summarize.func<- match.arg(summarize.func)
  	res<-apply(exprs(eset), 1, match.fun(summarize.func))
  	res<-as.numeric(res)

	mat<-exprs(eset)	
	colnames(mat)<-paste(header, " ", pData(eset)[,col_id], "uM",sep = "")
		
	res<-cbind(fData(eset), score = res, mat)
	res<-res[order(res$score, decreasing = TRUE),, drop = FALSE]
	colnames(res)[colnames(res) %in% "score"]<- "Summary Score"
	return(res)
}

get_genecard_link<-function(genesymbol){
  sprintf('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s&keywords=%s" target="_blank" class="btn btn-primary">%s</a>',
    genesymbol, genesymbol, genesymbol)
}

get_de<-function(input, tab, 
  eset, annot_prof, 
  match_id = "sig_id", 
  col_id = uniq_id, 
  header = "ModZScore", 
  landmark = FALSE, do.scorecutoff = TRUE, scorecutoff = c(-2, 2), 
  do.nmarkers = TRUE, nmarkers = c(100, 100),
  summarize.func = c("mean", "median", "max", "min", "Q1", "Q3"), 
  landmark_id = "Landmark Gene", 
  landmark_positive = "Yes"){
  i<-get_BUID(input, tab)
  pData(eset)<-annot_prof[match(colnames(eset), annot_prof[, match_id]),]
  eset<-eset[, eset$BUID %in% i]
  
  if(landmark)
    eset<-eset[fData(eset)[, landmark_id] %in% landmark_positive,]
  
  mat<-exprs(eset)
  fdat<-fData(eset)
  pdat<-pData(eset)

  colnames(mat)<-paste(header, " ", pdat[, col_id], "uM", sep = "")
  res<-summarize_eset(mat, summarize.func, do.scorecutoff, scorecutoff,
    do.nmarkers, nmarkers)

  res.ind<-res$inds
  res.scores<-res$scores
  direction<-sapply(res.scores, function(i){
  	if(i > 0) return("Up")
  	else return("Down")
  	})
  tab<-cbind(fdat[res.ind,, drop = FALSE], Direction = direction,
  	SummaryScore=res.scores, mat[res.ind,, drop = FALSE])
  colnames(tab)[colnames(tab) %in% "SummaryScore"]<-"Summary Score"

  tab$"Gene Symbol"<-sapply(as.character(tab$"Gene Symbol"), get_genecard_link)

  return(tab)
}

get_chemical_description<-function(i){
	buid<-get_BUID(input = i, tab = annot_chem)
	res<-annot_chem[annot_chem$BUID %in% buid, ]
	return(res[1,])
}

add_padding<-function(){
	##insert empty space
	fluidRow(column(width = 1, 
		offset = 0, 
		style='padding:10px;')
	)
}

selectInputWithTooltip<-function(inputId, label, choices, bId, helptext, ...){
  selectInput(inputId, tags$span(label,  
            tipify(bsButton(bId, "?", style = "inverse", size = "extra-small"), 
              helptext)),
             choices, ...)}

get_geneset_link<-function(geneset){
  sprintf('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/%s" target="_blank" class="btn btn-primary">%s</a>',
    geneset, geneset)

}

get_morpheus_link<-function(url = "toy.gct", 
	domain = "https://carcinogenome.org/data", tas = 0){

	url<-sprintf('{
	"dataset" : "%s/%s",
	"columns" : [ { "field" : "id", "display" : ["color"] },
				  { "field" : "Chemical Name", "display" : ["color"] },
				  { "field" : "CAS", "display" : ["color"] },
				  { "field" : "dose (uM)", "display" : ["color"] },
				  { "field" : "Carcinogenicity", "display" : ["color"] },
				  { "field" : "Genotoxicity", "display" : ["color"] },
			      { "field" : "TAS", "display" : ["color"] } ],
	"rows" : [ { "field" : "genesets", "display" : ["text"] }],
	"columnFilter" : {
		"filters" : [{
		"field" : "TAS",
		"type" : "range",
		"min" : %f,
		"max" : 1.0
		}]
		}
	}', domain, url, tas)

  url<-URLencode(URL = url)
  url<-paste("https://software.broadinstitute.org/morpheus/?json=", url, sep = "")
  return(url)
}


get_heatmap_gct<-function(ds, dsmap, method, domain){
  dsname<-dsmap[[ds]]
  url<-paste(dsname, "_", method, ".gct", sep = "")
}

dsmap <- list(Hallmark="gsscores_h.all.v5.0",
    C2="gsscores_c2.cp.reactome.v5.0", 
    NURSA="gsscores_nursa_consensome_Cbyfdrvalue_0.01.gmt")

connmap<-list(PCL = "pcl", PERT = "pert")
names(connmap)<-c("Perturbagen Classes", "Perturbagens")

gsmethods<-c("gsva", "ssgsea", "zscore", "gsproj")

#tooltip texts

helptextgsname<-HTML(paste("Hallmark: MSigDB Hallmark Pathways (v5.0)",
              "C2: MSigDB C2 reactome Pathways (v5.0)",
              "NURSA: Nuclear Receptor Signaling Atlas, consensome data for human", 
               sep="<br/>"))

helptextgsmethod<-HTML(paste(
              "gsva, ssgea, zscore: from R Bioconductor package GSVA",
              "gsproj: GeneSetProjection for R package montilab:CBMRtools", 
               sep="<br/>"))

chemicals<-get_ids_pdat(pdat = dat[["Chemical Annotation"]])
annot_chem<-dat[["Chemical Annotation"]]
annot_prof<-dat[["Profile Annotation"]]

mybutton<-"color: #98978b; background-color: #F8F5F0; padding: 7px 3px; margin: 1px; text-align:center; 
display:inline-block; font-size:10px"
sumrules<-c("max", "median", "mean", "min", "Q1", "Q3")
markers<-c("Genes", "Gene Sets", "CMap Connectivity")
markers_hm<-c("Genes (Landmark)", "Gene Sets", "CMap Connectivity")

Q1<-function(x){ 
	quantile(x, 0.25, na.rm = T)
}

Q3<-function(x){
	quantile(x, 0.75, na.rm = T)
}

get_de_by_gene_table<-function(input, eset, annot_prof, 
	match_id = "sig_id",  
	header = "mod Z-scores",
	tas = 0){

	if(input == "" | is.null(input) | is.na(input) | !(input %in% rownames(eset))) return(NULL)

	pData(eset)<-annot_prof[match(colnames(eset), annot_prof[, match_id]),]

	eset<-eset[, eset$TAS >= tas]
	

    rowid<-which(rownames(eset) %in% input)[1]
    x<-as.numeric(exprs(eset)[rowid,])
    ord<-order(x, decreasing = TRUE)
    pdat<-pData(eset)

    df<-cbind(value = x, pdat)
    df<-df[ord,]

    colnames(df)[colnames(df) %in% "value"]<-header
    return(df)
}

get_de_by_gene_hist<-function(input, eset, annot_prof, 
	match_id = "sig_id", 
	col_id = NA, 
	col_colors = c("grey", "green", "orange", "grey"), 
	col_names = c("N/A", "-", "+"), 
	header = "mod Z-scores",
	tas = 0, 
	plot = "Density"){

	if(input == "" | is.null(input) | is.na(input) | !(input %in% rownames(eset))) return(NULL)

	pData(eset)<-annot_prof[match(colnames(eset), annot_prof[, match_id]),]

	eset<-eset[, eset$TAS >= tas]
	

    rowid<-which(rownames(eset) %in% input)[1]
    x<-as.numeric(exprs(eset)[rowid,])

    plot_wrapper<-function( plot = "Density", ...){
    	if(plot %in% "Density") {
    		res<-ggplot(df, aes_string(x = "x", fill = "cols"))+ geom_density(position = "identity", alpha = 0.5, ...)
    	} else {
    		
    		res<-ggplot(df, aes_string(x = "cols", y= "x", fill = "cols"))+ 
    		geom_boxplot(
    		  position = "identity", 
    			width = 0.2,
    			alpha = 0.5, 
    			outlier.fill = NULL,
    			outlier.alpha = NULL,
    			...)
    	
    	}
    	return(res)
    }
    
    if(is.na(col_id)){
    	p.title<-paste("Distribution of ", header, " across profiles for ", input, " (Overall)", sep = "")
    	background<-as.numeric(exprs(eset))
	    df<-rbind(data.frame(x = x, cols = "query"),
	    	data.frame(x = background, cols = "background"))

	    df$cols <- factor(df$cols, levels = c("background", "query"))

	    col_vec<-c("grey", "red")
	    names(col_vec)<-c("background", "query")
	    p<-plot_wrapper(plot)+
	    xlab(header) + 
	    ylab("Count")+ 
	     scale_fill_manual(name = "Overall", 
	     	values = col_vec,
	     	breaks = names(col_vec),
	     	labels = names(col_vec)
	     	)+
	    ggtitle(p.title)+
	    theme_bw()+
	    theme(plot.title = element_text(hjust = 0.5))
	}
    if(!is.na(col_id)){
    	p.title<-paste("Distribution of ", header, " across profiles for ", input, " (by ", col_id, ")", sep = "")
    	cols<-pData(eset)[, col_id]
	    cols<-as.character(cols)
	    cols_match<-col_colors
	    names(cols_match)<-col_names
	    df<-data.frame(x = x, cols= cols)

	    df$cols <- factor(df$cols, levels = col_names)

	    p<-plot_wrapper(plot)+
	    scale_fill_manual(
	    	name = col_id,
	    	values = cols_match, 
	    	breaks = names(cols_match), 
	    	labels = names(cols_match))+
	    xlab(header) + 
	    ylab("Density")+ 
	    ggtitle(p.title)+
	    theme_bw()+
	    theme(plot.title = element_text(hjust = 0.5))
	 	
    }
    return(p)
}

app<-shinyApp(

	ui = shinyUI(
		
		fluidPage(
			theme = shinytheme("sandstone"),
			#override css here
			#turn off slider highlighting
			#adjust button container width
  			tags$style(HTML(".irs-bar {background: none; border-top: none; border-bottom: none; border-left:none;}
   			.irs-bar-edge {background: none; border-top: none; border-bottom: none; border-left:none;}
  			.col-sm-1 {width: 40px; margin: 0px; padding: 0px 0px;}")),
	
			#tags for the entire page
			tags$head(includeScript("google-analytics.js")),
			#added padding to to navbarPage position = fixed-top
			tags$style(type="text/css", "body {padding-top: 70px;}"),

			navbarPage(dat[["title"]], position = c("fixed-top"),

				tabPanel("About",
					titlePanel(dat[["title"]]),
					      fluidRow(
					        column(11,
					          includeMarkdown(dat[["about page"]])
					        ),
					        img(src="../logo.png", align = "left", width = 600)
					      )
					),

				tabPanel("Annotation",
					radioButtons("annot_table_selection", 
						label = "Annotation Type",
						choices = c("Chemicals", "Samples"), 
						selected = "Chemicals", 
						inline = TRUE),
					conditionalPanel(
						condition = "input.annot_table_selection == 'Chemicals'",
						DT::dataTableOutput("t1") 
						),
					conditionalPanel(
						condition = "input.annot_table_selection == 'Samples'",
						DT::dataTableOutput("t2") 
						)
					),

				tabPanel("Chemical Explorer",
					fluidRow(
				    column(3,selectizeInput("chem", 
				    	"Select chemical:",   
				     	options = list(
    						placeholder = 'Please select an option below',
    						onInitialize = I('function() { this.setValue(""); }')
							),
				     	choices = chemicals,
				     	selected = ""
				     	)),
				    column(9,
				    	conditionalPanel("input.chem !=''", 
				    		DT::dataTableOutput("t3")))
				    ),
				    add_padding(),
	            	conditionalPanel("input.chem != ''",
	            		tabsetPanel(
	            	  		tabPanel("Gene Expression", 
	            	  			fluidRow(
	            	  			column(11,
	            	  				de_opts()),
	            	  			column(1,
	            	  			actionButton("de_hide", "Hide",
    								style=mybutton)),
	            	  			actionButton("de_show", "Show",
    								style=mybutton)),	                        	  			
	            	  			DT::dataTableOutput("t4") 
	            	  			),
	            	  		tabPanel("Gene Set Enrichment",
	            	  			column(10,
	            	  				gs_opts()),
	            	  			DT::dataTableOutput("t5")
	            	  			),
	            	  		tabPanel("Connectivity", 
	            	  			column(10,
	            	  				conn_opts()),
	            	  			DT::dataTableOutput("t6")
	            	  			)
	            	  		)	            		
	            		)     
				),	
				tabPanel("Marker Explorer", 
 					
					fluidRow(
 					column(3,selectizeInput("marker", 
				    	"Select marker set:",   
				     	options = list(
    						placeholder = 'Please select an option below',
    						onInitialize = I('function() { this.setValue(""); }')
							),
				     	choices = markers,
				     	selected = "Genes"
				     	)),

 					column(3,sliderInput("marker_tas", "TAS range", min = 0,  max = maxTAS, 
	  				value = 0.2, step = 0.1)),

	  				column(3, 
	  					radioButtons("marker_view", "Plot Type", 
	  						choices = c("Density", "Boxplot"), selected = "Density"))

 					),
				    conditionalPanel("input.marker == 'Genes'",
				    		 					
				    selectizeInput("marker_gene", 
				    	"Select a gene:",   
				     	options = list(
    						placeholder = 'Please select an option below',
    						onInitialize = I('function() { this.setValue(""); }')
							),
				     	choices = "",
				     	selected = ""
				     	
				     	),

				    conditionalPanel("input.marker_gene != false && input.marker_gene != ''",
				    	conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            tags$div("Loading...",id="loadmessage")),	
						conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",						
							downloadLink("t7_download_png", "Download png"),
							downloadLink("t7_download_pdf", "pdf"),
					     	plotOutput("t7_1"),
					     	plotOutput("t7_2"),
					     	plotOutput("t7_3"),
					     	column(8, align = "left", h2(textOutput("t7_text"))),
					     	DT::dataTableOutput("t7_table")					     
				     		)
				    	)

				    ),
				    conditionalPanel("input.marker == 'Gene Sets'",

					marker_gs_opts(),
					selectizeInput("marker_gs",
						"Select a gene set:",
 				     	options = list(
    						placeholder = 'Please select an option below',
    						onInitialize = I('function() { this.setValue(""); }')
							),
				     	choices = "",
				     	selected = ""
						),
					
					conditionalPanel("input.marker_gs != ''",
				    	conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            tags$div("Loading...",id="loadmessage")),	
						conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
							downloadLink("t8_download_pdf", "Download pdf"),
							downloadLink("t8_download_png", "png"),
					     	plotOutput("t8_1"),
					     	plotOutput("t8_2"),
					     	plotOutput("t8_3"),
					     	column(8, align = "left", h2(textOutput("t8_text"))),
						    DT::dataTableOutput("t8_table")
						    )
				     	)						

				    ),

				    conditionalPanel("input.marker == 'CMap Connectivity'",
				    marker_conn_opts(),

				    selectizeInput("marker_conn",
						"Select a CMap hit:",
 				     	options = list(
    						placeholder = 'Please select an option below',
    						onInitialize = I('function() { this.setValue(""); }')
							),
				     	choices = "",
				     	selected = ""
						),

				    conditionalPanel("input.marker_conn != ''",
				    	conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            tags$div("Loading...",id="loadmessage")),	
						conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
							downloadLink("t9_download_pdf", "Download pdf"),
							downloadLink("t9_download_png", "png"),
					     	plotOutput("t9_1"),
					     	plotOutput("t9_2"),
					     	plotOutput("t9_3"),
					     	column(8, align = "left", h2(textOutput("t9_text"))),
						    DT::dataTableOutput("t9_table")
						    )				     	
				     	)				    
				    )		
				),


				tabPanel("Heatmap Explorer", 
 					
					fluidRow(
 					column(3,selectizeInput("marker_hm", 
				    	"Select marker set:",   
				     	options = list(
    						placeholder = 'Please select an option below',
    						onInitialize = I('function() { this.setValue(""); }')
							),
				     	choices = markers_hm,
				     	selected = "Genes (Landmark)"
				     	)),

 					column(3,sliderInput("marker_tas_hm", "TAS range", min = 0,  max = maxTAS, 
	  				value = 0.2, step = 0.01)),
	  				##warning
	  				textOutput("warning")
 					),
				    conditionalPanel("input.marker_hm == 'Genes (Landmark)'",	
				    	conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            tags$div("Loading...",id="loadmessage")),	
				     	conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
							uiOutput("t10_combined"))
				    ),
				    conditionalPanel("input.marker_hm == 'Gene Sets'",
						marker_gs_opts(inputId = "marker_gsname_hm",
							bId = "Bgsname_marker_hm",
							inputId2 = "marker_gsmethod_hm"),	
						
						conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            tags$div("Loading...",id="loadmessage")),	
						conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
							uiOutput("morpheus_result_link"),
							uiOutput("t11_combined"))	
						
				    ),

				    conditionalPanel("input.marker_hm == 'CMap Connectivity'",
				    	marker_conn_opts(inputId = "marker_conn_name_hm", hm = TRUE),	
				    	conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            tags$div("Loading...",id="loadmessage")),
						conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
							uiOutput("t12_combined"))
				    )		
				),

   			tabPanel(HTML("</a></li><li><a href=\"../\">Home"))

			)
		)
	),

	server = shinyServer(function(input, output, session){
		get_gs_call<-function(){
			get_gsproj(input = input$chem, 
				gslist = dat[["Gene Set Enrichment"]], 
				tab = annot_chem[, c("Chemical Name", "BUID", "CAS")], 
				gsname = input$gsname, 
				gsmethod = input$gsmethod,
				annot_prof = dat[["Profile Annotation"]], 
				match_id = "sig_id", 
				col_id = uniq_id, 
				header = "GS Score",
				summarize.func = input$summarize_gs)
		}
		get_de_call<-function(){
		    get_de(input$chem, 
				tab=annot_chem[, c("Chemical Name", "BUID", "CAS")], 
		        eset = dat[["Gene Expression"]], 
		        annot_prof = annot_prof,
		        match_id = "sig_id",
		        col_id = uniq_id,
		        landmark = input$landmark_de, 
		        do.scorecutoff = "score" %in% input$filterbyinput_de, 
		        scorecutoff = c(input$range_de[1], input$range_de[2]), 
		        do.nmarkers = "number" %in% input$filterbyinput_de, 
		        nmarkers = c(input$numberthresleft_de, input$numberthresright_de),
		        summarize.func = input$summarizefunc_de)

		}

		get_conn_call<-function(){
			get_connectivity(input$chem, 
				connlist = dat[["Connectivity"]], 
				tab = annot_chem[, c("Chemical Name", "BUID", "CAS")], 
				conn_name = input$conn_name,
				annot_prof = annot_prof,
				match_id = "sig_id",
				col_id = uniq_id,
				header = "Connectivity Score",
				summarize.func = input$summarizefunc_conn)
		}

		plot_heatmap_static<-function(plot_name, eset, tas){

			col_legend<-list(Carcinogenicity = list(col_breaks = c("+", "-", "N/A"), 
			col_values = sapply(c("orange", "green", "grey"), to.hex),
			col_labels = c("+", "-", "N/A")),
			Genotoxicity = list(col_breaks = c("+", "-", "N/A"), 
			col_values = sapply(c("purple", "pink", "grey"), to.hex),
			col_labels = c("+", "-", "N/A")))

			hmcolors<-function(...) scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, ...)

			eset<-eset[, pData(eset)[, "TAS"]>= tas]
			ind.remove<-apply(exprs(eset), 1, function(i){
				any(is.nan(i))
				})

			eset<-eset[!ind.remove, ]

			hc<-clust_eset(eset)

			ncols<-ncol(eset)
			nrows<-nrow(eset)

			h<-min(nrows*15+200, 3000)
			w<-max(min(ncols*15+50, 3000), 400)
			
			xsize<-4
			if (ncols<500) 
			xsize<-5
			if (ncols<100)
			xsize<-6
			if(ncols<10)
			xsize<-10
			if(ncols > 1000)
			xsize<-0 
			
			ysize<-4
			if (nrows<500) 
			ysize<-5
			if (nrows<100)
			ysize<-6
			if(nrows<10)
			ysize<-10
			if(nrows > 1000)
			ysize<-0 

			col_lab<-c("Carcinogenicity", "Genotoxicity")

			ph2<-3*length(col_lab)
			ph3<-nrows+20
			ph1<-0.1*(ph2+ph3)

			p.heights<-c(ph1, ph2, ph3)

			#render plot
			p<-ggheat.continuous.single(eset = eset, 
			hc = hc$hc, 
			hr = hc$hr, 
			hmcolors = hmcolors,
			hmtitle = "",
			col_lab = col_lab, 
			col_legend = col_legend,
			ylabstr = "",
			fout = NA, 
			p.heights = p.heights,
			xsize = xsize,
			ysize = ysize, 
			ysizelab = 7)

			output[[plot_name]]<-renderPlot({
				do.call(grid.arrange, p)
			}, width = w, height = h
			)

			#render download link
			output[[paste0(plot_name, "_download_pdf")]]<-downloadHandler(
	      		filename = paste0("carcinogenome_download_", "heatmap", ".pdf"),
	      		content = function(file){ 			
	      			ggsave(plot = do.call(grid.arrange, p), 
	      				filename = file, 
	      				device = "pdf", 
	      				width = w/61, height = h/61, 
	      				dpi = 300
	      				)
	      		})
			output[[paste0(plot_name, "_download_png")]]<-downloadHandler(
	      		filename = paste0("carcinogenome_download_", "heatmap", ".png"),
	      		content = function(file){ 			
	      			ggsave(plot = do.call(grid.arrange, p), 
	      				filename = file, 
	      				device = "png", 
	      				width = w/61, height = h/61, 
	      				dpi = 300
	      				)
	      		})
			#render combined
			output[[paste0(plot_name, "_combined")]] <- renderUI({	
				fluidRow(			
				downloadLink(paste0(plot_name, "_download_pdf"), "Download pdf"),
				downloadLink(paste0(plot_name, "_download_png"), "png"),
	      		plotOutput(plot_name, height = h, width = w)
	      		)
	    	})
		}

		plot_hist<-function(plot_name, eset, header, markerid, plot = "Density"){

			
			p1<-get_de_by_gene_hist(markerid, eset,
	      			annot_prof = dat[["Profile Annotation"]], 
	      			match_id = "sig_id",
	      			header = header,
	      			tas = input$marker_tas, 
	      			plot = plot)

			p2<-get_de_by_gene_hist(markerid, eset, 
				annot_prof = dat[["Profile Annotation"]], 
				match_id = "sig_id", 
				col_id = "Carcinogenicity", 
				col_colors = c("grey","green", "orange"), 
				col_names = c("N/A","-", "+"),
				header = header,
				tas = input$marker_tas, 
				plot = plot)
			p3<-get_de_by_gene_hist(markerid, eset, 
				annot_prof = dat[["Profile Annotation"]], 
				match_id = "sig_id", 
				col_id = "Genotoxicity", 
				col_colors = c("grey","pink", "purple"), 
				col_names = c("N/A", "-", "+"),
				header = header,
				tas = input$marker_tas,
				plot = plot)

			w<-1000
			h<-300

			if(plot == "Density"){
		
				output[[paste0(plot_name, "_", 1)]]<-renderPlot(p1,
		      		width = w, height = h
		      		)

		      	output[[paste0(plot_name, "_", 2)]]<-renderPlot(p2,
					width = w, height = h
		      		)

		      	output[[paste0(plot_name, "_", 3)]]<-renderPlot(p3,
					width = w, height = h
		      		)

		      	output[[paste0(plot_name, "_download_pdf")]]<-downloadHandler(
	      		filename = paste0("carcinogenome_download_", header, "_", markerid, ".pdf"),
	      		content = function(file){
	      			ggsave(file, device = "pdf", grid.arrange(p1, p2, p3),
	      				width = w/130, height = h/35, 
	      				dpi = 300)
	      		})

		      	output[[paste0(plot_name, "_download_png")]]<-downloadHandler(
		      		filename = paste0("carcinogenome_download_", header, "_", markerid, ".png"),
		      		content = function(file){
		      			ggsave(file, device = "png", grid.arrange(p1, p2, p3),
		      				width = w/130, height = h/35, 
		      				dpi = 300)
		      	})

			} else if (plot == "Boxplot"){
			
				w<-600
				h<-400
				
				output[[paste0(plot_name, "_", 1)]]<-renderPlot(p1,
		      		width = w, height = h
		      		)

		      	output[[paste0(plot_name, "_", 2)]]<-renderPlot(p2,
					width = w, height = h
		      		)

		      	output[[paste0(plot_name, "_", 3)]]<-renderPlot(p3,
					width = w, height = h
		      		)

		      	output[[paste0(plot_name, "_download_pdf")]]<-downloadHandler(
	      		filename = paste0("carcinogenome_download_", header, "_", markerid, ".pdf"),
	      		content = function(file){
	      			ggsave(file, device = "pdf", grid.arrange(p1, p2, p3),
	      				width = w/70, height = h/35, units = "in",
	      				dpi = 150)
	      		})

		      	output[[paste0(plot_name, "_download_png")]]<-downloadHandler(
		      		filename = paste0("carcinogenome_download_", header, "_", markerid, ".png"),
		      		content = function(file){
		      			ggsave(file, device = "png", grid.arrange(p1, p2, p3),
		      				width = w/70, height = h/35, units = "in", 
		      				dpi = 150)
		      	})
	      	}  	
		}

		plot_table<-function(plot_name, eset, header, markerid){
			tab<-get_de_by_gene_table(markerid, eset, 
				annot_prof = dat[["Profile Annotation"]], 
				match_id = "sig_id",  
				header = header,
				tas = input$marker_tas)

			output[[plot_name]]<-customTable(tab)
			
		}

		plot_text<-function(plot_name, header){
			output[[plot_name]]<-renderText(paste0("Table of profiles ranked by ", header))
		}

		#Annotation Tables
		output$t1<-customTable(dat[["Chemical Annotation"]])
		output$t2<-customTable(dat[["Profile Annotation"]])

		#Chemical Explorer Tables
		output$t3<-DT::renderDataTable(
			get_chemical_description(input$chem),
			options = list(dom = ''), 
			rownames = FALSE
			)

        observeEvent(c(input$landmark_de,
        	input$summarizefunc_de,
        	input$filterbyinput_de,
        	input$range_de,
        	input$numberthresleft_de, 
        	input$numberthresright_de, 
        	input$chem), {
	      		output$t4<-customTable(
		      		get_de_call()
	      			)
            }, once = FALSE)

           observeEvent(c(input$restore), {
            	updateCheckboxInput(session, inputId = "landmark_de", value = defaults[["landmark_de"]])
            	updateSelectInput(session, inputId = "summarizefunc_de", selected = defaults[["summarizefunc_de"]])
				updateCheckboxGroupInput(session, inputId = "filterbyinput_de", selected = defaults[["filterbyinput_de"]])
				updateSelectInput(session, inputId = "range_de", selected = defaults[["range_de"]])
				updateSliderInput(session, inputId = "numberthresleft_de", value = defaults[["numberthresleft_de"]])
				updateSliderInput(session, inputId = "numberthresright_de", value = defaults[["numberthresright_de"]])
           })

        observeEvent(input$de_hide, ({
    	updateCollapse(session, "de_opt_panel", close = "Options")
   			}))
        observeEvent(input$de_show, ({
    	updateCollapse(session, "de_opt_panel", open = "Options")
   			}))

      	observeEvent(c(input$chem,
        	input$gsname, 
        	input$gsmethod, 
        	input$summarize_gs),
        	output$t5<-customTable(get_gs_call())
        	)

      	observeEvent(c(input$chem,
        	input$conn_name, 
        	input$summarizefunc_conn),
        	output$t6<-customTable(get_conn_call())
        	)

      	#Marker Explorer plots
		observeEvent(c(input$marker_gene, input$marker_view, input$marker_tas), {

			es<-get_de_eset(annot_prof = dat[["Profile Annotation"]], 
				match_id = "sig_id")
			es<-es[, pData(es)[, "TAS"]>input$marker_tas]
	    	header<-"mod Z-scores"
	    	plot_hist(plot_name = "t7",
			eset = es,
			header = header, 
			markerid = input$marker_gene, plot = input$marker_view)

	    	plot_text(plot_name = "t7_text", header = header)

	    	plot_table(plot_name = "t7_table",
			eset = es,
			header = header, 
			markerid = input$marker_gene)


            }, once = FALSE)

		observeEvent(c(input$marker_gs, input$marker_view, input$marker_tas), {
			es<-get_gs_eset(gslist = dat[["Gene Set Enrichment"]], 
			 	gsname = input$marker_gsname, 
			 	gsmethod = input$marker_gsmethod, 
			 	annot_prof = dat[["Profile Annotation"]],
			 	match_id = "sig_id")
			es<-es[, pData(es)[, "TAS"]>input$marker_tas]
			header<-"Gene Set Scores"

	    	plot_hist(plot_name = "t8",
			eset = es, 
			header = header, 
			markerid = input$marker_gs, 
			plot = input$marker_view)

	    	plot_text(plot_name = "t8_text", header = header)

			plot_table(plot_name = "t8_table",
			eset = es,
			header = header, 
			markerid = input$marker_gs)

            }, once = FALSE)

		observeEvent(c(input$marker_conn, input$marker_view, input$marker_tas), {
			es<-get_conn_eset(connlist = dat[["Connectivity"]], 
 		 		conn_name = input$marker_conn_name, 
 		 		annot_prof = dat[["Profile Annotation"]],
			 	match_id = "sig_id")
			es<-es[, pData(es)[, "TAS"]>input$marker_tas]
			header<- "Connectivity Score (Percentile)"
	    	
	    	plot_hist(plot_name = "t9",
			eset = es, 
			header = header, 
			markerid = input$marker_conn, 
			plot = input$marker_view)

	    	plot_text(plot_name = "t9_text", header = header)

			plot_table(plot_name = "t9_table",
			eset = es,
			header = header, 
			markerid = input$marker_conn)

            }, once = FALSE)
		

		observe({
    		updateSelectizeInput(session, 
    			'marker_gene', 
    			choices = sort(rownames(dat[["Gene Expression"]]))
    		)
  		})

        observe({
    		updateSelectizeInput(session, 
    			'marker_gs', 
    			choices = 
 			sort(rownames(get_gs_eset(gslist = dat[["Gene Set Enrichment"]], 
 				gsname = input$marker_gsname, 
 				gsmethod = input$marker_gsmethod, 
 				annot_prof = dat[["Profile Annotation"]],
				match_id = "sig_id"))
    			))
  		})
  		
  		observe({
    		updateSelectizeInput(session, 
    			'marker_conn', 
    			choices = 
 			rownames(get_conn_eset(connlist = dat[["Connectivity"]], 
 				conn_name = input$marker_conn_name, 
 				annot_prof = dat[["Profile Annotation"]],
				match_id = "sig_id")
    			))
  		})

  		##Heatmap marker Exploper
  		output$warning<-renderText("Warning: TAS < 0.2 is slow to load!")

		observeEvent(c(input$marker_hm, 
			input$marker_tas_hm
			), {
			if(input$marker_hm == "Genes (Landmark)"){
				es<-get_de_eset(annot_prof = dat[["Profile Annotation"]], 
					match_id = "sig_id", hm = TRUE)
				plot_heatmap_static("t10", es, input$marker_tas_hm)
				}
            }, once = FALSE)

  		##Heatmap marker Exploper
		observeEvent(c(input$marker_hm, 
			input$marker_tas_hm, 
			input$marker_gsname_hm, 
			input$marker_gsmethod_hm
			), {
			if(input$marker_hm == "Gene Sets"){
				
				output$morpheus_result_link <- renderUI({
				  
				  htmlfile <- paste0("JSON/", sub(" .*", "", dat$title), "/", dsmap[[input$marker_gsname_hm]], "_", input$marker_gsmethod_hm, ".html")
				  
				  a(href=htmlfile, target="_blank", alt="heatmap", width = "100%", height="auto", span(icon("fas fa-arrow-circle-right"), "Open interactive heatmap in Morpheus"))
				  
				 	# paste(c('<a target="_blank" href="',
				 	# 	get_morpheus_link(url =
				 	# 		get_heatmap_gct(ds=input$marker_gsname_hm, 
				 	# 			dsmap = dsmap, 
				 	# 		method = input$marker_gsmethod_hm), 
				 	# 		domain = domain, tas = input$marker_tas_hm),
				 	# 	'">', 'Open interactive heatmap in Morpheus', '</a>'), sep = "")
				 	
				 	})

				es<-get_gs_eset(gslist = dat[["Gene Set Enrichment"]], 
				 	gsname = input$marker_gsname_hm, 
				 	gsmethod = input$marker_gsmethod_hm, 
				 	annot_prof = dat[["Profile Annotation"]],
				 	match_id = "sig_id")
				plot_heatmap_static("t11", es, input$marker_tas_hm)
				} 
            }, once = FALSE)
		
		observeEvent(c(input$marker_hm, 
			input$marker_tas_hm, 
			input$marker_conn_name_hm
			), {
			if (input$marker_hm == "CMap Connectivity"){
				es<-get_conn_eset(connlist = dat[["Connectivity"]], 
	 		 		conn_name = input$marker_conn_name_hm, 
	 		 		annot_prof = dat[["Profile Annotation"]],
				 	match_id = "sig_id")
				plot_heatmap_static("t12", es, input$marker_tas_hm)
				}
            }, once = FALSE)

		})	
)