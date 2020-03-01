## Functions for enhancedGraphics app with RCy3

#' EnhanceGraphics Pie Chart
#'
#' @description Prepares syntax for enhancedGraphics pie charts.
#' @param columns List of column names
#' @param colors Columns-matched list of pairs of colors to associate with minimum and
#' maximum values per column. Up to 4 brewer palette color pairs are 
#' used by default.
#' @param range List of two values: global minimum and maximum to be used for color 
#' gradients. Min and max across all columns are used by default.
#' @param portions: Columns-matched list of values to define pie slices. Pie is divided
#' into equal slices by default.
#' @param arcstart: Position to start pie slices. 0 is 9 o'clock position. Default is
#' 90 or 6 o'clock position.
#' @param scale: Size of graphic relative to node. Default is 1.0.
#' @param position: Position of graphic relative to node. Default is center.
#' @param showlabels: Display labels on pie slices. Default is FALSE.
#' @param labellist: Columns-matched list of labels for slices. Default is column names.
#' @param labelsize: Font size for pie slice labels. Default is 12.
#' @param nodelist: Which nodes to apply enhanded Graphics. Default is all nodes by SUID. 
#' If a list is provided, then it match node names. If a string prefixed
#' with "column:" is provided, then it will match non-null values for 
#' that column name.
#' @param newcolumn: Name of column to add enhancedGraphics syntax for style mapping.
#' @param style: Name of style to add mapping to. Default is "default". 
#' @param network: Name or SUID of network to apply enhancedGraphics.
#' @param base.url: Ignore unless you need to specify a custom domain, port or version
#' of CyREST.
#' @examples
#   egPie(c("HEKScore","JurkatScore"), range=c(0,1))
#   egPie(c("HEKScore","JurkatScore"), range=c(0,1), nodelist="column:UniProt", style="AP-MS Style")
#   egPie(c("HEKScore","JurkatScore"), c("white", "blue", "white", "purple"), range=c(0,1), nodelist="column:UniProt", style="AP-MS Style")
#
egPie <- function(columns, colors=NULL, range=NULL, portions=NULL, arcstart=90,
                  scale=1.0, position="center", showlabels=FALSE, 
                  labellist=NULL, labelsize=12, 
                  nodelist=NULL, newcolumn="egPie", style="default", 
                  network=NULL, base.url = 'http://localhost:1234/v1') {
  
  require(RCy3)
  
  net.colnames<-getTableColumnNames(network = network, base.url = base.url)
  not.col<-setdiff(columns, net.colnames)
  if (length(not.col) > 0)
    stop(sprintf("ERROR: %s column is not found in node table.", not.col))
  
  net.cols<-getTableColumns(columns = columns, 
                            network = network, base.url = base.url)
  col.min<-unlist(lapply(columns, function(x){
    if (is.null(range))
      min(net.cols[,x], na.rm=TRUE)
    else
      min(range)
  }))
  col.max<-unlist(lapply(columns, function(x){
    if (is.null(range))
      max(net.cols[,x], na.rm=TRUE)
    else
      max(range)
  }))
  
  if(is.null(colors)){
    brewer.cols1<-c('#ffffff','#3182bd', '#ffffff','#e6550d',
                    '#ffffff','#756bb1', '#ffffff','#31a354')
    brewer.cols2<-c('#4393C3','#D6604D', '#4393C3','#D6604D',
                    '#4393C3','#D6604D', '#4393C3','#D6604D')
    if(min(col.min) < 0)
      colors<-brewer.cols2[1:(2*length(columns))]
    else
      colors<-brewer.cols1[1:(2*length(columns))]
  } 

  colors1<-colors[c(TRUE,FALSE)]
  colors2<-colors[c(FALSE,TRUE)]
  
  if (length(colors) != length(columns)*2)
    stop(sprintf("ERROR: List of colors must be twice the 
                 length of list of columns."))
  
  colors<-unlist(lapply(columns, function(x){
    if (min(col.min) < 0)
      paste0('down:',colors1[match(x,columns)],
             ',zero:#ffffff',
             ',up:', colors2[match(x,columns)])
    else
      paste0('down:',colors1[match(x,columns)],
             ',zero:',colors1[match(x,columns)],
             ',up:',colors2[match(x,columns)])
  }))
  
  if(is.null(range))
    range<-c(min(col.min),max(col.max))

  if(is.null(portions))
    portions<-rep(round(1/length(columns),3),length(columns))  
  
  if(is.null(labellist))
    labellist<-columns
  
  if(showlabels){
    labellist<-paste0(' labellist="',paste(labellist, collapse = ','),'"')
    labelsize<-paste0(' labelsize="',labelsize,'"')
  }  else {
    labellist<-''
    labelsize<-''
  }
  
  eg.str<- paste0('piechart:',
                  ' attributelist="',paste(columns, collapse = ','),'"',
                  ' colorlist="', paste(colors, collapse = ';') ,'"',
                  ' range="',paste(range, collapse = ','),'"',
                  ' valuelist="',paste(portions, collapse = ','),'"',
                  ' arcstart="',arcstart,'"',
                  ' showlabels="',showlabels,'"',
                  toString(labellist),
                  toString(labelsize),
                  ' position="',paste(position, collapse = ','),'"',
                  ' scale="',formatC(scale,digits=1,format="f"),'"')
    
  # Populate column and map to style
  nodelist.df <- NULL
  tkc <- 'SUID'
  if(is.null(nodelist)){ #working with all nodes
    nodelist.df <- net.cols
  } else if (startsWith(nodelist, "column:")){ #working with column name
    nodelist<-strsplit(nodelist, "column:")[[1]][2]
    nodelist.df<-getTableColumns(columns=c(nodelist, "SUID"), network = network, 
                                 base.url = base.url)
    nodelist.df<-nodelist.df[which(nodelist.df[nodelist]!=""),]
  } else { #working with list of node names
    nodelist.df<-data.frame(nodelist)
    row.names(nodelist.df)<-nodelist
    tkc<-'name'
  }
    
  egpie.data <- cbind(nodelist.df, eg.str)
  egpie.data <- egpie.data['eg.str']
  names(egpie.data)<-newcolumn
  cyrestPOST(paste('networks',getNetworkSuid(),'tables','defaultnode','columns',sep = '/' ),
             body=list(name=newcolumn,
                       type="String"),
             base.url=base.url)
  loadTableData(egpie.data, table.key.column = tkc, 
                network = network, base.url = base.url)
  egpie.map <- mapVisualProperty('NODE_CUSTOMGRAPHICS_1', newcolumn, 'p',
                                 network = network, base.url = base.url)
  updateStyleMapping(style, egpie.map, base.url = base.url)
}


