library(dplyr)
library(ROCR)
library(glmnet)
library(biglasso)
library(readr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(plotly)
library(tableHTML)
library(htmlwidgets)

#' Evaluation Function 
#' 
#' This function evaluates the performance of at least one PRS model on the same 
#' outcome by utilizing various evaluation methods, depending on whether the outcome is binary
#' or continuous. 
#' 
#' @param phenotype matrix containing the phenotype information
#' @param ... dataframes that contain one set of polygenic risk scores each. 
#' Please provide at least one set of polygenic risk scores. You can provide the 
#' name of each set of polygenic risk score by providing a name for the parameter 
#' in this format: \code{model_evaluation(phenotype=phenotype, 'Pre-trained PRS'=df1, 'ABPRS'=df2)}. 
#' @param binary boolean that determines whether the data is binary or not 
#' @param filename the name of the exported file
#' @return An html file with the filename will appear in the current directory,
#' which contains different plots and data depending on whether the outcome is binary
#' or continuous. 
#' 
#' For Binary Outcomes:
#' - a downloadable performance score table with AUC scores for each set of PRSs
#' - an interactable performance score comparison bar plot
#' - an interactable polygenic risk score distribution density plot
#' - an interactable percentage of cases vs polygenic risk score graph
#' - an interactable odds ratio plot
#' - the legend, ggplot code, and downloadable dataframe that produced the above plots
#' 
#' 
#' For Continuous Outcomes: 
#' - a downloadable performance score table with the Mean-Squared Error and \eqn{R^2}
#' for each set of PRSs
#' - an interactable performance score comparison bar plot
#' - an interactable mean phenotype vs polygenic risk score percentile graph
#' - an interactable mean difference plot
#' - the legend, ggplot code, downloadable dataframe that produced the above plots
#' 
#' @export
#' 
model_evaluation <- function(phenotype, ..., binary, bin=10){
  
  #Retrieve PRS Scores
  all_prs <- list(...)
  all_prs <- lapply(all_prs, function(x) as.data.frame(x))
  # Ensure at least one PRS score is provided
  if (length(all_prs) == 0) {
    stop("At least one PRS score must be provided.")
  }
  
  #Retrieve Names 
  call <- match.call()
  call_list <- as.list(call)
  all_names <- names(call_list)[-c(1,2)] 
  all_names <- all_names[-c(length(all_names), length(all_names)-1)]
  if(any(nchar(all_names)==0)){
    all_names <- paste0("PRS", 1:length(all_prs))
  }
  
  style <- "<style>
  body { font-family: Arial, sans-serif; }
  div {max-width:800px; margin-inline: auto; margin-top:10px;}
  pre { background-color: #f4f4f4; padding: 10px; border: 1px solid #ddd; }
  iframe { display: block; margin: 0 auto; border:none;}
  button { background-color: lightblue; border: none; color: black; 
  padding: 10px 20px; text-align: center; text-decoration: none; 
  display: inline-block; font-size: 16px; margin: 4px 2px; cursor: pointer; }
  button:hover { background-color: #2c7bb6; }
  table { width: 50%; border: 1px solid;}
  th { background-color: lightblue; border: 1px solid; padding: 8px;}
  td {border: 1px solid; padding: 8px;}
</style>"
  
  if(binary){
    html_binary(phenotype, all_prs, all_names, style, bin=bin)
  }else{
    html_continuous(phenotype, all_prs, all_names, style, bin=bin)
  }
  
  return("Done")
}

html_binary<- function(phenotype, all_prs, all_names, style, bin){
  
  n_model <- length(all_names)
  
  # PLOT 1: Performance Scores Comparison
  all_auc <- sapply(all_prs, function(prs) AUC_Score(prs, phenotype))
  PerformanceScores <- data.frame(Model=all_names, AUC=all_auc)
  ylab <- "AUC"
  BarPlot <- Plot_Score(PerformanceScores, ylab)
  barplot_code <- "
   plot <- ggplot(data=Scores, aes(x=reorder(Model, +Score), y=Score, fill=Model)) +
    geom_bar(stat=\"identity\", width=0.5) +
    geom_text(aes(label=round(Score, 2)), size=3.5, vjust=-0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    labs(x=\"Model\", y=ylab, fill=\"Model\") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = \"black\")) +
    scale_fill_manual(values = setNames(c(\"#ef8a62\", \"#67a9cf\"), c(name1, name2)))"
  saveWidget(BarPlot, "BarPlot.html", selfcontained = TRUE)
  PerformanceScores[,2]<- round(PerformanceScores[,2], 3)
  scores_csv <- convert_df_to_string(PerformanceScores)
  
  # PLOT 2: Density Plots
  Phenotype<-as.factor(rep(phenotype,n_model))
  Standardized_PRS<-unlist(lapply(all_prs, scale))
  Method<-rep(all_names, each=length(phenotype))
  PRSTable<-data.frame(Phenotype,Standardized_PRS,Method)
  DensityPlot <- Density_Plot(PRSTable, all_names)
  densityplot_code <- "
    StandardizedPRS$Method <- factor(StandardizedPRS$Method, levels = c(name1, name2))
    plot<-ggplot(StandardizedPRS,aes(x=PRS, fill=Phenotype))+
      geom_density(alpha=0.5)+
      scale_fill_manual(values = c(\"#67a9cf\", \"#ef8a62\"))+
    labs(x=\"Standardized PRS\", y=\"Density\") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x = element_line(colour = \"black\")) + 
    facet_grid(cols=vars(Method),scales=\"free_x\")"
  saveWidget(DensityPlot, "DensityPlot.html", selfcontained = TRUE)
  prs_csv <- convert_df_to_string(PRSTable)
  
  # PLOT 3: Percentage of Cases vs. PRS Percentile
  prev_data = Prevalence_Data(phenotype, all_prs, all_names, bin)
  prev_csv <- convert_df_to_string(prev_data)
  PrevalencePlot = Prevalence_Plot(prev_data)
  prevalence_code <- "
  plot <- ggplot(PrevalenceData, aes(x = Percentile, y = Prevalence, color=Model)) +
    geom_point(size = 3) +
    labs(x = \"Risk Score Percentile\", y = \"Case Prevalence\",
         title = \"Prevalence vs. Risk Score Percentile\") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = \"black\"))"
  saveWidget(PrevalencePlot, "PrevalencePlot.html", selfcontained = TRUE)
  
  # PLOT 4: Odds Ratio Plot
  #TBA
  
  cat(
    paste(

"<!DOCTYPE html>
<html>
<head>
<title>Evaluation</title>",

style,

"<script>

var csv1 = \"", scores_csv, "\"
var filename1 = \"PerformanceScores.csv\"

var csv2 = \"", prs_csv, "\"
var filename2 = \"StandardizedPRS.csv\"

var csv3 = \"", prev_csv, "\"
var filename3 = \"PrevalenceData.csv\"

function download_csv_file(csv, filename) {
    var hiddenElement = document.createElement('a');
    hiddenElement.href = 'data:text/csv;charset=utf-8,' + encodeURI(csv);
    hiddenElement.target = '_blank';
    
    //provide the name for the CSV file to be downloaded
    hiddenElement.download = filename;
    hiddenElement.click();
}
</script>

</head>

<body>

<h1><center>Evaluation Results</center></h1>

<h2>Performance Score Table</h2>",
tableHTML(PerformanceScores),
"<button onclick=\"download_csv_file(csv1, filename1)\">Download PerformanceScores.csv</button>",

"<h2>Performance Score Comparisons</h2>",
"<iframe src='BarPlot.html' width='800' height='500'></iframe>",
"<div>This barplot compares the AUC performance score between different 
sets of PRSs derived from different models. The models are ordered by increasing 
performance score. </div>",
"<pre><code>", barplot_code, "</code></pre>",

"<h2>Polygenic Risk Score Distributions</h2>
<iframe src='DensityPlot.html' width='800' height='500'></iframe>",
"<div>The figure shows the density curves of the PRS of control (0) and case (1) 
phenotype for each set of PRS. The PRSs are standardized with a mean of 0 and 
standard deviation of 1 for each set. The goal of this figure is to see how well
the polygenic risk scores can distinguish between case and control. </div>",
"<button onclick=\"download_csv_file(csv2, filename2)\">Download StandardizedPRS.csv</button>",
"<pre><code>", densityplot_code, "</code></pre>",

"<h2>Phenotype vs PRS Percentile</h2>
 <iframe src='PrevalencePlot.html' width='800' height='500'></iframe>",
"<div>The figure plots percentage of cases(prevalence) against the risk score percentile
for PRSs derived from different models. For each model, 15 quantiles are plotted 
in the graph. A model that performs better should have a higher prevalence in the
higher risk score percentiles and a lower prevalence in the lower percentiles.  
</div>",
"<button onclick=\"download_csv_file(csv3, filename3)\">Download PrevalenceData.csv</button>",
"<pre><code>", prevalence_code, "</code></pre>",

"</body>
</html>"), 

file = "Evaluation.html")
}

html_continuous <- function(phenotype, all_prs, all_names, style, bin){
  
  n_model <- length(all_names)
  
  # Plot1: Performance Score Plot
  all_mse <- sapply(all_prs, function(prs) MSE_Score(prs, phenotype))
  all_r2 <- sapply(all_prs, function(prs) RSquared_Score(prs, phenotype))
  PerformanceScores <- data.frame(Model=all_names, MSE=all_mse, `R-Squared`=all_r2)
  ylab = "MSE"
  BarPlot <- Plot_Score(PerformanceScores[,c(1:2)], ylab)
  barplot_code <- "
   plot <- ggplot(data=Scores, aes(x=reorder(Model, +Score), y=Score, fill=Model)) +
    geom_bar(stat=\"identity\", width=0.5) +
    geom_text(aes(label=round(Score, 2)), size=3.5, vjust=-0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    labs(x=\"Model\", y=ylab, fill=\"Model\") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = \"black\"))"
  saveWidget(BarPlot, "BarPlot.html", selfcontained = TRUE)
  PerformanceScores[,-1]<- round(PerformanceScores[,-1], 3)
  scores_csv <- convert_df_to_string(PerformanceScores)
  
  # Plot 2: Mean of Phenotype vs. Risk Score Percentile
  mean_data = Mean_Data(phenotype, all_prs, all_names, bin)
  mean_csv <- convert_df_to_string(mean_data)
  MeanPlot = Mean_Plot(mean_data)
  mean_code <- "
  plot <- ggplot(MeanData, aes(x = Percentile, y = Prevalence, color=Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, alpha=0.8) +
    labs(x = \"Risk Score Percentile\", y = \"Mean of Phenotype\",
         title = \"Mean of Phenotype vs. Risk Score Percentile\") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = \"black\"))"
  saveWidget(MeanPlot, "MeanPlot.html", selfcontained = TRUE)

  # Plot 3: Mean Difference Plot
  MeanDiff <- data.frame(TBA=c("TBA", "TBA"))
  
  cat(
    paste(
      
"<!DOCTYPE html>
<html>
<head>
<title>Evaluation</title>",
      
style,
      
"<script>

var csv1 = \"", scores_csv, "\"
var filename1 = \"PerformanceScores.csv\"

var csv2 = \"", mean_csv, "\"
var filename2 = \"MeanData.csv\"

function download_csv_file(csv, filename) {
    var hiddenElement = document.createElement('a');
    hiddenElement.href = 'data:text/csv;charset=utf-8,' + encodeURI(csv);
    hiddenElement.target = '_blank';
    
    //provide the name for the CSV file to be downloaded
    hiddenElement.download = filename;
    hiddenElement.click();
}
</script>

</head>

<body>

<h1><center>Evaluation Results</center></h1>

<h2>Performance Score Table</h2>",
tableHTML(PerformanceScores),
"<button onclick=\"download_csv_file(csv1, filename1)\">Download PerformanceScores.csv</button>",

"<h2>Performance Score Comparisons</h2>",
"<iframe src='BarPlot.html' width='800' height='500'></iframe>",
"<div>This barplot compares the mean squared error (MSE) performance scores between different 
sets of PRSs derived from different models. </div>",
"<pre><code>", barplot_code, "</code></pre>",

"<h2>Phenotype vs PRS Percentile</h2>
 <iframe src='MeanPlot.html' width='800' height='500'></iframe>",
"<div>This figure plots the mean of the continuous phenotype and +/- 
1 standard deviation in each polygenic risk score quantiles, for different 
PRS models. For each model,", bin, "quantiles are plotted in the graph. 
A model that performs better should have a higher mean in the
higher risk score percentiles and a lower mean in the lower percentiles.  
</div>",
"<button onclick=\"download_csv_file(csv2, filename2)\">Download MeanData.csv</button>",
"<pre><code>", mean_code, "</code></pre>",

"</body>
</html>"), 

file = "Evaluation.html")
}

convert_df_to_string <- function(df) {
  column_names <- colnames(df)
  header <- paste(column_names, collapse = ",")
  body <- apply(df, 1, function(row) paste(row, collapse = ",")) 
  result <- paste(header, paste(body, collapse = "\\n"), sep = "\\n")
  return(result)
}

Plot_Score <- function(PerformanceScores, ylab){
  colnames(PerformanceScores) = c("Model", "Score")
  PerformanceScores$Model <- factor(PerformanceScores$Model, levels = PerformanceScores$Model)
  
  plot <- ggplot(data=PerformanceScores, aes(x=Model, y=Score, fill=Model)) +
    geom_bar(stat="identity", width=0.5) +
    geom_text(aes(label=round(Score, 2)), size=3.5, vjust=-0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    labs(x="Model", y=ylab, fill="Model") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) #+
    #scale_fill_manual(values = setNames(c("#ef8a62", "#67a9cf"), c(name1, name2)))
  fig <- ggplotly(plot)
  return(fig)
}

#Calculate AUC score from data
AUC_Score<-function(prs, pheno){
  # Create prediction model, obtain predictions, calculate AUC score from prediction
  mod <- glm(pheno~., data=prs, family="binomial")
  prediction <- prediction(predict(mod, prs, type="response"), pheno)
  auc <- performance(prediction, measure="auc")@y.values[[1]]
  return(auc)
}

#Calculate MSE score from data
MSE_Score<-function(prs, pheno){
  mod <- lm(pheno~., data=prs)
  prediction <- (predict(mod, prs, type="response")- pheno)^2
  mse <- mean(prediction)
  return(mse)
}

#Calculate R Squared from data
RSquared_Score <- function(prs, pheno){
  mod <- lm(pheno~., data=prs)
  r2 <- summary(mod)$r.squared
  return(r2)
}

Density_Plot <- function(PRSTable, model_names){
  PRSTable$Method <- factor(PRSTable$Method, levels = model_names)
  plot<-ggplot(PRSTable,aes(x=Standardized_PRS, fill=Phenotype))+
    geom_density(alpha=0.5)+
    scale_fill_manual(values = c("#67a9cf", "#ef8a62"))+
    labs(x="Standardized PRS", y="Density") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line.x = element_line(colour = "black")) + 
    facet_grid(cols=vars(Method),scales="free_x")
  fig <- ggplotly(plot)
  return(fig)
}

Prevalence_Data <- function(pheno, all_prs, all_names, n){
  
  all_prevalence_data <- data.frame()
  
  # Loop over each PRS vector
  for (i in seq_along(all_prs)) {
    prs <- as.numeric(all_prs[[i]][,1])
    name <- all_names[i]
    
    percentile <- ecdf(as.matrix(prs))(as.matrix(prs)) * 100
    
    data <- data.frame(Percentile = percentile, Phenotype = pheno)
    
    # Calculate Prevalence
    prevalence <- data %>%
      mutate(Percentile = ntile(Percentile, n)/n) %>%
      group_by(Percentile) %>%
      summarize(Prevalence = mean(Phenotype)) %>%
      mutate(Model = name)
    
    # Combine the current prevalence data with the accumulated data
    all_prevalence_data <- rbind(all_prevalence_data, prevalence)
  }
  
  return(all_prevalence_data)
}

Prevalence_Plot <- function(PrevalenceData){
  
  plot <- ggplot(PrevalenceData, aes(x = Percentile, y = Prevalence, color=Model)) +
    geom_point(size = 3) +
    labs(x = "Risk Score Percentile", y = "Percentage of Cases",
         title = "Percentage of Cases vs. Risk Score Percentile")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) #+ 
    #scale_color_manual(values = setNames(c("#ef8a62", "#67a9cf"), c(name1, name2)))
  fig <- ggplotly(plot)
  
  return(fig)
}

Mean_Data <- function(pheno, all_prs, all_names, n){
  
  all_mean_data <- data.frame()
  
  # Loop over each PRS vector
  for (i in seq_along(all_prs)) {
    prs <- as.numeric(all_prs[[i]][,1])
    name <- all_names[i]
    
    percentile <- ecdf(as.matrix(prs))(as.matrix(prs)) * 100
    
    data <- data.frame(Percentile = percentile, Phenotype = pheno)
    
    # Calculate Mean
    mean <- data %>%
      mutate(Percentile = ntile(Percentile, n)/n) %>%
      group_by(Percentile) %>%
      summarize(Mean = mean(Phenotype), SD=sd(Phenotype)) %>%
      mutate(Model = name) 
    
    # Combine the current mean data with the accumulated data
    all_mean_data <- rbind(all_mean_data, mean)
  }
  
  return(all_mean_data)
}

Mean_Plot <- function(MeanData){
  
  plot <- ggplot(MeanData, aes(x = Percentile, y = Mean, color=Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0, alpha=0.8) +
    labs(x = "Risk Score Percentile", y = "Mean of Phenotype",
         title = "Mean of Phenotype vs. Risk Score Percentile")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  #scale_color_manual(values = setNames(c("#ef8a62", "#67a9cf"), c(name1, name2)))
  fig <- ggplotly(plot)
  
  return(fig)
}