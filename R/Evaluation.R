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

#' Evaluate Function 
#' @param pheno matrix containing the phenotype information
#' @param prs1 dataframe containing the first set of polygenic risk scores
#' @param prs2 dataframe containing the second set of polygenic risk scores
#' @param binary boolean that determines whether the data is binary or not 
#' @return A file named "Evaluation.html" would appear in the current directory, 
#' which contains:
#' - a downloadable performance score table 
#' - an interactable performance score comparison bar plot
#' - an interactable polygenic risk score distribution density plot
#' - an interactable case prevalence vs polygenic risk score graph
#' - the ggplot code and downloadable dataframe that produced the above plots
#' @export
#' 
Evaluate <- function(pheno, prs1, prs2, binary=TRUE){
  
  prs1 <- as.data.frame(prs1)
  prs2 <- as.data.frame(prs2)

  #Calculate Performance Score
  if(binary){
    mod1 <- glm(pheno~., data=prs1, family="binomial")
    mod2 <- glm(pheno~., data=prs2, family="binomial")
    score1 <- AUC_Score(prs1, pheno, mod1)
    score2 <- AUC_Score(prs2, pheno, mod2)
    Scores <- data.frame(Model=c("Model 1","Model 2"), AUC=c(score1, score2))
    ylab <- "AUC"
  }else{
    mod1 <- glm(pheno~., data=prs1, family="gaussian")
    mod2 <- glm(pheno~., data=prs2, family="gaussian")
    score_MSE1 <- MSE_Score(df1, pheno, mod1)
    score_MSE2 <- MSE_Score(df2, pheno, mod2)
    score_R1 <- summary(mod1)$r.squared
    score_R2 <- summary(mod2)$r.squared
    Scores <- data.frame(Model=c("Model 1","Model 2"), MSE=c(score_MSE1, score_MSE2), R_Squared=c(score_R1, score_R2))
    ylab = "MSE"
  }
  
  convert_df_to_string <- function(df) {
    column_names <- colnames(df)
    header <- paste(column_names, collapse = ",")
    body <- apply(df, 1, function(row) paste(row, collapse = ",")) 
    result <- paste(header, paste(body, collapse = "\\n"), sep = "\\n")
    return(result)
  }
  
  #BarPlot 
  BarPlot <- Plot_Score(Scores, ylab)
  saveWidget(BarPlot, "BarPlot.html", selfcontained = TRUE)
  
  #Change Score to 3 dp and convert to CSV string
  Scores[,2] = round(Scores[,2],3)
  scores_csv <- convert_df_to_string(Scores)
  
  #DensityPlot
  DensityPlot = Density_Plot(pheno, prs1, prs2)
  saveWidget(DensityPlot, "DensityPlot.html", selfcontained = TRUE)
  PRS <- data.frame(PRS1=prs1, PRS2=prs2)
  prs_csv <- convert_df_to_string(PRS)
  
  #PrevalencePlot
  prev_data = Prevalence_Data(pheno, prs1, prs2, 15)
  prev_csv <- convert_df_to_string(prev_data)
  PrevalencePlot = Prevalence_Plot(prev_data)
  saveWidget(PrevalencePlot, "PrevalencePlot.html", selfcontained = TRUE)

  # ggplot Code String
  barplot_code <- "
  p <- ggplot(data = df, aes(x = reorder(Model, +score), y = score, fill = Model))
  plot <- p +
    geom_bar(stat = 'identity', width = 0.5) +
    geom_text(aes(label = round(score, 2)), size = 3.5, vjust = -0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    labs(x = 'Model', y = ylab, fill = 'Model') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
    scale_fill_manual(values = c('#ef8a62', '#67a9cf'))"
  
  densityplot_code <- "
  df<-data.frame(Phenotype,PRS,Method)
  p<-ggplot(df,aes(x=PRS, fill=Phenotype))
  plot<- p +
    geom_density(alpha=0.5) +
    scale_fill_manual(values = c('#ef8a62', '#67a9cf')) +
    labs(x='Standardized PRS', y='Density') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line.x = element_line(colour = 'black')) + 
    facet_grid(cols=vars(Method),scales='free_x')"
  
  prevalence_code <- "
  plot <- ggplot(Prevalence_Data, aes(x = Percentile, y = Prevalence, color=Type)) +
    geom_point(size = 3) +
    labs(x = 'Risk Score Percentile', y = 'Case Prevalence',
    title = 'Case Prevalence vs. Risk Score Percentile') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = 'black')) + 
    scale_color_manual(values = c('PRS1' = '#ef8a62', 'PRS2' = '#67a9cf'))"
  
  cat(
    paste(
      
      "<!DOCTYPE html>
<html>
<head>
<title>Evaluation</title>

<style>

  body { font-family: Arial, sans-serif; }
  pre { background-color: #f4f4f4; padding: 10px; border: 1px solid #ddd; }
  iframe { display: block; margin: 0 auto; border:none;}
  
  button {
      background-color: lightblue; /* Blue background */
        border: none; /* Remove borders */
        color: black; /* Black text */
        padding: 10px 20px; /* Some padding */
        text-align: center; /* Centered text */
        text-decoration: none; /* Remove underline */
        display: inline-block; /* Inline block */
        font-size: 16px; /* Set font size */
        margin: 4px 2px; /* Some margins */
        cursor: pointer; /* Pointer cursor on hover */
    }
  button:hover {
    background-color: #2c7bb6; /* Darker blue on hover */
  }
  
  table {
    width: 50%;
    border: 1px solid;
  }
  
  th {
    background-color: lightblue; 
    border: 1px solid; 
    padding: 8px;
  }
  
  td {
    border: 1px solid; 
    padding: 8px;
  }
  
</style>

<script>

var csv1 = \"", scores_csv, "\"
var filename1 = \"Score.csv\"

var csv2 = \"", prs_csv, "\"
var filename2 = \"prs.csv\"

var csv3 = \"", prev_csv, "\"
var filename3 = \"Prevalence.csv\"

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
tableHTML(Scores),
"<button onclick=\"download_csv_file(csv1, filename1)\">Download Score.csv</button>",

"<h2>Performance Score Comparisons</h2>",
"<iframe src='BarPlot.html' width='800' height='500'></iframe>",
"<pre><code>", barplot_code, "</code></pre>",

"<h2>Polygenic Risk Score Distributions</h2>
 <iframe src='DensityPlot.html' width='800' height='500'></iframe>",
"<button onclick=\"download_csv_file(csv2, filename2)\">Download prs.csv</button>",
"<pre><code>", densityplot_code, "</code></pre>",

"<h2>Case Prevalence vs. Polygenic Risk Score</h2>
 <iframe src='PrevalencePlot.html' width='800' height='500'></iframe>",
"<button onclick=\"download_csv_file(csv3, filename3)\">Download Prevalence.csv</button>",
"<pre><code>", prevalence_code, "</code></pre>",

"</body>
</html>"), 

file = "Evaluation.html")
  
  return("Done")
}

#Calculate AUC score from data
AUC_Score<-function(prs, pheno, mod){
  # Obtain predictions
  prediction <- prediction(predict(mod, prs, type="response"), pheno)
  # Find auc score from prediction
  auc <- performance(prediction, measure="auc")@y.values[[1]]
  return(auc)
}

#Calculate MSE score from data
MSE_Score<-function(prs, pheno, mod){
  prediction <- (predict(mod, x, type="response")- y)^2
  mse <- mean(prediction)
  return(mse)
}

Plot_Score <- function(df, ylab){
  colnames(df) = c("Model", "Score")
  
  p<-ggplot(data=df, aes(x=reorder(Model, +Score), y=Score, fill=Model))
  
  plot <- p +
    geom_bar(stat="identity", width=0.5) +
    geom_text(aes(label=round(Score, 2)), size=3.5, vjust=-0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    labs(x="Model", y=ylab, fill="Model") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    #theme(element_text(margin=margin(b=50))) +
    scale_fill_manual(values = c("#ef8a62", "#67a9cf"))
  
  fig <- ggplotly(plot)
  
  return(fig)
}

Density_Plot <- function(pheno, prs1, prs2){
  Phenotype<-as.factor(c(pheno, pheno))
  PRS<-c(c(scale(prs1)),c(scale(prs2)))
  Method<-rep(c("PRS1","PRS2"), each=length(pheno))
  df<-data.frame(Phenotype,PRS,Method)
  p<-ggplot(df,aes(x=PRS, fill=Phenotype))
  plot<-p+geom_density(alpha=0.5)+
    scale_fill_manual(values = c("#ef8a62", "#67a9cf"))+
    labs(x="Standardized PRS", y="Density") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line.x = element_line(colour = "black")) + 
    facet_grid(cols=vars(Method),scales="free_x")
  fig <- ggplotly(plot)
  return(fig)
}

Prevalence_Data <- function(pheno, prs1, prs2, n){
  
  percentile1 <- ecdf(as.matrix(prs1))(as.matrix(prs1)) * 100
  percentile2 <- ecdf(as.matrix(prs2))(as.matrix(prs2)) * 100
  
  data <- data.frame(Percentile1 = percentile1, Percentile2 = percentile2, 
                     Phenotype = pheno)
  
  #Calculate Prevalance
  prevalence1 <- data %>%
    mutate(Percentile = ntile(Percentile1, n)/n) %>%
    group_by(Percentile) %>%
    summarize(Prevalence = mean(Phenotype))
  
  prevalence2 <- data %>%
    mutate(Percentile = ntile(Percentile2, n)/n) %>%
    group_by(Percentile) %>%
    summarize(Prevalence = mean(Phenotype))
  
  
  Prevalence_Data <- rbind(prevalence1, prevalence2)
  Type = rep(c("PRS1", "PRS2"), each=n)
  Prevalence_Data <- cbind(Prevalence_Data, Type=Type)
  return(Prevalence_Data)
}

Prevalence_Plot <- function(Prevalence_Data){
  
  plot <- ggplot(Prevalence_Data, aes(x = Percentile, y = Prevalence, color=Type)) +
    geom_point(size = 3) +
    labs(x = "Risk Score Percentile", y = "Case Prevalence",
         title = "Case Prevalence vs. Risk Score Percentile") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_color_manual(values = c("PRS1" = "#ef8a62", "PRS2" = "#67a9cf"))
  
  fig <- ggplotly(plot)
  
  return(fig)
}
