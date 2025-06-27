# Source additional R scripts
source("01_library.R")
source("02_plot_utils.R")
source("03_utils.R")
source("04_GSEA.R")
source("05_ML_list.R")
source("06_MLL.R")

# Define UI
ui <- dashboardPage(skin = "black",
  header = dashboardHeader(title = "GENOXA", titleWidth = 278),
  ## Sidebar content
  sidebar = dashboardSidebar(
  width = 278,
  useShinyjs(),
  sidebarMenu(
  id = "tabs",
  menuItem("About", tabName = "home", icon = icon("home")),
  menuItem("User Manual", tabName = "manual",icon = icon("book")),
  menuItem("Select Data", expandedName = "upload", icon = icon("file"),
           radioButtons("dataset_choice", "Select Data:", 
                        choices = c("TCGA-PRAD", "TCGA-BRCA", "TCGA-LUAD",
                                    "TCGA-KIRC", "TCGA-COAD", "TCGA-STAD",
                                    "TCGA-BLCA", "TCGA-LIHC"),
                        selected = "TCGA-PRAD"),  # Default selection is TCGA-PRAD
           numericInput("sgn","Select gene counts (Total 2500) for data analysis:", min = 1, max = 9999999, value = 2000)
  ),
  menuItem("Data Preprocessing:",tabName="upload", icon = icon("cogs"),
           numericInput("dp1","Filter raw reads by rowSum:", min = 0.01, max = 9999999, value = 5),
           numericInput("dp2","Filter raw reads by FilterByExpr:", min = 0.01, max = 9999999, value = 5),
           selectInput("dp3","Select Normalization method for calcNormFactors:", choices = list("TMM","RLE","upperquartile","none"), selected = "TMM")
  ),
  menuItem("Differential expression analysis:",tabName="upload", icon = icon("tasks"),
           selectInput("refg", "Select reference group for DE:", choices = NULL),
           numericInput("dex1","Cutoff value for p-values and adjusted p-values:",min = 0.01, max = 99999, value = 0.05),
           numericInput("dex2","Minimum absolute log2-fold change (log2FC) to use as threshold for differential expression:",min = 1, max = 9999, value = 1),
           selectInput("dex3","Select adjustment method:", choices = list("none", "bonferroni", "holm", "BH", "fdr", "BY"), selected = "BH")
  ),
  menuItem("GO over-representation:",tabName="upload", icon = icon("codepen"),
           organismDbOptions,
           keytypeOptions,
           selectInput("ontt", "Select subontology:",
                       choices = c("BP", "MF", "CC", "ALL"), selected = "ALL"),
           numericInput("minGSSize", "Minimum Gene Set Size:", value = 10, min = 1, max = 9999999),
           numericInput("maxGSSize", "Maximum Gene Set Size:", value = 500, min = 1, max = 9999999),
           numericInput("pvalueCutoff", "P-Value Cutoff:", value = 0.05, min = 0, max = 9999999),
           numericInput("qvalueCutoff", "q-Value Cutoff:", value = 0.02, min = 0, max = 9999999),
           selectInput("pAdjustMethod", "Adjust Method:",
                       choices = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), selected = "BH")
  ),
  menuItem("Gene Modeling using ML:",tabName="upload", icon = icon("react"),
           h5(tags$b("*Select parameters for split data:")),
           numericInput("gms1","Enter The size of the training data set:",min = 0.01, max = 999999, value = 0.80),
           numericInput("gms2","Seed:",min = 1, max = 999999, value = 8314),
           h5(tags$b("*Select parameters for train models:")),
           selectInput("gmgmt1", "Select resampling method to use:",
                       c("None" = "none",
                         "boot" = "boot",
                         "boot632" = "boot632",
                         "optimism_boot" = "optimism_boot",
                         "boot_all" = "boot_all",
                         "cv" = "cv",
                         "repeatedcv" = "repeatedcv",
                         "LOOCV" = "LOOCV",
                         "LGOCV" = "LGOCV",
                         "oob" = "oob",
                         "timeslice" = "timeslice",
                         "adaptive_cv" = "adaptive_cv",
                         "adaptive_boot" = "adaptive_boot",
                         "adaptive_LGOCV" = "adaptive_LGOCV"), selected = "repeatedcv"),
           conditionalPanel(condition = "input.gmgmt1 == 'repeatedcv'",
                            numericInput("gmgmt2","The number of complete sets of folds to compute:",min = 1, max = 999999, value = 3)),
           numericInput("gmgmt3","Number of resampling iterations:",min = 1, max = 999999, value = 10),
           numericInput("gmgmt4","Random number seed:",min = 1, max = 999999, value = 351),
           ML_model_options,
           h5(tags$b("**Select minimum two and maximum")),
           h5(tags$b("four ML algorithms for better data")),
           h5(tags$b("visualization."))
  ),
  menuItem("Results", tabName = "results", icon = icon("sliders")),
  tags$br(),
  actionButton('subm', 'Start Analysis'),
  actionButton("reset_button", "Reset")
 ),

 tags$script(HTML("
    $(document).on('shiny:inputchanged', function(event) {
      if (event.name === 'gmgmt5') {
        var selectedOptions = $('#gmgmt5').val();
        var selectedCount = selectedOptions ? selectedOptions.length : 0;
        
        // Clear previous modals
        $('.modal').remove();
        
        if (selectedCount < 2) {
          Shiny.setInputValue('showModal', 'min', {priority: 'event'});
        } else if (selectedCount > 4) {
          Shiny.setInputValue('showModal', 'max', {priority: 'event'});
        } else {
          Shiny.setInputValue('showModal', 'none', {priority: 'event'});
        }
      }
    });

    $(document).on('shiny:inputchanged', function(event) {
      if (event.name === 'subm') {
        var selectedOptions = $('#gmgmt5').val();
        var selectedCount = selectedOptions ? selectedOptions.length : 0;
        
        if (selectedCount < 2 || selectedCount > 4) {
          event.preventDefault(); // Prevent form submission
          if (selectedCount < 2) {
            Shiny.setInputValue('showModal', 'min', {priority: 'event'});
          } else if (selectedCount > 4) {
            Shiny.setInputValue('showModal', 'max', {priority: 'event'});
          }
        } else {
          Shiny.setInputValue('showModal', 'none', {priority: 'event'}); // Hide modal if valid
        }
      }
    });
  "))

),

## Body content
body = dashboardBody(

# Include the CSS files  
 tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),

# General container for body content
 div(class = "container",
 tabItems(
 # First tab content
 tabItem(tabName = "home",
  fluidRow(
  box(
  title = "Welcome to GENOXA",
  width = 12,
  solidHeader = TRUE,
  status = "success",
  h2(tags$strong("GENOXA Ver 0.2")),
  tags$ul(
  tags$li(p("GENOXA is an interactive web application designed for analyzing bulk RNA-Seq data. It provides traditional bioinformatics results using limma and edgeR packages.")), 
  tags$li(p("GENOXA provides predictive models using DE genes based on machine learning-based modeling.")),
  tags$li(p("GENOXA is no-code interface which allows users without programming experience to easily analyze bulk RNA-Seq data and obtain results quickly."))
  ),
  h3(tags$strong("Key Features")),
  tags$ul(
   tags$li(p("Traditional bioinformatics results for bulk RNA-Seq data using limma and edgeR.")),
   tags$li(p("Interactive and high-quality visualizations, including options to download plots in .png or .pdf formats with customizable sizes and resolutions.")),
   tags$li(p("Multiple color schemes and customization options for plots to enhance data presentation.")),
   tags$li(p("Comprehensive Gene Ontology (GO) over-representation analysis, including GO plots, GO network plots, KEGG plots, and other valuable visualizations.")),
   tags$li(p("Detailed reports of analysis parameters and methods, available for download in .html, .pdf, and .doc formats."))
  ),
  h3(tags$strong(actionLink("ABTR", "Results Provided", class = "collapsible-header"))),
  div(id = "ABTRR", style = "display: none;",  # Initially hidden
  h4(tags$strong("1. User Data Insights")),
  tags$ul(
  tags$li(p(tags$b("Raw Count Data Table:"), "Easily view and explore your raw count data.")),
  tags$li(p(tags$b("Meta Data Overview:"), "Understand your sample characteristics at a glance.")),
  tags$li(p(tags$b("Clinical/Annotation data:"), "View and explore your uploaded Clinical/Annotation data.")),
  tags$li(p(tags$b("Raw Read Count Plot:"), "Visualize the distribution of your raw read counts."))
  ),
  h4(tags$strong("2. Quality Control and Data Normalization")),
  tags$ul(
  tags$li(p(tags$b("Filter Raw Counts Plot:"), "Assess the quality of your filtered data.")),
  tags$li(p(tags$b("Before and After Filtering Plots:"), "Visualize the impact of your filtering process.")),
  tags$li(p(tags$b("Data Normalization Plot:"), "Ensure your data is ready for analysis after data normalization step."))
  ),
  h4(tags$strong("3. Exploratory Data Analysis")),
  tags$ul(
  tags$li(p(tags$b("PCA Plot:"), "Identify patterns and relationships within your data.")),
  tags$li(p(tags$b("MDS Plot:"), "Explore multidimensional relationships effectively.")),
  tags$li(p(tags$b("Voom Plot:"), "Visualize the mean-variance relationship.")),
  tags$li(p(tags$b("Mean-Variance Trend Plot:"), "Gain insights into your data variability."))
  ),
  h4(tags$strong("4. Differential Expression (DE) Results")),
  tags$ul(
  tags$li(p(tags$b("DE Table:"), "Access detailed results of differential expression analysis.")),
  tags$li(p(tags$b("Volcano Plot:"), "Easily identify significant and not significant genes.")),
  tags$li(p(tags$b("DE Heat Map Plot:"), "Visualize expression patterns across samples."))
  ),
  h4(tags$strong("5. Gene Modeling using ML Results")),
  tags$ul(
  tags$li(p(tags$b("Performance Plot:"), "Evaluate the effectiveness of your model.")),
  tags$li(p(tags$b("Variable Importance Plot:"), "Discover key predictors in your data.")),
  tags$li(p(tags$b("ROC Plot:"), "Assess model performance through receiver operating characteristics."))
  ),
  h4(tags$strong("6. Over-Representation Analysis (ORA) Results")),
  tags$ul(
  tags$li(p(tags$b("Bar Plot & Dot Plot:"), "Visualize gene set enrichment results.")),
  tags$li(p(tags$b("Gene-Concept Network Plot:"), "Explore relationships between genes and biological concepts.")),
  tags$li(p(tags$b("Heat Plot & Enrichment Map Plot:"), "Gain insights into enriched pathways and processes.")),
  tags$li(p(tags$b("UpSet Plot:"), "Identify intersections in your data.")),
  tags$li(p(tags$b("KEGG Results & Plot:"), "Explore pathway-level insights."))
  ),
  h4(tags$strong("7. Comprehensive Report Generation")),
  tags$ul(
  tags$li(p("Generate a detailed report documenting the analysis parameters and tools used in your data analysis."))
  )
  ),
  h3(tags$strong("Benefits")),
  tags$ul(
  tags$li(p("No-code software that significantly reduces the time and complexity of RNA-Seq data analysis with results using ML algorithms for DE genes.")),
  tags$li(p("Valuable for researchers, students, and other users of bulk RNA-Seq data seeking deeper insights through ML algorithms."))
  ),
  h3(tags$strong("GENOXA workflow")),
  tags$br(),
  tags$img(src="./wrk.jpg", width="800"),
  tags$br(),
  h3(tags$b("YouTube tutorial for GENOXA App:")),
  tags$ul(
  tags$li(p("Please click the following YouTube link to see how to use the GENOXA app.")),
  tags$li(tags$b(p(a(href="https://youtu.be/6nPyV5qtOUc", target="_blank", "YouTube tutorial of GENOXA App"))))
  ),
  h3(tags$strong("GENOXA tool developer")),
  tags$ul(
  tags$li(p(tags$b("Sagar Patel, Ph.D."),"|", (tags$b(a(href="https://www.linkedin.com/in/sgr308", target="_blank", "LinkedIn"))))),
  tags$li(p(tags$b("Email: sgr.bnf(at)gmail.com")))
  ),
  h3(tags$strong("News and Updates")),
  tags$ul(
  tags$li(tags$b("On May 8th, 2025, GENOXA Version 0.2 was released, featuring integration with few TCGA datasets from the GDC portal as part of the GDC Tool Challenge."))
         )
       )
     )
   ),
   
# End of first tab panel -Home
#Second tab item
tabItem(tabName = "manual",
 tabBox(title = "User manual for GENOXA", width = 12,
 tabPanel(title = "Select Data",
  wellPanel(
  h4(tags$b("TCGA data from GDC portal:")),
  tags$ul(
    tags$li(p(tags$b("The TCGA (The Cancer Genome Atlas) is a landmark cancer genomics program that provides comprehensive genomic data from thousands of cancer patients across various cancer types. Bulk RNA-Seq data for the TCGA project were downloaded from the GDC portal. The dataset falls under the Transcriptome Profiling category, utilizing the RNA-Seq experimental strategy and the STAR - Counts workflow. It includes Primary Tumor and Solid Tissue Normal sample types, and is publicly accessible with open access. Below are brief descriptions of the datasets available for you to explore. Each dataset includes molecular, clinical, and expression data for different cancer types, allowing for in-depth research and analysis."))),
  ),
  h4(tags$b("1. TCGA-PRAD (Prostate Adenocarcinoma) - Prostate Cancer")),
  tags$ul(
  tags$li(p(tags$b("The TCGA-PRAD dataset contains data for prostate adenocarcinoma, a common and aggressive form of prostate cancer in men. The dataset includes transcriptomic and clinical data that can help identify molecular markers associated with prostate cancer progression and therapeutic response.")))
  ),
  h4(tags$b("2. TCGA-BRCA (Breast Invasive Carcinoma) - Breast Cancer")),
  tags$ul(
    tags$li(p(tags$b("This dataset features invasive breast carcinoma data, one of the most frequently diagnosed cancers in women. The TCGA-BRCA dataset provides insights into the genetic and molecular alterations that drive breast cancer, with applications in targeted therapy and prognostic studies.")))
  ),
  h4(tags$b("3. TCGA-LUAD (Lung Adenocarcinoma) - Lung Cancer")),
  tags$ul(
    tags$li(p(tags$b("The TCGA-LUAD dataset focuses on lung adenocarcinoma, a subtype of non-small cell lung cancer (NSCLC). The dataset is vital for understanding the transcriptomic changes that lead to the development of lung cancer and for identifying novel therapeutic targets.")))
  ),
  h4(tags$b("4. TCGA-KIRC (Kidney Renal Clear Cell Carcinoma) - Kidney Cancer")),
  tags$ul(
    tags$li(p(tags$b("The TCGA-KIRC dataset provides data on clear cell renal carcinoma, the most common subtype of kidney cancer. This dataset helps to investigate the molecular mechanisms underlying kidney cancer development and response to therapies.")))
  ),
  h4(tags$b("5. TCGA-COAD (Colon Adenocarcinoma) - Colon Cancer")),
  tags$ul(
    tags$li(p(tags$b("TCGA-COAD contains data from colon adenocarcinoma, a common and often aggressive cancer of the colon. The dataset is useful for identifying key genetic alterations that promote tumorigenesis in the colon and may guide clinical decision-making.")))
  ),
  h4(tags$b("6. TCGA-STAD (Stomach Adenocarcinoma) - Stomach Cancer")),
  tags$ul(
    tags$li(p(tags$b("The TCGA-STAD dataset focuses on stomach adenocarcinoma, a common form of gastric cancer. The data is crucial for exploring genetic drivers of tumor growth and metastasis in the stomach, as well as for identifying potential biomarkers for early detection and treatment.")))
  ),
  h4(tags$b("7. TCGA-BLCA (Bladder Urothelial Carcinoma) - Bladder Cancer")),
  tags$ul(
    tags$li(p(tags$b("TCGA-BLCA includes data for urothelial carcinoma of the bladder, a major cause of morbidity and mortality. The dataset helps researchers understand the genetic landscape of bladder cancer and identify molecular targets for novel therapies.")))
  ),
  h4(tags$b("8. TCGA-LIHC (Liver Hepatocellular Carcinoma) - Liver Cancer")),
  tags$ul(
    tags$li(p(tags$b("The TCGA-LIHC dataset covers liver hepatocellular carcinoma, a leading cause of cancer-related deaths worldwide. The dataset provides insights into the molecular mechanisms driving liver cancer and the potential for personalized treatment approaches.")))
  )
)),
 tabPanel(title = "Data Preprocessing",
 wellPanel(
  h4(tags$b("Data Preprocessing parameters:")),
  tags$ol(
  tags$li(p(tags$b("Filter raw reads by rowSum:"), "The rowSum option in limma allows you to filter genes based on the total number of reads they have across all samples. Filtering out such genes before performing differential expression analysis can improve the quality of your results and reduce computational burden.")),  
  tags$li(p(tags$b("Filter raw reads by FilterByExpr:"), "The FilterByExpr option in limma is designed to help with this by filtering genes based on their expression levels across samples. This method helps in removing genes that are not expressed enough to be of interest, thereby improving the robustness of your differential expression analysis.")),
  tags$li(p(tags$b("Select Normalization method for calcNormFactors:"), "The calcNormFactors function in edgeR helps to normalize the raw counts, so they are comparable across different samples. This function provides several normalization methods to choose from.")),
  h5(tags$b("You can select any of the following normalization method for the calcNormFactors option that best fits your analysis needs:")),
  tags$ul(
  tags$li(p(tags$b("TMM (Trimmed Mean of M-values):"), "TMM normalization is designed to account for differences in library sizes and compositional biases between samples. It adjusts for both systematic and sample-specific effects by comparing the relative abundance of genes between samples. TMM is particularly useful when dealing with RNA-seq data that has significant differences in library sizes or composition between samples.")),
  tags$li(p(tags$b("RLE (Relative Log Expression):"), "RLE normalization adjusts for library size by scaling each sample to have the same median relative log expression. This method assumes that most genes are not differentially expressed and uses the median expression of all genes as a scaling factor. RLE is often used when the assumption that most genes are not differentially expressed holds true and is suitable for data with moderate to high sequencing depth.")),
  tags$li(p(tags$b("Upper Quartile:"), "This method normalizes data by scaling each sample based on the upper quartile of gene counts. It is less sensitive to extreme values compared to methods like TMM and RLE. The Upper Quartile method is useful in situations where you want to minimize the influence of highly expressed genes on the normalization process.")),
  tags$li(p(tags$b("None:"), "Selecting *none* means no normalization will be applied. This option is generally not recommended unless you have already normalized your data through other means or you want to perform analysis without any normalization adjustments. Use this option only if you are certain that no additional normalization is needed or if you have pre-processed the data with other normalization techniques."))
  )
  )
  )
 ),
 tabPanel(title = "Differential expression analysis",
 wellPanel(
 h4(tags$b("Differential expression (DE) analysis parameters:")),
 tags$ol(
 tags$li(p(tags$b("Select reference group for DE:"), "Choose the condition that will serve as the reference group for differential expression analysis. The reference group should represent a baseline or control condition, against which all other conditions will be compared. This reference group will be set using the relevel() function from the edgeR package, ensuring meaningful comparisons across your experimental conditions. For example, if you select *Control* as the reference group, it will be the baseline condition used to compare with other experimental conditions in your analysis.")),  
 tags$li(p(tags$b("Cutoff value for p-values and adjusted p-values:"), "Set the threshold for significance when evaluating p-values and adjusted p-values in your differential expression analysis. This cutoff value determines the maximum p-value or adjusted p-value considered significant.")),
 tags$ul(
 tags$li(p(tags$b("p-value cutoff:"), "Defines the threshold for the raw p-values before correction for multiple testing.")),
 tags$li(p(tags$b("Adjusted p-value cutoff:"), "Sets the threshold for the p-values after applying corrections for multiple comparisons, such as the False Discovery Rate (FDR)."))
 ),
 tags$li(p(tags$b("Minimum absolute log2-fold change (log2FC) to use as threshold for differential expression:"), "Set the minimum absolute log2-fold change (log2FC) to define the threshold for differential expression. This value determines the smallest magnitude of expression change required for a gene to be considered significantly different between conditions. For example, a threshold of 1 will include genes with log2FC values of ≥1 or ≤-1, highlighting those with substantial changes in expression.")),
 tags$li(p(tags$b("Select adjustment method:"), "Choose the method for adjusting p-values to account for multiple comparisons. This adjustment controls the false discovery rate and reduces the likelihood of type I errors.")),
 h5(tags$b("Choose the adjustment method that best fits your analysis needs from the following options:")),
 tags$ul(
 tags$li(p(tags$b("None:"), "No adjustment applied.")),
 tags$li(p(tags$b("Bonferroni:"), "A stringent method that adjusts p-values based on the number of tests.")),
 tags$li(p(tags$b("Holm:"), "A stepwise procedure that is less conservative than Bonferroni.")),
 tags$li(p(tags$b("BH (Benjamini-Hochberg):"), "Controls the false discovery rate, offering a balance between type I and type II errors.")),
 tags$li(p(tags$b("FDR (False Discovery Rate):"), "Similar to BH, it aims to control the proportion of false discoveries among the rejected hypotheses.")),
 tags$li(p(tags$b("BY (Benjamini-Yekutieli):"), "Adjusts for dependencies between tests, providing a more conservative control of the false discovery rate."))
 )
 )
 )
 ),
 tabPanel(title = "GO over-representation",
 wellPanel(
  h4(tags$b("Gene Ontology (GO) over-representation parameters:")),
  h4(actionLink("GOP1", "1. Select Organism Database:", class = "collapsible-header")),
  div(id = "GOP11", style = "display: none;",  # Initially hidden
  h5(tags$b("Choose the appropriate organism database for your analysis to ensure accurate gene annotations and biological insights. The available options are:")),
  tags$ul(
  tags$li(p(tags$b("Human (org.Hs.eg.db):"), "Database for human genes.")),
  tags$li(p(tags$b("Mouse (org.Mm.eg.db):"), "Database for mouse genes.")),
  tags$li(p(tags$b("Rat (org.Rn.eg.db):"), "Database for rat genes.")),
  tags$li(p(tags$b("Yeast (org.Sc.sgd.db):"), "Database for yeast genes.")),
  tags$li(p(tags$b("Fly (org.Dm.eg.db):"), "Database for fruit fly genes.")),
  tags$li(p(tags$b("Arabidopsis (org.At.tair.db):"), "Database for Arabidopsis genes.")),
  tags$li(p(tags$b("Zebrafish (org.Dr.eg.db):"), "Database for zebrafish genes.")),
  tags$li(p(tags$b("Bovine (org.Bt.eg.db):"), "Database for bovine genes.")),
  tags$li(p(tags$b("Worm (org.Ce.eg.db):"), "Database for worm genes.")),
  tags$li(p(tags$b("Chicken (org.Gg.eg.db):"), "Database for chicken genes.")),
  tags$li(p(tags$b("Canine (org.Cf.eg.db):"), "Database for dog genes.")),
  tags$li(p(tags$b("Pig (org.Ss.eg.db):"), "Database for pig genes.")),
  tags$li(p(tags$b("E coli strain K12 (org.EcK12.eg.db):"), "Database for E. coli K12 genes.")),
  tags$li(p(tags$b("Xenopus (org.Xl.eg.db):"), "Database for Xenopus genes.")),
  tags$li(p(tags$b("Chimp (org.Pt.eg.db):"), "Database for chimpanzee genes.")),
  tags$li(p(tags$b("Anopheles (org.Ag.eg.db):"), "Database for Anopheles mosquitoes.")),
  tags$li(p(tags$b("E coli strain Sakai (org.EcSakai.eg.db):"), "Database for E. coli Sakai genes"))
  )),
  h4(actionLink("GOP2", "2. Select Key Type:", class = "collapsible-header")),
  div(id = "GOP22", style = "display: none;",  # Initially hidden
  h5(tags$b("Select the identifier type for your gene data. This key type will dictate how gene identifiers are represented and interpreted in your analysis. Ensure that your gene names align with one of the following options and choose the appropriate one. Available options include:")),
  tags$ul(
  tags$li(p(tags$b("ACCNUM:"), "Accession number (e.g., NM_007294) - an accession number from GenBank.")),
  tags$li(p(tags$b("ALIAS:"), "Gene alias names (e.g., BRCA1) - a common name or alias for the gene.")),
  tags$li(p(tags$b("ENSEMBL:"), "Ensembl gene IDs (e.g., ENSG00000012048).")),
  tags$li(p(tags$b("ENSEMBLPROT:"), "Ensembl protein IDs (e.g., ENSP00000354819).")),
  tags$li(p(tags$b("ENSEMBLTRANS:"), "Ensembl transcript IDs (e.g., ENST00000357654).")),
  tags$li(p(tags$b("ENTREZID:"), "Entrez Gene IDs (e.g., 672).")),
  tags$li(p(tags$b("GENENAME:"), "Gene names (e.g., BRCA1).")),
  tags$li(p(tags$b("REFSEQ:"), "RefSeq gene IDs (e.g., NM_007294).")),
  tags$li(p(tags$b("SYMBOL:"), "Gene symbols (e.g., BRCA1) - official gene symbol.")),
  tags$li(p(tags$b("UNIGENE:"), "UniGene cluster IDs (e.g., Hs.5074).")),
  tags$li(p(tags$b("UNIPROT:"), "UniProt protein IDs (e.g., P38398)."))
  )),
  h4(actionLink("GOP3", "3. Select subontology:", class = "collapsible-header")),
  div(id = "GOP33", style = "display: none;",  # Initially hidden
  h5(tags$b("Choose the specific subontology for Gene Ontology (GO) analysis to focus on particular aspects of gene functions. The available options are:")),
  tags$ul(
  tags$li(p(tags$b("BP (Biological Process):"), "Encompasses processes and series of events at the cellular or organismal level.")),
  tags$li(p(tags$b("MF (Molecular Function):"), "Describes the activities performed by individual gene products at the molecular level.")),
  tags$li(p(tags$b("CC (Cellular Component):"), "Relates to the locations within the cell where gene products are active.")),
  tags$li(p(tags$b("ALL:"), "Includes all three subontologies (Biological Process, Molecular Function, and Cellular Component)"))
  )),
  h4(actionLink("GOP4", "4. Minimum Gene Set Size:", class = "collapsible-header")),
  div(id = "GOP44", style = "display: none;",  # Initially hidden
  tags$ul(tags$li(p("Specify the minimum number of genes required in a gene set for it to be considered in the analysis. This threshold helps filter out gene sets with too few genes, ensuring that only those with a sufficient number of genes are included. For example, setting a minimum gene set size of 10 means that only gene sets with 10 or more genes will be analyzed."))
  )),
  h4(actionLink("GOP5", "5. Maximum Gene Set Size:", class = "collapsible-header")),
  div(id = "GOP55", style = "display: none;",  # Initially hidden
  tags$ul(tags$li(p("Specify the maximum number of genes allowed in a gene set for it to be included in the analysis. This threshold helps filter out gene sets that are too large, focusing the analysis on more manageable and relevant gene sets. For example, setting a maximum gene set size of 500 means that only gene sets with 500 or fewer genes will be considered in the analysis."))
  )),
  h4(actionLink("GOP6", "6. P-Value Cutoff:", class = "collapsible-header")),
  div(id = "GOP66", style = "display: none;",  # Initially hidden
  tags$ul(tags$li(p("Specify the threshold p-value for determining statistical significance in Gene Ontology (GO) over-representation analysis. This cutoff defines the maximum p-value considered significant for identifying enriched GO terms. For example, setting a cutoff of 0.05 means that only GO terms with p-values less than or equal to 0.05 will be considered significantly enriched. Adjust this value to balance the sensitivity and specificity of detecting meaningful biological insights."))
  )), 
  h4(actionLink("GOP7", "7. q-Value Cutoff:", class = "collapsible-header")),
  div(id = "GOP77", style = "display: none;",  # Initially hidden
  tags$ul(tags$li(p("Specify the threshold q-value for determining statistical significance in Gene Ontology (GO) over-representation analysis. The q-value represents the false discovery rate (FDR) and adjusts p-values for multiple comparisons. Setting a cutoff, such as 0.05, means that only GO terms with q-values less than or equal to 0.05 will be considered significantly enriched. Adjust this value to manage the trade-off between detecting true enrichments and controlling the rate of false discoveries."))
  )),
  h4(actionLink("GOP8", "8. Adjust Method:", class = "collapsible-header")),
  div(id = "GOP88", style = "display: none;",  # Initially hidden
  h5("Choose the method for adjusting p-values in Gene Ontology (GO) over-representation analysis to account for multiple comparisons. This adjustment helps control for false positives and manage the false discovery rate. Available options include:"),
  tags$ul(
  tags$li(p(tags$b("None:"), "No adjustment applied.")),
  tags$li(p(tags$b("Holm:"), "A stepwise adjustment method that is less conservative than Bonferroni.")),
  tags$li(p(tags$b("Hochberg:"), "A step-up procedure that controls the family-wise error rate with more power than Holm.")),
  tags$li(p(tags$b("Hommel:"), "A step-down method that controls the family-wise error rate, offering a less conservative alternative to Bonferroni.")),
  tags$li(p(tags$b("Bonferroni:"), "A stringent adjustment method that divides p-values by the number of tests, reducing the chance of false positives.")),
  tags$li(p(tags$b("BH (Benjamini-Hochberg):"), "Controls the false discovery rate, balancing between type I and type II errors.")),
  tags$li(p(tags$b("FDR (False Discovery Rate):"), "Similar to BH, it aims to control the proportion of false discoveries among the rejected hypotheses.")),
  tags$li(p(tags$b("BY (Benjamini-Yekutieli):"), "Adjusts for dependencies between tests, providing a more conservative control of the false discovery rate."))
  ))
  )
  ),
tabPanel(title = "Gene modeling",
 wellPanel(
 h4(tags$b("Gene modeling parameters:")),
 h4(tags$b("*Select parameters for split data:")),
 h4(actionLink("GOP11", "1. Enter The size of the training data set:", class = "collapsible-header")),
 div(id = "GOP1111", style = "display: none;",  # Initially hidden
     tags$ul(tags$li(p("This option allows you to specify the number of samples used for training your machine learning model. It is essential for setting the model’s learning capacity and overall performance. Configuring this setting correctly ensures effective training and accurate predictions.")))
 ),
 h4(actionLink("GOP12", "2. Seed:", class = "collapsible-header")),
 div(id = "GOP1212", style = "display: none;",  # Initially hidden
     tags$ul(tags$li(p("This option allows you to set a starting point for random number generation used during model training or data processing. By specifying a seed value, you ensure that the random processes are reproducible, allowing for consistent results across different runs and facilitating easier debugging and comparison.")))
 ),
 h4(tags$b("*Select parameters for train models:")),
 h4(actionLink("GOP9", "3. Select resampling method to use:", class = "collapsible-header")),
 div(id = "GOP99", style = "display: none;",  # Initially hidden
 h5(tags$b("Choose the resampling method to use for training your model from the dropdown list. This option determines how your training data is split or resampled to validate the model’s performance.")),
 h5(tags$b("You can select any of the following resampling method that best fits your analysis needs:")),
 tags$ul(
  tags$li(p(tags$b("None (none):"), "No resampling method is applied. This option uses the entire dataset for training and evaluation without any form of cross-validation or bootstrapping.")),
  tags$li(p(tags$b("Bootstrapping (boot):"), "Resamples the dataset with replacement to create multiple training sets. This method helps estimate the variability of the model’s performance and provides insight into the stability of the model.")),
  tags$li(p(tags$b("Bootstrap with 0.632 correction (boot632):"), "An extension of bootstrapping that includes a correction factor to account for the bias in the bootstrap estimate. This method helps to provide a more accurate estimate of model performance.")),
  tags$li(p(tags$b("Optimism Bootstrapping (optimism_boot):"), "This method applies bootstrapping and adjusts for the *optimism* in performance estimates by comparing the bootstrap estimates to the original model. It helps in assessing how optimistic the initial performance estimate might be.")),
  tags$li(p(tags$b("Bootstrapping All (boot_all):"), "Uses bootstrapping on all available data, rather than a subset. This method provides a thorough estimate of model performance by resampling the entire dataset multiple times.")),
  tags$li(p(tags$b("Cross-Validation (cv):"), "Divides the dataset into k folds and trains the model k times, each time using a different fold as the validation set and the remaining folds as the training set. This method provides a robust estimate of model performance by averaging results over multiple folds.")),
  tags$li(p(tags$b("Repeated Cross-Validation (repeatedcv):"), "Repeats the k-fold cross-validation process multiple times with different random splits. This method enhances the robustness of performance estimates by reducing variance and providing a more comprehensive evaluation.")),
  tags$li(p(tags$b("Number of Complete Sets of Folds:"), "If you select *repeatedcv* as the resampling method, you will need to specify the number of complete sets of folds to compute. Enter a positive integer value to define how many times the cross-validation process should be repeated. This helps in improving the robustness and reliability of the model evaluation.")),
  tags$li(p(tags$b("Leave-One-Out Cross-Validation (LOOCV):"), "A special case of cross-validation where each data point is used once as a validation set while the remaining data points form the training set. This method is computationally intensive but provides a nearly unbiased estimate of model performance.")),
  tags$li(p(tags$b("Leave-Group-Out Cross-Validation (LGOCV):"), "Similar to LOOCV but used when the data is grouped. One or more groups are left out for validation, and the remaining groups are used for training. This method is useful when dealing with grouped data.")),
  tags$li(p(tags$b("Out-of-Bag (oob):"), "An estimation method used with ensemble methods like random forests. Each tree is trained on a bootstrap sample, and out-of-bag samples are used to estimate performance without needing a separate validation set.")),
  tags$li(p(tags$b("Timeslice (timeslice):"), "Suitable for time series data, this method involves splitting the data into contiguous time periods for training and validation. It respects the temporal order of the data, ensuring that training always occurs before validation.")),
  tags$li(p(tags$b("Adaptive Cross-Validation (adaptive_cv):"), "A variation of cross-validation that adapts the size of the validation set based on the performance of the model. This method can be more flexible and responsive to model performance changes.")),
  tags$li(p(tags$b("Adaptive Bootstrapping (adaptive_boot):"), "An adaptive version of bootstrapping that adjusts the resampling process based on the model’s performance. This method aims to improve the estimation accuracy by dynamically changing the resampling strategy.")),
  tags$li(p(tags$b("Adaptive Leave-Group-Out Cross-Validation (adaptive_LGOCV):"), "An adaptive version of LGOCV where the group size or the number of groups left out for validation adjusts based on model performance. This approach provides a more tailored evaluation process for grouped data."))
  )
 ),
 h4(actionLink("GOP13", "4. Number of resampling iterations:", class = "collapsible-header")),
 div(id = "GOP1313", style = "display: none;",  # Initially hidden
     tags$ul(tags$li(p("This option specifies how many times the app will train and validate your machine learning model on different data subsets using cross-validation. More iterations provide a more accurate estimate of your model's performance by averaging results across various splits, offering a clearer view of its generalization ability. However, this also increases computation time.")))
 ),
 h4(actionLink("GOP14", "5. Random number seed:", class = "collapsible-header")),
 div(id = "GOP1414", style = "display: none;",  # Initially hidden
     tags$ul(tags$li(p("This option sets the initial value for random number generation, ensuring that the results are reproducible. By using the same seed, you can consistently recreate the same random data splits and model outcomes, which is useful for comparing results and debugging. Different seeds will generate different results, so changing the seed will affect the random aspects of your model training and validation process.")))
 ),
 h4(actionLink("GOP10", "6. Select ML algorithms for prediction:", class = "collapsible-header")),
 div(id = "GOP1010", style = "display: none;",  # Initially hidden
  h5(tags$b("This option lets you choose which machine learning algorithms to use for making predictions. By selecting different algorithms, you can apply various models to your data, each with its own strengths and weaknesses. Experimenting with different algorithms helps you find the best approach for your specific problem, potentially improving prediction accuracy and performance.")),
  h5(tags$b(a(href="https://topepo.github.io/caret/available-models.html", target="_blank", "Click here for detailed guide on machine learning algorithms."))),
  h5(tags$b("Please select minimum two and maximum four ML algorithms for better data visualization.")),
  tags$ul(
   tags$li(p(tags$b("Bagged CART (treebag):"), "Uses bagging (Bootstrap Aggregating) with Classification and Regression Trees (CART). It builds multiple trees on random samples of the data and combines their results to improve accuracy and reduce overfitting.")),
   tags$li(p(tags$b("Bagged FDA using gCV Pruning (bagFDAGCV):"), "Applies bagging to Flexible Discriminant Analysis (FDA) with generalized cross-validation (gCV) pruning to enhance model stability and performance.")),
   tags$li(p(tags$b("Bagged Flexible Discriminant Analysis (bagFDA):"), "Uses bagging with Flexible Discriminant Analysis, combining multiple FDA models built on random data subsets to improve prediction robustness.")),
   tags$li(p(tags$b("Boosted Generalized Linear Model (glmboost):"), "Enhances a Generalized Linear Model (GLM) using boosting, which sequentially adjusts the model to correct errors from previous iterations, improving predictive accuracy.")),
   tags$li(p(tags$b("C5.0 (C5.0):"), "An advanced version of the C4.5 decision tree algorithm, known for creating accurate and efficient classification models using a tree-like structure.")),
   tags$li(p(tags$b("CART (rpart):"), "Classification and Regression Trees (CART) build decision trees for classification and regression tasks, providing a straightforward model that splits data based on feature values.")),
   tags$li(p(tags$b("Conditional Inference Random Forest (cforest):"), "A Random Forest algorithm that uses conditional inference trees, offering robust and unbiased predictions by avoiding overfitting and handling interactions between features effectively.")),
   tags$li(p(tags$b("eXtreme Gradient Boosting (xgbLinear):"), "A variant of eXtreme Gradient Boosting optimized for linear models, enhancing performance through gradient boosting techniques that focus on improving model accuracy.")),
   tags$li(p(tags$b("eXtreme Gradient Boosting (xgbTree):"), "Uses gradient boosting with tree-based models to improve predictive accuracy, incorporating advanced techniques to handle complex data relationships.")),
   tags$li(p(tags$b("Flexible Discriminant Analysis (fda):"), "A discriminant analysis method that can adapt to different data structures and relationships, useful for classification tasks with flexible model assumptions.")),
   tags$li(p(tags$b("Generalized Linear Model (glm):"), "A versatile model that extends linear regression to handle various types of response variables, including binary, count, and continuous outcomes.")),
   tags$li(p(tags$b("glmnet (glmnet):"), "A regularized version of the Generalized Linear Model that applies Lasso or Ridge regression to improve model performance and handle high-dimensional data.")),
   tags$li(p(tags$b("Parallel Random Forest (parRF):"), "A Random Forest algorithm designed to leverage parallel processing for faster computation and improved efficiency, particularly useful with large datasets.")),
   tags$li(p(tags$b("Partial Least Squares (kernelpls):"), "A technique that combines features of principal component analysis and regression to handle multicollinearity and high-dimensional data.")),
   tags$li(p(tags$b("Random Forest (rf):"), "An ensemble learning method that builds multiple decision trees and aggregates their predictions, known for its robustness and ability to handle various data types.")),
   tags$li(p(tags$b("Regularized Random Forest (RRF):"), "Enhances Random Forest with regularization techniques to improve model performance and reduce overfitting.")),
   tags$li(p(tags$b("Regularized Random Forest (RRFglobal):"), "A variant of RRF that applies global regularization, aiming to enhance the model's predictive accuracy and generalization.")),
   tags$li(p(tags$b("Rotation Forest (rotationForest):"), "An ensemble method that uses rotated feature spaces to build a set of decision trees, enhancing classification accuracy through feature transformation.")),
   tags$li(p(tags$b("Rotation Forest (rotationForestCp):"), "A variant of Rotation Forest that applies complexity pruning to the individual trees, aiming to improve model performance and interpretability.")),
   tags$li(p(tags$b("Single C5.0 Ruleset (C5.0Rules):"), "Generates a ruleset from the C5.0 algorithm, providing a simplified model with a set of easy-to-interpret classification rules.")),
   tags$li(p(tags$b("Single C5.0 Tree (C5.0Tree):"), "Constructs a single decision tree using the C5.0 algorithm, offering a clear and interpretable model for classification tasks.")),
   tags$li(p(tags$b("Sparse Distance Weighted Discrimination (sdwd):"), "A classification method that uses distance-weighted discrimination with sparsity constraints, effective for handling high-dimensional and sparse data."))
   )),
 h3(tags$b("Important Disclaimer on Machine Learning Results")),
 tags$ul(
   tags$li(p("Please note that the machine learning results generated by our app are based on advanced ML models. We emphasize that these results should be validated with additional methods before drawing any conclusions. Our company does not take responsibility for any claims or decisions made based on the ML results retrieved using our app."))
 )
 )
) # END of ML tabpanel
) # END of manual tabbox      
), # END of second tabpanel

# Third tab content
tabItem(tabName = "results",
 tabsetPanel(
 tabPanel(title = "User Data Information",
 tabBox(title = "User Data Information", width = 12,
 tabPanel(title = "Raw Read Counts Data",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of rows in data:"),
 tags$b(textOutput("rc_genes"))
 )),
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of columns in data:"),
 tags$b(textOutput("rc_samples"))
 )
 )),
 DT::dataTableOutput("raw_counts")
 ),
 tabPanel(title = "Sample and Group Information",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of rows in data:"),
 tags$b(textOutput("samp_rows"))
 )),
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of columns in data:"),
 tags$b(textOutput("samp_columns"))
 )
 )),
 DT::dataTableOutput("sample_info")
 ),
 tabPanel(title = "Raw Read Counts Plot",
  fluidRow(
  column(
  width = 6,
  wellPanel(
  # Sample Selection Inputs
  fluidRow(
  column(6, numericInput("rcs1", "First sample to plot:", min = 1, max = 999999, value = 1)),
  column(6, numericInput("rcs2", "Last sample to plot:", min = 1, max = 999999, value = 10))
  ),
  tags$b("Note: If you want to view all samples in the plot, replace the value 10 with the total number of samples shown below."),
  tags$b(textOutput("srange"))
  ),
  wellPanel(
  # Download Options
  tags$b("Note: The following options apply only to the downloaded plot."),
  tags$hr(),
  fluidRow(
  column(6, numericInput("rc_ph", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
  column(6, numericInput("rc_pw", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
  column(6, numericInput("rc_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
  column(6, radioButtons("rc_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
  ),
  downloadButton("rc_plotd11", "Download Plot")
  )
  ),
  column(
  width = 6,
  wellPanel(
  # Plot Customization Inputs
  radioButtons("rc_plotpe1", "Choose Plot Type:",
   choices = c("Boxplot" = "boxplot", "Column" = "col", "Violin" = "violin"), 
   inline = TRUE),
  radioButtons("rc_dtr1", "Choose Data Transformation:",
  choices = c("Count" = "Count", "Log2 Count" = "log2_Count"), 
  inline = TRUE, selected = "log2_Count"),
  radioButtons("rc_plt1", "Select Color Palette for plots:",
   choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
   inline = TRUE),
  fluidRow(
   column(6,
   checkboxInput("show_x_labels1", "Show X Axis Label", value = TRUE),
   checkboxInput("show_y_labels1", "Show Y Axis Label", value = TRUE)
   ),
   column(6,
   checkboxInput("show_plot_title1", "Show Plot Title", value = TRUE),
   checkboxInput("show_legend1", "Show Legend", value = TRUE)
   )
  ),
 tags$br(),
 # Define your conditional panels with input fields
 fluidRow(
  column(6,
  conditionalPanel(
   condition = "input.show_x_labels1",
   textInput("x_axis_label1", "X Axis Label", value = "Sample Name")
   ),
  conditionalPanel(
   condition = "input.show_y_labels1",
   textInput("y_axis_label1", "Y Axis Label", value = "Gene Counts")
   )),
  column(6,
   conditionalPanel(
    condition = "input.show_plot_title1",
    textInput("plot_title1", "Plot Title", value = "Raw Gene Count Data Plot")
    ),
   conditionalPanel(
    condition = "input.show_legend1",
    textInput("legend_title1", "Legend Title", value = "Sample_Type")
    )
   )
  )))
  ),
  wellPanel(
  plotOutput("rc_plot11", height = "600px", width = "920px")
  )
 ),
 tabPanel(title = "Clinical/Annotation data",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of rows in data:"),
 tags$b(textOutput("cli_rows"))
 )),
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of columns in data:"),
 tags$b(textOutput("cli_cols"))
 )
 )),
 DT::dataTableOutput("cl_data_table")
 )
) 
), # End of First tabpanel "User Data"
    
# Data pre-processing tabpanel - Filt_plo1
tabPanel(title = "QC and Data Normalization",
 tabBox(title = "QC and Data Normalization", width = 12,
# First QC plot panel
 tabPanel(title = "Filtered Read Counts Plot",
  fluidRow(
  column(
  width = 6,
  wellPanel(
  tags$b("Total number of filtered rows in data:"),
  tags$b(textOutput("filt_1_genes"))
  ),
  wellPanel(
  # Sample Selection Inputs
  fluidRow(
  column(6, numericInput("filt1_1", "First sample to plot:", min = 1, max = 999999, value = 1)),
  column(6, numericInput("filt1_2", "Last sample to plot:", min = 1, max = 999999, value = 10))
  ),
  tags$b("Note: If you want to view all samples in the plot, replace the value 10 with the total number of samples shown below."),
  tags$b(textOutput("filt_1_samples"))
  ),
  wellPanel(
  # Download Options
  tags$b("Note: The following options apply only to the downloaded plot."),
  tags$hr(),
  fluidRow(
  column(6, numericInput("filt1h_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
  column(6, numericInput("filt1w_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
  column(6, numericInput("filt1_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
  column(6, radioButtons("filt1_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
  ),
  downloadButton("filt1_plotd", "Download Plot")
  )
  ),
  column(
  width = 6,
  wellPanel(
  # Plot Customization Inputs
  radioButtons("filt1_plotpe", "Choose Plot Type:",
   choices = c("Boxplot" = "boxplot", "Column" = "col", "Violin" = "violin"), 
               inline = TRUE),
  radioButtons("filt1_dtr", "Choose Data Transformation:",
   choices = c("Count" = "Count", "Log2 Count" = "log2_Count"), 
               inline = TRUE, selected = "log2_Count"),
  radioButtons("filt1_plt", "Select Color Palette for plots:",
   choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
               inline = TRUE),
  fluidRow(
  column(6,
  checkboxInput("filt1_show_x_labels", "Show X Axis Label", value = TRUE),
  checkboxInput("filt1_show_y_labels", "Show Y Axis Label", value = TRUE)
  ),
  column(6,
  checkboxInput("filt1_show_plot_title", "Show Plot Title", value = TRUE),
  checkboxInput("filt1_show_legend", "Show Legend", value = TRUE)
  )
  ),
  tags$br(),
  # Define your conditional panels with input fields
  fluidRow(
  column(6,
  conditionalPanel(
    condition = "input.filt1_show_x_labels",
    textInput("filt1_x_axis_label", "X Axis Label", value = "Sample Name")
  ),
  conditionalPanel(
    condition = "input.filt1_show_y_labels",
    textInput("filt1_y_axis_label", "Y Axis Label", value = "Gene Counts")
  )),
  column(6,
  conditionalPanel(
    condition = "input.filt1_show_plot_title",
    textInput("filt1_plot_title", "Plot Title", value = "Filtered read counts plot")
    ),
  conditionalPanel(
    condition = "input.filt1_show_legend",
    textInput("filt1_legend_title", "Legend Title", value = "Sample_Type")
    )
   )
 )))
 ),
 wellPanel(
 plotOutput("filt1_plot", height = "600px", width = "920px")
 )
),
# Second QC plot panel
tabPanel(title = "FilterByExpr Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Total number of filtered rows (after FilterByExpr) in the data:"),
 tags$b(textOutput("filt2_genes"))
 ),
 wellPanel(
 # Sample Selection Inputs
 fluidRow(
 column(6, numericInput("filt2_1", "First sample to plot:", min = 1, max = 999999, value = 1)),
 column(6, numericInput("filt2_2", "Last sample to plot:", min = 1, max = 999999, value = 10))
 ),
 tags$b("Note: If you want to view all samples in the plot, replace the value 10 with the total number of samples shown below."),
 tags$b(textOutput("filt2_samples"))
 ),
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("filt2h_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("filt2w_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("filt2_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("filt2_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("filt2_plotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 # Plot Customization Inputs
 radioButtons("filt2_plotpe", "Choose Plot Type:",
  choices = c("Boxplot" = "boxplot", "Column" = "col", "Violin" = "violin"), 
  inline = TRUE),
 radioButtons("filt2_dtr", "Choose Data Transformation:",
  choices = c("Count" = "Count", "Log2 Count" = "log2_Count"), 
  inline = TRUE, selected = "log2_Count"),
 radioButtons("filt2_plt", "Select Color Palette for plots:",
  choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
  inline = TRUE),
 fluidRow(
  column(6,
  checkboxInput("filt2_show_x_labels", "Show X Axis Label", value = TRUE),
  checkboxInput("filt2_show_y_labels", "Show Y Axis Label", value = TRUE)
  ),
 column(6,
 checkboxInput("filt2_show_plot_title", "Show Plot Title", value = TRUE),
 checkboxInput("filt2_show_legend", "Show Legend", value = TRUE)
  )
 ),
 tags$br(),
# Define your conditional panels with input fields
 fluidRow(
  column(6,
  conditionalPanel(
   condition = "input.filt2_show_x_labels",
    textInput("filt2_x_axis_label", "X Axis Label", value = "Sample Name")
  ),
  conditionalPanel(
   condition = "input.filt2_show_y_labels",
    textInput("filt2_y_axis_label", "Y Axis Label", value = "Gene Counts")
  )),
  column(6,
  conditionalPanel(
   condition = "input.filt2_show_plot_title",
   textInput("filt2_plot_title", "Plot Title", value = "FilterByExpr plot")
   ),
  conditionalPanel(
   condition = "input.filt2_show_legend",
   textInput("filt2_legend_title", "Legend Title", value = "Sample_Type")
    )
   )
  )))
 ),
 wellPanel(
 plotOutput("filt2_plot", height = "600px", width = "920px")
 )
),
# Before and after filter panel
tabPanel(title = "Before and after filtering data",
 fluidRow(
  column(
  width = 6,
  wellPanel(
  fluidRow(
  column(6, numericInput("cmp_1", "First sample to plot:", min = 1, max = 999999, value = 1)),
  column(6, numericInput("cmp_2", "Last sample to plot:", min = 1, max = 999999, value = 10))
  ),
  tags$b("If you want to see all samples in the plot, replace the value 10 with the total number of samples shown below:"),
  tags$b(textOutput("cmp_smp"))
  ),
  wellPanel(
  tags$b("Total number of original rows in the data (left plot):"),
  tags$b(textOutput("rc_genes2")),
  tags$br(),
  tags$b("Total number of filtered rows (after FilterByExpr) in the data (right plot):"),
  tags$b(textOutput("filt22_genes"))
  ),
  wellPanel(
  # Additional Download Options
  tags$b("Note: The following options apply only to the downloaded plot."),
  tags$hr(),
  fluidRow(
  column(6, numericInput("cmp_ph", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
  column(6, numericInput("cmp_pw", "Width of the plot (inch):", min = 1, max = 5000, value = 12))
  ),
  fluidRow(
  column(6, numericInput("cmp_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
  column(6, radioButtons("cmp_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
  ),
  downloadButton("cmp_plotd", "Download Plot")
  )),
  column(
  width = 6,
  wellPanel(
  tags$b("Plot options for Raw Counts Plot:"),
  tags$hr(),
  fluidRow(
  column(5,
  checkboxInput("cmp_show_x_labels", "Show X Axis Label", value = TRUE),
  checkboxInput("cmp_show_y_labels", "Show Y Axis Label", value = TRUE),
  checkboxInput("cmp_show_plot_title", "Show Plot Title", value = TRUE)),
  column(7,
  conditionalPanel(
   condition = "input.cmp_show_x_labels",
   textInput("cmp_x_axis_label", "X Axis Label", value = "Log-cpm")
   ),
  conditionalPanel(
   condition = "input.cmp_show_y_labels",
   textInput("cmp_y_axis_label", "Y Axis Label", value = "Density")
   ),
  conditionalPanel(
   condition = "input.cmp_show_plot_title",
   textInput("cmp_plot_title", "Plot Title", value = "A. Raw Counts Plot")
   )
  )
  )),
  wellPanel(
  tags$b("Plot options for Filtered Counts (FilterByExpr) Plot:"),
  tags$hr(),
  fluidRow(
  column(5,
  checkboxInput("cmp1_show_x_labels", "Show X Axis Label", value = TRUE),
  checkboxInput("cmp1_show_y_labels", "Show Y Axis Label", value = TRUE),
  checkboxInput("cmp1_show_plot_title", "Show Plot Title", value = TRUE)),
  column(7,
  conditionalPanel(
   condition = "input.cmp1_show_x_labels",
   textInput("cmp1_x_axis_label", "X Axis Label", value = "Log-cpm")
   ),
  conditionalPanel(
   condition = "input.cmp1_show_y_labels",
   textInput("cmp1_y_axis_label", "Y Axis Label", value = "Density")
   ),
  conditionalPanel(
   condition = "input.cmp1_show_plot_title",
   textInput("cmp1_plot_title", "Plot Title", value = "B. Filtered Counts (FilterByExpr) Plot")
   )
  )))  
 )),
 tags$br(),             
 wellPanel(
 plotOutput("cmp_plot", height = "600px", width = "920px")
 )
), # End of Before and after filter panel
# Data Normalization
# Second QC plot panel
tabPanel(title = "Data Normalization Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
   tags$b("Total number of rows (after Data Normalization) in the data:"),
   tags$b(textOutput("DNORM_genes"))
 ),
 wellPanel(
  fluidRow(
  column(6, numericInput("DNORM_1", "First sample to plot:", min = 1, max = 999999, value = 1)),
  column(6, numericInput("DNORM_2", "Last sample to plot:", min = 1, max = 999999, value = 10))
  ),
  tags$b("Note: If you want to view all samples in the plot, replace the value 10 with the total number of samples shown below."),
  tags$b(textOutput("DNORM_samples"))
  ),
 wellPanel(
  # Download Options
  tags$b("Note: The following options apply only to the downloaded plot."),
  tags$hr(),
  fluidRow(
   column(6, numericInput("DNORMh_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
   column(6, numericInput("DNORMw_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
   column(6, numericInput("DNORM_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
   column(6, radioButtons("DNORM_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
   ),
  downloadButton("DNORM_plotd", "Download Plot")
 )
 ),
  column(
  width = 6,
  wellPanel(
  # Plot Customization Inputs
  radioButtons("DNORM_plotpe", "Choose Plot Type:",
                choices = c("Boxplot" = "boxplot", "Column" = "col", "Violin" = "violin"), 
                inline = TRUE),
  radioButtons("DNORM_plt", "Select Color Palette for plots:",
                choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
                inline = TRUE),
  fluidRow(
   column(6,
   checkboxInput("DNORM_show_x_labels", "Show X Axis Label", value = TRUE),
   checkboxInput("DNORM_show_y_labels", "Show Y Axis Label", value = TRUE)
   ),
   column(6,
   checkboxInput("DNORM_show_plot_title", "Show Plot Title", value = TRUE),
   checkboxInput("DNORM_show_legend", "Show Legend", value = TRUE)
   )
  ),
  tags$br(),
  # Define your conditional panels with input fields
  fluidRow(
   column(6,
   conditionalPanel(
    condition = "input.DNORM_show_x_labels",
    textInput("DNORM_x_axis_label", "X Axis Label", value = "Sample Name")
   ),
   conditionalPanel(
    condition = "input.DNORM_show_y_labels",
    textInput("DNORM_y_axis_label", "Y Axis Label", value = "Gene Counts")
   )),
   column(6,
    conditionalPanel(
     condition = "input.DNORM_show_plot_title",
     textInput("DNORM_plot_title", "Plot Title", value = "Data Normalization plot")
     ),
    conditionalPanel(
     condition = "input.DNORM_show_legend",
     textInput("DNORM_legend_title", "Legend Title", value = "Sample_Type")
     )
    )
  )))
  ),
  wellPanel(
   plotOutput("DNORM_plot", height = "600px", width = "920px")
  )
) # End of Data Normalization
) # End of QC tab box
), #End of "QC and Data Normalization"

# Data Exploratory panel
# Exploratory data analysis Panel    
tabPanel(title = "Exploratory data analysis",
 tabBox(title = "Exploratory data analysis", width = 12,
 tabPanel(title = "PCA Plot",
  fluidRow(
  column(
  width = 6,
  wellPanel(
   # Download Options
   tags$b("Note: The following options apply only to the downloaded plot."),
   tags$hr(),
   fluidRow(
   column(6, numericInput("PCAPh_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
   column(6, numericInput("PCAPw_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
   column(6, numericInput("PCAP_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
   column(6, radioButtons("PCAP_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
   ),
   downloadButton("pca_plottd", "Download Plot")
  )
  ),
  column(
  width = 6,
  wellPanel(
   radioButtons("PCAP_plt", "Select Color Palette for plots:",
    choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
    inline = TRUE),
   fluidRow(
    column(6,
    checkboxInput("PCAP_show_x_labels", "Show X Axis Label", value = TRUE),
    checkboxInput("PCAP_show_y_labels", "Show Y Axis Label", value = TRUE)
    ),
    column(6,
    checkboxInput("PCAP_show_plot_title", "Show Plot Title", value = TRUE),
    checkboxInput("PCAP_show_legend", "Show Legend", value = TRUE)
    )
    ),
   tags$br(),
   fluidRow(
   column(6,
    conditionalPanel(
     condition = "input.PCAP_show_x_labels",
     textInput("PCAP_x_axis_label", "X Axis Label", value = "PC1")
     ),
   conditionalPanel(
    condition = "input.PCAP_show_y_labels",
    textInput("PCAP_y_axis_label", "Y Axis Label", value = "PC2")
   )),
   column(6,
    conditionalPanel(
     condition = "input.PCAP_show_plot_title",
     textInput("PCAP_plot_title", "Plot Title", value = "PCA plot")
     ),
   conditionalPanel(
    condition = "input.PCAP_show_legend",
    textInput("PCAP_legend_title", "Legend Title", value = "Sample_Type")
   )
   )
  )))
 ),
 wellPanel(
  plotOutput("pca_plott", height = "600px", width = "920px")
  )
), # END of PCA plot panel
# MDS plot
tabPanel(title = "MDS Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("MDSP2h_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("MDSP2w_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("MDSP2_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("MDSP2_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("MDS_plottd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 radioButtons("MDSP2_plt", "Select Color Palette for plots:",
  choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
  inline = TRUE),
 fluidRow(
 column(6,
 checkboxInput("MDSP2_show_x_labels", "Show X Axis Label", value = TRUE),
 checkboxInput("MDSP2_show_y_labels", "Show Y Axis Label", value = TRUE)
 ),
 column(6,
 checkboxInput("MDSP2_show_plot_title", "Show Plot Title", value = TRUE),
 checkboxInput("MDSP2_show_legend", "Show Legend", value = TRUE)
 )
 ),
 tags$br(),
 fluidRow(
 column(6,
 conditionalPanel(
  condition = "input.MDSP2_show_x_labels",
  textInput("MDSP2_x_axis_label", "X Axis Label", value = "Leading logFC dim1")
  ),
 conditionalPanel(
  condition = "input.MDSP2_show_y_labels",
  textInput("MDSP2_y_axis_label", "Y Axis Label", value = "Leading logFC dim2")
 )),
 column(6,
 conditionalPanel(
  condition = "input.MDSP2_show_plot_title",
  textInput("MDSP2_plot_title", "Plot Title", value = "MDS plot")
 ),
 conditionalPanel(
  condition = "input.MDSP2_show_legend",
  textInput("MDSP2_legend_title", "Legend Title", value = "Sample_Type")
 ))
 )
 )
 )
 ),
 wellPanel(
  plotOutput("MDS_plott", height = "600px", width = "920px")
 )
), #END of MDS plot
# Voom plot
tabPanel(title = "Voom Plot",
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(3, numericInput("VMh_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(3, numericInput("VMw_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(3, numericInput("VM_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(3, radioButtons("VM_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("voom_plottd", "Download Plot")
 ),
 wellPanel(
   plotOutput("voom_plott", height = "600px", width = "920px")
 )
), # End of voom plot
# Mean-variance trend plot
tabPanel(title = "Mean-variance trend Plot",
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(3, numericInput("meanh_1", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(3, numericInput("meanw_1", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(3, numericInput("meanv_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(3, radioButtons("meanv_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("meanv_plottd", "Download Plot")
 ),
 wellPanel(
 plotOutput("meanv_plott", height = "600px", width = "920px")
 )
) # END of Mean-variance trend plot
) # END of Exploratory data analysis tab box
), # END of Exploratory data analysis Panel         
# DE genes panel
tabPanel(title = "DE results",
 tabBox(title = "DE results", width = 12,
 tabPanel(title = "DE table",
  infoBoxOutput("DE_box", width = 12),
  DT::dataTableOutput("DE_genes"),
  tags$hr(),
  wellPanel(
    # Download Options
    tags$b("Note: The following options apply only to the downloaded files."),
    tags$hr(),
    column(3, numericInput("DE_top", "Top hits:", min = 1, max = 999999, value = 10)),
    tags$br(),
    tags$hr(),
    tags$br(),
    tags$b("Note: By default, it downloads the top 10 hits. You can change '10' to any number to customize."),
    tags$br(),
    tags$br(),
    downloadButton("DE_resultsd", "Download all results"),
    downloadButton("Top_DE", "Download Top Hits results"),
    tags$br(),
    tags$br(),
    downloadButton("SIG_UPDWd", "Download Up and Down regulated genes results"),
    downloadButton("allgenes_UPd", "Download Up-regulated genes only"),
    downloadButton("allgenes_DWd", "Download Down-regulated genes only")
  )
), # END of DE table tabpanel
# Volcano Plot
  tabPanel(title = "Volcano Plot",
   fluidRow(
   column(
   width = 6,
   wellPanel(
   # Download Options
    tags$b("Note: The following options apply only to the downloaded plot."),
    tags$hr(),
    fluidRow(
    column(6, numericInput("vl_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
    column(6, numericInput("vl_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
    column(6, numericInput("vl_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
    column(6, radioButtons("vl_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
    ),
    downloadButton("volc_plotd", "Download Plot")
   )
   ),
   column(
   width = 6,
   wellPanel(
    # Plot Customization Inputs
    fluidRow(
    column(6, radioButtons("vl_1", "Criteria to denote significance:", choices = list("adjP", "P"), inline = TRUE, selected = "adjP"),
              helpText("*adjP for adjusted p-value and *P for p-value.")),
    column(6, checkboxInput("vl_2", "See a dotted line to indicate the lfc threshold in the plot.", TRUE)),
    column(6, checkboxInput("vl_3", "See a dotted line to indicate the p-value cutoff.", TRUE)),
    column(6, checkboxInput("vl_4", "Label Genes:", TRUE)),
    conditionalPanel(
     condition = "input.vl_4 == true",
     column(6, numericInput("vl_5", "Text size in the plot:", min = 1, max = 100000, value = 15)),
     column(6, numericInput("vl_6", "Top Hits to show in the plot:", min = 1, max = 100000, value = 10))
    )
    ),
   radioButtons("vl_7", "Select Color Palette for plots:",
                 choices = list("green", "red", "gray", "teal", "purple", "orange", "blue", "magenta"),
                 inline = TRUE, selected = "green"),
   textInput("vl_8", "Plot Title", value = "Volcano Plot")
   )
  )
  ),
  wellPanel(
  plotOutput("volc_plot", height = "600px", width = "920px")
 )
), # End of volcano plot
# Heat Map
tabPanel(title = "DE Heat Map Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 # Sample Selection Inputs
 fluidRow(
 column(6, numericInput("ht_1x", "First DE gene to plot:", min = 1, max = 999999, value = 1)),
 column(6, numericInput("ht_2x", "Last DE gene to plot:", min = 1, max = 999999, value = 10))
 ),
 tags$b("Note: If you want to view all genes in the plot, replace the value 10 with the total number of genes shown below."),
 tags$b(textOutput("heat_gtx"))
 ),
 wellPanel(
 fluidRow(
 column(6, numericInput("ht_3x", "First Sample to plot:", min = 1, max = 999999, value = 1)),
 column(6, numericInput("ht_4x", "Last Sample to plot:", min = 1, max = 999999, value = 10))
 ),
 tags$b("Note: If you want to view all samples in the plot, replace the value 10 with the total number of samples shown below."),
 tags$b(textOutput("heat_clx"))
 ),
 wellPanel(
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$br(),
 tags$br(),
 fluidRow(
 column(6, numericInput("ht_hx", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("ht_wx", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("ht_dpix", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("ht_typex", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("heat_plotdx", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(column(6,
 textInput("ht_5x", "Plot Title", value = "Heat map Plot"))),
 fluidRow(column(6,
 numericInput("ht_14x", "Font Size:", min = 1, max = 5000, value = 10))),
 radioButtons("ht_6x", "Select Color Palette:", 
 choices = c("RdYlBu","BluYl", "Greens", "Oranges", "Reds", "YlOrBr", "YlOrRd", "BuPu", "GnBu", "OrRd"),
 inline = TRUE, selected = "RdYlBu"),
 radioButtons("ht_7x", "Scale data:", choices = list("row", "column","none"), inline = TRUE, selected = "row"),
 radioButtons("ht_8x", "Border Color:", choices = list("Black", "Grey","Brown", "Purple", "Blue"), inline = TRUE, selected = "Black"),
 checkboxInput("ht_9x", "Cluster Rows", TRUE),
 checkboxInput("ht_10x", "Cluster Columns", TRUE),
 checkboxInput("ht_11x", "Show Legend", TRUE),
 checkboxInput("ht_12x", "Show Column names", TRUE),
 checkboxInput("ht_13x", "Show Row names", TRUE),
 checkboxInput("ht_20x", "Show Clinical/Annotation data", TRUE)
 )
 )
 ),
 wellPanel(
 plotOutput("heat_plotx", height = "600px", width = "920px")
 )
)# End of heat map with cli
) # END od DE tabbox
),# END of DE tabpanel

# ML results
tabPanel(title = "Gene modeling",
 tabBox(title = "Gene modeling Results", width = 12,
 tabPanel(title = "Performance Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("MLPR_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("MLPR_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("MLPR_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("MLPR_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("MLPR_perf_plotd", "Download Plot"),
 downloadButton("model_results_matd", "Download Performance matrix")
 )
 ),
 column(
 width = 6,
 wellPanel(
 textInput("MLPR_ttl", "Plot Title", value = "Performance Plot"),
 radioButtons("MLPR_plt", "Select Color Palette for plots:",
  choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
  inline = TRUE),
 fluidRow(column(6,numericInput("MLPR_txt", "Text Size", min = 1, max = 5000, value = 10)))
 )
 )
 ),
 wellPanel(
  plotOutput("MLPR_perf_plot", height = "600px", width = "920px")
  )
),

tabPanel(title = "Variable Importance Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(5, numericInput("VIMP_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(5, numericInput("VIMP_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(5, numericInput("VIMP_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(5, radioButtons("VIMP_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("VIMP_plottd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(
   column(6,
          textInput("VIMP_ttl", "Plot Title", value = "Variable Importance Plot"),
          radioButtons("VIMP_ptype", "Choose Plot Type:",
                       choices = c("Bar" = "bar", "Lollipop" = "lollipop"), inline = TRUE, selected = "lollipop"),
          numericInput("VIMP_txt", "Text Size", min = 1, max = 5000, value = 10)),
   column(6, numericInput("VIMP_topimp","Top Importance:",min = 1, max = 100000, value = 10),
          numericInput("VIMP_rw","Number of rows for plots:",min = 1, max = 10000, value = 2),
          numericInput("VIMP_clm","Number of columns for plots:",min = 1, max = 10000, value = 2)
 )),
 radioButtons("VIMP_plt", "Select Color Palette for plots:",
              choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
              inline = TRUE),
  )
 )
 ),
 wellPanel(
  plotOutput("VIMP_plott", height = "600px", width = "920px")
 )
), # End of VarImp plot
# ROC plot
tabPanel(title = "ROC Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(5, numericInput("ROCML_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(5, numericInput("ROCML_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(5, numericInput("ROCML_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(5, radioButtons("ROCML_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("ROCML_plottd", "Download Plot"),
 downloadButton("ROCML_cnf_table", "Download Confusion matrix")
 )
 ),
 column(
  width = 6,
  wellPanel(
  fluidRow(
  column(6,
  textInput("ROCML_ttl", "Plot Title", value = "ROC Plot"),
  numericInput("ROCML_txt", "Text Size", min = 1, max = 5000, value = 10),
  checkboxInput("ROCML_mlt", "Multiple Plots", TRUE)
  ),
  ),
  radioButtons("ROCML_plt", "Select Color Palette for plots:",
                choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
                inline = TRUE)
  )
  )
  ),
  wellPanel(
  plotOutput("ROCML_plott", height = "600px", width = "920px")
  )
) # End of ROC plot
) # END of ML tabbox
), # END of ML tabpanel

# ORA 
tabPanel(title = "ORA results",
 tabBox(title = "ORA results", width = 12,
 tabPanel(title = "GO over-representation",
fluidRow(
  wellPanel(
  # DE gene Selection Inputs
  fluidRow(column(width = 4, numericInput("GOTOPH","Number of Genes to Select for GO analysis:", value = 10, min = 1))),
  tags$b("Note: By default, the top ten differentially expressed (DE) genes are used for Gene Ontology (GO) over-representation analysis."),
  tags$br(),
  tags$b("The following are the total DE genes identified, and you can replace the number 10 with any other number that suits your data analysis."),
  tags$b(textOutput("GOTOPH_DE")),
  tags$br(),
  tags$b("Please click the *Update Results* button to get the latest GO results."),
  tags$br(),
  actionButton("update_ORA1", "Update Results")
  ),
  DT::dataTableOutput("ORA_GO_table"),
  wellPanel(
  tags$b("Click the button below to download the GO over-representation results."),
  tags$br(),
  tags$br(),
  downloadButton("ORA_GO_tabled", "Download GO over-representation results")
  )
  )
), # END of GO results
                
tabPanel(title = "Bar Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 # Download Options
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("GOB_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("GOB_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("GOB_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("GOB_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("ORA_barplotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(
 column(6,
 textInput("GOB_pt", "Plot Title", value = "Bar Plot"),
 numericInput("GOB_ct", "Show GO category:", min = 1, max = 99999999, value = 10),
 numericInput("GOB_fz", "Font Size:", min = 1, max = 5000, value = 8)
 )
 )
 )
 )
 ),
 wellPanel(
 plotOutput("ORA_barplot", height = "600px", width = "920px")
 )
), # END of GO bar plot
# Dot plot
tabPanel(title = "Dot Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("DOT_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("DOT_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("DOT_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("DOT_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("ORA_dotplotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(
 column(6,
 textInput("DOT_pt", "Plot Title", value = "Dot Plot"),
 numericInput("DOT_ct", "Show GO category:", min = 1, max = 99999999, value = 10),
 numericInput("DOT_fz", "Font Size:", min = 1, max = 5000, value = 8)
    )
   )
  )
 )
 ),
 wellPanel(
  plotOutput("ORA_dotplot", height = "600px", width = "920px")
  )
), # END of GO Dot plot
#CNET Plot
tabPanel(title = "Gene-Concept Network Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("CNT_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("CNT_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("CNT_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("CNT_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("ORA_CNTplotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(
 column(7,
 textInput("CNT_ttl", "Plot Title", value = "Gene-Concept Network Plot"),
 numericInput("CNT_ct", "Show GO category:", min = 1, max = 99999, value = 5),
 checkboxInput("CNT_fz", "Coloring edge by enriched terms", FALSE),
 checkboxInput("CNT_pt", " Circular layout", TRUE)
 )
 )
 )
 )
 ),
 wellPanel(
 plotOutput("ORA_CNTplot", height = "800px", width = "920px")
 )
), # END of GO CNT plot
# Heat plot
tabPanel(title = "Heat Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("HGT_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("HGT_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("HGT_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("HGT_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("ORA_HGTplotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(
 column(7,
 textInput("HGT_ttl", "Plot Title", value = "Heat Plot"),
 numericInput("HGT_ct", "Show GO category:", min = 1, max = 99999999, value = 5)
 )
 )
 )
 )
 ),
 wellPanel(
 plotOutput("ORA_HGTplot", height = "600px", width = "920px")
 )
), # END of GO Heat plot
# Enrichment Map plot
tabPanel(title = "Enrichment Map Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
  column(6, numericInput("ENRT_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
  column(6, numericInput("ENRT_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
  column(6, numericInput("ENRT_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
  column(6, radioButtons("ENRT_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
  ),
 downloadButton("ORA_ENRTplotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
  fluidRow(
  column(7,
  textInput("ENRT_ttl", "Plot Title", value = "Enrichment Map plot"),
  numericInput("ENRT_ct", "Show GO category:", min = 1, max = 99999999, value = 30)
  )
 )
 )
 )
 ),
 wellPanel(
 plotOutput("ORA_ENRTplot", height = "600px", width = "920px")
 )
), # END of GO Enrichment Map plot
# UpSet Plot
tabPanel(title = "UpSet Plot",
 fluidRow(
 column(
 width = 6,
 wellPanel(
 tags$b("Note: The following options apply only to the downloaded plot."),
 tags$hr(),
 fluidRow(
 column(6, numericInput("UPS_h", "Height of the plot (inch):", min = 1, max = 5000, value = 8)),
 column(6, numericInput("UPS_w", "Width of the plot (inch):", min = 1, max = 5000, value = 12)),
 column(6, numericInput("UPS_dpi", "Plot resolution (dpi):", min = 1, max = 5000, value = 300)),
 column(6, radioButtons("UPS_type", "Select file type:", choices = list("png", "pdf"), inline = TRUE))
 ),
 downloadButton("ORA_UPSplotd", "Download Plot")
 )
 ),
 column(
 width = 6,
 wellPanel(
 fluidRow(
 column(9,
 textInput("UPS_ttl", "Plot Title", value = "UpSet Plot"),
 numericInput("UPS_ct", "Show GO category:", min = 1, max = 99999999, value = 10)
 )
 )
 )
 )
 ),
 wellPanel(
 plotOutput("ORA_UPSplot", height = "600px", width = "920px")
 )
), # END of GO UpSet Plot
# KEGG Results
tabPanel(title = "KEGG Results",
 DT::dataTableOutput("ORA_KEGG_table"),
 wellPanel(
 tags$b("Click the button below to download the KEGG results."),
 tags$br(),
 tags$br(),
 downloadButton("ORA_KEGG_tabled", "Download KEGG results")
 )
),
# KEGG plot
tabPanel(title = "KEGG Plot",
 wellPanel(
 fluidRow(
 column(6,
 selectizeInput("pathwayIds",
 label = "Select KEGG Pathway ID",
 choices = NULL,  # Choices will be populated dynamically
 multiple = FALSE,
 options = list(
 placeholder = "Start typing pathway ID"
 )
 )
 )),
 fluidRow(
 column(6,
 actionButton("generatePathview", "Generate Pathview", class = "btn btn-info")
 )
 ),
 fluidRow(
 # Conditional panel to display download buttons when plots are available
 conditionalPanel(
  condition = "output.pathviewPlotsAvailable",
  column(3,
  downloadButton("downloadPathviewPng", "Download .png", class = "btn btn-warning", style = "margin: 7px;")
  ),
  column(3,
  downloadButton("downloadPathviewPdf", "Download .pdf", class = "btn btn-warning", style = "margin: 7px;")
  )
  )
 ),
 tags$div(class = "clearBoth")
 ),
 wellPanel(
 withSpinner(
  imageOutput(outputId = "pathview_plot", inline = T, height = "600px", width = "920px"),
  type = 8,
  color = "#e87b07"  # Customize spinner color to match the theme
  ),
 style = "text-align: center; padding: 5px; overflow: auto; max-height: 650px; max-width: 940px;"
 )  
 ) # END of pathview
 ) # END of ORA tabbox
), 
# END of ORA tabpanel

# Report
tabPanel(title = "Report",
wellPanel(
  h3("1. The report includes methods, results, and data analysis parameters."),
  tags$br(),
  h3("2. Choose your preferred file format and click the download button to get your report."),
  tags$br(),
  radioButtons("RMD_format", "Select Report Format:",
            choices = c("PDF" = "pdf_document",
                        "Word" = "word_document",
                        "HTML" = "html_document"), inline = TRUE),
tags$br(),
downloadButton("RMD_Report", "Download Data analysis Report")
)
)

) # End of tabsetpanel  
) # End of "results" tabitem
# END of main code
 )
 ) 
 )
)

# Define server logic
server <- function(input, output, session) {

# Set global options for maximum request size (999999 MB)
options(shiny.maxRequestSize = 999999 * 1024^2)

# GO options
observeEvent(input$ABTR, {toggle("ABTRR")})
observeEvent(input$GOP1, {toggle("GOP11")})
observeEvent(input$GOP2, {toggle("GOP22")})
observeEvent(input$GOP3, {toggle("GOP33")})
observeEvent(input$GOP4, {toggle("GOP44")})
observeEvent(input$GOP5, {toggle("GOP55")})
observeEvent(input$GOP6, {toggle("GOP66")})
observeEvent(input$GOP7, {toggle("GOP77")})
observeEvent(input$GOP8, {toggle("GOP88")})
observeEvent(input$GOP9, {toggle("GOP99")})
observeEvent(input$GOP10, {toggle("GOP1010")})
observeEvent(input$GOP11, {toggle("GOP1111")})
observeEvent(input$GOP12, {toggle("GOP1212")})
observeEvent(input$GOP13, {toggle("GOP1313")})
observeEvent(input$GOP14, {toggle("GOP1414")})
  
# ML options  
observe({
 selectedOptions <- input$gmgmt5
 selectedCount <- length(selectedOptions)
    
    # Manage modals based on selection count
    if (selectedCount < 2) {
      session$sendInputMessage("showModal", "min")
    } else if (selectedCount > 4) {
      session$sendInputMessage("showModal", "max")
    } else {
      session$sendInputMessage("showModal", "none")
    }
})
  
  observe({
    selectedOptions <- input$gmgmt5
    selectedCount <- length(selectedOptions)
    
    # Manage modals based on selection count
    if (selectedCount < 2) {
      session$sendInputMessage("showModal", "min")
    } else if (selectedCount > 4) {
      session$sendInputMessage("showModal", "max")
    } else {
      session$sendInputMessage("showModal", "none")
    }
  })
  
  observe({
    if (is.null(input$showModal)) return()
    
    if (input$showModal == 'min') {
      showModal(modalDialog(
        title = "ML model options!",
        "Please select a minimum of two ML models.",
        easyClose = FALSE,
        footer = modalButton("OK")
      ))
    } else if (input$showModal == 'max') {
      showModal(modalDialog(
        title = "ML model options!",
        "Please select a maximum of four ML models.",
        easyClose = FALSE,
        footer = modalButton("OK")
      ))
    } else if (input$showModal == 'none') {
      removeModal()
    }
  })
  
  # Observe submit button
  observeEvent(input$subm, {
    selectedOptions <- input$gmgmt5
    selectedCount <- length(selectedOptions)
    
    if (selectedCount < 2 || selectedCount > 4) {
      # Trigger modal and prevent submission
      if (selectedCount < 2) {
        showModal(modalDialog(
          title = "ML model options!",
          "Please select a minimum of two ML models.",
          easyClose = FALSE,
          footer = modalButton("OK")
        ))
      } else if (selectedCount > 4) {
        showModal(modalDialog(
          title = "ML model options!",
          "Please select a maximum of four ML models.",
          easyClose = FALSE,
          footer = modalButton("OK")
        ))
      }
    } 
  })
  
# To select group for DE
observe({
# Based on the selected dataset, generate the appropriate file path
dataset_prefix <- input$dataset_choice  # Get the selected dataset prefix
file2_path <- paste0("www/", dataset_prefix, "-metadata.txt")  # Construct the file path
  
# Read the metadata file from the 'www' folder
refgg <- read.table(file2_path, header = TRUE, sep = "\t", check.names = FALSE)
  
# Update selectInput with the second column of the metadata as choices
updateSelectInput(session, "refg",
                  choices = refgg[, 2],  # The second column as choices
                  selected = refgg[, 2][1])  # Set default selection to the first element
})

# Main observe
observeEvent(input$subm, {
 shinyalert(
 title = '<h2 style="text-align: center; color: #000000;">Please wait patiently :)</h2>',
 text = '<h3 style="text-align: center; color: #000000;">Your Data analysis has started...</h3><h3 style="text-align: center; font-size: 18px; color: #000000;">Wait until table and plots appear on the screen.</h3>',
 type = "info",
 html = TRUE,
 closeOnClickOutside = TRUE,
 closeOnEsc = TRUE
)

# Register the parallel backend
numCores11 <- parallelly::availableCores()
cl22 <- makePSOCKcluster(numCores11)
registerDoParallel(cl22)

# Update to the "results" tab
updateTabsetPanel(session, "tabs", selected = "results")

# Load the data
# Determine dataset prefix based on user selection
dataset_prefix <- input$dataset_choice  # No need for a reactive here
  
# Load the raw_gene_counts directly from the 'www/' folder
raw_gene_counts <- read.table(paste0("www/", dataset_prefix, "-counts.txt"), header = TRUE, sep = "\t", check.names = FALSE)
raw_gene_counts <- raw_gene_counts[1:input$sgn,]

# Load the sample metadata directly from the 'www/' folder
exp_design <- read.table(paste0("www/", dataset_prefix, "-metadata.txt"), header = TRUE, sep = "\t", check.names = FALSE)
  
# Load the clinical/annotation data directly from the 'www/' folder
clin_data <- read.table(paste0("www/", dataset_prefix, "-clinical_or_annotation.txt"), header = TRUE, sep = "\t", check.names = FALSE)
  
# Create a data frame to count rows and columns of uploaded data
# Gene counts file
r1 <- as.data.frame(raw_gene_counts)
output$rc_genes1 <- renderText({ nrow(r1) })

output$rc_samples1 <- renderText({ ncol(r1) })
output$srange1 <- renderText({ ncol(r1) })
output$cmp_smp <- renderText({ ncol(r1) })
    
output$rc_genes <- renderText({ nrow(r1) })
output$rc_samples <- renderText({ ncol(r1) })
output$srange <- renderText({ ncol(r1) })

output$rc_genes2 <- renderText({ nrow(r1) })

# Sample and Group file
samp1 <- as.data.frame(exp_design)
output$samp_rows1 <- renderText({ nrow(samp1) })
output$samp_columns1 <- renderText({ ncol(samp1) })
    
output$samp_rows <- renderText({ nrow(samp1) })
output$samp_columns <- renderText({ ncol(samp1) })
    
# View uploaded Data
output$raw_counts1 <- DT::renderDataTable({ raw_gene_counts }, options = list(scrollX = TRUE))
output$raw_counts11 <- DT::renderDataTable({ raw_gene_counts }, options = list(scrollX = TRUE))
output$sample_info1 <- DT::renderDataTable({ exp_design }, options = list(scrollX = TRUE))
output$raw_counts <- DT::renderDataTable({ raw_gene_counts }, options = list(scrollX = TRUE))
output$sample_info <- DT::renderDataTable({ exp_design }, options = list(scrollX = TRUE))
    
# Prepare data for Raw count data plot
raw_cnt_plot <- reactive({
  # Get the range of rows
  rcs_range <- input$rcs11:input$rcs21
  # Get user inputs
  count_type <- input$rc_dtr
  plot_type <- input$rc_plotpe
  selected_palette <- input$rc_plt
      
  # Define a function to get the color palette
  color_palette_func <- function(n) {
  set_col(palette = selected_palette, n = n)
  }
      
  # Call the plot function with inputs from the UI
  plot_raw_gene_counts(
   raw_gene_counts = raw_gene_counts, 
   exp_design = exp_design, 
   rcs_range = rcs_range,
   count_type = input$rc_dtr, 
   plot_type = input$rc_plotpe, 
   color_palette_func = color_palette_func,
   show_x_labels = input$show_x_labels, 
   show_y_labels = input$show_y_labels,
   x_axis_label = input$x_axis_label, 
   y_axis_label = input$y_axis_label, 
   plot_title = input$plot_title,
   show_plot_title = input$show_plot_title,
   legend_title = input$legend_title, 
   show_legend = input$show_legend
  )
})
    
# Plot
output$rc_plot1 <- renderPlot({ raw_cnt_plot() })
    
# Download plot
output$rc_plotd1 <- downloadHandler(
 filename = function() {
 paste("Raw_read_counts_plot", input$rc_type1, sep = ".")
 },
 content = function(file) {
    if (input$rc_type1 == "png") {
    ggsave(file, plot = raw_cnt_plot(), device = "png", width = input$rc_pw1, height = input$rc_ph1, units = "in", dpi = input$rc_dpi1, bg = "white")
    } else if (input$rc_type1 == "pdf") {
    ggsave(file, plot = raw_cnt_plot(), device = "pdf", width = input$rc_pw1, height = input$rc_ph1)
    } else {
    stop("Unsupported file type selected.")
    }
    }
 )

# Clinical/Annotation data
# Sample and Group file
clin_data_info <- as.data.frame(clin_data)
output$cli_rows <- renderText({ nrow(clin_data_info) })
output$cli_cols <- renderText({ ncol(clin_data_info) })

# View uploaded Data
output$cl_data_table <- DT::renderDataTable({clin_data},
                        options = list(scrollX = TRUE))

# Plot
# Prepare data for Raw count data plot
raw_cnt_plot1 <- reactive({
  # Get the range of rows
  rcs_range1 <- input$rcs1:input$rcs2
  
  # Get user inputs
  count_type <- input$rc_dtr1
  plot_type <- input$rc_plotpe1
  selected_palette <- input$rc_plt1
  
  # Define a function to get the color palette
  color_palette_func1 <- function(n) {
    set_col(palette = selected_palette, n = n)
  }
  
  # Call the plot function with inputs from the UI
  plot_raw_gene_counts(
    raw_gene_counts = raw_gene_counts, 
    exp_design = exp_design, 
    rcs_range = rcs_range1,
    count_type = input$rc_dtr1, 
    plot_type = input$rc_plotpe1, 
    color_palette_func = color_palette_func1,
    show_x_labels = input$show_x_labels1, 
    show_y_labels = input$show_y_labels1,
    x_axis_label = input$x_axis_label1, 
    y_axis_label = input$y_axis_label1, 
    plot_title = input$plot_title1,
    show_plot_title = input$show_plot_title1,
    legend_title = input$legend_title1, 
    show_legend = input$show_legend1
  )
})

# Plot
output$rc_plot11 <- renderPlot({ raw_cnt_plot1() })

# Download plot
output$rc_plotd11 <- downloadHandler(
  filename = function() {
    paste("Raw_read_counts_plot", input$rc_type, sep = ".")
  },
  content = function(file) {
    if (input$rc_type == "png") {
      ggsave(file, plot = raw_cnt_plot1(), device = "png", width = input$rc_pw, height = input$rc_ph, units = "in", dpi = input$rc_dpi, bg = "white")
    } else if (input$rc_type == "pdf") {
      ggsave(file, plot = raw_cnt_plot1(), device = "pdf", width = input$rc_pw, height = input$rc_ph)
    } else {
      stop("Unsupported file type selected.")
    }
  }
)

# Data Filtering/Pre-processing
# Filter low count genes by total row counts
Filt = rowSums(raw_gene_counts) >= input$dp1
Filt_counts_1 = raw_gene_counts[Filt,]
rm (Filt)

# QC Plot 1 - After total row count filtering plot
output$filt_1_genes <- renderText({ nrow(Filt_counts_1) })
output$filt_1_samples <- renderText({ ncol(Filt_counts_1) })

# Prepare data for Raw count data plot
filt_1_plot <- reactive({
  # Get the range of rows
  filt_1_rcs_range <- input$filt1_1:input$filt1_2
  
  # Get user inputs
  count_type <- input$filt1_dtr
  plot_type <- input$filt1_plotpe
  selected_palette <- input$filt1_plt
  
  # Define a function to get the color palette
  filt1_color_palette_func <- function(n) {
    set_col(palette = selected_palette, n = n)
  }
  
  # Call the plot function with inputs from the UI
  plot_raw_gene_counts(
    raw_gene_counts = Filt_counts_1, 
    exp_design = exp_design, 
    rcs_range = filt_1_rcs_range,
    count_type = input$filt1_dtr, 
    plot_type = input$filt1_plotpe, 
    color_palette_func = filt1_color_palette_func,
    show_x_labels = input$filt1_show_x_labels, 
    show_y_labels = input$filt1_show_y_labels,
    x_axis_label = input$filt1_x_axis_label, 
    y_axis_label = input$filt1_y_axis_label, 
    plot_title = input$filt1_plot_title,
    show_plot_title = input$filt1_show_plot_title,
    legend_title = input$filt1_legend_title, 
    show_legend = input$filt1_show_legend
  )
})

# Plot
output$filt1_plot <- renderPlot({ filt_1_plot() })

# Download plot
output$filt1_plotd <- downloadHandler(
  filename = function() {
    paste("Filtered_counts_plot", input$filt1_type, sep = ".")
  },
  content = function(file) {
    if (input$filt1_type == "png") {
      ggsave(file, plot = filt_1_plot(), device = "png", width = input$filt1w_1, height = input$filt1h_1, units = "in", dpi = input$filt1_dpi, bg = "white")
    } else if (input$filt1_type == "pdf") {
      ggsave(file, plot = filt_1_plot(), device = "pdf", width = input$filt1w_1, height = input$filt1h_1)
    } else {
      stop("Unsupported file type selected.")
    }
  }
)
# End of Filtered_read_counts_plot

# Prepare group
group = factor(exp_design[,2])
group = relevel(group, ref= input$refg)

# Design
design <- model.matrix(~ group)

# Creating a DGEList object
dgeF = DGEList(Filt_counts_1)

# FILTER COUNTS
# Data Filtering by library size
keep = filterByExpr(dgeF,design, min.count= input$dp2) # defining which genes to keep
dgeF = dgeF[keep,]
rm(keep)

# FilterByExpr plot
# QC Plot2 - FilterByExpr filtering plot
filt2c <- as.data.frame(dgeF$counts)

output$filt2_genes <- renderText({ nrow(filt2c) })
output$filt2_samples <- renderText({ ncol(filt2c) })
output$filt22_genes <- renderText({ nrow(filt2c) })

# Plot
# Prepare data for FilterByExpr plot
filt_2_plot <- reactive({
  # Get the range of rows
  filt_2_rcs_range <- input$filt2_1:input$filt2_2
  
  # Get user inputs
  count_type <- input$filt2_dtr
  plot_type <- input$filt2_plotpe
  selected_palette <- input$filt2_plt
  
  # Define a function to get the color palette
  filt2_color_palette_func <- function(n) {
    set_col(palette = selected_palette, n = n)
  }
  
  # Call the plot function with inputs from the UI
  plot_raw_gene_counts(
    raw_gene_counts = filt2c, 
    exp_design = exp_design, 
    rcs_range = filt_2_rcs_range,
    count_type = input$filt2_dtr, 
    plot_type = input$filt2_plotpe, 
    color_palette_func = filt2_color_palette_func,
    show_x_labels = input$filt2_show_x_labels, 
    show_y_labels = input$filt2_show_y_labels,
    x_axis_label = input$filt2_x_axis_label, 
    y_axis_label = input$filt2_y_axis_label, 
    plot_title = input$filt2_plot_title,
    show_plot_title = input$filt2_show_plot_title,
    legend_title = input$filt2_legend_title, 
    show_legend = input$filt2_show_legend
  )
})

# Plot
output$filt2_plot <- renderPlot({ filt_2_plot() })

# Download plot
output$filt2_plotd <- downloadHandler(
  filename = function() {
    paste("FilterByExpr_plot", input$filt2_type, sep = ".")
  },
  content = function(file) {
    if (input$filt2_type == "png") {
      ggsave(file, plot = filt_2_plot(), device = "png", width = input$filt2w_1, height = input$filt2h_1, units = "in", dpi = input$filt2_dpi, bg = "white")
    } else if (input$filt2_type == "pdf") {
      ggsave(file, plot = filt_2_plot(), device = "pdf", width = input$filt2w_1, height = input$filt2h_1)
    } else {
      stop("Unsupported file type selected.")
    }
  }
 )
# End of FilterByExpr plot

# Before and after filtering
cmp_filt <- reactive({
  # Calculate log-cpm values
  lcpm <- cpm(raw_gene_counts, log=TRUE)
  lcpm2 <- cpm(filt2c, log=TRUE)
  
  nsamples <- ncol(raw_gene_counts[input$cmp_1:input$cmp_2])
  
  if (nsamples <= 0) {
    stop("Invalid number of samples. Please check the input.")
  }
  
  # Generate color palette with Viridis
  color_palette <- viridis(n = nsamples)
  sample_names <- colnames(raw_gene_counts)[input$cmp_1:input$cmp_2]
  
  # Prepare data for plotting
  plot_data_raw <- lcpm %>% 
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "Sample", values_to = "Log_cpm") %>%
    mutate(Sample = factor(Sample, levels = sample_names))
  
  plot_data_filtered <- lcpm2 %>% 
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "Sample", values_to = "Log_cpm") %>%
    mutate(Sample = factor(Sample, levels = sample_names))
  
  # Create ggplot objects for raw and filtered data
  p1 <- ggplot(plot_data_raw, aes(x=Log_cpm, color=Sample)) +
    geom_density(size=1) +
    scale_color_viridis(discrete = TRUE) +
    labs(
      title = if (!is.null(input$cmp_show_plot_title) && input$cmp_show_plot_title) input$cmp_plot_title else NULL,
      x = if (!is.null(input$cmp_show_x_labels) && input$cmp_show_x_labels) input$cmp_x_axis_label else NULL,
      y = if (!is.null(input$cmp_show_y_labels) && input$cmp_show_y_labels) input$cmp_y_axis_label else NULL
    ) +
    extra_theme() +
    theme(
      plot.title = if (!is.null(input$cmp_show_plot_title) && input$cmp_show_plot_title) element_text(color = "black", size = 12, face = "bold", hjust = 0.5) else element_blank(),
      axis.title.x = if (!is.null(input$cmp_show_x_labels) && input$cmp_show_x_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.title.y = if (!is.null(input$cmp_show_y_labels) && input$cmp_show_y_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.text.x = if (!is.null(input$cmp_show_x_labels) && input$cmp_show_x_labels) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.text.y = if (!is.null(input$cmp_show_y_labels) && input$cmp_show_y_labels) element_text(size = 8) else element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  p2 <- ggplot(plot_data_filtered, aes(x=Log_cpm, color=Sample)) +
    geom_density(size=1) +
    scale_color_viridis(discrete = TRUE) +
    labs(
      title = if (!is.null(input$cmp1_show_plot_title) && input$cmp1_show_plot_title) input$cmp1_plot_title else NULL,
      x = if (!is.null(input$cmp1_show_x_labels) && input$cmp1_show_x_labels) input$cmp1_x_axis_label else NULL,
      y = if (!is.null(input$cmp1_show_y_labels) && input$cmp1_show_y_labels) input$cmp1_y_axis_label else NULL
    ) +
    extra_theme() +
    theme(
      plot.title = if (!is.null(input$cmp1_show_plot_title) && input$cmp1_show_plot_title) element_text(color = "black", size = 12, face = "bold", hjust = 0.5) else element_blank(),
      axis.title.x = if (!is.null(input$cmp1_show_x_labels) && input$cmp1_show_x_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.title.y = if (!is.null(input$cmp1_show_y_labels) && input$cmp1_show_y_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.text.x = if (!is.null(input$cmp1_show_x_labels) && input$cmp1_show_x_labels) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.text.y = if (!is.null(input$cmp1_show_y_labels) && input$cmp1_show_y_labels) element_text(size = 8) else element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  # Arrange plots side by side
  grid.arrange(p1, p2, ncol=2)
})

# Plot
output$cmp_plot <- renderPlot({ cmp_filt() })

# Download plot
output$cmp_plotd <- downloadHandler(
  filename = function() {
    paste("Before_and_after_data_filtering_plot", input$cmp_type, sep = ".")
  },
  content = function(file) {
    if (input$cmp_type == "png") {
      ggsave(file, plot = cmp_filt(), device = "png", width = input$cmp_pw, height = input$cmp_ph, units = "in", dpi = input$cmp_dpi, bg = "white")
    } else if (input$cmp_type == "pdf") {
      ggsave(file, plot = cmp_filt(), device = "pdf", width = input$cmp_pw, height = input$cmp_ph)
    } else {
      stop("Unsupported file type selected.")
    }
  }
) # End of before and after filtering 

# Data Normalization
dgeF = calcNormFactors(dgeF,method = input$dp3)

# Voom plot
vm = voom(dgeF,design,plot= FALSE)

# For plot only
normpl <- as.data.frame (vm$E)
output$DNORM_genes <- renderText({ nrow(normpl) })
output$DNORM_samples <- renderText({ ncol(normpl) })

# Plot
# Prepare data for Data_Normalization_plot
DNORM_plot <- reactive({
  # Get the range of rows
  DNORM_rcs_range <- input$DNORM_1:input$DNORM_2
  
  # Get user inputs
  count_type <- input$DNORM_dtr
  plot_type <- input$DNORM_plotpe
  selected_palette <- input$DNORM_plt
  
  # Define a function to get the color palette
  DNORM_color_palette_func <- function(n) {
    set_col(palette = selected_palette, n = n)
  }
  
  # Call the plot function with inputs from the UI
  plot_raw_gene_counts(
    raw_gene_counts = normpl, 
    exp_design = exp_design, 
    rcs_range = DNORM_rcs_range,
    count_type = "Count", 
    plot_type = input$DNORM_plotpe, 
    color_palette_func = DNORM_color_palette_func,
    show_x_labels = input$DNORM_show_x_labels, 
    show_y_labels = input$DNORM_show_y_labels,
    x_axis_label = input$DNORM_x_axis_label, 
    y_axis_label = input$DNORM_y_axis_label, 
    plot_title = input$DNORM_plot_title,
    show_plot_title = input$DNORM_show_plot_title,
    legend_title = input$DNORM_legend_title, 
    show_legend = input$DNORM_show_legend
  )
})

# Plot
output$DNORM_plot <- renderPlot({ DNORM_plot() })

# Download plot
output$DNORM_plotd <- downloadHandler(
  filename = function() {
    paste("Data_Normalization_plot", input$DNORM_type, sep = ".")
  },
  content = function(file) {
    if (input$DNORM_type == "png") {
      ggsave(file, plot = DNORM_plot(), device = "png", width = input$DNORMw_1, height = input$DNORMh_1, units = "in", dpi = input$DNORM_dpi, bg = "white")
    } else if (input$DNORM_type == "pdf") {
      ggsave(file, plot = DNORM_plot(), device = "pdf", width = input$DNORMw_1, height = input$DNORMh_1)
    } else {
      stop("Unsupported file type selected.")
    }
  }
) # End of Data_Normalization_plot
# PCA plot
#PCA plot data prep
vm$targets$grp <- samp1[,2]

# Define the PCA plot as a reactive expression
# Reactive expression for PCA plot
pca_plotr <- reactive({
  # Ensure vm is available in the global or server context
  group <- factor(vm$targets$grp)
  
  # Perform PCA
  pca <- prcomp(t(vm$E))
  
  # Create a data frame for ggplot
  pca_data <- data.frame(PC1 = pca$x[,1], 
                         PC2 = pca$x[,2], 
                         Group = group)
  
  # Number of levels in the group variable
  n_groups <- length(levels(group))
  
  # Select color palette based on user input
  palette_name <- input$PCAP_plt
  colors <- set_col(palette_name, n_groups) # Function that returns color values for palettes
  
  # Generate the PCA plot using ggplot2
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    extra_theme() +
    theme(
      plot.title = if (!is.null(input$PCAP_show_plot_title) && input$PCAP_show_plot_title) element_text(color = "black", size = 12, face = "bold", hjust = 0.5) else element_blank(),
      axis.title.x = if (!is.null(input$PCAP_show_x_labels) && input$PCAP_show_x_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.title.y = if (!is.null(input$PCAP_show_y_labels) && input$PCAP_show_y_labels) element_text(color = "black", size = 10, face = "bold") else element_blank(),
      axis.text.x = if (!is.null(input$PCAP_show_x_labels) && input$PCAP_show_x_labels) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.text.y = if (!is.null(input$PCAP_show_y_labels) && input$PCAP_show_y_labels) element_text(size = 8) else element_blank(),
      legend.title = if (!is.null(input$PCAP_show_legend) && input$PCAP_show_legend) element_text(size = 10, face = "bold") else element_blank(),
      legend.text = if (!is.null(input$PCAP_show_legend) && input$PCAP_show_legend) element_text(size = 8) else element_blank(),
      legend.position = if (!is.null(input$PCAP_show_legend) && input$PCAP_show_legend) "bottom" else "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
    
  # Conditional labels and titles
  if (input$PCAP_show_x_labels) {
    pca_plot <- pca_plot + labs(x = input$PCAP_x_axis_label)
  }
  if (input$PCAP_show_y_labels) {
    pca_plot <- pca_plot + labs(y = input$PCAP_y_axis_label)
  }
  if (input$PCAP_show_plot_title) {
    pca_plot <- pca_plot + labs(title = input$PCAP_plot_title)
  }
  if (input$PCAP_show_legend) {
    pca_plot <- pca_plot + labs(color = input$PCAP_legend_title)
  }
  
  return(pca_plot)
})

# PCA plot
output$pca_plott <- renderPlot({pca_plotr()})

# Download plot
output$pca_plottd <- downloadHandler(
  filename = function() {
    paste("PCA_plot", input$PCAP_type, sep = ".")
  },
  content = function(file) {
    if (input$PCAP_type == "png") {
      ggsave(file, plot = pca_plotr(), device = "png", width = input$PCAPw_1, height = input$PCAPh_1, units = "in", dpi = input$PCAP_dpi, bg = "white")
    } else if (input$PCAP_type == "pdf") {
      ggsave(file, plot = pca_plotr(), device = "pdf", width = input$PCAPw_1, height = input$PCAPh_1)
    } else {
      stop("Unsupported file type selected.")
    }
  }
) # End of PCA plot
# MDS plot
MDS_plot <- reactive({
  palette_choice <- input$MDSP2_plt
  x_axis_label <- if (input$MDSP2_show_x_labels) input$MDSP2_x_axis_label else NULL
  y_axis_label <- if (input$MDSP2_show_y_labels) input$MDSP2_y_axis_label else NULL
  plot_title <- if (input$MDSP2_show_plot_title) input$MDSP2_plot_title else NULL
  legend_title <- if (input$MDSP2_show_legend) input$MDSP2_legend_title else NULL
  
  num_groups <- length(levels(group))
  colors <- set_col(palette_choice, num_groups)
  
  # Create a temporary plot object to capture the plot
  plot_object <- function() {
    # Increase right margin to ensure the legend fits
    par(mar = c(6, 5, 4, 14) + 0.1)  # Right margin adjusted further
    
    plotMDS(dgeF, col = colors[as.numeric(group)], main = plot_title, pch = 19,
            xlab = x_axis_label, ylab = y_axis_label)
    
    if (!is.null(legend_title)) {
      # Move legend further to the right
      legend("topright", inset = c(-0.2, 0), legend = levels(group), pch = 19, 
             col = colors, title = legend_title, xpd = TRUE, cex = 0.8)  # Adjust cex if needed
    }
  }
  
  plot_object
})

# MDS plot
output$MDS_plott <- renderPlot({ MDS_plot()()})

# Output for downloading plot
output$MDS_plottd <- downloadHandler(
  filename = function() {
    paste("MDS_plot", input$MDSP2_type, sep = ".")
  },
  content = function(file) {
    if (input$MDSP2_type == "png") {
      png(file, width = input$MDSP2w_1, height = input$MDSP2h_1, units = "in", res = input$MDSP2_dpi)
    } else if (input$MDSP2_type == "pdf") {
      pdf(file, width = input$MDSP2w_1, height = input$MDSP2h_1)
    } else {
      stop("Unsupported file type selected.")
    }
    
    # Render the plot
    print(MDS_plot()())  # Call the reactive plot function
    
    # Close the device
    dev.off()
  }
)
# Voom Reactive
voom_p <- reactive({tx = voom(dgeF,design,plot=TRUE)})

# Voom plot
output$voom_plott <- renderPlot({voom_p()})

# Voom download
output$voom_plottd <- downloadHandler(
  filename = function() {
    paste("Voom_plot", input$VM_type, sep = ".")
  },
  content = function(file) {
    # Determine the appropriate device to use
    if (input$VM_type == "png") {
      png(file, width = input$VMw_1, height = input$VMh_1, units = "in", res = input$VM_dpi)
    } else if (input$VM_type == "pdf") {
      pdf(file, width = input$VMw_1, height = input$VMh_1)
    } else {
      stop("Unsupported file type selected.")
    }
    
    # Render the plot to the file
    tx = voom(dgeF,design,plot=TRUE)
    
    # Close the device
    dev.off()
  }
) # End of voom plot
# Mean-variance trend plot
# Fit model to data given design
fit = lmFit(vm, design)
fit = eBayes(fit)

# Mean-variance plot
# Reactive
meanv_p <- reactive({plotSA(fit, main="Final model: Mean-variance trend")})

# plot
output$meanv_plott <- renderPlot({meanv_p()})

# download
output$meanv_plottd <- downloadHandler(
  filename = function() {
    paste("Mean_variance_trend_plot", input$meanv_type, sep = ".")
  },
  content = function(file) {
    # Determine the appropriate device to use
    if (input$meanv_type == "png") {
      png(file, width = input$meanw_1, height = input$meanh_1, units = "in", res = input$meanv_dpi)
    } else if (input$meanv_type == "pdf") {
      pdf(file, width = input$meanw_1, height = input$meanh_1)
    } else {
      stop("Unsupported file type selected.")
    }
    
    # Render the plot to the file
    plotSA(fit, main="Final model: Mean-variance trend")
    
    # Close the device
    dev.off()
  }
) # End of Mean-variance trend plot
# DE table

# Define a reactive value to hold counts
reactive_counts <- reactiveValues()

# DE genes
DE_table <- reactive({
  cutoff <- input$dex1
  lfc <- input$dex2
  adj_method_DE <- input$dex3
  
  # Make a list of DE results based on provided criteria
  de_results <- limma::decideTests(fit, lfc = lfc, adjust.method = adj_method_DE)
  
  # Extract topTable results
  allgenes <- limma::topTable(fit, coef = 2, adjust.method = adj_method_DE, n = Inf)
  
  # Rearrange order of columns
  allgenesF <- allgenes[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  
  # Extract genes with absolute logFC > lfc
  allgenesFF <- allgenesF[abs(allgenesF$logFC) > lfc, ]
  
  # Extract significant DE genes and order from smallest to largest p-values
  allgenesFF <- allgenesFF[allgenesFF$adj.P.Val < cutoff, ]
  allgenesFF <- allgenesFF[order(allgenesFF$P.Value, allgenesFF$adj.P.Val), ]
  
  # Up and Down regulatory genes
  SIG_UPDW <- allgenesFF %>% 
    mutate(Gene_Expression = case_when(adj.P.Val < cutoff & logFC > lfc ~ "Up-regulated",
                                       adj.P.Val < cutoff & -logFC > lfc ~ "Down-regulated"))
  
  # Download Up and down regulatory genes
  output$SIG_UPDWd <- downloadHandler(
    filename = function() {
      paste("Up_and_Down_regulated_genes", ".csv", sep = "")
    },
    content = function(file) {
      # Write the current DE_table data to a CSV file
      write.csv(SIG_UPDW, file, row.names = TRUE)
    }
  )
  
  # Up regulated genes
  allgenes_UP <- allgenesFF[allgenesFF$logFC > lfc, ]

  # Download Up and down regulatory genes
  output$allgenes_UPd <- downloadHandler(
    filename = function() {
      paste("Up_regulated_genes", ".csv", sep = "")
    },
    content = function(file) {
      # Write the current DE_table data to a CSV file
      write.csv(allgenes_UP, file, row.names = TRUE)
    }
  )
  
  # Down regulated genes
  allgenes_DW <- allgenesFF[-(allgenesFF$logFC) > lfc, ]
  
  # Download Up and down regulatory genes
  output$allgenes_DWd <- downloadHandler(
    filename = function() {
      paste("Down_regulated_genes", ".csv", sep = "")
    },
    content = function(file) {
      # Write the current DE_table data to a CSV file
      write.csv(allgenes_DW, file, row.names = TRUE)
    }
  )
  
  # Update reactive values with counts
  reactive_counts$total_genes <- nrow(allgenesFF)
  reactive_counts$upregulated_genes <- nrow(allgenesFF[allgenesFF$logFC > lfc, ])
  reactive_counts$downregulated_genes <- nrow(allgenesFF[allgenesFF$logFC < -lfc, ])
  
  # Render information box with the count of significant DE genes
  output$DE_box <- renderInfoBox({
    total_genes <- nrow(allgenesFF)
    upregulated_genes <- nrow(allgenes_UP)
    downregulated_genes <- nrow(allgenes_DW)
    
    infoBox(
      title = "",
      value = paste0(
        "Total ", total_genes, " significantly differentially expressed genes were identified."), 
        paste0(
        "Of these, ", upregulated_genes, " are upregulated genes and ", downregulated_genes, " are downregulated genes."
      ),
      color = "purple",
      fill = TRUE
    )
  })
 
 return(allgenesFF)  # Ensure the result is returned
})

# Expose reactive counts as a reactive expression
counts <- reactive({
  reactive_counts
})

# DE genes table output
output$DE_genes <- DT::renderDataTable({
  DE_table()  # Call the reactive expression
}, rownames = TRUE, options = list(scrollX = TRUE))

# Download results
output$DE_resultsd <- downloadHandler(
  filename = function() {
    paste("Differential_expression_analysis_all_data", ".csv", sep = "")
  },
  content = function(file) {
    # Write the current DE_table data to a CSV file
    write.csv(DE_table(), file, row.names = TRUE)
  }
)

# Top DE
output$Top_DE <- downloadHandler(
  filename = function() {
    paste("Top_Differential_expression_analysis", ".csv", sep = "")
  },
  content = function(file) {
    # Write the current DE_table data to a CSV file
    write.csv(DE_table()[1:input$DE_top, ], file, row.names = TRUE)
  }
)

# Volcano plot

# Example set_col function
volcano_plot <- function(fit,
                         adj_method = "BH",
                         sig = "adjP",
                         cutoff = 0.05,
                         lfc = 1,
                         line_fc = TRUE,
                         line_p = TRUE,
                         palette = "viridis",
                         text_size = 10,
                         label_top = FALSE,
                         n_top = 10,
                         title = NULL) {
  
  # Extract the required data from the fit object to make the volcano plot
  res_de <- limma::topTable(fit, coef = 2, adjust.method = adj_method, n = Inf)
  
  # Identify significant genes
  res_de <- res_de %>%
    dplyr::mutate(Significant = if (sig == "P") {
      P.Value < cutoff & abs(logFC) > lfc
    } else {
      adj.P.Val < cutoff & abs(logFC) > lfc
    }) %>%
    dplyr::mutate(Label = ifelse(Significant, rownames(res_de), ""))
  
  # Determine the y-axis based on significance type
  res_de$y_value <- if (sig == "P") {
    -log10(res_de$P.Value)
  } else {
    -log10(res_de$adj.P.Val)
  }
  
  # Generate the color palette using the custom set_col function
  n_colors <- length(unique(res_de$Significant))
  colors <- set_col2(palette = input$vl_7, n_colors)
  
  # Create the volcano plot
  de_volcanoplot <- ggplot2::ggplot(res_de, aes(x = logFC, y = y_value, color = Significant, label = Label)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +  # Fixed size of points
    ggplot2::xlab(expression("log"[2] * " fold change")) +
    ggplot2::ylab(if (sig == "P") expression("-log"[10] * "(P-value)") else expression("-log"[10] * "(adj.P-value)")) +
    ggplot2::scale_color_manual(values = c("FALSE" = "tan4", "TRUE" = colors[2]), labels = c("Not Significant", "Significant")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(color = "black", size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(color = "black", size = 10, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "right",
      panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
      panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
      panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  # Add plot title if specified
  if (!is.null(title)) {
    de_volcanoplot <- de_volcanoplot + ggplot2::ggtitle(title)
  }
  
  # Add vertical and horizontal lines if required
  if (line_fc) {
    de_volcanoplot <- de_volcanoplot +
      ggplot2::geom_vline(xintercept = c(-lfc, lfc), color = "grey60", linetype = 2, size = 0.5, alpha = 0.8)
  }
  
  if (line_p) {
    de_volcanoplot <- de_volcanoplot +
      ggplot2::geom_hline(yintercept = -log10(cutoff), color = "grey60", linetype = 2, size = 0.5, alpha = 0.8)
  }
  
  # Add labels for top significant points if required
  if (label_top) {
    top_data <- res_de[order(-res_de$y_value)[1:n_top], ]
    de_volcanoplot <- de_volcanoplot +
      ggrepel::geom_text_repel(
        data = top_data,
        aes(label = Label),
        size = text_size / 4,
        color = "blue"
      )
  }
  
  return(de_volcanoplot)
}

#Volcano Plot
vol_plott <- reactive({
  volcano_plot(fit, adj_method = input$dex3, cutoff = input$dex1, lfc = input$dex2, n_top = input$vl_6, sig = input$vl_1, line_fc = input$vl_2, line_p = input$vl_3, label_top = input$vl_4, text_size = input$vl_5, palette = input$vl_7, title = input$vl_8)
  })
# Plot
output$volc_plot <- renderPlot({ vol_plott ()})
# Download plot
output$volc_plotd <- downloadHandler(
  filename = function() {
    paste("Volcano_plot", input$vl_type, sep = ".")
  },
  content = function(file) {
    if (input$vl_type == "png") {
      ggsave(file, plot = vol_plott (), device = "png", width = input$vl_w, height = input$vl_h, units = "in", dpi = input$vl_dpi, bg = "white")
    } else if (input$vl_type == "pdf") {
      ggsave(file, plot = vol_plott (), device = "pdf", width = input$vl_w, height = input$vl_h)
    } else {
      stop("Unsupported file type selected.")
    }
  }
) # End of volcano plot
# Heat map 
# Heatmap with cli
get_palette1x <- reactive({
  palette_namex <- input$ht_6x
  n_colorsx <- 50  # You can adjust this or make it dynamic
  set_col3(palette_namex, n_colorsx)
})
heat_px <- reactive({
  clin_data <- clin_data
  pht_normx <- data.frame(vm$E)
  heat_dex <- pht_normx[rownames(DE_table())[input$ht_1x:input$ht_2x],]
  rownames(clin_data) <- colnames(heat_dex)
  output$heat_gtx <- renderText({ nrow(DE_table()) })
  output$heat_clx <- renderText({ ncol(heat_dex) })
  if (input$ht_20x == TRUE){
    pheatmap(heat_dex[input$ht_3x:input$ht_4x], main = input$ht_5x, color = get_palette1x(), scale = input$ht_7x, border_color = input$ht_8x, cluster_rows=input$ht_9x, cluster_cols=input$ht_10x, legend = input$ht_11x, show_colnames = input$ht_12x,  show_rownames = input$ht_13x, fontsize = input$ht_14x, annotation_col= clin_data)
  }else
    pheatmap(heat_dex[input$ht_3x:input$ht_4x], main = input$ht_5x, color = get_palette1x(), scale = input$ht_7x, border_color = input$ht_8x, cluster_rows=input$ht_9x, cluster_cols=input$ht_10x, legend = input$ht_11x, show_colnames = input$ht_12x,  show_rownames = input$ht_13x, fontsize = input$ht_14x)
})

# Plot
output$heat_plotx <- renderPlot({ heat_px ()})

# Download plot
output$heat_plotdx <- downloadHandler(
  filename = function() {
    paste("Heatmap_plot", input$ht_typex, sep = ".")
  },
  content = function(file) {
    if (input$ht_typex == "png") {
      ggsave(file, plot = heat_px (), device = "png", width = input$ht_wx, height = input$ht_hx, units = "in", dpi = input$ht_dpix, bg = "white")
    } else if (input$ht_typex == "pdf") {
      ggsave(file, plot = heat_px (), device = "pdf", width = input$ht_wx, height = input$ht_hx)
    } else {
      stop("Unsupported file type selected.")
    }
  }
)


### ML
## ----------- Pre-process gene expression data for modeling ----------------
  # Get Sample name & group or condition information
  ML_group <- exp_design
  
  # Get Normalized genes (E) data to compare with fit2 DE genes
  ML_norm <- data.frame(vm$E)
  
  #Get norm log DE of top fit2 DEG
  DE_gene_list <- ML_norm[rownames(DE_table()),]
  
  # Transpose data frame
  topint_trans <- as.data.frame(t(DE_gene_list))
  
  # Remove sample names.
  rownames(topint_trans) <- NULL
  
  # Add a new column with the group or condition information.
  topint_trans$group <- ML_group$Group_Type
  topint_trans_1 <- topint_trans
  
  # Convert group names to R compatible names
  topint_trans_1$group <- make.names(topint_trans_1$group)
  
  # Convert group to a factor (important for varimp calculations)
  topint_trans_1$group <- factor(topint_trans_1$group)
  
  model_df <- topint_trans_1
  
  ## Split the data frame into training and test data sets
  split_df <- split_data(model_df = model_df, train_size = input$gms1, seed = input$gms2)
  
  # Create a model_list object by training models on the training data set using a custom list of algorithms.Don't forget to fix the random seed for reproducibility.
  ## Fit models using a user-specified list of ML algorithms.
  model_list <- train_models(
    split_df,
    resample_method = input$gmgmt1,
    resample_iterations = input$gmgmt3,
    num_repeats = input$gmgmt2,
    algorithm_list = input$gmgmt5,
    seed = input$gmgmt4)
  
  # Get model results
  # Model performance plot- DOT PLOT
  model_results_mat <- reactive({extract_and_sort_metrics(model_list)})
  
  output$model_results_matd <- downloadHandler(
    filename = function() {
      paste("Model_performance_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(model_results_mat(), file, row.names = FALSE)
    }
  )
  
  # Model performance plot- DOT PLOT
  MLPR_plot <- reactive({
    perf_plot(model_list, ggtitle = input$MLPR_ttl, text_size = input$MLPR_txt, palette = input$MLPR_plt)
  })
  
  # Performance Plots
  output$MLPR_perf_plot <- renderPlot({MLPR_plot ()})
  
  # Download plot
  output$MLPR_perf_plotd <- downloadHandler(
    filename = function() {
      paste("ML_performance_plot", input$MLPR_type, sep = ".")
    },
    content = function(file) {
      if (input$MLPR_type == "png") {
        ggsave(file, plot = MLPR_plot (), device = "png", width = input$MLPR_w, height = input$MLPR_h, units = "in", dpi = input$MLPR_dpi, bg = "white")
      } else if (input$MLPR_type == "pdf") {
        ggsave(file, plot = MLPR_plot (), device = "pdf", width = input$MLPR_w, height = input$MLPR_h)
      } else {
        stop("Unsupported file type selected.")
      }
    }
  ) # End of performance_plot
  
  # Variable_importance_plot
  varimp_plotF <- function(model_list,
                           ...,
                           type = "lollipop",
                           text_size = 10,
                           palette = "viridis",
                           n_row,
                           n_col) {
    
    # Binding for global variable
    genes <- Importance <- NULL
    
    n_row <- input$VIMP_rw
    n_col <- input$VIMP_clm
    
    # Calculate variable importance with VarImp for each ML algorithm-based
    # model in the list
    vimp <- lapply(
      model_list,
      function(x) {
        tryCatch(caret::varImp(
          x,
          ...
        ),
        error = function(e) NULL
        )
      }
    )
    
    # Drop Null items from the list
    vimp <- Filter(Negate(is.null), vimp)
    
    # Make necessary changes to 'importance data frames' in the list before
    # plotting and select top 10 variables
    plot_imp_data <- lapply(
      vimp,
      function(x) {
        x <- x$importance[1]
        x$genes <- rownames(x)
        rownames(x) <- NULL
        colnames(x) <- c("Importance", "genes")
        x$genes <- gsub("`|\\\\", "", x$genes)
        x$Importance[is.nan(x$Importance)] <- NA
        x <- x[rowSums(is.na(x)) == 0, ]
        x <- head(x, input$VIMP_topimp)  # Select top 10
        x
      }
    )
    
    # Drop empty data frames before proceeding
    plot_imp_data <- Filter(function(x) dim(x)[1] > 0, plot_imp_data)
    
    # Make plots based on type
    if (type == "bar") {
      # Make bar plots
      vi_plots <- lapply(seq_along(plot_imp_data), function(i) {
        t <- plot_imp_data[[i]]
        ggplot2::ggplot(
          data = t,
          aes(
            x = reorder(genes, Importance),
            y = Importance,
            fill = Importance
          )
        ) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::coord_flip() +
          viridis::scale_fill_viridis(
            option = palette,
            direction = -1
          ) +
          xlab("") +
          extra_facet_theme() +
          ggplot2::theme(
            plot.title = element_text(
              size = text_size * 1.2,
              face = "bold",
              hjust = 0.5
            ),
            plot.subtitle = element_text(
              size = 10,
              face = "bold",
              hjust = 0.5
            ),
            plot.margin = ggplot2::margin(10, 10, 10, 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
            legend.key.width = grid::unit(0.8, "cm"),
            legend.key.height = grid::unit(0.2, "cm"),
            legend.title = element_text(size = text_size * 0.7),
            legend.text = element_text(size = text_size * 0.5),
            axis.text.y = element_text(
              size = text_size * 0.8,
              face = "bold"
            ),
            panel.grid.major.x = element_line(
              size = 0.1,
              color = "grey80"
            )
          ) +
          ggplot2::guides(fill = guide_colorbar(
            title.position = "top"
          )) +
          ggtitle(input$VIMP_ttl) +
          labs(subtitle = model_list[[i]]$modelInfo$label)
      })
    } else {
      # Make lollipop plots
      vi_plots <- lapply(seq_along(plot_imp_data), function(i) {
        t <- plot_imp_data[[i]]
        ggplot2::ggplot(
          data = t,
          aes(
            x = reorder(genes, Importance),
            y = Importance,
            color = Importance
          )
        ) +
          ggplot2::geom_segment(aes(
            xend = genes,
            y = 0,
            yend = Importance
          ),
          lwd = text_size * 0.2
          ) +
          ggplot2::geom_point(aes(color = Importance),
                              size = text_size,
                              show.legend = FALSE
          ) +
          ggplot2::coord_flip() +
          ggplot2::geom_label(
            label = round(t$Importance, digits = 1),
            fill = NA,
            label.size = NA,
            colour = "white",
            size = text_size * 0.3
          ) +
          ggplot2::xlab("") +
          viridis::scale_color_viridis(
            option = palette,
            direction = -1
          ) +
          extra_facet_theme() +
          ggplot2::theme(
            plot.title = element_text(
              size = text_size * 1.2,
              face = "bold",
              hjust = 0.5
            ),
            plot.subtitle = element_text(
              size = 10,
              face = "bold",
              hjust = 0.5
            ),
            plot.margin = ggplot2::margin(10, 10, 10, 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
            legend.key.width = grid::unit(0.8, "cm"),
            legend.key.height = grid::unit(0.2, "cm"),
            legend.title = element_text(size = text_size * 0.7),
            legend.text = element_text(size = text_size * 0.5),
            axis.text.y = element_text(
              size = text_size * 0.8,
              face = "bold"
            ),
            axis.text.x = element_text(
              size = text_size * 0.9
            ),
            panel.grid.major.x = element_line(
              size = 0.1,
              color = "grey80"
            )
          ) +
          ggplot2::guides(colour = guide_colourbar(
            title = "Importance",
            title.position = "top",
            title.hjust = 0.5
          )) +
          ggtitle(input$VIMP_ttl) +
          labs(subtitle = model_list[[i]]$modelInfo$label)
      })
    }
    
    grid.arrange(grobs = vi_plots, nrow = n_row, ncol = n_col, newpage = TRUE)
  }
  
  # Var Imp Plots
  VIMP_plot <- reactive({
    varimp_plotF(model_list, type = input$VIMP_ptype, text_size = input$VIMP_txt, palette = input$VIMP_plt)
  })
  
  # Performance Plots
  output$VIMP_plott <- renderPlot({VIMP_plot ()})
  
  # Download plot
  output$VIMP_plottd <- downloadHandler(
    filename = function() {
      paste("ML_Variable_importance_plot", input$VIMP_type, sep = ".")
    },
    content = function(file) {
      if (input$VIMP_type == "png") {
        ggsave(file, plot = VIMP_plot (), device = "png", width = input$VIMP_w, height = input$VIMP_h, units = "in", dpi = input$VIMP_dpi, bg = "white")
      } else if (input$VIMP_type == "pdf") {
        ggsave(file, plot = VIMP_plot (), device = "pdf", width = input$VIMP_w, height = input$VIMP_h)
      } else {
        stop("Unsupported file type selected.")
      }
    }
  ) # End of Variable_importance_plot
  
  # Test models - prob
  prob_list <- test_models(model_list, split_df, type = "prob")
  
  # Confusion matrix
  cnf_mat <- reactive({
    test_models(model_list, split_df, type = "raw", save_confusionmatrix = FALSE)
  })
  
  output$ROCML_cnf_table <- downloadHandler(
    filename = function() {
      paste("ML_Confusion_matrix", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cnf_mat(), file, row.names = FALSE)
    }
  )
  
  # ROC PLOT
  # Reactive value to trigger update
  update_palette <- reactiveVal(NULL)
  # Reactive 
  ROCML_plot <- reactive({
    # Update palette to force reactivity
    update_palette(input$ROCML_plt)
    roc_plot(prob_list, split_df, title = input$ROCML_ttl, text_size = input$ROCML_txt, palette = input$ROCML_plt, multiple_plots = input$ROCML_mlt)
  })
  
  # ROC Plots
  output$ROCML_plott <- renderPlot({ROCML_plot ()})
  
  # Download plot
  output$ROCML_plottd <- downloadHandler(
    filename = function() {
      paste("ML_ROC_plot", input$ROCML_type, sep = ".")
    },
    content = function(file) {
      if (input$ROCML_type == "png") {
        ggsave(file, plot = ROCML_plot (), device = "png", width = input$ROCML_w, height = input$ROCML_h, units = "in", dpi = input$ROCML_dpi, bg = "white")
      } else if (input$ROCML_type == "pdf") {
        ggsave(file, plot = ROCML_plot (), device = "pdf", width = input$ROCML_w, height = input$ROCML_h)
      } else {
        stop("Unsupported file type selected.")
      }
    }
  ) # End of ROC plot
# END of ML option

# ORA
# Dynamically load the chosen organism database
library(input$organismDb, character.only = TRUE)
annDb <- get(input$organismDb)
keytypes <- keytypes(annDb)
updateSelectInput(session, "keytype", choices = keytypes, selected = input$keytype)

# Reactive values
myValues <- reactiveValues()

# ORA for DE only
output$GOTOPH_DE <- renderText({ nrow(DE_table()) })

# Observer to handle button clicks and perform actions ORA for DE only
observeEvent(input$update_ORA1, {
  req(input$GOTOPH)
  G11 <- gsub("\\.[0-9]*$", "", rownames(DE_table()))
  G1 <- G11[1:input$GOTOPH]
  req(G1)
  
  shinyalert(
    title = '',
    text = '<h3 style="text-align: center; color: #000000;">The GO analysis is in progress. Please wait...</h3>',
    type = "info",
    html = TRUE,
    closeOnClickOutside = TRUE,
    closeOnEsc = TRUE
  )
  
  # Perform GO enrichment analysis
  GS <- enrichGO(
    G1, 
    OrgDb = input$organismDb, 
    keyType = input$keytype, 
    ont = input$ontt, 
    pvalueCutoff = input$pvalueCutoff,
    qvalueCutoff = input$qvalueCutoff,
    pAdjustMethod = input$pAdjustMethod,
    minGSSize = input$minGSSize, 
    maxGSSize = input$maxGSSize
  )
  
  # GO Outputs
  output$ORA_GO_table <- DT::renderDataTable({
    as.data.frame(GS)
  }, rownames = TRUE, options = list(scrollX = TRUE))
  # GO download
  output$ORA_GO_tabled <- downloadHandler(
    filename = function() {
      paste("GO_over_representation_results", ".csv", sep = "")
    },
    content = function(file) {
      # Write the current DE_table data to a CSV file
      write.csv(as.data.frame(GS), file, row.names = TRUE)
    }
  )
  
  # Bar plot
  ORA_barp <- reactive({
    barplot(GS, showCategory = input$GOB_ct, font.size = input$GOB_fz, title = input$GOB_pt) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
        axis.title.x = ggplot2::element_text(color = "black", size = 10, face = "bold"),
        panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
        panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
        panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
  })
  
  # GO bar plot
  output$ORA_barplot <- renderPlot({ORA_barp()})
  # GO bar plot download
  output$ORA_barplotd <- downloadHandler(
    filename = function() {
      paste("GO_bar_plot", input$GOB_type, sep = ".")
    },
    content = function(file) {
      if (input$GOB_type == "png") {
        png(file, width = input$GOB_w, height = input$GOB_h, units = "in", res = input$GOB_dpi)
      } else if (input$GOB_type == "pdf") {
        pdf(file, width = input$GOB_w, height = input$GOB_h)
      } else {
        stop("Unsupported file type selected.")
      }
      
      # Render the plot
      print(ORA_barp())  # Call the reactive plot function
      
      # Close the device
      dev.off()
    }
  )
  # DOT PLOT
  ORA_dotp <- reactive({
    # Initialize the plot based on input$ontt
    if (input$ontt == "ALL") {
      D1 <- dotplot(GS, 
                    showCategory = input$DOT_ct, 
                    font.size = input$DOT_fz, 
                    title = input$DOT_pt, 
                    split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scale = "free")
    } else {
      D1 <- dotplot(GS, 
                    showCategory = input$DOT_ct, 
                    font.size = input$DOT_fz, 
                    title = input$DOT_pt)
    }
    
    # Add ggplot2 theme customizations
    D1 + ggplot2::theme(
      plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(color = "black", size = 10, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
      panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
      panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  })
  
  #GO dot plot
  output$ORA_dotplot <- renderPlot(ORA_dotp())
  
  # GO dot plot download
  output$ORA_dotplotd <- downloadHandler(
    filename = function() {
      paste("GO_dot_plot", input$DOT_type, sep = ".")
    },
    content = function(file) {
      if (input$DOT_type == "png") {
        png(file, width = input$DOT_w, height = input$DOT_h, units = "in", res = input$DOT_dpi)
      } else if (input$DOT_type == "pdf") {
        pdf(file, width = input$DOT_w, height = input$DOT_h)
      } else {
        stop("Unsupported file type selected.")
      }
      
      # Render the plot
      print(ORA_dotp())  # Call the reactive plot function
      
      # Close the device
      dev.off()
    }
  ) # End of GO dot plot
  # CNET plot
  ORA_CNTp <- reactive({
    GS2 <- simplify(GS)
    cnetplot(GS2, color.params = list(foldChange = GS), showCategory = input$CNT_ct, colorEdge = input$CNT_fz, circular = input$CNT_pt) +
      ggtitle(input$CNT_ttl) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
        panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
        panel.grid.minor = ggplot2::element_blank()  # Remove minor grid lines
      )
  })
  #GO CNT plot
  output$ORA_CNTplot <- renderPlot(ORA_CNTp())
  
  # GO CNT plot download
  output$ORA_CNTplotd <- downloadHandler(
    filename = function() {
      paste("GO_Gene_Concept_Network_plot", input$CNT_type, sep = ".")
    },
    content = function(file) {
      if (input$CNT_type == "png") {
        png(file, width = input$CNT_w, height = input$CNT_h, units = "in", res = input$CNT_dpi)
      } else if (input$CNT_type == "pdf") {
        pdf(file, width = input$CNT_w, height = input$CNT_h)
      } else {
        stop("Unsupported file type selected.")
      }
      
      # Render the plot
      print(ORA_CNTp())  # Call the reactive plot function
      
      # Close the device
      dev.off()
    }
  ) # End of GO CNT plot
  # Heat plot
  ORA_HGTp <- reactive({
    GS22 <- simplify(GS)
    enrichplot::heatplot(GS, showCategory = input$HGT_ct) +
      ggtitle(input$HGT_ttl) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
        panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
        panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
  })
  #GO HGT plot
  output$ORA_HGTplot <- renderPlot(ORA_HGTp())
  
  # GO HGT plot download
  output$ORA_HGTplotd <- downloadHandler(
    filename = function() {
      paste("GO_Heat_plot", input$HGT_type, sep = ".")
    },
    content = function(file) {
      if (input$HGT_type == "png") {
        png(file, width = input$HGT_w, height = input$HGT_h, units = "in", res = input$HGT_dpi)
      } else if (input$HGT_type == "pdf") {
        pdf(file, width = input$HGT_w, height = input$HGT_h)
      } else {
        stop("Unsupported file type selected.")
      }
      
      # Render the plot
      print(ORA_HGTp())  # Call the reactive plot function
      
      # Close the device
      dev.off()
    }
  ) # End of GO Heat plot
  # Enrichment Map plot
  ORA_ENRTp <- reactive({
    GSF <- pairwise_termsim(GS)
    emapplot(GSF, showCategory = input$ENRT_ct) +
      ggtitle(input$ENRT_ttl) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
        panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
      )
  })
  #GO ENRT plot
  output$ORA_ENRTplot <- renderPlot(ORA_ENRTp())
  
  # GO ENRT plot download
  output$ORA_ENRTplotd <- downloadHandler(
    filename = function() {
      paste("GO_Enrichment_Map_plot", input$ENRT_type, sep = ".")
    },
    content = function(file) {
      if (input$ENRT_type == "png") {
        png(file, width = input$ENRT_w, height = input$ENRT_h, units = "in", res = input$ENRT_dpi)
      } else if (input$ENRT_type == "pdf") {
        pdf(file, width = input$ENRT_w, height = input$ENRT_h)
      } else {
        stop("Unsupported file type selected.")
      }
      
      # Render the plot
      print(ORA_ENRTp())  # Call the reactive plot function
      
      # Close the device
      dev.off()
    }
  ) # End of GO Enrichment Map plot
  # UpSet Plot
  ORA_UPSp <- reactive({
    GS1FF <- simplify(GS)
    upsetplot(GS1FF, n = input$UPS_ct) +
      ggtitle(input$UPS_ttl) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white"),  # Set panel background to white
        panel.grid.major = ggplot2::element_blank(),  # Remove major grid lines
        panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
  })
  #GO UPS plot
  output$ORA_UPSplot <- renderPlot(ORA_UPSp())
  
  # GO UPS plot download
  output$ORA_UPSplotd <- downloadHandler(
    filename = function() {
      paste("GO_UpSet_plot", input$UPS_type, sep = ".")
    },
    content = function(file) {
      if (input$UPS_type == "png") {
        png(file, width = input$UPS_w, height = input$UPS_h, units = "in", res = input$UPS_dpi)
      } else if (input$UPS_type == "pdf") {
        pdf(file, width = input$UPS_w, height = input$UPS_h)
      } else {
        stop("Unsupported file type selected.")
      }
      
      # Render the plot
      print(ORA_UPSp())  # Call the reactive plot function
      
      # Close the device
      dev.off()
    }
  ) # End of GO UpSet Plot
  
  # KEGG Enrichment Analysis
  organismsDbKegg <- c(
    "org.Hs.eg.db" = "hsa", "org.Mm.eg.db" = "mmu", "org.Rn.eg.db" = "rno",
    "org.Sc.sgd.db" = "sce", "org.Dm.eg.db" = "dme", "org.At.tair.db" = "ath",
    "org.Dr.eg.db" = "dre", "org.Bt.eg.db" = "bta", "org.Ce.eg.db" = "cel",
    "org.Gg.eg.db" = "gga", "org.Cf.eg.db" = "cfa", "org.Ss.eg.db" = "ssc",
    "org.EcK12.eg.db" = "eck", "org.Xl.eg.db" = "xla",
    "org.Pt.eg.db" = "ptr", "org.Ag.eg.db" = "aga","org.EcSakai.eg.db" = "ecs"
  )
  
  K1 <- bitr(G1, fromType = input$keytype, toType = "ENTREZID", OrgDb = input$organismDb)
  
  # Reactive expression to load KEGG results
  KEGG_results <- reactive({
    req(input$organismDb)  # Ensure input$organismDb is available
    
    kegg_genes <- K1$ENTREZID
    kegg_result <- enrichKEGG(
      gene = kegg_genes,
      organism = organismsDbKegg[[input$organismDb]],
      pvalueCutoff = input$pvalueCutoff,
      qvalueCutoff = input$qvalueCutoff,
      keyType = "kegg",
      minGSSize = input$minGSSize,
      maxGSSize = input$maxGSSize
    )
    
    # Debugging output
    print("KEGG Results:")
    print(head(kegg_result$ID))  # Print the first few pathway IDs
    kegg_result
  })
  
  # KEGG Results
  output$ORA_KEGG_table <- DT::renderDataTable({
    as.data.frame(KEGG_results())
  }, rownames = TRUE, options = list(scrollX = TRUE))
  # KEGG results download
  output$ORA_KEGG_tabled <- downloadHandler(
    filename = function() {
      paste("KEGG_results", ".csv", sep = "")
    },
    content = function(file) {
      # Write the current DE_table data to a CSV file
      write.csv(as.data.frame(KEGG_results()), file, row.names = TRUE)
    }
  )
  
  # Update pathway ID choices dynamically
  observe({
    kegg_result <- KEGG_results()
    if (!is.null(kegg_result) && length(kegg_result$ID) > 0) {
      available_pathways <- as.character(kegg_result$ID)
      print("Updating pathway selector with IDs:")
      print(available_pathways)
      
      updateSelectizeInput(session, "pathwayIds", choices = available_pathways, selected = NULL)
    } else {
      updateSelectizeInput(session, "pathwayIds", choices = NULL, selected = NULL)
      output$pathview_message <- renderUI({
        h3("No pathways available for the selected organism.")
      })
      output$pathview_plot <- renderImage({
        list(
          src = "",
          filetype = "image/png",
          alt = "No Pathview Plot Available"
        )
      })
    }
  })
  
  # Pathview
  # Reactive expression to generate Pathview plot
  pathviewReactive <- eventReactive(input$generatePathview, {
    req(input$pathwayIds)  # Ensure pathway ID is selected
    
    withProgress(message = "Plotting Pathview ...", {
      isolate({
        tryCatch({
          gene_data <- K1
          gene_list <- gene_data$ENTREZID
          names(gene_list) <- gene_data$SYMBOL
          
          # Define file paths in the current working directory
          png_file <- paste0(input$pathwayIds, ".pathview.png")
          pdf_file <- paste0(input$pathwayIds, ".pathview.pdf")
          temp_image <- tempfile(fileext = ".png")  # Temporary file to store image
          
          setProgress(value = 0.3, detail = paste0("Generating Pathview for ID ", input$pathwayIds, " ..."))
          
          # Plot and save PNG
          pathview(gene.data = gene_list, pathway.id = input$pathwayIds, species = organismsDbKegg[[input$organismDb]])
          if (file.exists(png_file)) {
            file.copy(png_file, temp_image, overwrite = TRUE)
          } else {
            showNotification("PNG file not generated.")
            return(NULL)
          }
          
          setProgress(value = 0.7, detail = paste0("Generating PDF for ID ", input$pathwayIds, " ..."))
          
          # Plot and save PDF
          pathview(gene.data = gene_list, pathway.id = input$pathwayIds, species = organismsDbKegg[[input$organismDb]], kegg.native = FALSE)
          if (!file.exists(pdf_file)) {
            showNotification("PDF file not generated.")
          }
          
          # Return the image list for rendering
          return(list(
            src = temp_image,          # Path to the temporary image file
            contentType = "image/png", # Content type for the image
            alt = "Pathview image",    # Alt text for the image
            deleteFile = TRUE          # Request Shiny to delete the file after serving
          ))
        }, error = function(e) {
          showNotification(paste("Error generating Pathview plot:", e$message))
          return(NULL)
        })
      })
    })
  })
  
  # Render Pathview plot
  output$pathview_plot <- renderImage({
    req(pathviewReactive())  # Ensure the reactive expression is available
    pathviewReactive()
  }, deleteFile = TRUE)
  
  # Reactive expression to check if plots are available
  output$pathviewPlotsAvailable <- reactive({
    !is.null(pathviewReactive())
  })
  
  outputOptions(output, "pathviewPlotsAvailable", suspendWhenHidden = FALSE)
  
  # Download handlers for Pathview plot files
  output$downloadPathviewPng <- downloadHandler(
    filename = function() {
      paste0(input$pathwayIds, ".pathview.png")
    },
    content = function(file) {
      png_file <- paste0(input$pathwayIds, ".pathview.png")
      if (file.exists(png_file)) {
        file.copy(png_file, file, overwrite = TRUE)
      } else {
        showNotification("PNG file not found.")
      }
      # Ensure file is deleted after download
      unlink(png_file)
    },
    contentType = "image/png"
  )
  
  output$downloadPathviewPdf <- downloadHandler(
    filename = function() {
      paste0(input$pathwayIds, ".pathview.pdf")
    },
    content = function(file) {
      pdf_file <- paste0(input$pathwayIds, ".pathview.pdf")
      if (file.exists(pdf_file)) {
        file.copy(pdf_file, file, overwrite = TRUE)
      } else {
        showNotification("PDF file not found.")
      }
      # Ensure file is deleted after download
      unlink(pdf_file)
    },
    contentType = "application/pdf"
  )
  # END of pathview

# Simulate fetching data (replace with your actual data fetching code)
  Sys.sleep(30)  # Simulate a delay of 3 seconds
  # Once data is fetched, close the alert
  closeAlert()
})

# RMARKDOWN_REPORT
total_genes <- reactive({ counts()$total_genes })
upregulated_genes <- reactive({ counts()$upregulated_genes })
downregulated_genes <- reactive({ counts()$downregulated_genes })

# RMARKDOWN_ML
output$RMD_Report <- downloadHandler(
  filename = function() {
    paste("GENOXA_report_", Sys.Date(), ".", switch(input$RMD_format,
                                                      pdf_document = "pdf",
                                                      word_document = "docx",
                                                      html_document = "html"), sep = "")
  },
  content = function(file) {
    # Define the path to the temporary Rmd file
    tempReport <- file.path(tempdir(), "R1.Rmd")
    file.copy("R2.Rmd", tempReport, overwrite = TRUE)
    
    # Set up parameters to pass to the Rmd document
    params <- list(
      dp1 = input$dp1,
      dp2 = input$dp2,
      dp3 = input$dp3,
      refg = input$refg,
      dex1 = input$dex1,
      dex2 = input$dex2,
      dex3 = input$dex3,
      organismDb = input$organismDb,
      keytype = input$keytype,
      ontt = input$ontt,
      minGSSize = input$minGSSize,
      maxGSSize = input$maxGSSize,
      pvalueCutoff = input$pvalueCutoff,
      qvalueCutoff = input$qvalueCutoff,
      pAdjustMethod = input$pAdjustMethod,
      r1 = r1,
      samp1 = samp1,
      Filt_counts_1 = Filt_counts_1,
      filt2c = filt2c,
      total_genes= total_genes(),
      upregulated_genes = upregulated_genes(),
      downregulated_genes = downregulated_genes(),
      gms1 = input$gms1,
      gms2 = input$gms2,
      gmgmt1 = input$gmgmt1,
      gmgmt3 = input$gmgmt3,
      gmgmt4 = input$gmgmt4,
      gmgmt5 = input$gmgmt5
    )
    
    # Render the R Markdown document with the selected output format
    tryCatch({
      rmarkdown::render(
        input = tempReport,
        output_file = file,
        output_format = switch(input$RMD_format,
                               pdf_document = rmarkdown::pdf_document(),
                               word_document = rmarkdown::word_document(),
                               html_document = rmarkdown::html_document()),
        params = params,
        envir = new.env(parent = globalenv())
      )
    }, error = function(e) {
      message("Error rendering report: ", e$message)
    })
  }
)
# END OF RMARKDOWN_with_ML

# END of main code

# Simulate fetching data (replace with your actual data fetching code)
Sys.sleep(50)  # Simulate a delay of 3 seconds
# Once data is fetched, close the alert
closeAlert()

# Stop the cluster
stopCluster(cl22)

}) # End of Observe

# Reset app on button click
observeEvent(input$reset_button, {
  session$reload()
})

} #End of server code

# Run the application
shinyApp(ui, server, options = list(launch.browser = TRUE))

