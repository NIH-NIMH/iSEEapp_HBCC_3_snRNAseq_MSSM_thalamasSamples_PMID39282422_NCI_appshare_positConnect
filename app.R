# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
################################################################################
# Install in sequence
# R version 4.5.2

# Check package versions
# BiocManager::version()
# packageVersion("iSEE")
# packageVersion("SparseArray")
# packageVersion("S4Arrays")
# packageVersion("SingleCellExperiment")
# packageVersion("SummarizedExperiment")
# 
# # LOAD PACKAGES
library(iSEE); #BiocManager::install("iSEE",ask = FALSE,update = FALSE, force = TRUE)
library(ggplot2) #install.packages("ggplot2")
library(memoise) #install.packages("memoise")
library(htmltools) #install.packages("htmltools")
library(base64enc) #install.packages("base64enc")
################################################################################
# Set working directory: 
# install.packages("this.path")
library(this.path)
setwd(this.path::here())
################################################################################
# Create RENV file
#install.packages("renv")
# library(renv)
# renv::snapshot()
# Select "y" after this
# Then select "1"
################################################################################

################################################################################
# IMPORTANT: Use this script for any scRNA seq study.
# Requirements: 
# 1. Update variable names e.g. nCounts_RNA, nFeature_RNA, cclass, ctypes
#nCounts_RNA = "nCount_RNA"
#nFeature_RNA = "nFeature_RNA"
#cclass = "cclass"
#ctypes = "ctypes"
# 3. Update specific names of the genes that need plot in feature plot and heatmap
defaultGene = "GAD1"
heatMapGenes = "TCF7L2\nSLC17A6\nDRD2\nGAD1\nRELN\nSOX14\nGFAP\nFGFR3\nST18\nMOG\nOLIG1\nPDGFRA\nDNAH3\nFOXJ1\nTTR\nHTR2C\nC3\nCX3CR1\nPTPRC\nVWF\nPDGFRB\nDCN"
# 4. Update the Title PanelId = <number> and Footer PanelId = <number>
appNum = 3L
# 5. Update publication title and Pubmed ID
pubRef = "https://pubmed.ncbi.nlm.nih.gov/39282422/"
studyNam = "A conserved spectrum of cell types across the human mediodorsal and paraventricular thalamus (MSSM samples)"
# 6. Update and read "tour.txt" file
tour <- read.delim("tour.txt", sep=";", stringsAsFactors = FALSE, row.names = NULL)
################################################################################

################################################################################
# Install in sequence
# R version 4.5.2

# Check package versions
# BiocManager::version()
# packageVersion("iSEE")
# packageVersion("SparseArray")
# packageVersion("S4Arrays")
# packageVersion("SingleCellExperiment")
# packageVersion("SummarizedExperiment")
# 
# # LOAD PACKAGES
library(iSEE); #BiocManager::install("iSEE",ask = FALSE,update = FALSE, force = TRUE)
library(ggplot2) #install.packages("ggplot2")
library(memoise) #install.packages("memoise")
library(htmltools) #install.packages("htmltools")
library(base64enc) #install.packages("base64enc")
################################################################################
# Create RENV file
#install.packages("renv")
#library(renv)
#renv::snapshot()
# Select "y" after this
# Then select "1"
################################################################################
################################################################################
# To retrieve an option: "getOption('timeout')" and to set an option
options(timeout=600)
# Speed up app launching time by using Memoise with cached SCE
# 1) Memoised loader for your SCE
get_sce <- memoise(function() {
  # OR read input from local 
  sce <- readRDS("sce_mssm.RDS")
  sce
})
# 2) Use the cached SCE to build the app
sce_small <- get_sce()
################################################################################
################################################################################
# Specify number of colors for each cell type
library(RColorBrewer)
numColor = as.numeric(length(unique(sce_small$ctypes)))
n <- numColor
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sample(col_vector, n)
names(col_vector) <- as.vector(unique(sce_small$ctypes))
############################################################################################
# Create Custom plot for title
############################################################################################
# Set variables
logo_path_t =  "nimh-logo.png"      # optional: path to a PNG logo
img_t <- png::readPNG("nimh-logo.png")
grob_t <- grid::rasterGrob(img_t, interpolate = TRUE)
# Create a function to use with custom plot
TITLE_FUN <- function(se, rows, columns) {
  
  ggplot() +
    geom_rect(
      aes(xmin = 0, xmax = 100, ymin = 0, ymax = 10),
      fill = "white",
      color = NA
    ) +
    coord_cartesian(
      xlim = c(0, 100),
      ylim = c(0, 10),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
}

TitlePanel <- iSEE::createCustomPlot(
  FUN       = TITLE_FUN,
  restrict  = NULL,
  className = "TitlePanel",
  fullName  = "iSEE App ID = HBCC"
)

# Subclass one panel type (e.g., ReducedDimensionPlot) and hide its Data box
setClass("Title", contains = "TitlePanel")

setMethod(".hideInterface", "TitlePanel",
          function(x, field) {
            if (field == "DataBoxOpen" | field == "SelectionBoxOpen") return(TRUE)  # hide the whole Data parameters box
            callNextMethod()
          }
)

# Add custom output method to include hyperlink below the plot
setMethod(".defineOutput", "TitlePanel", function(x) {
  
  tagList(
    
    # ============================================================
    # HEADER
    # ============================================================
    
    tags$div(
      style = "
    display: flex;
    align-items: center;
    gap: 20px;
  ",
      
      # Clickable NIMH logo
      tags$a(
        href = "https://www.nimh.nih.gov/",
        target = "_blank",
        title = "National Institute of Mental Health",
        
        tags$img(
          src = paste0("data:image/png;base64,", logo1_base64),
          height = "75px",
          alt = "National Institute of Mental Health logo",
          style = "display: block;"
        )
      ),
      
      # Clickable HBCC title
      tags$a(
        href = "https://www.nimh.nih.gov/research/research-conducted-at-nimh/research-areas/research-support-services/hbcc",
        target = "_blank",
        title = "Human Brain Collection Core",
        
        "Human Brain Collection Core Data Portal at NIMH",
        
        style = "
      color: royalblue4;
      font-size: 52px;
      font-weight: bold;
      text-decoration: none;
    "
      )
    ),
    
    # ============================================================
    # STUDY NAME
    # ============================================================
    
    tags$div(
      style = "
        padding: 10px 20px;
        text-align: left;
        font-size: 30px;
        color: royalblue4;
      ",
      
      tags$span(
        tags$strong("Study Name:"),
        " ",
        tags$a(
          href = pubRef,
          target = "_blank",
          studyNam,
          style = "
            color: royalblue4;
            font-size: 28px;
            text-decoration: underline;
          "
        )
      )
    ),
    
    # ============================================================
    # PMID
    # ============================================================
    
    tags$div(
      style = "
        padding: 10px 20px;
        text-align: left;
        font-size: 30px;
        color: royalblue4;
      ",
      
      tags$span(
        tags$strong("PMID:"),
        " ",
        tags$a(
          href = pubRef,
          target = "_blank",
          "39282422",
          style = "
            color: royalblue4;
            font-size: 28px;
            text-decoration: underline;
          "
        )
      )
    )
  )
})
############################################################################################
# Create list of all the panels by starting with an empty list
initial <- list()
initial[["TitlePanel"]] <- new("TitlePanel", PanelWidth = 12L, PanelHeight = 400L, PanelId = appNum)
################################################################################
# Settings for Reduced dimension plot 1
################################################################################
initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "UMAP", XAxis = 1L, 
                                          YAxis = 2L, FacetRowByColData = "Barcode", FacetColumnByColData = "Barcode", 
                                          ColorByColumnData = "cclass", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Column data", 
                                          SizeByColumnData = "sum", FacetRowBy = "None", FacetColumnBy = "None", 
                                          ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                          ColorByFeatureName = defaultGene, ColorByFeatureSource = "---", 
                                          ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "---", 
                                          ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                          VisualChoices = c("Color", "Shape"), ContourAdd = FALSE, 
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                          CustomLabelsText = "---", FontSize = 1, 
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                          LabelCenters = FALSE, LabelCentersBy = "Barcode", LabelCentersColor = "#000000", 
                                          VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version", 
                                                                                                             "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L), 
                                          PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                          RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = TRUE, 
                                          SelectionHistory = list())
################################################################################
# Settings for Feature assay plot 1
################################################################################
initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
                                      XAxisColumnData = "cclass", XAxisFeatureName = defaultGene,
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = defaultGene, YAxisFeatureSource = "---",
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "Barcode",
                                      FacetColumnByColData = "Barcode", ColorByColumnData = "cclass",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "sample_pseudonym", SizeByColumnData = "sum", FacetRowBy = "None",
                                      FacetColumnBy = "None", ColorBy = "Column data", ColorByDefaultColor = "#000000",
                                      ColorByFeatureName = defaultGene, ColorByFeatureSource = "---",
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "{{cellone}}",
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
                                      PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                      CustomLabels = FALSE, CustomLabelsText = "{{cellone}}",
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Barcode",
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 4L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 600L,
                                      PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                      ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                      ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                      ColumnSelectionRestrict = TRUE, SelectionHistory = list())
################################################################################
# Settings for Complex heatmap 1
################################################################################
initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = heatMapGenes, 
                                        ClusterRows = TRUE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("cclass", "ctypes"), RowData = character(0), 
                                        CustomBounds = FALSE, LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE, 
                                        AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow", 
                                        ShowDimNames = "Rows", LegendPosition = "Right", LegendDirection = "Horizontal", 
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10, 
                                        ShowColumnSelection = FALSE, OrderColumnSelection = TRUE, 
                                        VersionInfo = list(iSEE = structure(list(c(2L, 6L, 0L)), class = c("package_version", "numeric_version"))), 
                                        PanelId = 1L, PanelHeight = 600L, PanelWidth = 12L, 
                                        SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                        SelectionHistory = list())
################################################################################
# Settings for Column data plot 1
################################################################################
initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nCount_RNA", 
                                    XAxisColumnData = "ctypes", FacetRowByColData = "sample_pseudonym", 
                                    FacetColumnByColData = "sample_pseudonym", ColorByColumnData = "ctypes", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "sample_pseudonym", SizeByColumnData = "nCount_RNA", 
                                    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
                                    ColorByDefaultColor = "#000000", ColorByFeatureName = "---", 
                                    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                    ColorBySampleName = "---", ColorBySampleSource = "---", 
                                    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
                                    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                    CustomLabelsText = "---", FontSize = 1, 
                                    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                    LabelCenters = FALSE, LabelCentersBy = "sample_pseudonym", LabelCentersColor = "#000000", 
                                    VersionInfo = list(iSEE = structure(list(c(2L, 6L, 0L)), class = c("package_version", 
                                                                                                       "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 12L, 
                                    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())
################################################################################
# Settings for Row data table 1
################################################################################
# initial[["RowDataTable1"]] <- new("RowDataTable", Selected = defaultGene, Search = "", SearchColumns = c("",
#                                                                                                       "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
#                                                                                                       "", "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0),
#                                   VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
#                                                                                                      "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 600L,
#                                   PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
#                                   ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
#                                   ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
#                                   ColumnSelectionRestrict = FALSE, SelectionHistory = list())
############################################################################################
# Create Custom plot for FOOTER
############################################################################################
# Create a function to use with custom plot
FOOTER_FUN <- function(se, rows, columns) {
  
  ggplot() +
    geom_rect(
      aes(xmin = 0, xmax = 100, ymin = 0, ymax = 10),
      fill = "#0B3D5B",
      color = NA
    ) +
    coord_cartesian(
      xlim = c(0, 100),
      ylim = c(0, 10),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
}
FooterPanel <- iSEE::createCustomPlot(
  FUN       = FOOTER_FUN,
  restrict  = NULL,
  className = "FooterPanel",
  fullName  = "iSEE App ID = HBCC"
)

# Subclass one panel type (e.g., ReducedDimensionPlot) and hide its Data box
setClass("Footer", contains = "FooterPanel")
setMethod(".hideInterface", "FooterPanel",
          function(x, field) {
            if (field == "DataBoxOpen" | field == "SelectionBoxOpen") return(TRUE)  # hide the whole Data parameters box
            callNextMethod()
          }
)

# Encode image into base64 encoding
socialMedia_base64 <- base64encode("socialMedia.png")

setMethod(".defineOutput", "FooterPanel", function(x) {
  
  tagList(
    
    # ============================================================
    # MAIN BLUE FOOTER
    # ============================================================
    
    tags$div(
      style = "
        background-color: #0B3D5B;
        color: white;
        padding: 30px 35px;
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        font-family: 'Open Sans', sans-serif;
      ",
      
      # --------------------------
      # LEFT SIDE
      # --------------------------
      
      tags$div(
        style = "text-align: left;",
        
        tags$a(
          href = "https://www.nimh.nih.gov/",
          target = "_blank",
          style = "
            color: white;
            text-decoration: none;
            font-size: 26px;
            font-weight: bold;
          ",
          "National Institute of Mental Health"
        ),
        
        tags$br(),
        
        tags$a(
          href = "https://www.nih.gov/",
          target = "_blank",
          style = "
            color: white;
            text-decoration: none;
            font-size: 20px;
          ",
          "at the National Institutes of Health"
        )
      ),
      
      
      # --------------------------
      # RIGHT SIDE
      # --------------------------
      
      tags$div(
        style = "
          text-align: right;
          font-size: 18px;
          line-height: 1.7;
        ",
        
        tags$div(
          style = "
            font-size: 24px;
            font-weight: bold;
          ",
          "Contact Us"
        ),
        
        tags$a(
          href = "mailto:nimhinfo@nih.gov",
          style = "
            color: white;
            text-decoration: underline;
          ",
          "nimhinfo@nih.gov"
        ),
        
        tags$br(),
        
        tags$a(
          href = "https://www.hhs.gov/",
          target = "_blank",
          style = "color: white; text-decoration: none;",
          "U.S. Department of Health and Human Services"
        ),
        
        tags$br(),
        
        tags$a(
          href = "https://www.nih.gov/",
          target = "_blank",
          style = "color: white; text-decoration: none;",
          "National Institutes of Health"
        ),
        
        tags$br(),
        
        tags$a(
          href = "https://www.nimh.nih.gov/",
          target = "_blank",
          style = "color: white; text-decoration: none;",
          "National Institute of Mental Health"
        ),
        
        tags$br(),
        
        tags$a(
          href = "https://www.usa.gov/",
          target = "_blank",
          style = "color: white; text-decoration: none;",
          "USA.gov"
        )
      )
    ),
    
    
    # ============================================================
    # FOLLOW US
    # ============================================================
    
    tags$div(
      style = "
        padding: 10px;
        text-align: left;
      ",
      
      tags$span(
        "Follow Us",
        style = "
          color: #0B3D5B;
          font-weight: bold;
          font-size: 26px;
          margin-right: 10px;
          vertical-align: middle;
        "
      ),
      
      tags$a(
        href = "https://www.nih.gov/news-events/nih-social-media",
        target = "_blank",
        
        tags$img(
          src = paste0(
            "data:image/png;base64,",
            socialMedia_base64
          ),
          height = "26px",
          alt = "NIH social media"
        )
      )
    )
  )
})

### Panel ends  
initial[["FooterPanel"]] <- new("FooterPanel", PanelWidth = 12L, PanelHeight = 800L, PanelId = appNum)
################################################################################

################################################################################
# # Overriding the default panel colors in the package
iSEEOptions$set(panel.color=c(FeatureAssayPlot="#3565AA", RowDataPlot="#3565AA", ColumnDataPlot="#3565AA", ComplexHeatmapPlot="#3565AA", RowDataTable="#3565AA", TitlePanel="#3565AA", BannerPanel1 = "white", FooterPanel = "#3565AA"))

sce_small <- registerAppOptions(sce_small, color.maxlevels = n)
################################################################################

# Add US flag along with text title
# Encode  images
flag_base64 <- base64encode("flag.png")
logo1_base64 <- base64encode("nimh-logo.png")
link_title <- "An official website of the United States government"

# Define custom CSS for iSEE app title
custom_css <- "
  .main-header .logo {
    text-align: left;
    padding-left: 75px; /* Adjust padding as needed */;
  }
"
################################################################################
# Run the iSEE wrapper function to launch the app
iSEE(
  sce_small,
  tour = tour,
  appTitle = tags$div(
    
    # US flag
    tags$img(
      src = paste0(
        "data:image/png;base64,",
        flag_base64
      ),
      height = "15px",
      alt = "USA flag",
      style = "
      vertical-align: middle;
      margin-right: 8px;
    "
    ),
    
    # Government website text
    tags$span(
      "An official website of the United States government",
      style = "
      font-size: 15px;
      font-family: 'Open Sans', sans-serif;
      color: black;
      margin-right: 15px;
      vertical-align: middle;
    "
    ),
    
    # CSS
    tags$head(
      tags$style(HTML(custom_css))
    )
  ),
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    ctypes = function(n) {col_vector[!grepl("drop", names(col_vector))]}))
)
################################################################################
# Code Ends
################################################################################
