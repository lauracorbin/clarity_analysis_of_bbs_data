# modified bubble heatmap functions

MyformatData <- function(data,
                       rowName = "rowName",
                       colName = "colName",
                       plotGroup = "plotGroup",
                       colorValue = "colorValue",
                       sizeValue = "sizeValue",
                       rowOrder = NULL,
                       colOrder = NULL,
                       listOrder = NULL,
                       nightingale = TRUE) {
  data <- as.data.frame(data)
  print(paste0("nightingale set of: ", nightingale))
  # Split data into plot groups
  data <- split(data, data[, plotGroup])
  print(summary(data))
  # Cast matrices
  reformatCast <- function(dat, val) {
    if (nrow(dat) == 1) {
      m <- matrix(data = as.numeric(dat[val]),
                  nrow = 1,
                  ncol = 1,
                  dimnames = list(dat[, rowName], dat[, colName]))
    } else if (nrow(unique(dat[colName])) == 1) {
      m <- matrix(data = as.numeric(as.matrix(dat[val])),
                  nrow = nrow(dat),
                  ncol = 1,
                  dimnames = list(dat[, rowName], dat[, colName][1]))
    } else {
      m <- as.matrix(
        reshape::cast(dat,
                      formula = as.formula(paste(rowName, " ~ ", colName)),
                      value = val))
    }
  }
  colorList <- lapply(data, reformatCast, val = colorValue)
  print(summary(colorList))
  sizeList <- lapply(data, reformatCast, val = sizeValue)
  print(summary(sizeList))
  # Sort list of plot matrices#
    # listOrder <- c("Absolute Measures", "Ratio to Total Lipids", "Aggregate",
    #                "Fatty Acids", "Mean Diameter", "Other Lipids",
    #                "Apolipoproteins", "Amino Acids", "BCAAs", "Glycolysis",
    #                "Ketones", "Fluid Balance", "Inflammation")
      listOrder <- c("Absolute Measures", "Aggregate",
                     "Fatty Acids", "Mean Diameter", "Other Lipids",
                     "Apolipoproteins", "Amino Acids", "BCAAs", "Glycolysis",
                     "Ketones", "Fluid Balance", "Inflammation")
  if (!(is.null(listOrder))) {
    colorList <- colorList[listOrder]
    sizeList <- sizeList[listOrder]
  }
  # Sort rows and columns of plot matrices
  matSort <- function(x, listIn, order, RC) {
    if (x %in% names(order)) {
      if (RC == "C") {
        m <- as.matrix(rbind(listIn[[x]][, colOrder[[x]]]))
        row.names(m) <- row.names(listIn[[x]])
        return(m)
      } else if (RC == "R") {
        #m <- as.matrix(cbind(listIn[[x]][rowOrder[[x]], ]))
        if(length(listIn[[x]]) < length(rowOrder[[x]])) {
          m <- as.matrix(cbind(listIn[[x]]))
          } else {
            m <- as.matrix(cbind(listIn[[x]][rowOrder[[x]], ]))
            }
        colnames(m) <- colnames(listIn[[x]])
        return(m)
      }
    } else {
      return(listIn[[x]])
    }
  }

  colOrder <- list(
    `Absolute Measures` = c("CM & XL VLDL", "VL VLDL", "L VLDL", "M VLDL",
                            "S VLDL", "VS VLDL", "IDL", "L LDL", "M LDL",
                            "S LDL", "VL HDL", "L HDL", "M HDL", "S HDL"),
    # `Ratio to Total Lipids` = c("CM & XL VLDL", "VL VLDL", "L VLDL", "M VLDL",
    #                             "S VLDL", "VS VLDL", "IDL", "L LDL", "M LDL",
    #                             "S LDL", "VL HDL", "L HDL", "M HDL", "S HDL"),
    #`Aggregate` = c("Total", "VLDL", "LDL", "HDL", "Clinical LDL", "Non-HDL", "Remnant"),
    `Aggregate` = c("VLDL", "LDL", "HDL", "Clinical LDL", "Non-HDL", "Remnant"),
    #`Fatty Acids` = c("Total FA", "o Unsaturation", "Sat FA", "PUFA", "MUFA", "n-6", "22:6 DHA",
     #                 "n-3","18:2 LA", "PUFA:MUFA", "n-6:n-3")
    `Fatty Acids` = c("o Unsaturation", "Sat FA", "PUFA", "MUFA", "n-6", "22:6 DHA",
                      "n-3","18:2 LA")
  )
  rowOrder <- list(
    `Absolute Measures` = c("Particle Concentration", "Total Lipids",
                            "Total Cholesterol", "Free Cholesterol",
                            "Cholesterol Esters", "Triglycerides",
                            "Phospholipids"),
    `Aggregate` = c("Particle Concentration", "Total Lipids", "Total Cholesterol",
                    "Free Cholesterol", "Cholesterol Esters", "Triglycerides",
                    "Phospholipids"),
    #`Apolipoproteins` = c("ApoA1", "ApoB", "ApoB:ApoA1"),
    `Apolipoproteins` = c("ApoA1", "ApoB"),
    `Mean Diameter` = c("VLDL", "HDL", "LDL"),
    #`Ratio to Total Lipids` = c("Total Cholesterol", "Free Cholesterol",
     #                           "Cholesterol Esters", "Triglycerides",
      #                          "Phospholipids"),
    `Fatty Acids` = c("Concentration", "Ratio to Total*"),
    `Amino Acids` = c("Ala", "Gln", "Gly", "His", "Phe", "Tyr"),
    #`BCAAs` = c("Ile", "Leu", "Val", "Total"),
    `BCAAs` = c("Ile", "Leu", "Val"),
    `Fluid Balance` = c("Albumin", "Creatinine"),
    `Glycolysis` = c("Glucose", "Lactate", "Citrate",
                     "Pyruvate", "Glycerol"),
    `Ketones` = c("3OH-butyrate", "Acetate",
                  "Acetoacetate", "Acetone"),
    #`Other Lipids` = c("PG", "Trig:PG", "Phosphatidylcholine",
     #                  "Cholines", "Sphingomyelins")
    `Other Lipids` = c("PG", "Phosphatidylcholine",
                       "Cholines", "Sphingomyelins")
  )
if (!is.null(colOrder)) {
  colorList <- lapply(setNames(names(colorList), names(colorList)), matSort,
                      listIn = colorList, order = colOrder, RC = "C")
  sizeList <- lapply(setNames(names(sizeList), names(sizeList)), matSort,
                     listIn = sizeList, order = colOrder, RC = "C")
}
if (!is.null(rowOrder)) {
  colorList <- lapply(setNames(names(colorList), names(colorList)), matSort,
                      listIn = colorList, order = rowOrder, RC = "R")
  sizeList <- lapply(setNames(names(sizeList), names(sizeList)), matSort,
                     listIn = sizeList, order = rowOrder, RC = "R")
}
return(list(colorList = colorList, sizeList = sizeList))
}


###############################################################################
###############################################################################
###############################################################################


MymultiPlotInput <- function(colorList,
                           sizeList = colorList,
                           nightingale = TRUE,
                           treeNames = names(colorList),
                           context = gpar(cex = 0.8),
                           colorSeq = c("#053061", "#2166AC", "#4393C3",
                                        "#92C5DE", "#D1E5F0", "#FDDBC7",
                                        "#F4A582", "#D6604D", "#B2182B",
                                        "#67001F"),
                           unitBase = unit(0.5, "cm"),
                           diameter = 0.8,
                           colorLim = c(NA, NA),
                           sizeLim = c(NA, NA),
                           showTopLabels = TRUE,
                           showLeftLabels = TRUE,
                           leftLabelsTitles = FALSE,
                           showRowBracket = FALSE,
                           rowTitles = FALSE,
                           showColBracket = FALSE,
                           colTitles = FALSE,
                           plotTitles = FALSE,
                           xTitles = FALSE,
                           yTitles = FALSE,
                           showColorLegend = TRUE,
                           showBubbleLegend = TRUE,
                           colorBreaks = NULL,
                           sizeBreaks = NULL,
                           legendHeight = nrow(colorList[[1]]),
                           legendTitles = c(expression("-log"[10]*"p"), "r")) {
  #### Prepare Data Matrices ####
  # If min or max limits are preset for color or size then#
  #edit any values outside limits#
  if (!all(is.na(c(colorLim, sizeLim)))) {
    applyLim <- function(x, min, max) {
      if (!is.na(min)) x[x < min] <- min
      if (!is.na(max)) x[x > max] <- max
      return(x)
    }
    colorList <- lapply(colorList, applyLim, colorLim[1], colorLim[2])
    sizeList <- lapply(sizeList, applyLim, sizeLim[1], sizeLim[2])
    for(i in 1:length(colorList)){
      storage.mode(colorList[[i]]) <= "numeric"
      storage.mode(sizeList[[i]]) <= "numeric"
    }}
  #If limits not preset then set range using min/max across all plots#
  if(any(is.na(c(colorLim, sizeLim)))) {
    if(is.na(c(colorLim[1]))){colorLim[1] <- min(do.call(c, colorList)[is.finite(as.numeric(do.call(c, colorList)))])}
    if(is.na(c(colorLim[2]))){colorLim[2] <- max(do.call(c, colorList)[is.finite(as.numeric(do.call(c, colorList)))])}
    if(is.na(c(sizeLim[1]))){sizeLim[1] <- min(do.call(c, sizeList)[is.finite(as.numeric(do.call(c, sizeList)))])}
    if(is.na(c(sizeLim[2]))){sizeLim[2] <- max(do.call(c, sizeList)[is.finite(as.numeric(do.call(c, sizeList)))])}}
  #### Set Legend Breakpoints####
  #If bubble legend labels are not supplied and#
  #bubble legend is required then set size breakpoints#
  if (is.null(sizeBreaks) & showBubbleLegend != FALSE) {
    sizeBreaks <- pretty(c(sizeLim[1], sizeLim[2]), legendHeight)
    if (sizeBreaks[1] > sizeLim[1]) sizeBreaks <-
        c(sizeBreaks[1] - (sizeBreaks[2] - sizeBreaks[1]), sizeBreaks)
    if (tail(sizeBreaks, 1) < sizeLim[2]) sizeBreaks <-
        c(sizeBreaks, tail(sizeBreaks, 1) + (sizeBreaks[2] - sizeBreaks[1]))
    sizeLim <- c(sizeBreaks[1], tail(sizeBreaks, 1))
  }
  #If color legend labels are not supplied and#
  #color legend is required then set color breakpoints#
  #(Attempt to match legend sizes)#
  if (is.null(colorBreaks) & showColorLegend != FALSE) {
    if (!(is.null(sizeBreaks))) legendHeight <- length(sizeBreaks)
    colorBreaks <- pretty(c(colorLim[1], colorLim[2]), legendHeight)
    if (colorBreaks[1] > colorLim[1]) colorBreaks <-
      c(round(colorBreaks[1] - (colorBreaks[2] - colorBreaks[1]),digits=1), round(colorBreaks,digits=1))
    if (tail(colorBreaks, 1) < colorLim[2]) colorBreaks <-
      c(round(colorBreaks,digits=1), round(tail(sizeBreaks, 1) + (colorBreaks[2] - colorBreaks[1]),digits=1))
    colorLim <- c(colorBreaks[1], tail(colorBreaks, 1))
  }
  #### Set color palette####
  colorSeq <- grDevices::colorRampPalette(colorSeq)(200)
  #### Set Variable Vectors For Each Plot####
  if (nightingale == TRUE) {
    plotTitles <- c(rep.int(FALSE, length(colorList)))
    xTitles <- rep.int(FALSE, length(colorList))
    yTitles <- rep.int(FALSE, length(colorList))
    showTopLabels <- c(TRUE, TRUE, TRUE, TRUE, rep.int(FALSE, length(colorList) - 4))
    showLeftLabels <- c(TRUE, TRUE, FALSE, rep.int(TRUE, length(colorList) - 3))
    leftLabelsTitles <- names(colorList)
    #leftLabelsTitles[3] <- FALSE
    showRowBracket <- rep.int(FALSE, length(colorList))
    rowTitles <- rep.int(FALSE, length(colorList))
    showColBracket <- rep.int(FALSE, length(colorList))
    colTitles <- rep.int(FALSE, length(colorList))
    showBubbleLegend <- TRUE
    showColorLegend <- TRUE
    showTopLabels[names(colorList) == "Aggregate"] <- TRUE
  } else {
    varList <- c("plotTitles", "xTitles", "yTitles", "showTopLabels",
                 "showLeftLabels", "leftLabelsTitles", "showRowBracket",
                 "rowTitles", "showColBracket", "colTitles",
                 "showBubbleLegend", "showColorLegend")
    for (x in varList) {
      if (length(get(x)) != (0 || 1 || length(colorList))) {
        stop(paste("length", x, "must equal 0, 1 or length(colorList)"))
      }
    }
  }
  if (length(showColorLegend) == 1 & showColorLegend == TRUE) {
    sCLOverall <- TRUE
    showColorLegend <- FALSE
  }
  if (length(showBubbleLegend) == 1 & showBubbleLegend == TRUE) {
    sBLOverall <- TRUE
    showBubbleLegend <- FALSE
  }
  vectorList <- list(showTopLabels = showTopLabels,
                     showLeftLabels = showLeftLabels,
                     showRowBracket = showRowBracket,
                     showColBracket = showColBracket,
                     rowTitle = rowTitles,
                     colTitle = colTitles,
                     xTitle = xTitles,
                     yTitle = yTitles,
                     leftLabelsTitle = leftLabelsTitles,
                     plotTitle = plotTitles,
                     showBubbleLegend = showBubbleLegend,
                     showColorLegend = showColorLegend)
  #### Combine variables to list for output ####

  inputA <- c(list(f = list, colorMat = colorList,
                   sizeMat = sizeList, treeName = treeNames), vectorList)
  inputA <- do.call(Map, inputA)
  inputB <- c(list(colorSeq = colorSeq, unitBase = unitBase,
                   colorLim = colorLim, sizeLim = sizeLim,
                   colorBreaks = colorBreaks, sizeBreaks = sizeBreaks,
                   legendTitles = legendTitles, legendHeight = legendHeight,
                   context = context,
                   diameter = diameter))
  if (sCLOverall == TRUE || sBLOverall == TRUE) {
    legend <- c(list(colorSeq = colorSeq, unitBase = unitBase,
                     diameter = diameter, context = context,
                     colorBreaks = round(colorBreaks,digits=1), sizeBreaks = round(sizeBreaks,digits=1),
                     legendTitles = legendTitles, legendHeight = legendHeight,
                     showColorLegend = sCLOverall,
                     showBubbleLegend = sBLOverall))
  }
  return(list(trees = lapply(inputA, function(a, b) {
     c(a, b)
   }, b = inputB), legends = legend))
}


###############################################################################
###############################################################################
###############################################################################

MymetabFigurePlot <- function(treeList) {

  #Get unit size
  unitBase <- treeList[[1]]$children$GridSq$y[1] - treeList[[1]]$children$GridSq$y[2]

  # List layout groupings#
  pG <- list(
    pH = c("Absolute Measures", "Aggregate", "Apolipoproteins"),
    pW = c("Absolute Measures", "Other Lipids", "Glycolysis"),
    p1 = "Absolute Measures",
    p2 = c("Aggregate", "Apolipoproteins", "Mean Diameter"),
    p2B = c("Apolipoproteins", "Mean Diameter"),
    #p3 = "Ratio to Total Lipids",
    p4 = "Fatty Acids",
    p5A = c("Other Lipids", "Glycolysis"),
    p5B = c("Ketones", "Fluid Balance", "Inflammation"),
    p5C = c("Amino Acids", "BCAAs"),
    p6 = "Legends")

  #Split Fatty Acids
  if("PUFA:MUFA" %in% treeList$`Fatty Acids`$children$TopLabels$label){
    len <- length(treeList$`Fatty Acids`$children$TopLabels$x)
    treeList$`Fatty Acids`$children$TopLabels$x[(len-1):len] <-
      treeList$`Fatty Acids`$children$TopLabels$x[(len-1):len] + 0.5*unitBase
    len <- length(treeList$`Fatty Acids`$children$GridSq$x)
    treeList$`Fatty Acids`$children$GridSq$x[(len-1):len] <-
      treeList$`Fatty Acids`$children$GridSq$x[(len-1):len] + 0.5*unitBase
    len <- length(treeList$`Fatty Acids`$children$Bubbles$x)
    treeList$`Fatty Acids`$children$Bubbles$x[(len-1):len] <-
      treeList$`Fatty Acids`$children$Bubbles$x[(len-1):len] + 0.5*unitBase
    treeList$`Fatty Acids` <- getLayoutSizes(treeList$`Fatty Acids`, unitBase)}

  ##### Align widths#####
  treeList[pG$p1]$`Absolute Measures`$childrenvp$parent$layout$widths[4] <- unit(1, "null")
  treeList[pG$pW] <- matchLayoutWidths(treeList[pG$pW], 1, 2)
  #treeList[c(pG$p1, pG$p3)] <- matchLayoutWidths(treeList[c(pG$p1, pG$p3)], 1, 5)
  treeList[c(pG$p1)] <- matchLayoutWidths(treeList[c(pG$p1)], 1, 5)
  treeList[c(pG$p4, pG$p2B)] <- matchLayoutWidths(treeList[c(pG$p4, pG$p2B)], 4, 4)
  treeList[pG$p2B] <- matchLayoutWidths(treeList[pG$p2B], 1, 2)
  treeList[pG$p5B] <- matchLayoutWidths(treeList[pG$p5B], 1, 2)
  treeList[pG$p5C] <- matchLayoutWidths(treeList[pG$p5C], 1, 2)
  treeList[pG$p2[1]][[1]]$childrenvp$parent$layout$widths[4] <- unit(0, "cm")
  treeList[pG$p4][[1]]$childrenvp$parent$layout$widths[1] <- unit(1, "null")

  ##### Align heights#####
  treeList[pG$pH] <- matchLayoutHeights(treeList[pG$pH], 1, 3)
  treeList[pG$pH[1:2]] <- matchLayoutHeights(treeList[pG$pH[1:2]], 1, 5)
  #treeList[pG$p4][[1]]$childrenvp$parent$layout$heights[1] <- unit(1, "null")
  #treeList[pG$p3][[1]]$childrenvp$parent$layout$heights[1] <- unit(1, "null")

  # Combine apolipoproteins & mean diameter (panel 2B)
  p2B <- arrangePanels("Panel2B", subTrees = treeList[pG$p2B],
                       nrow = 2, ncol = 1, just = c("right", "top"),
                       spacer = unit(0, "null"))
  p2B <- gTree(name="Panel2B", children=do.call(gList, p2B$grobList), childrenvp=p2B$vp)


  # Combine nine miscellaneous plots into one tree (panel 5)#
  treeList[pG$p5C[2]][[1]]$children$LeftLabelsTitle$gp <- gpar(fontface="italic")
  p5A <- arrangePanels("Misc1", subTrees = treeList[pG$p5A],
                       nrow = 2, ncol = 1, just = c(0, 0),
                       spacer = unit(1, "null"))
  p5B <- arrangePanels("Misc2", subTrees = treeList[pG$p5B],
                       nrow = 3, ncol = 1, just = c(0, 0),
                       spacer = unit(1, "null"))
  p5C <- arrangePanels("Misc3", subTrees = treeList[pG$p5C],
                       nrow = 2, ncol = 1, just = c(0, 0),
                       spacer = unit(1, "null"))
  p5 <- combinePanels("Misc", panels = list(p5A, p5B, p5C),
                      nrow = 1, ncol = 3, just = c(0, 0),
                      spacer = unit(1, "null"))

  # Arrange final panels
  metabPanel <-
    arrangePanels("metabTree",
                  subTrees = c(treeList[pG$p1],
                               treeList[pG$p2[1]],
                               LP_Other = list(p2B),
                               #treeList[pG$p3],
                               #list(nullGrob()), 
                               Misc = list(p5),
                               treeList[pG$p6],
                               list(nullGrob()),
                               treeList[pG$p4],
                               list(nullGrob()),
                               list(nullGrob())
                               ),
                  nrow = 3, ncol = 3, spacer = unitBase)

  # Fine adjustments to alignment
  #metabPanel$grobList[pG$p4][[1]]$childrenvp$parent$layout.pos.col <- as.integer(c(1, 2))
  #metabPanel$grobList[pG$p4][[1]]$childrenvp$parent$valid.pos.col <- as.integer(c(1, 2))
  metabPanel$grobList[pG$p6][[1]]$childrenvp$parent$layout.pos.col <- as.integer(c(4, 6))
  metabPanel$grobList[pG$p6][[1]]$childrenvp$parent$valid.pos.col <- as.integer(c(4, 6))
  metabPanel$vp$layout$heights[3] <- 0.5*unitBase
  metabPanel$vp$layout$widths[4] <-
    sum(metabPanel$grobList[pG$p2[1]][[1]]$childrenvp$parent$layout$widths)
  metabPanel$vp$layout$widths[6] <-
    sum(metabPanel$grobList$LP_Other$childrenvp$layout$widths)
  metabPanel$vp$layout$widths[5] <- max(
    sum(metabPanel$grobList[pG$p4][[1]]$childrenvp$parent$layout$widths) -
      (sum(metabPanel$grobList[pG$p2[1]][[1]]$childrenvp$parent$layout$widths) +
         sum(metabPanel$grobList$LP_Other$childrenvp$layout$widths)),
    unit(1, "char"))
  metabTree <- gTree(name="metabTree", children=do.call(gList,
                                                        metabPanel$grobList), childrenvp=metabPanel$vp)
  # ----- selectively remove LeftLabels above the Aggregate panel -----
  remove_selected_leftlabels <- function(grob, target_sections = c("Aggregate")) {
    # If the current grob name matches one of the target sections
    if (!is.null(grob$name) && any(target_sections %in% grob$name)) {
      if (!is.null(grob$children)) {
        leftlabs <- names(grob$children)[grepl("LeftLabels", names(grob$children))]
        if (length(leftlabs) > 0) {
          for (nm in leftlabs) grob$children[[nm]] <- grid::nullGrob()
        }
      }
    }
    
    # Recurse into children
    if (!is.null(grob$children)) {
      for (nm in names(grob$children)) {
        grob$children[[nm]] <- remove_selected_leftlabels(grob$children[[nm]], target_sections)
      }
    }
    if (!is.null(grob$grobList)) {
      for (i in seq_along(grob$grobList)) {
        grob$grobList[[i]] <- remove_selected_leftlabels(grob$grobList[[i]], target_sections)
      }
    }
    grob
  }
  
  remove_selected_toplabels <- function(grob, target_sections = c("MeanDiameter")) {
    # If the current grob name matches one of the target sections
    if (!is.null(grob$name) && any(target_sections %in% grob$name)) {
      if (!is.null(grob$children)) {
        toplabs <- names(grob$children)[grepl("TopLabels", names(grob$children))]
        if (length(toplabs) > 0) {
          for (nm in toplabs) grob$children[[nm]] <- grid::nullGrob()
        }
      }
    }
    
    # Recurse into children
    if (!is.null(grob$children)) {
      for (nm in names(grob$children)) {
        grob$children[[nm]] <- remove_selected_toplabels(grob$children[[nm]], target_sections)
      }
    }
    if (!is.null(grob$grobList)) {
      for (i in seq_along(grob$grobList)) {
        grob$grobList[[i]] <- remove_selected_toplabels(grob$grobList[[i]], target_sections)
      }
    }
    grob
  }
  
  # Remove TopLabels only from panels that sit above Aggregate
  #metabTree <- remove_selected_leftlabels(metabTree,
  #                                        target_sections = c("Aggregate"))
  metabTree <- remove_selected_toplabels(metabTree,
                                         target_sections = c("MeanDiameter"))
  return(metabTree)
}


