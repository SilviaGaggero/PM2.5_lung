# ============================================================
# Utility functions for analysis notebooks
# This file is meant to be sourced from an R Project
# ============================================================

# ---- Packages ----
# These functions assume the following packages are installed
# biomaRt, openxlsx, htmltools, base64enc, ggplot2, DESeq2


# ============================================================
# Convert Ensembl IDs to human gene names
# ============================================================
get_human_gene_names = function(gene.set, mirror = "", path_to_mart) {
  
  # Load local mart file (fallback)
  human.gene.names = read.table(path_to_mart,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE)
  
  # Load biomaRT if not already loaded
  if (!"package:biomaRt" %in% search()) {
    suppressPackageStartupMessages(library(biomaRt))
  }
  
  # Try BioMart first
  tryCatch({
    biomartCacheClear()
    
    ensembl = useEnsembl(biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      mirror = mirror)
    
    # Test connection
    testError = tryCatch(listEnsemblArchives(),
      error = function(e) e)
    
    #  If connection to ensembl works, get gene names instead of using local file
    if (!inherits(testError, "error")) {
      cat("Connecting to BioMart")
      
      human.gene.names = getBM(attributes = c("ensembl_gene_id", 
                                                       "hgnc_symbol"),
        values = gene.set,
        mart = ensembl)
    }
  }, error = function(e) {
    cat("Could not connect to BioMart, using local mart file\n\n")
  })
  
  cat("Converting Ensembl IDs to gene names\n\n")
  n.ensembl.id = length(gene.set)
  cat(paste0("Number of Ensembl IDs in the dataset: ", n.ensembl.id,"\n\n"))
  # Switch ensembl id to gene name
  gene.names = gene.set
  for (i in 1:n.ensembl.id) {
    if (gene.set[i] %in% human.gene.names$ensembl_gene_id) {
      idx = which(human.gene.names$ensembl_gene_id == gene.set[i])
      # Handle multiple matches by taking the first non-empty gene name
      selected_name = human.gene.names$hgnc_symbol[idx][1]
      if (!(selected_name %in% gene.names)) {
        gene.names[i] = selected_name
      }
    }
  }
  n.gene.names = n.ensembl.id - length(which(gene.set %in% gene.names))
  n.duplicated = n.ensembl.id - length(unique(gene.names))
  cat(paste0("Number of corresponding gene names: ", n.gene.names,"\n\n"))
  cat(paste0("Number of duplicated gene names: ", n.duplicated,"\n\n"))
  return(gene.names)
}

# ============================================================
# find GO names for human genes
# ============================================================
get_human_GO_terms = function(df, gene_col = "gene") {
  
  pkgs = c("org.Hs.eg.db", "GO.db", "AnnotationDbi", "dplyr", "purrr")
  
  missing_pkgs = pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(missing_pkgs) > 0) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
  }
  
  for(pkg in pkgs) {
    if(!paste0("package:", pkg) %in% search())
      library(pkg, character.only = TRUE)
  }
  
  # Split HGNC vs Ensembl
  is_ensembl = grepl("^ENS", df[[gene_col]])
  hgnc_genes = df[[gene_col]][!is_ensembl]
  ensembl_genes = df[[gene_col]][is_ensembl]
  
  # Map to Entrez
  entrez_hgnc = mapIds(org.Hs.eg.db,
                       keys = hgnc_genes,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  entrez_ensembl = mapIds(org.Hs.eg.db,
                          keys = ensembl_genes,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  
  entrez_ids = c(entrez_hgnc, entrez_ensembl)
  entrez_ids = entrez_ids[!is.na(entrez_ids)]
  
  if(length(entrez_ids) == 0) {
    warning("No genes could be mapped to Entrez IDs.")
    df$GO_CC = NA
    df$GO_BP = NA
    df$GO_MF = NA
    return(df)
  }
  
  # Gene → Entrez lookup
  gene_entrez_lookup = data.frame(gene = names(entrez_ids),
                                  ENTREZID = as.character(entrez_ids),
                                  stringsAsFactors = FALSE)
  
  # Fetch GO IDs
  go_list = AnnotationDbi::select(org.Hs.eg.db,
                                  keys = unique(gene_entrez_lookup$ENTREZID),
                                  columns = c("GO", "ONTOLOGY"),
                                  keytype = "ENTREZID") %>%
    filter(!is.na(GO))
  
  # Fetch GO term names
  go_term_lookup = AnnotationDbi::select(GO.db,
                                         keys = unique(as.character(go_list$GO)),
                                         columns = "TERM",
                                         keytype = "GOID")
  
  # Join GO IDs → GO names
  go_list = go_list %>%
    left_join(go_term_lookup, by = c("GO" = "GOID"))
  
  # Collapse by ontology (TERM names only)
  go_CC = go_list %>%
    filter(ONTOLOGY == "CC") %>%
    group_by(ENTREZID) %>%
    summarise(GO_CC = paste(unique(TERM), collapse = ";"), .groups = "drop")
  
  go_BP = go_list %>%
    filter(ONTOLOGY == "BP") %>%
    group_by(ENTREZID) %>%
    summarise(GO_BP = paste(unique(TERM), collapse = ";"), .groups = "drop")
  
  go_MF = go_list %>%
    filter(ONTOLOGY == "MF") %>%
    group_by(ENTREZID) %>%
    summarise(GO_MF = paste(unique(TERM), collapse = ";"), .groups = "drop")
  
  go_terms = reduce(list(go_CC, go_BP, go_MF),
                    full_join,
                    by = "ENTREZID")
  
  # Merge back to original dataframe
  df = df %>%
    left_join(gene_entrez_lookup, by = setNames("gene", gene_col)) %>%
    left_join(go_terms, by = "ENTREZID") %>%
    dplyr::select(-ENTREZID)
  
  return(df)
}

# ============================================================
# HTML download buttons for tables (TXT / XLSX)
# ============================================================
table_download = function(data,
                           base_filename = "table",
                           txt = TRUE,
                           xlsx = TRUE) {
  
  buttons = list()
  
  # ---- TXT button ----
  if (txt) {
    txt_content = paste(capture.output(write.table(data, 
                                                   sep = "\t", 
                                                   row.names = FALSE)),
      collapse = "\n")
    
    buttons = c(buttons,
      list(htmltools::tags$a(href = paste0("data:text/plain;charset=utf-8,",
            utils::URLencode(txt_content)),
          download = paste0(base_filename, ".txt"),
          class = "btn btn-default",
          style = "margin-right: 6px;",
          "Download TXT")))
  }
  
  # ---- XLSX button ----
  if (xlsx) {
    xlsx_path = tempfile(fileext = ".xlsx")
    
    wb = openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Data")
    openxlsx::writeData(wb, "Data", data)
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    
    buttons = c(buttons,
      list(htmltools::tags$a(href = xlsx_path,
          download = paste0(base_filename, ".xlsx"),
          class = "btn btn-default",
          "Download XLSX")))
    }
  
  htmltools::tags$div(style = "margin-bottom: 10px;",
    buttons )
}


# ============================================================
# Embed ggplot as static image (Notebook / HTML safe)
# ============================================================
embed_plot = function(plot, width = 6, height = 4, dpi = 200) {
  
  # Create a temporary file for the plot
  temp_file = tempfile(fileext = ".png")
  
  ggplot2::ggsave(filename = temp_file,
                  plot = plot,
                  width = width,
                  height = height,
                  dpi = dpi )
  
  # Convert to base64 and embed
  img_base64 = base64enc::base64encode(temp_file)
  
  tab = htmltools::tagList(htmltools::tags$img(src = paste0("data:image/png;base64,", 
                                                            img_base64),
                                               style = "max-width: 100%;"))
  
  print(tab)
  unlink(temp_file)
}


# ============================================================
# Extract results from DESeq2 object as data frame
# ============================================================
export_DE_results = function(df, contrast, fc_limit=1, padj_limit=0.05) {
  
  df = DESeq2::results(df, contrast= contrast)
  df = df[order(df$padj),]
  
  # Convert to dataframe and remove NA padj
  df=as.data.frame(df)
  df=na.omit(df)
  
  # Add gene name column
  df$gene = rownames(df)
  rownames(df) = NULL
  
  # add labels for differentially expressed genes
  df$diff.expr. = "NO"
  df$diff.expr.[df$log2FoldChange > fc_limit & df$padj < padj_limit] = "UP"
  df$diff.expr.[df$log2FoldChange < -fc_limit & df$padj < padj_limit] = "DOWN"
  
  return(df)
}

# ============================================================
# P-value histogram plot
# ============================================================
pval_hist = function(df, main = "") {
  
  # Load ggplot2 if not already loaded
  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(library(ggplot2))
  }
  
  df = as.data.frame(df)
  
  # Remove rows with NA adjusted p-values
  if (any(is.na(df$padj))) {
    df = df[!is.na(df$padj), , drop = FALSE]
  }
  
  p = ggplot(df, ggplot2::aes(x = pvalue)) +
    geom_histogram(breaks = seq(0, 1, 0.02),
                   fill = "black",       
                   color = "black") +
    ggtitle(main) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, NA),
                       labels = scales::comma) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, NA))
  
  return(p)
}


# ============================================================
# Volcano plot
# ============================================================
volcano.plot = function(df, lfc = 1, padj = 0.05, main = "", 
                        imp.genes = c(), top_n = 20,
                        point_stroke_color = c("DOWN" = "#00344A", 
                                               "UP" = "#B53737", 
                                               "NO" = "#B2B2B2"),
                        point_fill_colors = c("DOWN" = "#628392", 
                                              "UP" = "#D38383", 
                                              "NO" = "#D1D1D1"),
                        bg_down_color = "#F4F9FF",
                        bg_up_color = "#FFF4F4") {
  
  # Load ggplot2 if not already loaded
  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(library(ggplot2))
  }
  
  # Compute y-axis value
  df$neglog10padj = -log10(df$padj)
  
  # Determine important genes to label
  if (!is.null(imp.genes)) {
    imp.genes = intersect(imp.genes, df$gene[df$diff.expr. != "NO"])
  } else {
    top_up = df %>%
      filter(diff.expr. == "UP") %>%
      arrange(desc(log2FoldChange)) %>%
      slice_head(n = top_n) %>%
      pull(gene)
    
    top_down = df %>%
      filter(diff.expr. == "DOWN") %>%
      arrange(log2FoldChange) %>%
      slice_head(n = top_n) %>%
      pull(gene)
    
    low_padj = df %>%
      filter(diff.expr. %in% c("UP", "DOWN")) %>%
      mutate(score = abs(log2FoldChange) / padj) %>%
      group_by(diff.expr.) %>%
      slice_max(score, n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      pull(gene)
    
    imp.genes = unique(c(top_up, top_down, low_padj))
    imp.genes = imp.genes[!grepl("ENSG", imp.genes)]
  }
  
  # Count UP / DOWN genes 
  n_up = sum(df$diff.expr. == "UP", na.rm = TRUE)
  n_down = sum(df$diff.expr. == "DOWN", na.rm = TRUE)
  
  # Axis lims for annotation placement 
  xmax = max(abs(df$log2FoldChange), na.rm = TRUE)
  ymax = max(df$neglog10padj, na.rm = TRUE)
  
  # Volcano plot
  p = ggplot(data=df, aes(text=gene, 
                             x=log2FoldChange,
                             y=neglog10padj)) +
    # background regions 
    annotate("rect", xmin = -Inf, xmax = -lfc, ymin = -log10(padj), ymax = Inf,
             fill = bg_down_color) +
    annotate("rect", xmin =  lfc, xmax = Inf, ymin = -log10(padj), ymax = Inf,
             fill = bg_up_color) +
    # Base scatter
    geom_point(aes(color = diff.expr., fill = diff.expr.),
               size = 1.5, alpha = 0.8, shape = 21, stroke = 0.3) +
    scale_color_manual(values = point_stroke_color) +
    scale_fill_manual(values = point_fill_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))+
    geom_vline(xintercept=c(-lfc, lfc), color="#B2B2B2",
               linetype="dashed", linewidth=0.3) +
    geom_hline(yintercept=-log10(padj), color="#B2B2B2",
               linetype="dashed", linewidth=0.3) +
    labs(title = main,
         subtitle = paste(n_up, "UP genes &", n_down, "DOWN genes"),
         x = expression(log[2]("fold change")),
         y = expression(-log[10]("padj"))) + 
    theme_classic() +
    
    # reset scales so labeled genes get their own
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    
    # Highlight and label important genes
      geom_point(data = subset(df, gene %in% imp.genes),
               aes(color = diff.expr., fill = diff.expr.),
               size = 2, shape = 21, stroke = 0.4) +
    ggrepel::geom_label_repel(data = subset(df, gene %in% imp.genes),
                              aes(label = gene),
                              size = 3.4,
                              max.overlaps = 30,
                              min.segment.length = 0.1,
                              force = 1,
                              box.padding = 0.5,
                              label.size = NA,
                              label.padding = 0.1,
                              fill = alpha(c("white"), 0.5)) +     
    scale_color_manual(values = point_stroke_color, guide = "none") +
    scale_fill_manual(values = point_stroke_color, guide = "none")

  return(p)
  
}


# ============================================================
# MA plot
# ============================================================
MA_plot = function(df, lfc = 1, padj = 0.05, imp.genes = NULL, 
                   top_n = 10, main = "",
                   point_stroke_color = c("DOWN" = "#00344A", 
                                          "UP"   = "#B53737", 
                                          "NO"   = "#B2B2B2"),
                   point_fill_colors  = c("DOWN" = "#628392", 
                                          "UP"   = "#D38383", 
                                          "NO"   = "#D1D1D1")) {
  
  # Convert to data frame and remove NAs
  df = as.data.frame(df)
  df = na.omit(df)
  
  # Add diff expression labels
  df$diff.expr. = "NO"
  df$diff.expr.[df$log2FoldChange > lfc & df$padj < padj] = "UP"
  df$diff.expr.[df$log2FoldChange < -lfc & df$padj < padj] = "DOWN"
  
  # Add gene name column if missing
  if (is.null(df$gene)) df$gene = rownames(df)
  
  # Count UP / DOWN genes
  n_up = sum(df$diff.expr. == "UP", na.rm = TRUE)
  n_down = sum(df$diff.expr. == "DOWN", na.rm = TRUE)
  
  # Determine important genes to label
  if (is.null(imp.genes)) {
    top_up = df %>%
      filter(diff.expr. == "UP") %>%
      arrange(desc(log2FoldChange)) %>%
      slice_head(n = top_n) %>%
      pull(gene)
    
    top_down = df %>%
      filter(diff.expr. == "DOWN") %>%
      arrange(log2FoldChange) %>%
      slice_head(n = top_n) %>%
      pull(gene)
    
    low_padj = df %>%
      filter(diff.expr. %in% c("UP", "DOWN")) %>%
      mutate(score = abs(log2FoldChange) / padj) %>%
      group_by(diff.expr.) %>%
      slice_max(score, n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      pull(gene)
    
    imp.genes = unique(c(top_up, top_down, low_padj))
    imp.genes = imp.genes[!grepl("ENSG", imp.genes)]
  } else {
    imp.genes = intersect(imp.genes, df$gene[df$diff.expr. != "NO"])
  }
  
  # Load ggplot2 if not already loaded
  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(library(ggplot2))
  }
  
  # Base MA plot
  p = ggplot(df, aes(x = log2(baseMean), y = log2FoldChange, text = gene)) +
    geom_point(aes(color = diff.expr.), alpha = 0.6) +
    scale_color_manual(values = point_stroke_color) +
    geom_hline(yintercept = c(-lfc, lfc), color = "#B2B2B2", 
               linetype = "dashed", 
               linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))+
    labs(title = main,
         subtitle = paste(n_up, "UP genes &", n_down, "DOWN genes"),
         x = expression(log[2]("mean expression")),
         y = expression(log[2]("fold change"))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    
    # Highlight and label important genes
    geom_point(data = subset(df, gene %in% imp.genes),
               aes(color = diff.expr., fill = diff.expr.),
               size = 2.5, shape = 21, stroke = 0.4) +
    ggrepel::geom_label_repel(data = subset(df, gene %in% imp.genes),
                              aes(label = gene),
                              size = 3.4,
                              max.overlaps = 30,
                              min.segment.length = 0.1,
                              force = 1,
                              box.padding = 0.5,
                              label.size = NA,
                              label.padding = 0.1,
                              fill = alpha("white", 0.5)) +
    scale_color_manual(values = point_stroke_color, guide = "none") +
    scale_fill_manual(values = point_stroke_color, guide = "none")
  
  return(p)
}

# ============================================================
# Significance brackets
# ============================================================

add_sig_brackets = function(sig_df, group_levels, y_base, step = NULL,
                             bracket_width = 0.08, text_size = 2.5,
                             sig_threshold = 0.05) {
  
  if (is.null(sig_df) || nrow(sig_df) == 0) return(list())
  
  sig_df = sig_df[sig_df$pvalue < sig_threshold, , drop = FALSE]
  if (nrow(sig_df) == 0) return(list())
  
  # Sort by span so wider brackets stack on top
  sig_df$idx1 = match(sig_df$group1, group_levels)
  sig_df$idx2 = match(sig_df$group2, group_levels)
  sig_df$span = abs(sig_df$idx2 - sig_df$idx1)
  sig_df = sig_df[order(sig_df$span), ]
  
  if (is.null(step)) step = max(y_base * 0.08, 0.5)
  
  layers = list()
  for (i in seq_len(nrow(sig_df))) {
    row = sig_df[i, ]
    y_pos = y_base + i * step
    
    # Format p-value text
    if (row$pvalue < 0.001) {
      pval_text = "p < 0.001"
    } else {
      pval_text = sprintf("p = %.3f", row$pvalue)
    }
    
    layers = c(layers, list(
      # Left tick
      annotate("segment", x = row$idx1 - bracket_width,
               xend = row$idx1 - bracket_width,
               y = y_pos - step * 0.3, yend = y_pos, linewidth = 0.4),
      # Right tick
      annotate("segment", x = row$idx2 + bracket_width,
               xend = row$idx2 + bracket_width,
               y = y_pos - step * 0.3, yend = y_pos, linewidth = 0.4),
      # Horizontal line
      annotate("segment", x = row$idx1 - bracket_width,
               xend = row$idx2 + bracket_width,
               y = y_pos, yend = y_pos, linewidth = 0.4),
      # P-value text
      annotate("text", x = (row$idx1 + row$idx2) / 2,
               y = y_pos + step * 0.15, label = pval_text,
               size = text_size, vjust = 0, hjust = 0.5)))
  }
  return(layers)
}
