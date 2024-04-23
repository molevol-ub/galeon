#### Script to estimate the probability of finding two or more gene copies in a specific distance

args = commandArgs(trailingOnly = TRUE)

#### COMMAND LINE INPUT ARGUMENTS ####
number_genes = as.numeric(args[1])
genome_size = as.numeric(args[2])
g = args[3]
outdir = args[4]

#### OUTPUT FILE ####
outfile = "g_estimation.table.txt"
outfile = paste(outdir, outfile, sep="/")

## Example data:
## Edit the variables and run all code
# number_genes <- 3357 # Specify here the number of genes from your gene family (i.e.: 411)
# genome_size <- 3592 # (Mb) Specify here the genome size in Mb (i.e: 1365.69 Mb here)
# g <- 100 # (kb) Specify the maximum distance between two copies of a given family to consider that they are clustered in a unit of kb (i.e.: 100 kb)
## Several g values can be provided as such: "100,200,300"

#### FUNCTIONS ####
# 1 #
read_g_argum <- function(i_g_argument) {
  
  # Function to read the input "g argument" and detect whether it is just one value or a list of values
  
  if (grepl(",", i_g_argument) == FALSE) {
    out_g_arg <- as.numeric(i_g_argument)
    return(out_g_arg)
    
  } else if ((grepl(",", i_g_argument) == TRUE)) {
    temp <- as.list(strsplit(i_g_argument, ","))[[1]]
    out_g_arg <- as.numeric(unlist(temp))
    
    return(out_g_arg)
  }
}

# 2 #
estimate_g <- function(i_number_genes, i_genome_size, i_gval) {
  
  # Function to estimate several statistics, based on input Gene number, Genome size and g value.
  
  # Code
  mb_per_gene <- i_genome_size/i_number_genes
  mb_per_gene <- round(mb_per_gene, 2)
  # cat ("We would expect to find one gene each", mb_per_gene, "Mb \n")
  genes_per_mb <- i_number_genes/i_genome_size
  # cat ("- ", genes_per_mb, "genes are expected to be found each Mb \n")
  genes_per_kb <- (i_number_genes/1000)/i_genome_size
  # cat ("- ", genes_per_kb, "genes are expected to be found each kb \n")
  genes_per_g <- (i_number_genes/(1000/i_gval))/i_genome_size
  # cat ("- ", genes_per_g, "genes are expected to be found each", i_gval,"kb, specified as g \n")

  p <- 1 - ppois(1, lambda = genes_per_g)
  # cat ("The probability of finding by chance two (or more) genes in a", i_gval, "kb stretch is p =", p, " (Poisson distribution, λ =", genes_per_g, ") ### \n " )
  out_info <- c(i_gval, mb_per_gene, genes_per_mb, genes_per_kb, genes_per_g, p, genes_per_g)
  return(out_info)
}

## EXECUTION ##

# Read g value argument and parse it
g <- read_g_argum(g)

# Estimate parameters

if (length(g) == 1) { # If only one g value is provided
  
  info_list = estimate_g(number_genes, genome_size, g)
  info_df <- data.frame(matrix(ncol=7, nrow=0))
  info_df <- rbind(info_df, info_list)
  colnames(info_df) <- c("g value", "Exp. 1 gene each X Mb", "Exp. genes / Mb", "Exp. Genes / kb", "Exp. Genes/g value",  "P(X>=2) / g value", "Poisson's λ")
  
  # 
  
  info = paste(c("g value: ", g, " --- P(X>=2) = ", info_list[6]), sep = "", collapse = "")
  print(info)
  
} else if (length(g) > 1) { # If several g value are provided
  
  info_list = list()
  
  for (i in 1:length(g)) {
    gvalue = g[i]

    temp = estimate_g(number_genes, genome_size, gvalue)
    info_list[[i]] <- temp
    
    info = paste(c(i, ") g value: ", gvalue, " --- P(X>=2) = ", temp[6]), sep = "", collapse = "")
    print(info)
  
  }
  
  info_df <- as.data.frame(do.call(rbind, info_list))
  colnames(info_df) <- c("g value", "Exp. 1 gene each X Mb", "Exp. genes / Mb", "Exp. Genes / kb", "Exp. Genes/g value",  "P(X>=2) / g value ", "Poisson's λ")

} else {
  print("Unknown error, check the inputs")
}

if (file.exists(outfile) == FALSE) {
  write.table(info_df, outfile, sep = "\t", row.names = FALSE) 
} else if (file.exists(outfile) == TRUE) {
  write.table(info_df, outfile, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
}

temp = paste("Script output table: ", outfile, sep="")
print(temp)

info_df