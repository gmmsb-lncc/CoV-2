library(ggmsa)
library(RColorBrewer)
library(pals)

msa_heatmap <- function(sequences, pos_start, pos_end) {
  my_pal <- colorRampPalette(rev(brewer.pal(n = 9, name = "Reds")))
  my_custom <- data.frame(
    names = c(LETTERS[1:26], "-"), 
    color = my_pal(27), 
    stringsAsFactors = FALSE
  )
  
  pals::pal.bands(my_custom$color)
  
  res_msa <- ggmsa(sequences, 
                   pos_start, 
                   pos_end, 
                   custom_color = my_custom,
                   by_conservation = TRUE,
                   border = "white",
                   seq_name = TRUE,
                   show.legend = TRUE) + 
    geom_seqlogo(color = "Chemistry_AA") +
    geom_msaBar()
  
  return(res_msa)
}

# Input and parameters
sequences = "./canonical_consensus/VARIANTS_canonical_consensus_MSA.fasta"
pos_start = 415
pos_end   = 502

# Execute function
msa_heatmap(sequences, pos_start, pos_end)
