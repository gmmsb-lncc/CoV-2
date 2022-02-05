library(ggmsa)
library(RColorBrewer)
library(pals)


# heatmap range
sequences = "./canonical_consensus/VARIANTS_canonical_consensus_MSA.fasta"
pos_start = 415
pos_end   = 502

msa_heatmap <- function(sequences, pos_start, pos_end){
  
  my_pal <- colorRampPalette(rev(brewer.pal(n = 9, name = "Reds")))
  my_cutstom <- data.frame(names = c(LETTERS[1:26],"-"), 
                           color = my_pal(27), 
                           stringsAsFactors = FALSE)
  head(my_cutstom)
  
  pals::pal.bands(my_cutstom$color)
  
  res_msa <- ggmsa(sequences, 
                  pos_start, 
                  pos_end, 
        #char_width = 0,
                  font = NULL,
                  custom_color = my_cutstom,
                  by_conservation = TRUE,
                  border = "white",
                  seq_name = TRUE,
                  show.legend = TRUE) + 
                  geom_seqlogo(color = "Chemistry_AA") +
                  geom_msaBar()
  
  return(res_msa)  
}

msa_heatmap(sequences, pos_start, pos_end)

