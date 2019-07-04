## This script can be used to generate a legend for the maternal SNP density

# 
# SNPs_density$alpha[SNPs_density$MAT_Freq < 0.25] <- 0.3 + (0.15-SNPs_density$MAT_Freq[SNPs_density$MAT_Freq < 0.25])/0.2
# SNPs_density$alpha[SNPs_density$MAT_Freq > 0.5] <- 0.4 + (SNPs_density$MAT_Freq[SNPs_density$MAT_Freq > 0.5]-0.4)/0.6
# SNPs_density$alpha[SNPs_density$alpha > 1] <- 1
# SNPs_density$alpha[SNPs_density$alpha < 0] <- 0
# 
# SNPs_density$color[SNPs_density$MAT_Freq < 0.25] <- rgb(255,128,0, SNPs_density$alpha[SNPs_density$MAT_Freq < 0.25]*255, maxColorValue = 255)
# SNPs_density$color[SNPs_density$MAT_Freq > 0.5] <- rgb(127,56,236,  SNPs_density$alpha[SNPs_density$MAT_Freq > 0.5]*255, maxColorValue = 255)

library(ggplot2)
output_folder <- "~/hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Results/Embryo_karyograms_1000000/"

pdf(paste(output_folder, "Legend_PAT_MAT.pdf", sep = ""), width = 5, height = 1)
par(mar=c(2, 4, 2, 2))
plot.new()  

xlim <- c(0,1000)
ylim <- c(0,2)

plot.window(xlim=xlim, ylim=ylim)

title("Maternal SNP density")

colors <- data.frame(x = seq(1,1000, 1), color = "white", alpha = 1, stringsAsFactors = F)
colors$alpha[colors$x < 250] <-  0.3 + (150-colors$x[colors$x < 250])/250

colors$alpha[colors$x > 500] <-  0.3 + (colors$x[colors$x > 500]-500)/500

colors$alpha[colors$alpha > 1] <- 1
colors$alpha[colors$alpha < 0] <- 0

colors$color[colors$x < 250] <-  rgb(255,128,0, colors$alpha[colors$x < 250]*255, maxColorValue = 255)
colors$color[colors$x > 500] <-  rgb(127,56,236, colors$alpha[colors$x > 500]*255, maxColorValue = 255)


rect(xleft = (colors$x-1), xright = colors$x, ybottom = 0, ytop = 2, col = colors$color, border = NA)
rect(xleft = 0, xright = 1000, ybottom = 0, ytop = 2, col = NA, border = "black")
axis(1, line = 0, at =  seq(0,1000, 100), col = NA, col.ticks =  "black", labels = F)

axis(1, line = 0, labels = seq(0,1, 0.2), at =  seq(0,1000, 200), col = NA, col.ticks = "black")

dev.off()


## Legend for copy number states

dummy_data <- data.frame(CN_state = c("0", "1", "2", "3", "4", "5+"), colors = c("gray", "#E41A1C", "chartreuse3", "#2A6AFF", "#1034A6", "#111E6C"), score = 1, stringsAsFactors = F)
dummy_data
ggplot(dummy_data, aes(x = colors, fill = CN_state)) + geom_bar() + scale_fill_manual(values = dummy_data$colors) +
  theme_bw(base_size = 8) + 
  labs(fill = "Copy\nNumber") +
  theme(strip.text.x = element_text(size = 7, face = "bold"), 
        panel.border = element_rect(colour = "black"),
        legend.key.size = unit(0.6, "line"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title =  element_text(size = 7),
        panel.grid.major.x = element_blank(),
        text=element_text(family="Arial"))
ggsave(filename = paste(output_folder, "Legend_CN_State.pdf", sep = ""), width = 70, height = 55, units = "mm", dpi = 300, device = cairo_pdf)



