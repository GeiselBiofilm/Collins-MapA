library(ggplot2);library(reshape2)

substrates.unique <- read.csv("Data/All_lapg_targets.csv", stringsAsFactors = FALSE, header = TRUE)
no.substrates2 <- read.csv("Data/all_non_lapg_targets.csv", stringsAsFactors = FALSE, header = TRUE)

all.pseudomonas.genomes <- read.table("Data/all_pseudomonas_genomes.txt", stringsAsFactors = FALSE, sep = '\t')$V1


protein.count <- data.frame(table(unlist(substrates.unique$organism)))
non.lapA.but.lapDG.orgs <- data.frame(table(unlist(no.substrates2$organism)))[,1]


multiples <- protein.count[protein.count$Var1[which(protein.count$Freq >1)],]
singles <- protein.count[protein.count$Var1[which(protein.count$Freq == 1)],]

multiple.rtx <- substrates.unique[which(substrates.unique$organism == multiples$Var1[1]),]
single.rtx <- substrates.unique[which(substrates.unique$organism == singles$Var1[1]),]


for ( i in 2: length(multiples$Var1)){
  multiple.rtx <- rbind(multiple.rtx, substrates.unique[which(substrates.unique$organism == multiples$Var1[i]),])
  
}

for ( i in 2: length(singles$Var1)){
  single.rtx <- rbind(single.rtx, substrates.unique[which(substrates.unique$organism == singles$Var1[i]),])
  
}

pseudomonas.multi <- multiple.rtx[ grepl( "Pseudomonas" , multiple.rtx$organism ), ]
pseudomonas.single <- single.rtx[ grepl( "Pseudomonas" , single.rtx$organism ), ]
pseudomonas.none <- no.substrates2[ grepl( "Pseudomonas" , no.substrates2$organism ), ]

pseudomonas.multi$number <- "Multiple"
pseudomonas.single$number <- "One"
pseudomonas.none$number <- "None"

pseudomonas.all <- rbind(pseudomonas.multi[,c(1,6)], pseudomonas.single[,c(1,6)], pseudomonas.none[,c(1,7)])

pseudomonas.nolapDG <- all.pseudomonas.genomes[!all.pseudomonas.genomes %in% pseudomonas.all$organism]

Shewanella.multi <- multiple.rtx[ grepl( "Shewanella" , multiple.rtx$organism ), ]
Shewanella.single <- single.rtx[ grepl( "Shewanella" , single.rtx$organism ), ]

Vibrio.multi <- multiple.rtx[ grepl( "Vibrio" , multiple.rtx$organism ), ]
Vibrio.single <- single.rtx[ grepl( "Vibrio" , single.rtx$organism ), ]

##### Output files used for Figure 11 #####

write.csv(pseudomonas.all, "Output/Pseudomonas_LapA_counts.csv", row.names = FALSE)

write.table(pseudomonas.nolapDG, "Output/Pseudomonas_no_lapDG.txt", row.names = FALSE, col.names = FALSE)


##### Plotting #####

ggplot(data = pseudomonas.multi, aes(protein.length, fill = 'multiple LapA-like protein')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,16500), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +  
  geom_density(alpha=0.25)+
  geom_density(data = pseudomonas.single, aes(protein.length, fill = 'Single LapA-like protein'), alpha=0.25) + 
  labs(x = "Protein Length", y = "Density", fill = "") +
  ggsave("Plots/Pseudomonas_substrate_sizes.png", width = 14, height = 6, units = "cm", dpi = 500)

ggplot(data = Shewanella.multi, aes(protein.length, fill = 'multiple LapA-like protein')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,16500), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.0003, len = 4), labels = function(x) format(x, scientific = TRUE)) + 
  geom_density(alpha=0.25)+
  geom_density(data = Shewanella.single, aes(protein.length, fill = 'Single LapA-like protein'), alpha=0.25) + 
  labs(x = "Protein Length", y = "Density", fill = "")+ 
  ggsave("Plots/Shewanella_substrate_sizes.png", width = 14, height = 6, units = "cm", dpi = 500)

ggplot(data = Vibrio.multi, aes(protein.length, fill = 'multiple LapA-like protein')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,16500), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +  
  geom_density(alpha=0.25)+
  geom_density(data = Vibrio.single, aes(protein.length, fill = 'Single LapA-like protein'), alpha=0.25) + 
  labs(x = "Protein Length", y = "Density", fill = "")+ 
  ggsave("Plots/Vibrio_substrate_sizes.png", width = 14, height = 6, units = "cm", dpi = 500)