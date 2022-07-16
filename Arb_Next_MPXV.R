library(ape)
library(aplot)
library(cowplot)
library(ggforce)
library(ggtree)
library(ggfun)
library(tidytree)
library(tidyverse)
library(treeio)


#IQtree tree
tree<-read.iqtree("masked.fasta.treefile")
metadata<-read_tsv("metadata.tsv")

tree_lab<-as_tibble(tree)
tree<-groupClade(tree,425)

#hMPXV1 tree

general <-tree %>% 
    ggtree(aes(fill=group), color='gray50', size=0.9)%<+% metadata +
      theme_tree()+
      geom_hilight(node = 636, color='red', fill='gray',alpha=0.2)+
      geom_hilight(node=425, fill="steelblue",alpha =0.3, to.bottom = TRUE)+
      scale_fill_manual(values=c("gray50","cadetblue3"))+
      geom_tippoint(aes(fill=group), color="black",shape=21, size=1)+
      geom_tippoint(aes(
        subset=(country=='Mexico')),size=1.7, fill='red', color='black',shape=23)+
      geom_tippoint(aes(
        subset=(region=='Europe')), fill='gold', color='black',shape=21, size=1)+
      coord_cartesian(clip = 'off')+
      theme(legend.position = 'none')+
      geom_treescale(x=0.00030, y=6, color='black',offset = 5, fontsize = 2)+
      annotate(geom="text", x=0.00025, y=75, label="B.1",
                color="cadetblue3")

#Node 641

node_NL<-tree_subset(tree, 641, levels_back = 1)

subset<-ggtree(node_NL, color='gray50')%<+% metadata +
  theme_tree()+
  geom_tippoint(size=1, color='black',shape=21 )+
  geom_tippoint(aes(
    subset=(country=='Germany')), fill='gold', color='black',shape=21, size=1)+
  geom_tippoint(aes(
    subset=(country=='Portugal')), fill='gold', color='black',shape=21, size=1)+
  geom_tippoint(
    aes(subset=(grepl('ON911481', label, fixed=TRUE)==TRUE)),
        size=1.5, fill='red', color='black',shape=23)+
  coord_cartesian(clip = 'off')+
  geom_nodelab(aes(x=branch, label=label), size=2, hjust = -1.5, vjust=0.6) +
  geom_cladelab(node=6, label="NLE-UANL-001/Mexico/2022", align=TRUE,  
                textcolor='red',offset = 0, fontsize=3)+
  geom_cladelab(node=1, label="PT0031/Portugal/2022", align=FALSE,  
                offset = 0.0000004, fontsize=2.5)+
  geom_cladelab(node=2, label="RKI133/Germany/2022", align=FALSE,  
                offset = 0.0000004, fontsize=2.5)+
  geom_cladelab(node=3, label="RKI043/Germany/2022", align=FALSE,  
                offset = 0.0000001, fontsize=2.5)+
  geom_cladelab(node=4, label="RKI049/Germany/2022", align=FALSE,  
                offset = 0.0000001, fontsize=2.5)+
  geom_cladelab(node=5, label="RKI158/Germany/2022", align=FALSE,  
                offset = 0.0000001, fontsize=2.5)+
  geom_cladelab(node=7, label="RKI165/Germany/2022", align=FALSE,  
                offset = 0.0000001, fontsize=2.5)+
  geom_cladelab(node=8, label="RKI157/Germany/2022", align=FALSE,  
                offset = 0.0000001, fontsize=2.5)+
  hexpand(.4, direction = 1)+ 
  geom_treescale(x=0.00001, y=-0.2, color='black',
                fontsize = 2, offset = 0.25)


#Figure composition

plot_list(general, subset, tag_levels = "A",widths = c(0.4,0.7) )
ggsave('compose.tiff', width = 8, height = 6, units = "in",limitsize = FALSE, dpi = 600)


#Nextstrain

tree_nxt<-read.newick("tree.nwk")
metadata<-read_tsv("metadata.tsv")

tree_lab_nxt<-as_tibble(tree_nxt)
#tree_nxt<-groupClade(tree_nxt,425)

nxt <-tree_nxt %>% 
  ggtree(aes(fill=country), size=0.9)%<+% metadata +
  theme_tree()+
  geom_tippoint(color="black",shape=21, size=1)+
  geom_tippoint(aes(
    subset=(country=='Mexico')),size=2, fill='red', color='black',shape=23)+
  geom_tippoint(aes(
    subset=(region=='Europe')), fill='gold', color='black',shape=21, size=1)+
  coord_cartesian(clip = 'off')+
  geom_treescale(x=0.00030, y=6, color='black',offset = 5, fontsize = 2)
  
#node 444 nextstrain

node_NL_nxt<-tree_subset(tree_nxt, 444, levels_back = 1)

subset_nxt<-node_NL_nxt %>% 
  ggtree(aes(fill=country))%<+% metadata +
  theme_tree()+
  geom_tiplab(offset = 0.0000005)+
  geom_tippoint(size=1,shape=21 )+
  geom_tippoint(
    aes(subset=(grepl('ON911481', label, fixed=TRUE)==TRUE)),
    size=1.5, fill='red', color='black',shape=23)+
  coord_cartesian(clip = 'off')+
  theme(legend.position = 'none')+
  hexpand(.5, direction = 1)+ 
  geom_treescale(x=0.00001, y=-0.2, color='black',
                 fontsize = 2, offset = 0.25)

#Figure composition sup

plot_list(nxt, subset_nxt, tag_levels = "A",widths = c(0.7,0.3) )
ggsave('compose_nxt.tiff', width = 10, height = 6, units = "in",limitsize = FALSE, dpi = 600)




################### Enfoque de clado #####
set.seed(2019-08-05)
x <- rtree(30)
nn <- tidytree::offspring(x, 43, self_include=TRUE)
ggtree(x) + ggforce::facet_zoom(xy = node %in% nn)



