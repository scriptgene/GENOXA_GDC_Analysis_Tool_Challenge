# Define selectInput options for organism database
organismDbOptions <- selectInput("organismDb", "Select Organism Database:",
                                 choices = c(
                                   "Human (org.Hs.eg.db)" = "org.Hs.eg.db", 
                                   "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db", 
                                   "Rat (org.Rn.eg.db)" = "org.Rn.eg.db",
                                   "Yeast (org.Sc.sgd.db)" = "org.Sc.sgd.db", 
                                   "Fly (org.Dm.eg.db)" = "org.Dm.eg.db", 
                                   "Arabidopsis (org.At.tair.db)" = "org.At.tair.db",
                                   "Zebrafish (org.Dr.eg.db)" = "org.Dr.eg.db", 
                                   "Bovine (org.Bt.eg.db)" = "org.Bt.eg.db", 
                                   "Worm (org.Ce.eg.db)" = "org.Ce.eg.db",
                                   "Chicken (org.Gg.eg.db)" = "org.Gg.eg.db", 
                                   "Canine (org.Cf.eg.db)" = "org.Cf.eg.db", 
                                   "Pig (org.Ss.eg.db)" = "org.Ss.eg.db",
                                   "E coli strain K12 (org.EcK12.eg.db)" = "org.EcK12.eg.db", 
                                   "Xenopus (org.Xl.eg.db)" = "org.Xl.eg.db",
                                   "Chimp (org.Pt.eg.db)" = "org.Pt.eg.db", 
                                   "Anopheles (org.Ag.eg.db)" = "org.Ag.eg.db", 
                                   "E coli strain Sakai (org.EcSakai.eg.db)" = "org.EcSakai.eg.db"
                                 ), selected = "org.Hs.eg.db"
                          )

# Define selectInput options for key type
keytypeOptions <- selectInput("keytype", "Select Key Type:",
                              choices = c(
                                "ACCNUM" = "ACCNUM",
                                "ALIAS" = "ALIAS",
                                "ENSEMBL" = "ENSEMBL",
                                "ENSEMBLPROT" = "ENSEMBLPROT",
                                "ENSEMBLTRANS" = "ENSEMBLTRANS",
                                "ENTREZID" = "ENTREZID",
                                "GENENAME" = "GENENAME",
                                "REFSEQ" = "REFSEQ",
                                "SYMBOL" = "SYMBOL",
                                "UNIGENE" = "UNIGENE",
                                "UNIPROT" = "UNIPROT"
                              ),
                              selected = "ENSEMBL"
                        )
