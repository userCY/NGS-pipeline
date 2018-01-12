# merge
NHEJ_overlap12 <- merge(NHEJ_beta_0.3_1,NHEJ_beta_0.3_2,by = 'Gene')
NHEJ_overlap13 <- merge(NHEJ_beta_0.3_1,NHEJ_beta_0.3_3,by = 'Gene')
NHEJ_overlap23 <- merge(NHEJ_beta_0.3_2,NHEJ_beta_0.3_3,by = 'Gene')
NHEJ_overlap123 <- merge(NHEJ_overlap12,NHEJ_beta_0.3_3,by = 'Gene')
NHEJ_overlap45 <- merge(NHEJ_beta_0.3_4,NHEJ_beta_0.3_5,by = 'Gene')
NHEJ_overlap46 <- merge(NHEJ_beta_0.3_4,NHEJ_beta_0.3_6,by = 'Gene')
NHEJ_overlap56 <- merge(NHEJ_beta_0.3_5,NHEJ_beta_0.3_6,by = 'Gene')
NHEJ_overlap456 <- merge(NHEJ_overlap45,NHEJ_beta_0.3_6,by = 'Gene')
