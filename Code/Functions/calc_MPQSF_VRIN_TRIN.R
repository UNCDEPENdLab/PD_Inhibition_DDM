calc_MPQSF_VRIN_TRIN <- function(df){
# browser()
    df <- df %>% mutate_at(paste0("MPQ_",seq(1, 155,1)), funs(recode(., `1`=0, `2`=1)))

VRIN <- df %>%
  mutate(VRINbf1 = abs(MPQ_5 - MPQ_112),
         VRINbf2 = abs(MPQ_20 - MPQ_139), 
         VRINbf3 = abs(MPQ_44 - MPQ_140), 
         VRINbf4 = abs(MPQ_2 - MPQ_110), 
         VRINbf5 = abs(MPQ_88 - MPQ_135), 
         VRINbf6 = abs(MPQ_51 - MPQ_45), 
         VRINbf7 = abs(MPQ_6 - MPQ_90),
         VRINbf8 = abs(MPQ_76 - MPQ_146),
         VRINbf9 = abs(MPQ_38 - MPQ_109),
         VRINbf10 = abs(MPQ_30 - MPQ_91),
         VRINbf11 = abs(MPQ_82 - MPQ_142),
         VRINbf12 = abs(MPQ_33 - MPQ_152),
         VRINbf13 = abs(MPQ_50 - MPQ_62),
         VRINbf14 = abs(MPQ_28 - MPQ_136),
         VRINbf15 = abs(MPQ_52 - MPQ_111),
         VRINbf16 = abs(MPQ_65 - MPQ_148),
         VRINbf17 = abs(MPQ_71 - MPQ_130),
         VRINbf18 = abs(MPQ_85 - MPQ_144),
         VRINbf19 = abs(MPQ_132 - MPQ_53),
         VRINbf20 = abs(MPQ_15 - MPQ_117),
         VRINbf21 = abs(MPQ_10 - MPQ_29)
  ) %>% select(subject, starts_with("VRIN")) %>% mutate("VRIN" = rowSums(select(., -subject))) %>% 
  select(subject, VRIN) %>% mutate("Z_VRIN" = (VRIN-5.3511)/2.11416,
"T_VRIN" = Z_VRIN*10 + 50)

TRIN <- df %>% mutate(
  TRINbft1 = ifelse(MPQ_28 == 1 & MPQ_112 == 1, 1,0),
  TRINbft2 = ifelse(MPQ_52 == 1 & MPQ_64 == 1,1,0),
  TRINbft3 = ifelse(MPQ_70 == 1 & MPQ_118 == 1,1,0),
  TRINbft4 = ifelse(MPQ_99 == 1 & MPQ_146 == 1,1,0),
  TRINbft5 = ifelse(MPQ_5 == 1 & MPQ_136 == 1,1,0),
  TRINbft6 = ifelse(MPQ_20 == 1 & MPQ_79 == 1,1,0),
  TRINbft7 = ifelse(MPQ_15 == 1 & MPQ_98 == 1,1,0),
  TRINbft8 = ifelse(MPQ_110 == 1 & MPQ_145 == 1,1,0),
  
  TRINbff1 = ifelse(MPQ_51 == 0 & MPQ_63 == 0,1,0),
  TRINbff2 = ifelse(MPQ_40 == 0 & MPQ_148 == 0,1,0),
  TRINbff3 = ifelse(MPQ_69 == 0 & MPQ_46 == 0,1,0),
  TRINbff4 = ifelse(MPQ_80 == 0 & MPQ_140 == 0,1,0),
  TRINbff5 = ifelse(MPQ_5 == 0 & MPQ_136 == 0,1,0),
  TRINbff6 = ifelse(MPQ_20 == 0 & MPQ_79 == 0,1,0),
  TRINbff7 = ifelse(MPQ_15 == 0 & MPQ_98 == 0,1,0),
  TRINbff8 = ifelse(MPQ_110 == 0 & MPQ_145 == 0,1,0),
  TRINbff9 = ifelse(MPQ_132 == 0 & MPQ_60 == 0,1,0),
  TRINbff10 = ifelse(MPQ_86 == 0 & MPQ_9 == 0,1,0),
  TRINbff11 = ifelse(MPQ_30 == 0 & MPQ_73 == 0,1,0),
  TRINbff12 = ifelse(MPQ_32 == 0 & MPQ_104 == 0,1,0),
) %>% select(subject, starts_with("TRIN"))

TRINbft <- TRIN %>% select(subject, starts_with("TRINbft")) %>% mutate("TRINbft" = rowSums(select(., -subject))) %>% 
  select(subject, TRINbft) 
TRINbff <- TRIN %>% select(subject, starts_with("TRINbff")) %>% mutate("TRINbff" = rowSums(select(., -subject))) %>% 
  select(subject, TRINbff) 

TRIN <- left_join(TRINbft, TRINbff, by = "subject") %>% mutate("TRIN" = 12+TRINbft - TRINbff,
                                                               "Z_TRIN" = (TRIN - 11.8296)/1.49445,
                                                               "T_TRIN" = (Z_TRIN*10)+50)

valids <- left_join(VRIN, select(TRIN, subject, TRIN, Z_TRIN, T_TRIN), by = "subject")

return(valids)

}


