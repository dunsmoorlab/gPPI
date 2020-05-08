require(lme4)
require(ez)
setwd("C:/Users/ACH/Dropbox (LewPeaLab)/Gus/gPPI")

lev <- read.csv('../Desktop/level_df.csv')
lmod <- ezANOVA(lev,dv=.(rsa),wid=.(subject),within=.(encode_phase,level,roi),between=.(group))

c_lev <- subset(lev, group %in% 'control')
p_lev <- subset(lev, group %in% 'ptsd')

c_lev_mod <- ezANOVA(c_lev,dv=.(rsa),wid=.(subject),within=.(encode_phase,level,roi)) 
p_lev_mod <- ezANOVA(p_lev,dv=.(rsa),wid=.(subject),within=.(encode_phase,level,roi))


ws <- read.csv('../Desktop/memory_df.csv')
wmod <- ezANOVA(ws,dv=.(rsa),wid=.(subject),within=.(encode_phase,memory_phase,roi),between=.(group))

c_ws <- subset(ws, group %in% 'control')
p_ws <- subset(ws, group %in% 'ptsd')

c_ws_mod <- ezANOVA(c_ws,dv=.(rsa),wid=.(subject),within=.(encode_phase,memory_phase,roi)) 
p_ws_mod <- ezANOVA(p_ws,dv=.(rsa),wid=.(subject),within=.(encode_phase,memory_phase,roi)) 


us <- read.csv('us_rsa.csv')
usa <- ezANOVA(us,dv=.(rsa),wid=.(subject),within=.(roi,US),between=.(group))

mem <- read.csv('mem_rsa.csv')
usa <- ezANOVA(mem,dv=.(rsa),wid=.(subject),within=.(roi,high_confidence_accuracy,encode_phase),between=.(group))

ers <- read.csv('ers_cs_comp.csv')
erA <- ezANOVA(ers,dv=.(rsa),wid=.(subject),within=.(roi,encode_phase),between=.(group))

fdf <- read.csv('foi_phase.csv')
fA <- ezANOVA(fdf,dv=.(rsa),wid=.(subject),within=.(roi,encode_phase),between=.(group))

