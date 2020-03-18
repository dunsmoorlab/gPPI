require(lme4)
require(ez)

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
