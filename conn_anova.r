require(ez)


setwd("C:/Users/ACH/Documents/gPPI")

c <- read.csv('pfc_roi_conn.csv')
str(c)
c$subject <- as.factor(c$subject)

hc_control_ext_acq <- c[with(c,group=='healthy' & seed %in% c('hc_head','hc_body','hc_tail') & cope=='ext_acq'),]
am_control_ext_acq <- c[with(c,group=='healthy' & seed %in% c('amyg_cem','amyg_bla') & cope=='ext_acq'),]
hc_ptsd_ext_acq <- c[with(c,group=='ptsd' & seed %in% c('hc_head','hc_body','hc_tail') & cope=='ext_acq'),]
am_ptsd_ext_acq <- c[with(c,group=='ptsd' & seed %in% c('amyg_cem','amyg_bla') & cope=='ext_acq'),]

hc_control_acq <- c[with(c,group=='healthy' & seed %in% c('hc_head','hc_body','hc_tail') & cope=='acq_csp_csm'),]
am_control_acq <- c[with(c,group=='healthy' & seed %in% c('amyg_cem','amyg_bla') & cope=='acq_csp_csm'),]
hc_ptsd_acq <- c[with(c,group=='ptsd' & seed %in% c('hc_head','hc_body','hc_tail') & cope=='acq_csp_csm'),]
am_ptsd_acq <- c[with(c,group=='ptsd' & seed %in% c('amyg_cem','amyg_bla') & cope=='acq_csp_csm'),]


hc_control_ext <- c[with(c,group=='healthy' & seed %in% c('hc_head','hc_body','hc_tail') & cope=='ext_csp_csm'),]
am_control_ext <- c[with(c,group=='healthy' & seed %in% c('amyg_cem','amyg_bla') & cope=='ext_csp_csm'),]
hc_ptsd_ext <- c[with(c,group=='ptsd' & seed %in% c('hc_head','hc_body','hc_tail') & cope=='ext_csp_csm'),]
am_ptsd_ext <- c[with(c,group=='ptsd' & seed %in% c('amyg_cem','amyg_bla') & cope=='ext_csp_csm'),]


ezANOVA(data=am_ptsd_ext,dv=.(conn),wid=.(subject),within=.(seed,target),type=2)


