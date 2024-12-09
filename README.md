# From eggs to adulthood: sustained effects of early developmental temperature and corticosterone exposure on physiology and body size in an Australian lizard

[![DOI](https://zenodo.org/badge/489925392.svg)](https://doi.org/10.5281/zenodo.14343654)

# Citation
This paper is now published and can be cited as follows:

Crino, O. L., Wild, K.H., Friesen, C., Leibold, D., Laven, N., Peardon, A.Y., Recio, P., Salin., K., and Noble, D.W. 2024. From eggs to adulthood: sustained effects of early developmental temperature and corticosterone exposure on physiology and body size in an Australian lizard. Journal of Experimental Biology, 227, jeb249234. doi: 10. 1242/jeb.249234

# Abstract
Developing animals are increasingly exposed to elevated temperatures as global temperatures rise due to climate change. Vertebrates can be affected by elevated temperatures during development directly, and indirectly through maternal effects (e.g., exposure to prenatal glucocorticoid hormones). Past studies have examined how elevated temperatures and glucocorticoid exposure during development independently affect vertebrates. However, exposure to elevated temperatures and prenatal corticosterone could have interactive effects on developing animals that affect physiology and life history traits across life. We tested interactions between incubation temperature and prenatal corticosterone exposure in the delicate skink (Lampropholis delicata). We treated eggs with high or low dose corticosterone treatments and incubated eggs at 23°C (cool) or 28°C (warm). We measured the effects of these treatments on development time, body size, and survival from hatching to adulthood and on adult hormone levels and mitochondrial respiration. We found no evidence for interactive effects of incubation temperature and prenatal corticosterone exposure on phenotype. However, incubation temperature and corticosterone treatment each independently decreased body size at hatching and these effects were sustained into the juvenile period and adulthood. Lizards exposed to low doses of corticosterone during development had elevated levels of baseline corticosterone as adults. Additionally, lizards incubated at cool temperatures had higher levels of baseline corticosterone and more efficient mitochondria as adults compared to lizards incubated at warm temperatures. Our results show that developmental conditions can have sustained effects on morphological and physiological traits in oviparous lizards but suggest that incubation temperature and prenatal corticosterone do not have interactive effects.

# Scripts
Lampro_Cort_Temp_Analysis_final.R contains the code needed for all analysis and figures. Running this script will provide the final models used in the manuscript that are found in the 'models' folder. 

The .rmd file 'Tables_SI.rmd' provide the code needed to construct all supplementary tables. These are knited to a word document.

# Data files
All raw data used for scripts can be found in folder "data".In this folder you will find raw data used for analysis,script, and figures. In this folder you will find 3 .csv files of raw data that are used in the Lampro_Cort_Temp_Analysis_final.R script. The 'Yolk_hormone.csv' provides the hormone treatment effects on yolk corticosterone levels that was tested on a subset of individuals. The 'Cort_hatch_juv.csv' includes the morphometric data from hathling and juvenile individuals. The 'Cort_adult.csv' includes the adult data for all individuals including final morphometrics, hormone levels, and mitochondrial bioenergetics (basal, leak, OXPHOS, and RCR) for all individuals. The two .csv files 'Cort_hatch_juv.csv' and 'Cort_adult.csv' are combined in final analysis 'Lampro_Cort_Temp_Analysis_final.R' to provide how individuals grew over time. Finally, the .csv file 'interaction_table.csv' is the interaction table (S1) in the supplementary. 

# Final files
All saved models can be found in folder "models" and are saved as .rds files. There you will find each saved model that describes models of interest that was saved as an output from the analysis script 'Lampro_Cort_Temp_Analysis_final.R'. 

# Figure Files
All figures (and supplementary figures) used in manuscript can be found in "figures" folder.

# Manuscript
Final manuscript and supplementary information can be found in the folder 'MS'
 