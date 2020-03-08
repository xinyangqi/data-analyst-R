# Stats 506, Fall 2019
# This script contains solutions for problem set 4
#
# This script creates four beautiful tables for results of problem set 4,
# which are output by Stata.
#
# Author : Xinyang Qi
# Updated : Nov 18 , 2019


library(tidyverse)
library(readxl)

# q1

q1=read_excel("ps4_1.xls") %>%
  rename(measures=...1)
knitr::kable(q1,caption = "Condition relative effect for curvature measures.(95% CI)")

# q2

decode_division = function(x) {
  # Decodes numeric codes for DIVISION
  # Args: 
  #   x: a single numeric code for DIVISION
  #
  # Output : the labels
  
  # Throw an error if x isn't numeric
  if(!is.numeric(x)) stop('decode_states expects numeric input indexed from 1!')
  switch(x,
         "New England","Middle Atlantic","East North Central","West North Central",
         "South Atlantic","East South Central","West South Central",
         "Mountain North","Mountain South","Pacific")
}
decode_all_division = function(x){
  # Vectorizes decode_division above
  #
  # Args: 
  #  x: a vector of integer-valued reportable_domains
  #
  # Returns: A vector of labels or a "list" if some are unmatched.
  sapply(x,decode_division)
}

q2=read.csv("ps4_q2.csv")
q2=mutate(q2,division=decode_all_division(division))
knitr::kable(q2,caption = "Internet disparity between urban and rural.(95% CI)")

# q3(c)

q3_c=read_excel("ps4_3c.xls") %>%
  rename(measures=...1)
knitr::kable(q3_c,caption = "Odd ratio and margin effect of the first model .(95% CI)")

# q3(d)

q3_d=read_excel("ps4_3d.xls") %>%
  rename(measures=...1)
knitr::kable(q3_d,caption = "Odd ratio and margin effect of the second model .(95% CI)")
