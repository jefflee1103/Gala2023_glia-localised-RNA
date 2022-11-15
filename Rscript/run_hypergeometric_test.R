run_hypergeometric_test <- function(overlap, group_a, group_b, total){
  phyper(
    overlap - 1, 
    group_a, 
    total - group_a,
    group_b,
    lower.tail = FALSE
  )
}
