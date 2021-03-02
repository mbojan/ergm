f <- function(e) {
  o <- capture.output(e)
  cat(paste0("## ", o), sep="\n")
}



data("faux.mesa.high")
frm <- faux.mesa.high ~ nodefactor("Grade") + nodefactor("Race") +
  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
  nodematch("Race",diff=TRUE, levels=-c(1,4)) + nodematch("Sex",diff=FALSE) + edges
frm.full <- faux.mesa.high ~ edges + nodefactor("Grade") + nodefactor("Race") +
  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
  nodematch("Race",diff=TRUE, keep=-c(1,4)) + nodematch("Sex",diff=FALSE) +
  gwdegree(1.0,fixed=TRUE) + gwesp(1.0,fixed=TRUE) + gwdsp(1.0,fixed=TRUE) 
fit <- ergm(frm)
em <- ergm_model(fit$formula, nw = fit$network)
str(em$terms)
em.full <- ergm_model(frm.full, nw = faux.mesa.high)
str(em.full$terms)
term_names <- structure(
  lapply(em$terms, function(x) list(coef.names = x$coef.names, )),
  names = vapply(em$terms, "[[", character(1), "name")
)
to_drop <- lapply()


R <- cov2cor(vcov(fit))

vif_ergm_new(fit)
vif_ergm_new(fit, drop_terms = "edges")

