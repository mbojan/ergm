# vif_ergm.ergm_model() -----------------------------------------------------


data("faux.mesa.high")
frm <- faux.mesa.high ~ nodefactor("Grade") + edges + nodefactor("Sex")
em <- ergm_model(frm)
coefnames <- unlist(sapply(em$terms, "[[", "coef.names"))
fake_vcov <- matrix(0.2, 7, 7)
diag(fake_vcov) <- 2
rownames(fake_vcov) <- colnames(fake_vcov) <- coefnames



test_that("vif_ergm.ergm_model() just works for a formula with nodefactor terms", {
  expect_silent(
    r <- vif_ergm(ergm_model(frm), v = fake_vcov)
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 7)
  expect_true(nrow(r$vif_term) == 3)
})

test_that("vif_ergm.ergm_model() just works dropping term edges", {
  expect_silent(
    r <- vif_ergm(em, v = fake_vcov, drop_terms = "edges")
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 6)
  expect_true(nrow(r$vif_term) == 2)
})

test_that("vif_ergm.ergm_model() just works dropping coef edges", {
  expect_silent(
    r <- vif_ergm(em, v = fake_vcov, drop_coef = "edges")
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 6)
  expect_true(nrow(r$vif_term) == 2)
})

test_that("vif_ergm.ergm_model() just works dropping coef edges and term nodefactor('Sex')", {
  expect_silent(
    r <- vif_ergm(em, v = fake_vcov, drop_coef = "edges", drop_term='nodefactor("Sex")')
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 5)
  expect_true(nrow(r$vif_term) == 1)
})




# vif_ergm.ergm() ---------------------------------------------------------

fit <- ergm(frm)

test_that("vif_ergm.ergm() just works for a formula with nodefactor terms", {
  expect_silent(
    r <- vif_ergm(fit)
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 7)
  expect_true(nrow(r$vif_term) == 3)
})

test_that("vif_ergm.ergm_model() just works dropping term edges", {
  expect_silent(
    r <- vif_ergm(fit, drop_terms = "edges")
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 6)
  expect_true(nrow(r$vif_term) == 2)
})

test_that("vif_ergm.ergm_model() just works dropping coef edges", {
  expect_silent(
    r <- vif_ergm(fit, drop_coef = "edges")
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 6)
  expect_true(nrow(r$vif_term) == 2)
})

test_that("vif_ergm.ergm_model() just works dropping coef edges and term nodefactor('Sex')", {
  expect_silent(
    r <- vif_ergm(fit, drop_coef = "edges", drop_term='nodefactor("Sex")')
  )
  expect_s3_class(r, "vif_ergm")
  expect_true(nrow(r$vif_coef) == 5)
  expect_true(nrow(r$vif_term) == 1)
})





test_that("vif_ergm(complexModel)", {
  fit <- ergm(faux.mesa.high ~ nodefactor("Grade") + 
                nodefactor("Race") +
                nodefactor("Sex") + 
                nodematch("Grade",diff=TRUE) +
                nodematch("Race",diff=TRUE, levels=-c(1,4)) +
                nodematch("Sex",diff=FALSE) + 
                edges)
  expect_silent(r <- vif_ergm(fit, drop_coef="edges"))
})