## Tests on plot.amce
## 10/14/16 Elisha Cohen
results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` +
                  `Country of Origin` + Job + `Job Experience` + `Job Plans` +
                  `Reason for Application` + `Prior Entry`, data=immigrationconjoint,
                cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)

cat("\nTEST PLOT OBJECT")
test_that("Plot object gives error with incorrect theme",{
  expect_error(plot(results, plot.theme = fake), "object 'fake' not found")
})


test_that("Incorrect CI returns error",{
  expect_output(plot(results, ci=1.5),"Invalid confidence interval -- Defaulting to 95%")
})

test_that("Valid plot display option",{
  expect_error(plot(results, plot.display = "other"),"Error-- plot.display must be once of:  all, unconditional, interaction")
})