# Interactive example: Adjust prior and likelihood
prior <- 0.01#"base rate" of disease
likelihood <- 0.95
specificity <- 0.90

# Marginal probability of positive test
p_positive <- prior * likelihood + (1 - prior) * (1 - specificity)

posterior <- (prior * likelihood) / p_positive

cat("Posterior Probability:", round(posterior, 4))