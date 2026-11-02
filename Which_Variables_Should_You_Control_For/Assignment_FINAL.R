
library(ggplot2)
library(ggdag)
library(DiagrammeR)
library(iWeigReg)
library(Matching)
library(dagitty)
library(drgee)

# Set seed, same as source article

set.seed(28359)

# Data generating process

N <- 500 

simulated_data_list <- list()

# 1000 simulated data sets with sample size 500

for (i in 1:1000) {
  
  # Define parameters
  
  n <- 10000

  beta_YT <- 0.5 # Treatment effect

  # Observed confounder is also drawn from the standard normal distribution N
  
  U <- rnorm(n, mean = 0, sd = 1)

  # All continuous covariates with no incoming arrows were 
  # drawn from the standard normal distribution N

  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  
  # X3 is a function of U and an error term
  
  X3_raw <- U + rnorm(n, mean = 0, sd = 1)
  
  # X3 was standardized afterward to have zero mean and unit variance.
  X3 <- as.numeric(scale(X3_raw))
  
  X4 <- rnorm(n, mean = 0, sd = 1)
  
  # X5 is a function of X4 and an error term
  
  X5_raw <- X4 + rnorm(n, mean = 0, sd = 1)
  
  X5 <- as.numeric(scale(X5_raw))
  
  X6 <- rnorm(n, mean = 0, sd = 1)
  X7 <- rnorm(n, mean = 0, sd = 1)
  
  # X8 is a function of X2, X6, and an error term
  
  X8_raw <- (X2) + (X6) + rnorm(n, mean = 0, sd = 1)
  
  X8 <- as.numeric(scale(X8_raw))
  
  X9 <- rnorm(n, mean = 0, sd = 1)

  # Calculate the dependent variables

  # Generate treatment 
  
  prob_T <- 1 / (1 + exp(-(X1 + X2 + X3 + X4)))
  
  Tr <- rbinom(n, 1, prob_T)
  
  Y <- (U) + (X5) + (X6) + (X7) +
    (beta_YT * Tr) + rnorm(n, mean = 0, sd = 1)

  simulated_data_list[[i]] <- data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9,
                                         Tr, Y = Y)
  
  

}
# We subset the simulated data sets to the first 500 rows.
simulated_data_list <- lapply(simulated_data_list, function(x){return(x[1:N, ])})

library(ggplot2)
library(ggforce)
library(patchwork)

# We create the DAG structure for the data generating process as described in the source article. 

nodes <- data.frame(
  name = c("W", "Y", "U", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"),
  x    = c(2.0, 6.0, 4.8, 0.2, 1.2, 3.5, 2.8, 5.2, 5.5, 7.2, 3.5, 7.5),
  y    = c(3.0, 3.0, 4.2, 3.0, 1.2, 4.2, 5.5, 5.5, 1.2, 4.5, 1.2, 1.8),
  type = c("treatment", "outcome", "unobserved", rep("covariate", 9)),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c("X1","X2","X2","X4","X4","X3","U", "U", "X5","X6","X6","X7","W"),
  to   = c("W", "W", "X8","W", "X5","W", "X3","Y", "Y", "Y", "X8","Y", "Y"),
  stringsAsFactors = FALSE
)
edges$dashed <- ifelse(edges$from == "U", TRUE, FALSE)

edges <- merge(edges, nodes[, c("name","x","y")], by.x = "from", by.y = "name")
names(edges)[names(edges) %in% c("x","y")] <- c("x1","y1")
edges <- merge(edges, nodes[, c("name","x","y")], by.x = "to", by.y = "name")
names(edges)[names(edges) %in% c("x","y")] <- c("x2","y2")

shorten_arrow <- function(x1, y1, x2, y2, r1 = 0.41, r2 = 0.41) {
  dx <- x2 - x1; dy <- y2 - y1
  d  <- sqrt(dx^2 + dy^2)
  ux <- dx / d;   uy <- dy / d
  data.frame(
    xs = x1 + ux * r1, ys = y1 + uy * r1,
    xe = x2 - ux * r2, ye = y2 - uy * r2
  )
}

arrow_coords <- with(edges, shorten_arrow(x1, y1, x2, y2))
edges <- cbind(edges, arrow_coords)

target_sets <- list(
  list(
    title   = expression(bold(X) * "*: Full set {" * X[1] * " ... " * X[9] * "}"),
    members = c("X1","X2","X3","X4","X5","X6","X7","X8","X9"),
    color   = "#0F6E56"
  ),
  list(
    title   = expression(bold(X)[Y] * ": Outcome-oriented {" * X[3] * ", " * X[5] * ", " * X[6] * ", " * X[7] * "}"),
    members = c("X3", "X5", "X6", "X7"),
    color   = "#0F6E56"
  ),
  list(
    title   = expression(bold(X)[T] * ": Treatment-oriented {" * X[1] * ", " * X[2] * ", " * X[3] * ", " * X[4] * "}"),
    members = c("X1", "X2", "X3", "X4"),
    color   = "#0F6E56"
  ),
  list(
    title   = expression(bold(X)[min] * ": Minimal set {" * X[3] * ", " * X[5] * "}"),
    members = c("X3", "X5"),
    color   = "#0F6E56"
  )
)

make_dag_panel <- function(set_info) {
  members <- set_info$members
  col     <- set_info$color

  cov_nodes <- nodes[nodes$type == "covariate" & nodes$name %in% members, ]
  cov_nodes$fill   <- col
  cov_nodes$border <- col
  cov_nodes$tcol   <- "white"
  cov_nodes$r      <- 0.30

  visible <- c("W", "Y", "U", members)

  edges_vis <- edges[edges$from %in% visible & edges$to %in% visible, ]
  edges_solid  <- edges_vis[!edges_vis$dashed, ]
  edges_dashed <- edges_vis[edges_vis$dashed, ]
  edges_solid$lw <- ifelse(edges_solid$from == "W", 0.9, 0.5)

  special <- nodes[nodes$type %in% c("treatment", "outcome", "unobserved"), ]

  p <- ggplot() +
    coord_fixed(xlim = c(-0.6, 8.4), ylim = c(-0.2, 6.8)) +
    {if (nrow(edges_solid) > 0)
      geom_segment(
        data = edges_solid,
        aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = lw),
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        color = "#A0A098", alpha = 0.5
      )
    } +
    scale_linewidth_identity() +
    {if (nrow(edges_dashed) > 0)
      geom_segment(
        data = edges_dashed,
        aes(x = xs, y = ys, xend = xe, yend = ye),
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        linewidth = 0.5, color = "#A0A098", alpha = 0.35,
        linetype = "dashed"
      )
    } +
    {if (nrow(cov_nodes) > 0)
      geom_circle(
        data = cov_nodes,
        aes(x0 = x, y0 = y, r = r, fill = fill, color = border)
      )
    } +
    scale_fill_identity() +
    scale_color_identity() +
    scale_alpha_identity() +
    {if (nrow(cov_nodes) > 0)
      geom_text(
        data = cov_nodes,
        aes(x = x, y = y, label = name, colour = tcol),
        size = 3.2, fontface = "italic"
      )
    } +
    geom_circle(
      data = special[special$name == "W", ],
      aes(x0 = x, y0 = y, r = 0.38),
      fill = "#7F77DD", color = "#534AB7", linewidth = 0.8
    ) +
    geom_text(
      data = special[special$name == "W", ],
      aes(x = x, y = y), label = "W",
      color = "white", size = 5, fontface = "bold.italic"
    ) +
    geom_circle(
      data = special[special$name == "Y", ],
      aes(x0 = x, y0 = y, r = 0.38),
      fill = "#D85A30", color = "#993C1D", linewidth = 0.8
    ) +
    geom_text(
      data = special[special$name == "Y", ],
      aes(x = x, y = y), label = "Y",
      color = "white", size = 5, fontface = "bold.italic"
    ) +
    geom_circle(
      data = special[special$name == "U", ],
      aes(x0 = x, y0 = y, r = 0.28),
      fill = NA, color = "#999999", linewidth = 0.5, linetype = "dashed"
    ) +
    geom_text(
      data = special[special$name == "U", ],
      aes(x = x, y = y), label = "U",
      color = "#999999", size = 3.2, fontface = "italic"
    ) +
    ggtitle(set_info$title) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 5)),
      plot.margin = margin(8, 8, 8, 8)
    )

  return(p)
}

panels <- lapply(target_sets, make_dag_panel)

combined <- (panels[[1]] | panels[[2]]) /
            (panels[[3]] | panels[[4]]) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(
        size = 16, face = "bold", hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  )

print(combined)

#ggsave(
#  filename = "vector_dag.pdf", 
#  plot = combined, 
#  width = 10, 
#  height = 8,
#  bg = "transparent")

################## LINEAR REGRESSION #########################################################################

# Create function to extract ACE from linear regression model

get_ace <- function(data, formula) {
  model <- lm(formula, data = data)
  return(coef(model)["Tr"]) # We extract the coefficient for the treatment variable "Tr" which represents the estimated average causal effect (ACE) if the model is
  # correctly specified. 
}

# We then apply this function to each target set with our simulated data
# set to get the ACE estimates 

# Minimal set

ace_min <- sapply(simulated_data_list, get_ace, formula = Y ~ Tr + X3 + X5)

# Full set

ace_full <- sapply(simulated_data_list, get_ace,
                   formula = Y ~ Tr + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9)

# Treatment oriented
ace_treat <- sapply(simulated_data_list, get_ace, formula = Y ~ Tr + X1 + X2 + X3 + X4)

# Outcome oriented
ace_out <- sapply(simulated_data_list, get_ace, formula = Y ~ Tr + X3 + X5 + X6 + X7)

# For later analysis

# Minimal set + X1

ace_min2 <- sapply(simulated_data_list, get_ace, formula = Y ~ Tr + X3 + X5 + X1)

# Outcome oriented + X1
ace_out2 <- sapply(simulated_data_list, get_ace, formula = Y ~ Tr + X3 + X5 + X6 + X7 + X1)

###### MATCHING ########################################################################################

# We now create another function for extracting the ACE using matching

get_ace_matching <- function(obs_df, variables) {
  out <- Match(Y = obs_df$Y, Tr = obs_df$Tr, X = obs_df[, variables], estimand = "ATE", 
               Weight = 1) # Takes the Mahalanobis distance between the covariates 
                           # to perform matching and extracts the ACE estimate
  return(out$est)
}

# We use sapply again to apply the function to our data for the different target sets
# to get the ACE estimates for the different target sets 

# Full set
ace_matching_full <- sapply(simulated_data_list, get_ace_matching,
                            variables = c("X1", "X2", "X3", "X4", "X5", "X6", 
                                          "X7", "X8", "X9"))

# Minimal set 
ace_matching_min <- sapply(simulated_data_list, get_ace_matching, variables = c("X3", "X5"))

# Treatment oriented
ace_matching_treat <- sapply(simulated_data_list, get_ace_matching,
                            variables = c("X1", "X2", "X3", "X4"))

# Outcome oriented
ace_matching_out <- sapply(simulated_data_list, get_ace_matching, 
                            variables = c("X3", "X5", "X6", "X7"))

######################## PS MATCHING ###############################################

# For propensity score matching, we first need to estimate the propensity scores using a logistic regression model. 
# These scores are then used to perform the PS score matching and extract the ACE. 

get_ace_ps <- function(data_frame, variables) {
  # We use paste variables so we dont have to write out the formula for each set
  formula_str <- paste(variables, collapse = " + ")
  # The formula is then used to fit a logistic regression model to estimate the propensity scores
  ps_model <- glm(as.formula(paste("Tr ~", formula_str)), data = data_frame, family = binomial("logit"))
  # The predicted propensity scores are extracted from the fitted model
  ps_scores <- predict(ps_model, type = "response") 
  # We use the estimated propensity scores to perform PS matching and extract the ACE using the Match function from the Matching package
  fit <- Match(Y = data_frame$Y, Tr = data_frame$Tr, X = ps_scores, estimand = "ATE", Weight = 1)
  # Finally, we return the estimated ACE from the matching procedure
  return(fit$est)
}

# Full set
ace_ps_full <- sapply(simulated_data_list, get_ace_ps,
                            variables = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"))

# Minimal set 
ace_ps_min <- sapply(simulated_data_list, get_ace_ps, variables = c("X3", "X5"))

# Treatment oriented
ace_ps_treat <- sapply(simulated_data_list, get_ace_ps,
                             variables = c("X1", "X2", "X3", "X4"))

# Outcome oriented
ace_ps_out <- sapply(simulated_data_list, get_ace_ps, 
                           variables = c("X3", "X5", "X6", "X7"))

########################## DOUBLY ROBUST ########################################################################################################

# Now we create a function for the doubly robust estimator which combines
# an outcome regression model and propensity score weighting to get the ACE estimates for the different target sets.

get_ace_dr <- function(data_frame, variables) {
  # Again we use paste so we don't have to write out the formula for each set
  formula_str <- paste(variables, collapse = " + ")
  # Propensity score model and outcome regression models require the same formula structure, so we create them using the same formula string
  form_Tr <- as.formula(paste("Tr ~", formula_str))
  form_Y <- as.formula(paste("Y ~", formula_str))
  
  # Propensity score model, what is the probability of treatment based on the covariates, which can be estimated using a logistic regression model
  ppi.glm <- glm(form_Tr, data = data_frame, family = binomial)
  # Extract the predicted propensity scores from the fitted model
  ppi.hat <- ppi.glm$fitted
  
  # Then we look at solely the treated and solely the control group to estimate the potential outcomes for each group using linear regression
  eta1.glm <- glm(form_Y, subset = (Tr == 1), data = data_frame)
  # Then we use predict to get the potential outcome predictions for all individuals in the data frame based on the fitted model for the treated group. For 
  # the observations not receiving treatment, this is their counterfactual outcome if they had received treatment, and for the people receiving treatment, this is their predicted outcome based on the model pattern for the treated group
  eta1.hat <- predict(eta1.glm, newdata = data_frame)
  
  # We now look at the untreated to estimate the potential outcomes for the control group using linear regression
  eta0.glm <- glm(form_Y, subset = (Tr == 0), data = data_frame)
  # Then again we use predict to get the predicted potential outcomes for all individuals in the data frame based on the fitted model for the untreated group. For the people receiving treatment, this is their counterfactual outcome if they had not received treatment, and for the people not receiving treatment, this is their predicted outcome for the model pattern 
  # for the untreated group
  eta0.hat <- predict(eta0.glm, newdata = data_frame)
  
  # Finally we combine all parts into the doubly robust estimator
  # We send the observed outcomes, treatment indicators, predicted propensity scores, and the predicted potential outcomes for both the treated and untreated groups 
  # to the ate.clik function from the drgee package to get the ACE estimate
  out.lik <- ate.clik(data_frame$Y, data_frame$Tr, ppi.hat, 
                      g0 = cbind(1, eta0.hat), 
                      g1 = cbind(1, eta1.hat))
  
  return(as.numeric(out.lik$diff))
}

# Full set
ace_dr_full <- sapply(simulated_data_list, get_ace_dr,
                      variables = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"))
# Minimal set
ace_dr_min <- sapply(simulated_data_list, get_ace_dr,
                     variables = c("X3", "X5"))
# Treatment oriented
ace_dr_treat <- sapply(simulated_data_list, get_ace_dr,
                       variables = c("X1", "X2", "X3", "X4"))
# Outcome oriented
ace_dr_out <- sapply(simulated_data_list, get_ace_dr,
                     variables = c("X3", "X5", "X6", "X7"))

######## COMBINING DATA ##########

# We combine all the ACE estimates from the different methods and variable sets into a single data frame for easier comparison and plotting.

n_sim <- length(simulated_data_list)
sets  <- c("X*", "Xmin", "XY", "XT")

all_estimates <- do.call(rbind, list(
  data.frame(ACE = c(ace_full, ace_min, ace_out, ace_treat),
             Method    = factor(rep(sets, each = n_sim), levels = sets),
             Estimator = "Linear Regression"),
  data.frame(ACE = c(ace_matching_full, ace_matching_min, ace_matching_out, ace_matching_treat),
             Method    = factor(rep(sets, each = n_sim), levels = sets),
             Estimator = "Matching"),
  data.frame(ACE = c(ace_ps_full, ace_ps_min, ace_ps_out, ace_ps_treat),
             Method    = factor(rep(sets, each = n_sim), levels = sets),
             Estimator = "PS Matching"),
  data.frame(ACE = c(ace_dr_full, ace_dr_min, ace_dr_out, ace_dr_treat),
             Method    = factor(rep(sets, each = n_sim), levels = sets),
             Estimator = "Doubly Robust")
))

######## PLOTTING ##########
png("plot_output.png", width = 15, height = 5, units = "in", res = 600)

# Define the specific order
target_order <- c("Linear Regression", "Matching", "PS Matching", "Doubly Robust")

# Convert Estimator to a factor with the specified levels
all_estimates$Estimator <- factor(all_estimates$Estimator, levels = target_order)

ggplot(all_estimates, aes(x = Method, y = ACE, fill = Method)) +
  geom_boxplot(outlier.size = 0.8, lwd = 0.5, alpha = 0.8, width = 0.6) + 
  facet_wrap(~Estimator, ncol = 4) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkblue", size = 0.5) +
  coord_cartesian(ylim = c(-1, 1.7)) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2")) +
  labs(x = NULL, y = "ACE") +
  theme_bw() + 
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    aspect.ratio = 1 
  )

dev.off()

###########################################################################################
# SUMMARY DATAFRAME — RMSE, Bias, Empirical SE (all 16 combinations)
###########################################################################################

compute_metrics <- function(estimates, true_val = 0.5) {
  bias  <- mean(estimates) - true_val
  emp_se <- sd(estimates)
  rmse  <- sqrt(mean((estimates - true_val)^2))
  c(Bias = bias, Emp_SE = emp_se, RMSE = rmse)
}

results_list <- list(
  list(method = "Linear Regression", set = "X*", ests = ace_full),
  list(method = "Linear Regression", set = "Xmin", ests = ace_min),
  list(method = "Linear Regression", set = "XT", ests = ace_treat),
  list(method = "Linear Regression", set = "XY", ests = ace_out),
  list(method = "Linear Regression", set = "XY + X1", ests = ace_out2),
  list(method = "Linear Regression", set = "Xmin + X1", ests = ace_min2),
  
  list(method = "Matching", set = "X*", ests = ace_matching_full),
  list(method = "Matching", set = "Xmin", ests = ace_matching_min),
  list(method = "Matching", set = "XT", ests = ace_matching_treat),
  list(method = "Matching", set = "XY", ests = ace_matching_out),
  
  list(method = "PS Matching", set = "X*", ests = ace_ps_full),
  list(method = "PS Matching", set = "Xmin", ests = ace_ps_min),
  list(method = "PS Matching", set = "XT", ests = ace_ps_treat),
  list(method = "PS Matching", set = "XY", ests = ace_ps_out),
  
  list(method = "Doubly Robust", set = "X*", ests = ace_dr_full),
  list(method = "Doubly Robust", set = "Xmin", ests = ace_dr_min),
  list(method = "Doubly Robust", set = "XT", ests = ace_dr_treat),
  list(method = "Doubly Robust", set = "XY", ests = ace_dr_out)
)

metrics_df <- do.call(rbind, lapply(results_list, function(x) {
  m <- compute_metrics(x$ests)
  data.frame(
    Method        = x$method,
    Variable_Set  = x$set,
    Bias          = round(m["Bias"], 4),
    Emp_SE        = round(m["Emp_SE"], 4),
    RMSE          = round(m["RMSE"], 4),
    stringsAsFactors = FALSE
  )
}))

rownames(metrics_df) <- NULL
print(metrics_df)

library(knitr)
library(kableExtra)

my_format <- if (is_latex_output()) "latex" else "html"

metrics_df %>%
  kable(
    format = my_format,
    booktabs = TRUE,
    col.names = c("Method", "Variable Set", "Bias", "Empirical SE", "RMSE"),
    align = "llccc",
    caption = "Estimation Performance Metrics by Method and Variable Set",
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    bootstrap_options = c("striped", "hover"), 
    full_width = FALSE
  ) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
  row_spec(0, bold = TRUE)

