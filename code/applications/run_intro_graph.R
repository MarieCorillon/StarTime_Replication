## 1. Data Preparation 
### Create the plots needed for the introduction

# 1. Silently load the data by suppressing messages, warnings, and column specifications
Variances_10min <- suppressWarnings(suppressMessages(
  read_csv("data/financial/Variances_10min.csv", show_col_types = FALSE, name_repair = "minimal")
))
variances_df <- Variances_10min

stock="AAPL"
max_lag=20

# 1. Get log-transformed realized volatility
log_rv <- log(variances_df[[stock]])
raw_block <- as.numeric(log_rv)

# 2. Standardize
mu <- mean(raw_block)
sigma <- sd(raw_block)
if (sigma == 0) sigma <- 1
scaled_block <- (raw_block - mu) / sigma

# 3. Create the lags
mat <- stats::embed(scaled_block, max_lag + 1)
y_all <- mat[, 1]
X_all <- mat[, -1, drop=FALSE]

# 4. HAR matrices
rv_lag1 <- X_all[, 1]
rv_week <- rowMeans(X_all[, 1:5, drop=FALSE])
rv_month <- rowMeans(X_all[, 1:20, drop=FALSE])
X_har <- cbind(rv_lag1, rv_week, rv_month)

# 5. HAR estimation
har_res_raw <- ols(X_har, y_all)
beta1 <- har_res_raw[1] + (har_res_raw[2]/5) + (har_res_raw[3]/20)
beta2 <- (har_res_raw[2]/5) + (har_res_raw[3]/20)
beta3 <- (har_res_raw[3]/20)
har_res <- c(beta1, rep(beta2,4), rep(beta3, 15))

# as.numeric(X_har[1, ]) %*% as.numeric(har_res_raw)

cvfit <- cv.glmnet(X_all, y_all, alpha=1)
lasso_res_raw <- glmnet(X_all, y_all, alpha = 1, lambda = cvfit$lambda.min)
lasso_res <- as.numeric(lasso_res_raw$beta)

P_list <- list(c(20, 4, 1))
Ki_list <- list(c(5, 4, 1))

# 6. Silently run StarTime estimations using capture.output to swallow prints/cats
invisible(capture.output(suppressMessages(suppressWarnings({
  post_startime_raw <- startime_ic(X_all, y_all, P_list, Ki_list, post=TRUE)
}))))
post_startime_res <- post_startime_raw$beta

invisible(capture.output(suppressMessages(suppressWarnings({
  simple_startime_raw <- startime_ic(X_all, y_all, P_list, Ki_list, post=FALSE)
}))))
simple_startime_res <- simple_startime_raw$beta


ts_data <- data.frame(
  Date = variances_df[[1]], # Using your Date column for the x-axis
  LogRV = raw_block
)

results_df <- data.frame(
  Lag = 1:20,
  HAR = har_res,
  Lasso = as.numeric(lasso_res), 
  Post_STARTIME = as.numeric(post_startime_res), 
  Simple_STARTIME = as.numeric(simple_startime_res)
)

plot_data <- results_df %>%
  pivot_longer(cols = -Lag, names_to = "Estimator_Raw", values_to = "Coefficient") %>%
  mutate(
    Estimator = factor(Estimator_Raw, 
                       levels = c("HAR", "Lasso", "Post_STARTIME", "Simple_STARTIME"),
                       labels = c("HAR", "Lasso", "Post StarTime", "Simple StarTime")),
    Color_Group = as.factor(round(Coefficient, 3)) 
  )

## 2. Your Original Custom Color Scheme
user_colors <- c("#648FFF", "#DC267F","#785EF0", "#FE6100", "#FFB000", "#4ECDC4", "#FF6F61")
n_groups <- length(unique(plot_data$Color_Group))
my_custom_ramp <- colorRampPalette(user_colors)(n_groups)

## 3. Build Top Plot (Time Series)
p_ts <- ggplot(ts_data, aes(x = Date, y = LogRV)) +
  geom_line(color = "#2C3E50", linewidth = 0.5) +
  labs(title = "Realized Variance: AAPL", y = "Log(RV)", x = NULL) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

## 4. Build Bottom Plot (Coefficients)
p_coef <- ggplot(plot_data, aes(x = Lag, y = Coefficient, fill = Color_Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~Estimator, ncol = 2, strip.position = "bottom") + 
  scale_fill_manual(values = my_custom_ramp) + 
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  labs(y = "Coefficient Magnitude", x = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "white", color = NA),
    strip.placement = "outside", 
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

## 5. Combine using Patchwork 
final_plot <- p_ts / p_coef + plot_layout(heights = c(1, 1.8))

out_path <- file.path(FIGURES_DIR, "applications/financial/intro_graph.pdf")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

# Save the plot
suppressMessages(
  ggsave(out_path, final_plot, width = 10, height = 8, dpi = 300)
)

# Print the final success message!
message("   ✓ Introductory graph generated and saved to: ", out_path)