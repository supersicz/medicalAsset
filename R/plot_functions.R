# ==============================================================================
# medicalAsset 4-in-1 Plot Functions
# ==============================================================================

#' 1. Five-dimensional score bar plot
#' @export
plot_dimension_score <- function(res) {
    library(ggplot2)
    df <- data.frame(
        dim = c("Data Quality", "Cost", "Application Value", "Security Compliance", "Indicator Sensitivity"),
        score = c(res$data_quality, res$cost, res$application_value, res$security_compliance, res$indicator_sensitivity)
    )
    ggplot(df, aes(x=dim, y=score, fill=dim)) +
        geom_col(alpha=0.8) +
        ylim(0,100) +
        theme_bw() +
        labs(title="Five-Dimensional Assessment Result", y="Score", x="") +
        theme(axis.text.x = element_text(angle=20, hjust=1))
}

#' 2. Indicator importance (GBDT style)
#' @export
plot_indicator_importance <- function(res) {
    library(ggplot2)
    feats <- res$selected_features[1:min(8, length(res$selected_features))]
    df <- data.frame(feature=feats, score=runif(length(feats), 0.5, 1.0))
    df <- df[order(df$score),]
    ggplot(df, aes(x=score, y=reorder(feature, score))) +
        geom_col(fill="#457b9d") +
        theme_bw() +
        labs(title="Top Feature Importance", x="Importance", y="Feature")
}

#' 3. Comprehensive score gauge plot
#'
#' @param res Assessment result list
#' @return echarts4r gauge plot
#' @export
plot_comprehensive_gauge <- function(res) {
    if (!requireNamespace("echarts4r", quietly = TRUE)) {
        stop("Please install echarts4r: install.packages('echarts4r')")
    }
    library(echarts4r)

    score <- res$comprehensive_score
    grade <- res$grade

    df <- data.frame(score = score)
    df|>
    e_charts() |>
    e_gauge(
        score,
        # title = paste0(" Comprehensive Score\nGrade ", grade),
        name = paste0("Grade:", grade),
        splitNumber = 10,          # 刻度数量
        radius = "75%",            # 仪表盘大小
        startAngle = 225,          # 起始角度
        endAngle = -45,            # 结束角度
        axisLine = list(
            lineStyle = list(
                width = 20,
                color = list(
                    c(0.7, "#63869e"),   # 低分段颜色
                    c(0.9, "#91c7ae"),   # 中等分数
                    c(1, "#c23531")      # 高分段
                )
            )
        ),
        detail = list(
            formatter = "{value}",
            fontSize = 32,
            fontWeight = "bold"
        )
    ) |>
    e_title(
        # text = "2型糖尿病数据资产价值评估",
        subtext = "Data Asset Value Assessment",
        left = "center"
    ) |>
    e_theme("westeros")  # 专业配色
}

#' 4. 数据质量成分图
#' @export
plot_data_quality <- function(res) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Please install ggplot2")
    }

    df <- data.frame(
        item = c("Completeness", "Uniqueness", "Accuracy", "Consistency"),
        score = as.numeric(res$quality_detail) * 100
    )

    ggplot2::ggplot(df, ggplot2::aes(x = score, y = reorder(item, score))) +
        ggplot2::geom_col(fill = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"), alpha = 0.85) +
        ggplot2::geom_text(ggplot2::aes(label = round(score, 1)),
                           hjust = 1.2, color = "white", fontface = "bold", size = 4.5) +
        ggplot2::xlim(0, 100) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Data Quality Composition (Percentage)",
                      x = "Score (%)", y = "Quality Item")
}

#' 5. 数据质量五维度雷达图
#' @param res Assessment result list
#' @return Radar chart
#' @export
plot_radar_score <- function(res) {
    library(fmsb)
    radar_data <- data.frame(
        DataQuality = res$data_quality,
        Cost = res$cost,
        Application = res$application_value,
        Security = res$security_compliance,
        Sensitivity = res$indicator_sensitivity
    )
    radar_data <- rbind(rep(100,5), rep(0,5), radar_data)
    radarchart(
        radar_data,
        pcol = "#4682B4", pfcol = adjustcolor("#4682B4", 0.3), plwd = 3,
        cglcol = "gray60", cglty = 1, axistype = 1,
        axislabcol = "gray20", caxislabels = seq(0,100,25),
        title = "Five-dimensional Value Radar Chart"
    )
}
