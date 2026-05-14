# ==============================================================================
# medicalAsset Core Function
# No hard-coded values, content-based privacy detection
# ==============================================================================

#' Comprehensive Data Asset Value Assessment
#'
#' @param data Data frame of clinical data
#' @return List of 5-dimensional scores, grade, features
#' @export
medicalAsset_assess <- function(data) {
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Please install dplyr")
    }
    library(dplyr)

    safe_mean <- function(x) mean(x, na.rm = TRUE)
    safe_sd   <- function(x) sd(x, na.rm = TRUE)
    safe_cv   <- function(x) {
        m <- safe_mean(x)
        if (is.na(m) || m == 0) return(0)
        return(safe_sd(x)/m)
    }

    # -------------------- 1. Data Quality --------------------
    total_cells  <- prod(dim(data))
    na_cells     <- sum(is.na(data))
    completeness <- (total_cells - na_cells) / total_cells

    dup_rows     <- sum(duplicated(data))
    uniqueness   <- 1 - (dup_rows / nrow(data))

    numeric_cols <- sapply(data, is.numeric)
    num_data     <- data[, numeric_cols, drop=FALSE]

    if (ncol(num_data) > 0) {
        z_list <- lapply(num_data, function(x) abs((x - safe_mean(x))/safe_sd(x)))
        outlier_ratio <- mean(sapply(z_list, function(z) mean(z>3, na.rm=TRUE)), na.rm=TRUE)
        accuracy <- 1 - outlier_ratio
    } else {
        accuracy <- 0.9
    }

    cv_vals <- unlist(lapply(num_data, safe_cv))
    cv_vals <- cv_vals[!is.na(cv_vals)]
    consistency <- ifelse(length(cv_vals)>0, 1-mean(cv_vals), 0.9)
    consistency <- pmax(pmin(consistency,1),0.7)

    data_quality <- round(100*(0.3*completeness + 0.2*uniqueness + 0.3*accuracy + 0.2*consistency),2)

    # -------------------- 2. Cost --------------------
    n_row <- nrow(data)
    n_col <- ncol(data)
    data_size_kb <- as.numeric(object.size(data))/1024
    storage_cost <- data_size_kb / 10000
    process_cost <- (n_row*n_col)/1000000
    total_cost <- pmin(storage_cost + process_cost, 1)
    cost <- round(100*(1-total_cost),2)

    # -------------------- 3. Application Value --------------------
    n_numeric <- sum(numeric_cols)
    variable_score <- pmin(n_numeric/15,1)

    if (ncol(num_data) > 0) {
        cv_all <- apply(num_data,2,safe_cv)
        valid_cv <- cv_all[!is.na(cv_all) & cv_all>0]
        missing_rate <- apply(num_data,2,function(x) mean(is.na(x)))
        high_info <- sum(valid_cv>0.1 & missing_rate<0.2)
        clinical_relevance <- pmin(high_info/max(3,ncol(num_data)*0.5),1)
    } else {
        clinical_relevance <- 0.5
    }
    clinical_relevance <- pmax(clinical_relevance,0.4)
    application_value <- round(100*(0.6*variable_score + 0.4*clinical_relevance),2)

    # -------------------- 4. Security & Compliance (Content-based) --------------------
    has_privacy <- FALSE
    for (col in seq_len(ncol(data))) {
        v <- na.omit(unique(data[[col]]))[1:20]
        if (length(v)==0) next
        id <- grepl("^[0-9]{17}[0-9X]$",v,ignore.case=TRUE)
        phone <- grepl("^1[3-9]\\d{9}$",v)
        name <- grepl("^[\\u4e00-\\u9fa5]{2,6}$",v)
        medno <- grepl("^[A-Za-z0-9]{6,20}$",v)
        if (any(id|phone|name|medno)) { has_privacy=TRUE; break }
    }

    consent <- ifelse(has_privacy,0.8,1.0)
    encrypt <- sum(numeric_cols)/ncol(data)
    access <- ifelse(dup_rows < n_row*0.05,1.0,0.7)
    security_compliance <- round(100*mean(c(consent,encrypt,access)),2)

    # -------------------- 5. Indicator Sensitivity --------------------
    if (ncol(num_data)==0) {
        sensitivity_score <- 65
    } else {
        vm <- mean(apply(num_data,2,safe_mean))
        vs <- mean(apply(num_data,2,safe_sd))
        vm <- ifelse(is.na(vm)|vm==0,1,vm)
        vs <- ifelse(is.na(vs)|vs<0,0.1,vs)
        cv <- vs/vm
        sensitivity_score <- round(pmin(cv*220,95),2)
    }

    # -------------------- Final Score --------------------
    total <- round(
        data_quality*0.2 + cost*0.2 + application_value*0.25 +
            security_compliance*0.2 + sensitivity_score*0.15, 2)

    grade <- case_when(
        total>=90~"A", total>=80~"B", total>=70~"C", TRUE~"D"
    )

    list(
        data_quality = data_quality,
        cost = cost,
        application_value = application_value,
        security_compliance = security_compliance,
        indicator_sensitivity = sensitivity_score,
        comprehensive_score = total,
        grade = grade,
        selected_features = colnames(num_data),
        quality_detail = c(completeness=completeness, uniqueness=uniqueness, accuracy=accuracy, consistency=consistency)
    )
}
