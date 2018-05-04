library("oro.nifti")
library("parallel")
library("waldi")
library("brglm2")
library("dplyr")

n_cores <- 8

## Data are from
## https://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/software/bsglmm/manual/data_demo.dat.tar.gz
## See also https://github.com/nicholst/BSGLMM/releases/tag/v0.3
path <- "."

## Read the data
patient_characteristics <- read.table(paste0(path, "/", "lesion data/data_demo.dat"), header = TRUE)
n_patients <- nrow(patient_characteristics)
patient_data <- as.list(numeric(n_patients <- nrow(patient_characteristics)))
names(patient_data) <- patient_characteristics$image_file
for (file_name in patient_characteristics$image_file) {
    current_path <- paste0(path, "/", "lesion data/images/", file_name)
    patient_data[[file_name]] <- readNIfTI(current_path)
    cat("reding file", file_name, "done\n")
}

## Ensure that patient_data are ordered according to patient_characteristics
patient_data <- patient_data[patient_characteristics$image_file]

## Voxel-by-voxel organisation of the data
## number of logistic regressions we need to carry out:
n_voxels <- prod(dims <- dim(patient_data[[1]]@.Data))
## Greedy but OK at these resolutions
voxel_id <- expand.grid(d1 = seq.int(dims[1]),
                        d2 = seq.int(dims[2]),
                        d3 = seq.int(dims[3]))

lesions <- t(sapply(patient_data, function(x) c(x@.Data)))
rownames(lesions) <- patient_characteristics$patient


## get white matter (check in TN webpage)
white_matter <- readNIfTI(paste0(path, "/", "lesion data/images/", "avg152T1_white.nii.gz"))
non_white_matter_voxels <- which(white_matter == 0, arr.ind = TRUE)

## Identify unique lesion configurations of the 50 patients per voxel
## Again greedy but OK
lesions_string <- apply(lesions, 2, paste, collapse = "")
counts <- table(lesions_string)
unique_indicators <- match(names(counts), lesions_string)
array_indices <- match(lesions_string, names(counts))
unique_lesions <- lesions[, unique_indicators]
colnames(unique_lesions) <- names(counts)
n_unique_lesions <- length(counts)
## Must all be TRUE
## all.equal(names(counts), lesions_string[unique_indicators])
## all(apply(unique_lesions, 2, paste, collapse = "") == colnames(unique_lesions))
## all.equal(unique_lesions[, array_indices], lesions)

## Assumming that type 1 is relapsing-remitting MS, and type 2 is secondary progressive MS
base_formula <- voxel_lesion ~ type2 + age + sex + DD + EDSS + PASAT

## Massive univariate regression
link <- "probit"
fits <- mclapply(seq.int(n_unique_lesions), function(j) {
    current_data <- data.frame(voxel_lesion = unique_lesions[, j], patient_characteristics)
    fit_br <- glm(base_formula, family = binomial(link), data = current_data,
                  method = "brglm_fit", type = "AS_mean")
    fit_mbr <- glm(base_formula, family = binomial(link), data = current_data,
                   method = "brglm_fit", type = "AS_median")
    fit_ml <- glm(base_formula, family = binomial(link), data = current_data,
                  maxit = 100)
    coefs_ml <- coef(fit_ml)
    coefnames <- names(coefs_ml)[-1]
    mod <- glm(base_formula, family = binomial(link), data = current_data,
                   method = "detect_separation")
    is_inf <-  is.infinite(mod$beta)[coefnames]
    nvar <- length(coefnames)
    z_ml <- coef(summary(fit_ml))[coefnames, "z value"]
    z_br <- coef(summary(fit_br))[coefnames, "z value"]
    z_mbr <- coef(summary(fit_mbr))[coefnames, "z value"]
    corz_ml <- try(waldi(fit_ml, null = 0, what = NULL)[coefnames], silent = TRUE)
    corz_br <- try(waldi(fit_br, null = 0, what = NULL)[coefnames], silent = TRUE)
    corz_mbr <- try(waldi(fit_mbr, null = 0, what = NULL)[coefnames], silent = TRUE)
    if (inherits(corz_ml, "try-error")) {
        corz_ml <- rep(NA, nvar)
    }
    if (inherits(corz_br, "try-error")) {
        corz_br <- rep(NA, nvar)
    }
    if (inherits(corz_mbr, "try-error")) {
        corz_mbr <- rep(NA, nvar)
    }
    ## Give z_ml and corz_ml value zero if infinite estimate
    corz_ml[is_inf] <- 0
    z_ml[is_inf] <- 0
    r <- sapply(coefnames, function(variable) {
        f <- paste(deparse(base_formula), "-", variable)
        W <- anova(update(fit_ml, formula = f), fit_ml)$Deviance[2]
        out <- ifelse(W < 0, -Inf, sign(coefs_ml[variable]) * sqrt(W)) ## treat negative W values as infinities
        unname(out)
    })[coefnames]
    scores_br <- na.omit(fit_br$grad)
    scores_mbr <- na.omit(fit_mbr$grad)
    status_br <- all(abs(scores_br) < 1e-05)
    status_mbr <- all(abs(scores_mbr) < 1e-05)
    out <- data.frame(statistic = rep(c("z_ml", "z_br", "z_mbr", "corz_ml", "corz_br", "corz_mbr", "r"), each = nvar),
                      value = c(z_ml, z_br, z_mbr, corz_ml, corz_br, corz_mbr, r),
                      parameter = coefnames,
                      infinite = is_inf,
                      converged_br = status_br,
                      converged_mbr = status_mbr,
                      voxel = j,
                      voxel_id = paste(unique_lesions[, j], collapse = ""))
    if (j %% 100 == 0) cat(j, "\n")
    out
}, mc.cores = n_cores)

## Produce fits matrix
df_counts <- data.frame(voxel_id = names(counts), count = as.vector(unclass(counts)))
fits_mat <- left_join(do.call("rbind", fits), df_counts, by = "voxel_id")

save(lesions, fits_mat, white_matter, array_indices,
     file = paste(path, "results/brains_case_study.rda", sep = "/"))
