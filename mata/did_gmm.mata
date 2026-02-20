*! did_gmm.mata - GMM optimal weighting and double DID estimation
*!
*! Implements the generalized method of moments (GMM) framework for combining
*! the standard DID estimator and the sequential DID estimator into the double
*! difference-in-differences (double DID) estimator.
*!
*! GMM Objective:
*!   tau_dDID = argmin (tau - tau_DID, tau - tau_sDID)' W (tau - tau_DID, tau - tau_sDID)
*!
*! Optimal Weight Computation:
*!   W = Sigma^{-1}                             (precision matrix)
*!   w_DID  = (W[1,1] + W[1,2]) / sum(W)        (weight for DID estimator)
*!   w_sDID = (W[2,2] + W[1,2]) / sum(W)        (weight for sequential DID)
*!   tau_dDID = w_DID * tau_DID + w_sDID * tau_sDID
*!   Var(tau_dDID) = 1 / sum(W)                 (asymptotic variance)

version 16.0

mata:
mata set matastrict on

// ============================================================================
// GMM OPTIMAL WEIGHTING AND DOUBLE DID ESTIMATION
// ============================================================================
// This module implements the GMM framework for double DID estimation:
//   1. Optimal weight matrix computation from bootstrap variance-covariance
//   2. Double DID point estimation combining DID and sequential DID
//   3. Staggered adoption design extensions with time-weighted aggregation
// ============================================================================

// ----------------------------------------------------------------------------
// DATA STRUCTURES
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * struct gmm_weights - GMM Optimal Weight Structure
 *
 * Stores the precision matrix and derived optimal weights for combining
 * the DID and sequential DID estimators via GMM.
 *---------------------------------------------------------------------------*/
struct gmm_weights {
    real matrix W                // W = Sigma^{-1}, precision matrix (2x2)
    real matrix vcov             // Sigma, bootstrap VCOV of (tau_DID, tau_sDID)
    real scalar w_did            // Optimal weight for DID estimator
    real scalar w_sdid           // Optimal weight for sequential DID estimator
}

/*---------------------------------------------------------------------------
 * struct ddid_result - Double DID Estimation Result
 *
 * Stores complete results from double DID estimation including point
 * estimates, variance, confidence intervals, and GMM weights.
 *---------------------------------------------------------------------------*/
struct ddid_result {
    real scalar estimate         // tau_dDID, double DID point estimate
    real scalar variance         // Var(tau_dDID), asymptotic or bootstrap
    real scalar std_error        // SE(tau_dDID) = sqrt(variance)
    real scalar ci_low           // Confidence interval lower bound
    real scalar ci_high          // Confidence interval upper bound
    real scalar w_did            // GMM weight for DID estimator
    real scalar w_sdid           // GMM weight for sequential DID estimator
    real scalar tau_did          // tau_DID, standard DID point estimate
    real scalar tau_sdid         // tau_sDID, sequential DID point estimate
    real scalar var_did          // Var(tau_DID) from bootstrap
    real scalar var_sdid         // Var(tau_sDID) from bootstrap
    real scalar lead             // Lead value (post-treatment periods ahead)
}

/*---------------------------------------------------------------------------
 * struct sa_ddid_result - Staggered Adoption Double DID Result
 *
 * Stores complete results for staggered adoption double DID estimation
 * across multiple lead values, including GMM weights and precision matrices.
 *---------------------------------------------------------------------------*/
struct sa_ddid_result {
    real colvector estimate      // tau_dDID, double DID estimates (n_lead x 1)
    real colvector variance      // Var(tau_dDID) from bootstrap (n_lead x 1)
    real colvector std_error     // SE(tau_dDID) = sqrt(variance) (n_lead x 1)
    real colvector ci_low        // CI lower bounds (n_lead x 1)
    real colvector ci_high       // CI upper bounds (n_lead x 1)
    real colvector w_did         // GMM weights for DID (n_lead x 1)
    real colvector w_sdid        // GMM weights for sequential DID (n_lead x 1)
    real colvector tau_did       // tau_DID estimates (n_lead x 1)
    real colvector tau_sdid      // tau_sDID estimates (n_lead x 1)
    real colvector var_did       // Var(tau_DID) from bootstrap (n_lead x 1)
    real colvector var_sdid      // Var(tau_sDID) from bootstrap (n_lead x 1)
    pointer vector W_matrices    // Precision matrices W, one per lead
    pointer vector VCOV_matrices // Bootstrap VCOV matrices Sigma, one per lead
}

/*---------------------------------------------------------------------------
 * struct sa_ddid_var_result - Staggered Adoption Variance Result
 *
 * Stores bootstrap variance and percentile confidence intervals for
 * double DID, DID, and sequential DID estimators at a specific lead.
 *---------------------------------------------------------------------------*/
struct sa_ddid_var_result {
    real scalar var              // Bootstrap variance of double DID
    real rowvector ci_low        // Percentile CI lower bounds (dDID, DID, sDID)
    real rowvector ci_high       // Percentile CI upper bounds (dDID, DID, sDID)
}

/*---------------------------------------------------------------------------
 * struct sa_point - Staggered Adoption Point Estimates
 *
 * Stores time-weighted DID and sequential DID estimates for each lead
 * value in staggered adoption designs.
 *---------------------------------------------------------------------------*/
struct sa_point {
    real rowvector DID           // Time-weighted tau_DID for each lead (1 x n_lead)
    real rowvector sDID          // Time-weighted tau_sDID for each lead (1 x n_lead)
}

// ----------------------------------------------------------------------------
// GMM WEIGHT COMPUTATION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * compute_weights() - Compute Optimal GMM Weight Matrix
 *
 * Computes the optimal weight matrix W and derived scalar weights from
 * the bootstrap variance-covariance matrix of (tau_DID, tau_sDID).
 *
 * Arguments:
 *   vcov : real matrix (2x2) - Bootstrap VCOV matrix Sigma
 *
 * Returns:
 *   struct gmm_weights - Contains W, vcov, w_did, w_sdid
 *
 * GMM Weight Computation:
 *   W = Sigma^{-1}                                (precision matrix)
 *   w_DID  = (W[1,1] + W[1,2]) / sum(W)           (optimal DID weight)
 *   w_sDID = (W[2,2] + W[1,2]) / sum(W)           (optimal sDID weight)
 *
 * Closed-form equivalence:
 *   w_DID  = (Var(sDID) - Cov) / (Var(DID) + Var(sDID) - 2*Cov)
 *   w_sDID = (Var(DID) - Cov) / (Var(DID) + Var(sDID) - 2*Cov)
 *
 * Property: w_DID + w_sDID = 1.0 (convex combination)
 *---------------------------------------------------------------------------*/
struct gmm_weights scalar compute_weights(real matrix vcov)
{
    struct gmm_weights scalar weights
    real matrix W
    real scalar sum_W, success
    real rowvector eigs
    real scalar cond_num
    
    // Initialize result structure with input VCOV and missing weights
    weights.vcov = vcov
    weights.W = J(2, 2, .)
    weights.w_did = .
    weights.w_sdid = .
    
    // -------------------------------------------------------------------------
    // Validate input dimensions
    // -------------------------------------------------------------------------
    if (rows(vcov) != 2 || cols(vcov) != 2) {
        errprintf("Error: VCOV must be a 2x2 matrix\n")
        return(weights)
    }
    
    // -------------------------------------------------------------------------
    // Validate variance positivity (diagonal elements must be non-negative)
    // -------------------------------------------------------------------------
    if (vcov[1,1] < 0 || vcov[2,2] < 0) {
        errprintf("Error: Negative variance detected\n")
        errprintf("  Var(DID) = %g, Var(sDID) = %g\n", vcov[1,1], vcov[2,2])
        return(weights)
    }
    
    // -------------------------------------------------------------------------
    // Assess numerical stability via condition number
    // A high condition number indicates near-singularity
    // -------------------------------------------------------------------------
    eigs = symeigenvalues(vcov)
    if (min(eigs) > 0) {
        cond_num = max(eigs) / min(eigs)
        if (cond_num > 1e10) {
            printf("{txt}Warning: VCOV matrix is near-singular (cond=%g), results may be unreliable\n", cond_num)
        }
        else if (cond_num > 1e8) {
            printf("{txt}Note: VCOV matrix condition number is moderately large (cond=%g)\n", cond_num)
        }
    }
    
    // -------------------------------------------------------------------------
    // Compute precision matrix W = Sigma^{-1}
    // Tolerance 1e-10 used for numerical stability in matrix inversion
    // -------------------------------------------------------------------------
    W = safe_invert(vcov, 1e-10, &success)
    
    // -------------------------------------------------------------------------
    // Handle singular matrix case
    // -------------------------------------------------------------------------
    if (success == 0 || missing(W[1,1])) {
        errprintf("{err}Error: Variance-covariance matrix is singular, cannot compute GMM weights\n")
        errprintf("{err}       This may be caused by:\n")
        errprintf("{err}       - Insufficient bootstrap samples\n")
        errprintf("{err}       - Collinear data\n")
        errprintf("{err}       - Extreme outliers\n")
        weights.W = J(2, 2, .)
        weights.w_did = .
        weights.w_sdid = .
        return(weights)
    }
    
    // -------------------------------------------------------------------------
    // Compute weight sum: sum(W) = W[1,1] + W[2,2] + 2*W[1,2]
    // -------------------------------------------------------------------------
    sum_W = sum(W)
    
    // -------------------------------------------------------------------------
    // Compute optimal GMM weights
    // -------------------------------------------------------------------------
    weights.W = W
    if (sum_W <= 1e-10 | missing(sum_W)) {
        errprintf("Warning: sum(W) is zero or near-zero, using equal weights\n")
        weights.w_did = 0.5
        weights.w_sdid = 0.5
    }
    else {
        weights.w_did = (W[1,1] + W[1,2]) / sum_W
        weights.w_sdid = (W[2,2] + W[1,2]) / sum_W
    }
    
    // -------------------------------------------------------------------------
    // Verify convex combination property: w_DID + w_sDID = 1
    // -------------------------------------------------------------------------
    if (abs(weights.w_did + weights.w_sdid - 1.0) > 1e-10) {
        printf("{txt}Warning: GMM weights do not sum to 1.0 (sum = %18.15f)\n", 
               weights.w_did + weights.w_sdid)
    }
    
    // -------------------------------------------------------------------------
    // Check for weights outside [0,1] (occurs with high positive correlation)
    // -------------------------------------------------------------------------
    if (weights.w_did < 0 | weights.w_did > 1 | weights.w_sdid < 0 | weights.w_sdid > 1) {
        printf("{txt}Warning: GMM weights are outside [0,1] range (w_did=%g, w_sdid=%g)\n", 
               weights.w_did, weights.w_sdid)
        printf("{txt}         This may indicate high positive correlation between DID and sDID estimates\n")
    }
    
    return(weights)
}

// ----------------------------------------------------------------------------
// DOUBLE DID ESTIMATION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * compute_double_did() - Compute Double DID Point Estimate and Inference
 *
 * Combines the DID and sequential DID estimators using optimal GMM weights
 * to produce the double DID estimator with variance and confidence intervals.
 *
 * Arguments:
 *   tau_did   : real scalar - Standard DID point estimate
 *   tau_sdid  : real scalar - Sequential DID point estimate
 *   weights   : struct gmm_weights - Optimal weights from compute_weights()
 *   boot_est  : real matrix (B x 2) - Bootstrap estimates [tau_DID, tau_sDID] (optional)
 *   se_boot   : real scalar - Inference method: 0 = asymptotic, 1 = bootstrap
 *   level     : real scalar - Confidence level in percent (default: 95)
 *
 * Returns:
 *   struct ddid_result - Complete double DID results
 *
 * Point Estimation:
 *   tau_dDID = w_DID * tau_DID + w_sDID * tau_sDID
 *
 * Asymptotic Inference (se_boot = 0):
 *   Var(tau_dDID) = 1 / sum(W)
 *   CI: tau_dDID +/- z_{alpha/2} * sqrt(Var)
 *
 * Bootstrap Inference (se_boot = 1):
 *   tau_dDID^{(b)} = w_DID * tau_DID^{(b)} + w_sDID * tau_sDID^{(b)}
 *   Var(tau_dDID) = sample variance of tau_dDID^{(b)}
 *   CI: percentile method at (alpha/2, 1-alpha/2)
 *---------------------------------------------------------------------------*/
struct ddid_result scalar compute_double_did(real scalar tau_did,
                                             real scalar tau_sdid,
                                             struct gmm_weights scalar weights,
                                             | real matrix boot_est,
                                             real scalar se_boot,
                                             real scalar level)
{
    struct ddid_result scalar result
    real scalar tau_ddid, var_ddid, var_did, var_sdid
    real scalar alpha, z, n_boot
    real colvector boot_ddid
    real scalar w_did, w_sdid
    
    // -------------------------------------------------------------------------
    // Set default parameters
    // -------------------------------------------------------------------------
    if (args() < 5) se_boot = 0
    if (args() < 6) level = 95
    
    // Validate confidence level
    if (level <= 0 || level >= 100) {
        printf("{txt}Warning: Invalid confidence level %g%%, using 95%%\n", level)
        level = 95
    }
    
    // -------------------------------------------------------------------------
    // Extract weights
    // -------------------------------------------------------------------------
    w_did = weights.w_did
    w_sdid = weights.w_sdid
    
    // Return missing result if weights are invalid
    if (missing(w_did) || missing(w_sdid)) {
        result.estimate = .
        result.variance = .
        result.std_error = .
        result.ci_low = .
        result.ci_high = .
        result.w_did = .
        result.w_sdid = .
        result.tau_did = tau_did
        result.tau_sdid = tau_sdid
        result.var_did = .
        result.var_sdid = .
        result.lead = .
        return(result)
    }
    
    // -------------------------------------------------------------------------
    // Compute double DID point estimate as weighted combination
    // -------------------------------------------------------------------------
    tau_ddid = w_did * tau_did + w_sdid * tau_sdid
    
    // -------------------------------------------------------------------------
    // Extract component variances from bootstrap VCOV diagonal
    // -------------------------------------------------------------------------
    var_did = weights.vcov[1,1]
    var_sdid = weights.vcov[2,2]
    
    // -------------------------------------------------------------------------
    // Compute critical value for confidence interval
    // -------------------------------------------------------------------------
    alpha = 1 - level / 100
    z = invnormal(1 - alpha / 2)
    
    // -------------------------------------------------------------------------
    // Select inference method based on se_boot flag and bootstrap availability
    // -------------------------------------------------------------------------
    n_boot = (args() >= 4 && rows(boot_est) > 0) ? rows(boot_est) : 0
    
    if (se_boot == 0 || n_boot < 2) {
        // ---------------------------------------------------------------------
        // Asymptotic inference: variance derived from precision matrix
        // ---------------------------------------------------------------------
        {
            real scalar sum_W
            sum_W = sum(weights.W)
            if (sum_W <= 1e-10 | missing(sum_W)) {
                var_ddid = .
                result.ci_low = .
                result.ci_high = .
            }
            else {
                var_ddid = 1 / sum_W
                result.ci_low = tau_ddid - z * sqrt(var_ddid)
                result.ci_high = tau_ddid + z * sqrt(var_ddid)
            }
        }
    } 
    else {
        // ---------------------------------------------------------------------
        // Bootstrap inference: variance and CI from bootstrap distribution
        // ---------------------------------------------------------------------
        
        // Apply GMM weights to bootstrap estimates
        boot_ddid = w_did * boot_est[., 1] + w_sdid * boot_est[., 2]
        
        // Exclude missing values before variance computation
        real colvector boot_ddid_valid
        real scalar n_valid
        boot_ddid_valid = select(boot_ddid, boot_ddid :< .)
        n_valid = rows(boot_ddid_valid)
        
        // Sample variance with Bessel correction (B-1 denominator)
        if (n_valid >= 2) {
            var_ddid = variance(boot_ddid_valid)
        }
        else {
            var_ddid = .
        }
        
        // Percentile confidence interval from bootstrap distribution
        result.ci_low = quantile_sorted(boot_ddid, alpha / 2)
        result.ci_high = quantile_sorted(boot_ddid, 1 - alpha / 2)
    }
    
    // -------------------------------------------------------------------------
    // Populate result structure with all computed values
    // -------------------------------------------------------------------------
    result.estimate = tau_ddid
    result.variance = var_ddid
    result.std_error = sqrt(var_ddid)
    result.w_did = w_did
    result.w_sdid = w_sdid
    result.tau_did = tau_did
    result.tau_sdid = tau_sdid
    result.var_did = var_did
    result.var_sdid = var_sdid
    result.lead = .
    
    return(result)
}

// ----------------------------------------------------------------------------
// STAGGERED ADOPTION DESIGN FUNCTIONS
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * sa_to_ddid() - Staggered Adoption Double DID Estimation
 *
 * Computes the double DID estimator for staggered adoption designs across
 * multiple lead values using GMM optimal weighting at each lead.
 *
 * Arguments:
 *   point_est : struct sa_point - Time-weighted DID and sDID estimates
 *   boot_est  : pointer vector - Bootstrap sa_point structures (B elements)
 *   lead      : real rowvector - Lead values, e.g., (0) or (0, 1, 2)
 *   level     : real scalar - Confidence level in percent (default: 95)
 *
 * Returns:
 *   struct sa_ddid_result - Complete results for all lead values
 *
 * Algorithm (for each lead value l):
 *   1. Compute bootstrap VCOV: Sigma_l = sa_calc_cov(boot_est, l)
 *   2. Compute precision matrix: W_l = Sigma_l^{-1}
 *   3. Compute GMM weights:
 *        w_DID  = (W[1,1] + W[1,2]) / sum(W)
 *        w_sDID = (W[2,2] + W[1,2]) / sum(W)
 *   4. Compute double DID: tau_dDID = w_DID * tau_DID + w_sDID * tau_sDID
 *   5. Compute bootstrap variance and percentile CI via sa_calc_ddid_var()
 *---------------------------------------------------------------------------*/
struct sa_ddid_result scalar sa_to_ddid(struct sa_point scalar point_est,
                                        pointer(struct sa_point scalar) vector boot_est,
                                        real rowvector lead,
                                        real scalar level)
{
    struct sa_ddid_result scalar result
    struct sa_ddid_var_result scalar var_result
    real matrix VC, W
    real scalar n_lead, n_boot, ll, lead_idx
    real scalar w_did, w_sdid, sum_W
    real scalar tau_did, tau_sdid, tau_ddid
    real scalar var_did, var_sdid, var_ddid
    real scalar success
    
    // -------------------------------------------------------------------------
    // Set default parameters and validate inputs
    // -------------------------------------------------------------------------
    if (args() < 4) level = 95
    
    n_lead = cols(lead)
    n_boot = length(boot_est)
    
    // -------------------------------------------------------------------------
    // Initialize result structure
    // -------------------------------------------------------------------------
    result.estimate = J(n_lead, 1, .)
    result.variance = J(n_lead, 1, .)
    result.std_error = J(n_lead, 1, .)
    result.ci_low = J(n_lead, 1, .)
    result.ci_high = J(n_lead, 1, .)
    result.w_did = J(n_lead, 1, .)
    result.w_sdid = J(n_lead, 1, .)
    result.tau_did = J(n_lead, 1, .)
    result.tau_sdid = J(n_lead, 1, .)
    result.var_did = J(n_lead, 1, .)
    result.var_sdid = J(n_lead, 1, .)
    result.W_matrices = J(n_lead, 1, NULL)
    result.VCOV_matrices = J(n_lead, 1, NULL)
    
    // -------------------------------------------------------------------------
    // Validate bootstrap sample availability
    // -------------------------------------------------------------------------
    if (n_boot == 0) {
        errprintf("sa_to_ddid(): No Bootstrap samples\n")
        return(result)
    }
    
    if (n_boot < 2) {
        errprintf("sa_to_ddid(): Cannot compute variance with single Bootstrap sample\n")
        return(result)
    }
    
    // -------------------------------------------------------------------------
    // Loop over each lead value
    // -------------------------------------------------------------------------
    for (ll = 1; ll <= n_lead; ll++) {
        
        // Lead index for bootstrap access (1-based)
        lead_idx = ll
        
        // ---------------------------------------------------------------------
        // Get point estimates for this lead
        // ---------------------------------------------------------------------
        tau_did = point_est.DID[ll]
        tau_sdid = point_est.sDID[ll]
        
        result.tau_did[ll] = tau_did
        result.tau_sdid[ll] = tau_sdid
        
        // ---------------------------------------------------------------------
        // Compute bootstrap VCOV matrix for this lead
        // ---------------------------------------------------------------------
        result.VCOV_matrices[ll] = &(sa_calc_cov(boot_est, lead_idx))
        VC = *result.VCOV_matrices[ll]
        
        // Handle invalid VCOV (missing values indicate bootstrap failure)
        if (missing(VC[1,1]) || missing(VC[2,2])) {
            errprintf("{err}Warning: Failed to compute VCOV for lead %g\n", lead[ll])
            errprintf("{err}         Bootstrap VCOV contains missing values\n")
            errprintf("{err}         This may be caused by insufficient valid bootstrap iterations\n")
            result.estimate[ll] = .
            result.variance[ll] = .
            result.std_error[ll] = .
            result.ci_low[ll] = .
            result.ci_high[ll] = .
            result.w_did[ll] = .
            result.w_sdid[ll] = .
            result.var_did[ll] = .
            result.var_sdid[ll] = .
            continue
        }
        
        // Store variances from VCOV diagonal
        var_did = VC[1,1]
        var_sdid = VC[2,2]
        result.var_did[ll] = var_did
        result.var_sdid[ll] = var_sdid
        
        // ---------------------------------------------------------------------
        // Compute precision matrix W = Sigma^{-1} for GMM weighting
        // ---------------------------------------------------------------------
        result.W_matrices[ll] = &(safe_invert(VC, 1e-10, &success))
        W = *result.W_matrices[ll]
        
        // Handle singular precision matrix
        if (success == 0 || missing(W[1,1])) {
            errprintf("{err}Warning: Variance-covariance matrix is singular for lead %g\n", lead[ll])
            errprintf("{err}         Cannot compute GMM weights for SA design\n")
            errprintf("{err}         This may be caused by:\n")
            errprintf("{err}         - Insufficient bootstrap samples for this lead\n")
            errprintf("{err}         - Collinear data at this lead value\n")
            errprintf("{err}         - Insufficient variation in treatment timing\n")
            result.estimate[ll] = .
            result.variance[ll] = .
            result.std_error[ll] = .
            result.ci_low[ll] = .
            result.ci_high[ll] = .
            result.w_did[ll] = .
            result.w_sdid[ll] = .
            continue
        }
        
        // ---------------------------------------------------------------------
        // Compute GMM weights
        // w_did = (W[1,1] + W[1,2]) / sum(W)
        // w_sdid = (W[2,2] + W[1,2]) / sum(W)
        // ---------------------------------------------------------------------
        sum_W = sum(W)
        if (sum_W <= 1e-10 | missing(sum_W)) {
            errprintf("Warning: sum(W) is zero or near-zero for lead %g, using equal weights\n", lead[ll])
            w_did = 0.5
            w_sdid = 0.5
        }
        else {
            w_did = (W[1,1] + W[1,2]) / sum_W
            w_sdid = (W[2,2] + W[1,2]) / sum_W
        }
        
        // Verify convex combination property
        if (abs(w_did + w_sdid - 1.0) > 1e-10) {
            printf("{txt}Warning: GMM weights do not sum to 1.0 for lead %g (sum = %18.15f)\n", 
                   lead[ll], w_did + w_sdid)
        }
        
        // Warn if weights are outside [0,1] range
        if (w_did < 0 | w_did > 1 | w_sdid < 0 | w_sdid > 1) {
            printf("{txt}Warning: GMM weights outside [0,1] for lead %g (w_did=%g, w_sdid=%g)\n", 
                   lead[ll], w_did, w_sdid)
            printf("{txt}         This may indicate high positive correlation between DID and sDID estimates\n")
        }
        
        result.w_did[ll] = w_did
        result.w_sdid[ll] = w_sdid
        
        // ---------------------------------------------------------------------
        // Compute SA-Double-DID point estimate
        // tau_dDID = w_did * tau_DID + w_sdid * tau_sDID
        // ---------------------------------------------------------------------
        tau_ddid = w_did * tau_did + w_sdid * tau_sdid
        result.estimate[ll] = tau_ddid
        
        // ---------------------------------------------------------------------
        // Compute bootstrap variance and CI
        // ---------------------------------------------------------------------
        var_result = sa_calc_ddid_var(boot_est, lead_idx, w_did, w_sdid, level)
        
        // Store variance and standard error
        var_ddid = var_result.var
        result.variance[ll] = var_ddid
        result.std_error[ll] = sqrt(var_ddid)
        
        // Store CI bounds
        result.ci_low[ll] = var_result.ci_low[1]
        result.ci_high[ll] = var_result.ci_high[1]
    }
    
    return(result)
}

/*---------------------------------------------------------------------------
 * sa_calc_cov() - Compute Bootstrap VCOV for Staggered Adoption
 *
 * Extracts bootstrap DID and sequential DID estimates for a specific lead
 * index and computes the 2x2 variance-covariance matrix Sigma.
 *
 * Arguments:
 *   boot_est : pointer vector - Bootstrap sa_point structures (B elements)
 *   lead_idx : real scalar - Lead index (1-based)
 *
 * Returns:
 *   real matrix (2x2) - Bootstrap VCOV matrix Sigma:
 *     [Var(tau_DID),          Cov(tau_DID, tau_sDID)]
 *     [Cov(tau_DID, tau_sDID), Var(tau_sDID)        ]
 *---------------------------------------------------------------------------*/
real matrix sa_calc_cov(pointer(struct sa_point scalar) vector boot_est, real scalar lead_idx)
{
    real matrix combined, combined_valid, vcov
    real scalar n_boot, b, n_valid
    real colvector valid_idx
    struct sa_point scalar pt
    
    n_boot = length(boot_est)
    
    // -------------------------------------------------------------------------
    // Require at least 2 bootstrap samples for variance estimation
    // -------------------------------------------------------------------------
    if (n_boot == 0) {
        return(J(2, 2, .))
    }
    
    if (n_boot < 2) {
        return(J(2, 2, .))
    }
    
    // -------------------------------------------------------------------------
    // Extract bootstrap estimates for specified lead index
    // -------------------------------------------------------------------------
    combined = J(n_boot, 2, .)
    valid_idx = J(n_boot, 1, 0)
    
    for (b = 1; b <= n_boot; b++) {
        // Skip null pointers from failed bootstrap iterations
        if (boot_est[b] == NULL) {
            continue
        }
        
        // Dereference pointer to access bootstrap estimates
        pt = *boot_est[b]
        
        // Validate lead_idx bounds before vector access
        if (lead_idx < 1 || lead_idx > cols(pt.DID) || lead_idx > cols(pt.sDID)) {
            continue
        }
        
        // Extract bootstrap estimates for this lead
        combined[b, 1] = pt.DID[lead_idx]
        combined[b, 2] = pt.sDID[lead_idx]
        
        // Mark as valid only if both estimates are non-missing
        if (!missing(combined[b, 1]) && !missing(combined[b, 2])) {
            valid_idx[b] = 1
        }
    }
    
    // -------------------------------------------------------------------------
    // Exclude missing entries before VCOV computation
    // -------------------------------------------------------------------------
    n_valid = sum(valid_idx)
    
    if (n_valid < 2) {
        return(J(2, 2, .))
    }
    
    combined_valid = select(combined, valid_idx)
    
    // -------------------------------------------------------------------------
    // Compute sample VCOV with Bessel correction (n-1 denominator)
    // -------------------------------------------------------------------------
    vcov = compute_vcov(combined_valid)
    
    return(vcov)
}

/*---------------------------------------------------------------------------
 * sa_calc_ddid_var() - Compute Bootstrap Variance and CI for SA Double DID
 *
 * Computes bootstrap variance and percentile confidence intervals for the
 * double DID estimator and its component estimators (DID and sequential DID).
 *
 * Arguments:
 *   boot_est : pointer vector - Bootstrap sa_point structures (B elements)
 *   lead_idx : real scalar - Lead index (1-based)
 *   w_did    : real scalar - GMM weight for DID estimator
 *   w_sdid   : real scalar - GMM weight for sequential DID estimator
 *   level    : real scalar - Confidence level in percent (default: 95)
 *
 * Returns:
 *   struct sa_ddid_var_result:
 *     var     : Bootstrap variance of SA-Double-DID
 *     ci_low  : (1x3) lower CI bounds for (dDID, DID, sDID)
 *     ci_high : (1x3) upper CI bounds for (dDID, DID, sDID)
 *---------------------------------------------------------------------------*/
struct sa_ddid_var_result scalar sa_calc_ddid_var(pointer(struct sa_point scalar) vector boot_est,
                                                   real scalar lead_idx,
                                                   real scalar w_did,
                                                   real scalar w_sdid,
                                                   real scalar level)
{
    struct sa_ddid_var_result scalar result
    struct sa_point scalar pt
    real colvector boot_ddid, boot_did, boot_sdid
    real scalar n_boot, b, alpha
    
    // -------------------------------------------------------------------------
    // Set default parameters
    // -------------------------------------------------------------------------
    if (args() < 5) level = 95
    
    // Initialize result with missing values
    result.var = .
    result.ci_low = (., ., .)
    result.ci_high = (., ., .)
    
    n_boot = length(boot_est)
    
    // -------------------------------------------------------------------------
    // Require at least 2 bootstrap samples for variance estimation
    // -------------------------------------------------------------------------
    if (n_boot == 0) {
        return(result)
    }
    
    if (n_boot < 2) {
        return(result)
    }
    
    // -------------------------------------------------------------------------
    // Extract bootstrap estimates and apply GMM weights
    // -------------------------------------------------------------------------
    boot_ddid = J(n_boot, 1, .)
    boot_did = J(n_boot, 1, .)
    boot_sdid = J(n_boot, 1, .)
    
    for (b = 1; b <= n_boot; b++) {
        // Skip null pointers from failed bootstrap iterations
        if (boot_est[b] == NULL) {
            continue
        }
        
        // Dereference pointer to access bootstrap estimates
        pt = *boot_est[b]
        
        // Validate lead_idx bounds before vector access
        if (lead_idx < 1 || lead_idx > cols(pt.DID) || lead_idx > cols(pt.sDID)) {
            continue
        }
        
        // Extract bootstrap estimates for this lead
        boot_did[b] = pt.DID[lead_idx]
        boot_sdid[b] = pt.sDID[lead_idx]
        
        // Apply GMM weights to compute bootstrap double DID estimate
        boot_ddid[b] = w_did * boot_did[b] + w_sdid * boot_sdid[b]
    }
    
    // -------------------------------------------------------------------------
    // Exclude missing values before variance computation
    // -------------------------------------------------------------------------
    real colvector boot_ddid_valid, boot_did_valid, boot_sdid_valid
    real scalar n_valid
    
    boot_ddid_valid = select(boot_ddid, boot_ddid :< .)
    n_valid = rows(boot_ddid_valid)
    
    // -------------------------------------------------------------------------
    // Compute sample variance with Bessel correction (n-1 denominator)
    // -------------------------------------------------------------------------
    if (n_valid >= 2) {
        result.var = variance(boot_ddid_valid)
    }
    
    // -------------------------------------------------------------------------
    // Compute percentile confidence intervals from bootstrap distribution
    // -------------------------------------------------------------------------
    alpha = 1 - level / 100
    
    // Double DID percentile CI
    result.ci_low[1] = quantile_sorted(boot_ddid, alpha / 2)
    result.ci_high[1] = quantile_sorted(boot_ddid, 1 - alpha / 2)
    
    // DID percentile CI
    result.ci_low[2] = quantile_sorted(boot_did, alpha / 2)
    result.ci_high[2] = quantile_sorted(boot_did, 1 - alpha / 2)
    
    // Sequential DID percentile CI
    result.ci_low[3] = quantile_sorted(boot_sdid, alpha / 2)
    result.ci_high[3] = quantile_sorted(boot_sdid, 1 - alpha / 2)
    
    return(result)
}

// ----------------------------------------------------------------------------
// MAIN ESTIMATION FUNCTION FOR STANDARD DID
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _did_std_main() - Main Estimation Orchestrator for Standard DID Design
 *
 * Coordinates the complete double DID estimation workflow for the standard
 * (non-staggered) difference-in-differences design.
 *
 * Workflow:
 *   1. Initialize result storage matrices
 *   2. Run cluster bootstrap for variance estimation
 *   3. For each lead value:
 *      a. Compute DID and sequential DID point estimates
 *      b. Compute GMM optimal weights from bootstrap VCOV
 *      c. Compute double DID point estimate and inference
 *   4. Store results in global variables for Stata retrieval
 *
 * Arguments:
 *   lead    : real rowvector - Lead values, e.g., (0) or (0, 1, 2)
 *   n_boot  : real scalar - Number of bootstrap iterations
 *   se_boot : real scalar - Inference method: 0 = asymptotic, 1 = bootstrap
 *   level   : real scalar - Confidence level in percent, e.g., 95
 *
 * Returns:
 *   real scalar - Return code: 0 = success, non-zero = error
 *     1 = Bootstrap failed (insufficient successful iterations)
 *     2 = VCOV computation failed (singular matrix)
 *
 * Side Effects:
 *   Populates global result variables: _did_b, _did_V, _did_estimates,
 *   _did_lead_values, _did_weights, _did_W, _did_vcov_gmm, _did_n_boot_success
 *---------------------------------------------------------------------------*/
real scalar _did_std_main(real rowvector lead, real scalar n_boot, 
                          real scalar se_boot, real scalar level)
{
    // External data structures populated by data preparation
    external struct did_data scalar did_dat
    external struct did_option scalar did_opt
    
    // External result variables for Stata retrieval
    external real rowvector _did_b
    external real matrix _did_V
    external real matrix _did_estimates
    external real rowvector _did_lead_values
    external real matrix _did_weights
    external real matrix _did_W
    external real matrix _did_vcov_gmm
    external real scalar _did_n_boot_success
    
    struct boot_result scalar boot_res
    struct gmm_weights scalar weights
    struct ddid_result scalar ddid_res
    
    real scalar n_lead, l, row, alpha, z
    real rowvector point_est
    real scalar se_did, se_sdid
    real matrix boot_est_l
    real scalar tau_did, tau_sdid, tau_ddid
    real scalar var_did, var_sdid, var_ddid
    real scalar ci_lo_did, ci_hi_did, ci_lo_sdid, ci_hi_sdid
    real scalar ci_lo_ddid, ci_hi_ddid
    real scalar w_did, w_sdid
    
    // -------------------------------------------------------------------------
    // Initialize result storage matrices
    // -------------------------------------------------------------------------
    n_lead = cols(lead)
    
    _did_b = J(1, 3 * n_lead, .)           // Coefficient vector: [dDID, DID, sDID] per lead
    _did_V = J(3 * n_lead, 3 * n_lead, 0)  // Variance-covariance matrix
    _did_estimates = J(3 * n_lead, 6, .)   // Full results table
    _did_lead_values = lead                 // Lead values for reference
    _did_weights = J(n_lead, 2, .)          // GMM weights (w_DID, w_sDID) per lead
    _did_W = J(n_lead, 4, .)                // Precision matrices (flattened) per lead
    _did_vcov_gmm = J(n_lead, 4, .)         // VCOV matrices (flattened) per lead
    _did_n_boot_success = n_boot            // Number of successful bootstrap iterations
    
    // Compute critical value for confidence intervals
    alpha = 1 - level / 100
    z = invnormal(1 - alpha / 2)
    
    // -------------------------------------------------------------------------
    // Run bootstrap for all lead values
    // -------------------------------------------------------------------------
    boot_res = did_boot_std(did_dat, lead, n_boot, did_opt.seed)
    _did_n_boot_success = boot_res.n_successful
    
    // Require minimum bootstrap successes for variance estimation
    if (boot_res.n_successful < 2) {
        errprintf("Error: Bootstrap failed - insufficient successful iterations\n")
        errprintf("       Only %g of %g iterations succeeded\n", 
                  boot_res.n_successful, n_boot)
        return(1)
    }
    
    // -------------------------------------------------------------------------
    // Compute estimates for each lead value
    // -------------------------------------------------------------------------
    row = 1
    for (l = 1; l <= n_lead; l++) {
        
        // ---------------------------------------------------------------------
        // Compute DID and sequential DID point estimates
        // ---------------------------------------------------------------------
        point_est = did_fit(
            did_dat.outcome,
            did_dat.outcome_delta,
            did_dat.Gi,
            did_dat.It,
            did_dat.covariates,
            did_dat.id_time_std,
            lead[l]
        )
        
        tau_did = point_est[1]
        tau_sdid = point_est[2]
        
        // ---------------------------------------------------------------------
        // Extract bootstrap VCOV and compute component standard errors
        // ---------------------------------------------------------------------
        if (boot_res.vcov[l] != NULL) {
            // Store flattened VCOV for Stata retrieval
            _did_vcov_gmm[l, .] = vec(*boot_res.vcov[l])'
            
            // Extract variances from VCOV diagonal
            var_did = (*boot_res.vcov[l])[1, 1]
            var_sdid = (*boot_res.vcov[l])[2, 2]
            
            // Compute component standard errors
            se_did = sqrt(var_did)
            se_sdid = sqrt(var_sdid)
            
            // -----------------------------------------------------------------
            // Compute GMM optimal weights from precision matrix
            // -----------------------------------------------------------------
            weights = compute_weights(*boot_res.vcov[l])
            w_did = weights.w_did
            w_sdid = weights.w_sdid
            
            // Store weights and precision matrix for Stata retrieval
            _did_weights[l, .] = (w_did, w_sdid)
            _did_W[l, .] = vec(weights.W)'
            
            // -----------------------------------------------------------------
            // Compute double DID point estimate and inference
            // -----------------------------------------------------------------
            boot_est_l = boot_res.estimates[., (2*l-1)..(2*l)]
            
            ddid_res = compute_double_did(tau_did, tau_sdid, weights, 
                                          boot_est_l, se_boot, level)
            
            tau_ddid = ddid_res.estimate
            var_ddid = ddid_res.variance
            ci_lo_ddid = ddid_res.ci_low
            ci_hi_ddid = ddid_res.ci_high
            
            // -----------------------------------------------------------------
            // Compute confidence intervals for component estimators
            // -----------------------------------------------------------------
            if (se_boot) {
                // Bootstrap percentile confidence intervals
                ci_lo_did = quantile_sorted(boot_est_l[., 1], alpha / 2)
                ci_hi_did = quantile_sorted(boot_est_l[., 1], 1 - alpha / 2)
                ci_lo_sdid = quantile_sorted(boot_est_l[., 2], alpha / 2)
                ci_hi_sdid = quantile_sorted(boot_est_l[., 2], 1 - alpha / 2)
            }
            else {
                // Asymptotic Wald confidence intervals
                ci_lo_did = tau_did - z * se_did
                ci_hi_did = tau_did + z * se_did
                ci_lo_sdid = tau_sdid - z * se_sdid
                ci_hi_sdid = tau_sdid + z * se_sdid
            }
        }
        else {
            // VCOV computation failed
            errprintf("{err}Error: Bootstrap VCOV computation failed for lead %g\n", lead[l])
            errprintf("{err}       VCOV pointer is NULL - no valid bootstrap covariance available\n")
            errprintf("{err}       This may be caused by:\n")
            errprintf("{err}       - Insufficient successful bootstrap iterations\n")
            errprintf("{err}       - All bootstrap estimates are missing for this lead\n")
            return(2)
        }
        
        // ---------------------------------------------------------------------
        // Populate Stata result matrices
        // ---------------------------------------------------------------------
        
        // Coefficient vector e(b): [dDID, DID, sDID] for each lead
        _did_b[1, 3*(l-1)+1] = tau_ddid
        _did_b[1, 3*(l-1)+2] = tau_did
        _did_b[1, 3*(l-1)+3] = tau_sdid
        
        // Variance matrix e(V): diagonal elements
        _did_V[3*(l-1)+1, 3*(l-1)+1] = var_ddid
        _did_V[3*(l-1)+2, 3*(l-1)+2] = var_did
        _did_V[3*(l-1)+3, 3*(l-1)+3] = var_sdid
        
        // Off-diagonal: DID-sDID covariance from bootstrap
        if (boot_res.vcov[l] != NULL) {
            _did_V[3*(l-1)+2, 3*(l-1)+3] = (*boot_res.vcov[l])[1, 2]
            _did_V[3*(l-1)+3, 3*(l-1)+2] = (*boot_res.vcov[l])[2, 1]
        }
        
        // Results table e(estimates): [lead, estimate, SE, CI_lo, CI_hi, weight]
        // Row order per lead: double DID, DID, sequential DID
        
        // Double DID result row
        _did_estimates[row, 1] = lead[l]
        _did_estimates[row, 2] = tau_ddid
        _did_estimates[row, 3] = sqrt(var_ddid)
        _did_estimates[row, 4] = ci_lo_ddid
        _did_estimates[row, 5] = ci_hi_ddid
        _did_estimates[row, 6] = .
        row = row + 1
        
        // DID result row
        _did_estimates[row, 1] = lead[l]
        _did_estimates[row, 2] = tau_did
        _did_estimates[row, 3] = se_did
        _did_estimates[row, 4] = ci_lo_did
        _did_estimates[row, 5] = ci_hi_did
        _did_estimates[row, 6] = w_did
        row = row + 1
        
        // Sequential DID result row
        _did_estimates[row, 1] = lead[l]
        _did_estimates[row, 2] = tau_sdid
        _did_estimates[row, 3] = se_sdid
        _did_estimates[row, 4] = ci_lo_sdid
        _did_estimates[row, 5] = ci_hi_sdid
        _did_estimates[row, 6] = w_sdid
        row = row + 1
    }
    
    return(0)
}

// ----------------------------------------------------------------------------
// MODULE VERIFICATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _did_gmm_loaded() - Module Load Verification
 *
 * Prints confirmation message when module is successfully loaded.
 *---------------------------------------------------------------------------*/
void _did_gmm_loaded()
{
    printf("{txt}did_gmm.mata loaded successfully\n")
}

end
