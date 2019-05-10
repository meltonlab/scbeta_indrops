import logging
from typing import *
import numpy as np
import scipy as sp
import pandas as pd

import loompy
import cytograph as cg

from sklearn.svm import SVR

from scbeta_scrnaseq import utils
import scbeta_scrnaseq.cytograph_inmem_utils as cgm

def vgam_regression(ds):

    from rpy2.robjects.packages import importr
    from rpy2.robjects import r, pandas2ri
    import rpy2.robjects as robjects
    pandas2ri.activate()

    ds.vals = sp.sparse.csr_matrix(ds.layers[""][:, :])

    # Regress theta (1/dispersion) for each gene
    totals_normalizer = cgm.CountsNormalizer(mean_center=False, 
                                    normalize_variance=False)
    totals_normalizer.fit(ds.vals)

    mu = totals_normalizer.mu
    sd = totals_normalizer.sd
    disp = (sd**2 - mu)/(mu**2)

    ok = (mu > 0) & (disp > 0) & (sd > 0)

    log2_mu = np.log2(mu[ok])
    log2_disp = np.log2(disp[ok])

    svr_gamma = 1000. / len(log2_mu)
    clf = SVR(gamma=svr_gamma)
    clf.fit(log2_mu[:, None], log2_disp)
    fitted_fun = clf.predict
    log2_disp_est = fitted_fun(log2_mu[:, None])

    # Adjust estimates to make them monotonic
    log2_mu_at_disp_min = log2_mu[log2_disp_est.argmin()]
    log2_mu_at_disp_max = log2_mu[log2_disp_est.argmax()]

    log2_disp_est[log2_mu > log2_mu_at_disp_min] = log2_disp_est.min()
    log2_disp_est[log2_mu < log2_mu_at_disp_max] = log2_disp_est.max()


    theta_est = np.ones_like(mu)*np.nan
    theta_est[ok] = 1/(2**log2_disp_est)
    ds.ra.NegBinomThetaEstimate = theta_est

    
    regress_genes = np.where((~np.isnan(ds.ra.NegBinomThetaEstimate)) & (ds.ra._Regress > 0))[0]
    valid_cells = np.where(ds.ca._Valid > 0)[0]
    
    # Check if we  model should include branches
    use_branch_model = 'PseudotimeBranch' in ds.ca.keys()

    # Prepare output layers
    n_pred_points = 100
    n_genes = len(ds.ra._Valid)

    if use_branch_model:
        import itertools
        _branch, _pseudotime = list(zip(*itertools.product([0, 1], np.linspace(0, 1, n_pred_points))))
        prediction_input_data = pd.DataFrame({
                'branch' : _branch,
                'pseudotime' : _pseudotime,
            })

        lrt_pseudotime_pvals = np.zeros((n_genes, ))*np.nan
        lrt_pseudotime_chisq = np.zeros((n_genes, ))*np.nan
        lrt_branch_pvals = np.zeros((n_genes, ))*np.nan
        lrt_branch_chisq = np.zeros((n_genes, ))*np.nan

        models_log_lik = np.zeros((n_genes, 3))*np.nan

        pred__pseudotime = np.zeros((n_genes, 2*n_pred_points))*np.nan
        pred__branch = np.zeros((n_genes, 2*n_pred_points))*np.nan

    else:
        prediction_input_data = pd.DataFrame({
                'pseudotime' : np.linspace(0, 1, n_pred_points)
            })

        lrt_pseudotime_pvals = np.zeros((n_genes, ))*np.nan
        lrt_pseudotime_chisq = np.zeros((n_genes, ))*np.nan

        models_log_lik = np.zeros((n_genes, 2))*np.nan

        pred__pseudotime = np.zeros((n_genes, n_pred_points))*np.nan


    # Begin regressions
    vgam = importr('VGAM')
    negbinomial_size = robjects.r('negbinomial.size')
    ns_df = 3
    gi_done = 0
    for gi in regress_genes:
        try:
            g_theta = ds.ra.NegBinomThetaEstimate[gi]

            g_data = pd.DataFrame({
                    'counts' : totals_normalizer.transform(ds.vals, genes=[gi], cells=valid_cells).reshape(-1).round(),
                    'pseudotime': ds.ca.PseudotimeRank[valid_cells], 
                })
            if use_branch_model:
                g_data['branch'] = ds.ca.PseudotimeBranch[valid_cells]

            # Fit the models
            null_model = vgam.vglm(robjects.Formula('counts ~ 1'), 
                  data = g_data,
                  family = negbinomial_size(g_theta),
                  epsilon = 1e-1)

            pseudotime_model = vgam.vglm(robjects.Formula(f'counts ~ sm.ns(pseudotime, df={ns_df})'), 
                  data = g_data,
                  family = negbinomial_size(g_theta),
                  epsilon = 1e-1)

            if use_branch_model:
                branch_model = vgam.vglm(robjects.Formula(f'counts ~ sm.ns(pseudotime, df={ns_df})*branch'), 
                  data = g_data,
                  family = negbinomial_size(g_theta),
                  epsilon = 1e-1)

            # Use log-likelihood ratio test to compute p-values
            pseudotime_lrt = vgam.lrtest_vglm(pseudotime_model, null_model)
            lrt_pseudotime_pvals[gi] = np.array(pseudotime_lrt.slots['Body'])[:, 1][4]
            lrt_pseudotime_chisq[gi] = np.array(pseudotime_lrt.slots['Body'])[:, 1][3]

            models_log_lik[gi, :2] = np.array(pseudotime_lrt.slots['Body'])[1, :].ravel()

            if use_branch_model:
                branch_lrt = vgam.lrtest_vglm(branch_model, pseudotime_model)
                lrt_branch_pvals[gi] = np.array(branch_lrt.slots['Body'])[:, 1][4]
                lrt_branch_chisq[gi] = np.array(branch_lrt.slots['Body'])[:, 1][3]
                
                models_log_lik[gi, :2] = np.array(branch_lrt.slots['Body'])[1, :].ravel()
                models_log_lik[gi, 1:3] = np.array(pseudotime_lrt.slots['Body'])[1, :].ravel()

            # Predict expression values on values spanning the input interval
            pred__pseudotime[gi] = np.array(vgam.predict(pseudotime_model, newdata=prediction_input_data)).reshape(-1)
            if use_branch_model:
                pred__branch[gi] = np.array(vgam.predict(branch_model, newdata=prediction_input_data)).reshape(-1)

            gi_done += 1
        except Exception as ex:
            print(ex)

        if (gi_done % 250) == 0:
            print(f'Regression completed for {gi_done} genes')

    ds.ra.LRT_Pseudotime_pval = lrt_pseudotime_pvals
    ds.ra.ModelLogLik = models_log_lik 
    ds.ra.Pred__pseudotime = pred__pseudotime
    if use_branch_model:
        ds.ra.LRT_Branch_pval = lrt_branch_pvals
        ds.ra.Pred__branch = pred__branch

def annotate_vgam_ds(pdt_ds, pv=0.1, min_pval=10**-300, fdr_alpha=0.001, n_smooth=5):
    branched = 'PseudotimeBranch' in pdt_ds.ca.keys()

    _regressed = np.where(~np.isnan(pdt_ds.ra.LRT_Pseudotime_pval))[0]

    #Adjust predicted expression using pseudovalue.
    pred_expr = np.exp(pdt_ds.ra.Pred__pseudotime)[_regressed]
    pred_expr[pred_expr < pv] = pv

    log2_pred_expr = np.log2(pred_expr)
    log2_start_end_fc = log2_pred_expr[:, -n_smooth:].mean(1) - log2_pred_expr[:, :n_smooth].mean(1)

    log2_min_max_fc = np.sign(log2_start_end_fc) * (log2_pred_expr.max(1) - log2_pred_expr.min(1))

    pseudotime_pval = pdt_ds.ra.LRT_Pseudotime_pval[_regressed] + min_pval

    from statsmodels.sandbox.stats.multicomp import multipletests
    _, pseudotime_qval, _, _, = multipletests(pseudotime_pval, fdr_alpha, method="fdr_bh")

    _reg = np.zeros_like(pdt_ds.ra._Valid)
    _reg[_regressed] = 1
    pdt_ds.ra['_Regressed'] = _reg

    _fc = np.zeros_like(pdt_ds.ra.LRT_Pseudotime_pval)*np.nan
    _fc[_regressed] = log2_start_end_fc
    pdt_ds.ra['Pred__pseudotime__log2fc__start_end'] = _fc

    _fc = np.zeros_like(pdt_ds.ra.LRT_Pseudotime_pval)*np.nan
    _fc[_regressed] = log2_min_max_fc
    pdt_ds.ra['Pred__pseudotime__log2fc__min_max'] = _fc

    _qval = np.zeros_like(pdt_ds.ra.LRT_Pseudotime_pval)*np.nan
    _qval[_regressed] = pseudotime_qval
    pdt_ds.ra['LRT_Pseudotime_qval'] = _qval

    _pred = np.zeros_like(pdt_ds.ra.Pred__pseudotime)*np.nan
    _pred[_regressed] = pred_expr
    pdt_ds.ra['Pred__pseudotime__tpm'] = _pred

    if branched:

        n_pred_points = 100
        pred_br = np.exp(pdt_ds.ra.Pred__branch[_regressed])
        pred_br[pred_br < pv] = pv

        log2_pred_br = np.log2(pred_br)
        br_fc = log2_pred_br[:,(n_pred_points-n_smooth-1):(n_pred_points-1)].mean(1) - log2_pred_br[:,-n_smooth:].mean(1)

        br_pval = pdt_ds.ra.LRT_Branch_pval[_regressed] + min_pval
        _, br_qval, _, _, = multipletests(br_pval, fdr_alpha, method="fdr_bh")

        _fc = np.zeros_like(pdt_ds.ra.LRT_Branch_pval)*np.nan
        _fc[_regressed] = br_fc
        pdt_ds.ra['Pred__branch__log2fc'] = _fc

        _qval = np.zeros_like(pdt_ds.ra.LRT_Branch_pval)*np.nan
        _qval[_regressed] = br_qval
        pdt_ds.ra['LRT_branch_qval'] = _qval

        _pred = np.zeros_like(pdt_ds.ra.Pred__branch)*np.nan
        _pred[_regressed] = pred_br
        pdt_ds.ra['Pred__branch__tpm'] = _pred


