  # print('running will old long code.')
        # baseline_counts_vs_cgc_status_data = preliminary_data[["BASELINE", ref_gene_column]].copy()
        # our_method_counts_vs_cgc_status_data = preliminary_data[["OUR_METHOD", ref_gene_column]].copy()
        # elaspic_cov_vs_cgc_status_data = preliminary_data[["ELASPIC_COVERAGE", ref_gene_column]].copy()
        #
        # baseline_counts_vs_cgc_status_tcga_data = preliminary_data[["BASELINE", ref_gene_column_cohort]].copy()
        # our_method_counts_vs_cgc_status_tcga_data = preliminary_data[["OUR_METHOD", ref_gene_column_cohort]].copy()
        # elaspic_cov_vs_cgc_status_tcga_data = preliminary_data[["ELASPIC_COVERAGE", ref_gene_column_cohort]].copy()
        #
        # fpr_baseline, tpr_baseline, ths_baseline = get_fpr_tpr_ths(
        #     np.array(baseline_counts_vs_cgc_status_data["CancerMine_STATUS"]),  # FIXME: ref_gene_column
        #     np.array(baseline_counts_vs_cgc_status_data["BASELINE"]))
        #
        # fpr_our_method, tpr_our_method, ths_our_method = get_fpr_tpr_ths(
        #     np.array(our_method_counts_vs_cgc_status_data["CancerMine_STATUS"]),
        #     np.array(our_method_counts_vs_cgc_status_data["OUR_METHOD"]))
        #
        # fpr_elaspic_cov, tpr_elaspic_cov, ths_elaspic_cov = get_fpr_tpr_ths(
        #     np.array(elaspic_cov_vs_cgc_status_data["CancerMine_STATUS"]),
        #     np.array(elaspic_cov_vs_cgc_status_data["ELASPIC_COVERAGE"]))
        #
        # fpr_baseline_brca, tpr_baseline_brca, ths_baseline_brca = get_fpr_tpr_ths(
        #     np.array(baseline_counts_vs_cgc_status_tcga_data["CancerMine_STATUS (BRCA)"]),
        #     np.array(baseline_counts_vs_cgc_status_tcga_data["BASELINE"]))
        #
        # fpr_our_method_brca, tpr_our_method_brca, ths_our_method_brca = get_fpr_tpr_ths(
        #     np.array(our_method_counts_vs_cgc_status_tcga_data["CancerMine_STATUS (BRCA)"]),
        #     np.array(our_method_counts_vs_cgc_status_tcga_data["OUR_METHOD"]))
        #
        # fpr_elaspic_cov_brca, tpr_elaspic_cov_brca, ths_elaspic_cov_brca = get_fpr_tpr_ths(
        #     np.array(elaspic_cov_vs_cgc_status_tcga_data["CancerMine_STATUS (BRCA)"]),
        #     np.array(elaspic_cov_vs_cgc_status_tcga_data["ELASPIC_COVERAGE"]))
        #
        # roc_auc_baseline = metrics.auc(fpr_baseline, tpr_baseline)
        # roc_auc_our_method = metrics.auc(fpr_our_method, tpr_our_method)
        # roc_auc_elaspic_cov = metrics.auc(fpr_elaspic_cov, tpr_elaspic_cov)
        #
        # roc_auc_baseline_tcga = metrics.auc(fpr_baseline_brca, tpr_baseline_brca)
        # roc_auc_our_method_tcga = metrics.auc(fpr_our_method_brca, tpr_our_method_brca)
        # roc_auc_elaspic_cov_tcga = metrics.auc(fpr_elaspic_cov_brca, tpr_elaspic_cov_brca)
        #
        # # Plotting CancerMine_STATUS
        # plt.figure(figsize=(7, 5))
        # sns.set(style='white', font_scale=1.50)
        #
        # plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ STATUS$ vs Various Columns')
        # plt.plot(fpr_baseline, tpr_baseline, label='baseline (%0.3f)' % roc_auc_baseline)
        # plt.plot(fpr_our_method, tpr_our_method, label='our_method (%0.3f)' % roc_auc_our_method)
        # plt.plot(fpr_elaspic_cov, tpr_elaspic_cov, label='elaspic_cov (%0.3f)' % roc_auc_elaspic_cov)
        #
        # plt.legend(loc='lower right')
        # plt.plot([0, 1], [0, 1], 'k--')
        # plt.xlim([0, 1])
        # plt.ylim([0, 1])
        # plt.ylabel('True Positive Rate\n(Sensitivity)')
        # plt.xlabel('False Positive Rate\n(1-Specificity)')
        # plt.show()
        #
        # # Plotting CancerMine_STATUS_brca
        # plt.figure(figsize=(7, 5))
        # sns.set(style='white', font_scale=1.50)
        #
        # plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ BRCA\ STATUS$ vs Various Columns')
        # plt.plot(fpr_baseline_brca, tpr_baseline_brca, label='baseline BRCA  (%0.3f)' % roc_auc_baseline_tcga)
        # plt.plot(fpr_our_method_brca, tpr_our_method_brca, label='our_method BRCA  (%0.3f)' % roc_auc_our_method_tcga)
        # plt.plot(fpr_elaspic_cov_brca, tpr_elaspic_cov_brca,
        #          label='elaspic_cov BRCA  (%0.3f)' % roc_auc_elaspic_cov_tcga)
        #
        # plt.legend(loc='lower right')
        # plt.plot([0, 1], [0, 1], 'k--')
        # plt.xlim([0, 1])
        # plt.ylim([0, 1])
        # plt.ylabel('True Positive Rate\n(Sensitivity)')
        # plt.xlabel('False Positive Rate\n(1-Specificity)')
        # plt.show()
        #
        # # # ROC with `BASELINE` Non-zero Columns # # # # # # # # # # # # # # #
        # preliminary_data_baseline_nonzero = preliminary_data[preliminary_data['BASELINE'] != 0].copy()
        # print(preliminary_data_baseline_nonzero.shape)
        #
        # baseline_counts_vs_cgc_status_bnz_data = preliminary_data_baseline_nonzero[["BASELINE", ref_gene_column]].copy()
        # our_method_counts_vs_cgc_status_bnz_data = preliminary_data_baseline_nonzero[["OUR_METHOD", ref_gene_column]].copy()
        # elaspic_cov_vs_cgc_status_bnz_data = preliminary_data_baseline_nonzero[["ELASPIC_COVERAGE", ref_gene_column]].copy()
        #
        # baseline_counts_vs_cgc_status_brca_bnz_data = preliminary_data_baseline_nonzero[
        #     ["BASELINE", ref_gene_column_cohort]].copy()
        # our_method_counts_vs_cgc_status_brca_bnz_data = preliminary_data_baseline_nonzero[
        #     ["OUR_METHOD", ref_gene_column_cohort]].copy()
        # elaspic_cov_vs_cgc_status_brca_bnz_data = preliminary_data_baseline_nonzero[
        #     ["ELASPIC_COVERAGE", ref_gene_column_cohort]].copy()
        #
        # fpr_baseline_bnz, tpr_baseline_bnz, ths_baseline_bnz = get_fpr_tpr_ths(
        #     np.array(baseline_counts_vs_cgc_status_bnz_data["CancerMine_STATUS"]),
        #     np.array(baseline_counts_vs_cgc_status_bnz_data["BASELINE"]))
        #
        # fpr_our_method_bnz, tpr_our_method_bnz, ths_our_method_bnz = get_fpr_tpr_ths(
        #     np.array(our_method_counts_vs_cgc_status_bnz_data["CancerMine_STATUS"]),
        #     np.array(our_method_counts_vs_cgc_status_bnz_data["OUR_METHOD"]))
        #
        # fpr_elaspic_cov_bnz, tpr_elaspic_cov_bnz, ths_elaspic_cov_bnz = get_fpr_tpr_ths(
        #     np.array(elaspic_cov_vs_cgc_status_bnz_data["CancerMine_STATUS"]),
        #     np.array(elaspic_cov_vs_cgc_status_bnz_data["ELASPIC_COVERAGE"]))
        #
        # fpr_baseline_bnz_brca, tpr_baseline_bnz_brca, ths_baseline_bnz_brca = get_fpr_tpr_ths(
        #     np.array(baseline_counts_vs_cgc_status_brca_bnz_data["CancerMine_STATUS (BRCA)"]),
        #     np.array(baseline_counts_vs_cgc_status_brca_bnz_data["BASELINE"]))
        #
        # fpr_our_method_bnz_brca, tpr_our_method_bnz_brca, ths_our_method_bnz_brca = get_fpr_tpr_ths(
        #     np.array(our_method_counts_vs_cgc_status_brca_bnz_data["CancerMine_STATUS (BRCA)"]),
        #     np.array(our_method_counts_vs_cgc_status_brca_bnz_data["OUR_METHOD"]))
        #
        # fpr_elaspic_cov_bnz_brca, tpr_elaspic_cov_bnz_brca, ths_elaspic_cov_bnz_brca = get_fpr_tpr_ths(
        #     np.array(elaspic_cov_vs_cgc_status_brca_bnz_data["CancerMine_STATUS (BRCA)"]),
        #     np.array(elaspic_cov_vs_cgc_status_brca_bnz_data["ELASPIC_COVERAGE"]))
        #
        # roc_auc_baseline_bnz = metrics.auc(fpr_baseline_bnz, tpr_baseline_bnz)
        # roc_auc_our_method_bnz = metrics.auc(fpr_our_method_bnz, tpr_our_method_bnz)
        # roc_auc_elaspic_cov_bnz = metrics.auc(fpr_elaspic_cov_bnz, tpr_elaspic_cov_bnz)
        #
        # roc_auc_baseline_bnz_brca = metrics.auc(fpr_baseline_bnz_brca, tpr_baseline_bnz_brca)
        # roc_auc_our_method_bnz_brca = metrics.auc(fpr_our_method_bnz_brca, tpr_our_method_bnz_brca)
        # roc_auc_elaspic_cov_bnz_brca = metrics.auc(fpr_elaspic_cov_bnz_brca, tpr_elaspic_cov_bnz_brca)
        #
        # # Plotting BNZ CancerMine
        # plt.figure(figsize=(7, 5))
        # sns.set(style='white', font_scale=1.50)
        #
        # plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ STATUS$ vs Various Columns')
        # plt.plot(fpr_baseline_bnz, tpr_baseline_bnz, label='baseline bnz (%0.3f)' % roc_auc_baseline_bnz)
        # plt.plot(fpr_our_method_bnz, tpr_our_method_bnz, label='our_method bnz (%0.3f)' % roc_auc_our_method_bnz)
        # plt.plot(fpr_elaspic_cov_bnz, tpr_elaspic_cov_bnz, label='elaspic_cov bnz (%0.3f)' % roc_auc_elaspic_cov_bnz)
        #
        # plt.legend(loc='lower right')
        # plt.plot([0, 1], [0, 1], 'k--')
        # plt.xlim([0, 1])
        # plt.ylim([0, 1])
        # plt.ylabel('True Positive Rate\n(Sensitivity)')
        # plt.xlabel('False Positive Rate\n(1-Specificity)')
        # plt.show()
        #
        # # Plotting BNZ CancerMine_STATUS_brca
        # plt.figure(figsize=(7, 5))
        # sns.set(style='white', font_scale=1.25)
        #
        # plt.title('Receiver Operating Characteristic (ROC)\n$CancerMine\ BRCA\ STATUS$ vs Various Columns')
        # plt.plot(fpr_baseline_bnz_brca, tpr_baseline_bnz_brca, label='baseline (%0.3f)' % roc_auc_baseline_bnz_brca)
        # plt.plot(fpr_our_method_bnz_brca, tpr_our_method_bnz_brca,
        #          label='disruptive interactions (%0.3f)' % roc_auc_our_method_bnz_brca)
        # plt.plot(fpr_elaspic_cov_bnz_brca, tpr_elaspic_cov_bnz_brca,
        #          label='coverage (%0.3f)' % roc_auc_elaspic_cov_bnz_brca)
        #
        # plt.legend(loc='lower right')
        # plt.plot([0, 1], [0, 1], 'k--')
        # plt.xlim([0, 1])
        # plt.ylim([0, 1])
        # plt.ylabel('$True Positive Rate\n(Sensitivity)$')
        # plt.xlabel('False Positive Rate\n(1-Specificity)')
        # plt.show()
        #
        # # displayers
        # display(pd.DataFrame({
        #     "roc_auc_baseline": [roc_auc_baseline],
        #     "roc_auc_our_method": [roc_auc_our_method],
        #     "roc_auc_elaspic_cov": [roc_auc_elaspic_cov],
        # }, index=['AUC']).T)
        #
        # display(pd.DataFrame({
        #     "roc_auc_baseline_brca": [roc_auc_baseline_tcga],
        #     "roc_auc_our_method_brca": [roc_auc_our_method_tcga],
        #     "roc_auc_elaspic_cov_brca": [roc_auc_elaspic_cov_tcga],
        # }, index=['AUC']).T)
        #
        # display(pd.DataFrame({
        #     "roc_auc_baseline_bnz": [roc_auc_baseline_bnz],
        #     "roc_auc_our_method_bnz": [roc_auc_our_method_bnz],
        #     "roc_auc_elaspic_cov_bnz": [roc_auc_elaspic_cov_bnz],
        # }, index=['AUC']).T)
        #
        # display(pd.DataFrame({
        #     "roc_auc_baseline_bnz_brca": [roc_auc_baseline_bnz_brca],
        #     "roc_auc_our_method_bnz_brca": [roc_auc_our_method_bnz_brca],
        #     "roc_auc_elaspic_cov_bnz_brca": [roc_auc_elaspic_cov_bnz_brca],
        # }, index=['AUC']).T)



# predator_analysis = PredatorAnalysis(
#     tcga_code="brca",
#     snv_path=SNV_BRCA_PATH,
#     prediction_data_path=PREDICTION_BRCA_REDUCED_PATH,
#     elaspic_core_path=BRCA_CORE_PATH,
#     elaspic_interface_path=BRCA_INTERFACE_PATH,
#     reference_data="cancermine",
#     reference_data_path=CANCER_MINE_BREAST_PATH,
# )