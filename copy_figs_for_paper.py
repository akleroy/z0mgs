import os

in_dir = '../plots/'
out_dir = '../plots/for_survey_paper/'
astroph_dir = '../plots/for_astroph/'

# Define a dictionary with plots for the paper.

plots = {
#   plot_distance_unc.pro
    'fig1a':'edd_vs_leda.eps',
    'fig1b':'edd_vs_cf3.eps',
    'fig1c':'dist_hist.eps',
    'fig1d':'dist_sketch.eps',
#   motivate_selection.pro
    'fig2':'size_vs_dist.png',
    'fig3':'ston_vs_absb.png',
    'fig4a':'select_wedge.png',
    'fig4b':'select_cdf.png',
    'fig5a':'morph_vs_dist.png',
#   plot_detection_vs_absb
    'fig5b':'absb_vs_dist.png',
#   plot_distribution.pro
    'fig6a':'distrib_dist_xy.png',
    'fig6b':'distrib_dist_xz.png',
#   plot_unwise_noise.pro
    'fig7a':'unwise_noise_band1.eps',
    'fig7b':'unwise_noise_band2.eps',
    'fig7c':'unwise_noise_band3.eps',
    'fig7d':'unwise_noise_band4.eps',
#   plot_galex_noise.pro
    'fig7e':'galex_noise_NUV.eps',
    'fig7f':'galex_noise_FUV.eps',
#   plot_unwise_stats_vs_b.pro
    'fig8a':'unwise_noise_lat_band1.eps',
    'fig8b':'unwise_noise_lat_band2.eps',
    'fig8c':'unwise_noise_lat_band3.eps',
    'fig8d':'unwise_noise_lat_band4.eps',
#   plot_galex_stats_vs_b.pro
    'fig8e':'galex_noise_lat_NUV.eps',
    'fig8f':'galex_noise_lat_FUV.eps',
#   plot_galex_time.pro
    'fig9a':'galex_time.eps',
#   plot_galex_noise_scaling.pro
    'fig9b':'galex_noise_vs_t.eps',
#   plot_galex_extinction.pro
    'fig10':'galex_extinction.eps',
#   plot_nearby_3color.pro
    'fig11':'gallery_3color_1.png',
    'fig12':'gallery_3color_2.png',
    'fig13':'gallery_3color_3.png',
#   plot_nga_comp.pro (res_str='gauss7p5' and not)
    'fig14a':'nga_comp_fuv_gauss7p5.eps',
    'fig14b':'nga_comp_fuv_gauss15.eps',
    'fig14c':'nga_comp_nuv_gauss7p5.eps',
    'fig14d':'nga_comp_nuv_gauss15.eps',
#   plot_s4g_comp.pro (res_str='gauss7p5' and not)
    'fig14e':'s4g_comp_wise1_gauss7p5.eps',
    'fig14f':'s4g_comp_wise1_gauss15.eps',
    'fig14g':'s4g_comp_wise2_gauss7p5.eps',
    'fig14h':'s4g_comp_wise2_gauss15.eps',    
#   plot_intens_vs_mag.pro
    'fig15a':'mag_vs_intens_wise_gauss7p5.png',
    'fig15b':'mag_vs_intens_wise_gauss15.png',
    'fig15c':'mag_vs_intens_galex_gauss7p5.png',
    'fig15d':'mag_vs_intens_galex_gauss15.png',
#   plot_beam_from_stars.pro
    'fig16a':'beam_from_stars_gauss7p5.eps',
    'fig16b':'beam_from_stars_gauss15.eps',
#   plot_photometry_clarkdale.pro
#   plot_photometry_s4g.pro
    'fig17a':'z0mgs_phot_checkwise1.eps',
    'fig17b':'s4g_integrated_wise1.eps',
    'fig17c':'z0mgs_phot_checkwise2.eps',
    'fig17d':'s4g_integrated_wise2.eps',  
    'fig18a':'z0mgs_phot_checkwise3.eps',
    'fig18b':'z0mgs_phot_checkwise4.eps',
    'fig18c':'z0mgs_phot_checknuv.eps',
    'fig18d':'z0mgs_phot_checkfuv.eps',
#   plot_mainseq.pro
    'fig19a':'z0mgs_mainseq_contours.eps',
    'fig19b':'z0mgs+gswlc_mainseq_contours.eps',
    'fig19c':'z0mgs_mainseq.png',
    'fig19d':'markup_mainseq.eps',
#   plot_samples.pro
    'fig20a':'nuvw1_vs_w3w1.eps',
    'fig20b':'band_vs_band.eps',
#   Appendix starts here
#   build_gswlc_grid.pro, plot_gswlc_grids.pro
    'fig21a':'mtol_grid_ssfrmstar.eps',
    'fig21b':'cfuvw4_grid_ssfrmstar.eps',
    'fig21c':'cfuvw3_grid_ssfrmstar.eps',
    'fig21d':'cjustw4_grid_ssfrmstar.eps',
#   plot_gswlc_mtol.pro
    'fig22a':'gswlc_mtol_hist.eps',
    'fig22b':'gswlc_mtol_ssfr.eps',
    'fig23a':'gswlc_mtol_wise1.eps',
    'fig23b':'gswlc_mtol_mstar.eps',
    'fig23c':'gswlc_mtol_w2w1.eps',
    'fig23d':'gswlc_mtol_nuvw1.eps',
    'fig23e':'gswlc_mtol_w4w1.eps',
    'fig23f':'gswlc_mtol_ssfrlike.eps',
#   plot_gswlc_mtolresids.pro
    'fig24a':'mtol_fix_estimate.eps',
    'fig24b':'mtol_w4w1_estimate.eps',
    'fig24c':'mtol_ssfrlikefuvw4_estimate.eps',
    'fig24d':'mtol_ssfr_estimate.eps',
#   plot_gswlc_sfruv.pro
    'fig25a':'gswlc_cfuv.eps',
    'fig25b':'cfuv_grid_ssfrmstar.eps',
#   plot_gswlc_sfr_calibs.pro
    'fig26':'gswlc_cwise.eps',
#   plot_gswlc_sfr.pro
    'fig27a':'gswlc_cfuvw4_ssfr.eps',
    'fig27b':'gswlc_cfuvw4_nuvw1.eps',
    }

# Copy those plots to the output directory for the paper

for key in plots.keys():
    this_fig_name = plots[key]
    line = 'cp '
    line += in_dir+this_fig_name+' '
    line += out_dir+this_fig_name+' '
    os.system(line)
