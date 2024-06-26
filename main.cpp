
#include "sofa/sofa.h"
#include "sofa/sofam.h"
#include <fstream>
#include "Catalog.h"
#include "random"
#include "map"
#include "GaiaTest.h"

extern std::vector<int> query_disc(int Nside, double ra, double dec, double range, bool r);

auto queryThreshold(std::string gaia_dir,double ra, double dec, double raidus, const std::string outputpath = "/dev/shm/", const std::string outname = "results.csv") {
    ZQQ_AT_CHINA_VO::GaiaTest c("GAIAdr3", "no", 20000000,gaia_dir);
    auto r = c.QueryThreshold(ra, dec, raidus, 20000000);
    std::string content = "";
    content.append(c.header_line);
    content.append("\n");
    for (auto &l: r) {
        content.append(l);
        content.append("\n");
    }
    std::remove((outputpath + outname).c_str());
    c.write_to_file(outputpath , outname, content);
}



auto matchCatalog(std::string csvfile, std::string raname = "ra", std::string decname = "dec",
                  const std::string outputpath = "/dev/shm/", const std::string outname = "results.csv") {

    auto start = clock();
    ZQQ_AT_CHINA_VO::CsvRead c(csvfile, "dynamic", 20000000);
    c.HaveProperMotion = false;
    c.setRaDecNameInHeader(raname, decname);
    c.inputSingleCsv(csvfile);
    ZQQ_AT_CHINA_VO::Log() << "catalog read" << c.getCatlogSize();

    ZQQ_AT_CHINA_VO::GaiaTest gaia("GAIAdr3", "no", 20000000);
    auto res = gaia.MatchCatalogClosest(c);
    auto header = "dis," + c.header_line + "," + gaia.header_line;
    std::string content = "";
    content.append(header);
    content.append("\n");
    for (auto &l: res) {
        content.append(l);
        content.append("\n");
    }
    std::remove((outputpath + outname).c_str());
    c.write_to_file(outputpath, outname, content);
}

auto split(std::string inputDir, std::string outputDir, size_t catalogSize, size_t memoryLimit) {

    auto start = clock();
    ZQQ_AT_CHINA_VO::CsvRead c("GAIAdr3", "dynamic", memoryLimit);//20000000 @ 32G memory
    c.HaveProperMotion = false;
    c.setRaDecNameInHeader("ra", "dec");
    c.setThreshold(catalogSize);
    c.temp_dir = outputDir;
    c.setSplitDir(c.temp_dir);
    c.inputdir(outputDir);
    c.split();
    int InputTotRead = c.getCatlogSize();
    ZQQ_AT_CHINA_VO::Log() << "Catalog read" << InputTotRead << " " << c.CatalogTotRead;
    auto endclock = clock();
    double endtime = (double) (endclock - start) / CLOCKS_PER_SEC;
    ZQQ_AT_CHINA_VO::Log() << "Total time: " << endtime << "s";
    ZQQ_AT_CHINA_VO::Log() << "Merge  read" << c.CatalogTotRead;
    ZQQ_AT_CHINA_VO::Log() << "tot split" << c.split_tot_write;
    return c;
}

auto merge(std::string inputDir, std::string outputDir, size_t catalogSize, size_t memoryLimit) {
    auto start = clock();
    ZQQ_AT_CHINA_VO::CsvRead c("GAIAdr3", "dynamic", memoryLimit);
    c.HaveProperMotion = false;
    c.setRaLocDefault();
    c.header_line = "solution_id,designation,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,parallax_over_error,pm,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_bad_obs_al,astrometric_gof_al,astrometric_chi2_al,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_params_solved,astrometric_primary_flag,nu_eff_used_in_astrometry,pseudocolour,pseudocolour_error,ra_pseudocolour_corr,dec_pseudocolour_corr,parallax_pseudocolour_corr,pmra_pseudocolour_corr,pmdec_pseudocolour_corr,astrometric_matched_transits,visibility_periods_used,astrometric_sigma5d_max,matched_transits,new_matched_transits,matched_transits_removed,ipd_gof_harmonic_amplitude,ipd_gof_harmonic_phase,ipd_frac_multi_peak,ipd_frac_odd_win,ruwe,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,duplicated_source,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_flux_over_error,phot_g_mean_mag,phot_bp_n_obs,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_bp_mean_flux_over_error,phot_bp_mean_mag,phot_rp_n_obs,phot_rp_mean_flux,phot_rp_mean_flux_error,phot_rp_mean_flux_over_error,phot_rp_mean_mag,phot_bp_rp_excess_factor,phot_bp_n_contaminated_transits,phot_bp_n_blended_transits,phot_rp_n_contaminated_transits,phot_rp_n_blended_transits,phot_proc_mode,bp_rp,bp_g,g_rp,radial_velocity,radial_velocity_error,rv_method_used,rv_nb_transits,rv_nb_deblended_transits,rv_visibility_periods_used,rv_expected_sig_to_noise,rv_renormalised_gof,rv_chisq_pvalue,rv_time_duration,rv_amplitude_robust,rv_template_teff,rv_template_logg,rv_template_fe_h,rv_atm_param_origin,vbroad,vbroad_error,vbroad_nb_transits,grvs_mag,grvs_mag_error,grvs_mag_nb_transits,rvs_spec_sig_to_noise,phot_variable_flag,l,b,ecl_lon,ecl_lat,in_qso_candidates,in_galaxy_candidates,non_single_star,has_xp_continuous,has_xp_sampled,has_rvs,has_epoch_photometry,has_epoch_rv,has_mcmc_gspphot,has_mcmc_msc,in_andromeda_survey,classprob_dsc_combmod_quasar,classprob_dsc_combmod_galaxy,classprob_dsc_combmod_star,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,mh_gspphot,mh_gspphot_lower,mh_gspphot_upper,distance_gspphot,distance_gspphot_lower,distance_gspphot_upper,azero_gspphot,azero_gspphot_lower,azero_gspphot_upper,ag_gspphot,ag_gspphot_lower,ag_gspphot_upper,ebpminrp_gspphot,ebpminrp_gspphot_lower,ebpminrp_gspphot_upper,libname_gspphot";

    c.setThreshold(3000);
    c.final_split_dir = outputDir;
    c.temp_dir = inputDir;
    c.MergedDir = outputDir;
    c.setSplitDir(c.temp_dir);
    c.mergeCatalog("");
    int InputTotRead = c.getCatlogSize();
    ZQQ_AT_CHINA_VO::Log() << "catalog read" << InputTotRead << " " << c.CatalogTotRead;
    auto endclock = clock();
    double endtime = (double) (endclock - start) / CLOCKS_PER_SEC;
    ZQQ_AT_CHINA_VO::Log() << "Total time: " << endtime << "s";
    ZQQ_AT_CHINA_VO::Log() << "Merge  read" << c.CatalogTotRead;
    ZQQ_AT_CHINA_VO::Log() << "tot split" << c.split_tot_write;
    return c;
}

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Error: No command provided. See manual" << std::endl;
        std::cerr << "Usage: GaiaCatalogQuery [command] [options]" << std::endl;
        return 1;
    }

    std::string command = argv[1];

    if (command == "split") {
        if (argc < 4) {
            std::cerr << "Error: Not enough arguments for split-temp command. See manual" << std::endl;
            std::cerr << "Usage: GaiaCatalogQuery split <inputDir> <outputDir> [catalogSize] [memoryLimit]"
                      << std::endl;
            std::cerr
                    << "Example: GaiaCatalogQuery split /home/user/Gaia_dr3/gaia_source_uncompressed/ /home/user/gaia_temp/ 1000 20000000"
                    << std::endl;

            return 1;
        }
        std::string inputDir = argv[2];
        std::string outputDir = argv[3];
        size_t catalogSize = (argc > 4) ? std::stoul(argv[4]) : 1000;
        size_t memoryLimit = (argc > 5) ? std::stoul(argv[5]) : 20000000;
        split(inputDir, outputDir, catalogSize, memoryLimit);
    } else if (command == "merge") {
        if (argc < 4) {
            std::cerr << "Error: Not enough arguments for split-temp command. See manual" << std::endl;
            std::cerr << "Usage: GaiaCatalogQuery merge <inputDir> <outputDir> [catalogSize] [memoryLimit]"
                      << std::endl;
            std::cerr << "Example: GaiaCatalogQuery merge /home/user/gaia_temp/ /home/user/gaia_final/ 1000 20000000"
                      << std::endl;

            return 1;
        }
        std::string inputDir = argv[2];
        std::string outputDir = argv[3];
        size_t catalogSize = (argc > 4) ? std::stoul(argv[4]) : 1000;
        size_t memoryLimit = (argc > 5) ? std::stoul(argv[5]) : 20000000;
        merge(inputDir, outputDir, catalogSize, memoryLimit);
    } else if (command == "search") {
        if (argc < 5) {
            std::cerr << "Error: Not enough arguments for search command." << std::endl;
            std::cerr << "Usage: program search <catalog_folder> <ra> <dec> <radius>" << std::endl;
            return 1;
        }
        std::string gaia_dir= argv[2];
        double ra = std::stod(argv[3]);
        double dec = std::stod(argv[4]);
        double radius = std::stod(argv[5]);
        queryThreshold(gaia_dir,ra, dec, radius);

    } else if (command == "cross-match") {
        std::string inputCatalog= argv[2];

        //matchCatalog();
    } else {
        std::cerr << "Error: Unknown command '" << command << "'." << std::endl;
        std::cerr << "Usage: program [split-unordered|split-temporary|search|cross-identify]" << std::endl;
        return 1;
    }

    return 0;

}
