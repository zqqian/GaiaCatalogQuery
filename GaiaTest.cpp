//
// Created by zqqia on 2022/08/15.
//

#include "GaiaTest.h"
//#include "match.pb.h"
//#include "healpix/healpix.h"
#include <thread>

#include <mutex>
#include <queue>

#include <condition_variable>
//void ZQQ_AT_CHINA_VO::GaiaTest::inputPixList(const std::string &fileName) {
//#if defined(_WIN32) || defined(_WIN64)
//
//    if ((_access(fileName.c_str(), 0)) == -1) {
//        //cout << "Unable to access " << name << endl;
//        return;
//    }
//#endif
//
//#if  defined(__GNUC__) || defined(__linux__)
//    if ((access(fileName.c_str(), 0)) == -1) {
//        std::cout << "Unable to access " << fileName << std::endl;
//        return;
//    }
//#endif
//    CrossMatch::PixList p;
//    std::fstream input(fileName, std::ios::in | std::ios::binary);
//    p.ParseFromIstream(&input);
//    Log() << fileName << p.pix_size();
//
//
////int max=0;
//    for (auto &pix: p.pix()) {
//        int pix2 = (pix >> 6);
//        // v[pix2]++;
//        //l->setId(getCatlogSize());
//        treeAdd(pix2, 1);
//        // Log() << uniq;
////        std::shared_ptr<Catalog::line> l = std::make_shared<Catalog::line>();
////        l->setPix(pix2);
////        addLine(l);
//    }
//}

int ZQQ_AT_CHINA_VO::GaiaTest::wirteCatalog(int nowlevel, int ipix, std::string Dir) {
//Log()<<"nowlevel "<<nowlevel<<" ipix "<<ipix;
    if (nowlevel < currentTreeLevel) {
        wirteCatalog(nowlevel + 1, ipix << 2, Dir);
        wirteCatalog(nowlevel + 1, (ipix << 2) + 1, Dir);
        wirteCatalog(nowlevel + 1, (ipix << 2) + 2, Dir);
        wirteCatalog(nowlevel + 1, (ipix << 2) + 3, Dir);

    } else {

        treeAdd(ipix, -1 * Tree[calUNIQ(currentTreeLevel, ipix, currentTreeLevel)]);

    }
    return 0;
}

ZQQ_AT_CHINA_VO::GaiaTest::GaiaTest(std::string name, std::string strindexType, long long int MEMORY_SIZE_LIMIT)
        : CsvRead(name, strindexType, MEMORY_SIZE_LIMIT) {
    header_line = "solution_id,designation,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,parallax_over_error,pm,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_bad_obs_al,astrometric_gof_al,astrometric_chi2_al,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_params_solved,astrometric_primary_flag,nu_eff_used_in_astrometry,pseudocolour,pseudocolour_error,ra_pseudocolour_corr,dec_pseudocolour_corr,parallax_pseudocolour_corr,pmra_pseudocolour_corr,pmdec_pseudocolour_corr,astrometric_matched_transits,visibility_periods_used,astrometric_sigma5d_max,matched_transits,new_matched_transits,matched_transits_removed,ipd_gof_harmonic_amplitude,ipd_gof_harmonic_phase,ipd_frac_multi_peak,ipd_frac_odd_win,ruwe,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,duplicated_source,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_flux_over_error,phot_g_mean_mag,phot_bp_n_obs,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_bp_mean_flux_over_error,phot_bp_mean_mag,phot_rp_n_obs,phot_rp_mean_flux,phot_rp_mean_flux_error,phot_rp_mean_flux_over_error,phot_rp_mean_mag,phot_bp_rp_excess_factor,phot_bp_n_contaminated_transits,phot_bp_n_blended_transits,phot_rp_n_contaminated_transits,phot_rp_n_blended_transits,phot_proc_mode,bp_rp,bp_g,g_rp,radial_velocity,radial_velocity_error,rv_method_used,rv_nb_transits,rv_nb_deblended_transits,rv_visibility_periods_used,rv_expected_sig_to_noise,rv_renormalised_gof,rv_chisq_pvalue,rv_time_duration,rv_amplitude_robust,rv_template_teff,rv_template_logg,rv_template_fe_h,rv_atm_param_origin,vbroad,vbroad_error,vbroad_nb_transits,grvs_mag,grvs_mag_error,grvs_mag_nb_transits,rvs_spec_sig_to_noise,phot_variable_flag,l,b,ecl_lon,ecl_lat,in_qso_candidates,in_galaxy_candidates,non_single_star,has_xp_continuous,has_xp_sampled,has_rvs,has_epoch_photometry,has_epoch_rv,has_mcmc_gspphot,has_mcmc_msc,in_andromeda_survey,classprob_dsc_combmod_quasar,classprob_dsc_combmod_galaxy,classprob_dsc_combmod_star,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,mh_gspphot,mh_gspphot_lower,mh_gspphot_upper,distance_gspphot,distance_gspphot_lower,distance_gspphot_upper,azero_gspphot,azero_gspphot_lower,azero_gspphot_upper,ag_gspphot,ag_gspphot_lower,ag_gspphot_upper,ebpminrp_gspphot,ebpminrp_gspphot_lower,ebpminrp_gspphot_upper,libname_gspphot";
}

//void ZQQ_AT_CHINA_VO::GaiaTest::writeTree() {
//    CrossMatch::PixList p;
//    for (auto &node: Tree) {
//        p.add_pix(node);
//    }
//    std::fstream output(final_split_dir + std::to_string(GaiaPbfTimes) + "Tree.pbf",
//                        std::ios::out | std::ios::trunc | std::ios::binary);
//    if (!p.SerializeToOstream(&output)) {
//        Log() << Log::ERR << "Failed to write";
//    }
//    google::protobuf::ShutdownProtobufLibrary();
//}

//void ZQQ_AT_CHINA_VO::GaiaTest::inputTree(const std::string &fileName) {
//#if defined(_WIN32) || defined(_WIN64)
//
//    if ((_access(fileName.c_str(), 0)) == -1) {
//        //cout << "Unable to access " << name << endl;
//        return;
//    }
//#endif
//
//#if  defined(__GNUC__) || defined(__linux__)
//    if ((access(fileName.c_str(), 0)) == -1) {
//        std::cout << "Unable to access " << fileName << std::endl;
//        return;
//    }
//#endif
//    CrossMatch::PixList p;
//    std::fstream input(fileName, std::ios::in | std::ios::binary);
//    p.ParseFromIstream(&input);
//    Log() << fileName << p.pix_size();
//
//
////int max=0;
//    int now = 0;
//    for (auto &pix: p.pix()) {
//        Tree[now++] += pix;
//        // v[pix2]++;
//        //l->setId(getCatlogSize());
//
//        // Log() << uniq;
////        std::shared_ptr<Catalog::line> l = std::make_shared<Catalog::line>();
////        l->setPix(pix2);
////        addLine(l);
//    }
//}

//int ZQQ_AT_CHINA_VO::GaiaTest::calPixTot(const std::string &fileName, int uniq, int level) {
//#if defined(_WIN32) || defined(_WIN64)
//
//    if ((_access(fileName.c_str(), 0)) == -1) {
//        //cout << "Unable to access " << name << endl;
//        return 0;
//    }
//#endif
//
//#if  defined(__GNUC__) || defined(__linux__)
//    if ((access(fileName.c_str(), 0)) == -1) {
//        std::cout << "Unable to access " << fileName << std::endl;
//        return 0;
//    }
//#endif
//    CrossMatch::Catalog p;
//    std::fstream input(fileName, std::ios::in | std::ios::binary);
//    p.ParseFromIstream(&input);
//    int tot = 0;
//    for (auto &i: p.lines()) {
//        if (eq2pix_nest(NSIDE, i.ra(), i.dec()) == i.pix() && calUNIQ(level, i.pix(), nside2Level(NSIDE)) == uniq) {
//            tot++;
//            HASH *= (((long long) (i.pix())) * 19260817);
//            HASH %= (long long) (1e9 + 7);
//        } else {
//            Log() << fileName << p.tot_lines() << " not match!!!  e:" << tot;
//
//        }
//    }
//    if (tot == p.tot_lines()) {
//        return tot;
//    } else {
//        Log() << fileName << p.tot_lines() << " not match!!! real exist:" << tot;
//        return -1 * tot;
//    }
//    //Log()<<fileName<<p.tot_lines();
////    p.ParsePartialFromIstream()
//    return p.tot_lines();
//}

//void ZQQ_AT_CHINA_VO::GaiaTest::loadUniqlist(std::string filename) {
//    uniqlist.clear();
//    // auto filelist = get_file_list("/home/smart/sdb/gaia/");
//    // std::string filename="/home/sysadmin/tempPBF/3000/";
//    ;
//
//    auto filelist = get_file_list(filename);
//    long long tot = 0;
//    ZQQ_AT_CHINA_VO::Log() << "file num" << filelist.size();
//    //sort(l.begin(),l.end());
//    std::vector<int> filelistint;
//    for (auto &file: filelist) {
//        //Log() << file;
//        std::string s;
//        if (file[filename.size() + 2] == '/')
//            s = file.substr(filename.size() + 3, file.size() - 4).c_str();
//        else
//            s = file.substr(filename.size() + 4, file.size() - 4).c_str();
//        int fileid = atoi(s.c_str());
//        Log() << file<<" "<<fileid;
////        filelistint.push_back(fileid);
//        uniqlist[fileid] = file;
//    }
//}
void ZQQ_AT_CHINA_VO::GaiaTest::loadUniqlist(const std::string &filepath, const std::string &filelist) {
    uniqlist.clear();

    std::ifstream file(filelist);
    if (!file.is_open()) {
        Log() << ("Failed to open filelist.txt");
        return;
    }

    std::vector<int> filelistint;
    std::string line;
    while (std::getline(file, line)) {
        int fileid = std::stoi(line);
        //Log() <<("Read file id: " + std::to_string(fileid));
        uniqlist[fileid] =
                filepath + std::to_string(fileid) + ".csv"; // Replace "placeholder" with actual file path if needed
    }
    file.close();

    Log() << ("Total files: " + std::to_string(uniqlist.size()));
}

std::vector<std::string> ZQQ_AT_CHINA_VO::GaiaTest::NeedtoReadFileList(std::vector<int> pixlist) {
    std::map<int, std::string> list;
    list.clear();
    for (auto &i: pixlist) {
        int level = nside2Level(NSIDE);
        while (level >= 0) {
            int uniq = calUNIQ(level, i, nside2Level(NSIDE));
            if (!uniqlist[uniq].empty()) {
                if (list[uniq].empty()) {
                    list[uniq] = uniqlist[uniq];
                }
                break;
            } else {
                level--;
            }
        }
    }
    std::vector<std::string> filelist;
    filelist.clear();
    for (auto &i: list) {
        ZQQ_AT_CHINA_VO::Log() << "add uniq" << i.first << i.second;
        filelist.push_back(i.second);
    }
    return filelist;
}

extern std::vector<int> query_disc(int Nside, double ra, double dec, double range, bool r);

//CrossMatch::Results
//ZQQ_AT_CHINA_VO::GaiaTest::QueryClosestReturnPBF(double ra, double dec, double radiusInDegree, int limit) {
//    if (uniqlist.empty()) {
//        CrossMatch::Results r;
//        r.set_status(-1);
//        r.set_match_num(0);
//        r.set_error("Catalog File list empty!");
//        ZQQ_AT_CHINA_VO::Log() << "Catalog File list empty!";
//        return r;
//    }
//    if (!checkRaDecLegal(ra, dec)) {
//        CrossMatch::Results r;
//        r.set_status(-4);
//        r.set_match_num(0);
//        r.set_error("RA DEC error ! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec));
//        ZQQ_AT_CHINA_VO::Log() << "RA DEC error! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
//        return r;
//    }
//
//    ZQQ_AT_CHINA_VO::Log() << "uniq file list size:" << uniqlist.size();
//    double query_ra = ra;
//    double query_dec = dec;
//    double R = radiusInDegree * 3.14159265359 / 360.0;
//    auto pixlist = query_disc(NSIDE, query_ra, query_dec, R, false);
//    auto fl = NeedtoReadFileList(pixlist);
//    if (fl.size() != 0) {
//        ZQQ_AT_CHINA_VO::GaiaTest d("GAIAdr2", "no", 200000000);
//
//        for (auto &f: fl) {
//            d.input(f);
//        }
//        d.getCatlogSize();
//        ZQQ_AT_CHINA_VO::CatalogCalMaster cal;
//        cal.addCatalog(d);
//        cal.buildAll();
//        CatalogCalMaster::kd_query_result res = cal.QueryClosetWithReturn(query_ra, query_dec);
//        ZQQ_AT_CHINA_VO::Log() << std::to_string(res.dis);
//        CrossMatch::Results r;
//        if (res.dis >= 360.0 || res.dis < 0) {
//            r.set_status(-3);
//            r.set_match_num(0);
//            r.set_error("ERROR in query ! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec));
//            ZQQ_AT_CHINA_VO::Log() << "ERROR in query! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
//            return r;
//        }
//        r.set_status(1);
//        r.set_match_num(1);
//        r.set_tot_read(d.getCatlogSize());
//        r.set_error("Success");
//        auto m = r.add_details();
//        m->set_tot_matched_num(1);
//        m->set_ra(ra);
//        m->set_dec(dec);
//
//        auto line = m->add_lines();
//
//        line->add_line(res.Line->content);
//        line->set_pix(res.Line->getPix());
//        line->set_ra(res.Line->getRa());
//        line->set_dec(res.Line->getDec());
//        return r;
//
//    } else {
//        CrossMatch::Results r;
//        r.set_status(-2);
//        r.set_match_num(0);
//        r.set_error("Query File empty! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec));
//        ZQQ_AT_CHINA_VO::Log() << "Query File empty! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
//        return r;
//    }
//
//    CrossMatch::Results r;
//    r.set_status(-3);
//    r.set_match_num(0);
//    r.set_error("??? IMPOSSIBLE");
//    ZQQ_AT_CHINA_VO::Log() << "IMPOSSIBLE QUERY";
//    return r;
//}
//
//CrossMatch::Results
//ZQQ_AT_CHINA_VO::GaiaTest::QueryThresholdReturnPBF(double ra, double dec, double radiusInDegree, int limit) {
//    if (uniqlist.empty()) {
//        CrossMatch::Results r;
//        r.set_status(-1);
//        r.set_match_num(0);
//        r.set_error("Catalog File list empty!");
//        ZQQ_AT_CHINA_VO::Log() << "Catalog File list empty!";
//        return r;
//    }
//    if (!checkRaDecLegal(ra, dec) || radiusInDegree > 360) {
//        CrossMatch::Results r;
//        r.set_status(-4);
//        r.set_match_num(0);
//        r.set_error("parameter error ! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec) + " radius: " +
//                    std::to_string(radiusInDegree));
//        ZQQ_AT_CHINA_VO::Log() << "parameter error! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
//        return r;
//    }
//    ZQQ_AT_CHINA_VO::Log() << "uniq file list size:" << uniqlist.size();
//    double query_ra = ra;
//    double query_dec = dec;
//    double R = radiusInDegree * 3.14159265359 / 360.0;
//    auto pixlist = query_disc(NSIDE, query_ra, query_dec, R, false);
//    auto fl = NeedtoReadFileList(pixlist);
//    if (fl.size() != 0) {
//        ZQQ_AT_CHINA_VO::GaiaTest d("GAIAdr2", "no", 200000000);
//
//        for (auto &f: fl) {
//            d.input(f);
//        }
//        d.getCatlogSize();
//        ZQQ_AT_CHINA_VO::CatalogCalMaster cal;
//        cal.addCatalog(d);
//        cal.setSearchRadius(radiusInDegree);
//        cal.buildAll();
//        auto res = cal.QueryThresholdWithReturn(ra, dec, limit);
//        ZQQ_AT_CHINA_VO::Log() << std::to_string(res.size());
//        CrossMatch::Results r;
//        r.set_status(2);
//        r.set_match_num(1);
//        r.set_tot_read(d.getCatlogSize());
//        r.set_error("Success");
//        auto m = r.add_details();
//        m->set_tot_matched_num(res.size());
//        m->set_ra(ra);
//        m->set_dec(dec);
//        for (auto &i: res) {
//            auto line = m->add_lines();
//            line->add_line(i.Line->content);
//            line->set_pix(i.Line->getPix());
//            line->set_ra(i.Line->getRa());
//            line->set_dec(i.Line->getDec());
//        }
//
//        return r;
//
//    } else {
//        CrossMatch::Results r;
//        r.set_status(-2);
//        r.set_match_num(0);
//        r.set_error("Query File empty! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec));
//        ZQQ_AT_CHINA_VO::Log() << "Query File empty! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
//        return r;
//    }
//
//    CrossMatch::Results r;
//    r.set_status(-3);
//    r.set_match_num(0);
//    r.set_error("??? IMPOSSIBLE");
//    ZQQ_AT_CHINA_VO::Log() << "IMPOSSIBLE QUERY";
//    return r;
//    return CrossMatch::Results();
//}

std::vector<std::string>
ZQQ_AT_CHINA_VO::GaiaTest::QueryThreshold(double ra, double dec, double radiusInDegree, int limit) {
    if (uniqlist.empty()) {

        ZQQ_AT_CHINA_VO::Log() << "Catalog File list empty!";
        return {};
    }
    if (!checkRaDecLegal(ra, dec) || radiusInDegree > 360) {
//        CrossMatch::Results r;

        ZQQ_AT_CHINA_VO::Log() << "parameter error! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
        return {};
    }
    ZQQ_AT_CHINA_VO::Log() << "uniq file list size:" << uniqlist.size();
    double query_ra = ra;
    double query_dec = dec;
    double R = radiusInDegree * 3.14159265359 / 180.0;
    auto pixlist = query_disc(NSIDE, query_ra, query_dec, R, false);
    auto fl = NeedtoReadFileList(pixlist);

    if (fl.size() != 0) {
        ZQQ_AT_CHINA_VO::GaiaTest d("GAIAdr3", "no", 350000000);
        d.setRaLocDefault();

        for (auto &f: fl) {
            ZQQ_AT_CHINA_VO::Log() << f;
            d.input(f);
        }
        ZQQ_AT_CHINA_VO::Log() << "catalog size" << d.getCatlogSize();
        ZQQ_AT_CHINA_VO::CatalogCalMaster cal;
        cal.addCatalog(d);
        cal.setSearchRadius(radiusInDegree);
        cal.buildAll();
        auto res = cal.QueryThresholdWithReturn(ra, dec, limit);
        ZQQ_AT_CHINA_VO::Log() << std::to_string(res.size());
        //ZQQ_AT_CHINA_VO::Log()<<res.dis<<res.Line->content;
        std::vector<std::string> ret;
        for (auto it = res.rbegin(); it != res.rend(); ++it) {
            ret.push_back(it->Line->content);
        }

        return ret;


    } else {

        ZQQ_AT_CHINA_VO::Log() << "Query File empty! Check ra :" + std::to_string(ra) + " dec:" + std::to_string(dec);
        return {};
    }


    ZQQ_AT_CHINA_VO::Log() << "IMPOSSIBLE QUERY";

}

void ZQQ_AT_CHINA_VO::GaiaTest::getHeader(std::string filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        Log() << ("Failed to open header");
        return;
    }

    std::string line;
    std::getline(file, line);
    header = line;
    file.close();

}

std::vector<std::string> ZQQ_AT_CHINA_VO::GaiaTest::MatchCatalogClosest(ZQQ_AT_CHINA_VO::Catalog c) {
    double R = 0.1 * 3.14159265359 / 180.0;
    std::vector<int> pixlist_full;
    for (auto &i: c.getCatalog()) {
        auto ra = i->getRa();
        auto dec = i->getDec();
        auto pix = i->getPix();
        auto pixlist = query_disc(NSIDE, ra, dec, R, false);
        pixlist_full.insert(pixlist_full.end(), pixlist.begin(), pixlist.end());
    }
    std::sort(pixlist_full.begin(), pixlist_full.end());
    // 使用unique去重，unique返回一个去重后的尾部迭代器
    auto last = std::unique(pixlist_full.begin(), pixlist_full.end());
    // 删除重复元素
    pixlist_full.erase(last, pixlist_full.end());

    auto fl = NeedtoReadFileList(pixlist_full);

    if (fl.size() != 0) {
        ZQQ_AT_CHINA_VO::GaiaTest d("GAIAdr3", "no", 350000000);
        d.setRaLocDefault();

        for (auto &f: fl) {
            ZQQ_AT_CHINA_VO::Log() << f;
            d.input(f);
        }
        ZQQ_AT_CHINA_VO::Log() << "Gaia catalog size" << d.getCatlogSize();
        ZQQ_AT_CHINA_VO::CatalogCalMaster cal;
        cal.addCatalog(d);
        cal.buildAll();
        ZQQ_AT_CHINA_VO::Log() << "Tree built, start match";

//        std::vector<std::string> ret;
//        int tot=0;
//        for (auto &i: c.getCatalog()) {
//
//            auto ra = i->getRa();
//            auto dec = i->getDec();
//            auto res=cal.QueryClosetWithReturn(ra,dec);
//            std::string result=std::to_string(res.dis)+","+i->content+","+res.Line->content;
//            ret.push_back(result);
//            tot++;
//            if(tot%1000==0){
//                ZQQ_AT_CHINA_VO::Log() <<"matching "<<tot<<"/"<<c.getCatlogSize();
//            }
//        }
//        ZQQ_AT_CHINA_VO::Log() <<"matching "<<tot<<"/"<<c.getCatlogSize();
        auto ret = MatchCatalogClosestMT(c, cal);
        return ret;


    } else {

        ZQQ_AT_CHINA_VO::Log() << "Query File empty! ";
        return {};
    }

    return std::vector<std::string>();
}

void
worker(std::queue<std::shared_ptr<ZQQ_AT_CHINA_VO::Catalog::line>> &taskQueue, ZQQ_AT_CHINA_VO::CatalogCalMaster &cal,
       std::vector<std::string> &ret, int &tot, std::mutex &queueMtx, std::mutex &retMtx, std::condition_variable &cv,
       bool &done) {
    std::vector<std::string> rets;
    std::queue<std::shared_ptr<ZQQ_AT_CHINA_VO::Catalog::line>> thread_taskQueue;
    while (true) {
        std::shared_ptr<ZQQ_AT_CHINA_VO::Catalog::line> task = nullptr;
        {
            std::unique_lock<std::mutex> lock(queueMtx);
            cv.wait(lock, [&taskQueue, &done] { return !taskQueue.empty() || done; });
            while (!taskQueue.empty() && thread_taskQueue.size() < 300) {
                task = taskQueue.front();
                taskQueue.pop();
                thread_taskQueue.push(task);
            }

        }
//        ZQQ_AT_CHINA_VO::Log() << "thread_taskQueue size " << thread_taskQueue.size() << "/" << tot ;
        while (!thread_taskQueue.empty()) {
            task = thread_taskQueue.front();
            thread_taskQueue.pop();
            auto ra = task->getRa();
            auto dec = task->getDec();
            auto res = cal.QueryClosetWithReturn(ra, dec);
            double dis = res.dis * 3600;
            std::string result = std::to_string(dis) + "," + task->content + "," + res.Line->content;

            {
                rets.push_back(result);

            }
        }
        if (taskQueue.empty() && done&&thread_taskQueue.empty()) {
            break;
        }

    }
    {
        std::lock_guard<std::mutex> lock(retMtx);
        ret.insert(ret.end(), rets.begin(), rets.end());
        ZQQ_AT_CHINA_VO::Log() << "Total matching " << ret.size() << "/" << tot;
        rets.clear();
    }


}

std::vector<std::string>
ZQQ_AT_CHINA_VO::GaiaTest::MatchCatalogClosestMT(ZQQ_AT_CHINA_VO::Catalog c, ZQQ_AT_CHINA_VO::CatalogCalMaster cal) {
    std::vector<std::string> ret;
    int tot = 0;
    std::queue<std::shared_ptr<Catalog::line>> taskQueue;
    std::mutex queueMtx, retMtx;
    std::condition_variable cv;
    bool done = false;
    unsigned int numCores = std::thread::hardware_concurrency();
    if (numCores == 0) {
        numCores = 1; // 如果无法检测到核心数，默认为1
    }
    ZQQ_AT_CHINA_VO::Log() << "Detected CPU Cores " << numCores;
    //numCores=4;
    std::vector<std::thread> threads;
    threads.reserve(numCores); // 预留空间
    // Fill the task queue
    for (auto &item: c.getCatalog()) {
        taskQueue.push(item);
    }
    tot = c.getCatlogSize();

    for (int i = 0; i < numCores; ++i) {
        threads.emplace_back(worker, std::ref(taskQueue), std::ref(cal), std::ref(ret), std::ref(tot),
                             std::ref(queueMtx), std::ref(retMtx), std::ref(cv), std::ref(done));
    }

    {
        std::lock_guard<std::mutex> lock(queueMtx);
        done = true;
    }
    cv.notify_all();

    for (auto &t: threads) {
        t.join();
    }

    return ret;
}

ZQQ_AT_CHINA_VO::GaiaTest::GaiaTest(std::string name, std::string strindexType, long long int MEMORY_SIZE_LIMIT,
                                    std::string filepath, std::string filelist) {
    header_line = "solution_id,designation,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,parallax_over_error,pm,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_bad_obs_al,astrometric_gof_al,astrometric_chi2_al,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_params_solved,astrometric_primary_flag,nu_eff_used_in_astrometry,pseudocolour,pseudocolour_error,ra_pseudocolour_corr,dec_pseudocolour_corr,parallax_pseudocolour_corr,pmra_pseudocolour_corr,pmdec_pseudocolour_corr,astrometric_matched_transits,visibility_periods_used,astrometric_sigma5d_max,matched_transits,new_matched_transits,matched_transits_removed,ipd_gof_harmonic_amplitude,ipd_gof_harmonic_phase,ipd_frac_multi_peak,ipd_frac_odd_win,ruwe,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,duplicated_source,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_flux_over_error,phot_g_mean_mag,phot_bp_n_obs,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_bp_mean_flux_over_error,phot_bp_mean_mag,phot_rp_n_obs,phot_rp_mean_flux,phot_rp_mean_flux_error,phot_rp_mean_flux_over_error,phot_rp_mean_mag,phot_bp_rp_excess_factor,phot_bp_n_contaminated_transits,phot_bp_n_blended_transits,phot_rp_n_contaminated_transits,phot_rp_n_blended_transits,phot_proc_mode,bp_rp,bp_g,g_rp,radial_velocity,radial_velocity_error,rv_method_used,rv_nb_transits,rv_nb_deblended_transits,rv_visibility_periods_used,rv_expected_sig_to_noise,rv_renormalised_gof,rv_chisq_pvalue,rv_time_duration,rv_amplitude_robust,rv_template_teff,rv_template_logg,rv_template_fe_h,rv_atm_param_origin,vbroad,vbroad_error,vbroad_nb_transits,grvs_mag,grvs_mag_error,grvs_mag_nb_transits,rvs_spec_sig_to_noise,phot_variable_flag,l,b,ecl_lon,ecl_lat,in_qso_candidates,in_galaxy_candidates,non_single_star,has_xp_continuous,has_xp_sampled,has_rvs,has_epoch_photometry,has_epoch_rv,has_mcmc_gspphot,has_mcmc_msc,in_andromeda_survey,classprob_dsc_combmod_quasar,classprob_dsc_combmod_galaxy,classprob_dsc_combmod_star,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,mh_gspphot,mh_gspphot_lower,mh_gspphot_upper,distance_gspphot,distance_gspphot_lower,distance_gspphot_upper,azero_gspphot,azero_gspphot_lower,azero_gspphot_upper,ag_gspphot,ag_gspphot_lower,ag_gspphot_upper,ebpminrp_gspphot,ebpminrp_gspphot_lower,ebpminrp_gspphot_upper,libname_gspphot";
    loadUniqlist(filepath,filelist);
}
