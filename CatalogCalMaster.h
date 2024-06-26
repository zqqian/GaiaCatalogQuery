//
// Created by zqqia on 2021/12/14.
//
#pragma once
#ifndef CROSSMATCHING_CATALOGCALMASTER_H
#define CROSSMATCHING_CATALOGCALMASTER_H
#include "Catalog.h"
#include "Log.h"
#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#include  <io.h>
#endif
#if  defined(__linux__)
#include <dirent.h>
#include  <sys/io.h>
#include <unistd.h>
#endif
#include <sys/stat.h>

#include  <cstdio>
#include  <cstdlib>

#include "atomic"
#include <queue>
long eq2pix_nest(long nside, double ra, double dec);

namespace ZQQ_AT_CHINA_VO {
    class CatalogCalMaster {

    public:
        CatalogCalMaster() = default;
        CatalogCalMaster(Catalog c);
        ~CatalogCalMaster() = default;
        void addCatalog(Catalog c);
        // virtual void input(std::string name) = 0;

        virtual void output() ;
        //  std::string outputMatchResultsPBF();

        static bool checkRaDecLegal(double ra,double dec){
            if(!(ra>0&&ra<360))
                return false;
            if(!(dec>-90&&dec<90))
                return false;
            return true;
        }
        void query();

        void findAllinThreshold();

        void findClosest();

        void findClosestViolence();

        void setSearchRadius(double r);

        double getSearchRadius();
        void singleFileQueryThreshold();
        //debug only:
        long long quicktimes{0};
        long long havertimes = 0;

        void buildAll();
//         struct line {
//            int uniq;
//            double ra, dec;
//            std::string content;
//            std::vector<std::string >content_details;
//
//        };
//        struct catalogFile {
//            std::string name;
//            std::vector<std::string >header_vector;
//            std::vector<line> catalog;
//
//        };
        std::vector<Catalog> UserVector;
        void singleFileQueryClosest(std::vector<std::shared_ptr<Catalog::line>> l);
        bool queryClosest(double ra, double dec);

        //   std::vector<double>Time_cal;//temp

        struct kd_query_result {
            int id;
            double dis;
            std::shared_ptr<Catalog::line> Line;
            bool operator<(const kd_query_result &k) const {
                return dis < k.dis;
            }
        };
        kd_query_result QueryClosetWithReturn(double ra,double dec);
        std::vector<kd_query_result> QueryThresholdWithReturn(double ra,double dec,int lineLimit);

    protected:


        void build(std::vector<std::shared_ptr<Catalog::line>> BuildTreeVector);

        // std::vector<line> BuildTreeVector;
        void printRes();

        void printTot();
        long long totMatched();

        double SearchThreshold = 0.01 * (1.0 / 3600);

        void queryKdClosest(double ra, double dec, int now, kd_query_result &closest);

    private:
        double min_dis;
        enum match_type {
            ALL_IN_THRESHOLD,
            CLOSEST,
            VIOLENCE
        };

        struct wait_to_build_kd {
            double coordinate[2];
            int id;
            std::shared_ptr<Catalog::line> Line;
        };
        struct matched{
            double dis;
            std::shared_ptr<Catalog::line> l;
        };
        struct matched_line{
            Catalog::line original;
            std::vector<matched>match_details;
        };
        struct matched_file{
            std::string filename;
            std::vector<matched_line> line_details;
        };

        std::vector<matched_file> res;

        std::vector<matched_file>results;
        std::vector<wait_to_build_kd> nodelist;

        struct tree_node {
            int fa{};
            int deep{};
            int is_vaild;
            double coordinate[2]{};
            int file_id{};
            int lson{}, rson{};
            std::shared_ptr<Catalog::line> Line;
            tree_node() : is_vaild(0) {}
        };


        std::vector<tree_node> tree_node_list;

        void queryAllCatalogs(match_type type);

        void tree_build(int fa, int l, int r, int deep, int type);

        void singleFileQueryAllinThreshold(std::vector<std::shared_ptr<Catalog::line>> l, matched_file &f);

        void queryKdThreshold(double ra, double dec, int now, std::priority_queue <kd_query_result> &kp);
    };

}


#endif //CROSSMATCHING_CATALOGCALMASTER_H
