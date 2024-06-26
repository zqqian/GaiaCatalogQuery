//
// Created by zqqia on 6/3/2022.
//

#ifndef CROSSMATCHING_CATALOG_H
#define CROSSMATCHING_CATALOG_H

#include "queue"
#include <iostream>
#include <string>
#include "vector"
#include "Log.h"
#include "climits"
#include <memory>

#if  defined(__GNUC__) || defined(__linux__)

#include "dirent.h"

#endif

#if defined(__linux__)
#include "sys/io.h"
#include <unistd.h>
#include <sys/stat.h>
#include <fstream>


#endif
#if defined(_WIN32) || defined(_WIN64)

#include "io.h"
#include <direct.h>
#include <fstream>

#endif

#define NSIDE 1024
namespace ZQQ_AT_CHINA_VO {
    class Catalog {
    public:
        Catalog();

        Catalog(std::string name);

        Catalog(std::string name, double epoch);

        Catalog(std::string name, std::string strindexType);

        Catalog(std::string name, std::string strindexType, long long MEMORY_SIZE_LIMIT);

        Catalog(std::string name, std::string strindexType, double epoch);


        Catalog(std::string name, std::string strindexType, double epoch, long long MEMORY_SIZE_LIMIT);

        class line {
        public:
            line(double ra, double dec);

            line(const int pix, const double ra, const double dec, const std::string content, const long long id) : pix(
                    pix), ra(ra),
                                                                                                                    dec(dec),
                                                                                                                    content(content),
                                                                                                                    id(id) {}

            line(const int &pix, const double &ra, const double &dec, const std::string &content,
                 const std::vector<std::string> &content_details) : pix(pix), ra(ra), dec(dec), content(content),
                                                                    content_details(content_details) {}

            line(const line &l);

            line();

            std::string content;
            std::vector<std::string> content_details;

            inline void setRa(const double &e) {
                if (!(e >= 0 && e < 360)) {
                    Log() << "ra value error!" << e << Log::ERROR;
                }
                ra = e;
            }

            inline void setDec(const double &e) {
                if (!(e > -90 && e < 90)) {
                    Log() << "dec value error!" << e << Log::ERROR;
                }
                dec = e;
            }

            inline double getRa() const {
                return ra;
            }

            inline double getDec() const {
                return dec;
            }

//            inline void setPmr(const double &e) {
//                properMotion.pmr = e;
//            }

//            inline void setPmd(const double &e) {
//                properMotion.pmd = e;
//            }

            inline double getPmr() const {

                return properMotion.pmr;
            }

            inline double getPmd() const {
                return properMotion.pmd;
            }

//            inline void setPx(const double &e) {
//                properMotion.px = e;
//            }

//            inline void setRv(const double &e) {
//                properMotion.rv = e;
//            }

            inline double getPx() const {
                return properMotion.px;
            }

            inline double getRv() const {
                return properMotion.rv;
            }

            inline void setPix(const int &e) {
                pix = e;

            }

            inline int getPix() const {
                return pix;
            }

            bool setProperMotion(const double &pmr, const double &pmd, const double &px, const double &rv) {
                if (pmr != 0 || pmd != 0 || px != 0 || rv != 0) {
                    properMotion.pmr = pmr;
                    properMotion.pmd = pmd;
                    properMotion.px = px;
                    properMotion.rv = rv;
                    properMotion.haveProperMotion = true;
                    return true;
                }
                return false;
            }

            long long getId() const {
                return id;
            }

            void setId(long long id) {
                this->id = id;
            }

        private:
            double ra, dec;
            int pix;
            long long id;
            //double pmr, pmd, px, rv;
            struct pm {
                double pmr, pmd, px, rv;
                bool haveProperMotion = false;
            } properMotion;
        };


        std::string name;
        std::string header_line;
        std::vector<std::string> header_vector;
        int split_tot_write = 0;
        long long MEMORY_SIZE_LIMIT;

        std::vector<int> Tree;
        std::vector<int> pix_list_tot;
        std::vector<std::vector<int> > node2pixid;

        struct node {
            int pix;
            int val;
            int level;
        };
        std::vector<node> split_list;


        std::string MergedDir;

        static int nside2Level(int nside);

        static void UNIQ2LevelPix(int uniq, int &level, int &pix);

        void setMemoryLimit(long long limit);

        long long getMemoryLimit() const;

        // void addLine(line &l);

        void addLine(const int &pix, const double &ra, const double &dec, const std::string &content);

        void addLine(std::shared_ptr<Catalog::line> l);

        void addIndex(const int &pix, const int &id);

        virtual int wirteCatalog(int nowlevel, int ipix, std::string Dir);

        void treeBuild();

        void clear();

        void split();

        line &getLineById(const int &id);

        inline void setEpoch(const double &e) {
            haveEpoch = true;
            epoch = e;
        }

        int calUNIQ(int k, int npix, int regional);

        inline double getEpoch() const {
            return epoch;
        }

        void setCatalog(std::vector<std::shared_ptr<Catalog::line >> c);

        void setThreshold(int t);

        int epochConvert(double targetEpoch);

        virtual void input(const std::string &fileName);

        virtual void inputDir(const std::string &dirName);

        virtual void output(const std::string &fileName);

        virtual void output();

        virtual void outputDir(const std::string &dirName);

        virtual void readCatalogSingleFileFromTemp(int uniq);

        static bool checkRaDecLegal(double ra, double dec) {
            if (!(ra > 0 && ra < 360))
                return false;
            if (!(dec > -90 && dec < 90))
                return false;
            return true;
        }


        std::vector<std::shared_ptr<Catalog::line >> &getCatalog() {
            return catalog;
        }

        std::vector<std::shared_ptr<Catalog::line >> getPixbyUniq(int);

        std::vector<std::shared_ptr<Catalog::line >> getPixbyLevelPix(int, int);

        void writeHeader(std::string name = "header.csv");

        long long getCatlogSize() const {
            return catalog.size();
        }


        void treeBigger();

        void setSplitDir(const std::string &dir) {
            SplitedDir = dir;
            Log() << "set split dir to " << dir;
            dir_check(SplitedDir);
        }

        virtual void write_to_file(std::string filepath, std::string filename, std::string content);

        static std::string read_line(std::string filename);


        void mergeTemp(std::vector<int> &uniqList);


        std::string getSplitedDir() {
            return SplitedDir;
        }


    private:
        double epoch;
        enum {
            LEFT,
            MID,
            RIGHT
        } SPLIT_TYPE = MID;//DO  NOT CHANGE THIS BEFORE YOU KNOW WHAT YOU ARE DOING
        bool haveEpoch;
        enum {
            NO_INDEX,
            DYNAMIC_INDEX,
            STATIC_INDEX
        } indexType;
        enum {
            NO_TREE,
            STATIC_TREE,
            DYNAMIC_TREE
        } TreeType;
        long long totInMemory = 0;
        long long totCatalog = 0;
        int singleTreeNodeMaxElements = 100;
        std::vector<std::shared_ptr<Catalog::line>> catalog;

        std::string SplitedDir;

        int check(int now_level, int ipix);

        int checkRIGHT(int now_level, int ipix);


        bool checkTreeNeedBig(int uniq);

        int mergeTempCheck(std::vector<int> &uniqList, int nowLevel, int ipix);

        void dir_check(std::string &dir);

    protected:
        bool SPLIT_FUNCTION = true;

        int findCatalogtoWrite(std::string &content, int nowlevel, int ipix);

        int findCatalogtoWrite(std::vector<std::shared_ptr<Catalog::line>> &v, int nowlevel, int ipix);

        void treeAdd(int pix, int val = 1);

        int treeLevelMax;
        //current
        int currentTreeLevel;
    };
}


#endif //CROSSMATCHING_CATALOG_H
