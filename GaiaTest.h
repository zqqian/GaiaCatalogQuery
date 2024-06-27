//
// Created by zqqia on 2022/08/15.
//

#ifndef CROSSMATCHING_GAIATEST_H
#define CROSSMATCHING_GAIATEST_H
#include "map"
#include "CsvRead.h"
namespace ZQQ_AT_CHINA_VO {
    class GaiaTest : public CsvRead{
    public:

        GaiaTest();
        GaiaTest(std::string name, std::string strindexType, long long MEMORY_SIZE_LIMIT);
        GaiaTest(std::string name, std::string strindexType, long long MEMORY_SIZE_LIMIT,std::string filepath,std::string filelist);
        std::vector<int>v;
        std::map<int,std::string>uniqlist;
        std::string header;
        void getHeader(std::string filename);
        long long HASH=1;
        void loadUniqlist(const std::string& filepath,const std::string& filelist);
        virtual  int wirteCatalog(int nowlevel, int ipix, std::string Dir);
        std::vector<std::string>NeedtoReadFileList(std::vector<int>pixlist);
        std::vector<std::string>QueryThreshold(double ra,double dec,double radiusInDegree,int limit=100);
        std::vector<std::string>MatchCatalogClosest(Catalog c);
        std::vector<std::string> MatchCatalogClosestMT(Catalog c, CatalogCalMaster cal);
    };

}
#endif //CROSSMATCHING_GAIATEST_H
