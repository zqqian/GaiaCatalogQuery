//
// Created by zqqia on 2022/04/25.
//

#include "Catalog.h"

void ZQQ_AT_CHINA_VO::Catalog::mergeTemp(std::vector<int> &uniqList) {
    for (int i = 0; i < 12; i++) {
        mergeTempCheck(uniqList, 0, i);
    }
}

int ZQQ_AT_CHINA_VO::Catalog::mergeTempCheck(std::vector<int> &uniqList, int nowLevel, int ipix) {
    if (nowLevel > treeLevelMax) return 0;
    auto uniq = calUNIQ(nowLevel, ipix, nowLevel);
    if (uniq > 16 * NSIDE * NSIDE || uniq < 0) {
        Log()<<Log::FATAL<<"???? Unexcepted ERROR uniq error uniq:"<<uniq<<"nowlevel:"<<nowLevel<<"ipix:"<<ipix;
        return 0;
    }

    readCatalogSingleFileFromTemp(uniq);
    if(TreeType==DYNAMIC_TREE){
        while (nowLevel >= currentTreeLevel && currentTreeLevel < treeLevelMax) {
            treeBigger();
        }
    }

    if (nowLevel < treeLevelMax) {
        mergeTempCheck(uniqList, nowLevel + 1, ipix << 2);
        mergeTempCheck(uniqList, nowLevel + 1, (ipix << 2) + 1);
        mergeTempCheck(uniqList, nowLevel + 1, (ipix << 2) + 2);
        mergeTempCheck(uniqList, nowLevel + 1, (ipix << 2) + 3);
    }
    if (Tree[uniq] > singleTreeNodeMaxElements || nowLevel == 0) {
        // split_list.push_back({uniq, Tree[uniq], nowLevel});
        int needWirte = Tree[uniq];
        int write_line = wirteCatalog(nowLevel, ipix, MergedDir);

        if (write_line != needWirte) {
            Log() << Log::FATAL << "write_line!=Tree[uniq] - tot" << write_line << " " << needWirte;
        }
        return 0;
    } else {
        return 0;
    }

}
