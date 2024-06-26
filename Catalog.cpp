//
// Created by zqqia on 6/3/2022.
//

#include "Catalog.h"
#include "healpix/healpix.h"
#include "sofa/sofa.h"

ZQQ_AT_CHINA_VO::Catalog::line::line(const double ra, const double dec) : ra(ra), dec(dec) {
    if (checkRaDecLegal(ra, dec)) {
        setRa(ra);
        setDec(dec);
        pix = eq2pix_nest(NSIDE, ra, dec);
        this->setPix(pix);
    } else {
        throw std::invalid_argument("ra or dec is not legal");
    }

}

ZQQ_AT_CHINA_VO::Catalog::line::line() {

}

ZQQ_AT_CHINA_VO::Catalog::line::line(const ZQQ_AT_CHINA_VO::Catalog::line &l) {
    this->setRa(l.getRa());
    this->setDec(l.getDec());
    this->setPix(l.getPix());
    this->setId(l.getId());
    this->setProperMotion(l.getPmr(), l.getPmd(), l.getPx(), l.getRv());
}

/*
 * return value:
 * 0 : convert success
 *
 * */
int ZQQ_AT_CHINA_VO::Catalog::epochConvert(double targetEpoch) {
    if (!haveEpoch)return 1;
    if (targetEpoch == epoch) {
        return 0;
    }
    Log() << Log::DEBUG << "start convert epoch";
    double sourceDjm0;
    double sourceDjm;
    iauEpb2jd(epoch, &sourceDjm0, &sourceDjm);
    double targetDjm0;
    double targetDjm;
    iauEpb2jd(targetEpoch, &targetDjm0, &targetDjm);
    Log() << Log::DEBUG << "sources epoch" << epoch << sourceDjm0 << sourceDjm << "target" << targetEpoch << targetDjm0
          << targetDjm;
    for (auto &l: catalog) {
        //int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
        //              double px1, double rv1,
        //              double ep1a, double ep1b, double ep2a, double ep2b,
        //              double *ra2, double *dec2, double *pmr2, double *pmd2,
        //              double *px2, double *rv2);
        double ra2, dec2, pmr2, pmd2, px2, rv2;
        int res = iauPmsafe(l->getRa(), l->getDec(), l->getPmr(), l->getPmd(), l->getPx(), l->getRv(), sourceDjm0,
                            sourceDjm,
                            targetDjm0, targetDjm, &ra2, &dec2, &pmr2, &pmd2, &px2, &rv2);
        if (res == 0) {
            l->setRa(ra2);
            l->setDec(dec2);
            l->setProperMotion(pmr2, pmd2, px2, rv2);

        } else {
            Log() << Log::ERROR << "epoch convert failed! code:" << res << " content:" << l->content;

        }
    }
    epoch = targetEpoch;
    return 0;
}

void ZQQ_AT_CHINA_VO::Catalog::input(const string &fileName) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
}

void ZQQ_AT_CHINA_VO::Catalog::inputDir(const string &dirName) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
}

void ZQQ_AT_CHINA_VO::Catalog::output(const string &fileName) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
}

void ZQQ_AT_CHINA_VO::Catalog::outputDir(const string &dirName) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
}

void ZQQ_AT_CHINA_VO::Catalog::output() {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
}

void ZQQ_AT_CHINA_VO::Catalog::setMemoryLimit(long long int limit) {
    MEMORY_SIZE_LIMIT = limit;
}

long long ZQQ_AT_CHINA_VO::Catalog::getMemoryLimit() const {
    return MEMORY_SIZE_LIMIT;
}

void ZQQ_AT_CHINA_VO::Catalog::addLine(std::shared_ptr<Catalog::line> l) {

        l->setId(getCatlogSize());
        addIndex(l->getPix(), l->getId());


    catalog.push_back(l);

    if (!Tree.empty()&&Tree[0] > MEMORY_SIZE_LIMIT) {
        Log() << Log::ERROR << "catalog size is too large, please set a larger memory limit";
        split();
        clear();
        treeBuild();
        return;
    }
}

void ZQQ_AT_CHINA_VO::Catalog::addLine(const int &pix, const double &ra, const double &dec, const string &content) {
    std::shared_ptr<Catalog::line> l;
    l->setPix(pix);
    l->setRa(ra);
    l->setDec(dec);
    l->content = content;
    addLine(l);


}

void ZQQ_AT_CHINA_VO::Catalog::split() {
    split_list.clear();
    // split_tot_write = 0;
    if (indexType == NO_INDEX) {
        Log() << Log::ERROR << "no index type";
        return;
    } else if (indexType == DYNAMIC_INDEX||indexType == STATIC_INDEX) {
        for (int i = 0; i < 12; i++) {
            int ans;
            if (SPLIT_TYPE == RIGHT) {
                ans = checkRIGHT(0, i);
            } else {
                ans = check(0, i);
            }

            if (Tree[4 + i] - ans != 0) {

                split_list.push_back({i + 4, Tree[4 + i] - ans, 0});
                wirteCatalog(0, i, SplitedDir);
            }
        }
    } else if (indexType == STATIC_INDEX) {
        Log() << "Static index split start";
        int start = calUNIQ(currentTreeLevel, 0, currentTreeLevel);
        int end = calUNIQ(currentTreeLevel + 1, 0, currentTreeLevel + 1);
        for (int i = start; i < end; i++) {
            if (Tree[i] != 0) {
                //    Log()<<"split index:"<<i<<" value:"<<Tree[i];

                split_list.push_back({i, Tree[i], currentTreeLevel});
                wirteCatalog(currentTreeLevel, i - start, MergedDir);
            }
        }
    }
}


void ZQQ_AT_CHINA_VO::Catalog::clear() {
    catalog.clear();
    totInMemory = 0;
    totCatalog = 0;
    pix_list_tot.clear();
}

int ZQQ_AT_CHINA_VO::Catalog::calUNIQ(int k, int npix, int regional) {
    if (regional == k) {
        return 4 * pow(4, k) + npix;
    }
    if (regional > k) {
        return calUNIQ(k, npix >> 2, regional - 1);
    }
    if (regional < k) {
        return calUNIQ(k, npix << 2, regional + 1);
    }

}


void ZQQ_AT_CHINA_VO::Catalog::addIndex(const int &pix, const int &id) {
    if (indexType == NO_INDEX) {
        return;
    }
    //Log() << Log::DEBUG << "add index:" << pix;
    if (pix < 0) {
        Log() << Log::ERROR << "pix is less than 0";
        return;
    }
    if (pix > 12 * NSIDE * NSIDE) {
        Log() << Log::ERROR << "pix is larger than 12*NSIDE*NSIDE";
        return;
    }
    node2pixid[calUNIQ(treeLevelMax, pix, treeLevelMax)].push_back(id);
    pix_list_tot.push_back(pix);
    totInMemory++;
    treeAdd(pix, 1);

}

void ZQQ_AT_CHINA_VO::Catalog::treeAdd(int pix, int val) {
    int now_level = nside2Level(NSIDE);
    int now_pix = pix;
    int pix_uniq = calUNIQ(currentTreeLevel, now_pix, now_level);
    while (now_level >= 0) {
        if (now_level <= currentTreeLevel) {
            int uniq = calUNIQ(now_level, now_pix, now_level);
            if (Tree.size() > uniq) {
                Tree[uniq] += val;
            }
        }
        now_pix = now_pix >> 2;
        now_level--;
    }
    Tree[0] += val;
    if (TreeType == DYNAMIC_TREE && checkTreeNeedBig(pix_uniq)) {
        treeBigger();
    }

}

ZQQ_AT_CHINA_VO::Catalog::line &ZQQ_AT_CHINA_VO::Catalog::getLineById(const int &id) {
    return *catalog[id];
}

int ZQQ_AT_CHINA_VO::Catalog::check(int now_level, int ipix) {
    auto uniq = calUNIQ(now_level, ipix, now_level);
    int tot = 0;

    if (Tree[uniq] > singleTreeNodeMaxElements && now_level < currentTreeLevel) {

        tot += check(now_level + 1, ipix << 2);
        tot += check(now_level + 1, (ipix << 2) + 1);
        tot += check(now_level + 1, (ipix << 2) + 2);
        tot += check(now_level + 1, (ipix << 2) + 3);
    }

    if (SPLIT_TYPE == LEFT) {
        if (Tree[uniq] - tot > 0) {
            split_list.push_back({uniq, Tree[uniq] - tot, now_level});
        }
        return Tree[uniq];
    }
    else if (SPLIT_TYPE == RIGHT) {
        Log() << Log::FATAL << "SPLIT_TYPE == RIGHT not implemented";
//        if (Tree[uniq] < singleTreeNodeMaxElements) {
//            return 0;
//        } else {
//            if (Tree[uniq] - tot > 0) {
//                split_list.push_back({uniq, Tree[uniq]});
//            }
//            return Tree[uniq];
//        }


    } else if (SPLIT_TYPE == MID) {

        if (Tree[uniq] - tot > singleTreeNodeMaxElements) {
            split_list.push_back({uniq, Tree[uniq] - tot, now_level});
//            if (write_line != Tree[uniq] - tot) {
//                Log() << Log::FATAL << "write_line!=Tree[uniq] - tot" << write_line << " " << Tree[uniq] - tot;
//            }
            wirteCatalog(now_level, ipix, SplitedDir);
            return 0;
        } else {
            return 0;
        }
    }
}

void ZQQ_AT_CHINA_VO::Catalog::treeBigger() {
    if (currentTreeLevel >= treeLevelMax) {
        //  Log() << Log::ERROR << "currentTreeLevel is larger than treeLevelMax,stop bigger" << currentTreeLevel << " "
        //       << treeLevelMax;
        return;
    }
    Log() << Log::NOTIFY << "tree need big" << "current size:" << Tree.size() << " Tree level:" << currentTreeLevel;

    currentTreeLevel++;
    // tree level 0 = 16
    Tree.resize(pow(4, currentTreeLevel + 2));
    for (int i = 0; i < 12 * pow(4, currentTreeLevel); i++) {
        Tree.push_back(0);
    }
    for (auto &i: pix_list_tot) {
        Tree[calUNIQ(currentTreeLevel, i, nside2Level(NSIDE))]++;
    }
    Log() << "big success ,tree size:" << Tree.size() << "tree level:" << currentTreeLevel;
}

int ZQQ_AT_CHINA_VO::Catalog::checkRIGHT(int now_level, int ipix) {
    auto uniq = calUNIQ(now_level, ipix, now_level);
    int tot = 0;
    if (Tree[uniq] > singleTreeNodeMaxElements) {
        if (now_level < currentTreeLevel &&
            Tree[calUNIQ(now_level + 1, ipix << 2, now_level + 1)] > singleTreeNodeMaxElements
            && Tree[calUNIQ(now_level + 1, (ipix << 2) + 1, now_level + 1)] > singleTreeNodeMaxElements
            && Tree[calUNIQ(now_level + 1, (ipix << 2) + 2, now_level + 1)] > singleTreeNodeMaxElements
            && Tree[calUNIQ(now_level + 1, (ipix << 2) + 3, now_level + 1)] > singleTreeNodeMaxElements
                ) {

            tot += checkRIGHT(now_level + 1, ipix << 2);
            tot += checkRIGHT(now_level + 1, (ipix << 2) + 1);
            tot += checkRIGHT(now_level + 1, (ipix << 2) + 2);
            tot += checkRIGHT(now_level + 1, (ipix << 2) + 3);
            return Tree[uniq];
        } else {
            split_list.push_back({uniq, Tree[uniq], now_level});
            return Tree[uniq];
        }
    } else {
        return 0;
    }
}

int ZQQ_AT_CHINA_VO::Catalog::wirteCatalog(int nowlevel, int ipix, std::string Dir) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
    return 0;
}

int ZQQ_AT_CHINA_VO::Catalog::findCatalogtoWrite(std::string &content, int nowlevel, int ipix) {
    int uniq = calUNIQ(nowlevel, ipix, nowlevel);
    int tot = 0;

    for (auto &i: node2pixid[uniq]) {

        content += catalog[i]->content;

        content += "\n";

        split_tot_write++;
        tot++;
        treeAdd(catalog[i]->getPix(), -1);
        catalog[i].reset();
        catalog[i] = NULL;
    }
    node2pixid[uniq].clear();
    if (nowlevel < treeLevelMax) {
        tot += findCatalogtoWrite(content, nowlevel + 1, ipix << 2);
        tot += findCatalogtoWrite(content, nowlevel + 1, (ipix << 2) + 1);
        tot += findCatalogtoWrite(content, nowlevel + 1, (ipix << 2) + 2);
        tot += findCatalogtoWrite(content, nowlevel + 1, (ipix << 2) + 3);

    }
    return tot;

}
int ZQQ_AT_CHINA_VO::Catalog::findCatalogtoWrite(vector<std::shared_ptr<Catalog::line>> &v,int nowlevel,int ipix) {
    int uniq = calUNIQ(nowlevel, ipix, nowlevel);
    int tot = 0;

    for (auto &i: node2pixid[uniq]) {
        v.push_back(catalog[i]);
        split_tot_write++;
        tot++;
        treeAdd(catalog[i]->getPix(), -1);
        catalog[i].reset();
       catalog[i] = NULL;
    }
    node2pixid[uniq].clear();
    if (nowlevel < treeLevelMax) {
        tot += findCatalogtoWrite(v, nowlevel + 1, ipix << 2);
        tot += findCatalogtoWrite(v, nowlevel + 1, (ipix << 2) + 1);
        tot += findCatalogtoWrite(v, nowlevel + 1, (ipix << 2) + 2);
        tot += findCatalogtoWrite(v, nowlevel + 1, (ipix << 2) + 3);

    }
    return tot;
}
void ZQQ_AT_CHINA_VO::Catalog::writeHeader(string name) {
    if (!header_line.empty() && !name.empty()) {
        write_to_file(SplitedDir, name, header_line);
    } else {
        Log() << Log::ERROR << "header line or header file name is empty";
    }

}


void ZQQ_AT_CHINA_VO::Catalog::readCatalogSingleFileFromTemp(int uniq) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
}

int ZQQ_AT_CHINA_VO::Catalog::nside2Level(int nside) {
    switch (nside) {
        case 1:
            return 0;
        case 2:
            return 1;
        case 4:
            return 2;
        case 8:
            return 3;
        case 16:
            return 4;
        case 32:
            return 5;
        case 64:
            return 6;
        case 128:
            return 7;
        case 256:
            return 8;
        case 512:
            return 9;
        case 1024:
            return 10;
        case 2048:
            return 11;
        case 4096:
            return 12;
        case 8192:
            return 13;
        case 16384:
            return 14;
        case 32768:
            return 15;
        case 65536:
            return 16;
        case 131072:
            return 17;
        case 262144:
            return 18;
        case 524288:
            return 19;
        case 1048576:
            return 20;
        case 2097152:
            return 21;
        case 4194304:
            return 22;
        case 8388608:
            return 23;
        case 16777216:
            return 24;
        case 33554432:
            return 25;
        case 67108864:
            return 26;
        default:
            return -1;
    }
}

void ZQQ_AT_CHINA_VO::Catalog::treeBuild() {
    Tree.clear();
    if (indexType == NO_INDEX) {
        TreeType = NO_TREE;
        currentTreeLevel = 0;
        return;
    } else if (indexType == STATIC_INDEX) {
        TreeType = STATIC_TREE;
        currentTreeLevel = nside2Level(NSIDE);
        treeLevelMax = currentTreeLevel;
        Tree.reserve(16 * NSIDE * NSIDE);
        for (int i = 0; i < 16 * NSIDE * NSIDE; i++) {
            Tree.push_back(0);
        }
    } else if (indexType == DYNAMIC_INDEX) {
        currentTreeLevel = 0;
        TreeType = DYNAMIC_TREE;
        treeLevelMax = nside2Level(NSIDE);
        for (int i = 0; i < 16; i++) {
            Tree.push_back(0);
        }
    }
    node2pixid.clear();
    node2pixid.reserve(16 * NSIDE * NSIDE);
    for (int i = 0; i < 16 * NSIDE * NSIDE; i++) {

        // Log()<<i;

        node2pixid.push_back(std::vector<int>());
    }
    Log() << "treelevel: " << currentTreeLevel;

}

void ZQQ_AT_CHINA_VO::Catalog::write_to_file(std::string filepath, std::string filename, std::string content) {
    Log() << Log::ERROR << " virtual function! readCatalogSingleFileFromTemp ERROR";
    //    // auto startTime = std::chrono::high_resolution_clock::now();
//    if (filepath.size() == 0)return;
//    if (filepath[filepath.size() - 1] != '/') {
//        filepath.append("/");
//    }
//#if defined(_WIN32) || defined(_WIN64)
//    if (_access(filepath.c_str(), 00) == -1) {
//        _mkdir(filepath.c_str());
//    }
//#endif
//#if   defined(__linux__)
//    if (access(filepath.c_str(), 0) == -1) {
//            mkdir(filepath.c_str(),S_IRWXU);
//        }
//
//#endif
//    std::ofstream outFile;
//    const size_t bufsize = 256 * 1024;
//    char buf[bufsize];
//    outFile.rdbuf()->pubsetbuf(buf, bufsize);
//    outFile.open((filepath + filename).c_str(), std::ios_base::app);
//    outFile.write(content.c_str(), content.size());
//
//    outFile.close();
//    //auto endTime = std::chrono::high_resolution_clock::now();
//
//    Log() << "Write to " << (filepath + filename) << " size " << content.size()
//          << "\n";// " cost time "<< std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms \n";
}

std::string ZQQ_AT_CHINA_VO::Catalog::read_line(std::string filename) {
    FILE *fp;
    char str[10000];

    /* 打开用于读取的文件 */
    fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        std::string err = "error occured when open " + filename;
        perror(err.c_str());
        exit(-1);
    }
    fgets(str, 10000, fp);

    fclose(fp);
    return str;
}

void ZQQ_AT_CHINA_VO::Catalog::dir_check(string &dir) {
    Log() << "Checking dir " << dir << "\n";

    if (dir[dir.size() - 1] != '/') {
        dir += '/';
    }
#if defined(_WIN32) || defined(_WIN64)
    if (_access(dir.c_str(), 00) == -1) {
        Log() << dir << " not found,now creating." << "\n";
        if (_mkdir(dir.c_str()) == 0) {
            Log() << dir << " created successfully." << "\n";
        } else {
            Log() << "ERROR:Problem creating path " << dir << " , PLEASE REFERENCE MANUAL" << "\n";
            exit(-1);
        }
    } else {
        Log() << dir << " found successfully." << "\n";

    }
#endif
#if   defined(__linux__)
    if (access(dir.c_str(), 0) == -1) {
           if( mkdir(dir.c_str(),S_IRWXU)==0){
               Log()<<dir<<" created successfully."<<"\n";
           }else{
               Log()<<"ERROR:Problem creating path "<<dir<<" , PLEASE REFERENCE MANUAL"<<"\n";
                       exit(-1);

           }

        }else{
        Log()<<dir<<" found successfully."<<"\n";
    }
#endif
    std::string s = "TEST TEST TEST ";
    srand((int) time(0));

    s = s + dir;
    s = s + std::to_string(rand());
    Log() << "Write contents: \"" << s << "\"\n";
    write_to_file(dir, "test.csv", s);

    auto read = read_line(dir + "/test.csv");
    Log() << "Read contents: \"" << read << "\"\n";
    if (read == s) {
        Log() << dir << " Write and Read test successfully.\n";
    } else {
        Log() << "ERROR write contents in" << dir << "\n";
        exit(-1);
    }
    if (remove((dir + "/test.csv").c_str()) == 0) {
        Log() << (dir + "/test.csv") << " remove successfully\n";
    } else {
        Log() << "ERROR remove" << (dir + "/test.csv") << "\n";
        exit(-1);
    }

}

ZQQ_AT_CHINA_VO::Catalog::Catalog(std::string name, std::string strindexType, double epoch,
                                  long long int MEMORY_SIZE_LIMIT) : epoch(epoch),
                                                                     haveEpoch(
                                                                             true),
                                                                     name(name),
                                                                     MEMORY_SIZE_LIMIT(
                                                                             MEMORY_SIZE_LIMIT),
                                                                     indexType(
                                                                             NO_INDEX) {
    totInMemory = 0;
    if (strindexType == "static") {
        indexType = STATIC_INDEX;
    } else if (strindexType == "dynamic") {
        indexType = DYNAMIC_INDEX;
    } else {
        indexType = NO_INDEX;
    }
    treeBuild();
    Log() << "Catalog(std::string name, std::string strindexType, double epoch, long long MEMORY_SIZE_LIMIT)";
    Log() << "epoch: " << epoch;
    Log() << "haveEpoch: " << haveEpoch;
    Log() << "name: " << name;
    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
    Log() << "indexType: " << indexType;


}

ZQQ_AT_CHINA_VO::Catalog::Catalog(std::string name, std::string strindexType, double epoch) : epoch(epoch),
                                                                                              haveEpoch(true),
                                                                                              name(name),
                                                                                              MEMORY_SIZE_LIMIT(
                                                                                                      LLONG_MAX),
                                                                                              indexType(NO_INDEX) {
    totInMemory = 0;
    if (strindexType == "static") {
        indexType = STATIC_INDEX;
    } else if (strindexType == "dynamic") {
        indexType = DYNAMIC_INDEX;
    } else {
        indexType = NO_INDEX;
    }
    treeBuild();
    Log() << "Catalog(std::string name, std::string strindexType, double epoch)";
    Log() << "epoch: " << epoch;
    Log() << "haveEpoch: " << haveEpoch;
    Log() << "name: " << name;
    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
    Log() << "indexType: " << indexType;

}

ZQQ_AT_CHINA_VO::Catalog::Catalog(std::string name, std::string strindexType, long long int MEMORY_SIZE_LIMIT) : epoch(
        2000.0),
                                                                                                                 haveEpoch(
                                                                                                                         true),
                                                                                                                 name(name),
                                                                                                                 MEMORY_SIZE_LIMIT(
                                                                                                                         MEMORY_SIZE_LIMIT),
                                                                                                                 indexType(
                                                                                                                         NO_INDEX) {
    totInMemory = 0;
    if (strindexType == "static") {
        indexType = STATIC_INDEX;
    } else if (strindexType == "dynamic") {
        indexType = DYNAMIC_INDEX;
    } else {
        indexType = NO_INDEX;
    }
//    Log() << "Catalog(std::string name, std::string strindexType, long long MEMORY_SIZE_LIMIT)";
//    Log() << "epoch: " << epoch;
//    Log() << "haveEpoch: " << haveEpoch;
//    Log() << "name: " << name;
//    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
//    Log() << "indexType: " << indexType;
    treeBuild();


}

ZQQ_AT_CHINA_VO::Catalog::Catalog(std::string name, std::string strindexType) : epoch(2000.0), haveEpoch(false),
                                                                                name(name),
                                                                                MEMORY_SIZE_LIMIT(LLONG_MAX),
                                                                                indexType(NO_INDEX) {
    totInMemory = 0;
    if (strindexType == "static") {
        indexType = STATIC_INDEX;
    } else if (strindexType == "dynamic") {
        indexType = DYNAMIC_INDEX;
    } else {
        indexType = NO_INDEX;
    }
    treeBuild();
    Log() << "Catalog(std::string name, std::string strindexType)";
    Log() << "epoch: " << epoch;
    Log() << "haveEpoch: " << haveEpoch;
    Log() << "name: " << name;
    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
    Log() << "indexType: " << indexType;
}

ZQQ_AT_CHINA_VO::Catalog::Catalog(std::string name, double epoch) : epoch(epoch), haveEpoch(true), name(name),
                                                                    MEMORY_SIZE_LIMIT(LLONG_MAX),
                                                                    indexType(NO_INDEX) {
    totInMemory = 0;
    treeBuild();
    Log() << "Catalog(std::string name, double epoch)";
    Log() << "epoch: " << epoch;
    Log() << "haveEpoch: " << haveEpoch;
    Log() << "name: " << name;
    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
    Log() << "indexType: " << indexType;
}

ZQQ_AT_CHINA_VO::Catalog::Catalog(std::string name) : epoch(2000.0), haveEpoch(false), name(name),
                                                      MEMORY_SIZE_LIMIT(LLONG_MAX),
                                                      indexType(NO_INDEX) {
    totInMemory = 0;
    totCatalog = 0;
    treeBuild();
    Log() << "Catalog(std::string name)";
    Log() << "epoch: " << epoch;
    Log() << "haveEpoch: " << haveEpoch;
    Log() << "name: " << name;
    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
    Log() << "indexType: " << indexType;
}

ZQQ_AT_CHINA_VO::Catalog::Catalog() : epoch(2000.0), haveEpoch(false), name("default catalog"),
                                      MEMORY_SIZE_LIMIT(LLONG_MAX),
                                      indexType(NO_INDEX) {
    totInMemory = 0;
    totCatalog = 0;
    treeBuild();
//    Log() << "Catalog()";
//    Log() << "epoch: " << epoch;
//    Log() << "haveEpoch: " << haveEpoch;
//    Log() << "name: " << name;
//    Log() << "MEMORY_SIZE_LIMIT: " << MEMORY_SIZE_LIMIT;
//    Log() << "indexType: " << indexType;
}

bool ZQQ_AT_CHINA_VO::Catalog::checkTreeNeedBig(int uniq) {
    if (TreeType == STATIC_TREE) {
        return false;
    }
    if (currentTreeLevel >= treeLevelMax) {
        return false;
    }
    if ((4 << (currentTreeLevel)) * 12 * singleTreeNodeMaxElements < totInMemory) {
        return true;
    }
    if (Tree[uniq] > singleTreeNodeMaxElements) {
        return true;
    }
    return false;
}

std::vector<std::shared_ptr<ZQQ_AT_CHINA_VO::Catalog::line>> ZQQ_AT_CHINA_VO::Catalog::getPixbyUniq(int uniq) {
    if(TreeType == NO_TREE){
        return std::vector<std::shared_ptr<line>>();
    }
    int k = log2((uniq) / 4) / 2;
    int p = uniq - 4 * (1 << (2 * k));
    int tot = 1;
    while (k < treeLevelMax) {
        p = p << 2;
        tot *= 4;
        k++;
    }
    std::vector<std::shared_ptr<line>> pix;

    for (int i = 0; i < tot; i++) {
       int pp=calUNIQ(treeLevelMax,i+p,treeLevelMax);
        //Log() << "pp: " << pp;
        if (!node2pixid[pp].empty())
            for (auto p: node2pixid[pp]) {
                pix.push_back(catalog[p]);
            }
    }

    return pix;
}

//!!!!deprecated
std::vector<std::shared_ptr<ZQQ_AT_CHINA_VO::Catalog::line>>
ZQQ_AT_CHINA_VO::Catalog::getPixbyLevelPix(int level, int pix) {
    std::vector<std::shared_ptr<line>> pixx;
    int tot = 1;
    while (level < treeLevelMax) {
        pix = pix << 2;
        tot *= 4;
        level++;
    }
    for (int i = pix; i < tot + pix; i++) {
        if (!node2pixid[calUNIQ(treeLevelMax, i, treeLevelMax)].empty())
            for (auto p: node2pixid[calUNIQ(treeLevelMax, i, treeLevelMax)]) {
                pixx.push_back(catalog[p]);
            }

    }
    return pixx;
}

void ZQQ_AT_CHINA_VO::Catalog::setCatalog(std::vector<std::shared_ptr<Catalog::line >> c) {
catalog=c;
    totInMemory = c.size();
    totCatalog = c.size();
}

void ZQQ_AT_CHINA_VO::Catalog::setThreshold(int t) {
    singleTreeNodeMaxElements=t;
}

void ZQQ_AT_CHINA_VO::Catalog::UNIQ2LevelPix(int uniq, int &level, int &pix) {
    int k = log2((uniq) / 4) / 2;
    int p = uniq - 4 * (1 << (2 * k));
    level = k;
    pix=p;
}









