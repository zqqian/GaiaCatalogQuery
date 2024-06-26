//
// Created by zqqia on 2021/12/23.
//


#if defined(__linux__)
#include "sys/io.h"
#include <unistd.h>
#endif
#if defined(_WIN32) || defined(_WIN64)

#include "io.h"
#include <direct.h>

#endif


#include <fstream>
#include <string>
#include "CsvRead.h"
#include "cmath"
//#include "healpix/healpix.cpp"

#if defined(_WIN32) || defined(_WIN64)

/*
 * SUFFIX IS NOT SUPPORT AT LINUX
 * */
std::vector<std::string> ZQQ_AT_CHINA_VO::CsvRead::get_file_list(std::string dir, std::string suffix) {
    std::vector<std::string> allPath;


    std::string dir2 = dir + "*.*";
    intptr_t handle;
    _finddata_t findData;

    handle = _findfirst(dir2.c_str(), &findData);
    if (handle == -1) {
        return allPath;
    }
    do {
        if (findData.attrib & _A_SUBDIR) {
            if (strcmp(findData.name, ".") == 0 || strcmp(findData.name, "..") == 0)
                continue;
            std::string dirNew = dir + findData.name + "/";
            std::vector<std::string> tempPath = get_file_list(dirNew);
            allPath.insert(allPath.end(), tempPath.begin(), tempPath.end());
        } else {
            std::string s;
            s = findData.name;
            if (s.size() < suffix.size() + 1)continue;
            if (suffix == "") {
                std::string filePath = dir + findData.name;
                allPath.push_back(filePath);
            } else {
                if (s.substr(s.size() - (suffix.size() + 1), suffix.size() + 1) == "." + suffix) {
                    std::string filePath = dir + findData.name;
                    allPath.push_back(filePath);
                }
            }

        }
    } while (_findnext(handle, &findData) == 0);
    _findclose(handle);
    return allPath;

//    for (const auto & entry : fs::directory_iterator(dir))
//        allPath.push_back(entry.path().string());
    return allPath;
}

#endif

#if  defined(__linux__)
/*
 * SUFFIX IS NOT SUPPORT AT LINUX
 * */
std::vector<std::string> ZQQ_AT_CHINA_VO::CsvRead::get_file_list(std::string path,std::string suffix){

    DIR *dir; struct dirent *diread;
    std::vector<std::string> files;

    if ((dir = opendir(path.c_str())) != nullptr) {
        while ((diread = readdir(dir)) != nullptr) {
            std::string filename=diread->d_name;
            if(filename==".."||filename==".")continue;
            if(filename.find(".")!=filename.npos){
                files.push_back(path+filename);
            }else{
                auto f= get_file_list(path + filename + "/");
                files.insert(files.end(),f.begin(),f.end());
            }
        }
        closedir (dir);
    }
    return files;

}
#endif

std::istream &operator>>(std::istream &str, ZQQ_AT_CHINA_VO::CSVRow &data) {
    //  data.readNextRow(str);
    return str;
}

void ZQQ_AT_CHINA_VO::CsvRead::input(std::string name) {
    inputSingleCsv(name, false);
    std::cout << this->getCatlogSize() << " lines \n";
}

void ZQQ_AT_CHINA_VO::CsvRead::output() {
//    printTot();


}

void ZQQ_AT_CHINA_VO::CsvRead::inputdir(std::string name, std::string suffix) {
    if (name[name.size() - 1] != '/') {
        name += '/';
    }

    auto filelist = get_file_list(name, suffix);
    std::cout << "find " << filelist.size() << " in dir " << name << "\n";
    for (auto i: filelist) {
        std::cout << "now read " << i << " \n";
        CsvRead::inputSingleCsv(i, FIRST_LINE_HEADER, unlimitLineNumber);
       Log() << i<<" read " << CatalogTotRead<< " lines\n";
    }

}

/*
 * NOTICE:
 * line Limit also include header_vector line!
 * */
void ZQQ_AT_CHINA_VO::CsvRead::inputSingleCsv(std::string name, int headerLine, int lineLimit) {
   // Log() << "inputSingleCsv" << name;
    //    Catalog f;
    //  f.name = name;
    this->name = name;

    long long now_line = 0;//line number is start from 1

    FILE *fp;
    fp = fopen(name.c_str(), "r");

    if (fp == NULL) {
        std::string err = "error occured when open " + name;
        perror(err.c_str());
        return;
    }
    while (row.readNextRow(fp)) {
        // cout<<row.rawString()<<endl;
        if (row.rawString().size() == 0) continue;
        if (row.rawString()[0] == '#') continue;

        now_line++;
        if (CatalogTotRead % 1000000 == 0) {
            Log() << "read line " << CatalogTotRead << Log::WARNING;
        }
        if (lineLimit != unlimitLineNumber && now_line > lineLimit)break;
        if (headerLine == NO_HEADER) setRaLocDefault();
        if (now_line < headerLine) continue;

        if (now_line == headerLine) {
            header_line = row.rawString();
            this->header_vector = headerDeal();
            continue;
        }

       // Catalog::line l;
        std::shared_ptr<Catalog::line> l=std::make_shared<Catalog::line>();


        bool flag = lineDeal(*l);//check  this line
        CatalogTotRead++;
//        {//GAIA　test
//            GAIApixList.push_back(l->getPix());
//        if(GAIApixList.size()>100000000){
//            GaiaPbfTimes++;
//            writeGAIA();
//            GAIApixList.clear();
//        }
//
//        }
        if (SAVE_TO_CATALOG && flag){

            //std::unique_ptr<Catalog::line> lp(new Catalog::line(l));

            this->addLine(l);
        }



    }


   // Log() << "read " << (this->getCatlogSize()) << " lines." << Log::WARNING;
    if (SAVE_TO_CATALOG) {//UserVector.push_back(f);}





    }
    fclose(fp);
}

void ZQQ_AT_CHINA_VO::CsvRead::setRaLocDefault() {
    //raLoc = 105;
    //decLoc = 106;
    raLoc=5;
    decLoc=7;
}

//input a header, program will find ra dec automatically
void ZQQ_AT_CHINA_VO::CsvRead::setRcDecLocFromHeader(std::vector<std::string> header) {


//    int now = 0;
////    for (auto i: header_vector) {
////          cout<<i<<" "<<raname<<std::endl;
////        if (i == raname)
////            raLoc = now;
////        if (i==decname)
////            decLoc = now;
////        now++;
////    }
    raLoc = -1;
    decLoc = -1;
    for (int i = 0; i < header.size(); i++) {
        if (header[i] == raname) {
            raLoc = i;
        }
        if (header[i] == decname) {
            decLoc = i;
        }
    }
    if ((raLoc == -1 || decLoc == -1) && header.size() == 109) {
        //  setRaLocDefault();
    }
}

void ZQQ_AT_CHINA_VO::CsvRead::setRaDecNameInHeader(std::string ra_name, std::string dec_name) {
    if (ra_name.size() > 0)
        raname = ra_name;
    else
        Log() << "ra  name cannot be empty!" << Log::ERROR;
    if (dec_name.size() > 0)
        decname = dec_name;
    else
        Log() << "dec name cannot be empty!\n" << Log::ERROR;
}

ZQQ_AT_CHINA_VO::CsvRead::CsvRead() : raLoc(-1), decLoc(-1), raname("ra"), decname("dec"), row(true), HaveProperMotion(
        false), SAVE_TO_CATALOG(true) {
}

bool ZQQ_AT_CHINA_VO::CsvRead::lineDeal(ZQQ_AT_CHINA_VO::Catalog::line &l) {
    bool flag = RaDecDeal(l);



    return flag;
}

std::vector<std::string> ZQQ_AT_CHINA_VO::CsvRead::headerDeal() {

    std::vector<std::string> h;
    for (int i = 0; i < row.size(); i++) {
        h.push_back(std::string(row[i]));
//       Log() << std::string(row[i]) << Log::WARNING;
    }
    if (RaDecType == READ_FROM_HEADER) {
        setRcDecLocFromHeader(h);
        if (raLoc == -1 || decLoc == -1) {
            Log()<<header_line;
            Log() << "error can not find ra dec location!" << Log::FATAL;
        } else {
           // Log() << "ra " << raLoc << " dec " << decLoc;
        }
    }

    return h;
}

bool ZQQ_AT_CHINA_VO::CsvRead::RaDecDeal(ZQQ_AT_CHINA_VO::Catalog::line &l) {
    double ra, dec;
    if (!checkLocLegal(raLoc) || !checkLocLegal(decLoc)) {
        Log() << " raloc or decloc error!" << Log::FATAL;
    }
    ra = atof(std::string(row[raLoc]).c_str());
    dec = atof(std::string(row[decLoc]).c_str());
    if (!checkRaDecLegal(ra, dec)) {
        Log() << "ra dec not legal" << ra << dec << row.rawString();
        return false;
    }
    l.content = row.rawString();
    l.setRa(ra);
    l.setDec(dec);
    l.setPix(eq2pix_nest(NSIDE, ra, dec));
  //  double mag = atof(std::string(row[19]).c_str());
   // double bright = pow(10, mag * (-0.4));
    //pix_count[eq2pix_nest(NSIDE, ra, dec)]+=bright*100000;
    //pix_count.push_back(eq2pix_nest(NSIDE, ra, dec));
    //   cout<<l.pix<<" "<<l.ra <<" "<<l.dec<<" "<<endl;
    return true;
}

void ZQQ_AT_CHINA_VO::CsvRead::PmDeal(ZQQ_AT_CHINA_VO::Catalog::line &l) {
    double pmr, pmd, px, rv;
    if (!(checkLocLegal(pmLoc.pmdLoc) && checkLocLegal(pmLoc.pmrLoc) && checkLocLegal(pmLoc.pxLoc) &&
          checkLocLegal(pmLoc.rvLoc))) {
        Log() << " proper motion location error!" << Log::FATAL;

    }
    pmr = atof(std::string(row[pmLoc.pmrLoc]).c_str());
    pmd = atof(std::string(row[pmLoc.pmdLoc]).c_str());
    px = atof(std::string(row[pmLoc.pxLoc]).c_str());
    rv = atof(std::string(row[pmLoc.rvLoc]).c_str());
    l.setProperMotion(pmr, pmd, px, rv);

}

bool ZQQ_AT_CHINA_VO::CsvRead::checkLocLegal(int loc) {
    if (loc >= 0 && loc < row.dataSize()) {
        return true;
    }
    return false;
}

ZQQ_AT_CHINA_VO::CsvRead::CsvRead(std::string catalogname, int headerLine, int lineLimit) : Catalog(catalogname) {

}

ZQQ_AT_CHINA_VO::CsvRead::CsvRead(std::string catalogname) : Catalog(catalogname) {

}

ZQQ_AT_CHINA_VO::CsvRead::CsvRead(std::string Catalogname, std::string strindexType, long long int MEMORY_SIZE_LIMIT)
        : Catalog(Catalogname, strindexType, MEMORY_SIZE_LIMIT) {

}

ZQQ_AT_CHINA_VO::CsvRead::CsvRead(std::string Catalogname, double epoch) : Catalog(Catalogname, epoch) {

}

#if defined(_WIN32) || defined(_WIN64)

std::vector<std::string> FileList(std::string dir) {

    std::vector<std::string> allPath;


    std::string dir2 = dir + "*.*";
    intptr_t handle;
    _finddata_t findData;

    handle = _findfirst(dir2.c_str(), &findData);
    if (handle == -1) {
        return allPath;
    }
    do {
        if (findData.attrib & _A_SUBDIR) {
            if (strcmp(findData.name, ".") == 0 || strcmp(findData.name, "..") == 0)
                continue;
            std::string dirNew = dir + findData.name + "/";
            std::vector<std::string> tempPath = FileList(dirNew);
            allPath.insert(allPath.end(), tempPath.begin(), tempPath.end());
        } else {
            std::string s;
            s = findData.name;
            if (s.size() < 4)continue;
            if (s.substr(s.size() - 4, 4) == ".csv") {
                std::string filePath = dir + findData.name;
                allPath.push_back(filePath);
            }
        }
    } while (_findnext(handle, &findData) == 0);
    _findclose(handle);
    return allPath;

//    for (const auto & entry : fs::directory_iterator(dir))
//        allPath.push_back(entry.path().string());
    return allPath;


}

#endif

#if  defined(__linux__)
std::vector<std::string> FileList(std::string path){

    DIR *dir; struct dirent *diread;
    std::vector<std::string> files;

    if ((dir = opendir(path.c_str())) != nullptr) {
        while ((diread = readdir(dir)) != nullptr) {
            std::string filename=diread->d_name;
            if(filename==".."||filename==".")continue;
            if(filename.find(".")!=filename.npos){
                files.push_back(path+filename);
            }else{
                auto f= FileList(path + filename + "/");
                files.insert(files.end(),f.begin(),f.end());
            }
        }
        closedir (dir);
    }
    return files;

}
#endif


void ZQQ_AT_CHINA_VO::CsvRead::mergeCatalog(std::string dir) {
    setSplitDir(final_split_dir);
    auto filelist = get_file_list(dir);
    Log() << "filelist size:" << filelist.size();
    std::vector<int> filelistint;
    for (auto &file: filelist) {
        Log() << file;
        int fileid = atoi(file.substr(dir.size(), file.size() - 4).c_str());
        Log() << fileid;
        filelistint.push_back(fileid);
    }
    std::sort(filelistint.begin(), filelistint.end());
    mergeTemp(filelistint);
}

inline bool check_if_file_exists(const std::string &name) {
#if defined(_WIN32) || defined(_WIN64)
    return (_access(name.c_str(), 00) != -1);
#endif
#if   defined(__linux__)
    return ( access( name.c_str(), 0 ) != -1 );
#endif


}

void ZQQ_AT_CHINA_VO::CsvRead::readCatalogSingleFileFromTemp(int uniq) {
    if (!check_if_file_exists(temp_dir + std::to_string(uniq) + ".csv")) {
        return;
    }
    inputSingleCsv(temp_dir + std::to_string(uniq) + ".csv", NO_HEADER, -1);
    remove((temp_dir + std::to_string(uniq) + ".csv").c_str());


}

int ZQQ_AT_CHINA_VO::CsvRead::wirteCatalog(int nowlevel, int ipix, std::string Dir) {
    if (!SPLIT_FUNCTION) {
        //Log()<<Log::NOTIFY<<"split function is false,stop split";
        return 0;
    }
    int tot = 0;
    std::string content = "";
    tot += findCatalogtoWrite(content, nowlevel, ipix);
    int uniq = calUNIQ(nowlevel, ipix, nowlevel);
    // Log() << Log::NOTIFY << "write catalog:" << content << "uniq:" << uniq;
    write_to_file(Dir, std::to_string(uniq) + ".csv", content);
    return tot;
}

void ZQQ_AT_CHINA_VO::CsvRead::writeGAIA() {
Log()<<"NO";
}

