//
// Created by zqqia on 1/23/2022.
//

#include "CatalogMerge.h"
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

#include <string>



void ZQQ_AT_CHINA_VO::CatalogMerge::output() {
std::cout<<UserVector[0].getCatlogSize()<<std::endl;
    write_to_file("C:/Archive/","merged.csv",UserVector[0].getCatalog());
}





void ZQQ_AT_CHINA_VO::CatalogMerge::merge() {
    if(UserVector.size()<2){
        std::cout << "file not enough!" << "/n";
        return;
    }

    for(int i=1;i<UserVector.size();i++){
        std::cout<<"start build catalog size:"<<UserVector[0].getCatlogSize()<<std::endl;
        build(UserVector[0].getCatalog());
        std::cout<<" build  success"<<std::endl;
        for(auto l:UserVector[i].getCatalog()){
            kd_query_result closest;
            closest.dis = 360.0;
            closest.id = -1;
            queryKdClosest(l->getRa(), l->getDec(), 1, closest);
            if (closest.id != -1 && closest.dis < SearchThreshold){
                auto &r=UserVector[0].getLineById(closest.id);
                r.setRa((r.getRa()+l->getRa())/2);
                r.setDec((r.getDec()+l->getDec())/2);
                r.content_details.push_back(l->content);
            }else{
                UserVector[0].addLine(l);

            }

        }
    UserVector[i].clear();
    }
    for(auto &row:UserVector[0].getCatalog()){
        CSVRow c;
        c.readString(row->content);
        double magg= atof(std::string(c[5]).c_str());
        for(auto attach:row->content_details){
            c.readString(attach);
            magg+=atof(std::string(c[5]).c_str());
        }
        magg/=(row->content_details.size() + 1);
        row->content.append(","+ std::to_string(magg));
    }
}

void ZQQ_AT_CHINA_VO::CatalogMerge::write_to_file(std::string filepath, std::string filename, std::vector<std::shared_ptr<Catalog::line>> content) {
    if(filepath.size()==0)return;
    if(filepath[filepath.size()-1]!='/'){
        filepath.append("/");
    }
#if defined(_WIN32) || defined(_WIN64)
    if (_access(filepath.c_str(), 00) == -1) {
        _mkdir(filepath.c_str());
    }
#endif
#if   defined(__linux__)
    if (access(filepath.c_str(), 0) == -1) {
            mkdir(filepath.c_str(),S_IRWXU);
        }

#endif
    std::ofstream outFile;
    const size_t bufsize = 256 * 1024;
    char buf[bufsize];
    outFile.rdbuf()->pubsetbuf(buf, bufsize);
    outFile.open((filepath  + filename).c_str(), std::ios_base::app);
    for(auto row:content){
        outFile.write(row->content.c_str(), row->content.size());
    }
    outFile.close();
    //auto endTime = std::chrono::high_resolution_clock::now();

    //std::cout << "Write to " << (filepath  + filename) << " size " << content.size() <<"\n";// " cost time "<< std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms \n";


}
